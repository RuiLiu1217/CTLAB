#include "projModelUtilities.h"
#include "DD3_GPU_Back.h"
#include "CTLAB.h"
#include "Geometry.h"
#include <cassert>
#include "TextureObjectProvider.h"

#define BACK_BLKX 64 ///< Define the block size along x dimension for backprojection
#define BACK_BLKY 4  ///< Define the block size along y dimension for backprojection
#define BACK_BLKZ 1  ///< Define the block size along z dimension for backprojection

static INLINE __device__ double bilerp(int2 v0, int2 v1, int2 v2, int2 v3, double t1, double t2)
{
	double v0_ = __hiloint2double(v0.y, v0.x);
	double v1_ = __hiloint2double(v1.y, v1.x);
	double v2_ = __hiloint2double(v2.y, v2.x);
	double v3_ = __hiloint2double(v3.y, v3.x);

	double vv0 = v0_ * (1.0 - t1) + v1_ * t1;
	double vv1 = v2_ * (1.0 - t1) + v3_ * t1;
	return vv0 * (1 - t2) + vv1 * t2;
}

/// \brief BackProjectionMethod with different model/implementation
enum BackProjectionMethod{
	_BRANCHLESS, ///< Branchless DD model
	_PSEUDODD,  ///< Pseudo DD model
	_ZLINEBRANCHLESS,  ///< Branchless DD model with Z-line implementation
	_VOLUMERENDERING  ///< Pixel driven model
};

#ifndef CALDETPARAS
#define CALDETPARAS
/// \brief Calculate the size of the detector fan angle, detector cell height and the center of the index according to the x,y coordinates of the detector cells at initial position
/// \param xds x coordinate
/// \param yds y coordinate
/// \param zds z coordinate
/// \param x0 source position x
/// \param y0 source position y
/// \param z0 source position z
/// \param DNU detector cell number along channel direction
/// \param DNV detector cell number along bench moving direction
static float4 calDetParas(float* xds, float* yds, float* zds, float x0, float y0, float z0, int DNU, int DNV)
{
	float* bxds = new float[DNU + 1];
	float* byds = new float[DNU + 1];
	float* bzds = new float[DNV + 1];
	DD3Boundaries(DNU + 1, xds, bxds);
	DD3Boundaries(DNU + 1, yds, byds);
	DD3Boundaries(DNV + 1, zds, bzds);

	float ddv = (bzds[DNV] - bzds[0]) / DNV;
	float detCtrIdxV = (-(bzds[0] - z0) / ddv) - 0.5;
	float2 dir = normalize(make_float2(-x0, -y0));
	float2 dirL = normalize(make_float2(bxds[0] - x0, byds[0] - y0));
	float2 dirR = normalize(make_float2(bxds[DNU] - x0, byds[DNU] - y0));
	float dbeta = asin(dirL.x * dirR.y - dirL.y * dirR.x) / DNU;
	assert(dbeta != 0);
	float minBeta = asin(dir.x * dirL.y - dir.y * dirL.x);
	float detCtrIdxU = -minBeta / dbeta - 0.5;
	delete [] bxds;
	delete [] byds;
	delete [] bzds;
	return make_float4(detCtrIdxU, detCtrIdxV, dbeta, ddv);
}
#endif


/// \brief Add two columns and rows of all 0 borders on the left and top of each projection views
/// \param prjIn Input of the projection data
/// \param prjOut Output of the projection data
/// \param DNU detector cell number along channel direction
/// \param DNV detector cell number along bench moving direction
/// \param PN number of views
__global__ void addTwoSidedZeroBoarder(float* prjIn, float* prjOut,
	const int DNU, const int DNV, const int PN)
{
	int idv = threadIdx.x + blockIdx.x * blockDim.x;
	int idu = threadIdx.y + blockIdx.y * blockDim.y;
	int pn = threadIdx.z + blockIdx.z * blockDim.z;
	if (idu < DNU && idv < DNV && pn < PN)
	{
		int inIdx = (pn * DNU + idu) * DNV + idv;
		int outIdx = (pn * (DNU + 2) + (idu + 1)) * (DNV + 2) + idv + 1;
		prjOut[outIdx] = prjIn[inIdx];
	}
}

/// \brief Add one column and row of all 0 borders on the left and top of each projection views
/// \param prj_in Input of the projection data
/// \param prj_out Output of the projection data
/// \param DNU detector cell number along channel direction
/// \param DNV detector cell number along bench moving direction
/// \param PN number of views
template<typename T>
__global__ void addOneSidedZeroBoarder(const T* prj_in, T* prj_out, int DNU, int DNV, int PN)
{
	int idv = threadIdx.x + blockIdx.x * blockDim.x;
	int idu = threadIdx.y + blockIdx.y * blockDim.y;
	int pn = threadIdx.z + blockIdx.z * blockDim.z;
	if (idu < DNU && idv < DNV && pn < PN)
	{
		int i = (pn * DNU + idu) * DNV + idv;
		int ni = (pn * (DNU + 1) + (idu + 1)) * (DNV + 1) + idv + 1;
		prj_out[ni] = prj_in[i];
	}
}


/// \brief Integral image along vertical direction
/// \param prj Input/Output of the projection data
/// \param ZN image dimension along bench moving direction
/// \param N usually equals (XN x YN)
template<typename T>
__global__ void verticalIntegral2(T* prj, int ZN, int N) {
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < N) {
		int currentHead = idx * ZN;
		for (int ii = 1; ii < ZN; ++ii)	{
			prj[currentHead + ii] += prj[currentHead + ii - 1];
		}
	}
}

/// \brief Integral image along horizontal direction
/// \param prj Input/Output of the projection data
/// \param DNU detector cell number along channel direction
/// \param DNV detector cell number along bench moving direction
/// \param PN number of views
__global__ void heorizontalIntegral2(float* prj, int DNU, int DNV, int PN)
{
	int idv = threadIdx.x + blockIdx.x * blockDim.x;
	int pIdx = threadIdx.y + blockIdx.y * blockDim.y;
	if (idv < DNV && pIdx < PN)
	{
		int headPrt = pIdx * DNU * DNV + idv;
		for (int ii = 1; ii < DNU; ++ii)
		{
			prj[headPrt + ii * DNV] = prj[headPrt + ii * DNV] + prj[headPrt + (ii - 1) * DNV];
		}
	}
}

/// \brief Generate the summed area table for the projection data
/// \param hprj projection data stored in host memory
/// \param DNU detector cell number along channel direction
/// \param DNV detector cell number along bench moving direction
/// \param PN number of views
static thrust::device_vector<float> genSAT_of_Projection(
	float* hprj,
	int DNU, int DNV, int PN)
{
	const int siz = DNU * DNV * PN;
	const int nsiz = (DNU + 1) * (DNV + 1) * PN;
	thrust::device_vector<float> prjSAT(nsiz, 0);
	thrust::device_vector<float> prj(hprj, hprj + siz);
	dim3 copyBlk(64, 16, 1);
	dim3 copyGid(
		(DNV + copyBlk.x - 1) / copyBlk.x,
		(DNU + copyBlk.y - 1) / copyBlk.y,
		(PN + copyBlk.z - 1) / copyBlk.z);

	addOneSidedZeroBoarder << <copyGid, copyBlk >> >(
		thrust::raw_pointer_cast(&prj[0]),
		thrust::raw_pointer_cast(&prjSAT[0]),
		DNU, DNV, PN);
	const int nDNU = DNU + 1;
	const int nDNV = DNV + 1;

	copyBlk.x = 512;
	copyBlk.y = 1;
	copyBlk.z = 1;
	copyGid.x = (nDNU * PN + copyBlk.x - 1) / copyBlk.x;
	copyGid.y = 1;
	copyGid.z = 1;
	verticalIntegral2 << <copyGid, copyBlk >> >(
		thrust::raw_pointer_cast(&prjSAT[0]),
		nDNV, nDNU * PN);
	copyBlk.x = 64;
	copyBlk.y = 16;
	copyBlk.z = 1;
	copyGid.x = (nDNV + copyBlk.x - 1) / copyBlk.x;
	copyGid.y = (PN + copyBlk.y - 1) / copyBlk.y;
	copyGid.z = 1;


	heorizontalIntegral2 << <copyGid, copyBlk >> >(
		thrust::raw_pointer_cast(&prjSAT[0]),
		nDNU, nDNV, PN);

	return prjSAT;
}


template < BackProjectionMethod METHOD >
__global__ void DD3_gpu_back_ker(
	cudaTextureObject_t prjTexObj,
	float* vol,
	const byte* __restrict__ msk,
	const float3* __restrict__ cossinT,
	float3 s,
	float S2D,
	float3 curvox,
	float dx, float dz,
	float dbeta, float ddv,
	float2 detCntIdx,
	int3 VN,
	int PN, int squared)
{}


template < BackProjectionMethod METHOD >
__global__ void DD3_panel_gpu_back_ker(
	cudaTextureObject_t prjTexObj,
	float* vol,
	const byte* __restrict__ msk,
	const float3* __restrict__ cossinT,
	float3 s,
	float S2D,
	float3 curvox,
	float dx, float dz,
	float dbeta, float ddv,
	float2 detCntIdx,
	int3 VN,
	int PN, int squared)
{}


template<>
__global__ void DD3_panel_gpu_back_ker<_BRANCHLESS>
	(cudaTextureObject_t prjTexObj,
	float* vol,
	const byte* __restrict__ msk,
	const float3* __restrict__ cossinT,
	float3 s,
	float S2D,
	float3 curvox,
	float dx, float dz,
	float dbeta, //Detector size in channel direction
	float ddv, //detector size in Z direction
	float2 detCntIdx,
	int3 VN,
	int PN, int squared)
{
	int3 id;
	id.z = threadIdx.x + __umul24(blockIdx.x, blockDim.x);
	id.x = threadIdx.y + __umul24(blockIdx.y, blockDim.y);
	id.y = threadIdx.z + __umul24(blockIdx.z, blockDim.z);
	if (id.x < VN.x && id.y < VN.y && id.z < VN.z)
	{
		if (msk[id.y * VN.x + id.x] != 1)
			return;
		curvox = (id - curvox) * make_float3(dx, dx, dz);
		float3 cursour;
		float idxL, idxR, idxU, idxD;
		float cosVal;
		float summ = 0;

		float3 cossin;
		assert(s.x != 0 || s.y != 0);
		float inv_sid = 1.0 / sqrtf(s.x * s.x + s.y * s.y);
		float3 dir;
		float l_square;
		float l;
		float alpha;
		float S2D2 = S2D;
		S2D = S2D / ddv;
		dbeta = 1.0 / dbeta;
		dz = dz * 0.5;
		float2 dirL;
		float2 dirR;
		for (int angIdx = 0; angIdx < PN; ++angIdx)
		{
			cossin = cossinT[angIdx];
			cursour = make_float3(
				s.x * cossin.x - s.y * cossin.y,
				s.x * cossin.y + s.y * cossin.x,
				s.z + cossin.z);

			dir = curvox - cursour;
			l_square = dir.x * dir.x + dir.y * dir.y;
			l = rsqrtf(l_square);
			idxU = (dir.z + dz) * S2D * l + detCntIdx.y + 1;
			idxD = (dir.z - dz) * S2D * l + detCntIdx.y + 1;

			if (fabsf(cursour.x) > fabsf(cursour.y))
			{
				ddv = dir.x;
				dirL = normalize(make_float2(dir.x, dir.y - 0.5 * dx));
				dirR = normalize(make_float2(dir.x, dir.y + 0.5 * dx));
			}
			else
			{
				ddv = dir.y;
				dirL = normalize(make_float2(dir.x + 0.5 * dx, dir.y));
				dirR = normalize(make_float2(dir.x - 0.5 * dx, dir.y));
			}
			cosVal = dx / ddv * sqrtf(l_square + dir.z * dir.z);

			//! TODO: Test the correctness of this method.
			alpha = asinf((cursour.y * dirL.x - cursour.x * dirL.y) * inv_sid);
			idxL = (tanf(alpha) * S2D2 * dbeta) + detCntIdx.x + 1;
			alpha = asinf((cursour.y * dirR.x - cursour.x * dirR.y) * inv_sid);
			idxR = (tanf(alpha) * S2D2 * dbeta) + detCntIdx.x + 1;
			//summ += idxL;

			summ +=
				(-tex3D<float>(prjTexObj, idxD, idxR, angIdx + 0.5)
				- tex3D<float>(prjTexObj, idxU, idxL, angIdx + 0.5)
				+ tex3D<float>(prjTexObj, idxD, idxL, angIdx + 0.5)
				+ tex3D<float>(prjTexObj, idxU, idxR, angIdx + 0.5)) * cosVal;
		}
		__syncthreads();
		vol[__umul24((__umul24(id.y, VN.x) + id.x), VN.z) + id.z] = summ;
	}
}

template<>
__global__ void DD3_gpu_back_ker<_BRANCHLESS>(
	cudaTextureObject_t prjTexObj,
	float* vol,
	const byte* __restrict__ msk,
	const float3* __restrict__ cossinT,
	float3 s,
	float S2D,
	float3 curvox,
	float dx, float dz,
	float dbeta, float ddv,
	float2 detCntIdx,
	int3 VN,
	int PN, int squared)
{
	int3 id;
	id.z = threadIdx.x + __umul24(blockIdx.x, blockDim.x);
	id.x = threadIdx.y + __umul24(blockIdx.y, blockDim.y);
	id.y = threadIdx.z + __umul24(blockIdx.z, blockDim.z);
	if (id.x < VN.x && id.y < VN.y && id.z < VN.z)
	{
		if (msk[id.y * VN.x + id.x] != 1)
			return;
		curvox = (id - curvox) * make_float3(dx, dx, dz);
		float3 cursour;
		float idxL, idxR, idxU, idxD;
		float cosVal;
		float summ = 0;

		float3 cossin;
		assert(s.x != 0 || s.y != 0);
		float inv_sid = 1.0 / sqrtf(s.x * s.x + s.y * s.y);
		float3 dir;
		float l_square;
		float l;
		float alpha;
		float deltaAlpha;
		S2D = S2D / ddv;
		dbeta = 1.0 / dbeta;
		dz = dz * 0.5;
		for (int angIdx = 0; angIdx < PN; ++angIdx)
		{
			cossin = cossinT[angIdx];
			cursour = make_float3(
				s.x * cossin.x - s.y * cossin.y,
				s.x * cossin.y + s.y * cossin.x,
				s.z + cossin.z);

			dir = curvox - cursour;
			l_square = dir.x * dir.x + dir.y * dir.y;
			assert(l_square != 0);
			l = rsqrtf(l_square);
			idxU = (dir.z + dz) * S2D * l + detCntIdx.y + 1;
			idxD = (dir.z - dz) * S2D * l + detCntIdx.y + 1;

			alpha = asinf((cursour.y * dir.x - cursour.x * dir.y) * inv_sid * l);
			if (fabsf(cursour.x) > fabsf(cursour.y))
			{
				ddv = dir.x;
			}
			else
			{
				ddv = dir.y;
			}
			
			deltaAlpha = ddv / l_square * dx * 0.5;
			cosVal = dx / ddv * sqrtf(l_square + dir.z * dir.z);
			idxL = (alpha - deltaAlpha) * dbeta + detCntIdx.x + 1;
			idxR = (alpha + deltaAlpha) * dbeta + detCntIdx.x + 1;

			summ +=
				(-tex3D<float>(prjTexObj, idxD, idxR, angIdx + 0.5)
				- tex3D<float>(prjTexObj, idxU, idxL, angIdx + 0.5)
				+ tex3D<float>(prjTexObj, idxD, idxL, angIdx + 0.5)
				+ tex3D<float>(prjTexObj, idxU, idxR, angIdx + 0.5)) * cosVal;
		}
		__syncthreads();
		vol[__umul24((__umul24(id.y, VN.x) + id.x), VN.z) + id.z] = summ;
	}
}


__global__ void DD3_gpu_back_branchless_ker(
	cudaTextureObject_t prjTexObj,
	float* vol,
	const byte* __restrict__ msk,
	const float3* __restrict__ cossinT,
	float3 s,
	float S2D,
	float3 curvox,
	float dx, float dz,
	float dbeta, float ddv,
	float2 detCntIdx,
	int3 VN,
	int PN, int squared)
{
	int3 id;
	id.z = threadIdx.x + __umul24(blockIdx.x, blockDim.x);
	id.x = threadIdx.y + __umul24(blockIdx.y, blockDim.y);
	id.y = threadIdx.z + __umul24(blockIdx.z, blockDim.z);
	if (id.x < VN.x && id.y < VN.y && id.z < VN.z)
	{
		if (msk[id.y * VN.x + id.x] != 1)
			return;
		curvox = (id - curvox) * make_float3(dx, dx, dz);
		float3 cursour;
		float idxL, idxR, idxU, idxD;
		float cosVal;
		float summ = 0;

		float3 cossin;
		assert(s.x != 0 || s.y != 0);
		float inv_sid = 1.0 / sqrtf(s.x * s.x + s.y * s.y);
		float3 dir;
		float l_square;
		float l;
		float alpha;
		float deltaAlpha;
		S2D = S2D / ddv;
		dbeta = 1.0 / dbeta;
		dz = dz * 0.5;
		for (int angIdx = 0; angIdx < PN; ++angIdx)
		{
			cossin = cossinT[angIdx];
			cursour = make_float3(
				s.x * cossin.x - s.y * cossin.y,
				s.x * cossin.y + s.y * cossin.x,
				s.z + cossin.z);

			dir = curvox - cursour;
			l_square = dir.x * dir.x + dir.y * dir.y;
			assert(l_square != 0);
			l = rsqrtf(l_square);
			idxU = (dir.z + dz) * S2D * l + detCntIdx.y + 1;
			idxD = (dir.z - dz) * S2D * l + detCntIdx.y + 1;

			alpha = asinf((cursour.y * dir.x - cursour.x * dir.y) * inv_sid * l);
			if (fabsf(cursour.x) > fabsf(cursour.y))
			{
				ddv = dir.x;
			}
			else
			{
				ddv = dir.y;
			}
			deltaAlpha = ddv / l_square * dx * 0.5;
			cosVal = dx / ddv * sqrtf(l_square + dir.z * dir.z);
			idxL = (alpha - deltaAlpha) * dbeta + detCntIdx.x + 1;
			idxR = (alpha + deltaAlpha) * dbeta + detCntIdx.x + 1;

			summ +=
				(-tex3D<float>(prjTexObj, idxD, idxR, angIdx + 0.5)
				- tex3D<float>(prjTexObj, idxU, idxL, angIdx + 0.5)
				+ tex3D<float>(prjTexObj, idxD, idxL, angIdx + 0.5)
				+ tex3D<float>(prjTexObj, idxU, idxR, angIdx + 0.5)) * cosVal;
		}
		__syncthreads();
		vol[__umul24((__umul24(id.y, VN.x) + id.x), VN.z) + id.z] = summ;
	}
}

__global__ void DD3_gpu_back_volumerendering_ker(
	cudaTextureObject_t texObj,
	float* vol,
	const byte* __restrict__ msk,
	const float3* __restrict__ cossinZT,
	float3 s,
	float S2D,
	float3 objCntIdx,
	float dx, float dz, float dbeta, float ddv,
	float2 detCntIdx, int3 VN, int PN, int squared)
{
	int3 id;
	id.z = __mul24(blockDim.x, blockIdx.x) + threadIdx.x;
	id.x = __mul24(blockDim.y, blockIdx.y) + threadIdx.y;
	id.y = __mul24(blockDim.z, blockIdx.z) + threadIdx.z;
	if (id.x < VN.x && id.y < VN.y && id.z < VN.z)
	{
		if (msk[id.y * VN.x + id.x] != 1)
			return;
		float3 curVox = make_float3(
			(id.x - objCntIdx.x) * dx,
			(id.y - objCntIdx.y) * dx,
			(id.z - objCntIdx.z) * dz);
		float3 boxmin = curVox - make_float3(dx, dx, dz) * 0.5;
		float3 boxmax = curVox + make_float3(dx, dx, dz) * 0.5;
		float3 dir;
		float3 cursour;
		float tnear, tfar;
		float intersectLength;
		float summ = 0;
		float sid = sqrtf(s.x * s.x + s.y * s.y);
		float idxZ;
		float idxXY;
		float3 cossinT;
		float2 ds;
		float l;
		float alpha;
		for (int angIdx = 0; angIdx < PN; ++angIdx)
		{
			cossinT = cossinZT[angIdx];
			cursour = make_float3(
				s.x * cossinT.x - s.y * cossinT.y,
				s.x * cossinT.y + s.y * cossinT.x,
				s.z + cossinT.z);
			dir = curVox - cursour;
			ds = make_float2(-cursour.x, -cursour.y);
			l = sqrtf(dir.x * dir.x + dir.y * dir.y);
			assert(l != 0);
			idxZ = dir.z * S2D / l / ddv + detCntIdx.y + 0.5;
			alpha = asinf((ds.x * dir.y - ds.y * dir.x) / (l * sid));
			dir = normalize(dir);
			intersectLength = intersectBox(cursour, dir, boxmin, boxmax, &tnear, &tfar);
			intersectLength *= (tfar - tnear);
			idxXY = alpha / dbeta + detCntIdx.x + 0.5;
			summ += tex3D<float>(texObj, idxZ, idxXY, angIdx + 0.5f) * intersectLength;

		}
		vol[(id.y + VN.x + id.x) * VN.z + id.z] = summ;
	}
}




template<>
__global__ void DD3_gpu_back_ker<_VOLUMERENDERING>(
	cudaTextureObject_t texObj,
	float* vol,
	const byte* __restrict__ msk,
	const float3* __restrict__ cossinZT,
	float3 s,
	float S2D,
	float3 objCntIdx,
	float dx, float dz, float dbeta, float ddv,
	float2 detCntIdx, int3 VN, int PN, int squared)
{
	int3 id;
	id.z = __mul24(blockDim.x, blockIdx.x) + threadIdx.x;
	id.x = __mul24(blockDim.y, blockIdx.y) + threadIdx.y;
	id.y = __mul24(blockDim.z, blockIdx.z) + threadIdx.z;
	if (id.x < VN.x && id.y < VN.y && id.z < VN.z)
	{
		if (msk[id.y * VN.x + id.x] != 1)
			return;
		float3 curVox = make_float3(
			(id.x - objCntIdx.x) * dx,
			(id.y - objCntIdx.y) * dx,
			(id.z - objCntIdx.z) * dz);
		float3 boxmin = curVox - make_float3(dx, dx, dz) * 0.5;
		float3 boxmax = curVox + make_float3(dx, dx, dz) * 0.5;
		float3 dir;
		float3 cursour;
		float tnear, tfar;
		float intersectLength;
		float summ = 0;
		float sid = sqrtf(s.x * s.x + s.y * s.y);
		float idxZ;
		float idxXY;
		float3 cossinT;
		float2 ds;
		float l;
		float alpha;
		for (int angIdx = 0; angIdx < PN; ++angIdx)
		{
			cossinT = cossinZT[angIdx];
			cursour = make_float3(
				s.x * cossinT.x - s.y * cossinT.y,
				s.x * cossinT.y + s.y * cossinT.x,
				s.z + cossinT.z);
			dir = curVox - cursour;
			ds = make_float2(-cursour.x, -cursour.y);
			l = sqrtf(dir.x * dir.x + dir.y * dir.y);
			assert(l != 0);
			idxZ = dir.z * S2D / l / ddv + detCntIdx.y + 0.5;
			alpha = asinf((ds.x * dir.y - ds.y * dir.x) / (l * sid));
			dir = normalize(dir);
			intersectLength = intersectBox(cursour, dir, boxmin, boxmax, &tnear, &tfar);
			intersectLength *= (tfar - tnear);
			idxXY = alpha / dbeta + detCntIdx.x + 0.5;
			summ += tex3D<float>(texObj, idxZ, idxXY, angIdx + 0.5f) * intersectLength;

		}
		vol[(id.y + VN.x + id.x) * VN.z + id.z] = summ;
	}
}




static
void DD3_gpu_back_volumerendering(float x0, float y0, float z0,
	int DNU, int DNV, float* xds,
	float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN,
	float* hvol, float* hprj,
	float dx, float dz, byte* mask, int squared, int gpunum)
{
	cudaSetDevice(gpunum);
	cudaDeviceReset();
	const int TOTVON = XN * YN * ZN;
	const float objCntIdxX = (XN - 1.0) * 0.5 - imgXCenter / dx;
	const float objCntIdxY = (YN - 1.0) * 0.5 - imgYCenter / dx;
	const float objCntIdxZ = (ZN - 1.0) * 0.5 - imgZCenter / dz;

	thrust::device_vector<float> vol(TOTVON, 0);
	thrust::device_vector<float3> cursour(PN);
	thrust::device_vector<float3> cossinZT(PN);
	thrust::device_vector<float2> dirsour(PN);

	thrust::device_vector<float> angs(hangs, hangs + PN);
	thrust::device_vector<float> zPos(hzPos, hzPos + PN);
	thrust::device_vector<byte> d_msk(mask, mask + XN * YN);

	thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(angs.begin(), zPos.begin())),
		thrust::make_zip_iterator(thrust::make_tuple(angs.end(), zPos.end())),
		thrust::make_zip_iterator(thrust::make_tuple(
		cossinZT.begin(), cursour.begin(), dirsour.begin())),
		CTMBIR::ConstantForBackProjection3(x0, y0, z0));

	float* bxds = new float[DNU + 1];
	float* byds = new float[DNU + 1];
	float* bzds = new float[DNV + 1];
	DD3Boundaries(DNU + 1, xds, bxds);
	DD3Boundaries(DNU + 1, yds, byds);
	DD3Boundaries(DNV + 1, zds, bzds);

	float ddv = (bzds[DNV] - bzds[0]) / DNV;
	float detCtrIdxV = (-(bzds[0] - z0) / ddv) - 0.5;

	float2 dir = normalize(make_float2(-x0, -y0));
	float2 dirL = normalize(make_float2(bxds[0] - x0, byds[0] - y0));
	float2 dirR = normalize(make_float2(bxds[DNU] - x0, byds[DNU] - y0));
	float dbeta = asin(dirL.x * dirR.y - dirL.y * dirR.x) / DNU;
	float minBeta = asin(dir.x * dirL.y - dir.y * dirL.x);
	float detCtrIdxU = -minBeta / dbeta - 0.5;
	const float S2D = hypotf(xds[0] - x0, yds[0] - y0);

	cudaExtent prjSize;
	prjSize.width = DNV;
	prjSize.height = DNU;
	prjSize.depth = PN;

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	cudaArray *d_prjArray;
	cudaMalloc3DArray(&d_prjArray, &channelDesc, prjSize);
	cudaMemcpy3DParms copyParams = {0};
	copyParams.srcPtr = make_cudaPitchedPtr(
		(void*) hprj, prjSize.width * sizeof(float),
		prjSize.width, prjSize.height);
	copyParams.dstArray = d_prjArray;
	copyParams.extent = prjSize;
	copyParams.kind = cudaMemcpyHostToDevice;
	cudaMemcpy3D(&copyParams);
	cudaResourceDesc resDesc;
	memset(&resDesc, 0, sizeof(resDesc));
	resDesc.resType = cudaResourceTypeArray;
	resDesc.res.array.array = d_prjArray;
	cudaTextureDesc texDesc;
	memset(&texDesc, 0, sizeof(texDesc));
	texDesc.addressMode[0] = cudaAddressModeBorder;
	texDesc.addressMode[1] = cudaAddressModeBorder;
	texDesc.addressMode[2] = cudaAddressModeBorder;
	texDesc.filterMode = cudaFilterModeLinear;
	texDesc.readMode = cudaReadModeElementType;
	cudaTextureObject_t texObj;
	texDesc.normalizedCoords = false;
	cudaCreateTextureObject(&texObj, &resDesc, &texDesc, nullptr);

	dim3 blk(BACK_BLKX, BACK_BLKY, BACK_BLKZ);
	dim3 gid(
		(ZN + blk.x - 1) / blk.x,
		(XN + blk.y - 1) / blk.y,
		(ZN + blk.z - 1) / blk.z);

	DD3_gpu_back_volumerendering_ker << <gid, blk >> >
		(texObj,
		thrust::raw_pointer_cast(&vol[0]),
		thrust::raw_pointer_cast(&d_msk[0]),
		thrust::raw_pointer_cast(&cossinZT[0]),
		make_float3(x0, y0, z0),
		S2D,
		make_float3(objCntIdxX, objCntIdxY, objCntIdxZ),
		dx, dz, dbeta, ddv,
		make_float2(detCtrIdxU, detCtrIdxV),
		make_int3(XN, YN, ZN), PN, squared);


	thrust::copy(vol.begin(), vol.end(), hvol);
	delete [] bxds;
	delete [] byds;
	delete [] bzds;
}



__global__ void DD3_gpu_back_pseudodistancedriven_ker(
	cudaTextureObject_t texObj,
	float* vol,
	const byte* __restrict__ msk,
	const float3* __restrict__ cossinZT,
	float3 s,
	float S2D,
	float3 objCntIdx,
	float dx, float dz, float dbeta, float ddv,
	float2 detCntIdx,
	int3 VN, int PN, int squared)
{
	int k = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;
	int i = __mul24(blockIdx.y, blockDim.y) + threadIdx.y;
	int j = __mul24(blockIdx.z, blockDim.z) + threadIdx.z;
	if (i < VN.x && j < VN.y && k < VN.z)
	{
		if (msk[j * VN.x + i] != 1)
			return;
		float3 curVox = make_float3(
			(i - objCntIdx.x) * dx,
			(j - objCntIdx.y) * dx,
			(k - objCntIdx.z) * dz);

		float3 dir;
		float3 cursour;
		float invsid = rsqrtf(s.x * s.x + s.y * s.y);
		float invl;
		float idxZ;
		float idxXY;
		float alpha;
		float cosVal;
		float3 cossinT;
		float summ = 0;
		S2D = S2D / ddv;
		dbeta = 1.0 / dbeta;
		for (int angIdx = 0; angIdx != PN; ++angIdx)
		{
			cossinT = cossinZT[angIdx];
			cursour = make_float3(
				s.x * cossinT.x - s.y * cossinT.y,
				s.x * cossinT.y + s.y * cossinT.x,
				s.z + cossinT.z);

			dir = curVox - cursour;
			ddv = dir.x * dir.x + dir.y * dir.y;
			invl = rsqrtf(ddv);
			idxZ = dir.z * S2D * invl + detCntIdx.y + 0.5;
			alpha = asinf((cursour.y * dir.x - cursour.x * dir.y) * invl * invsid);
			if (fabsf(cursour.x) >= fabsf(cursour.y))
			{
				assert(dir.x != 0);
				cosVal = fabsf(1.0 / dir.x);
			}
			else
			{
				assert(dir.y != 0);
				cosVal = fabsf(1.0 / dir.y);
			}
			cosVal *= (dx * sqrtf(ddv + dir.z * dir.z));
			idxXY = alpha * dbeta + detCntIdx.x + 0.5;
			summ += tex3D<float>(texObj, idxZ, idxXY, angIdx + 0.5f) * cosVal;
		}
		__syncthreads();
		vol[(j * VN.x + i) * VN.z + k] = summ;
	}
}









template<>
__global__ void DD3_gpu_back_ker<_PSEUDODD>(
	cudaTextureObject_t texObj,
	float* vol,
	const byte* __restrict__ msk,
	const float3* __restrict__ cossinZT,
	float3 s,
	float S2D,
	float3 objCntIdx,
	float dx, float dz, float dbeta, float ddv,
	float2 detCntIdx,
	int3 VN, int PN, int squared)
{
	int k = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;
	int i = __mul24(blockIdx.y, blockDim.y) + threadIdx.y;
	int j = __mul24(blockIdx.z, blockDim.z) + threadIdx.z;
	if (i < VN.x && j < VN.y && k < VN.z)
	{
		if (msk[j * VN.x + i] != 1)
			return;
		float3 curVox = make_float3(
			(i - objCntIdx.x) * dx,
			(j - objCntIdx.y) * dx,
			(k - objCntIdx.z) * dz);

		float3 dir;
		float3 cursour;
		float invsid = rsqrtf(s.x * s.x + s.y * s.y);
		float invl;
		float idxZ;
		float idxXY;
		float alpha;
		float cosVal;
		float3 cossinT;
		float summ = 0;
		S2D = S2D / ddv;
		dbeta = 1.0 / dbeta;
		for (int angIdx = 0; angIdx != PN; ++angIdx)
		{
			cossinT = cossinZT[angIdx];
			cursour = make_float3(
				s.x * cossinT.x - s.y * cossinT.y,
				s.x * cossinT.y + s.y * cossinT.x,
				s.z + cossinT.z);

			dir = curVox - cursour;
			ddv = dir.x * dir.x + dir.y * dir.y;
			invl = rsqrtf(ddv);
			idxZ = dir.z * S2D * invl + detCntIdx.y + 0.5;
			alpha = asinf((cursour.y * dir.x - cursour.x * dir.y) * invl * invsid);
			if (fabsf(cursour.x) >= fabsf(cursour.y))
			{
				assert(dir.x != 0);
				cosVal = fabsf(1.0 / dir.x);
			}
			else
			{
				assert(dir.y != 0);
				cosVal = fabsf(1.0 / dir.y);
			}
			cosVal *= (dx * sqrtf(ddv + dir.z * dir.z));
			idxXY = alpha * dbeta + detCntIdx.x + 0.5;
			summ += tex3D<float>(texObj, idxZ, idxXY, angIdx + 0.5f) * cosVal;
		}
		__syncthreads();
		vol[(j * VN.x + i) * VN.z + k] = summ;
	}
}




static
void DD3_gpu_back_pseudodistancedriven(
	float x0, float y0, float z0,
	int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN,
	float* hvol, float* hprj,
	float dx, float dz,
	byte* mask, int squared, int gpunum)
{
	cudaSetDevice(gpunum);
	cudaDeviceReset();
	const int TOTVON = XN  * YN * ZN;
	const float objCntIdxX = (XN - 1.0) * 0.5 - imgXCenter / dx;
	const float objCntIdxY = (YN - 1.0) * 0.5 - imgYCenter / dx;
	const float objCntIdxZ = (ZN - 1.0) * 0.5 - imgZCenter / dz;

	thrust::device_vector<byte> d_msk(mask, mask + XN * YN);
	thrust::device_vector<float> vol(TOTVON, 0);
	thrust::device_vector<float3> cursour(PN);
	thrust::device_vector<float2> dirsour(PN);
	thrust::device_vector<float3> cossinZT(PN);

	thrust::device_vector<float> angs(hangs, hangs + PN);
	thrust::device_vector<float> zPos(hzPos, hzPos + PN);

	thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(angs.begin(), zPos.begin())),
		thrust::make_zip_iterator(thrust::make_tuple(angs.end(), zPos.end())),
		thrust::make_zip_iterator(thrust::make_tuple(cossinZT.begin(), cursour.begin(), dirsour.begin())),
		CTMBIR::ConstantForBackProjection3(x0, y0, z0));

	float* bxds = new float[DNU + 1];
	float* byds = new float[DNU + 1];
	float* bzds = new float[DNV + 1];
	DD3Boundaries(DNU + 1, xds, bxds);
	DD3Boundaries(DNU + 1, yds, byds);
	DD3Boundaries(DNV + 1, zds, bzds);

	float ddv = (bzds[DNV] - bzds[0]) / DNV;
	float detCtrIdxV = (-(bzds[0] - z0) / ddv) - 0.5;

	float2 dir = normalize(make_float2(-x0, -y0));
	float2 dirL = normalize(make_float2(bxds[0] - x0, byds[0] - y0));
	float2 dirR = normalize(make_float2(bxds[DNU] - x0, byds[DNU] - y0));
	float dbeta = asin(dirL.x * dirR.y - dirL.y * dirR.x) / DNU;
	assert(dbeta != 0);
	float minBeta = asin(dir.x * dirL.y - dir.y * dirL.x);
	float detCtrIdxU = -minBeta / dbeta - 0.5;
	const float S2D = hypotf(xds[0] - x0, yds[0] - y0);

	cudaExtent prjSize;
	prjSize.width = DNV;
	prjSize.height = DNU;
	prjSize.depth = PN;

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	cudaArray *d_prjArray;
	cudaMalloc3DArray(&d_prjArray, &channelDesc, prjSize);
	cudaMemcpy3DParms copyParams = { 0 };
	copyParams.srcPtr = make_cudaPitchedPtr(
		(void*) hprj, prjSize.width * sizeof(float),
		prjSize.width, prjSize.height);
	copyParams.dstArray = d_prjArray;
	copyParams.extent = prjSize;
	copyParams.kind = cudaMemcpyHostToDevice;
	cudaMemcpy3D(&copyParams);
	cudaResourceDesc resDesc;
	memset(&resDesc, 0, sizeof(resDesc));
	resDesc.resType = cudaResourceTypeArray;
	resDesc.res.array.array = d_prjArray;
	cudaTextureDesc texDesc;
	memset(&texDesc, 0, sizeof(texDesc));
	texDesc.addressMode[0] = cudaAddressModeBorder;
	texDesc.addressMode[1] = cudaAddressModeBorder;
	texDesc.addressMode[2] = cudaAddressModeBorder;
	texDesc.filterMode = cudaFilterModeLinear;
	texDesc.readMode = cudaReadModeElementType;
	cudaTextureObject_t texObj;
	texDesc.normalizedCoords = false;
	cudaCreateTextureObject(&texObj, &resDesc, &texDesc, nullptr);

	dim3 blk(BACK_BLKX, BACK_BLKY, BACK_BLKZ);
	dim3 gid(
		(ZN + blk.x - 1) / blk.x,
		(XN + blk.y - 1) / blk.y,
		(YN + blk.z - 1) / blk.z);

	DD3_gpu_back_pseudodistancedriven_ker << <gid, blk >> >(texObj,
		thrust::raw_pointer_cast(&vol[0]),
		thrust::raw_pointer_cast(&d_msk[0]),
		thrust::raw_pointer_cast(&cossinZT[0]),
		make_float3(x0, y0, z0),
		S2D,
		make_float3(objCntIdxX, objCntIdxY, objCntIdxZ),
		dx, dz, dbeta, ddv,
		make_float2(detCtrIdxU, detCtrIdxV),
		make_int3(XN, YN, ZN),
		PN, static_cast<int>(squared != 0));
	thrust::copy(vol.begin(), vol.end(), hvol);
	delete [] bxds;
	delete [] byds;
	delete [] bzds;
}



template<int LAYERS>
__global__ void DD3_gpu_back_zlinebasedbranchless_ker(
	cudaTextureObject_t prjTexObj,
	float* vol,
	const byte* __restrict__ msk,
	const float3* __restrict__ cossinZT,
	float3 s,
	float S2D,
	float3 objCntIdx,
	float dx, float dz, float dbeta, float ddv,
	float2 detCntIdx, int3 VN, int PN, int squared)
{
	int idx = threadIdx.x + __umul24(blockIdx.x, blockDim.x);
	int idy = threadIdx.y + __umul24(blockIdx.y, blockDim.y);
	__shared__ float summ[4][8][LAYERS + 1];
#pragma unroll
	for (int i = 0; i <= LAYERS; ++i)
	{
		summ[threadIdx.y][threadIdx.x][i] = 0;
	}
	__syncthreads();
	if (idx < VN.x && idy < VN.y)
	{
		if (msk[idy * VN.x + idx] != 1)
			return;
		float curang(0);
		float2 dirlft, dirrgh;
		float3 cursour;
		float idxL, idxR, idxD;
		float cosVal = 1.0;
		float2 curvox_xy = make_float2((idx - objCntIdx.x) * dx, (idy - objCntIdx.y) * dx);
		float2 dirxy;
		int LPs = VN.z / LAYERS;
		float dirZ;
		float minObj = 0;
		float s2vlength = 0;
		float3 cossinT;
		S2D = S2D / ddv;
		dbeta = 1.0 / dbeta;
		float invSID = rsqrtf(s.x * s.x + s.y * s.y);
		for (int lpIdx = 0; lpIdx != LPs; ++lpIdx)
		{
			minObj = (-objCntIdx.z + lpIdx * LAYERS) * dz;
			for (int angIdx = 0; angIdx < PN; ++angIdx)
			{
				cossinT = cossinZT[angIdx];
				cursour = make_float3(
					s.x * cossinT.x - s.y * cossinT.y,
					s.x * cossinT.y + s.y * cossinT.x,
					s.z + cossinT.z);
				dirxy.x = curvox_xy.x - cursour.x;
				dirxy.y = curvox_xy.y - cursour.y;
				s2vlength = hypotf(dirxy.x, dirxy.y);
				assert(s2vlength != 0);
				if (fabsf(cossinT.x) <= fabsf(cossinT.y))
				{
					dirlft = normalize(make_float2(dirxy.x, dirxy.y - 0.5 * dx));
					dirrgh = normalize(make_float2(dirxy.x, dirxy.y + 0.5 * dx));
					assert(dirxy.x != 0);
					cosVal = (dx * s2vlength / dirxy.x);
				}
				else
				{
					dirlft = normalize(make_float2(dirxy.x + 0.5f * dx, dirxy.y));
					dirrgh = normalize(make_float2(dirxy.x - 0.5f * dx, dirxy.y));
					assert(dirxy.y != 0);
					cosVal = (dx * s2vlength / dirxy.y);
				}
				idxL = asinf((cursour.y * dirlft.x - cursour.x * dirlft.y) * invSID) * dbeta + detCntIdx.x + 1;
				idxR = asinf((cursour.y * dirrgh.x - cursour.x * dirrgh.y) * invSID) * dbeta + detCntIdx.x + 1;

				curang = S2D / s2vlength;
#pragma unroll
				for (int idz = 0; idz <= LAYERS; ++idz)
				{
					dirZ = minObj + idz * dz - cursour.z;
					ddv = hypotf(dirZ, s2vlength) / s2vlength;
					idxD = (dirZ - 0.5 * dz) * curang + detCntIdx.y + 1;
					summ[threadIdx.y][threadIdx.x][idz] +=
						(tex3D<float>(prjTexObj, idxD, idxR, angIdx + 0.5) -
						tex3D<float>(prjTexObj, idxD, idxL, angIdx + 0.5)) * cosVal * ddv;
				}
			}
			__syncthreads();
			int vIdx = (idy * VN.x + idx) * VN.z + (lpIdx * LAYERS);
#pragma unroll
			for (int idz = 0; idz < LAYERS; ++idz)
			{
				vol[vIdx + idz] = summ[threadIdx.y][threadIdx.x][idz + 1] - summ[threadIdx.y][threadIdx.x][idz];
				summ[threadIdx.y][threadIdx.x][idz] = 0;
			}
			summ[threadIdx.y][threadIdx.x][LAYERS] = 0;
			__syncthreads();
		}
	}
}



template<>
__global__ void DD3_gpu_back_ker<_ZLINEBRANCHLESS>(
	cudaTextureObject_t prjTexObj,
	float* vol,
	const byte* __restrict__ msk,
	const float3* __restrict__ cossinZT,
	float3 s,
	float S2D,
	float3 objCntIdx,
	float dx, float dz, float dbeta, float ddv,
	float2 detCntIdx, int3 VN, int PN, int squared)
{
	int idx = threadIdx.x + __umul24(blockIdx.x, blockDim.x);
	int idy = threadIdx.y + __umul24(blockIdx.y, blockDim.y);
	__shared__ float summ[4][8][17];
#pragma unroll
	for (int i = 0; i <= 16; ++i)
	{
		summ[threadIdx.y][threadIdx.x][i] = 0;
	}
	__syncthreads();
	if (idx < VN.x && idy < VN.y)
	{
		if (msk[idy * VN.x + idx] != 1)
			return;
		float curang(0);
		float2 dirlft, dirrgh;
		float3 cursour;
		float idxL, idxR, idxD;
		float cosVal = 1.0;
		float2 curvox_xy = make_float2((idx - objCntIdx.x) * dx, (idy - objCntIdx.y) * dx);
		float2 dirxy;
		int LPs = VN.z >> 4;
		float dirZ;
		float minObj = 0;
		float s2vlength = 0;
		float3 cossinT;
		S2D = S2D / ddv;
		dbeta = 1.0 / dbeta;
		float invSID = rsqrtf(s.x * s.x + s.y * s.y);
		for (int lpIdx = 0; lpIdx != LPs; ++lpIdx)
		{
			minObj = (-objCntIdx.z + lpIdx * 16) * dz;
			for (int angIdx = 0; angIdx < PN; ++angIdx)
			{
				cossinT = cossinZT[angIdx];
				cursour = make_float3(
					s.x * cossinT.x - s.y * cossinT.y,
					s.x * cossinT.y + s.y * cossinT.x,
					s.z + cossinT.z);
				dirxy.x = curvox_xy.x - cursour.x;
				dirxy.y = curvox_xy.y - cursour.y;
				s2vlength = hypotf(dirxy.x, dirxy.y);
				assert(s2vlength != 0);
				if (fabsf(cossinT.x) <= fabsf(cossinT.y))
				{
					dirlft = normalize(make_float2(dirxy.x, dirxy.y - 0.5 * dx));
					dirrgh = normalize(make_float2(dirxy.x, dirxy.y + 0.5 * dx));
					assert(dirxy.x != 0);
					cosVal = (dx * s2vlength / dirxy.x);
				}
				else
				{
					dirlft = normalize(make_float2(dirxy.x + 0.5f * dx, dirxy.y));
					dirrgh = normalize(make_float2(dirxy.x - 0.5f * dx, dirxy.y));
					assert(dirxy.y != 0);
					cosVal = (dx * s2vlength / dirxy.y);
				}
				idxL = asinf((cursour.y * dirlft.x - cursour.x * dirlft.y) * invSID) * dbeta + detCntIdx.x + 1;
				idxR = asinf((cursour.y * dirrgh.x - cursour.x * dirrgh.y) * invSID) * dbeta + detCntIdx.x + 1;
				curang = S2D / s2vlength;
#pragma unroll
				for (int idz = 0; idz <= 16; ++idz)
				{
					dirZ = minObj + idz * dz - cursour.z;
					ddv = hypotf(dirZ, s2vlength) / s2vlength;
					idxD = (dirZ - 0.5 * dz) * curang + detCntIdx.y + 1;
					summ[threadIdx.y][threadIdx.x][idz] +=
						(tex3D<float>(prjTexObj, idxD, idxR, angIdx + 0.5) -
						tex3D<float>(prjTexObj, idxD, idxL, angIdx + 0.5)) * cosVal * ddv;
				}
			}
			__syncthreads();
			int vIdx = (idy * VN.x + idx) * VN.z + (lpIdx << 4);
#pragma unroll
			for (int idz = 0; idz < 16; ++idz)
			{
				vol[vIdx + idz] = summ[threadIdx.y][threadIdx.x][idz + 1] - summ[threadIdx.y][threadIdx.x][idz];
				summ[threadIdx.y][threadIdx.x][idz] = 0;
			}
			summ[threadIdx.y][threadIdx.x][16] = 0;
			__syncthreads();
		}
	}
}
static
void DD3_gpu_back_zlinebasedbranchless(
	float x0, float y0, float z0,
	int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN,
	float* hvol, float* hprj, float dx, float dz,
	byte* mask, int squared, int gpunum)
{
	cudaSetDevice(gpunum);
	cudaDeviceReset();
	const int TOTVON = XN * YN * ZN;
	const float objCntIdxX = (XN - 1.0) * 0.5 - imgXCenter / dx;
	const float objCntIdxY = (YN - 1.0) * 0.5 - imgYCenter / dx;
	const float objCntIdxZ = (ZN - 1.0) * 0.5 - imgZCenter / dz;

	thrust::device_vector<float> prjSAT(genSAT_of_Projection(hprj, DNU, DNV, PN));
	cudaExtent prjSize;
	prjSize.width = DNV + 1;
	prjSize.height = DNU + 1;
	prjSize.depth = PN;

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	cudaArray *d_prjSATarray;
	cudaMalloc3DArray(&d_prjSATarray, &channelDesc, prjSize);
	cudaMemcpy3DParms copyParams = { 0 };
	copyParams.srcPtr = make_cudaPitchedPtr(
		(void*) thrust::raw_pointer_cast(&prjSAT[0]), prjSize.width * sizeof(float),
		prjSize.width, prjSize.height);
	copyParams.dstArray = d_prjSATarray;
	copyParams.extent = prjSize;
	copyParams.kind = cudaMemcpyDeviceToDevice;
	cudaMemcpy3D(&copyParams);
	cudaResourceDesc resDesc;
	memset(&resDesc, 0, sizeof(resDesc));
	resDesc.resType = cudaResourceTypeArray;
	resDesc.res.array.array = d_prjSATarray;
	cudaTextureDesc texDesc;
	memset(&texDesc, 0, sizeof(texDesc));
	texDesc.addressMode[0] = cudaAddressModeClamp;
	texDesc.addressMode[1] = cudaAddressModeClamp;
	texDesc.addressMode[2] = cudaAddressModeClamp;
	texDesc.filterMode = cudaFilterModeLinear;
	texDesc.readMode = cudaReadModeElementType;
	cudaTextureObject_t texObj;
	texDesc.normalizedCoords = false;
	cudaCreateTextureObject(&texObj, &resDesc, &texDesc, nullptr);
	prjSAT.clear();

	float* bxds = new float[DNU + 1];
	float* byds = new float[DNU + 1];
	float* bzds = new float[DNV + 1];

	DD3Boundaries(DNU + 1, xds, bxds);
	DD3Boundaries(DNU + 1, yds, byds);
	DD3Boundaries(DNV + 1, zds, bzds);

	float ddv = (bzds[DNV] - bzds[0]) / DNV;
	float detCtrIdxV = (-(bzds[0] - z0) / ddv) - 0.5;

	float2 dir = normalize(make_float2(-x0, -y0));
	float2 dirL = normalize(make_float2(bxds[0] - x0, byds[0] - y0));
	float2 dirR = normalize(make_float2(bxds[DNU] - x0, byds[DNU] - y0));
	float dbeta = asin(dirL.x * dirR.y - dirL.y * dirR.x) / DNU;
	assert(dbeta != 0);
	float minBeta = asin(dir.x * dirL.y - dir.y * dirL.x);
	float detCtrIdxU = -minBeta / dbeta - 0.5;
	const float S2D = hypotf(xds[0] - x0, yds[0] - y0);

	thrust::device_vector<float> angs(hangs, hangs + PN);
	thrust::device_vector<float> zPos(hzPos, hzPos + PN);

	float3 sour(make_float3(x0, y0, z0));
	thrust::device_vector<float> vol(TOTVON, 0);
	thrust::device_vector<float3> cursour(PN);
	thrust::device_vector<float2> dirsour(PN);
	thrust::device_vector<float3> cossinZT(PN);


	thrust::transform(
		thrust::make_zip_iterator(thrust::make_tuple(angs.begin(), zPos.begin())),
		thrust::make_zip_iterator(thrust::make_tuple(angs.end(), zPos.end())),
		thrust::make_zip_iterator(thrust::make_tuple(cossinZT.begin(), cursour.begin(), dirsour.begin())),
		CTMBIR::ConstantForBackProjection3(x0, y0, z0));
	thrust::device_vector<byte> msk(mask, mask + XN * YN);
	dim3 blk_(8, 4);
	dim3 gid_(
		(XN + blk_.x - 1) / blk_.x,
		(YN + blk_.y - 1) / blk_.y);

	DD3_gpu_back_zlinebasedbranchless_ker<16> << <gid_, blk_ >> >(texObj,
		thrust::raw_pointer_cast(&vol[0]),
		thrust::raw_pointer_cast(&msk[0]),
		thrust::raw_pointer_cast(&cossinZT[0]),
		make_float3(x0, y0, z0),
		S2D,
		make_float3(objCntIdxX, objCntIdxY, objCntIdxZ),
		dx, dz, dbeta, ddv, make_float2(detCtrIdxU, detCtrIdxV),
		make_int3(XN, YN, ZN), PN, squared);

	thrust::copy(vol.begin(), vol.end(), hvol);


}


template<BackProjectionMethod METHOD>
static void DD3_gpu_back(
	float x0, float y0, float z0,
	int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN,
	float* hvol, float* hprj, float dx, float dz,
	byte* mask, int squared, int gpunum)
{
	cudaSetDevice(gpunum);
	cudaDeviceReset();
	//std::cout << gpunum << "\n";

	float3 objCntIdx = make_float3(
		(XN - 1.0) * 0.5 - imgXCenter / dx,
		(YN - 1.0) * 0.5 - imgYCenter / dx,
		(ZN - 1.0) * 0.5 - imgZCenter / dz);
	float3 sour = make_float3(x0, y0, z0);
	thrust::device_vector<byte> msk(mask, mask + XN * YN);
	thrust::device_vector<float> vol(XN * YN * ZN, 0);
	const float S2D = hypotf(xds[0] - x0, yds[0] - y0);

	thrust::device_vector<float3> cossinZT(PN);
	thrust::device_vector<float> angs(hangs, hangs + PN);
	thrust::device_vector<float> zPos(hzPos, hzPos + PN);
	thrust::transform(
		thrust::make_zip_iterator(thrust::make_tuple(angs.begin(), zPos.begin())),
		thrust::make_zip_iterator(thrust::make_tuple(angs.end(), zPos.end())),
		cossinZT.begin(),
		CTMBIR::ConstantForBackProjection4(x0, y0, z0));
	float4 detParas = calDetParas(xds, yds, zds, x0, y0, z0, DNU, DNV);
	cudaArray* d_prjArray = nullptr;
	cudaTextureObject_t texObj;
	dim3 blk;
	dim3 gid;
	thrust::device_vector<float> prjSAT;
	switch (METHOD)
	{
	case _PSEUDODD:
	case _VOLUMERENDERING:
		createTextureObject(texObj, d_prjArray, DNV, DNU, PN, hprj, cudaMemcpyHostToDevice,
			cudaAddressModeBorder, cudaFilterModeLinear, cudaReadModeElementType, false);
		blk.x = BACK_BLKX;
		blk.y = BACK_BLKY;
		blk.z = BACK_BLKZ;
		gid.x = (ZN + blk.x - 1) / blk.x;
		gid.y = (XN + blk.y - 1) / blk.y;
		gid.z = (YN + blk.z - 1) / blk.z;
		break;
	case _BRANCHLESS:
		prjSAT = genSAT_of_Projection(hprj, DNU, DNV, PN);
		createTextureObject(texObj, d_prjArray, DNV + 1, DNU + 1, PN,
			thrust::raw_pointer_cast(&prjSAT[0]),
			cudaMemcpyDeviceToDevice,
			cudaAddressModeClamp, cudaFilterModeLinear, cudaReadModeElementType, false);
		prjSAT.clear();
		blk.x = BACK_BLKX;
		blk.y = BACK_BLKY;
		blk.z = BACK_BLKZ;
		gid.x = (ZN + blk.x - 1) / blk.x;
		gid.y = (XN + blk.y - 1) / blk.y;
		gid.z = (YN + blk.z - 1) / blk.z;
		break;
	case _ZLINEBRANCHLESS:
		prjSAT = genSAT_of_Projection(hprj, DNU, DNV, PN);
		createTextureObject(texObj, d_prjArray, DNV + 1, DNU + 1, PN,
			thrust::raw_pointer_cast(&prjSAT[0]),
			cudaMemcpyDeviceToDevice,
			cudaAddressModeClamp, cudaFilterModeLinear, cudaReadModeElementType, false);
		prjSAT.clear();
		blk.x = 8;
		blk.y = 4;
		blk.z = 1;
		gid.x = (XN + blk.x - 1) / blk.x;
		gid.y = (YN + blk.y - 1) / blk.y;
		break;
	default:
		prjSAT = genSAT_of_Projection(hprj, DNU, DNV, PN);
		createTextureObject(texObj, d_prjArray, DNV + 1, DNU + 1, PN,
			thrust::raw_pointer_cast(&prjSAT[0]),
			cudaMemcpyDeviceToDevice,
			cudaAddressModeClamp, cudaFilterModeLinear, cudaReadModeElementType, false);
		prjSAT.clear();
		blk.x = BACK_BLKX;
		blk.y = BACK_BLKY;
		blk.z = BACK_BLKZ;
		gid.x = (ZN + blk.x - 1) / blk.x;
		gid.y = (XN + blk.y - 1) / blk.y;
		gid.z = (YN + blk.z - 1) / blk.z;
		break;
	}


	DD3_gpu_back_ker<METHOD> << <gid, blk >> >(texObj,
		thrust::raw_pointer_cast(&vol[0]),
		thrust::raw_pointer_cast(&msk[0]),
		thrust::raw_pointer_cast(&cossinZT[0]),
		make_float3(x0, y0, z0),
		S2D,
		make_float3(objCntIdx.x, objCntIdx.y, objCntIdx.z),
		dx, dz, detParas.z, detParas.w, make_float2(detParas.x, detParas.y),
		make_int3(XN, YN, ZN), PN, squared);
	copy(vol.begin(), vol.end(), hvol);

	destroyTextureObject(texObj, d_prjArray);

	vol.clear();
	msk.clear();
	angs.clear();
	zPos.clear();
	cossinZT.clear();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DOUBLE PRECISION BACKPROJECTION 

template<typename T1, typename T2>
static void DD3BoundariesDouble(int nrBoundaries,
	T1 *pCenters,
	T2 *pBoundaries)
{
	int i;
	if (nrBoundaries >= 3)
	{
		*pBoundaries++ = 1.5 * *pCenters - 0.5 * *(pCenters + 1);
		for (i = 1; i <= (nrBoundaries - 2); i++)
		{
			*pBoundaries++ = 0.5 * *pCenters + 0.5 * *(pCenters + 1);
			pCenters++;
		}
		*pBoundaries = 1.5 * *pCenters - 0.5 * *(pCenters - 1);
	}
	else
	{
		*pBoundaries = *pCenters - 0.5;
		*(pBoundaries + 1) = *pCenters + 0.5;
	}
}


//Calculate the detCtrIdx and detector cell size parameters.
static
double4 calDetParasDouble(float* xds, float* yds, float* zds, float x0, float y0, float z0, int DNU, int DNV)
{
	double* bxds = new double[DNU + 1];
	double* byds = new double[DNU + 1];
	double* bzds = new double[DNV + 1];

	DD3BoundariesDouble<float, double>(DNU + 1, xds, bxds);
	DD3BoundariesDouble<float, double>(DNU + 1, yds, byds);
	DD3BoundariesDouble<float, double>(DNV + 1, zds, bzds);

	// detector size in Z
	double ddv = (bzds[DNV] - bzds[0]) / DNV; // detector size in Z direction
	double detCtrIdxV = (-(bzds[0] - z0) / ddv) - 0.5; // detector center index in Z direction
	double2 dir = normalize(make_double2(-x0, -y0)); // source to origin vector (XY plane)
	double2 dirL = normalize(make_double2(bxds[0] - x0, byds[0] - y0)); // Left boundary direction vector
	double2 dirR = normalize(make_double2(bxds[DNU] - x0, byds[DNU] - y0)); // Right boundary direction vector
																			// the angular size of detector cell as seen by the source
	double dbeta = asin(dirL.x * dirR.y - dirL.y * dirR.x) / DNU; //detector size in channel direction
	double minBeta = asin(dir.x * dirL.y - dir.y * dirL.x); //the fan angle corresponding to the most left boundary
	double detCtrIdxU = -minBeta / dbeta - 0.5; //det center index in XY / channel direction
	delete[] bxds;
	delete[] byds;
	delete[] bzds;
	return make_double4(detCtrIdxU, detCtrIdxV, dbeta, ddv);
}



__global__ void horizontalIntegral2(double* prj, int DNU, int DNV, int PN)
{
	int idv = threadIdx.x + blockIdx.x * blockDim.x;
	int pIdx = threadIdx.y + blockIdx.y * blockDim.y;
	if (idv < DNV && pIdx < PN)
	{
		int headPtr = pIdx * DNU * DNV + idv;
		for (int ii = 1; ii < DNU; ++ii)
		{
			prj[headPtr + ii * DNV] = prj[headPtr + ii * DNV] + prj[headPtr + (ii - 1) * DNV];
		}
	}
}

struct TransFromDoubleToInt2 {
	__device__ int2 operator()(double a) const {
		return make_int2(__double2loint(a), __double2hiint(a));
	}
};
static
thrust::device_vector<int2> genSAT_of_Projection_Double(
	const float* const hprj, // Projection data in host memory
	int DNU, // Detector channel #
	int DNV, // Detector row #
	int PN)  // view #
{
	const int siz = DNU * DNV * PN;
	const int nsiz = (DNU + 1) * (DNV + 1) * PN;
	thrust::device_vector<double> prjSAT(nsiz, 0); // SAT in device memory
	thrust::device_vector<double> prj(hprj, hprj + siz); // projection in device memory
	dim3 copyBlk(64, 16, 1);
	dim3 copyGid(
		(DNV + copyBlk.x - 1) / copyBlk.x,
		(DNU + copyBlk.y - 1) / copyBlk.y,
		(PN + copyBlk.z - 1) / copyBlk.z);
	addOneSidedZeroBoarder << <copyGid, copyBlk >> >(
		thrust::raw_pointer_cast(&prj[0]),
		thrust::raw_pointer_cast(&prjSAT[0]),
		DNU, DNV, PN);
	const int nDNU = DNU + 1;
	const int nDNV = DNV + 1;

	//Generate SAT inplace
	copyBlk.x = 512;
	copyBlk.y = 1;
	copyBlk.z = 1;
	copyGid.x = (nDNU * PN + copyBlk.x - 1) / copyBlk.x;
	copyGid.y = 1;
	copyGid.z = 1;
	verticalIntegral2 << <copyGid, copyBlk >> >(thrust::raw_pointer_cast(&prjSAT[0]), nDNV, nDNU * PN);

	copyBlk.x = 64;
	copyBlk.y = 16;
	copyBlk.z = 1;
	copyGid.x = (nDNV + copyBlk.x - 1) / copyBlk.x;
	copyGid.y = (PN + copyBlk.y - 1) / copyBlk.y;
	copyGid.z = 1;
	horizontalIntegral2 << <copyGid, copyBlk >> >(thrust::raw_pointer_cast(&prjSAT[0]), nDNU, nDNV, PN);
	prj.clear();
	thrust::device_vector<int2> prjSATInt(nsiz);

	thrust::transform(prjSAT.begin(), prjSAT.end(), prjSATInt.begin(), TransFromDoubleToInt2());
	prjSAT.clear();
	return prjSATInt; // It has the device to device memory copy (TODO: avoid this by using reference parameters)
}


__global__ void DD3_gpu_back_ker_double(
	cudaTextureObject_t texObj,
	double* vol,
	byte* __restrict__ msk,
	double3* __restrict__ cossinZT,
	double3 s,
	double S2D,
	double3 curvox,
	double dx, double dz, double dbeta, double ddv,
	double2 detCntIdx, int3 VN, int PN, int squared)
{
	int3 id;
	id.z = threadIdx.x + __umul24(blockIdx.x, blockDim.x);  // BACK_BLKX
	id.x = threadIdx.y + __umul24(blockIdx.y, blockDim.y);  // BACK_BLKY
	id.y = threadIdx.z + __umul24(blockIdx.z, blockDim.z);  // BACK_BLKZ

	if (id.x < VN.x && id.y < VN.y && id.z < VN.z)
	{
		if (msk[id.y * VN.x + id.x] != 1)
			return;
		// Position of current pixel
		curvox = (id - curvox) * make_double3(dx, dx, dz);// make_float3((id.x - objCntIdx.x) * dx, (id.y - objCntIdx.y) * dx, (id.z - objCntIdx.z) * dz);
		double3 cursour; // src position (precomputed in global mem "cursours"
		double idxL, idxR, idxU, idxD; // detctor index corresponding to shadow of the current pixel
		double cosVal; // ray angle relative to the normal vector of detector face
		double summ = 0;


		double3 cossin;
		double inv_sid = 1.0 / sqrt(s.x * s.x + s.y * s.y);
		double3 dir;
		double l_square;
		double l;
		double alpha;

		double deltaAlpha;
		S2D = S2D / ddv;
		dbeta = 1.0 / dbeta;
		dz = dz * 0.5;

		int cidxU, fidxU, cidxD, fidxD, cidxL, fidxL, cidxR, fidxR;
		for (int angIdx = 0; angIdx < PN; ++angIdx)
		{
			cossin = cossinZT[angIdx];
			cursour = make_double3(
				s.x * cossin.x - s.y * cossin.y,
				s.x * cossin.y + s.y * cossin.x,
				s.z + cossin.z);

			dir = curvox - cursour;
			l_square = dir.x * dir.x + dir.y * dir.y;
			assert(l_square != 0);
			l = rsqrt(l_square);
			idxU = (dir.z + dz) * S2D * l + detCntIdx.y + 1; //0.5 offset Because we use the texture fetching
			idxD = (dir.z - dz) * S2D * l + detCntIdx.y + 1;
			alpha = asin((cursour.y * dir.x - cursour.x * dir.y) * inv_sid * l);

			if (fabs(cursour.x) > fabs(cursour.y))
			{
				ddv = dir.x;
			}
			else
			{
				ddv = dir.y;
			}
			deltaAlpha = ddv / l_square * dx * 0.5;
			cosVal = dx / ddv * sqrt(l_square + dir.z * dir.z);
			idxL = (alpha - deltaAlpha) * dbeta + detCntIdx.x + 1;
			idxR = (alpha + deltaAlpha) * dbeta + detCntIdx.x + 1;

			cidxU = ceil(idxU);
			fidxU = floor(idxU);
			cidxD = ceil(idxD);
			fidxD = floor(idxD);
			cidxL = ceil(idxL);
			fidxL = floor(idxL);
			cidxR = ceil(idxR);
			fidxR = floor(idxR);

			summ += (
				-(bilerp(tex3D<int2>(texObj, fidxD, fidxR, angIdx + 0.5),
					tex3D<int2>(texObj, fidxD, cidxR, angIdx + 0.5),
					tex3D<int2>(texObj, cidxD, fidxR, angIdx + 0.5),
					tex3D<int2>(texObj, cidxD, cidxR, angIdx + 0.5),
					idxR - fidxR, idxD - fidxD))
				- (bilerp(tex3D<int2>(texObj, fidxU, fidxL, angIdx + 0.5),
					tex3D<int2>(texObj, fidxU, cidxL, angIdx + 0.5),
					tex3D<int2>(texObj, cidxU, fidxL, angIdx + 0.5),
					tex3D<int2>(texObj, cidxU, cidxL, angIdx + 0.5),
					idxL - fidxL, idxU - fidxU))
				+ (bilerp(tex3D<int2>(texObj, fidxD, fidxL, angIdx + 0.5),
					tex3D<int2>(texObj, fidxD, cidxL, angIdx + 0.5),
					tex3D<int2>(texObj, cidxD, fidxL, angIdx + 0.5),
					tex3D<int2>(texObj, cidxD, cidxL, angIdx + 0.5),
					idxL - fidxL, idxD - fidxD))
				+ (bilerp(tex3D<int2>(texObj, fidxU, fidxR, angIdx + 0.5),
					tex3D<int2>(texObj, fidxU, cidxR, angIdx + 0.5),
					tex3D<int2>(texObj, cidxU, fidxR, angIdx + 0.5),
					tex3D<int2>(texObj, cidxU, cidxR, angIdx + 0.5),
					idxR - fidxR, idxU - fidxU))) * cosVal;
		}
		__syncthreads();
		vol[__umul24((__umul24(id.y, VN.x) + id.x), VN.z) + id.z] = summ;// * MSK[threadIdx.y][threadIdx.x];
	}
}

//Create texture object and corresponding cudaArray function
template<typename T>
static void createTextureObjectDouble(
	cudaTextureObject_t& texObj, //return: texture object pointing to the cudaArray
	cudaArray* d_prjArray, // return: cudaArray storing the data
	int Width, int Height, int Depth, // data size
	T* sourceData, // where is the data
	cudaMemcpyKind memcpyKind, // data from host or memory
	cudaTextureAddressMode addressMode, // how to address the texture (clamp, border ...)
	cudaTextureFilterMode textureFilterMode, // usually linear filtering (double --> int2 use pointer not linear interpolation)
	cudaTextureReadMode textureReadMode, // usually use element wise reading mode.
	bool isNormalized) // usually false
{
	cudaExtent prjSize;
	prjSize.width = Width;
	prjSize.height = Height;
	prjSize.depth = Depth;
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<T>();
	cudaMalloc3DArray(&d_prjArray, &channelDesc, prjSize);
	cudaMemcpy3DParms copyParams = { 0 };
	copyParams.srcPtr = make_cudaPitchedPtr(
		(void*)sourceData, prjSize.width * sizeof(T),
		prjSize.width, prjSize.height);
	copyParams.dstArray = d_prjArray;
	copyParams.extent = prjSize;
	copyParams.kind = memcpyKind;
	CUDA_CHECK_RETURN(cudaMemcpy3D(&copyParams));
	cudaResourceDesc resDesc;
	memset(&resDesc, 0, sizeof(resDesc));
	resDesc.resType = cudaResourceTypeArray;
	resDesc.res.array.array = d_prjArray;
	cudaTextureDesc texDesc;
	memset(&texDesc, 0, sizeof(texDesc));
	texDesc.addressMode[0] = addressMode;
	texDesc.addressMode[1] = addressMode;
	texDesc.addressMode[2] = addressMode;
	texDesc.filterMode = textureFilterMode;
	texDesc.readMode = textureReadMode;

	texDesc.normalizedCoords = isNormalized;
	CUDA_SAFE_CALL(cudaCreateTextureObject(&texObj, &resDesc, &texDesc, nullptr));
}

struct ConstantForBackProjectionDouble
{
	typedef thrust::tuple<double, double> InTuple;

	ConstantForBackProjectionDouble(const double _x0, const double _y0, const double _z0) :x0(_x0), y0(_y0), z0(_z0) {}

	double x0;
	double y0;
	double z0;

	__device__ double3 operator()(const InTuple& tp)
	{
		double curang = regularizeAngle(thrust::get<0>(tp));
		double zP = thrust::get<1>(tp);
		double cosT = cos(curang);
		double sinT = sin(curang);

		return make_double3(cosT, sinT, zP);
	}
};

static
void DD3_gpu_back_double(
	const float x0, const float y0, const float z0,
	const int DNU, const int DNV,
	const float *  xds, const float *  yds, const float *  zds,
	const float imgXCenter, const float imgYCenter, const float imgZCenter,
	const float *  hangs, const float *  hzPos, const int PN,
	const int XN, const int YN, const int ZN,
	float *hvol, const float *  hprj,
	const float dx, const float dz,
	const byte * const mask, const int squared, const int gpunum)
{
	CUDA_CHECK_RETURN(cudaSetDevice(gpunum));
	CUDA_CHECK_RETURN(cudaDeviceReset());
	double3 objCntIdx = make_double3((XN - 1.0) * 0.5 - imgXCenter / dx, (YN - 1.0) * 0.5 - imgYCenter / dx, (ZN - 1.0) * 0.5 - imgZCenter / dz);
	double3 sour = make_double3(x0, y0, z0);
	thrust::device_vector<byte> msk(mask, mask + XN * YN);
	thrust::device_vector<double> vol(XN * YN * ZN, 0);
	const double S2D = hypot(xds[0] - x0, yds[0] - y0);
	thrust::device_vector<double3> cossinZT(PN);
	thrust::device_vector<double> angs(hangs, hangs + PN);
	thrust::device_vector<double> zPos(hzPos, hzPos + PN);
	thrust::transform(
		thrust::make_zip_iterator(thrust::make_tuple(angs.begin(), zPos.begin())),
		thrust::make_zip_iterator(thrust::make_tuple(angs.end(), zPos.end())),
		cossinZT.begin(),
		ConstantForBackProjectionDouble(x0, y0, z0));
	//Calculate detCtrIdxU, detCtrIdxV, dbeta, ddv and bind them in a float4 datatype
	double4 detParas = calDetParasDouble(
		const_cast<float*>(xds),
		const_cast<float*>(yds),
		const_cast<float*>(zds), x0, y0, z0, DNU, DNV);

	cudaArray *d_prjArray = nullptr;
	cudaTextureObject_t texObj;
	dim3 blk;
	dim3 gid;
	thrust::device_vector<int2> prjSAT;

	prjSAT = genSAT_of_Projection_Double(const_cast<float*>(hprj), DNU, DNV, PN);
	createTextureObjectDouble<int2>(texObj, d_prjArray, DNV + 1, DNU + 1, PN,
		thrust::raw_pointer_cast(&prjSAT[0]), cudaMemcpyDeviceToDevice,
		cudaAddressModeClamp, cudaFilterModePoint, cudaReadModeElementType, false);
	prjSAT.clear();

	blk.x = BACK_BLKX;
	blk.y = BACK_BLKY;
	blk.z = BACK_BLKZ;
	gid.x = (ZN + blk.x - 1) / blk.x;
	gid.y = (XN + blk.y - 1) / blk.y;
	gid.z = (YN + blk.z - 1) / blk.z;

	DD3_gpu_back_ker_double << <gid, blk >> >(texObj,
		thrust::raw_pointer_cast(&vol[0]),
		thrust::raw_pointer_cast(&msk[0]),
		thrust::raw_pointer_cast(&cossinZT[0]),
		sour,
		S2D,
		objCntIdx,
		dx, dz, detParas.z, detParas.w,
		make_double2(detParas.x, detParas.y),
		make_int3(XN, YN, ZN), PN, static_cast<int>(squared != 0));

	thrust::copy(vol.begin(), vol.end(), hvol);
	destroyTextureObject(texObj, d_prjArray);

	vol.clear();
	msk.clear();
	angs.clear();
	zPos.clear();
	cossinZT.clear();
}


///////////////////////////////////////////////////////////////////////////////////////////
// Branches backprojection 
///////////////////////////////////////////////////////////////////////////////////////////
__device__ float crossProduct(const float2 a, const float2 b)
{
	return a.x * b.y - a.y * b.x;
}

__global__ void DD3_back_branches_ker(
	const float* proj, float* vol,
	float S2O, float detCntIdU, float detCntIdV,
	int XN, int YN, int ZN, int DNU, int DNV, int PN,
	float dbeta, float ddv, float dx, float dz,
	float S2D,
	float3* cossinZT,
	float x0, float y0, float z0,
	float objCntIdxX, float objCntIdxY, float objCntIdxZ, byte* mask)
{
	const int kk = threadIdx.x + blockIdx.x * blockDim.x;
	const int ii = threadIdx.y + blockIdx.y * blockDim.y;
	const int jj = threadIdx.z + blockIdx.z * blockDim.z;
	if (ii < XN && jj < YN && kk < ZN)
	{
		if (mask[jj * XN + ii] != 1)
			return;
		float summ = 0;

		float3 curvox = make_float3((ii - objCntIdxX) * dx, (jj - objCntIdxY) * dx, (kk - objCntIdxZ) * dz); // current voxel position
		for (int angIdx = 0; angIdx < PN; ++angIdx)
		{
			float3 cossinzt = cossinZT[angIdx];
			float3 cursour = make_float3(
				x0 * cossinzt.x - y0 * cossinzt.y,
				x0 * cossinzt.y + y0 * cossinzt.x,
				z0 + cossinzt.z);
			float3 dir = (cursour - curvox);
			float ldir = length(dir);
			//dir = normalize(dir);

			float2 dirSour = normalize(make_float2(-cursour.x, -cursour.y));

			float2 dir1;
			float2 dir2;
			float betaLft;
			float betaRgh;
			float w;

			if (fabsf(dir.x) < fabsf(dir.y))
			{
				dir1 = normalize(make_float2(curvox.x - 0.5 * dx - cursour.x, curvox.y - cursour.y));
				dir2 = normalize(make_float2(curvox.x + 0.5 * dx - cursour.x, curvox.y - cursour.y));
				w = fabsf(dir.y);
			}
			else
			{
				dir1 = normalize(make_float2(curvox.x - cursour.x, curvox.y - cursour.y - 0.5 * dx));
				dir2 = normalize(make_float2(curvox.x - cursour.x, curvox.y - cursour.y + 0.5 * dx));
				w = fabsf(dir.x);
			}
			assert(w != 0);

			betaLft = asinf(crossProduct(dirSour, dir1));
			betaRgh = asinf(crossProduct(dirSour, dir2));

			if (betaLft > betaRgh)
			{
				float tmp = betaLft;
				betaLft = betaRgh;
				betaRgh = tmp;
			}
			float betaWidth = betaRgh - betaLft;
			assert(betaWidth != 0);
			int idxL = floorf(betaLft / dbeta + detCntIdU) - 2;
			int idxR = ceilf(betaRgh / dbeta + detCntIdU) + 2;
			if (idxL > DNU - 1) { continue; }
			if (idxR < 0) { continue; }
			if (idxL < 0) { idxL = 0; }
			if (idxR > DNU - 1) { idxR = DNU - 1; }

			float sour2voxHorizontalDis = sqrtf((cursour.x - curvox.x) * (cursour.x - curvox.x) + (cursour.y - curvox.y) * (cursour.y - curvox.y));
			assert(sour2voxHorizontalDis != 0);
			float relvHeightUp = S2D / sour2voxHorizontalDis * (curvox.z + 0.5f * dz - cursour.z);
			float relvHeightDn = S2D / sour2voxHorizontalDis * (curvox.z - 0.5f * dz - cursour.z);
			
			int idxU = ceilf(relvHeightUp / ddv + detCntIdV) + 2;
			int idxD = floorf(relvHeightDn / ddv + detCntIdV) - 2;
			if (idxU < 0) { continue; }
			if (idxD > DNV - 1) { continue; }
			if (idxU > DNV - 1) { idxU = DNV - 1; }
			if (idxD < 0) { idxD = 0; }

			for (int uIdx = idxL; uIdx <= idxR; ++uIdx)
			{
				float lftBeta = (uIdx - detCntIdU - 0.5f) * dbeta;
				float rghBeta = (uIdx - detCntIdU + 0.5f) * dbeta;
				float betaIntersection = intersectLength<float>(betaLft, betaRgh, lftBeta, rghBeta);
				if (betaIntersection > 0) {
					for (int vIdx = idxD; vIdx <= idxU; ++vIdx)
					{
						float dnZ = (vIdx - detCntIdV - 0.5f) * ddv;
						float upZ = (vIdx - detCntIdV + 0.5f) * ddv;
						float zIntersection = intersectLength<float>(relvHeightDn, relvHeightUp, dnZ, upZ);
						if (zIntersection > 0)
						{
							summ += proj[(angIdx * DNU + uIdx) * DNV + vIdx] * dx / w * ldir * (betaIntersection * zIntersection) / betaWidth;// / (betaWidth * relvHeight);
						}
						else
						{
							summ += 0;
						}
					}

				}

			}
		} // End for angIdx
		vol[(jj * XN + ii) * ZN + kk] = summ;
	}
}

static
void DD3Back_branch_gpu(float x0, float y0, float z0, int DNU, int DNV,
	float* xds, float* yds, float* zds, float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN, int XN, int YN, int ZN, float* hvol, float* hprj,
	float dx, float dz, byte* mask, int gpunum, int squared)
{
	CUDA_CHECK_RETURN(cudaSetDevice(gpunum));
	CUDA_CHECK_RETURN(cudaDeviceReset());
	thrust::device_vector<float> proj(hprj, hprj + DNU * DNV * PN);
	thrust::device_vector<float> vol(hvol, hvol + XN * YN * ZN);
	thrust::device_vector<byte> dmask(mask, mask + XN * YN);

	thrust::host_vector<float3> cossinZT(PN);
	for (int i = 0; i != PN; ++i)
	{
		cossinZT[i].x = cosf(static_cast<float>(i) / static_cast<float>(PN) * 3.141592653589793 * 2.0f);
		cossinZT[i].y = sinf(static_cast<float>(i) / static_cast<float>(PN) * 3.141592653589793 * 2.0f);
		cossinZT[i].z = hzPos[i];
	}
	thrust::device_vector<float3> dcossinZT = cossinZT;

	const float objCntIdxX = (XN - 1.0f) * 0.5f;
	const float objCntIdxY = (YN - 1.0f) * 0.5f;
	const float objCntIdxZ = (ZN - 1.0f) * 0.5f;

	const float S2O = y0;
	const float S2D = sqrtf(powf(x0 - xds[0], 2.0f) + powf(y0 - yds[0], 2.0f));

	thrust::host_vector<float> bxds(DNU + 1, 0);
	thrust::host_vector<float> byds(DNU + 1, 0);
	thrust::host_vector<float> bzds(DNV + 1, 0);

	DD3Boundaries(DNU + 1, xds, &bxds[0]);
	DD3Boundaries(DNU + 1, yds, &byds[0]);
	DD3Boundaries(DNV + 1, zds, &bzds[0]);

	const float detWidth = sqrtf(powf(xds[1] - xds[0], 2.0f) + powf(yds[0] - y0, 2.0f));
	const float dbeta = atanf(detWidth *0.5f / S2D) * 2.0f;
	assert(byds[0] != y0);
	const float lftTheta = atanf(bxds[0] / (byds[0] - y0));
	const float detCntIdU = lftTheta / dbeta - 0.5f;
	const float ddv = zds[1] - zds[0];
	const float detCntIdV = -bzds[0] / ddv - 0.5f;

	dim3 blk(64, 4, 1);
	dim3 gid(
		(ZN + blk.x - 1) / blk.x,
		(XN + blk.y - 1) / blk.y,
		(YN + blk.z - 1) / blk.z);

	DD3_back_branches_ker << <gid, blk >> >(
		thrust::raw_pointer_cast(&proj[0]),
		thrust::raw_pointer_cast(&vol[0]),
		S2O,
		detCntIdU, detCntIdV,
		XN, YN, ZN, DNU, DNV, PN,
		dbeta, ddv, dx, dz, S2D,
		thrust::raw_pointer_cast(&dcossinZT[0]),
		x0, y0, z0, objCntIdxX, objCntIdxY, objCntIdxZ, thrust::raw_pointer_cast(&dmask[0]));
	thrust::copy(vol.begin(), vol.end(), hvol);
}




extern "C"
void DD3Back_gpu(
float x0, float y0, float z0,
int DNU, int DNV,
float* xds, float* yds, float* zds,
float imgXCenter, float imgYCenter, float imgZCenter,
float* hangs, float* hzPos, int PN,
int XN, int YN, int ZN,
float* hvol, float* hprj,
float dx, float dz,
byte* mask, int gpunum, int squared, int prjMode)
{
	switch (prjMode)
	{
	case 0:
		DD3_gpu_back<_BRANCHLESS>(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
			hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz, mask, squared, gpunum);
		break;
	case 1:
		DD3_gpu_back<_PSEUDODD>(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
			hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz, mask, squared, gpunum);
		break;
	case 2:
		DD3_gpu_back_double(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
			hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz, mask, squared, 0);
		break;
	case 3:
		DD3_gpu_back<_ZLINEBRANCHLESS>(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
			hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz, mask, squared, gpunum);
		break;
	case 4:
		DD3Back_branch_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
			hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz, mask, 0, 0);
		break;
	default:
		DD3_gpu_back<_BRANCHLESS>(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
			hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz, mask, squared, gpunum);
		break;
	}
}


extern "C"
void DD3BackHelical_3GPU(
	float x0, float y0, float z0,
	int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN,
	float* hvol, float* hprj,
	float dx, float dz,
	byte* mask, int methodId, int (&startVOL)[3])
{
	thrust::host_vector<float> h_angs(hangs,hangs+PN);
	thrust::host_vector<float> h_zPos(hzPos,hzPos+PN);

	int ObjZIdx_Start[3] = {startVOL[0], startVOL[1], startVOL[2]};
	int ObjZIdx_End[3] = {startVOL[1], startVOL[2], ZN};

	int prjIdx_Start[3] = {0, 0, 0};
	int prjIdx_End[3] = {0, 0, 0};

	float objCntIdxZ = (ZN - 1.0) / 2.0 - imgZCenter / dz; //object center in Z direction
	float detStpZ = zds[1] - zds[0]; // detector cell height
	assert(detStpZ != 0);
	float detCntIdxV = -zds[0] / detStpZ; // Detector center along Z direction

	int SZN[3] = {startVOL[1] - startVOL[0], startVOL[2] - startVOL[1], ZN - startVOL[2]};

	float** subVol = new float*[3];
	subVol[0] = new float[XN * YN * SZN[0]];
	subVol[1] = new float[XN * YN * SZN[1]];
	subVol[2] = new float[XN * YN * SZN[2]];

	float subImgZCenter[3];//
	int SPN[3];

	omp_set_num_threads(3);
#pragma omp parallel for
	for(int i = 0; i < 3; ++i)
	{
		getSubVolume<float>(hvol, XN * YN, ZN,
				ObjZIdx_Start[i], ObjZIdx_End[i],
				subVol[i]);
		getPrjIdxPair<float>(h_zPos, ObjZIdx_Start[i], ObjZIdx_End[i],
						objCntIdxZ, dz, ZN, detCntIdxV, detStpZ, DNV, prjIdx_Start[i], prjIdx_End[i]);
		SPN[i] = prjIdx_End[i] - prjIdx_Start[i];
		std::cout<<i<<" "<<prjIdx_Start[i]<<" "<<prjIdx_End[i]<<"\n";

		subImgZCenter[i] = ((ObjZIdx_End[i] + ObjZIdx_Start[i] - (ZN - 1.0)) * dz + imgZCenter * 2.0) / 2.0;

	}
	int prefixSPN[3] = {prjIdx_Start[0], prjIdx_Start[1], prjIdx_Start[2]};

//	//Not implemented yet.o
#pragma omp parallel for
	for(int i = 0; i < 3; ++i)
	{
		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
				imgXCenter, imgYCenter, subImgZCenter[i],
				hangs + prefixSPN[i] , hzPos + prefixSPN[i],
				SPN[i],	XN, YN, SZN[i], subVol[i],
				hprj + DNU * DNV * prefixSPN[i],dx,dz,mask,i,0,methodId);
	}

	//Gather the volumes together
	combineVolume<float>(hvol,XN,YN,ZN,subVol,SZN,3);
	delete[] subVol[0];
	delete[] subVol[1];
	delete[] subVol[2];
	delete[] subVol;


}


extern "C"
void DD3BackHelical_4GPU(
	float x0, float y0, float z0,
	int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN,
	float* hvol, float* hprj,
	float dx, float dz,
	byte* mask, int methodId, int (&startVOL)[4])
{
	thrust::host_vector<float> h_angs(hangs,hangs+PN);
	thrust::host_vector<float> h_zPos(hzPos,hzPos+PN);

	int ObjZIdx_Start[4] = {startVOL[0], startVOL[1], startVOL[2], startVOL[3]};
	int ObjZIdx_End[4] = {startVOL[1], startVOL[2], startVOL[3], ZN};

	int prjIdx_Start[4] = {0, 0, 0, 0};
	int prjIdx_End[4] = {0, 0, 0, 0};

	float objCntIdxZ = (ZN - 1.0) / 2.0 - imgZCenter / dz; //object center in Z direction
	float detStpZ = zds[1] - zds[0]; // detector cell height
	assert(detStpZ != 0);
	float detCntIdxV = -zds[0] / detStpZ; // Detector center along Z direction

	int SZN[4] = {startVOL[1] - startVOL[0], startVOL[2] - startVOL[1], startVOL[3] - startVOL[2], ZN - startVOL[3]};

	float** subVol = new float*[4];
	subVol[0] = new float[XN * YN * SZN[0]];
	subVol[1] = new float[XN * YN * SZN[1]];
	subVol[2] = new float[XN * YN * SZN[2]];
	subVol[3] = new float[XN * YN * SZN[3]];

	float subImgZCenter[4];//
	int SPN[4];

	omp_set_num_threads(4);
#pragma omp parallel for
	for(int i = 0; i < 4; ++i)
	{
		getSubVolume<float>(hvol, XN * YN, ZN,
				ObjZIdx_Start[i], ObjZIdx_End[i],
				subVol[i]);
		getPrjIdxPair<float>(h_zPos, ObjZIdx_Start[i], ObjZIdx_End[i],
						objCntIdxZ, dz, ZN, detCntIdxV, detStpZ, DNV, prjIdx_Start[i], prjIdx_End[i]);
		SPN[i] = prjIdx_End[i] - prjIdx_Start[i];
		std::cout<<i<<" "<<prjIdx_Start[i]<<" "<<prjIdx_End[i]<<"\n";

		subImgZCenter[i] = ((ObjZIdx_End[i] + ObjZIdx_Start[i] - (ZN - 1.0)) * dz + imgZCenter * 2.0) / 2.0;

	}
	int prefixSPN[4] = {prjIdx_Start[0], prjIdx_Start[1], prjIdx_Start[2], prjIdx_Start[3]};

//	//Not implemented yet.o
#pragma omp parallel for
	for(int i = 0; i < 4; ++i)
	{
		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
				imgXCenter, imgYCenter, subImgZCenter[i],
				hangs + prefixSPN[i] , hzPos + prefixSPN[i],
				SPN[i],	XN, YN, SZN[i], subVol[i],
				hprj + DNU * DNV * prefixSPN[i],dx,dz,mask,i,0,methodId);
	}

	//Gather the volumes together
	combineVolume<float>(hvol,XN,YN,ZN,subVol,SZN,4);
	delete[] subVol[0];
	delete[] subVol[1];
	delete[] subVol[2];
	delete[] subVol[3];
	delete[] subVol;


}


void DD3Back_branchless_sat2d_multiGPU(
	float x0, float y0, float z0,
	int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* h_angs, float* h_zPos, int PN,
	int XN, int YN, int ZN,
	float* hvol, float* hprj,
	float dx, float dz,
	byte* mask, int* startVOL, int gpuNum)
{
	const int nDNU = DNU + 1;
	const int nDNV = DNV + 1;

	thrust::host_vector<float> hangs(h_angs, h_angs + PN);
	thrust::host_vector<float> hzPos(h_zPos, h_zPos + PN);

	std::vector<int> ObjZIdx_Start(startVOL, startVOL + gpuNum);
	std::vector<int> ObjZIdx_End(ObjZIdx_Start.size());
	std::copy(ObjZIdx_Start.begin() + 1, ObjZIdx_Start.end(), ObjZIdx_End.begin());
	ObjZIdx_End[gpuNum - 1] = ZN;

	std::vector<int> prjIdx_Start(gpuNum);
	std::vector<int> prjIdx_End(gpuNum);

	const float objCntIdxZ = (ZN - 1.0f) * 0.5 - imgZCenter / dz;
	const float detStpZ = (zds[DNV - 1] - zds[0]) / (DNV - 1.0f); // detector cell height
	assert(detStpZ != 0);
	const float detCntIdxV = -zds[0] / detStpZ; // Detector Center along Z direction

	std::vector<int> SZN = ObjZIdx_End - ObjZIdx_Start; // sub volume slices number

	std::vector<float> subImgZCenter(gpuNum,0.0f);
	std::vector<int> SPN(gpuNum);

	const float objCntIdxX = (XN - 1.0f) * 0.5f - imgXCenter / dx;
	const float objCntIdxY = (YN - 1.0f) * 0.5f - imgYCenter / dx;

	std::vector<float3> sour(gpuNum);
	thrust::host_vector<thrust::device_vector<byte> > msk(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > vol(gpuNum);
	thrust::host_vector<thrust::device_vector<float3> > cossinZT(gpuNum);
	thrust::host_vector<cudaArray*> d_prjArray(gpuNum);
	thrust::host_vector<cudaTextureObject_t> texObj(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > prjSAT(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > prj(gpuNum);
	thrust::host_vector<cudaStream_t> stream(gpuNum);

	const float4 detParas = calDetParas(xds, yds, zds, x0, y0, z0, DNU, DNV);
	const float S2D = hypotf(xds[0] - x0, yds[0] - y0);

	// Pre calculate the cossin z positions
	thrust::device_vector<float3> COSSINZT(PN);
	thrust::device_vector<float> ANGS = hangs;
	thrust::device_vector<float> ZPOS = hzPos;
	thrust::transform(
		thrust::make_zip_iterator(thrust::make_tuple(ANGS.begin(), ZPOS.begin())),
		thrust::make_zip_iterator(thrust::make_tuple(ANGS.end(), ZPOS.end())),
		COSSINZT.begin(), CTMBIR::ConstantForBackProjection4(x0, y0, z0));

	dim3 copyBlk(64,16,1);
	thrust::host_vector<dim3> copyGid(gpuNum);
	dim3 blk(BACK_BLKX, BACK_BLKY, BACK_BLKZ);
	thrust::host_vector<dim3> gid(gpuNum);
	dim3 vertGenBlk(512,1,1);
	thrust::host_vector<dim3> vertGenGid(gpuNum);
	dim3 horzGenBlk(64,16,1);
	thrust::host_vector<dim3> horzGenGid(gpuNum);

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

	thrust::host_vector<thrust::host_vector<float> > subVol(gpuNum);

	std::vector<size_t> siz(gpuNum,0);
	std::vector<size_t> nsiz(gpuNum,0);

	omp_set_num_threads(gpuNum);
#pragma omp parallel for
	for(int i = 0; i < gpuNum; ++i)
	{
		// get projection view index pair
		getPrjIdxPair<float>(hzPos, ObjZIdx_Start[i], ObjZIdx_End[i],
						objCntIdxZ, dz, ZN, detCntIdxV, detStpZ, DNV,
						prjIdx_Start[i], prjIdx_End[i]);
		SPN[i] = prjIdx_End[i] - prjIdx_Start[i];
		//std::cout<<i<<" "<<prjIdx_Start[i]<<" "<<prjIdx_End[i]<<"\n";
		// Calculate the corresponding center position index of the sub volumes
		subImgZCenter[i] = -imgZCenter / dz + ZN * 0.5 - ObjZIdx_Start[i] - 0.5f; // index position

		CUDA_CHECK_RETURN(cudaSetDevice(i));
		CUDA_CHECK_RETURN(cudaStreamCreate(&stream[i]));

		// Generate the SAT for the projection data
		siz[i] = DNU * DNV * SPN[i];
		nsiz[i] = (DNU + 1) * (DNV + 1) * SPN[i];
		prjSAT[i].resize(nsiz[i]);
		prj[i].resize(siz[i]);
		thrust::copy(
				hprj + DNU * DNV * prjIdx_Start[i],
				hprj + DNU * DNV * prjIdx_End[i],
				prj[i].begin());

		copyGid[i].x = (DNV + copyBlk.x - 1) / copyBlk.x;
		copyGid[i].y = (DNU + copyBlk.y - 1) / copyBlk.y;
		copyGid[i].z = (SPN[i] + copyBlk.z - 1) / copyBlk.z;
		addOneSidedZeroBoarder<<<copyGid[i], copyBlk, 0, stream[i]>>>(
				thrust::raw_pointer_cast(&prj[i][0]),
				thrust::raw_pointer_cast(&prjSAT[i][0]),
				DNU, DNV, SPN[i]);
		//cudaStreamSynchronize(stream[i]);

		vertGenGid[i].x = (nDNU * SPN[i] + vertGenBlk.x - 1) / copyBlk.x;
		vertGenGid[i].y = 1;
		vertGenGid[i].z = 1;
		verticalIntegral2 << <vertGenGid[i], vertGenBlk, 0, stream[i] >> >(
			thrust::raw_pointer_cast(&prjSAT[i][0]), nDNV, nDNU * SPN[i]);

		horzGenGid[i].x = (nDNV + horzGenBlk.x - 1) / horzGenBlk.x;
		horzGenGid[i].y = (SPN[i] + horzGenBlk.y - 1) / horzGenBlk.y;
		horzGenGid[i].z = 1;
		heorizontalIntegral2 << <horzGenGid[i], horzGenBlk,0,stream[i] >> >(
			thrust::raw_pointer_cast(&prjSAT[i][0]), nDNU, nDNV, SPN[i]);

		prj[i].clear();

		cudaExtent prjSize;
		prjSize.width = DNV + 1;
		prjSize.height = DNU + 1;
		prjSize.depth = SPN[i];
		CUDA_CHECK_RETURN(cudaMalloc3DArray(&d_prjArray[i], &channelDesc, prjSize));

		cudaMemcpy3DParms copyParams = { 0 };
			copyParams.srcPtr = make_cudaPitchedPtr(
				(void*) thrust::raw_pointer_cast(&prjSAT[i][0]),
				prjSize.width * sizeof(float),
				prjSize.width, prjSize.height);
			copyParams.dstArray = d_prjArray[i];
			copyParams.extent = prjSize;
			copyParams.kind = cudaMemcpyDeviceToDevice;
		CUDA_CHECK_RETURN(cudaMemcpy3DAsync(&copyParams,stream[i]));

		cudaResourceDesc resDesc;
		memset(&resDesc, 0, sizeof(resDesc));
		resDesc.resType = cudaResourceTypeArray;
		resDesc.res.array.array = d_prjArray[i];
		cudaTextureDesc texDesc;
		memset(&texDesc, 0, sizeof(texDesc));
		texDesc.addressMode[0] = cudaAddressModeClamp;
		texDesc.addressMode[1] = cudaAddressModeClamp;
		texDesc.addressMode[2] = cudaAddressModeClamp;
		texDesc.filterMode = cudaFilterModeLinear;
		texDesc.readMode = cudaReadModeElementType;
		texDesc.normalizedCoords = false;
		CUDA_CHECK_RETURN(cudaCreateTextureObject(&texObj[i], &resDesc, &texDesc, nullptr));
		prjSAT[i].clear();
		// The part above are for branchless DD

		gid[i].x = (SZN[i] + blk.x - 1) / blk.x;
		gid[i].y = (XN + blk.y - 1) / blk.y;
		gid[i].z = (YN + blk.z - 1) / blk.z;

		vol[i].resize(XN * YN * SZN[i]);
		msk[i].resize(XN * YN);
		thrust::copy(mask, mask + XN * YN, msk[i].begin());

		cossinZT[i].resize(SPN[i]);
		thrust::copy(
				COSSINZT.begin() + prjIdx_Start[i],
				COSSINZT.begin() + prjIdx_End[i],
				cossinZT[i].begin());
	}
#pragma omp parallel for
	for(int i = 0; i < gpuNum; ++i)
	{

		CUDA_CHECK_RETURN(cudaSetDevice(i));
		DD3_gpu_back_ker<_BRANCHLESS> << <gid[i], blk, 0, stream[i] >> >(texObj[i],
				thrust::raw_pointer_cast(&vol[i][0]), thrust::raw_pointer_cast(&msk[i][0]),
				thrust::raw_pointer_cast(&cossinZT[i][0]), make_float3(x0, y0, z0), S2D,
				make_float3(objCntIdxX, objCntIdxY, subImgZCenter[i]), //  have to be changed
				dx, dz, detParas.z, detParas.w, make_float2(detParas.x, detParas.y),
				make_int3(XN, YN, SZN[i]), SPN[i], 0);
	}
#pragma omp barrier
#pragma omp parallel for
	for(int i = 0 ;i < gpuNum; ++i)
	{
		CUDA_CHECK_RETURN(cudaSetDevice(i));
		// copy the volume back.
		subVol[i].resize(XN * YN * SZN[i]);
		thrust::copy(vol[i].begin(), vol[i].end(), subVol[i].begin());

		vol[i].clear();
		msk[i].clear();
		cossinZT[i].clear();

		CUDA_CHECK_RETURN(cudaDestroyTextureObject(texObj[i]));
		CUDA_CHECK_RETURN(cudaFreeArray(d_prjArray[i]));
	}
	CUDA_CHECK_RETURN(cudaDeviceSynchronize());

	combineVolume<float>(hvol, XN, YN, ZN, subVol, &(SZN[0]), gpuNum);
}





void DD3Back_pseudo_multiGPU(
	float x0, float y0, float z0,
	int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* h_angs, float* h_zPos, int PN,
	int XN, int YN, int ZN,
	float* hvol, float* hprj,
	float dx, float dz,
	byte* mask, int* startVOL, int gpuNum)
{
	thrust::host_vector<float> hangs(h_angs, h_angs + PN);
	thrust::host_vector<float> hzPos(h_zPos, h_zPos + PN);

	std::vector<int> ObjZIdx_Start(startVOL, startVOL + gpuNum);
	std::vector<int> ObjZIdx_End(ObjZIdx_Start.size());
	std::copy(ObjZIdx_Start.begin() + 1, ObjZIdx_Start.end(), ObjZIdx_End.begin());
	ObjZIdx_End[gpuNum - 1] = ZN;

	std::vector<int> prjIdx_Start(gpuNum);
	std::vector<int> prjIdx_End(gpuNum);

	const float objCntIdxZ = (ZN - 1.0f) * 0.5 - imgZCenter / dz;
	const float detStpZ = (zds[DNV - 1] - zds[0]) / (DNV - 1.0f); // detector cell height
	const float detCntIdxV = -zds[0] / detStpZ; // Detector Center along Z direction

	std::vector<int> SZN = ObjZIdx_End - ObjZIdx_Start; // sub volume slices number

	std::vector<float> subImgZCenter(gpuNum,0.0f);
	std::vector<int> SPN(gpuNum);

	const float objCntIdxX = (XN - 1.0f) * 0.5f - imgXCenter / dx;
	const float objCntIdxY = (YN - 1.0f) * 0.5f - imgYCenter / dx;

	std::vector<float3> sour(gpuNum);
	thrust::host_vector<thrust::device_vector<byte> > msk(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > vol(gpuNum);
	thrust::host_vector<thrust::device_vector<float3> > cossinZT(gpuNum);
	thrust::host_vector<cudaArray*> d_prjArray(gpuNum);
	thrust::host_vector<cudaTextureObject_t> texObj(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > prj(gpuNum);
	thrust::host_vector<cudaStream_t> stream(gpuNum);

	const float4 detParas = calDetParas(xds, yds, zds, x0, y0, z0, DNU, DNV);
	const float S2D = hypotf(xds[0] - x0, yds[0] - y0);

	// Pre calculate the cossin z positions
	thrust::device_vector<float3> COSSINZT(PN);
	thrust::device_vector<float> ANGS = hangs;
	thrust::device_vector<float> ZPOS = hzPos;
	thrust::transform(
		thrust::make_zip_iterator(thrust::make_tuple(ANGS.begin(), ZPOS.begin())),
		thrust::make_zip_iterator(thrust::make_tuple(ANGS.end(), ZPOS.end())),
		COSSINZT.begin(), CTMBIR::ConstantForBackProjection4(x0, y0, z0));

	dim3 copyBlk(64,16,1);
	thrust::host_vector<dim3> copyGid(gpuNum);
	dim3 blk(BACK_BLKX, BACK_BLKY, BACK_BLKZ);
	thrust::host_vector<dim3> gid(gpuNum);
	dim3 vertGenBlk(512,1,1);
	thrust::host_vector<dim3> vertGenGid(gpuNum);
	dim3 horzGenBlk(64,16,1);
	thrust::host_vector<dim3> horzGenGid(gpuNum);

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

	thrust::host_vector<thrust::host_vector<float> > subVol(gpuNum);

	std::vector<size_t> siz(gpuNum,0);
	std::vector<size_t> nsiz(gpuNum,0);

	omp_set_num_threads(gpuNum);
#pragma omp parallel for
	for(int i = 0; i < gpuNum; ++i)
	{
		// get projection view index pair
		getPrjIdxPair<float>(hzPos, ObjZIdx_Start[i], ObjZIdx_End[i],
						objCntIdxZ, dz, ZN, detCntIdxV, detStpZ, DNV,
						prjIdx_Start[i], prjIdx_End[i]);
		SPN[i] = prjIdx_End[i] - prjIdx_Start[i];
		//std::cout<<i<<" "<<prjIdx_Start[i]<<" "<<prjIdx_End[i]<<"\n";
		// Calculate the corresponding center position index of the sub volumes
		subImgZCenter[i] = -imgZCenter / dz + ZN * 0.5 - ObjZIdx_Start[i] - 0.5f; // index position

		CUDA_CHECK_RETURN(cudaSetDevice(i));
		CUDA_CHECK_RETURN(cudaStreamCreate(&stream[i]));
		////////////////////////////////////////////////////////////////////////
		siz[i] = DNU * DNV * SPN[i];
		prj[i].resize(siz[i]);
		thrust::copy(
				hprj + DNU * DNV * prjIdx_Start[i],
				hprj + DNU * DNV * prjIdx_End[i],
				prj[i].begin());

		cudaExtent prjSize;
		prjSize.width = DNV;
		prjSize.height = DNU;
		prjSize.depth = SPN[i];
		CUDA_CHECK_RETURN(cudaMalloc3DArray(&d_prjArray[i], &channelDesc, prjSize));

		cudaMemcpy3DParms copyParams = { 0 };
			copyParams.srcPtr = make_cudaPitchedPtr(
				(void*) thrust::raw_pointer_cast(&prj[i][0]),
				prjSize.width * sizeof(float),
				prjSize.width, prjSize.height);
			copyParams.dstArray = d_prjArray[i];
			copyParams.extent = prjSize;
			copyParams.kind = cudaMemcpyDeviceToDevice;
		CUDA_CHECK_RETURN(cudaMemcpy3DAsync(&copyParams,stream[i]));

		cudaResourceDesc resDesc;
		memset(&resDesc, 0, sizeof(resDesc));
		resDesc.resType = cudaResourceTypeArray;
		resDesc.res.array.array = d_prjArray[i];
		cudaTextureDesc texDesc;
		memset(&texDesc, 0, sizeof(texDesc));
		texDesc.addressMode[0] = cudaAddressModeBorder;
		texDesc.addressMode[1] = cudaAddressModeBorder;
		texDesc.addressMode[2] = cudaAddressModeBorder;
		texDesc.filterMode = cudaFilterModeLinear;
		texDesc.readMode = cudaReadModeElementType;
		texDesc.normalizedCoords = false;
		CUDA_CHECK_RETURN(cudaCreateTextureObject(&texObj[i], &resDesc, &texDesc, nullptr));
		prj[i].clear();
		////////////////////////////////////////////////////////////////////////
		// Generate the SAT for the projection data
		// The part above are for branchless DD

		gid[i].x = (SZN[i] + blk.x - 1) / blk.x;
		gid[i].y = (XN + blk.y - 1) / blk.y;
		gid[i].z = (YN + blk.z - 1) / blk.z;

		vol[i].resize(XN * YN * SZN[i]);
		msk[i].resize(XN * YN);
		thrust::copy(mask, mask + XN * YN, msk[i].begin());

		cossinZT[i].resize(SPN[i]);
		thrust::copy(
				COSSINZT.begin() + prjIdx_Start[i],
				COSSINZT.begin() + prjIdx_End[i],
				cossinZT[i].begin());
	}
#pragma omp parallel for
	for(int i = 0; i < gpuNum; ++i)
	{
		CUDA_CHECK_RETURN(cudaSetDevice(i));
		DD3_gpu_back_ker<_PSEUDODD> << <gid[i], blk, 0, stream[i] >> >(texObj[i],
				thrust::raw_pointer_cast(&vol[i][0]), thrust::raw_pointer_cast(&msk[i][0]),
				thrust::raw_pointer_cast(&cossinZT[i][0]), make_float3(x0, y0, z0), S2D,
				make_float3(objCntIdxX, objCntIdxY, subImgZCenter[i]), //  have to be changed
				dx, dz, detParas.z, detParas.w, make_float2(detParas.x, detParas.y),
				make_int3(XN, YN, SZN[i]), SPN[i], 0);
	}
#pragma omp barrier
#pragma omp parallel for
	for (int i = 0; i < gpuNum; ++i)
	{
		// copy the volume back.
		subVol[i].resize(XN * YN * SZN[i]);
		thrust::copy(vol[i].begin(), vol[i].end(), subVol[i].begin());

		vol[i].clear();
		msk[i].clear();
		cossinZT[i].clear();

		CUDA_CHECK_RETURN(cudaDestroyTextureObject(texObj[i]));
		CUDA_CHECK_RETURN(cudaFreeArray(d_prjArray[i]));
	}
	CUDA_CHECK_RETURN(cudaDeviceSynchronize());

	combineVolume<float>(hvol, XN, YN, ZN, subVol, &(SZN[0]), gpuNum);
}




extern "C"
void DD3Back_multiGPU(
	float x0, float y0, float z0,
	int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN,
	float* hvol, float* hprj,
	float dx, float dz,
	byte* mask, int bakMode, int* startVOL, int gpuNum)
{
	switch(bakMode)
	{
	case 0: // Branchless backprojection
		DD3Back_branchless_sat2d_multiGPU(x0, y0, z0,
			DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
			hangs, hzPos, PN, XN, YN, ZN, hvol, hprj,
			dx, dz, mask, startVOL, gpuNum);
		break;
	case 1: // Volume Rendering backprojection
	case 2: // Pseudo DD backprojection
		DD3Back_pseudo_multiGPU(x0, y0, z0,
			DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
			hangs, hzPos, PN, XN, YN, ZN, hvol, hprj,
			dx, dz, mask, startVOL, gpuNum);
		break;
	case 3: // Z line branchless backprojection
	default: // branchless DD backprojection
		DD3Back_branchless_sat2d_multiGPU(x0, y0, z0,
			DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
			hangs, hzPos, PN, XN, YN, ZN, hvol, hprj,
			dx, dz, mask, startVOL, gpuNum);
		break;
	}
}





void CT::Back(std::vector<float>& hvol, std::vector<float>& hprj, Geometry geo, const std::string& backModel)
{
	std::vector<float> xds = geo.getXds();
	std::vector<float> yds = geo.getYds();
	std::vector<float> zds = geo.getZds();
	std::vector<float> hangs = geo.getAllAngles();
	std::vector<float> hzPos = geo.getAllZPoses();
	
	std::vector<unsigned char> mask(geo.getObjDimX() * geo.getObjDimY(), 1);
	DD3Back_gpu(0, geo.getSourToObj(), 0,
		geo.getDetNumWidth(), geo.getDetNumHeight(),
		&xds[0], &yds[0], &zds[0],
		geo.getObjCtrCoordX(), geo.getObjCtrCoordY(), geo.getObjCtrCoordZ(),
		&hangs[0], &hzPos[0], geo.getViewNumber(),
		geo.getObjDimX(), geo.getObjDimY(), geo.getObjDimZ(),
		&hvol[0], &hprj[0],
		geo.getVoxelSizeX(), geo.getVoxelSizeZ(),
		&mask[0], 0, 0, 0);
}