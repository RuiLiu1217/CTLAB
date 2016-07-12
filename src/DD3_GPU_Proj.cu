#include "utilities.cuh"
#include "DD3_GPU_Proj.h"

#define BLKX 32
#define BLKY 8
#define BLKZ 1


template<typename Ta, typename Tb>
__global__ void naive_copyToTwoVolumes(Ta* in_ZXY,
	Tb* out_ZXY, Tb* out_ZYX,
	int XN, int YN, int ZN)
{
	int idz = threadIdx.x + blockIdx.x * blockDim.x;
	int idx = threadIdx.y + blockIdx.y * blockDim.y;
	int idy = threadIdx.z + blockIdx.z * blockDim.z;
	if (idx < XN && idy < YN && idz < ZN)
	{
		int i = (idy * XN + idx) * ZN + idz;
		int ni = (idy * (XN + 1) + (idx + 1)) * (ZN + 1) + idz + 1;
		int nj = (idx * (YN + 1) + (idy + 1)) * (ZN + 1) + idz + 1;

		out_ZXY[ni] = in_ZXY[i];
		out_ZYX[nj] = in_ZXY[i];
	}
}

template<typename Ta, typename Tb>
__global__ void naive_herizontalIntegral(Ta* in, Tb* out, int N, int ZN)
{
	int zi = threadIdx.x + blockIdx.x * blockDim.x;
	if (zi < ZN)
	{
		out[zi] = in[zi];
		for (int i = 1; i < N; ++i)
		{
			out[i * ZN + zi] = out[(i - 1) * ZN + zi]
				+ in[i * ZN + zi];
		}
	}
}

template<typename Ta, typename Tb>
__global__ void naive_verticalIntegral(Ta* in, Tb* out, int N, int ZN)
{
	int xyi = threadIdx.x + blockIdx.x * blockDim.x;
	if (xyi < N)
	{
		out[xyi * ZN] = in[xyi * ZN];
		for (int ii = 1; ii < ZN; ++ii)
		{
			out[xyi * ZN + ii] = out[xyi * ZN + ii - 1]
				+ in[xyi * ZN + ii];
		}

	}
}


template<typename T>
__global__ void verticalIntegral(T* prj, int ZN, int N)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < N)
	{
		int currentHead = idx * ZN;
		for (int ii = 1; ii < ZN; ++ii)
		{
			prj[currentHead + ii] = prj[currentHead + ii] + prj[currentHead + ii - 1];
		}
	}
}

template<typename T>
__global__ void horizontalIntegral(T* prj, int DNU, int DNV, int PN)
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





__global__ void naive_vertialIntegral(double* in, int2* out, int N, int ZN)
{
	int xyi = threadIdx.x + blockIdx.x * blockDim.x;
	if (xyi < N)
	{
		double temp = in[xyi * ZN];
		out[xyi * ZN] = make_int2(__double2loint(temp), __double2hiint(temp));
		double temp2 = 0;
		for (int ii = 0; ii < ZN; ++ii)
		{
			temp2 = temp + in[xyi * ZN + ii];
			out[xyi * ZN + ii] = make_int2(__double2loint(temp2), __double2hiint(temp2));
			temp = temp2;
		}
	}
}



__global__ void verticalIntegral(float* prj, int ZN, int N)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < N)
	{
		int currentHead = idx * ZN;
		for (int ii = 1; ii < ZN; ++ii)
		{
			prj[currentHead + ii] = prj[currentHead + ii] + prj[currentHead + ii - 1];
		}
	}
}



__global__ void horizontalIntegral(float* prj, int DNU, int DNV, int PN)
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

void genSAT_fof_Volume(float* hvol,
	thrust::device_vector<float>&ZXY,
	thrust::device_vector<float>&ZYX,
	int XN, int YN, int ZN)
{
	const int siz = XN * YN * ZN;
	const int nsiz_ZXY = (ZN + 1) * (XN + 1) * YN;
	const int nsiz_ZYX = (ZN + 1) * (YN + 1) * XN;
	ZXY.resize(nsiz_ZXY);
	ZYX.resize(nsiz_ZYX);

	thrust::device_vector<float> vol(hvol, hvol + siz);

	dim3 blk(64, 16, 1);
	dim3 gid(
		(ZN + blk.x - 1) / blk.x,
		(XN + blk.y - 1) / blk.y,
		(YN + blk.z - 1) / blk.z);

	naive_copyToTwoVolumes << <gid, blk >> >(
		thrust::raw_pointer_cast(&vol[0]),
		thrust::raw_pointer_cast(&ZXY[0]),
		thrust::raw_pointer_cast(&ZYX[0]),
		XN, YN, ZN);

	vol.clear();
	const int nZN = ZN + 1;
	const int nXN = XN + 1;
	const int nYN = YN + 1;

	blk.x = 32;
	blk.y = 1;
	blk.z = 1;
	gid.x = (nXN * YN + blk.x - 1) / blk.x;
	gid.y = 1;
	gid.z = 1;
	verticalIntegral << <gid, blk >> >(
		thrust::raw_pointer_cast(&ZXY[0]),
		nZN, nXN * YN);

	blk.x = 64;
	blk.y = 16;
	blk.z = 1;
	gid.x = (nZN + blk.x - 1) / blk.x;
	gid.y = (YN + blk.y - 1) / blk.y;
	gid.z = 1;

	horizontalIntegral << <gid, blk >> >(
		thrust::raw_pointer_cast(&ZXY[0]),
		nXN, nZN, YN);

	blk.x = 32;
	blk.y = 1;
	blk.z = 1;
	gid.x = (nYN * XN + blk.x - 1) / blk.x;
	gid.y = 1;
	gid.z = 1;
	verticalIntegral << <gid, blk >> >(
		thrust::raw_pointer_cast(&ZYX[0]),
		nZN, nXN * YN);

	blk.x = 64;
	blk.y = 16;
	blk.z = 1;
	gid.x = (nZN + blk.x - 1) / blk.x;
	gid.y = (XN + blk.y - 1) / blk.y;
	gid.z = 1;

	horizontalIntegral << <gid, blk >> >(
		thrust::raw_pointer_cast(&ZYX[0]),
		nYN, nZN, XN);
}


void genSAT_fof_Volume_alreadyinGPU(
	const thrust::device_vector<float>& vol,
	thrust::device_vector<float>&ZXY,
	thrust::device_vector<float>&ZYX,
	int XN, int YN, int ZN)
{
	//const int siz = XN * YN * ZN;
	const int nsiz_ZXY = (ZN + 1) * (XN + 1) * YN;
	const int nsiz_ZYX = (ZN + 1) * (YN + 1) * XN;
	ZXY.resize(nsiz_ZXY);
	ZYX.resize(nsiz_ZYX);

	//thrust::device_vector<float> vol(hvol, hvol + siz);

	dim3 blk(64, 16, 1);
	dim3 gid(
		(ZN + blk.x - 1) / blk.x,
		(XN + blk.y - 1) / blk.y,
		(YN + blk.z - 1) / blk.z);

	naive_copyToTwoVolumes << <gid, blk >> >(
		thrust::raw_pointer_cast(&vol[0]),
		thrust::raw_pointer_cast(&ZXY[0]),
		thrust::raw_pointer_cast(&ZYX[0]),
		XN, YN, ZN);

	//vol.clear();
	const int nZN = ZN + 1;
	const int nXN = XN + 1;
	const int nYN = YN + 1;

	blk.x = 32;
	blk.y = 1;
	blk.z = 1;
	gid.x = (nXN * YN + blk.x - 1) / blk.x;
	gid.y = 1;
	gid.z = 1;
	verticalIntegral << <gid, blk >> >(
		thrust::raw_pointer_cast(&ZXY[0]),
		nZN, nXN * YN);

	blk.x = 64;
	blk.y = 16;
	blk.z = 1;
	gid.x = (nZN + blk.x - 1) / blk.x;
	gid.y = (YN + blk.y - 1) / blk.y;
	gid.z = 1;

	horizontalIntegral << <gid, blk >> >(
		thrust::raw_pointer_cast(&ZXY[0]),
		nXN, nZN, YN);

	blk.x = 32;
	blk.y = 1;
	blk.z = 1;
	gid.x = (nYN * XN + blk.x - 1) / blk.x;
	gid.y = 1;
	gid.z = 1;
	verticalIntegral << <gid, blk >> >(
		thrust::raw_pointer_cast(&ZYX[0]),
		nZN, nXN * YN);

	blk.x = 64;
	blk.y = 16;
	blk.z = 1;
	gid.x = (nZN + blk.x - 1) / blk.x;
	gid.y = (XN + blk.y - 1) / blk.y;
	gid.z = 1;

	horizontalIntegral << <gid, blk >> >(
		thrust::raw_pointer_cast(&ZYX[0]),
		nYN, nZN, XN);
}

__global__  void DD3_gpu_proj_branchless_sat2d_ker(
	cudaTextureObject_t volTex1,
	cudaTextureObject_t volTex2,
	float* proj,
	float3 s,
	const float3* __restrict__ cossinZT,
	const float* __restrict__ xds,
	const float* __restrict__ yds,
	const float* __restrict__ zds,
	const float* __restrict__ bxds,
	const float* __restrict__ byds,
	const float* __restrict__ bzds,
	float3 objCntIdx,
	float dx, float dz,
	int XN, int YN, int ZN,
	int DNU, int DNV, int PN)
{
	int detIdV = threadIdx.x + blockIdx.x * blockDim.x;
	int detIdU = threadIdx.y + blockIdx.y * blockDim.y;
	int angIdx = threadIdx.z + blockIdx.z * blockDim.z;
	__shared__ float _xds[BLKY];
	__shared__ float _yds[BLKY];
	_xds[threadIdx.y] = xds[detIdU];
	_yds[threadIdx.y] = yds[detIdU];
	__syncthreads();
	if (detIdU < DNU && detIdV < DNV && angIdx < PN)
	{
		float3 dir = cossinZT[angIdx];
		float3 cursour = make_float3(
			s.x * dir.x - s.y * dir.y,
			s.x * dir.y + s.y * dir.x,
			s.z + dir.z);
		s = cossinZT[angIdx];
		float summ = _xds[threadIdx.y] * s.x - _yds[threadIdx.y] * s.y;
		float obj = _xds[threadIdx.y] * s.y + _yds[threadIdx.y] * s.x;
		float realL = bxds[detIdU];
		float realR = byds[detIdU];
		float realU = bxds[detIdU + 1];
		float realD = byds[detIdU + 1];

		float2 curDetL = make_float2(
			realL * s.x - realR * s.y,
			realL * s.y + realR * s.x);
		float2 curDetR = make_float2(
			realU * s.x - realD * s.y,
			realU * s.y + realD * s.x);

		float4 curDet = make_float4(
			summ, obj, bzds[detIdV] + s.z,
			bzds[detIdV + 1] + s.z);

		dir = normalize(make_float3(summ, obj,
			zds[detIdV] + s.z) - cursour);

		summ = 0;
		obj = 0;
		float intersectLength, intersectHeight;
		float invdz = 1.0 / dz;
		float invdx = 1.0 / dx;

		float factL(1.0f);
		float factR(1.0f);
		float factU(1.0f);
		float factD(1.0f);

		float constVal = 0;
		if (fabsf(s.x) <= fabsf(s.y))
		{
			summ = 0;
			factL = (curDetL.y - cursour.y) / (curDetL.x - cursour.x);
			factR = (curDetR.y - cursour.y) / (curDetR.x - cursour.x);
			factU = (curDet.w - cursour.z) / (curDet.x - cursour.x);
			factD = (curDet.z - cursour.z) / (curDet.x - cursour.x);

			constVal = dx * dx * dz / fabsf(dir.x);
#pragma unroll
			for (int ii = 0; ii < XN; ++ii)
			{
				obj = (ii - objCntIdx.x) * dx;

				realL = (obj - curDetL.x) * factL + curDetL.y;
				realR = (obj - curDetR.x) * factR + curDetR.y;
				realU = (obj - curDet.x) * factU + curDet.w;
				realD = (obj - curDet.x) * factD + curDet.z;

				intersectLength = realR - realL;
				intersectHeight = realU - realD;


				realD = realD * invdz + objCntIdx.z + 1;
				realR = realR * invdx + objCntIdx.y + 1;
				realU = realU * invdz + objCntIdx.z + 1;
				realL = realL * invdx + objCntIdx.y + 1;

				summ +=
					(
					tex3D<float>(volTex2, realD, realL, ii + 0.5)
					+ tex3D<float>(volTex2, realU, realR, ii + 0.5)
					- tex3D<float>(volTex2, realU, realL, ii + 0.5)
					- tex3D<float>(volTex2, realD, realR, ii + 0.5)
					) / (intersectLength * intersectHeight);

			}
			__syncthreads();
			proj[(angIdx * DNU + detIdU) * DNV + detIdV] = summ * constVal;
		}
		else
		{
			summ = 0;
			factL = (curDetL.x - cursour.x) / (curDetL.y - cursour.y);
			factR = (curDetR.x - cursour.x) / (curDetR.y - cursour.y);
			factU = (curDet.w - cursour.z) / (curDet.y - cursour.y);
			factD = (curDet.z - cursour.z) / (curDet.y - cursour.y);

			constVal = dx * dx * dz / fabsf(dir.y);
#pragma unroll
			for (int jj = 0; jj < YN; ++jj)
			{
				obj = (jj - objCntIdx.y) * dx;

				realL = (obj - curDetL.y) * factL + curDetL.x;
				realR = (obj - curDetR.y) * factR + curDetR.x;
				realU = (obj - curDet.y) * factU + curDet.w;
				realD = (obj - curDet.y) * factD + curDet.z;

				intersectLength = realR - realL;
				intersectHeight = realU - realD;


				realD = realD * invdz + objCntIdx.z + 1;
				realR = realR * invdx + objCntIdx.x + 1;
				realU = realU * invdz + objCntIdx.z + 1;
				realL = realL * invdx + objCntIdx.x + 1;

				summ +=
					(
					tex3D<float>(volTex1, realD, realL, jj + 0.5)
					+ tex3D<float>(volTex1, realU, realR, jj + 0.5)
					- tex3D<float>(volTex1, realU, realL, jj + 0.5)
					- tex3D<float>(volTex1, realD, realR, jj + 0.5)
					) / (intersectLength * intersectHeight);

			}
			__syncthreads();
			proj[(angIdx * DNU + detIdU) * DNV + detIdV] = summ * constVal;
		}
	}
}


void DD3_gpu_proj_branchless_sat2d(
	float x0, float y0, float z0,
	int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN,
	float* vol, float* hprj, float dx, float dz, byte* mask, int gpunum)
{
	for (int i = 0; i != XN * YN; ++i)
	{
		byte v = mask[i];
		for (int z = 0; z != ZN; ++z)
		{
			vol[i * ZN + z] = vol[i * ZN + z] * v;
		}
	}

	CUDA_SAFE_CALL(cudaSetDevice(gpunum));
	cudaDeviceReset();

	float* bxds = new float[DNU + 1];
	float* byds = new float[DNU + 1];
	float* bzds = new float[DNV + 1];

	DD3Boundaries(DNU + 1, xds, bxds);
	DD3Boundaries(DNU + 1, yds, byds);
	DD3Boundaries(DNV + 1, zds, bzds);

//	cudaStream_t streams[4];
//	cudaStreamCreate(&streams[0]);
//	cudaStreamCreate(&streams[1]);
//	cudaStreamCreate(&streams[2]);
//	cudaStreamCreate(&streams[3]);

	float objCntIdxX = (XN - 1.0) * 0.5 - imgXCenter / dx;
	float objCntIdxY = (YN - 1.0) * 0.5 - imgYCenter / dx;
	float objCntIdxZ = (ZN - 1.0) * 0.5 - imgZCenter / dz;

	thrust::device_vector<float> SATZXY;
	thrust::device_vector<float> SATZYX;
	genSAT_fof_Volume(vol, SATZXY, SATZYX, XN, YN, ZN);

	cudaExtent volumeSize1;
	cudaExtent volumeSize2;
	volumeSize1.width = ZN + 1;
	volumeSize1.height = XN + 1;
	volumeSize1.depth = YN;

	volumeSize2.width = ZN + 1;
	volumeSize2.height = YN + 1;
	volumeSize2.depth = XN;

	cudaChannelFormatDesc channelDesc1 = cudaCreateChannelDesc<float>();
	cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<float>();

	cudaArray* d_volumeArray1;
	cudaArray* d_volumeArray2;

	cudaMalloc3DArray(&d_volumeArray1, &channelDesc1, volumeSize1);
	cudaMalloc3DArray(&d_volumeArray2, &channelDesc2, volumeSize2);

	cudaMemcpy3DParms copyParams1 = { 0 };
	copyParams1.srcPtr = make_cudaPitchedPtr((void*)
		thrust::raw_pointer_cast(&SATZXY[0]),
		volumeSize1.width * sizeof(float),
		volumeSize1.width, volumeSize1.height);
	copyParams1.dstArray = d_volumeArray1;
	copyParams1.extent = volumeSize1;
	copyParams1.kind = cudaMemcpyDeviceToDevice;

	cudaMemcpy3DParms copyParams2 = { 0 };
	copyParams2.srcPtr = make_cudaPitchedPtr((void*)
		thrust::raw_pointer_cast(&SATZYX[0]),
		volumeSize2.width * sizeof(float),
		volumeSize2.width, volumeSize2.height);
	copyParams2.dstArray = d_volumeArray2;
	copyParams2.extent = volumeSize2;
	copyParams2.kind = cudaMemcpyDeviceToDevice;

	cudaMemcpy3D(&copyParams1);
	cudaMemcpy3D(&copyParams2);

	SATZXY.clear();
	SATZYX.clear();

	cudaResourceDesc resDesc1;
	cudaResourceDesc resDesc2;
	memset(&resDesc1, 0, sizeof(resDesc1));
	memset(&resDesc2, 0, sizeof(resDesc2));

	resDesc1.resType = cudaResourceTypeArray;
	resDesc2.resType = cudaResourceTypeArray;

	resDesc1.res.array.array = d_volumeArray1;
	resDesc2.res.array.array = d_volumeArray2;

	cudaTextureDesc texDesc1;
	cudaTextureDesc texDesc2;

	memset(&texDesc1, 0, sizeof(texDesc1));
	memset(&texDesc2, 0, sizeof(texDesc2));

	texDesc1.addressMode[0] = cudaAddressModeClamp;
	texDesc1.addressMode[1] = cudaAddressModeClamp;
	texDesc1.addressMode[2] = cudaAddressModeClamp;

	texDesc2.addressMode[0] = cudaAddressModeClamp;
	texDesc2.addressMode[1] = cudaAddressModeClamp;
	texDesc2.addressMode[2] = cudaAddressModeClamp;

	texDesc1.filterMode = cudaFilterModeLinear;
	texDesc2.filterMode = cudaFilterModeLinear;

	texDesc1.readMode = cudaReadModeElementType;
	texDesc2.readMode = cudaReadModeElementType;

	texDesc1.normalizedCoords = false;
	texDesc2.normalizedCoords = false;

	cudaTextureObject_t texObj1 = 0;
	cudaTextureObject_t texObj2 = 0;

	cudaCreateTextureObject(&texObj1, &resDesc1, &texDesc1, nullptr);
	cudaCreateTextureObject(&texObj2, &resDesc2, &texDesc2, nullptr);

	thrust::device_vector<float> prj(DNU * DNV * PN, 0);
	thrust::device_vector<float> d_xds(xds, xds + DNU);
	thrust::device_vector<float> d_yds(yds, yds + DNU);
	thrust::device_vector<float> d_zds(zds, zds + DNV);

	thrust::device_vector<float> d_bxds(bxds, bxds + DNU + 1);
	thrust::device_vector<float> d_byds(byds, byds + DNU + 1);
	thrust::device_vector<float> d_bzds(bzds, bzds + DNV + 1);

	thrust::device_vector<float> angs(hangs, hangs + PN);
	thrust::device_vector<float> zPos(hzPos, hzPos + PN);

	thrust::device_vector<float3> cossinZT(PN);
	thrust::transform(
		thrust::make_zip_iterator(thrust::make_tuple(angs.begin(), zPos.begin())),
		thrust::make_zip_iterator(thrust::make_tuple(angs.end(), zPos.end())),
		cossinZT.begin(), CTMBIR::ConstantForBackProjection4(x0, y0, z0));

	dim3 blk(BLKX, BLKY, BLKZ);
	dim3 gid(
		(DNV + blk.x - 1) / blk.x,
		(DNU + blk.y - 1) / blk.y,
		(PN + blk.z - 1) / blk.z);

	DD3_gpu_proj_branchless_sat2d_ker << <gid, blk >> >(texObj1, texObj2,
		thrust::raw_pointer_cast(&prj[0]),
		make_float3(x0, y0, z0),
		thrust::raw_pointer_cast(&cossinZT[0]),
		thrust::raw_pointer_cast(&d_xds[0]),
		thrust::raw_pointer_cast(&d_yds[0]),
		thrust::raw_pointer_cast(&d_zds[0]),
		thrust::raw_pointer_cast(&d_bxds[0]),
		thrust::raw_pointer_cast(&d_byds[0]),
		thrust::raw_pointer_cast(&d_bzds[0]),
		make_float3(objCntIdxX, objCntIdxY, objCntIdxZ),
		dx, dz, XN, YN, ZN, DNU, DNV, PN);

	thrust::copy(prj.begin(), prj.end(), hprj);
	cudaDestroyTextureObject(texObj1);
	cudaDestroyTextureObject(texObj2);

	cudaFreeArray(d_volumeArray1);
	cudaFreeArray(d_volumeArray2);

	prj.clear();
	angs.clear();
	zPos.clear();

	d_xds.clear();
	d_yds.clear();
	d_zds.clear();
	d_bxds.clear();
	d_byds.clear();
	d_bzds.clear();
	cossinZT.clear();

	delete[] bxds;
	delete[] byds;
	delete[] bzds;

}


// Use the texture to store the volume, we try to optimize the device function, because it is too slow when call "calSiddonOneRayKer" function
// The latency is caused by L2 cache in device function
__global__ void DD3ProjSiddon_gpu_kernel(
		cudaTextureObject_t volTex,
		float* proj,
		float3 s,
		const float3* __restrict__ cossinZT, //cosine, sine, Z shift in (x,y,z)
		const float* __restrict__ xds,
		const float* __restrict__ yds,
		const float* __restrict__ zds,
		float3 objCntIdx,
		float dx, float dz,
		int XN, int YN, int ZN,
		int DNU, int DNV, int PN)
{
	int detIdV = threadIdx.x + blockIdx.x * blockDim.x;
	int detIdU = threadIdx.y + blockIdx.y * blockDim.y;
	int angIdx = threadIdx.z + blockIdx.z * blockDim.z;
	__shared__ float3 MINO;
	__shared__ float3 _cossinZT;
	__shared__ float3 curS;
	MINO = make_float3(
			-(objCntIdx.x + 0.5f) * dx,
			-(objCntIdx.y + 0.5f) * dx,
			-(objCntIdx.z + 0.5f) * dz);
	_cossinZT = cossinZT[angIdx];

	__syncthreads();
	curS = make_float3(
				s.x * _cossinZT.x - s.y * _cossinZT.y,
			    s.x * _cossinZT.y + s.y * _cossinZT.x,
				s.z + _cossinZT.z);
	__syncthreads();
	if(detIdV < DNV && detIdU < DNU && angIdx < PN)
	{
		float3 curD = make_float3(
				xds[detIdU] * _cossinZT.x - yds[detIdU] * _cossinZT.y,
			    xds[detIdU] * _cossinZT.y + yds[detIdU] * _cossinZT.x,
			    zds[detIdV] + _cossinZT.z);
		float totWeight = 0;
		float prj = calSiddonOneRayKer(
			curS.x, curS.y, curS.z, curD.x, curD.y, curD.z,
			MINO.x, MINO.y, MINO.z, dx, dx, dz, XN,YN,ZN,volTex,&totWeight);

		proj[(angIdx * DNU + detIdU) * DNV + detIdV] = prj;
	}

}


// Ray-driven algorithm
void DD3ProjSiddon_gpu(
float x0, float y0, float z0,
int DNU, int DNV,
float* xds, float* yds, float* zds,
float imgXCenter, float imgYCenter, float imgZCenter,
float* hangs, float* hzPos, int PN,
int XN, int YN, int ZN,
float* vol, float* hprj,
float dx, float dz,
byte* mask, int gpunum)
{
	for (int i = 0; i != XN * YN; ++i)
	{
		byte v = mask[i];
		for (int z = 0; z != ZN; ++z)
		{
			vol[i * ZN + z] = vol[i * ZN + z] * v;
		}
	}

	CUDA_SAFE_CALL(cudaSetDevice(gpunum));
	cudaDeviceReset();

	float objCntIdxX = (XN - 1.0f) * 0.5f - imgXCenter / dx;
	float objCntIdxY = (YN - 1.0f) * 0.5f - imgYCenter / dx;
	float objCntIdxZ = (ZN - 1.0f) * 0.5f - imgZCenter / dz;

	cudaExtent volumeSize;
	volumeSize.width = ZN;
	volumeSize.height = XN;
	volumeSize.depth = YN;

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	cudaArray* d_volumeArray;
	cudaMalloc3DArray(&d_volumeArray, &channelDesc, volumeSize);

	cudaMemcpy3DParms copyParams = { 0 };
	copyParams.srcPtr = make_cudaPitchedPtr((void*)vol, // Changed here
		volumeSize.width * sizeof(float),
		volumeSize.width, volumeSize.height);
	copyParams.dstArray = d_volumeArray;
	copyParams.extent = volumeSize;
	copyParams.kind = cudaMemcpyHostToDevice;

	cudaMemcpy3D(&copyParams);

	cudaResourceDesc resDesc;

	memset(&resDesc, 0, sizeof(resDesc));

	resDesc.resType = cudaResourceTypeArray;
	resDesc.res.array.array = d_volumeArray;

	cudaTextureDesc texDesc;

	memset(&texDesc, 0, sizeof(texDesc));

	texDesc.addressMode[0] = cudaAddressModeBorder;
	texDesc.addressMode[1] = cudaAddressModeBorder;
	texDesc.addressMode[2] = cudaAddressModeBorder;

	texDesc.filterMode = cudaFilterModeLinear;

	texDesc.readMode = cudaReadModeElementType;

	texDesc.normalizedCoords = false;

	cudaTextureObject_t texObj = 0;

	cudaCreateTextureObject(&texObj, &resDesc, &texDesc, nullptr);

	thrust::device_vector<float> prj(DNU * DNV * PN, 0);
	thrust::device_vector<float> d_xds(xds, xds + DNU);
	thrust::device_vector<float> d_yds(yds, yds + DNU);
	thrust::device_vector<float> d_zds(zds, zds + DNV);

	thrust::device_vector<float> angs(hangs, hangs + PN);
	thrust::device_vector<float> zPos(hzPos, hzPos + PN);

	thrust::device_vector<float3> cossinZT(PN);
	thrust::transform(
		thrust::make_zip_iterator(thrust::make_tuple(angs.begin(), zPos.begin())),
		thrust::make_zip_iterator(thrust::make_tuple(angs.end(), zPos.end())),
		cossinZT.begin(), CTMBIR::ConstantForBackProjection4(x0, y0, z0));

	dim3 blk(BLKX, BLKY, BLKZ);
	dim3 gid(
		(DNV + blk.x - 1) / blk.x,
		(DNU + blk.y - 1) / blk.y,
		(PN + blk.z - 1) / blk.z);

	DD3ProjSiddon_gpu_kernel<<<gid, blk>>>(texObj,
			thrust::raw_pointer_cast(&prj[0]),
			make_float3(x0, y0, z0),
			thrust::raw_pointer_cast(&cossinZT[0]),
			thrust::raw_pointer_cast(&d_xds[0]),
			thrust::raw_pointer_cast(&d_yds[0]),
			thrust::raw_pointer_cast(&d_zds[0]),
			make_float3(objCntIdxX, objCntIdxY, objCntIdxZ),
			dx, dz, XN, YN, ZN, DNU, DNV, PN);


	thrust::copy(prj.begin(), prj.end(), hprj);
	cudaDestroyTextureObject(texObj);

	cudaFreeArray(d_volumeArray);


	prj.clear();
	angs.clear();
	zPos.clear();

	d_xds.clear();
	d_yds.clear();
	d_zds.clear();
	cossinZT.clear();

}






void DD3_gpu_proj_branchless_sat2d_alreadyinGPU(
	float x0, float y0, float z0,
	int DNU, int DNV,
	const thrust::device_vector<float>& xds,
	const thrust::device_vector<float>& yds,
	const thrust::device_vector<float>& zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	const thrust::device_vector<float>& hangs,const thrust::device_vector<float>& hzPos, int PN,
	int XN, int YN, int ZN,
	const thrust::device_vector<float>& vol,
	thrust::device_vector<float>& hprj, float dx, float dz,
	const thrust::device_vector<byte>& mask, int gpunum)
{
	//cudaSetDevice(gpunum);
	//cudaDeviceReset();

	thrust::host_vector<float> hxds = xds;
	thrust::host_vector<float> hyds = yds;
	thrust::host_vector<float> hzds = zds;

	thrust::host_vector<float> bxds(DNU+1,0);
	thrust::host_vector<float> byds(DNU+1,0);
	thrust::host_vector<float> bzds(DNV+1,0);


	DD3Boundaries(DNU + 1, &hxds[0], &bxds[0]);
	DD3Boundaries(DNU + 1, &hyds[0], &byds[0]);
	DD3Boundaries(DNV + 1, &hzds[0], &bzds[0]);

	float objCntIdxX = (XN - 1.0) * 0.5 - imgXCenter / dx;
	float objCntIdxY = (YN - 1.0) * 0.5 - imgYCenter / dx;
	float objCntIdxZ = (ZN - 1.0) * 0.5 - imgZCenter / dz;

	thrust::device_vector<float> SATZXY;
	thrust::device_vector<float> SATZYX;
	genSAT_fof_Volume_alreadyinGPU(vol, SATZXY, SATZYX, XN, YN, ZN);

	cudaExtent volumeSize1;
	cudaExtent volumeSize2;
	volumeSize1.width = ZN + 1;
	volumeSize1.height = XN + 1;
	volumeSize1.depth = YN;

	volumeSize2.width = ZN + 1;
	volumeSize2.height = YN + 1;
	volumeSize2.depth = XN;

	cudaChannelFormatDesc channelDesc1 = cudaCreateChannelDesc<float>();
	cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<float>();

	cudaArray* d_volumeArray1;
	cudaArray* d_volumeArray2;

	cudaMalloc3DArray(&d_volumeArray1, &channelDesc1, volumeSize1);
	cudaMalloc3DArray(&d_volumeArray2, &channelDesc2, volumeSize2);

	cudaMemcpy3DParms copyParams1 = { 0 };
	copyParams1.srcPtr = make_cudaPitchedPtr((void*)
		thrust::raw_pointer_cast(&SATZXY[0]),
		volumeSize1.width * sizeof(float),
		volumeSize1.width, volumeSize1.height);
	copyParams1.dstArray = d_volumeArray1;
	copyParams1.extent = volumeSize1;
	copyParams1.kind = cudaMemcpyDeviceToDevice;

	cudaMemcpy3DParms copyParams2 = { 0 };
	copyParams2.srcPtr = make_cudaPitchedPtr((void*)
		thrust::raw_pointer_cast(&SATZYX[0]),
		volumeSize2.width * sizeof(float),
		volumeSize2.width, volumeSize2.height);
	copyParams2.dstArray = d_volumeArray2;
	copyParams2.extent = volumeSize2;
	copyParams2.kind = cudaMemcpyDeviceToDevice;

	cudaMemcpy3D(&copyParams1);
	cudaMemcpy3D(&copyParams2);

	SATZXY.clear();
	SATZYX.clear();

	cudaResourceDesc resDesc1;
	cudaResourceDesc resDesc2;
	memset(&resDesc1, 0, sizeof(resDesc1));
	memset(&resDesc2, 0, sizeof(resDesc2));

	resDesc1.resType = cudaResourceTypeArray;
	resDesc2.resType = cudaResourceTypeArray;

	resDesc1.res.array.array = d_volumeArray1;
	resDesc2.res.array.array = d_volumeArray2;

	cudaTextureDesc texDesc1;
	cudaTextureDesc texDesc2;

	memset(&texDesc1, 0, sizeof(texDesc1));
	memset(&texDesc2, 0, sizeof(texDesc2));

	texDesc1.addressMode[0] = cudaAddressModeClamp;
	texDesc1.addressMode[1] = cudaAddressModeClamp;
	texDesc1.addressMode[2] = cudaAddressModeClamp;

	texDesc2.addressMode[0] = cudaAddressModeClamp;
	texDesc2.addressMode[1] = cudaAddressModeClamp;
	texDesc2.addressMode[2] = cudaAddressModeClamp;

	texDesc1.filterMode = cudaFilterModeLinear;
	texDesc2.filterMode = cudaFilterModeLinear;

	texDesc1.readMode = cudaReadModeElementType;
	texDesc2.readMode = cudaReadModeElementType;

	texDesc1.normalizedCoords = false;
	texDesc2.normalizedCoords = false;

	cudaTextureObject_t texObj1 = 0;
	cudaTextureObject_t texObj2 = 0;

	cudaCreateTextureObject(&texObj1, &resDesc1, &texDesc1, nullptr);
	cudaCreateTextureObject(&texObj2, &resDesc2, &texDesc2, nullptr);

	thrust::device_vector<float> d_bxds = bxds;
	thrust::device_vector<float> d_byds = byds;
	thrust::device_vector<float> d_bzds = bzds;

	thrust::device_vector<float3> cossinZT(PN);
	thrust::transform(
		thrust::make_zip_iterator(thrust::make_tuple(hangs.begin(), hzPos.begin())),
		thrust::make_zip_iterator(thrust::make_tuple(hangs.end(), hzPos.end())),
		cossinZT.begin(), CTMBIR::ConstantForBackProjection4(x0, y0, z0));

	dim3 blk(BLKX, BLKY, BLKZ);
	dim3 gid(
		(DNV + blk.x - 1) / blk.x,
		(DNU + blk.y - 1) / blk.y,
		(PN + blk.z - 1) / blk.z);

	DD3_gpu_proj_branchless_sat2d_ker << <gid, blk >> >(texObj1, texObj2,
		thrust::raw_pointer_cast(&hprj[0]),
		make_float3(x0, y0, z0),
		thrust::raw_pointer_cast(&cossinZT[0]),
		thrust::raw_pointer_cast(&xds[0]),
		thrust::raw_pointer_cast(&yds[0]),
		thrust::raw_pointer_cast(&zds[0]),
		thrust::raw_pointer_cast(&d_bxds[0]),	thrust::raw_pointer_cast(&d_byds[0]),
		thrust::raw_pointer_cast(&d_bzds[0]),
		make_float3(objCntIdxX, objCntIdxY, objCntIdxZ),
		dx, dz, XN, YN, ZN, DNU, DNV, PN);


	cudaDestroyTextureObject(texObj1);
	cudaDestroyTextureObject(texObj2);

	cudaFreeArray(d_volumeArray1);
	cudaFreeArray(d_volumeArray2);


	d_bxds.clear();
	d_byds.clear();
	d_bzds.clear();
	cossinZT.clear();

	hxds.clear();
	hyds.clear();
	hzds.clear();
	bxds.clear();
	byds.clear();
	bzds.clear();

}



__global__ void DD3_gpu_proj_volumerendering_ker(
	cudaTextureObject_t volTex,
	float* proj,
	float x0, float y0, float z0,
	float* d_xds, float* d_yds, float* d_zds,
	float2* dirCent,
	float* cosTs,
	float* sinTs,
	float objCntIdxX, float objCntIdxY, float objCntIdxZ,
	float dx, float dz,
	int XN, int YN, int ZN,
	int DNU, int DNV, int PN,
	float* angs, float* zPos, float3* gcursour, float* gcosGamma)
{
	int detIdV = threadIdx.x + blockIdx.x * blockDim.x;
	int detIdU = threadIdx.y + blockIdx.y * blockDim.y;
	int angIdx = threadIdx.z + blockIdx.z * blockDim.z;

	__shared__ float cosGamma[BLKX];
	__shared__ float2 sdirCent[BLKZ][BLKY];
	__shared__ float2 normDir[BLKZ][BLKY];
	__shared__ float curang;
	__shared__ float cosT[BLKZ];
	__shared__ float sinT[BLKZ];
	cosGamma[threadIdx.x] = gcosGamma[detIdV];
	sdirCent[threadIdx.z][threadIdx.y] = dirCent[angIdx * DNU + detIdU];
	float l1 = sqrtf(sdirCent[threadIdx.z][threadIdx.y].x * sdirCent[threadIdx.z][threadIdx.y].x + sdirCent[threadIdx.z][threadIdx.y].y * sdirCent[threadIdx.z][threadIdx.y].y);
	normDir[threadIdx.z][threadIdx.y].x = sdirCent[threadIdx.z][threadIdx.y].x / l1;
	normDir[threadIdx.z][threadIdx.y].y = sdirCent[threadIdx.z][threadIdx.y].y / l1;

	curang = angs[angIdx];
	cosT[threadIdx.z] = cosTs[angIdx];
	sinT[threadIdx.z] = sinTs[angIdx];
	__syncthreads();
	if (detIdU < DNU && detIdV < DNV && angIdx < PN)
	{
		float zP = zPos[angIdx];
		float3 curSour = make_float3(
			x0 * cosT[threadIdx.z] - y0 * sinT[threadIdx.z],
			x0 * sinT[threadIdx.z] + y0 * cosT[threadIdx.z],
			z0 + zP);
		float initDetX = d_xds[detIdU];
		float initDetY = d_yds[detIdU];
		float3 curDet = make_float3(
			initDetX * cosT[threadIdx.z] - initDetY * sinT[threadIdx.z],
			initDetX * sinT[threadIdx.z] + initDetY * cosT[threadIdx.z],
			d_zds[detIdV] + zP);

		float3 dir = curDet - curSour;
		float summ = 0;
		float obj = 0;
		float real, realZ;
		float invdz = 1.0 / dz;
		float invdx = 1.0 / dx;

		float constVal(1.0);
		float fact1(1.0f), fact2(1.0f);
		if ((curang > PI_4 && curang <= PI_3_4) || (curang > PI_5_4 && curang <= PI_7_4))
		{
			summ = 0;
			fact1 = dir.y / dir.x;
			fact2 = dir.z / dir.x;
			constVal = dx / (abs(normDir[threadIdx.z][threadIdx.y].x) * cosGamma[threadIdx.x]);
			obj = (-objCntIdxX) * dx;
			real = fmaf(obj - curDet.x, fact1, curDet.y);
			realZ = fmaf(obj - curDet.x, fact2, curDet.z);

			fact1 = dx * fact1;
			fact2 = dx * fact2;
			for (int ii = 0; ii < XN; ++ii)
			{
				summ += tex3D<float>(volTex, fmaf(realZ, invdz, objCntIdxZ + 0.5), ii + 0.5, fmaf(real, invdx, objCntIdxY + 0.5));
				obj += dx;
				real += fact1;
				realZ += fact2;

			}
			__syncthreads();
			proj[(angIdx * DNU + detIdU) * DNV + detIdV] = summ * constVal;
		}
		else
		{
			summ = 0;
			fact1 = dir.x / dir.y;
			fact2 = dir.z / dir.y;
			constVal = dx / (abs(normDir[threadIdx.z][threadIdx.y].y) * cosGamma[threadIdx.x]);
			obj = (-objCntIdxX) * dx;
			real = fmaf(obj - curDet.y, fact1, curDet.x);
			realZ = fmaf(obj - curDet.y, fact2, curDet.z);

			fact1 = dx * fact1;
			fact2 = dx * fact2;
			for (int ii = 0; ii < YN; ++ii)
			{
				summ += tex3D<float>(volTex, fmaf(realZ, invdz, objCntIdxZ + 0.5), fmaf(real, invdx, objCntIdxX + 0.5), ii + 0.5);
				obj += dx;
				real += fact1;
				realZ += fact2;

			}
			__syncthreads();
			proj[(angIdx * DNU + detIdU) * DNV + detIdV] = summ * constVal;
		}

	}

}

void DD3_gpu_proj_volumerendering(
	float x0, float y0, float z0,
	int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN,
	float* vol, float* hprj, float dx, float dz,
	byte* mask, int gpunum)
{
	for (int ii = 0; ii != XN * YN; ++ii)
	{
		byte v = mask[ii];
		for (int jj = 0; jj != ZN; ++jj)
		{
			vol[ii * ZN + jj] = vol[ii * ZN + jj] * v;
		}
	}

	float* bxds = new float[DNU + 1];
	float* byds = new float[DNU + 1];
	float* bzds = new float[DNV + 1];
	DD3Boundaries(DNU + 1, xds, bxds);
	DD3Boundaries(DNU + 1, yds, byds);
	DD3Boundaries(DNV + 1, zds, bzds);

	cudaSetDevice(gpunum);
	cudaDeviceReset();
	cudaStream_t streams[4];
	cudaStreamCreate(&streams[0]);
	cudaStreamCreate(&streams[1]);
	cudaStreamCreate(&streams[2]);
	cudaStreamCreate(&streams[3]);

	float objCntIdxX = (XN - 1.0) * 0.5 - imgXCenter / dx;
	float objCntIdxY = (YN - 1.0) * 0.5 - imgYCenter / dx;
	float objCntIdxZ = (ZN - 1.0) * 0.5 - imgZCenter / dz;

	cudaExtent volumeSize;
	volumeSize.width = ZN;
	volumeSize.height = XN;
	volumeSize.depth = YN;

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	cudaArray* d_volumeArray;
	cudaMalloc3DArray(&d_volumeArray, &channelDesc, volumeSize);
	cudaMemcpy3DParms copyParams = { 0 };
	copyParams.srcPtr = make_cudaPitchedPtr((void*)vol, volumeSize.width * sizeof(float), volumeSize.width, volumeSize.height);
	copyParams.dstArray = d_volumeArray;
	copyParams.extent = volumeSize;
	copyParams.kind = cudaMemcpyHostToDevice;
	cudaMemcpy3D(&copyParams);

	cudaResourceDesc resDesc;
	memset(&resDesc, 0, sizeof(resDesc));
	resDesc.resType = cudaResourceTypeArray;
	resDesc.res.array.array = d_volumeArray;

	cudaTextureDesc texDesc;
	memset(&texDesc, 0, sizeof(texDesc));
	texDesc.addressMode[0] = cudaAddressModeBorder;
	texDesc.addressMode[1] = cudaAddressModeBorder;
	texDesc.addressMode[2] = cudaAddressModeBorder;
	texDesc.filterMode = cudaFilterModeLinear;
	texDesc.readMode = cudaReadModeElementType;
	texDesc.normalizedCoords = false;

	cudaTextureObject_t texObj = 0;
	cudaCreateTextureObject(&texObj, &resDesc, &texDesc, nullptr);

	thrust::device_vector<float> prj(DNU * DNV * PN, 0);
	thrust::device_vector<float> angs(hangs, hangs + PN);
	thrust::device_vector<float> zPos(hzPos, hzPos + PN);
	thrust::device_vector<float> d_xds(xds, xds + DNU);
	thrust::device_vector<float> d_yds(yds, yds + DNU);
	thrust::device_vector<float> d_zds(zds, zds + DNV);
	thrust::device_vector<float3> cursour(PN);
	thrust::device_vector<float2> dirCent(PN * DNU);
	thrust::device_vector<float> cosTs(PN);
	thrust::device_vector<float> sinTs(PN);
	thrust::device_vector<float> cosGamma(DNV);

	dim3 blkc(64, 16, 1);
	dim3 gidc(
		(DNV + blkc.x - 1) / blkc.x,
		(DNU + blkc.y - 1) / blkc.y,
		(PN + blkc.z - 1) / blkc.z);
	dim3 blkc2(64, 16);
	dim3 gidc2(
		(DNU + blkc2.x - 1) / blkc2.x,
		(PN + blkc2.y - 1) / blkc2.y);
	//Calculate the constant values;

	dim3 blk(BLKX, BLKY, BLKZ);
	dim3 gid(
		(DNV + blk.x - 1) / blk.x,
		(DNU + blk.y - 1) / blk.y,
		(PN + blk.z - 1) / blk.z);

	DD3_gpu_proj_volumerendering_ker << <gid, blk >> >
		(texObj,
		thrust::raw_pointer_cast(&prj[0]),
		x0, y0, z0,
		thrust::raw_pointer_cast(&d_xds[0]),
		thrust::raw_pointer_cast(&d_yds[0]),
		thrust::raw_pointer_cast(&d_zds[0]),
		thrust::raw_pointer_cast(&dirCent[0]),
		thrust::raw_pointer_cast(&cosTs[0]),
		thrust::raw_pointer_cast(&sinTs[0]),
		objCntIdxX, objCntIdxY, objCntIdxZ,
		dx, dz, XN, YN, ZN, DNU, DNV, PN,
		thrust::raw_pointer_cast(&angs[0]),
		thrust::raw_pointer_cast(&zPos[0]),
		thrust::raw_pointer_cast(&cursour[0]),
		thrust::raw_pointer_cast(&cosGamma[0]));

	thrust::copy(prj.begin(), prj.end(), hprj);

	cudaDestroyTextureObject(texObj);
	cudaFreeArray(d_volumeArray);
	prj.clear();
	angs.clear();
	cursour.clear();
	cosTs.clear();
	sinTs.clear();
	dirCent.clear();
	cosGamma.clear();
	d_xds.clear();
	d_yds.clear();
	d_zds.clear();

	delete[] bxds;
	delete[] byds;
	delete[] bzds;
}





__global__ void DD3_gpu_proj_pseudodistancedriven_ker(
	cudaTextureObject_t volTex,
	float* proj, float3 s,
	float* d_xds, float* d_yds, float* d_zds,
	float3* cossinT,
	float3 objCntIdx,
	float dx, float dz,
	int XN, int YN,
	int DNU, int DNV, int PN)
{
	int detIdV = threadIdx.x + blockIdx.x * blockDim.x;
	int detIdU = threadIdx.y + blockIdx.y * blockDim.y;
	int angIdx = threadIdx.z + blockIdx.z * blockDim.z;
	if (detIdV < DNV && detIdU < DNU && angIdx < PN)
	{
		float3 cossin = cossinT[angIdx];
		float3 cursour = make_float3(
			s.x * cossin.x - s.y * cossin.y,
			s.x * cossin.y + s.y * cossin.x,
			s.z + cossin.z);

		float summ = d_xds[detIdU];
		float obj = d_yds[detIdU];
		float idx = d_zds[detIdV];
		float3 curDet = make_float3(
			summ * cossin.x - obj * cossin.y,
			summ * cossin.y + obj * cossin.x,
			idx + cossin.z);

		float3 dir = normalize(curDet - cursour);
		summ = 0;
		obj = 0;
		float idxZ;
		if (fabsf(cossin.x) <= fabsf(cossin.y))
		{
			summ = 0;
			for (int ii = 0; ii < XN; ++ii)
			{
				obj = (ii - objCntIdx.x) * dx;
				idx = (obj - curDet.x) / dir.x * dir.y + curDet.y;
				idxZ = (obj - curDet.x) / dir.x * dir.z + curDet.z;

				idx = idx / dx + objCntIdx.y + 0.5;
				idxZ = idxZ / dz + objCntIdx.z + 0.5;
				summ += tex3D<float>(volTex, idxZ, ii + 0.5f, idx);
			}
			__syncthreads();
			proj[(angIdx * DNU + detIdU) * DNV + detIdV] = summ * dx / fabsf(dir.x);
		}
		else
		{
			summ = 0;
			for (int jj = 0; jj != YN; ++jj)
			{
				obj = (jj - objCntIdx.y) * dx;
				idx = (obj - curDet.y) / dir.y * dir.x + curDet.x;
				idxZ = (obj - curDet.y) / dir.y * dir.z + curDet.z;

				idx = idx / dx + objCntIdx.x + 0.5;
				idxZ = idxZ / dz + objCntIdx.z + 0.5;
				summ += tex3D<float>(volTex, idxZ, idx, jj + 0.5f);
			}
			__syncthreads();
			proj[(angIdx * DNU + detIdU) * DNV + detIdV] = summ * dx / fabsf(dir.y);
		}
	}
}

void DD3_gpu_proj_pseudodistancedriven(
	float x0, float y0, float z0,
	int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN, float* hvol, float* hprj,
	float dx, float dz, byte* mask, int gpunum)
{
	for (int ii = 0; ii != XN * YN; ++ii)
	{
		byte v = mask[ii];
		for (int jj = 0; jj != ZN; ++jj)
		{
			hvol[ii * ZN + jj] = hvol[ii * ZN + jj] * v;
		}
	}


	float* bxds = new float[DNU + 1];
	float* byds = new float[DNU + 1];
	float* bzds = new float[DNV + 1];
	DD3Boundaries(DNU + 1, xds, bxds);
	DD3Boundaries(DNU + 1, yds, byds);
	DD3Boundaries(DNV + 1, zds, bzds);

	cudaSetDevice(gpunum);
	cudaDeviceReset();

	const int TOTVN = XN * YN * ZN;
	float objCntIdxX = (XN - 1.0) * 0.5 - imgXCenter / dx;
	float objCntIdxY = (YN - 1.0) * 0.5 - imgYCenter / dx;
	float objCntIdxZ = (ZN - 1.0) * 0.5 - imgZCenter / dz;

	thrust::device_vector<float> vol(hvol, hvol + TOTVN);


	cudaExtent volumeSize;
	volumeSize.width = ZN;
	volumeSize.height = XN;
	volumeSize.depth = YN;

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

	cudaArray* d_volumeArray;
	cudaMalloc3DArray(&d_volumeArray, &channelDesc, volumeSize);
	cudaMemcpy3DParms copyParams = { 0 };
	copyParams.srcPtr = make_cudaPitchedPtr((void*)thrust::raw_pointer_cast(&vol[0]),
		volumeSize.width * sizeof(float),
		volumeSize.width, volumeSize.height);
	copyParams.dstArray = d_volumeArray;
	copyParams.extent = volumeSize;
	copyParams.kind = cudaMemcpyDeviceToDevice;
	cudaMemcpy3D(&copyParams);
	vol.clear();

	cudaResourceDesc resDesc;
	memset(&resDesc, 0, sizeof(resDesc));
	resDesc.resType = cudaResourceTypeArray;
	resDesc.res.array.array = d_volumeArray;
	cudaTextureDesc texDesc;
	memset(&texDesc, 0, sizeof(texDesc));
	texDesc.addressMode[0] = cudaAddressModeBorder;
	texDesc.addressMode[1] = cudaAddressModeBorder;
	texDesc.addressMode[2] = cudaAddressModeBorder;
	texDesc.filterMode = cudaFilterModeLinear;
	texDesc.readMode = cudaReadModeElementType;
	texDesc.normalizedCoords = false;

	cudaTextureObject_t texObj = 0;
	cudaCreateTextureObject(&texObj, &resDesc, &texDesc, nullptr);

	thrust::device_vector<float> prj(DNU * DNV * PN, 0);
	thrust::device_vector<float> angs(hangs, hangs + PN);
	thrust::device_vector<float> zPos(hzPos, hzPos + PN);
	thrust::device_vector<float> d_xds(xds, xds + DNU);
	thrust::device_vector<float> d_yds(yds, yds + DNU);
	thrust::device_vector<float> d_zds(zds, zds + DNV);

	thrust::device_vector<float3> cossinZT(PN);
	thrust::transform(
		thrust::make_zip_iterator(thrust::make_tuple(angs.begin(), zPos.begin())),
		thrust::make_zip_iterator(thrust::make_tuple(angs.end(), zPos.end())),
		cossinZT.begin(), CTMBIR::ConstantForBackProjection4(x0, y0, z0));


	dim3 blk(64, 16, 1);
	dim3 gid(
		(DNV + blk.x - 1) / blk.x,
		(DNU + blk.y - 1) / blk.y,
		(PN + blk.z - 1) / blk.z);

	DD3_gpu_proj_pseudodistancedriven_ker << <gid, blk >> >(
		texObj, thrust::raw_pointer_cast(&prj[0]),
		make_float3(x0, y0, z0),
		thrust::raw_pointer_cast(&d_xds[0]),
		thrust::raw_pointer_cast(&d_yds[0]),
		thrust::raw_pointer_cast(&d_zds[0]),
		thrust::raw_pointer_cast(&cossinZT[0]),
		make_float3(objCntIdxX, objCntIdxY, objCntIdxZ),
		dx, dz, XN, YN, DNU, DNV, PN);
	thrust::copy(prj.begin(), prj.end(), hprj);

	delete[] bxds;
	delete[] byds;
	delete[] bzds;


}




extern "C"
void DD3Proj_gpu_alreadyinGPU(
float x0, float y0, float z0,
int DNU, int DNV,
const thrust::device_vector<float>& xds,
const thrust::device_vector<float>& yds,
const thrust::device_vector<float>& zds,
float imgXCenter, float imgYCenter, float imgZCenter,
const thrust::device_vector<float>& hangs,
const thrust::device_vector<float>& hzPos, int PN,
int XN, int YN, int ZN,
const thrust::device_vector<float>& hvol,
thrust::device_vector<float>& hprj,
float dx, float dz,
const thrust::device_vector<byte>& mask, int gpunum, int prjMode)
{
	DD3_gpu_proj_branchless_sat2d_alreadyinGPU(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
					hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz, mask, gpunum);
}









__global__ void DDM3D_EA_helical_proj_GPU(float* proj, const float* vol,
	const float S2O, const float O2D, const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeV, const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN,
	const int DNU, const int DNV, const int PN,
	const float dbeta, const float ddv, const float dx, const float dy, const float dz,
	const float* angs,	const float* zShifts)
{
	const int detIdV = threadIdx.x + blockIdx.x * blockDim.x; //(Z,X,Y)
	const int detIdU = threadIdx.y + blockIdx.y * blockDim.y;
	const int angIdx = threadIdx.z + blockIdx.z * blockDim.z;
	if (detIdU < DNU && detIdV < DNV && angIdx < PN)
	{
		float curang = angs[angIdx];
		float cosT = cos(curang);
		float sinT = sin(curang);

		float cursourx = S2O * cosT;
		float cursoury = S2O * sinT;
		float cursourz = zShifts[angIdx]; //ÕâµãÓëcone beam ²»Í¬
		float summ = 0;
		float beta = (detIdU - detCntIdU) * dbeta;

		float cosBeta = cos(beta);
		float sinBeta = sin(beta);
		float initDetX = (-O2D) * cosBeta - 0 * sinBeta + S2O;
		float initDetY = (-O2D) * sinBeta + 0 * cosBeta;
		float initDetZ = (detIdV - detCntIdV) * ddv;

		beta = (detIdU - detCntIdU - 0.5) * dbeta;
		cosBeta = cos(beta);
		sinBeta = sin(beta);
		float initDetLX = (-O2D) * cosBeta - 0 * sinBeta + S2O;
		float initDetLY = (-O2D) * sinBeta + 0 * cosBeta;


		beta = (detIdU - detCntIdU + 0.5) * dbeta;
		cosBeta = cos(beta);
		sinBeta = sin(beta);
		float initDetRX = (-O2D) * cosBeta - 0 * sinBeta + S2O;
		float initDetRY = (-O2D) * sinBeta + 0 * cosBeta;

		float initDetDZ = (detIdV - detCntIdV - 0.5) * ddv;
		float initDetUZ = (detIdV - detCntIdV + 0.5) * ddv;

		float curDetLX = initDetLX * cosT - initDetLY * sinT;
		float curDetLY = initDetLX * sinT + initDetLY * cosT;

		float curDetRX = initDetRX * cosT - initDetRY * sinT;
		float curDetRY = initDetRX * sinT + initDetRY * cosT;

		float curDetDX = initDetX * cosT - initDetY * sinT;
		float curDetDY = initDetX * sinT + initDetY * cosT;
		float curDetDZ = initDetDZ + zShifts[angIdx]; //

		float curDetUX = curDetDX;
		float curDetUY = curDetDY;
		float curDetUZ = initDetUZ + zShifts[angIdx]; //

		float curDetX = initDetX * cosT - initDetY * sinT;
		float curDetY = initDetX * sinT + initDetY * cosT;

		float dirX = curDetX - cursourx;
		float dirY = curDetY - cursoury;

		if ((curang > PI * 0.25 && curang <= PI * 0.75) || (curang >= PI * 1.25 && curang < PI * 1.75))
		{
			float cosAlpha = abs(dirY / sqrt(dirY * dirY + dirX * dirX));
			float cosGamma = abs(sqrt((S2O + O2D) * (S2O + O2D) - initDetZ * initDetZ) / (S2O + O2D));


			float detPosLX = -cursoury * (curDetLX - cursourx) / (curDetLY - cursoury) + cursourx; //×ó±ßµãÔÚXOZÆœÃæÉÏµÄÍ¶Ó°;
			float detPosRX = -cursoury * (curDetRX - cursourx) / (curDetRY - cursoury) + cursourx;
			float detPosDZ = -cursoury * (curDetDZ - cursourz) / (curDetDY - cursoury) + cursourz;
			float detPosUZ = -cursoury * (curDetUZ - cursourz) / (curDetUY - cursoury) + cursourz;

			float detprojLength = abs(detPosLX - detPosRX);
			float detprojHeight = abs(detPosUZ - detPosDZ);

			//ŒÙÉè×ó±ßµÄÐ¡;
			if (detPosLX > detPosRX)
			{
				float tt = detPosLX;
				detPosLX = detPosRX;
				detPosRX = tt;
			}
			//ŒÙÉèÏÂ±ßµÄÐ¡;
			if (detPosDZ > detPosUZ)
			{
				float tt = detPosDZ;
				detPosDZ = detPosUZ;
				detPosUZ = tt;
			}


			for (size_t jj = 0; jj < YN; jj++)
			{
				float objY = (jj - YN / 2.0 + 0.5) * dy;
				float temp = (objY - cursoury) / (curDetLY - cursoury);

				float minX = temp * (curDetLX - cursourx) + cursourx;
				float maxX = temp * (curDetRX - cursourx) + cursourx;
				float minZ = temp * (curDetDZ - cursourz) + cursourz;
				float maxZ = temp * (curDetUZ - cursourz) + cursourz;
				if (minX > maxX)
				{
					float tt = minX;
					minX = maxX;
					maxX = tt;
					//std::swap(minX, maxX);
				}
				if (minZ > maxZ)
				{
					float tt = minZ;
					minZ = maxZ;
					maxZ = tt;
					//std::swap(minZ, maxZ);
				}
				int minXIdx = floor(minX / dx + XN / 2.0) - 2;
				int maxXIdx = ceil(maxX / dx + XN / 2.0) + 2;
				int minZIdx = floor(minZ / dz + ZN / 2.0) - 2;
				int maxZIdx = ceil(maxZ / dz + ZN / 2.0) + 2;
				if (maxXIdx < 0){ continue; }
				if (minXIdx > XN){ continue; }
				if (maxZIdx < 0){ continue; }
				if (minZIdx > ZN){ continue; }
				if (minXIdx < 0){ minXIdx = 0; }
				if (maxXIdx > XN){ maxXIdx = XN; }
				if (minZIdx < 0){ minZIdx = 0; }
				if (maxZIdx > ZN){ maxZIdx = ZN; }


				for (size_t ii = minXIdx; ii < maxXIdx; ii++)
				{
					float curminx = (cursourx - (ii - XN / 2.0) * dx) * cursoury / (objY - cursoury) + cursourx;
					float curmaxx = (cursourx - ((ii + 1) - XN / 2.0) * dx) * cursoury / (objY - cursoury) + cursourx;
					float intersectL = intersectLength<float>(detPosLX, detPosRX, curminx, curmaxx);
					if (intersectL > 0)
					{
						for (size_t kk = minZIdx; kk < maxZIdx; kk++)
						{

							float curminz = (cursourz - (kk - ZN / 2.0) * dz) * cursoury / (objY - cursoury) + cursourz;
							float curmaxz = (cursourz - ((kk + 1) - ZN / 2.0) * dz) * cursoury / (objY - cursoury) + cursourz;

							float intersectH = intersectLength<float>(detPosDZ, detPosUZ, curminz, curmaxz);
							if (intersectH > 0)
							{
								summ += vol[(jj * XN + ii) * ZN + kk] * (intersectL * intersectH) / (detprojLength * detprojHeight * cosAlpha * cosGamma) * dx;
							}
						}
					}
					else
					{
						continue;
					}
				}
			}
			proj[(angIdx * DNU + detIdU) * DNV + detIdV] = summ;
		}
		else
		{
			float cosAlpha = abs(dirX / sqrt(dirY * dirY + dirX * dirX));
			float cosGamma = abs(sqrt((S2O + O2D) * (S2O + O2D) - initDetZ * initDetZ) / (S2O + O2D));


			float detPosLY = -cursourx * (curDetLY - cursoury) / (curDetLX - cursourx) + cursoury; //×ó±ßµãÔÚXOZÆœÃæÉÏµÄÍ¶Ó°;
			float detPosRY = -cursourx * (curDetRY - cursoury) / (curDetRX - cursourx) + cursoury;
			float detPosDZ = -cursourx * (curDetDZ - cursourz) / (curDetDX - cursourx) + cursourz;
			float detPosUZ = -cursourx * (curDetUZ - cursourz) / (curDetUX - cursourx) + cursourz;

			float detprojLength = abs(detPosLY - detPosRY);
			float detprojHeight = abs(detPosUZ - detPosDZ);

			//
			if (detPosLY > detPosRY)
			{
				float tt = detPosLY;
				detPosLY = detPosRY;
				detPosRY = tt;
				//std::swap(detPosLY, detPosRY);
			}
			//
			if (detPosDZ > detPosUZ)
			{
				float tt = detPosDZ;
				detPosDZ = detPosUZ;
				detPosUZ = tt;
				//std::swap(detPosDZ, detPosUZ);
			}

			for (size_t ii = 0; ii < XN; ii++)
			{
				float objX = (ii - XN / 2.0 + 0.5) * dx;
				float temp = (objX - cursourx) / (curDetLX - cursourx);

				float minY = temp * (curDetLY - cursoury) + cursoury;
				float maxY = temp * (curDetRY - cursoury) + cursoury;
				float minZ = temp * (curDetDZ - cursourz) + cursourz;
				float maxZ = temp * (curDetUZ - cursourz) + cursourz;
				if (minY > maxY)
				{
					float tt = minY;
					minY = maxY;
					maxY = tt;
					//std::swap(minY, maxY);
				}
				if (minZ > maxZ)
				{
					float tt = minZ;
					minZ = maxZ;
					maxZ = tt;
					//std::swap(minZ, maxZ);
				}
				int minYIdx = floor(minY / dy + YN / 2.0) - 2;
				int maxYIdx = ceil(maxY / dy + YN / 2.0) + 2;
				int minZIdx = floor(minZ / dz + ZN / 2.0) - 2;
				int maxZIdx = ceil(maxZ / dz + ZN / 2.0) + 2;
				if (maxYIdx < 0){ continue; }
				if (minYIdx > XN){ continue; }
				if (maxZIdx < 0){ continue; }
				if (minZIdx > ZN){ continue; }
				if (minYIdx < 0){ minYIdx = 0; }
				if (maxYIdx > XN){ maxYIdx = YN; }
				if (minZIdx < 0){ minZIdx = 0; }
				if (maxZIdx > ZN){ maxZIdx = ZN; }


				for (size_t jj = minYIdx; jj < maxYIdx; jj++)
				{
					float curminy = (cursoury - (jj - YN / 2.0) * dy) * cursourx / (objX - cursourx) + cursoury;
					float curmaxy = (cursoury - ((jj + 1) - YN / 2.0) * dy) * cursourx / (objX - cursourx) + cursoury;
					float intersectL = intersectLength<float>(detPosLY, detPosRY, curminy, curmaxy);
					if (intersectL > 0)
					{
						for (size_t kk = minZIdx; kk < maxZIdx; kk++)
						{

							float curminz = (cursourz - (kk - ZN / 2.0) * dz) * cursourx / (objX - cursourx) + cursourz;
							float curmaxz = (cursourz - ((kk + 1) - ZN / 2.0) * dz) * cursourx / (objX - cursourx) + cursourz;

							float intersectH = intersectLength<float>(detPosDZ, detPosUZ, curminz, curmaxz);
							if (intersectH > 0)
							{
								summ += vol[(jj * XN + ii) * ZN + kk] * (intersectL * intersectH) / (detprojLength * detprojHeight * cosAlpha * cosGamma) * dx;
							}
						}
					}
					else
					{
						continue;
					}

				}

			}
			proj[(angIdx * DNU + detIdU) * DNV + detIdV] = summ;

		}

	}
}


void DD3_gpu_proj_bruteForce(
		float x0, float y0, float z0,
		int DNU, int DNV, float* xds, float* yds, float* zds,
		float imgXCenter, float imgYCenter, float imgZCenter,
		float* hangs, float* hzPos,
		int PN, int XN, int YN, int ZN,
		float* hvol,float* hprj,float dx,float dz,byte* mask,int gpunum)
{
//	dim3 blk(BLKX,BLKY,BLKZ);
//	dim3 gid(
//			(DNV + blk.x - 1) / blk.x,
//			(DNU + blk.y - 1) / blk.y,
//			(PN + blk.z - 1) / blk.z);
//
//	thrust::device_vector<float> prj(hprj, hprj + DNU * DNV * PN);
//	thrust::device_vector<float> vol(hvol, hvol + XN * YN * ZN);

//	thrust::device_vector<float> angs(hangs,hangs + PN);
//	thrust::device_vector<float> zShifts(hzPos, hzPos + PN);

//	const float S2O = sqrtf(x0 * x0 + y0 * y0);
//	const float S2D = sqrtf(powf(x0 - xds[0],2) + powf(y0 - yds[0],2));
//	const float O2D = S2D - S2O;

//	const float objSizeX = dx * XN;
//	const float objSizeY = objSizeX;
//	const float objSizeZ = dz * ZN;
//	const float detUSize = sqrtf(powf(xds[1] - xds[0],2) + powf(yds[1] - yds[0],2));
//	const float ddv = zds[1] - zds[0];
//	const float dbeta = atan(detUSize * 0.5 / S2D) * 2.0;
//	const float detArc = dbeta * DNU;
//	const float detSizeV = ddv * DNV;

//	DDM3D_EA_helical_proj_GPU<<<gid,blk>>>(thrust::raw_pointer_cast(&prj[0]), thrust::raw_pointer_cast(&vol[0]),
//		S2O, O2D, objSizeX, objSizeY, objSizeZ, detArc, detSizeV, const float detCntIdU, const float detCntIdV,
//		XN, YN, ZN, DNU, DNV, PN, dbeta, ddv, dx, dx, dz, thrust::raw_pointer_cast(&angs[0]), thrust::raw_pointer_cast(&zShifts[0]));
}





extern "C"
void DD3Proj_gpu(
	float x0, float y0, float z0,
	int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN,
	float* hvol, float* hprj,
	float dx, float dz,
	byte* mask, int gpunum, int prjMode)
{
	switch (prjMode)
	{
	case 0:
		DD3_gpu_proj_branchless_sat2d(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
			hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz, mask, gpunum);
		break;
	case 1:
		DD3_gpu_proj_volumerendering(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
			hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz, mask, gpunum);
		break;
	case 2:
	
	case 3:
		DD3_gpu_proj_pseudodistancedriven(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
			hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz, mask, gpunum);
		break;
	case 4: //Siddon's algorithm exactly which is slow because of the L2 cache in device function
		DD3ProjSiddon_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
				hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz, mask, gpunum);
		break;
	case 5: //Brute Force projection (Not finished)
		DD3_gpu_proj_bruteForce(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
				hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz, mask, gpunum);
		break;
	default:
		DD3_gpu_proj_branchless_sat2d(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
			hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz, mask, gpunum);
		break;
	}
}




//Use the split-collect method to do the projection
extern "C"
void DD3ProjHelical_3GPU(
	float x0, float y0, float z0,
	int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN,
	float* hvol, float* hprj,
	float dx, float dz,
	byte* mask, int methodId, int (&startPN)[3])
{
	thrust::host_vector<float> h_angs(hangs,hangs+PN);
	thrust::host_vector<float> h_zPos(hzPos,hzPos+PN);
	//Define the projection Volume Range
	int ObjIdx_Start[3];
	int ObjIdx_End[3];

	int PrjIdx_Start[3] = {startPN[0], startPN[1], startPN[2]}; // The start index of the projection angle
	int PrjIdx_End[3] = {startPN[1], startPN[2], PN}; // The end index of the projection angle (not included)
	int SPN[3] = {startPN[1] - startPN[0], startPN[2] - startPN[1], PN - startPN[2]}; // how many views are contained in each sub projection problems
	int prefixSPN[3] = {0, SPN[0], SPN[0] + SPN[1]}; // Prefix of the number of projection views

	int SZN[3] = {0, 0, 0}; // The slices number of each sub volume
	float detStpZ = zds[1] - zds[0]; // detector cell height
	float detCntIdxV = -zds[0] / detStpZ; // Detector center along Z direction
	float objCntIdxZ = (ZN - 1.0) / 2.0 - imgZCenter / dz; //object center in Z direction

	float** subVol = new float*[3];
	float subImgZCenter[3] = {0, 0, 0};

	omp_set_num_threads(3);
#pragma omp parallel for
	for(int i = 0; i < 3; ++i)  //The last one has problem!!!!!!!!!!
	{
		getVolZIdxPair<float>(h_zPos,PrjIdx_Start[i],PrjIdx_End[i],
					detCntIdxV, detStpZ, DNV, objCntIdxZ,	dz, ZN, ObjIdx_Start[i], ObjIdx_End[i]);
		std::cout<<i<<" "<<ObjIdx_Start[i]<<" "<<ObjIdx_End[i]<<"\n";
		SZN[i] = ObjIdx_End[i] - ObjIdx_Start[i];
		subVol[i] = new float[XN * YN * SZN[i]];
		//Divide the volume
		getSubVolume<float>(hvol, XN * YN, ZN, ObjIdx_Start[i], ObjIdx_End[i], subVol[i]);
		//Calculate the corresponding center position
		subImgZCenter[i] = ((ObjIdx_End[i] + ObjIdx_Start[i] - (ZN - 1.0)) * dz + imgZCenter * 2.0) / 2.0;

		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
				imgXCenter, imgYCenter, subImgZCenter[i], hangs + prefixSPN[i] , hzPos + prefixSPN[i], SPN[i],
				XN, YN, SZN[i], subVol[i], hprj + DNU * DNV * prefixSPN[i],dx,dz,mask,i, methodId);
	}

	delete[] subVol[0];
	delete[] subVol[1];
	delete[] subVol[2];
	delete[] subVol;
}



extern "C"
void DD3ProjHelical_4GPU(
	float x0, float y0, float z0,
	int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN,
	float* hvol, float* hprj,
	float dx, float dz,
	byte* mask, int methodId, int (&startPN)[4])
{
	thrust::host_vector<float> h_angs(hangs,hangs+PN);
	thrust::host_vector<float> h_zPos(hzPos,hzPos+PN);
	//Define the projection Volume Range
	int ObjIdx_Start[4];
	int ObjIdx_End[4];

	int PrjIdx_Start[4] = {startPN[0], startPN[1], startPN[2], startPN[3]}; // The start index of the projection angle
	int PrjIdx_End[4] = {startPN[1], startPN[2],startPN[3], PN}; // The end index of the projection angle (not included)
	int SPN[4] = {startPN[1] - startPN[0], startPN[2] - startPN[1], startPN[3] - startPN[2], PN - startPN[3]}; // how many views are contained in each sub projection problems
	int prefixSPN[4] = {0, SPN[0], SPN[0] + SPN[1], SPN[0] + SPN[1] + SPN[2]}; // Prefix of the number of projection views

	int SZN[4] = {0, 0, 0, 0}; // The slices number of each sub volume
	float detStpZ = zds[1] - zds[0]; // detector cell height
	float detCntIdxV = -zds[0] / detStpZ; // Detector center along Z direction
	float objCntIdxZ = (ZN - 1.0) / 2.0 - imgZCenter / dz; //object center in Z direction

	float** subVol = new float*[4];
	float subImgZCenter[4] = {0, 0, 0, 0};

	omp_set_num_threads(4);
#pragma omp parallel for
	for(int i = 0; i < 4; ++i)  //The last one has problem!!!!!!!!!!
	{
		getVolZIdxPair<float>(h_zPos,PrjIdx_Start[i],PrjIdx_End[i],
					detCntIdxV, detStpZ, DNV, objCntIdxZ,	dz, ZN, ObjIdx_Start[i], ObjIdx_End[i]);
		std::cout<<i<<" "<<ObjIdx_Start[i]<<" "<<ObjIdx_End[i]<<"\n";
		SZN[i] = ObjIdx_End[i] - ObjIdx_Start[i];
		subVol[i] = new float[XN * YN * SZN[i]];
		//Divide the volume
		getSubVolume<float>(hvol, XN * YN, ZN, ObjIdx_Start[i], ObjIdx_End[i], subVol[i]);
		//Calculate the corresponding center position
		subImgZCenter[i] = ((ObjIdx_End[i] + ObjIdx_Start[i] - (ZN - 1.0)) * dz + imgZCenter * 2.0) / 2.0;

		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
				imgXCenter, imgYCenter, subImgZCenter[i], hangs + prefixSPN[i] , hzPos + prefixSPN[i], SPN[i],
				XN, YN, SZN[i], subVol[i], hprj + DNU * DNV * prefixSPN[i],dx, dz, mask, i, methodId);
	}

	delete[] subVol[0];
	delete[] subVol[1];
	delete[] subVol[2];
	delete[] subVol[3];
	delete[] subVol;
}


