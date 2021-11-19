#include <algorithm>
#include "utilities.hpp"
#include "utilities.cuh"
#include "DD3_GPU_Proj.h"
#include "CTLAB.h"
#include "Geometry.h"
#include "TextureObjectProvider.h"
#include "GenerateSummedAreaTable.h"
#include "cudaCheckReturner.h"
#include "VolumeToTextures.h"
#define BLKX 32
#define BLKY 8
#define BLKZ 1

// masking the volume with provided mask. The mask size is XN * YN
// the volume is with size [ZN, XN, YN] ZN is the primary order.
template<typename T>
void maskingVolume(T* vol, byte* mask, int XN, int YN, int ZN) {
	//Pre compute mask.*vol;
	for (int ii = 0; ii != XN * YN; ++ii) {
		byte v = mask[ii];
		std::for_each(vol + ii * ZN, vol + (ii + 1) * ZN, [v](T& val) { val *= v; });
	}
}

template<typename Ta, typename Tb>
__global__ void naive_copyToTwoVolumes(Ta* in_ZXY,
	Tb* out_ZXY, Tb* out_ZYX,
	int XN, int YN, int ZN) {
	const int idz = threadIdx.x + blockIdx.x * blockDim.x; // Z dimension major
	const int idx = threadIdx.y + blockIdx.y * blockDim.y;
	const int idy = threadIdx.z + blockIdx.z * blockDim.z;
	if (idx < XN && idy < YN && idz < ZN) {
		int i = (idy * XN + idx) * ZN + idz;
		int ni = (idy * (XN + 1) + (idx + 1)) * (ZN + 1) + idz + 1;
		int nj = (idx * (YN + 1) + (idy + 1)) * (ZN + 1) + idz + 1;

		out_ZXY[ni] = in_ZXY[i];
		out_ZYX[nj] = in_ZXY[i];
	}
}

template<typename Ta, typename Tb>
__global__ void naive_herizontalIntegral(Ta* in, Tb* out, int N, int ZN) {
	const int zi = threadIdx.x + blockIdx.x * blockDim.x;
	if (zi < ZN) {
		out[zi] = in[zi];
		for (int i = 1; i < N; ++i) {
			out[i * ZN + zi] = out[(i - 1) * ZN + zi]
				+ in[i * ZN + zi];
		}
	}
}

template<typename Ta, typename Tb>
__global__ void naive_verticalIntegral(Ta* in, Tb* out, int N, int ZN) {
	int xyi = threadIdx.x + blockIdx.x * blockDim.x;
	if (xyi < N) {
		out[xyi * ZN] = in[xyi * ZN];
		for (int ii = 1; ii < ZN; ++ii)	{
			out[xyi * ZN + ii] = out[xyi * ZN + ii - 1]
				+ in[xyi * ZN + ii];
		}
	}
}

// template specialization for double precision into int2 datatype
template<>
__global__ void naive_verticalIntegral(double*  in, int2* out, int N, int ZN) {
	int xyi = threadIdx.x + blockIdx.x * blockDim.x;
	if (xyi < N) {
		double temp = in[xyi * ZN];
		out[xyi * ZN] = make_int2(__double2loint(temp), __double2hiint(temp));
		double temp2 = 0;
		for (int ii = 1; ii < ZN; ++ii) {
			temp2 = temp + in[xyi * ZN + ii];
			out[xyi * ZN + ii] = make_int2(__double2loint(temp2), __double2hiint(temp2));
			temp = temp2;
		}
	}
}

template<typename T>
__global__ void verticalIntegral(T* prj, int ZN, int N) {
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < N) {
		int currentHead = idx * ZN;
		for (int ii = 1; ii < ZN; ++ii) {
			prj[currentHead + ii] += prj[currentHead + ii - 1];
		}
	}
}

template<typename T>
__global__ void horizontalIntegral(T* prj, int DNU, int DNV, int PN) {
	int idv = threadIdx.x + blockIdx.x * blockDim.x;
	int pIdx = threadIdx.y + blockIdx.y * blockDim.y;
	if (idv < DNV && pIdx < PN)	{
		int headPtr = pIdx * DNU * DNV + idv;
		for (int ii = 1; ii < DNU; ++ii) {
			prj[headPtr + ii * DNV] += prj[headPtr + (ii - 1) * DNV];
		}
	}
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
			assert(curDetL.x != cursour.x);
			assert(curDetR.x != cursour.x);
			assert(curDet.x != cursour.x);
			factL = (curDetL.y - cursour.y) / (curDetL.x - cursour.x);
			factR = (curDetR.y - cursour.y) / (curDetR.x - cursour.x);
			factU = (curDet.w - cursour.z) / (curDet.x - cursour.x);
			factD = (curDet.z - cursour.z) / (curDet.x - cursour.x);
			assert(dir.x != 0);
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
				assert(intersectLength != 0 && intersectHeight != 0);

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
			assert(curDetL.y != cursour.y);
			assert(curDetR.y != cursour.y);
			assert(curDet.y != cursour.y);

			factL = (curDetL.x - cursour.x) / (curDetL.y - cursour.y);
			factR = (curDetR.x - cursour.x) / (curDetR.y - cursour.y);
			factU = (curDet.w - cursour.z) / (curDet.y - cursour.y);
			factD = (curDet.z - cursour.z) / (curDet.y - cursour.y);
			assert(dir.y != 0);
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
				assert(intersectLength != 0 && intersectHeight != 0);

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

thrust::device_vector<float3> getCosSinZPosFromAngleAndZPos(float* hangs, float* hzPos, int PN, const float x0, const float y0, const float z0) {
	thrust::device_vector<float> angs(hangs, hangs + PN);
	thrust::device_vector<float> zPos(hzPos, hzPos + PN);

	thrust::device_vector<float3> cossinZT(PN);
	thrust::transform(
		thrust::make_zip_iterator(thrust::make_tuple(angs.begin(), zPos.begin())),
		thrust::make_zip_iterator(thrust::make_tuple(angs.end(), zPos.end())),
		cossinZT.begin(), CTMBIR::ConstantForBackProjection4(x0, y0, z0));
	return cossinZT;
}

void DD3_gpu_proj_branchless_sat2d(
	float x0, float y0, float z0,
	int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN,
	float* vol, float* hprj, float dx, float dz, byte* mask, int gpunum) {

	maskingVolume(vol, mask, XN, YN, ZN);
	CUDA_CHECK_RETURN(cudaSetDevice(gpunum));
	CUDA_CHECK_RETURN(cudaDeviceReset());
	
	thrust::device_vector<float> d_bxds = getDD3Boundaries(xds, DNU);
	thrust::device_vector<float> d_byds = getDD3Boundaries(yds, DNU);
	thrust::device_vector<float> d_bzds = getDD3Boundaries(zds, DNV);

	float objCntIdxX = (XN - 1.0) * 0.5 - imgXCenter / dx;
	float objCntIdxY = (YN - 1.0) * 0.5 - imgYCenter / dx;
	float objCntIdxZ = (ZN - 1.0) * 0.5 - imgZCenter / dz;

	VolumeToTextures volumeToTexture(vol, XN, YN, ZN, cudaAddressModeClamp);

	thrust::device_vector<float> prj(DNU * DNV * PN, 0);
	thrust::device_vector<float> d_xds(xds, xds + DNU);
	thrust::device_vector<float> d_yds(yds, yds + DNU);
	thrust::device_vector<float> d_zds(zds, zds + DNV);

	thrust::device_vector<float3> cossinZT = getCosSinZPosFromAngleAndZPos(hangs, hzPos, PN, x0, y0, z0);
	
	dim3 blk(BLKX, BLKY, BLKZ);
	dim3 gid(
		(DNV + blk.x - 1) / blk.x,
		(DNU + blk.y - 1) / blk.y,
		(PN + blk.z - 1) / blk.z);

	DD3_gpu_proj_branchless_sat2d_ker << <gid, blk >> >(volumeToTexture.getTextureZXY(), volumeToTexture.getTextureZYX(),
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
	maskingVolume(vol, mask, XN, YN, ZN);

	CUDA_CHECK_RETURN(cudaSetDevice(gpunum));
	CUDA_CHECK_RETURN(cudaDeviceReset());

	float objCntIdxX = (XN - 1.0f) * 0.5f - imgXCenter / dx;
	float objCntIdxY = (YN - 1.0f) * 0.5f - imgYCenter / dx;
	float objCntIdxZ = (ZN - 1.0f) * 0.5f - imgZCenter / dz;

	cudaExtent volumeSize;
	volumeSize.width = ZN;
	volumeSize.height = XN;
	volumeSize.depth = YN;

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	cudaArray* d_volumeArray;
	CUDA_CHECK_RETURN(cudaMalloc3DArray(&d_volumeArray, &channelDesc, volumeSize));

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

	CUDA_CHECK_RETURN(cudaCreateTextureObject(&texObj, &resDesc, &texDesc, nullptr));

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
	CUDA_CHECK_RETURN(cudaDestroyTextureObject(texObj));

	CUDA_CHECK_RETURN(cudaFreeArray(d_volumeArray));
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
	cudaSetDevice(gpunum);
	cudaDeviceReset();

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
}

//
//
//__global__ void DD3_gpu_proj_volumerendering_ker(
//	cudaTextureObject_t volTex,
//	float* proj,
//	float x0, float y0, float z0,
//	float* d_xds, float* d_yds, float* d_zds,
//	float2* dirCent,
//	float* cosTs,
//	float* sinTs,
//	float objCntIdxX, float objCntIdxY, float objCntIdxZ,
//	float dx, float dz,
//	int XN, int YN, int ZN,
//	int DNU, int DNV, int PN,
//	float* angs, float* zPos, float3* gcursour, float* gcosGamma)
//{
//	int detIdV = threadIdx.x + blockIdx.x * blockDim.x;
//	int detIdU = threadIdx.y + blockIdx.y * blockDim.y;
//	int angIdx = threadIdx.z + blockIdx.z * blockDim.z;
//
//	__shared__ float cosGamma[BLKX];
//	__shared__ float2 sdirCent[BLKZ][BLKY];
//	__shared__ float2 normDir[BLKZ][BLKY];
//	__shared__ float curang;
//	__shared__ float cosT[BLKZ];
//	__shared__ float sinT[BLKZ];
//	cosGamma[threadIdx.x] = gcosGamma[detIdV];
//	sdirCent[threadIdx.z][threadIdx.y] = dirCent[angIdx * DNU + detIdU];
//	float l1 = sqrtf(sdirCent[threadIdx.z][threadIdx.y].x * sdirCent[threadIdx.z][threadIdx.y].x + sdirCent[threadIdx.z][threadIdx.y].y * sdirCent[threadIdx.z][threadIdx.y].y);
//	assert(l1 != 0);
//	normDir[threadIdx.z][threadIdx.y].x = sdirCent[threadIdx.z][threadIdx.y].x / l1;
//	normDir[threadIdx.z][threadIdx.y].y = sdirCent[threadIdx.z][threadIdx.y].y / l1;
//
//	curang = angs[angIdx];
//	cosT[threadIdx.z] = cosTs[angIdx];
//	sinT[threadIdx.z] = sinTs[angIdx];
//	__syncthreads();
//	if (detIdU < DNU && detIdV < DNV && angIdx < PN)
//	{
//		float zP = zPos[angIdx];
//		float3 curSour = make_float3(
//			x0 * cosT[threadIdx.z] - y0 * sinT[threadIdx.z],
//			x0 * sinT[threadIdx.z] + y0 * cosT[threadIdx.z],
//			z0 + zP);
//		float initDetX = d_xds[detIdU];
//		float initDetY = d_yds[detIdU];
//		float3 curDet = make_float3(
//			initDetX * cosT[threadIdx.z] - initDetY * sinT[threadIdx.z],
//			initDetX * sinT[threadIdx.z] + initDetY * cosT[threadIdx.z],
//			d_zds[detIdV] + zP);
//
//		float3 dir = curDet - curSour;
//		float summ = 0;
//		float obj = 0;
//		float real, realZ;
//		float invdz = 1.0 / dz;
//		float invdx = 1.0 / dx;
//
//		float constVal(1.0);
//		float fact1(1.0f), fact2(1.0f);
//		if ((curang > PI_4 && curang <= PI_3_4) || (curang > PI_5_4 && curang <= PI_7_4))
//		{
//			summ = 0;
//			assert(dir.x != 0);
//			fact1 = dir.y / dir.x;
//			fact2 = dir.z / dir.x;
//			assert(abs(normDir[threadIdx.z][threadIdx.y].x) * cosGamma[threadIdx.x] != 0);
//			constVal = dx / (abs(normDir[threadIdx.z][threadIdx.y].x) * cosGamma[threadIdx.x]);
//			obj = (-objCntIdxX) * dx;
//			real = fmaf(obj - curDet.x, fact1, curDet.y);
//			realZ = fmaf(obj - curDet.x, fact2, curDet.z);
//
//			fact1 = dx * fact1;
//			fact2 = dx * fact2;
//			for (int ii = 0; ii < XN; ++ii)
//			{
//				summ += tex3D<float>(volTex, fmaf(realZ, invdz, objCntIdxZ + 0.5), ii + 0.5, fmaf(real, invdx, objCntIdxY + 0.5));
//				obj += dx;
//				real += fact1;
//				realZ += fact2;
//
//			}
//			__syncthreads();
//			proj[(angIdx * DNU + detIdU) * DNV + detIdV] = summ * constVal;
//		}
//		else
//		{
//			summ = 0;
//			assert(dir.y != 0);
//			fact1 = dir.x / dir.y;
//			fact2 = dir.z / dir.y;
//			assert(abs(normDir[threadIdx.z][threadIdx.y].y) * cosGamma[threadIdx.x] != 0);
//			constVal = dx / (abs(normDir[threadIdx.z][threadIdx.y].y) * cosGamma[threadIdx.x]);
//			obj = (-objCntIdxX) * dx;
//			real = fmaf(obj - curDet.y, fact1, curDet.x);
//			realZ = fmaf(obj - curDet.y, fact2, curDet.z);
//
//			fact1 = dx * fact1;
//			fact2 = dx * fact2;
//			for (int ii = 0; ii < YN; ++ii)
//			{
//				summ += tex3D<float>(volTex, fmaf(realZ, invdz, objCntIdxZ + 0.5), fmaf(real, invdx, objCntIdxX + 0.5), ii + 0.5);
//				obj += dx;
//				real += fact1;
//				realZ += fact2;
//
//			}
//			__syncthreads();
//			proj[(angIdx * DNU + detIdU) * DNV + detIdV] = summ * constVal;
//		}
//
//	}
//
//}
//
//void DD3_gpu_proj_volumerendering(
//	float x0, float y0, float z0,
//	int DNU, int DNV,
//	float* xds, float* yds, float* zds,
//	float imgXCenter, float imgYCenter, float imgZCenter,
//	float* hangs, float* hzPos, int PN,
//	int XN, int YN, int ZN,
//	float* vol, float* hprj, float dx, float dz,
//	byte* mask, int gpunum)
//{
//	for (int ii = 0; ii != XN * YN; ++ii)
//	{
//		byte v = mask[ii];
//		for (int jj = 0; jj != ZN; ++jj)
//		{
//			vol[ii * ZN + jj] = vol[ii * ZN + jj] * v;
//		}
//	}
//
//	float* bxds = new float[DNU + 1];
//	float* byds = new float[DNU + 1];
//	float* bzds = new float[DNV + 1];
//	DD3Boundaries(DNU + 1, xds, bxds);
//	DD3Boundaries(DNU + 1, yds, byds);
//	DD3Boundaries(DNV + 1, zds, bzds);
//
//	cudaSetDevice(gpunum);
//	cudaDeviceReset();
//	cudaStream_t streams[4];
//	cudaStreamCreate(&streams[0]);
//	cudaStreamCreate(&streams[1]);
//	cudaStreamCreate(&streams[2]);
//	cudaStreamCreate(&streams[3]);
//
//	float objCntIdxX = (XN - 1.0) * 0.5 - imgXCenter / dx;
//	float objCntIdxY = (YN - 1.0) * 0.5 - imgYCenter / dx;
//	float objCntIdxZ = (ZN - 1.0) * 0.5 - imgZCenter / dz;
//
//	cudaExtent volumeSize;
//	volumeSize.width = ZN;
//	volumeSize.height = XN;
//	volumeSize.depth = YN;
//
//	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
//	cudaArray* d_volumeArray;
//	cudaMalloc3DArray(&d_volumeArray, &channelDesc, volumeSize);
//	cudaMemcpy3DParms copyParams = { 0 };
//	copyParams.srcPtr = make_cudaPitchedPtr((void*)vol, volumeSize.width * sizeof(float), volumeSize.width, volumeSize.height);
//	copyParams.dstArray = d_volumeArray;
//	copyParams.extent = volumeSize;
//	copyParams.kind = cudaMemcpyHostToDevice;
//	cudaMemcpy3D(&copyParams);
//
//	cudaResourceDesc resDesc;
//	memset(&resDesc, 0, sizeof(resDesc));
//	resDesc.resType = cudaResourceTypeArray;
//	resDesc.res.array.array = d_volumeArray;
//
//	cudaTextureDesc texDesc;
//	memset(&texDesc, 0, sizeof(texDesc));
//	texDesc.addressMode[0] = cudaAddressModeBorder;
//	texDesc.addressMode[1] = cudaAddressModeBorder;
//	texDesc.addressMode[2] = cudaAddressModeBorder;
//	texDesc.filterMode = cudaFilterModeLinear;
//	texDesc.readMode = cudaReadModeElementType;
//	texDesc.normalizedCoords = false;
//
//	cudaTextureObject_t texObj = 0;
//	cudaCreateTextureObject(&texObj, &resDesc, &texDesc, nullptr);
//
//	thrust::device_vector<float> prj(DNU * DNV * PN, 0);
//	thrust::device_vector<float> angs(hangs, hangs + PN);
//	thrust::device_vector<float> zPos(hzPos, hzPos + PN);
//	thrust::device_vector<float> d_xds(xds, xds + DNU);
//	thrust::device_vector<float> d_yds(yds, yds + DNU);
//	thrust::device_vector<float> d_zds(zds, zds + DNV);
//	thrust::device_vector<float3> cursour(PN);
//	thrust::device_vector<float2> dirCent(PN * DNU);
//	thrust::device_vector<float> cosTs(PN);
//	thrust::device_vector<float> sinTs(PN);
//	thrust::device_vector<float> cosGamma(DNV);
//
//	dim3 blkc(64, 16, 1);
//	dim3 gidc(
//		(DNV + blkc.x - 1) / blkc.x,
//		(DNU + blkc.y - 1) / blkc.y,
//		(PN + blkc.z - 1) / blkc.z);
//	dim3 blkc2(64, 16);
//	dim3 gidc2(
//		(DNU + blkc2.x - 1) / blkc2.x,
//		(PN + blkc2.y - 1) / blkc2.y);
//	//Calculate the constant values;
//
//	dim3 blk(BLKX, BLKY, BLKZ);
//	dim3 gid(
//		(DNV + blk.x - 1) / blk.x,
//		(DNU + blk.y - 1) / blk.y,
//		(PN + blk.z - 1) / blk.z);
//
//	DD3_gpu_proj_volumerendering_ker << <gid, blk >> >
//		(texObj,
//		thrust::raw_pointer_cast(&prj[0]),
//		x0, y0, z0,
//		thrust::raw_pointer_cast(&d_xds[0]),
//		thrust::raw_pointer_cast(&d_yds[0]),
//		thrust::raw_pointer_cast(&d_zds[0]),
//		thrust::raw_pointer_cast(&dirCent[0]),
//		thrust::raw_pointer_cast(&cosTs[0]),
//		thrust::raw_pointer_cast(&sinTs[0]),
//		objCntIdxX, objCntIdxY, objCntIdxZ,
//		dx, dz, XN, YN, ZN, DNU, DNV, PN,
//		thrust::raw_pointer_cast(&angs[0]),
//		thrust::raw_pointer_cast(&zPos[0]),
//		thrust::raw_pointer_cast(&cursour[0]),
//		thrust::raw_pointer_cast(&cosGamma[0]));
//
//	thrust::copy(prj.begin(), prj.end(), hprj);
//
//	cudaDestroyTextureObject(texObj);
//	cudaFreeArray(d_volumeArray);
//	prj.clear();
//	angs.clear();
//	cursour.clear();
//	cosTs.clear();
//	sinTs.clear();
//	dirCent.clear();
//	cosGamma.clear();
//	d_xds.clear();
//	d_yds.clear();
//	d_zds.clear();
//
//	delete[] bxds;
//	delete[] byds;
//	delete[] bzds;
//}
//
//


__global__ void DD3_gpu_proj_volumerendering_ker(cudaTextureObject_t texObj,
	float* prj, float3 s,
	const float3* __restrict cossinZT,
	const float* __restrict xds,
	const float* __restrict yds,
	const float* __restrict zds,
	float3 objCtrIdx,
	float dx, float dz, int XN, int YN, int ZN, int DNU, int DNV, int PN)
{
	int detIdV = threadIdx.x + blockIdx.x * blockDim.x;
	int detIdU = threadIdx.y + blockIdx.y * blockDim.y;
	int angIdx = threadIdx.z + blockIdx.z * blockDim.z;
	if (detIdV < DNV && detIdU < DNU && angIdx < PN)
	{
		float3 dir = cossinZT[angIdx];
		// Calculate current s position
		float3 cursour = make_float3(
			s.x * dir.x - s.y * dir.y,
			s.x * dir.y + s.y * dir.x,
			s.z + dir.z);
		// Calculate current detector position
		float3 curDet = make_float3(
			xds[detIdU] * dir.x - yds[detIdU] * dir.y,
			xds[detIdU] * dir.y + yds[detIdU] * dir.x,
			zds[detIdV] + dir.z);
		float leng = length(curDet - cursour);
		float lstp = sqrtf(dx * dx + dx * dx + dz * dz) / 4;
		int stpNum = ceil(leng / lstp);

		dir = (curDet - cursour) / leng;
		float3 curIdx;
		float summ = 0;

		for (int i = 0; i != stpNum; ++i) {
			cursour = cursour + dir * 0.25;
			curIdx.x = cursour.x / dx + objCtrIdx.x + 0.5; // 0.5 is used for texture object offset
			curIdx.y = cursour.y / dx + objCtrIdx.y + 0.5;
			curIdx.z = cursour.z / dz + objCtrIdx.z + 0.5;

			summ += tex3D<float>(texObj, curIdx.z, curIdx.x, curIdx.y) * lstp;
		}
		prj[(angIdx * DNU + detIdU) * DNV + detIdV] = summ;
	}
}

// Todo: Incorrect
void DD3_gpu_proj_volumerendering(float x0, float y0, float z0, int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN, int XN, int YN, int ZN,
	float* hvol, float* hprj, float dx, float dz, byte* mask, int gpunum)
{
	maskingVolume(hvol, mask, XN, YN, ZN);
	CUDA_CHECK_RETURN(cudaSetDevice(gpunum)); // choose the 2nd GPU for now (1st GPU used to be occupied on torb
	CUDA_CHECK_RETURN(cudaDeviceReset());

	thrust::device_vector<float> prj(DNU * DNV * PN, 0); // projection data
	thrust::device_vector<float> d_xds(xds, xds + DNU);
	thrust::device_vector<float> d_yds(yds, yds + DNU);
	thrust::device_vector<float> d_zds(zds, zds + DNV);

	thrust::device_vector<float> angs(hangs, hangs + PN);
	thrust::device_vector<float> zPos(hzPos, hzPos + PN);

	thrust::device_vector<float3> cossinZT(PN);
	thrust::transform(
		thrust::make_zip_iterator(thrust::make_tuple(angs.begin(), zPos.begin())),
		thrust::make_zip_iterator(thrust::make_tuple(angs.end(), zPos.end())),
		cossinZT.begin(), CTMBIR::ConstantForBackProjection<float>(x0, y0, z0));

	float3 objCtrIdx = make_float3(
		(XN - 1.0) * 0.5 - imgXCenter / dx,
		(YN - 1.0) * 0.5 - imgYCenter / dx,
		(ZN - 1.0) * 0.5 - imgZCenter / dz);

	//Configure BLOCKs for projection
	dim3 blk;
	dim3 gid;
	blk.x = BLKX; // det row index
	blk.y = BLKY; // det col index
	blk.z = BLKZ; // view index
	gid.x = (DNV + blk.x - 1) / blk.x;
	gid.y = (DNU + blk.y - 1) / blk.y;
	gid.z = (PN + blk.z - 1) / blk.z;

	/// Copy volumes to texture
	// Allocate CUDA array in device memory
	cudaTextureObject_t texObj;
	cudaArray* d_volumeArray = nullptr;

	createTextureObject<float>(texObj, d_volumeArray, ZN, XN, YN, hvol,
		cudaMemcpyHostToDevice, cudaAddressModeBorder,
		cudaFilterModeLinear, cudaReadModeElementType, false);

	DD3_gpu_proj_volumerendering_ker << <gid, blk >> >(texObj, getRawPtr(prj),
		make_float3(x0, y0, z0), getRawPtr(cossinZT),
		getRawPtr(d_xds), getRawPtr(d_yds), getRawPtr(d_zds),
		objCtrIdx, dx, dz, XN, YN, ZN, DNU, DNV, PN);

	thrust::copy(prj.begin(), prj.end(), hprj);
	destroyTextureObject(texObj, d_volumeArray);
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
			assert(dir.x != 0);
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
			assert(dir.y != 0);
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
	float dx, float dz, byte* mask, int gpunum) {
	maskingVolume(hvol, mask, XN, YN, ZN);

	CUDA_CHECK_RETURN(cudaSetDevice(gpunum));
	CUDA_CHECK_RETURN(cudaDeviceReset());
	
	thrust::device_vector<float> d_bxds = getDD3Boundaries(xds, DNU);
	thrust::device_vector<float> d_byds = getDD3Boundaries(yds, DNU);
	thrust::device_vector<float> d_bzds = getDD3Boundaries(zds, DNV);

	const int TOTVN = XN * YN * ZN;
	float objCntIdxX = (XN - 1.0) * 0.5 - imgXCenter / dx;
	float objCntIdxY = (YN - 1.0) * 0.5 - imgYCenter / dx;
	float objCntIdxZ = (ZN - 1.0) * 0.5 - imgZCenter / dz;

	VolumeToOneTexture volumeToOneTexture(hvol, XN, YN, ZN, cudaAddressModeBorder);

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
		volumeToOneTexture.getTexture(), 
		thrust::raw_pointer_cast(&prj[0]),
		make_float3(x0, y0, z0),
		thrust::raw_pointer_cast(&d_xds[0]),
		thrust::raw_pointer_cast(&d_yds[0]),
		thrust::raw_pointer_cast(&d_zds[0]),
		thrust::raw_pointer_cast(&cossinZT[0]),
		make_float3(objCntIdxX, objCntIdxY, objCntIdxZ),
		dx, dz, XN, YN, DNU, DNV, PN);
	thrust::copy(prj.begin(), prj.end(), hprj);
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
			assert(sqrt(dirY * dirY + dirX * dirX) != 0);
			float cosAlpha = abs(dirY / sqrt(dirY * dirY + dirX * dirX));
			float cosGamma = abs(sqrt((S2O + O2D) * (S2O + O2D) - initDetZ * initDetZ) / (S2O + O2D));

			assert(curDetLY != cursoury);
			assert(curDetRY != cursoury);
			assert(curDetDY != cursoury);
			assert(curDetUY != cursoury);

			float detPosLX = -cursoury * (curDetLX - cursourx) / (curDetLY - cursoury) + cursourx; //×ó±ßµãÔÚXOZÆœÃæÉÏµÄÍ¶Ó°;
			float detPosRX = -cursoury * (curDetRX - cursourx) / (curDetRY - cursoury) + cursourx;
			float detPosDZ = -cursoury * (curDetDZ - cursourz) / (curDetDY - cursoury) + cursourz;
			float detPosUZ = -cursoury * (curDetUZ - cursourz) / (curDetUY - cursoury) + cursourz;

			float detprojLength = abs(detPosLX - detPosRX);
			float detprojHeight = abs(detPosUZ - detPosDZ);
			assert(detprojLength != 0);
			assert(detprojHeight != 0);

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
				assert(temp != 0);
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
			assert(sqrt(dirY * dirY + dirX * dirX) != 0);
			float cosAlpha = abs(dirX / sqrt(dirY * dirY + dirX * dirX));
			float cosGamma = abs(sqrt((S2O + O2D) * (S2O + O2D) - initDetZ * initDetZ) / (S2O + O2D));

			assert(curDetLX != cursourx);
			assert(curDetRX != cursourx);
			assert(curDetDX != cursourx);
			assert(curDetUX != cursourx);

			float detPosLY = -cursourx * (curDetLY - cursoury) / (curDetLX - cursourx) + cursoury; //×ó±ßµãÔÚXOZÆœÃæÉÏµÄÍ¶Ó°;
			float detPosRY = -cursourx * (curDetRY - cursoury) / (curDetRX - cursourx) + cursoury;
			float detPosDZ = -cursourx * (curDetDZ - cursourz) / (curDetDX - cursourx) + cursourz;
			float detPosUZ = -cursourx * (curDetUZ - cursourz) / (curDetUX - cursourx) + cursourz;

			float detprojLength = abs(detPosLY - detPosRY);
			float detprojHeight = abs(detPosUZ - detPosDZ);
			assert(detprojLength != 0);
			assert(detprojHeight != 0);

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
				assert(temp != 0);
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


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Using DOUBLE PRECISION BRANCHLESS GPU
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// \brief Kernel function of Double precision based branchless DD with 2D SAT
__global__ void DD3_gpu_proj_doubleprecisionbranchless_ker(
	cudaTextureObject_t volTex1, // volume SAT in ZXY order
	cudaTextureObject_t volTex2, // volume SAT in ZYX order
	double* proj, // projection data
	double3 s, // initial source position
	const double3* __restrict cossinZT, // bind (cosine, sine, zshift)
	const double* __restrict xds,  //
	const double* __restrict yds,  //
	const double* __restrict zds,  // detector cells center positions
	const double* __restrict bxds,
	const double* __restrict byds,
	const double* __restrict bzds, // detector boundary positions
	double3 objCntIdx, // object center index
	double dx, double dz, // pixel size in xy plane and Z direction
	int XN, int YN, // pixel # in XY plane
	int DNU, int DNV, // detector cell # in xy plane and Z direction
	int PN) // view #
{
	int detIdV = threadIdx.x + blockIdx.x * blockDim.x;
	int detIdU = threadIdx.y + blockIdx.y * blockDim.y;
	int angIdx = threadIdx.z + blockIdx.z * blockDim.z;

	__shared__ double _xds[BLKY];
	__shared__ double _yds[BLKY];
	_xds[threadIdx.y] = xds[detIdU];
	_yds[threadIdx.y] = yds[detIdU];
	__syncthreads();

	if (detIdU < DNU && detIdV < DNV && angIdx < PN)
	{
		double3 dir = cossinZT[angIdx];
		double3 cursour = make_double3(
			s.x * dir.x - s.y * dir.y,
			s.x * dir.y + s.y * dir.x,
			s.z + dir.z);
		s = cossinZT[angIdx];
		double summ = _xds[threadIdx.y] * s.x - _yds[threadIdx.y] * s.y;
		double obj = _xds[threadIdx.y] * s.y + _yds[threadIdx.y] * s.x;
		double realL = bxds[detIdU];
		double realR = byds[detIdU];
		double realU = bxds[detIdU + 1];
		double realD = byds[detIdU + 1]; // intersection coordinates (mm); float2 is equv to (obj1,obj2) above
		double2 curDetL = make_double2(
			realL * s.x - realR * s.y,
			realL * s.y + realR * s.x);

		double2 curDetR = make_double2(
			realU * s.x - realD * s.y,
			realU * s.y + realD * s.x);
		double4 curDet = make_double4(summ, obj, bzds[detIdV] + s.z, bzds[detIdV + 1] + s.z); //(center x, center y, lower z, upper z)

		dir = normalize(make_double3(
			summ,
			obj,
			zds[detIdV] + s.z) - cursour);

		summ = 0; // to accumulate projection value
		obj = 0; // slice location (mm) along the ray tracing direction TODO: is this variable needed?

		double intersectLength, intersectHeight;
		double invdz = 1.0 / dz;
		double invdx = 1.0 / dx;


		double factL(1.0f); // dy/dx for (0,pi/4)
		double factR(1.0f);
		double factU(1.0f);
		double factD(1.0f);
		double constVal = 0;

		int crealD, crealR, crealU, crealL;
		int frealD, frealR, frealU, frealL;

		if (abs(s.x) <= abs(s.y))
		{
			summ = 0;
			// a few book keeping variables
			assert(curDetL.x != cursour.x);
			assert(curDetR.x != cursour.x);
			assert(curDet.x != cursour.x);

			factL = (curDetL.y - cursour.y) / (curDetL.x - cursour.x);
			factR = (curDetR.y - cursour.y) / (curDetR.x - cursour.x);
			factU = (curDet.w - cursour.z) / (curDet.x - cursour.x);
			factD = (curDet.z - cursour.z) / (curDet.x - cursour.x);

			assert(dir.x != 0);
			constVal = dx * dx * dz / (abs(dir.x));
#pragma unroll
			for (int ii = 0; ii < XN; ii++)
			{
				obj = (ii - objCntIdx.x) * dx;

				realL = (obj - curDetL.x) * factL + curDetL.y;
				realR = (obj - curDetR.x) * factR + curDetR.y;
				realU = (obj - curDet.x) * factU + curDet.w;
				realD = (obj - curDet.x) * factD + curDet.z;

				intersectLength = realR - realL;
				intersectHeight = realU - realD;
				assert(intersectLength != 0 && intersectHeight != 0);

				// 1D LUT to address inaccuracies in texture coordinates
				realD = realD * invdz + objCntIdx.z + 1;
				realR = realR * invdx + objCntIdx.y + 1;
				realU = realU * invdz + objCntIdx.z + 1;
				realL = realL * invdx + objCntIdx.y + 1;

				crealD = ceil(realD);
				crealR = ceil(realR);
				crealU = ceil(realU);
				crealL = ceil(realL);

				frealD = floor(realD);
				frealR = floor(realR);
				frealU = floor(realU);
				frealL = floor(realL);


				summ +=
					(bilerp(
						tex3D<int2>(volTex2, frealD, frealL, ii + 0.5),
						tex3D<int2>(volTex2, frealD, crealL, ii + 0.5),
						tex3D<int2>(volTex2, crealD, frealL, ii + 0.5),
						tex3D<int2>(volTex2, crealD, crealL, ii + 0.5),
						realL - frealL, realD - frealD) +
						bilerp(
							tex3D<int2>(volTex2, frealU, frealR, ii + 0.5),
							tex3D<int2>(volTex2, frealU, crealR, ii + 0.5),
							tex3D<int2>(volTex2, crealU, frealR, ii + 0.5),
							tex3D<int2>(volTex2, crealU, crealR, ii + 0.5),
							realR - frealR, realU - frealU) -
						bilerp(
							tex3D<int2>(volTex2, frealD, frealR, ii + 0.5),
							tex3D<int2>(volTex2, frealD, crealR, ii + 0.5),
							tex3D<int2>(volTex2, crealD, frealR, ii + 0.5),
							tex3D<int2>(volTex2, crealD, crealR, ii + 0.5),
							realR - frealR, realD - frealD) -
						bilerp(
							tex3D<int2>(volTex2, frealU, frealL, ii + 0.5),
							tex3D<int2>(volTex2, frealU, crealL, ii + 0.5),
							tex3D<int2>(volTex2, crealU, frealL, ii + 0.5),
							tex3D<int2>(volTex2, crealU, crealL, ii + 0.5),
							realL - frealL, realU - frealU)) / (intersectLength * intersectHeight);
			}
			__syncthreads();
			proj[(angIdx * DNU + detIdU) * DNV + detIdV] = summ * constVal;
		}
		else
		{
			summ = 0;
			assert(curDetL.y - cursour.y);
			assert(curDetR.y - cursour.y);
			assert(curDet.y - cursour.y);

			factL = (curDetL.x - cursour.x) / (curDetL.y - cursour.y);
			factR = (curDetR.x - cursour.x) / (curDetR.y - cursour.y);
			factU = (curDet.w - cursour.z) / (curDet.y - cursour.y);
			factD = (curDet.z - cursour.z) / (curDet.y - cursour.y);

			assert(dir.y);
			constVal = dx * dx * dz / (abs(dir.y));
#pragma unroll
			for (int jj = 0; jj < YN; jj++)
			{
				obj = (jj - objCntIdx.y) * dx;
				realL = (obj - curDetL.y) * factL + curDetL.x;
				realR = (obj - curDetR.y) * factR + curDetR.x;
				realU = (obj - curDet.y) * factU + curDet.w;
				realD = (obj - curDet.y) * factD + curDet.z;

				intersectLength = realR - realL;
				intersectHeight = realU - realD;
				assert(intersectLength != 0 && intersectHeight != 0);

				realD = realD * invdz + objCntIdx.z + 1;
				realR = realR * invdx + objCntIdx.x + 1;
				realU = realU * invdz + objCntIdx.z + 1;
				realL = realL * invdx + objCntIdx.x + 1;


				crealD = ceil(realD);
				crealR = ceil(realR);
				crealU = ceil(realU);
				crealL = ceil(realL);

				frealD = floor(realD);
				frealR = floor(realR);
				frealU = floor(realU);
				frealL = floor(realL);

				summ += (bilerp(
					tex3D<int2>(volTex1, frealD, frealL, jj + 0.5),
					tex3D<int2>(volTex1, frealD, crealL, jj + 0.5),
					tex3D<int2>(volTex1, crealD, frealL, jj + 0.5),
					tex3D<int2>(volTex1, crealD, crealL, jj + 0.5),
					realL - frealL, realD - frealD) +
					bilerp(
						tex3D<int2>(volTex1, frealU, frealR, jj + 0.5),
						tex3D<int2>(volTex1, frealU, crealR, jj + 0.5),
						tex3D<int2>(volTex1, crealU, frealR, jj + 0.5),
						tex3D<int2>(volTex1, crealU, crealR, jj + 0.5),
						realR - frealR, realU - frealU) -
					bilerp(
						tex3D<int2>(volTex1, frealD, frealR, jj + 0.5),
						tex3D<int2>(volTex1, frealD, crealR, jj + 0.5),
						tex3D<int2>(volTex1, crealD, frealR, jj + 0.5),
						tex3D<int2>(volTex1, crealD, crealR, jj + 0.5),
						realR - frealR, realD - frealD) -
					bilerp(
						tex3D<int2>(volTex1, frealU, frealL, jj + 0.5),
						tex3D<int2>(volTex1, frealU, crealL, jj + 0.5),
						tex3D<int2>(volTex1, crealU, frealL, jj + 0.5),
						tex3D<int2>(volTex1, crealU, crealL, jj + 0.5),
						realL - frealL, realU - frealU)) / (intersectLength * intersectHeight);
			}
			__syncthreads();
			proj[(angIdx * DNU + detIdU) * DNV + detIdV] = summ * constVal;
		}

	}
}

// \brief C interface of Double Precision Branchless DD with 2D SAT
void DD3_gpu_proj_doubleprecisionbranchless(
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
	//Pre compute mask.*vol;
	maskingVolume(vol, mask, XN, YN, ZN);

	float* bxds = new float[DNU + 1];
	float* byds = new float[DNU + 1];
	float* bzds = new float[DNV + 1];
	DD3Boundaries(DNU + 1, xds, bxds);
	DD3Boundaries(DNU + 1, yds, byds);
	DD3Boundaries(DNV + 1, zds, bzds);

	CUDA_CHECK_RETURN(cudaSetDevice(gpunum));
	CUDA_CHECK_RETURN(cudaDeviceReset());

	cudaStream_t streams[4];
	CUDA_CHECK_RETURN(cudaStreamCreate(&streams[0]));
	CUDA_CHECK_RETURN(cudaStreamCreate(&streams[1]));
	CUDA_CHECK_RETURN(cudaStreamCreate(&streams[2]));
	CUDA_CHECK_RETURN(cudaStreamCreate(&streams[3]));

	int TOTVN = XN * YN * ZN;
	double objCntIdxX = (XN - 1.0) * 0.5 - imgXCenter / dx;
	double objCntIdxY = (YN - 1.0) * 0.5 - imgYCenter / dx;
	double objCntIdxZ = (ZN - 1.0) * 0.5 - imgZCenter / dz;

	thrust::device_vector<float> in(vol, vol + TOTVN); // original img volume
	thrust::device_vector<double> in_ZXY((ZN + 1) * (XN + 1) * YN, 0); //
	thrust::device_vector<double> in_ZYX((ZN + 1) * (YN + 1) * XN, 0); // transposed img volume

	dim3 blk(64, 16, 1);
	dim3 gid(
		(ZN + blk.x - 1) / blk.x,
		(XN + blk.y - 1) / blk.y,
		(YN + blk.z - 1) / blk.z);
	// copy to original and transposed image volume with left- and top-side boarder padding to be consistent with SAT dimensions
	naive_copyToTwoVolumes<float, double> << <gid, blk >> >(
		thrust::raw_pointer_cast(&in[0]),
		thrust::raw_pointer_cast(&in_ZXY[0]),
		thrust::raw_pointer_cast(&in_ZYX[0]), XN, YN, ZN);
	in.clear();

	thrust::device_vector<double> in_ZXY_summ1((ZN + 1) * (XN + 1) * YN, 0);
	thrust::device_vector<int2> in_ZXY_summ((ZN + 1) * (XN + 1) * YN);


	blk.x = 64;							blk.y = 1;		blk.z = 1;
	gid.x = (ZN + blk.x) / blk.x;		gid.y = 1;		gid.z = 1;

	dim3 blk2(64);
	dim3 gid2((YN + blk2.x) / blk2.x);
	dim3 blk3(64);
	dim3 gid3((XN + blk3.x) / blk3.x);

	// compute SAT for the original img volume
	for (int jj = 0; jj != YN; ++jj)
	{
		// for each Y slice
		naive_herizontalIntegral << <gid, blk, 0, streams[0] >> >(
			thrust::raw_pointer_cast(&in_ZXY[0]) + jj * (ZN + 1) * (XN + 1),
			thrust::raw_pointer_cast(&in_ZXY_summ1[0]) + jj * (ZN + 1) * (XN + 1), XN + 1, ZN + 1);
		naive_verticalIntegral << <gid2, blk2, 0, streams[0] >> >(
			thrust::raw_pointer_cast(&in_ZXY_summ1[0]) + jj * (ZN + 1) * (XN + 1),
			thrust::raw_pointer_cast(&in_ZXY_summ[0]) + jj * (ZN + 1) * (XN + 1), XN + 1, ZN + 1);
	}
	in_ZXY.clear();
	in_ZXY_summ1.clear();


	cudaArray* d_volumeArray1 = nullptr;
	cudaTextureObject_t texObj1;

	createTextureObject<int2>(texObj1, d_volumeArray1, ZN + 1, XN + 1, YN,
		thrust::raw_pointer_cast(&in_ZXY_summ[0]),
		cudaMemcpyDeviceToDevice, cudaAddressModeClamp, cudaFilterModePoint,
		cudaReadModeElementType, false);
	in_ZXY_summ.clear();
	
	thrust::device_vector<double> in_ZYX_summ1((ZN + 1) * (YN + 1) * XN, 0); // SAT for the transposed img volume
	thrust::device_vector<int2> in_ZYX_summ((ZN + 1) * (YN + 1) * XN);
	// compute SAT for the transposed img volume
	for (int ii = 0; ii != XN; ++ii)
	{
		// for each X slice
		naive_herizontalIntegral << <gid, blk, 0, streams[1] >> >(
			thrust::raw_pointer_cast(&in_ZYX[0]) + ii * (ZN + 1) * (YN + 1),
			thrust::raw_pointer_cast(&in_ZYX_summ1[0]) + ii * (ZN + 1) * (YN + 1), YN + 1, ZN + 1);
		naive_verticalIntegral << <gid3, blk3, 0, streams[1] >> >(
			thrust::raw_pointer_cast(&in_ZYX_summ1[0]) + ii * (ZN + 1) * (YN + 1),
			thrust::raw_pointer_cast(&in_ZYX_summ[0]) + ii * (ZN + 1) * (YN + 1), YN + 1, ZN + 1);
	}
	in_ZYX.clear();
	in_ZYX_summ1.clear();

	cudaArray* d_volumeArray2 = nullptr;
	cudaTextureObject_t texObj2;

	createTextureObject<int2>(texObj2, d_volumeArray2, ZN + 1, YN + 1, XN,
		thrust::raw_pointer_cast(&in_ZYX_summ[0]),
		cudaMemcpyDeviceToDevice, cudaAddressModeClamp, cudaFilterModePoint,
		cudaReadModeElementType, false);
	in_ZYX_summ.clear();

	thrust::device_vector<double> prj(DNU * DNV * PN, 0);
	thrust::device_vector<double> angs(hangs, hangs + PN);
	thrust::device_vector<double> zPos(hzPos, hzPos + PN);
	thrust::device_vector<double> d_xds(xds, xds + DNU);
	thrust::device_vector<double> d_yds(yds, yds + DNU);
	thrust::device_vector<double> d_zds(zds, zds + DNV);
	thrust::device_vector<double> d_bxds(bxds, bxds + DNU + 1);
	thrust::device_vector<double> d_byds(byds, byds + DNU + 1);
	thrust::device_vector<double> d_bzds(bzds, bzds + DNV + 1);


	// constant values for DD calculation
	thrust::device_vector<double3> cossinZT(PN);
	thrust::transform(
		thrust::make_zip_iterator(thrust::make_tuple(angs.begin(), zPos.begin())),
		thrust::make_zip_iterator(thrust::make_tuple(angs.end(), zPos.end())),
		cossinZT.begin(), CTMBIR::ConstantForBackProjection<double>(x0, y0, z0));

	//precalculate all constant values in CUDA
	dim3 blkc(64, 16, 1);
	dim3 gidc(
		(DNV + blkc.x) / blkc.x,
		(DNU + blkc.y) / blkc.y,
		(PN + blkc.z - 1) / blkc.z);


	//Configure BLOCKs for projection
	blk.x = BLKX; // det row index
	blk.y = BLKY; // det col index
	blk.z = BLKZ; // view index
	gid.x = (DNV + blk.x - 1) / blk.x;
	gid.y = (DNU + blk.y - 1) / blk.y;
	gid.z = (PN + blk.z - 1) / blk.z;

	//Projection kernel
	DD3_gpu_proj_doubleprecisionbranchless_ker << <gid, blk >> >(texObj1, texObj2,
		thrust::raw_pointer_cast(&prj[0]), make_double3(x0, y0, z0),
		thrust::raw_pointer_cast(&cossinZT[0]),
		thrust::raw_pointer_cast(&d_xds[0]),
		thrust::raw_pointer_cast(&d_yds[0]),
		thrust::raw_pointer_cast(&d_zds[0]),
		thrust::raw_pointer_cast(&d_bxds[0]),
		thrust::raw_pointer_cast(&d_byds[0]),
		thrust::raw_pointer_cast(&d_bzds[0]),
		make_double3(objCntIdxX, objCntIdxY, objCntIdxZ),
		dx, dz, XN, YN, DNU, DNV, PN);
	thrust::copy(prj.begin(), prj.end(), hprj);

	CUDA_CHECK_RETURN(cudaDestroyTextureObject(texObj1));
	CUDA_CHECK_RETURN(cudaDestroyTextureObject(texObj2));

	destroyTextureObject(texObj1, d_volumeArray1);
	destroyTextureObject(texObj2, d_volumeArray2);
	prj.clear();
	angs.clear();
	zPos.clear();

	d_xds.clear();
	d_yds.clear();
	d_zds.clear();
	d_bxds.clear();
	d_byds.clear();
	d_bzds.clear();


	delete[] bxds;
	delete[] byds;
	delete[] bzds;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// END DOUBLE PRECISION BRANCHLESS 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DD3 BRANCHES PROJECTION
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void DD3Proj_branches_ker(float* proj, const float* vol,
	float x0, float y0, float z0,
	const float* xds, const float* yds,
	const float* bxds, const float* byds, const float* bzds,
	float3* cossinZT,
	float objCntIdxX, float objCntIdxY, float objCntIdxZ,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	float dx, float dz)
{
	const int detIdV = threadIdx.x + blockIdx.x * blockDim.x; //Cannot be changed
	const int detIdU = threadIdx.y + blockIdx.y * blockDim.y; // Cannot be changed
	const int angIdx = threadIdx.z + blockIdx.z * blockDim.z; // Cannot be changed
	if (detIdU < DNU && detIdV < DNV && angIdx < PN)
	{
		float cosT = cossinZT[angIdx].x;
		float sinT = cossinZT[angIdx].y;
		float zPos = cossinZT[angIdx].z;
		float3 cursour = make_float3(x0 * cosT - y0 * sinT, x0 * sinT + y0 * cosT, z0 + zPos); // Cannot be modified
		float summ = 0; // Cannot be modified

		float detPosLX = bxds[detIdU];
		float detPosRX = bxds[detIdU + 1];
		float detPosLY = byds[detIdU];
		float detPosRY = byds[detIdU + 1];
		float detPosDZ = bzds[detIdV];
		float detPosUZ = detPosDZ + dz;

		float curDirLX = detPosLX * cosT - detPosLY * sinT - cursour.x;
		float curDirLY = detPosLX * sinT + detPosLY * cosT - cursour.y;

		float curDirRX = detPosRX * cosT - detPosRY * sinT - cursour.x;
		float curDirRY = detPosRX * sinT + detPosRY * cosT - cursour.y;
		float curDirDZ = detPosDZ + zPos - cursour.z;
		float curDirUZ = detPosUZ + zPos - cursour.z;

		float dirX = xds[detIdU] * cosT - yds[detIdU] * sinT - cursour.x; // cannot be modified
		float dirY = xds[detIdU] * sinT + yds[detIdU] * cosT - cursour.y; // cannot be modified
		float dirZ = detPosDZ - 0.5 * dz - z0; // used cannot be modified

		if (fabsf(dirX) < fabsf(dirY))
		{
			assert(curDirLY != 0);
			assert(curDirRY != 0);
			assert(dirY != 0);
			detPosLX = -cursour.y * curDirLX / curDirLY + cursour.x;
			detPosRX = -cursour.y * curDirRX / curDirRY + cursour.x;
			detPosDZ = -cursour.y * curDirDZ / dirY + cursour.z;
			detPosUZ = -cursour.y * curDirUZ / dirY + cursour.z;
			assert(detPosLX != detPosRX);
			assert(detPosUZ != detPosDZ);
			dirZ = dx / fabsf(dirY / sqrtf(dirY * dirY + dirX * dirX + dirZ * dirZ) * (detPosLX - detPosRX) * (detPosUZ - detPosDZ));

			if (detPosLX > detPosRX)
			{
				zPos = detPosLX;
				detPosLX = detPosRX;
				detPosRX = zPos;
			}
			for (int jj = 0; jj < YN; jj++)
			{
				y0 = (jj - objCntIdxY) * dx - cursour.y;
				zPos = y0 / curDirLY;

				cosT = zPos * curDirLX + cursour.x;
				sinT = zPos * curDirRX + cursour.x;
				x0 = zPos * curDirDZ + cursour.z;
				z0 = zPos * curDirUZ + cursour.z;
				if (cosT > sinT)
				{
					zPos = cosT;
					cosT = sinT;
					sinT = zPos;
				}

				int minXIdx = floor(cosT / dx + objCntIdxX) - 1; // XN/2.0
				int maxXIdx = ceil(sinT / dx + objCntIdxX) + 1;
				int minZIdx = floor(x0 / dz + objCntIdxZ) - 1; // ZN/2.0
				int maxZIdx = ceil(z0 / dz + objCntIdxZ) + 1;
				if (maxXIdx < 0) { continue; }
				if (minXIdx > XN) { continue; }
				if (maxZIdx < 0) { continue; }
				if (minZIdx > ZN) { continue; }
				if (minXIdx < 0) { minXIdx = 0; }
				if (maxXIdx > XN) { maxXIdx = XN; }
				if (minZIdx < 0) { minZIdx = 0; }
				if (maxZIdx > ZN) { maxZIdx = ZN; }

				cosT = (cursour.x - (minXIdx - objCntIdxX) * dx) * cursour.y / y0 + cursour.x; // XN/2.0
				for (int ii = minXIdx; ii < maxXIdx; ii++)
				{
					sinT = cosT - dx * cursour.y / y0;
					dirX = intersectLength_device<float>(detPosLX, detPosRX, cosT, sinT);
					x0 = (cursour.z - (minZIdx - objCntIdxZ) * dz) * cursour.y / y0 + cursour.z;
					for (int kk = minZIdx; kk < maxZIdx; kk++)
					{
						z0 = x0 - dz * cursour.y / y0;
						dirY = intersectLength_device<float>(detPosDZ, detPosUZ, x0, z0);
						summ += vol[(jj * XN + ii) * ZN + kk] * (dirX * dirY);
						x0 = z0;
					}
					cosT = sinT;
				}
			}
			proj[(angIdx * DNU + detIdU) * DNV + detIdV] = summ * dirZ;
		}
		else
		{
			assert(curDirLX != 0);
			assert(curDirRX != 0);
			assert(dirX != 0);
			
			detPosLY = -cursour.x * curDirLY / curDirLX + cursour.y;
			detPosRY = -cursour.x * curDirRY / curDirRX + cursour.y;
			detPosDZ = -cursour.x * curDirDZ / dirX + cursour.z;
			detPosUZ = -cursour.x * curDirUZ / dirX + cursour.z;

			assert(detPosLY != detPosRY);
			assert(detPosUZ != detPosDZ);
			dirZ = dx / fabsf(dirX / sqrtf(dirY * dirY + dirX * dirX + dirZ * dirZ) * (detPosLY - detPosRY) * (detPosUZ - detPosDZ));

			if (detPosLY > detPosRY)
			{
				zPos = detPosLY;
				detPosLY = detPosRY;
				detPosRY = zPos;
			}

			for (int ii = 0; ii < XN; ii++)
			{
				x0 = (ii - objCntIdxX) * dx - cursour.x;
				zPos = x0 / curDirLX;

				cosT = zPos * curDirLY + cursour.y;
				sinT = zPos * curDirRY + cursour.y;
				y0 = zPos * curDirDZ + cursour.z;
				z0 = zPos * curDirUZ + cursour.z;
				if (cosT > sinT) { zPos = cosT; cosT = sinT; sinT = zPos; }

				int minYIdx = floor(cosT / dx + objCntIdxY) - 1;
				int maxYIdx = ceil(sinT / dx + objCntIdxY) + 1;
				int minZIdx = floor(y0 / dz + objCntIdxZ) - 1;
				int maxZIdx = ceil(z0 / dz + objCntIdxZ) + 1;
				if (maxYIdx < 0) { continue; }
				if (minYIdx > XN) { continue; }
				if (maxZIdx < 0) { continue; }
				if (minZIdx > ZN) { continue; }
				if (minYIdx < 0) { minYIdx = 0; }
				if (maxYIdx > XN) { maxYIdx = YN; }
				if (minZIdx < 0) { minZIdx = 0; }
				if (maxZIdx > ZN) { maxZIdx = ZN; }

				cosT = (cursour.y - (minYIdx - objCntIdxY) * dx) * cursour.x / x0 + cursour.y;
				for (int jj = minYIdx; jj < maxYIdx; jj++)
				{
					sinT = cosT - dx * cursour.x / x0;// (cursour.y - ((jj + 1) - YN / 2.0) * dx) * cursour.x / objdirX + cursour.y;
					dirX = intersectLength_device<float>(detPosLY, detPosRY, cosT, sinT);
					y0 = (cursour.z - (minZIdx - objCntIdxZ) * dz) * cursour.x / x0 + cursour.z;
					for (int kk = minZIdx; kk < maxZIdx; kk++)
					{
						z0 = y0 - dz * cursour.x / x0;
						dirY = intersectLength_device<float>(detPosDZ, detPosUZ, y0, z0);
						summ += vol[(jj * XN + ii) * ZN + kk] * (dirX * dirY);
						y0 = z0;
					}
					cosT = sinT;
				}
			}
			proj[(angIdx * DNU + detIdU) * DNV + detIdV] = summ * dirZ;
		}
	}
}

void DD3Proj_branches(float x0, float y0, float z0, int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN, float* hvol, float* hprj, float dx, float dz, byte* mask,
	int gpunum)
{
	CUDA_CHECK_RETURN(cudaSetDevice(gpunum));
	maskingVolume(hvol, mask, XN, YN, ZN);

	thrust::host_vector<float3> hcossinZT(PN);
	for (int i = 0; i != PN; ++i)
	{
		hcossinZT[i].x = cosf(hangs[i]);
		hcossinZT[i].y = sinf(hangs[i]);
		hcossinZT[i].z = hzPos[i];
	}
	thrust::device_vector<float3> cossinZT = hcossinZT;

	thrust::device_vector<float> dxds(xds, xds + DNU);
	thrust::device_vector<float> dyds(yds, yds + DNU);
	thrust::device_vector<float> dzds(zds, zds + DNU);

	thrust::device_vector<float> dbxds = getDD3Boundaries(xds, DNU);
	thrust::device_vector<float> dbyds = getDD3Boundaries(yds, DNU);
	thrust::device_vector<float> dbzds = getDD3Boundaries(zds, DNV);

	thrust::device_vector<float> dprj(hprj, hprj + DNU * DNV * PN);
	thrust::device_vector<float> dvol(hvol, hvol + XN * YN * ZN);

	dim3 blk(64, 4, 1);
	dim3 gid(
		(DNV + blk.x - 1) / blk.x,
		(DNU + blk.y - 1) / blk.y,
		(PN + blk.z - 1) / blk.z);
	float objCntIdxX = (static_cast<float>(XN) - 1.0f) * 0.5f - imgXCenter / dx;
	float objCntIdxY = (static_cast<float>(YN) - 1.0f) * 0.5f - imgYCenter / dx;
	float objCntIdxZ = (static_cast<float>(ZN) - 1.0f) * 0.5f - imgZCenter / dz;

	DD3Proj_branches_ker << <gid, blk >> >(
		thrust::raw_pointer_cast(&dprj[0]),
		thrust::raw_pointer_cast(&dvol[0]),
		x0, y0, z0,
		thrust::raw_pointer_cast(&dxds[0]), thrust::raw_pointer_cast(&dyds[0]),
		thrust::raw_pointer_cast(&dbxds[0]), thrust::raw_pointer_cast(&dbyds[0]), thrust::raw_pointer_cast(&dbzds[0]),
		thrust::raw_pointer_cast(&cossinZT[0]),
		objCntIdxX, objCntIdxY, objCntIdxZ,
		XN, YN, ZN, DNU, DNV, PN,
		dx, dz);
	thrust::copy(dprj.begin(), dprj.end(), hprj);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// END BRANCHES PROJECTION
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
	case 1: // Todo: incorrect
		DD3_gpu_proj_volumerendering(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
			hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz, mask, gpunum);
		break;
	case 2: // Todo: incorrect
		DD3_gpu_proj_doubleprecisionbranchless(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
			hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz, mask, gpunum);
	case 3:
		DD3_gpu_proj_pseudodistancedriven(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
			hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz, mask, gpunum);
		break;
	case 4: // Todo: incorrect
		DD3ProjSiddon_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
				hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz, mask, gpunum);
		break;
	case 5:
		DD3Proj_branches(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
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



void DD3_gpu_proj_pseudodistancedriven_multiGPU(
		float x0, float y0, float z0,
		int DNU, int DNV,
		float* xds, float* yds, float* zds,
		float imgXCenter, float imgYCenter, float imgZCenter,
		float* h_angs, float* h_zPos, int PN,
		int XN, int YN, int ZN,
		float* hvol, float* hprj,
		float dx, float dz,
		byte* mask, int* startPN, int gpuNum)
{
	thrust::host_vector<float> hangs(h_angs, h_angs + PN);
	thrust::host_vector<float> hzPos(h_zPos, h_zPos + PN);
	// Mask the volume
	for (int i = 0; i != XN * YN; ++i)
	{
		byte v = mask[i];
		for (int z = 0; z != ZN; ++z)
		{
			hvol[i * ZN + z] = hvol[i * ZN + z] * v;
		}
	}
	// Calculate the boundary positions

	const float objCntIdxX = (XN - 1.0) * 0.5 - imgXCenter / dx;
	const float objCntIdxY = (YN - 1.0) * 0.5 - imgYCenter / dx;
	const float objCntIdxZ = (ZN - 1.0) * 0.5 - imgZCenter / dz;

	// Divide the volume into sub volumes with overlaps according to the startPN
	std::vector<int> ObjIdx_Start(gpuNum, -1);
	std::vector<int> ObjIdx_End(gpuNum, -1);

	std::vector<int> PrjIdx_Start(startPN, startPN+gpuNum);
	std::vector<int> PrjIdx_End(gpuNum, 0);

	std::copy(PrjIdx_Start.begin()+1, PrjIdx_Start.end(), PrjIdx_End.begin());
	PrjIdx_End[gpuNum - 1] = PN;
	std::vector<int> SPN = PrjIdx_End - PrjIdx_Start;
	std::vector<int> prefixSPN = SPN;

	thrust::exclusive_scan(prefixSPN.begin(), prefixSPN.end(), prefixSPN.begin());
	//std::cout<<"prefixSPN are "<<prefixSPN[0]<<"  "<<prefixSPN[1]<<"  "<<prefixSPN[2]<<"\n";

	std::vector<int> SZN(gpuNum, 0); // The slices number of each sub volume
	const float detStpZ = (zds[DNV-1] - zds[0]) / (DNV - 1); // detector cell height
	assert(detStpZ != 0);
	const float detCntIdxV = -zds[0] / detStpZ; // Detector center along Z direction

	std::vector<std::vector<float> > subVol(gpuNum); // Used to store three sub volumes
	std::vector<float> subImgZCenter(gpuNum, 0); // the center of three sub volumes

	// Generate multiple streams;
	std::vector<cudaStream_t> stream(gpuNum);

	std::vector<int> siz(gpuNum, 0);
	std::vector<cudaExtent> volumeSize(gpuNum);
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	std::vector<cudaArray*> d_volumeArray(gpuNum);

	thrust::host_vector<thrust::device_vector<float> > d_vol(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > prj(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > d_xds(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > d_yds(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > d_zds(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > angs(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > zPos(gpuNum);
	thrust::host_vector<thrust::device_vector<float3> > cossinZT(gpuNum);

	dim3 blk(64, 16, 1);
	std::vector<dim3> gid(gpuNum);
	std::vector<cudaTextureObject_t> texObj(gpuNum);
	// First we define the main framework to see how it works.
	omp_set_num_threads(gpuNum);
#pragma omp parallel for
	for(int i = 0; i < gpuNum; ++i)
	{
		getVolZIdxPair<float>(hzPos, PrjIdx_Start[i], PrjIdx_End[i],
				detCntIdxV, detStpZ, DNV, objCntIdxZ, dz, ZN, ObjIdx_Start[i],
				ObjIdx_End[i]);
		//std::cout<<i<<" "<<ObjIdx_Start[i]<<" "<<ObjIdx_End[i]<<"\n";
		SZN[i] = ObjIdx_End[i] - ObjIdx_Start[i];
		subVol[i].resize(XN * YN * SZN[i]);
		// Divide the volume into multiple sets
		getSubVolume<float>(hvol, XN * YN, ZN, ObjIdx_Start[i], ObjIdx_End[i], &(subVol[i][0]));

		// NOTE: The explanation will be later:
		subImgZCenter[i] = -imgZCenter / dz + ZN * 0.5 - ObjIdx_Start[i] - 0.5f;

		CUDA_CHECK_RETURN(cudaSetDevice(i));
		// For each GPU generate two streams
		CUDA_CHECK_RETURN(cudaStreamCreate(&stream[i]));
		siz[i] = XN * YN * SZN[i];

		d_vol[i].resize(siz[i]);
		d_vol[i] = subVol[i];
		subVol[i].clear();

		volumeSize[i].width = SZN[i];
		volumeSize[i].height = XN;
		volumeSize[i].depth = YN;
		CUDA_CHECK_RETURN(cudaMalloc3DArray(&d_volumeArray[i], &channelDesc, volumeSize[i]));

		cudaMemcpy3DParms copyParams = { 0 };
		copyParams.srcPtr = make_cudaPitchedPtr((void*)
			thrust::raw_pointer_cast(&d_vol[i][0]),
			volumeSize[i].width * sizeof(float),
			volumeSize[i].width, volumeSize[i].height);
		copyParams.dstArray = d_volumeArray[i];
		copyParams.extent = volumeSize[i];
		copyParams.kind = cudaMemcpyDeviceToDevice;

		CUDA_CHECK_RETURN(cudaMemcpy3DAsync(&copyParams,stream[i]));
		d_vol[i].clear();


		cudaResourceDesc resDesc;
		memset(&resDesc, 0, sizeof(resDesc));
		resDesc.resType = cudaResourceTypeArray;
		resDesc.res.array.array = d_volumeArray[i];

		cudaTextureDesc texDesc;
		memset(&texDesc, 0, sizeof(texDesc));
		texDesc.addressMode[0] = cudaAddressModeBorder;
		texDesc.addressMode[1] = cudaAddressModeBorder;
		texDesc.addressMode[2] = cudaAddressModeBorder;
		texDesc.filterMode = cudaFilterModeLinear;
		texDesc.readMode = cudaReadModeElementType;
		texDesc.normalizedCoords = false;
		texObj[i] = 0;
		CUDA_CHECK_RETURN(cudaCreateTextureObject(&texObj[i], &resDesc, &texDesc, nullptr));


		prj[i].resize(DNU * DNV * SPN[i]); // Change here
		d_xds[i].resize(DNU);
		d_yds[i].resize(DNU);
		d_zds[i].resize(DNV);
		thrust::copy(xds,xds+DNU,d_xds[i].begin());
		thrust::copy(yds,yds+DNU,d_yds[i].begin());
		thrust::copy(zds,zds+DNV,d_zds[i].begin());

		angs[i].resize(SPN[i]);
		zPos[i].resize(SPN[i]);
		thrust::copy(hangs.begin() + PrjIdx_Start[i],
					 hangs.begin() + PrjIdx_Start[i] + SPN[i],
					 angs[i].begin());
		thrust::copy(hzPos.begin() + PrjIdx_Start[i],
					 hzPos.begin() + PrjIdx_Start[i] + SPN[i],
					 zPos[i].begin());
		cossinZT[i].resize(PN);

		thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(angs[i].begin(), zPos[i].begin())),
			thrust::make_zip_iterator(thrust::make_tuple(angs[i].end(), zPos[i].end())),
			cossinZT[i].begin(), CTMBIR::ConstantForBackProjection4(x0, y0, z0));
		angs[i].clear();
		zPos[i].clear();

		gid[i].x = (DNV + blk.x - 1) / blk.x;
		gid[i].y = (DNU + blk.y - 1) / blk.y;
		gid[i].z = (SPN[i] + blk.z - 1) / blk.z;
	}
#pragma omp parallel for
	for(int i = 0; i < gpuNum; ++i)
	{
		cudaSetDevice(i);
		DD3_gpu_proj_pseudodistancedriven_ker<< <gid[i], blk, 0, stream[i]>> >(
				texObj[i], thrust::raw_pointer_cast(&prj[i][0]),
			make_float3(x0, y0, z0),
			thrust::raw_pointer_cast(&d_xds[i][0]),
			thrust::raw_pointer_cast(&d_yds[i][0]),
			thrust::raw_pointer_cast(&d_zds[i][0]),
			thrust::raw_pointer_cast(&cossinZT[i][0]),
			make_float3(objCntIdxX, objCntIdxY, subImgZCenter[i]),
			dx, dz, XN, YN, DNU, DNV, SPN[i]);
	}
#pragma omp barrier
#pragma omp parallel for
	for(int i = 0; i < gpuNum; ++i)
	{
		cudaSetDevice(i);
		CUDA_CHECK_RETURN(cudaMemcpyAsync(hprj + DNU * DNV * prefixSPN[i],
				thrust::raw_pointer_cast(&prj[i][0]), sizeof(float) * DNU * DNV * SPN[i],
				cudaMemcpyDeviceToHost,stream[i]));
		d_xds[i].clear();
		d_yds[i].clear();
		d_zds[i].clear();
		cossinZT[i].clear();
		prj[i].clear();

		CUDA_CHECK_RETURN(cudaDestroyTextureObject(texObj[i]));
		CUDA_CHECK_RETURN(cudaFreeArray(d_volumeArray[i]));
	}
}


void DD3_gpu_proj_branchless_sat2d_multiGPU(
		float x0, float y0, float z0,
		int DNU, int DNV,
		float* xds, float* yds, float* zds,
		float imgXCenter, float imgYCenter, float imgZCenter,
		float* h_angs, float* h_zPos, int PN,
		int XN, int YN, int ZN,
		float* hvol, float* hprj,
		float dx, float dz,
		byte* mask, int* startPN, int gpuNum)
{
	thrust::host_vector<float> hangs(h_angs, h_angs+PN);
	thrust::host_vector<float> hzPos(h_zPos, h_zPos+PN);

	maskingVolume(hvol, mask, XN, YN, ZN);
	// Calculate the boundary positions
	std::vector<float> bxds(DNU + 1, 0.0f);
	std::vector<float> byds(DNU + 1, 0.0f);
	std::vector<float> bzds(DNV + 1, 0.0f);

	DD3Boundaries<float>(DNU + 1, xds, &(bxds[0]));
	DD3Boundaries<float>(DNU + 1, yds, &(byds[0]));
	DD3Boundaries<float>(DNV + 1, zds, &(bzds[0]));

	const float objCntIdxX = (XN - 1.0) * 0.5 - imgXCenter / dx;
	const float objCntIdxY = (YN - 1.0) * 0.5 - imgYCenter / dx;
	const float objCntIdxZ = (ZN - 1.0) * 0.5 - imgZCenter / dz;

	// Divide the volume into sub volumes with overlaps according to the startPN
	std::vector<int> ObjIdx_Start(gpuNum, -1);
	std::vector<int> ObjIdx_End(gpuNum, -1);

	std::vector<int> PrjIdx_Start(startPN, startPN+gpuNum);
	std::vector<int> PrjIdx_End(gpuNum, 0);

	std::copy(PrjIdx_Start.begin()+1, PrjIdx_Start.end(), PrjIdx_End.begin());
	PrjIdx_End[gpuNum - 1] = PN;
	std::vector<int> SPN = PrjIdx_End - PrjIdx_Start;
	std::vector<int> prefixSPN = SPN;
	thrust::exclusive_scan(prefixSPN.begin(), prefixSPN.end(), prefixSPN.begin());
	//std::cout<<"prefixSPN are "<<prefixSPN[0]<<"  "<<prefixSPN[1]<<"  "<<prefixSPN[2]<<"\n";

	std::vector<int> SZN(gpuNum, 0); // The slices number of each sub volume
	const float detStpZ = (zds[DNV-1] - zds[0]) / (DNV - 1); // detector cell height
	const float detCntIdxV = -zds[0] / detStpZ; // Detector center along Z direction

	std::vector<std::vector<float> > subVol(gpuNum); // Used to store three sub volumes
	std::vector<float> subImgZCenter(gpuNum, 0); // the center of three sub volumes

	// Generate multiple streams;
	std::vector<cudaStream_t> stream(gpuNum * 2);

	std::vector<int> siz(gpuNum, 0);
	std::vector<int> nsiz_ZXY(gpuNum, 0);
	std::vector<int> nsiz_ZYX(gpuNum, 0);
	std::vector<int> nZN(gpuNum,0);

	const int nXN = XN + 1;
	const int nYN = YN + 1;

	thrust::host_vector<thrust::device_vector<float> > d_vol(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > d_ZXY(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > d_ZYX(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > prj(gpuNum); // Change here
	thrust::host_vector<thrust::device_vector<float> > d_xds(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > d_yds(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > d_zds(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > d_bxds(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > d_byds(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > d_bzds(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > angs(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > zPos(gpuNum);
	thrust::host_vector<thrust::device_vector<float3> > cossinZT(gpuNum);

	// Copy to three volumes
	dim3 copyblk(64, 16, 1);
	std::vector<dim3> copygid(gpuNum);
	dim3 satblk1(32,1,1);
	dim3 satblk2(64,16,1);
	dim3 satgid1_1((nXN * YN + satblk1.x - 1) / satblk1.x, 1, 1);
	dim3 satgid1_2((nYN * XN + satblk1.x - 1) / satblk1.x, 1, 1);
	std::vector<dim3> satgid2_1(gpuNum);
	std::vector<dim3> satgid2_2(gpuNum);

	dim3 blk(BLKX, BLKY, BLKZ);
	std::vector<dim3> gid(gpuNum);

	std::vector<cudaExtent> volumeSize1(gpuNum);
	std::vector<cudaExtent> volumeSize2(gpuNum);

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

	std::vector<cudaArray*> d_volumeArray1(gpuNum);
	std::vector<cudaArray*> d_volumeArray2(gpuNum);

	std::vector<cudaTextureObject_t> texObj1(gpuNum);
	std::vector<cudaTextureObject_t> texObj2(gpuNum);

	omp_set_num_threads(gpuNum);
#pragma omp parallel for
	for(int i = 0; i < gpuNum; ++i)
	{
		getVolZIdxPair<float>(hzPos, PrjIdx_Start[i], PrjIdx_End[i],
				detCntIdxV, detStpZ, DNV, objCntIdxZ, dz, ZN, ObjIdx_Start[i],
				ObjIdx_End[i]);
		SZN[i] = ObjIdx_End[i] - ObjIdx_Start[i];
		subVol[i].resize(XN * YN * SZN[i]);
		// Divide the volume into multiple sets
		getSubVolume<float>(hvol, XN * YN, ZN, ObjIdx_Start[i], ObjIdx_End[i], &(subVol[i][0]));

		// NOTE: How it comes
		// We need to calculate the (ii - subImgZCenter[i]) * dz to define the
		// real physical position of the voxel.
		// Assume that the physical center of the whole volume is imgZCenter
		// The minimum lower position of the volume is imgZCenter - dz * N / 2;
		// Then the corresponding physical lower boundary position of ObjIdx_Start[i]
		// is --> imgZCenter - dz * N / 2 + ObjIdx_Start[i] * dz
		// while the corresponding physical center position of layer ObjIdx_Start[i]
		// is -->  imgZCenter - dz * N / 2 + ObjIdx_Start[i] * dz + 0.5 * dz
		// We need when ii==0 --> (ii - subImgZCenter[i]) * dz = imgZCenter - dz * N / 2 + ObjIdx_Start[i] * dz + 0.5 * dz
		// It means subImgZCenter[i] = -imgZCenter / dz + N / 2 - ObjIdx_Start[i] - 0.5;
		subImgZCenter[i] = -imgZCenter / dz + ZN * 0.5 - ObjIdx_Start[i] - 0.5f;


		CUDA_CHECK_RETURN(cudaSetDevice(i));
		// For each GPU generate two streams
		CUDA_CHECK_RETURN(cudaStreamCreate(&stream[i * 2]));
		CUDA_CHECK_RETURN(cudaStreamCreate(&stream[i * 2 + 1]));
		siz[i] = XN * YN * SZN[i];
		nZN[i] = SZN[i] + 1;
		nsiz_ZXY[i] = nZN[i] * nXN * YN;
		nsiz_ZYX[i] = nZN[i] * nYN * XN;

		d_ZXY[i].resize(nsiz_ZXY[i]);
		d_ZYX[i].resize(nsiz_ZYX[i]);
		d_vol[i].resize(siz[i]);
		d_vol[i] = subVol[i];
		subVol[i].clear();

		copygid[i].x = (SZN[i] + copyblk.x - 1) / copyblk.x;
		copygid[i].y = (XN + copyblk.y - 1) / copyblk.y;
		copygid[i].z = (YN + copyblk.z - 1) / copyblk.z;
		naive_copyToTwoVolumes << <copygid[i], copyblk, 0, stream[2 * i] >> >(
				thrust::raw_pointer_cast(&d_vol[i][0]),
				thrust::raw_pointer_cast(&d_ZXY[i][0]),
				thrust::raw_pointer_cast(&d_ZYX[i][0]),
				XN,YN,SZN[i]);
		CUDA_CHECK_RETURN(cudaStreamSynchronize(stream[2 * i]));
		CUDA_CHECK_RETURN(cudaStreamSynchronize(stream[2 * i + 1]));

		d_vol[i].clear();
		// Generate the SAT for two volumes
		satgid2_1[i].x = (nZN[i] + satblk2.x - 1) / satblk2.x;
		satgid2_1[i].y = (YN + satblk2.y - 1) / satblk2.y;
		satgid2_1[i].z = 1;

		satgid2_2[i].x = (nZN[i] + satblk2.x - 1) / satblk2.x;
		satgid2_2[i].y = (XN + satblk2.y - 1) / satblk2.y;
		satgid2_2[i].z = 1;

		verticalIntegral << <satgid1_1, satblk1, 0, stream[2 * i] >> >(
				thrust::raw_pointer_cast(&d_ZXY[i][0]), nZN[i], nXN * YN);
		horizontalIntegral << <satgid2_1[i], satblk2, 0, stream[2 * i] >> >(
				thrust::raw_pointer_cast(&d_ZXY[i][0]), nXN, nZN[i], YN);
		verticalIntegral << <satgid1_2, satblk1, 0, stream[2 * i + 1] >> >(
				thrust::raw_pointer_cast(&d_ZYX[i][0]), nZN[i], nYN * XN);
		horizontalIntegral << <satgid2_2[i], satblk2, 0, stream[2 * i + 1] >> >(
				thrust::raw_pointer_cast(&d_ZYX[i][0]), nYN, nZN[i], XN);

		//Bind to the texture;
		volumeSize1[i].width = nZN[i];
		volumeSize1[i].height = nXN;
		volumeSize1[i].depth = YN;

		volumeSize2[i].width = nZN[i];
		volumeSize2[i].height = nYN;
		volumeSize2[i].depth = XN;

		CUDA_CHECK_RETURN(cudaMalloc3DArray(&d_volumeArray1[i], &channelDesc, volumeSize1[i]));
		CUDA_CHECK_RETURN(cudaMalloc3DArray(&d_volumeArray2[i], &channelDesc, volumeSize2[i]));

		cudaMemcpy3DParms copyParams1 = { 0 };
		copyParams1.srcPtr = make_cudaPitchedPtr((void*)
			thrust::raw_pointer_cast(&d_ZXY[i][0]),
			volumeSize1[i].width * sizeof(float),
			volumeSize1[i].width, volumeSize1[i].height);
		copyParams1.dstArray = d_volumeArray1[i];
		copyParams1.extent = volumeSize1[i];
		copyParams1.kind = cudaMemcpyDeviceToDevice;

		cudaMemcpy3DParms copyParams2 = { 0 };
		copyParams2.srcPtr = make_cudaPitchedPtr((void*)
			thrust::raw_pointer_cast(&d_ZYX[i][0]),
			volumeSize2[i].width * sizeof(float),
			volumeSize2[i].width, volumeSize2[i].height);
		copyParams2.dstArray = d_volumeArray2[i];
		copyParams2.extent = volumeSize2[i];
		copyParams2.kind = cudaMemcpyDeviceToDevice;

		CUDA_CHECK_RETURN(cudaMemcpy3DAsync(&copyParams1,stream[2 * i]));
		CUDA_CHECK_RETURN(cudaMemcpy3DAsync(&copyParams2,stream[2 * i + 1]));

		d_ZXY[i].clear();
		d_ZYX[i].clear();

		cudaResourceDesc resDesc1;
		cudaResourceDesc resDesc2;
		memset(&resDesc1, 0, sizeof(resDesc1));
		memset(&resDesc2, 0, sizeof(resDesc2));
		resDesc1.resType = cudaResourceTypeArray;
		resDesc2.resType = cudaResourceTypeArray;
		resDesc1.res.array.array = d_volumeArray1[i];
		resDesc2.res.array.array = d_volumeArray2[i];
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
		texObj1[i] = 0;
		texObj2[i] = 0;
		CUDA_CHECK_RETURN(cudaCreateTextureObject(&texObj1[i], &resDesc1, &texDesc1, nullptr));
		CUDA_CHECK_RETURN(cudaCreateTextureObject(&texObj2[i], &resDesc2, &texDesc2, nullptr));


		prj[i].resize(DNU * DNV * SPN[i]); // Change here
		d_xds[i].resize(DNU);
		d_yds[i].resize(DNU);
		d_zds[i].resize(DNV);
		thrust::copy(xds,xds+DNU,d_xds[i].begin());
		thrust::copy(yds,yds+DNU,d_yds[i].begin());
		thrust::copy(zds,zds+DNV,d_zds[i].begin());
		d_bxds[i].resize(bxds.size());
		d_bxds[i] = bxds;
		d_byds[i].resize(byds.size());
		d_byds[i] = byds;
		d_bzds[i].resize(bzds.size());
		d_bzds[i] = bzds;

		angs[i].resize(SPN[i]);
		zPos[i].resize(SPN[i]);
		thrust::copy(hangs.begin() + PrjIdx_Start[i],
					 hangs.begin() + PrjIdx_Start[i] + SPN[i],
					 angs[i].begin());
		thrust::copy(hzPos.begin() + PrjIdx_Start[i],
					 hzPos.begin() + PrjIdx_Start[i] + SPN[i],
					 zPos[i].begin());
		cossinZT[i].resize(PN);

		thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(angs[i].begin(), zPos[i].begin())),
			thrust::make_zip_iterator(thrust::make_tuple(angs[i].end(), zPos[i].end())),
			cossinZT[i].begin(), CTMBIR::ConstantForBackProjection4(x0, y0, z0));
		angs[i].clear();
		zPos[i].clear();

		gid[i].x = (DNV + blk.x - 1) / blk.x;
		gid[i].y = (DNU + blk.y - 1) / blk.y;
		gid[i].z = (SPN[i] + blk.z - 1) / blk.z;

	}
#pragma omp parallel for
	for(int i = 0; i < gpuNum; ++i)
	{
		cudaSetDevice(i);
		DD3_gpu_proj_branchless_sat2d_ker << <gid[i], blk, 0, stream[i * 2]>> >(
				texObj1[i], texObj2[i],
			thrust::raw_pointer_cast(&prj[i][0]),
			make_float3(x0, y0, z0),
			thrust::raw_pointer_cast(&cossinZT[i][0]),
			thrust::raw_pointer_cast(&d_xds[i][0]),
			thrust::raw_pointer_cast(&d_yds[i][0]),
			thrust::raw_pointer_cast(&d_zds[i][0]),
			thrust::raw_pointer_cast(&d_bxds[i][0]),
			thrust::raw_pointer_cast(&d_byds[i][0]),
			thrust::raw_pointer_cast(&d_bzds[i][0]),
			make_float3(objCntIdxX, objCntIdxY, subImgZCenter[i]),
			dx, dz, XN, YN, ZN, DNU, DNV, SPN[i]);
	}
#pragma omp barrier
#pragma omp parallel for
	for(int i = 0; i < gpuNum; ++i)
	{
		cudaSetDevice(i);
		CUDA_CHECK_RETURN(cudaMemcpyAsync(hprj + DNU * DNV * prefixSPN[i],
				thrust::raw_pointer_cast(&prj[i][0]), sizeof(float) * DNU * DNV * SPN[i],
				cudaMemcpyDeviceToHost,stream[2*i]));
		d_xds[i].clear();
		d_yds[i].clear();
		d_zds[i].clear();
		d_bxds[i].clear();
		d_byds[i].clear();
		d_bzds[i].clear();
		cossinZT[i].clear();
		prj[i].clear();

		CUDA_CHECK_RETURN(cudaDestroyTextureObject(texObj1[i]));
		CUDA_CHECK_RETURN(cudaDestroyTextureObject(texObj2[i]));
		CUDA_CHECK_RETURN(cudaFreeArray(d_volumeArray1[i]));
		CUDA_CHECK_RETURN(cudaFreeArray(d_volumeArray2[i]));
	}
}



extern "C"
void DD3Proj_multiGPU(
	float x0, float y0, float z0,
	int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN,
	float* hvol, float* hprj,
	float dx, float dz,
	byte* mask, int prjMode, int* startPN, int gpuNum)
{
	switch(prjMode)
	{
	case 0: // Branchless DD model based multi-GPU projection
		DD3_gpu_proj_branchless_sat2d_multiGPU(x0, y0, z0, DNU, DNV,
				xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
				hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz,
				mask, startPN, gpuNum);
		break;
	case 1: // Volume rendering based multi-GPU projection
	case 2: // Pseudo DD based multi-GPUs projection
	case 3:
		DD3_gpu_proj_pseudodistancedriven_multiGPU(x0, y0, z0, DNU, DNV,
				xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
				hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz,
				mask, startPN, gpuNum);
		break;

	case 4: // Siddon's method multi-GPUs projection
	case 5: // Brute force DD projection in multi-GPUs projection
	default: // Branchless DD model based multi-GPUs projection
		DD3_gpu_proj_branchless_sat2d_multiGPU(x0, y0, z0, DNU, DNV,
				xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
				hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz,
				mask, startPN, gpuNum);
		break;

	}
}

void CT::Proj(std::vector<float>& hvol, std::vector<float>& hprj, Geometry geo, const std::string& projModel)
{
	std::vector<unsigned char> mask(geo.getObjDimX() * geo.getObjDimY(), 1);
	std::vector<float> xds = geo.getXds();
	std::vector<float> yds = geo.getYds();
	std::vector<float> zds = geo.getZds();
	std::vector<float> hangs = geo.getAllAngles();
	std::vector<float> hzPos = geo.getAllZPoses();

	DD3Proj_gpu(0, geo.getSourToObj(), 0,
		geo.getDetNumWidth(), geo.getDetNumHeight(),
		&xds[0], &yds[0], &zds[0],
		geo.getObjCtrCoordX(), geo.getObjCtrCoordY(), geo.getObjCtrCoordZ(),
		&hangs[0], &hzPos[0],
		geo.getViewNumber(),
		geo.getObjDimX(), geo.getObjDimY(), geo.getObjDimZ(),
		&hvol[0], &hprj[0],
		geo.getVoxelSizeX(), geo.getVoxelSizeZ(),
		&mask[0], 0, 0);
}