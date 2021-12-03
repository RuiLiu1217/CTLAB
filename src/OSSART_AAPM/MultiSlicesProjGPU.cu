/**
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 * Author : Rui Liu
 * Date   : Apr. 4, 2016
 */
#include "mex.h"
#include "matrix.h"
     
#include "utilities.cuh"
// CUDA thread block dimension
#define BLKX 32
#define BLKY 8
#define BLKZ 1

namespace DD2 {

	// Copy the volume from the original to
	template<typename Ta, typename Tb>
	__global__ void naive_copyToTwoVolumes(Ta* in_ZXY,
		Tb* out_ZXY, Tb* out_ZYX,
		int XN, int YN, int ZN) {

		int idz = threadIdx.x + blockIdx.x * blockDim.x;
		int idx = threadIdx.y + blockIdx.y * blockDim.y;
		int idy = threadIdx.z + blockIdx.z * blockDim.z;
		if (idx < XN && idy < YN && idz < ZN) {
			int i = (idy * XN + idx) * ZN + idz;
			int ni = (idy * (XN + 1) + (idx + 1)) * ZN + idz;
			int nj = (idx * (YN + 1) + (idy + 1)) * ZN + idz;
			out_ZXY[ni] = in_ZXY[i];
			out_ZYX[nj] = in_ZXY[i];
		}
	}

	__global__ void horizontalIntegral(float* prj, int DNU, int DNV, int PN) {
		int idv = threadIdx.x + blockIdx.x * blockDim.x;
		int pIdx = threadIdx.y + blockIdx.y * blockDim.y;
		if (idv < DNV && pIdx < PN) {
			int headPrt = pIdx * DNU * DNV + idv;
			for (int ii = 1; ii < DNU; ++ii) {
				prj[headPrt + ii * DNV] = prj[headPrt + ii * DNV] + prj[headPrt + (ii - 1) * DNV];
			}
		}
	}

	void genSAT_for_Volume_MultiSlice(float* hvol,
		float* ZXY,
		float* ZYX,
		int XN, int YN, int ZN) {

		const int siz = XN * YN * ZN;

		float* vol = nullptr;
		CUDA_CHECK_RETURN(cudaMalloc(&vol, sizeof(float) * siz));
		CUDA_CHECK_RETURN(cudaMemcpy(vol, hvol, sizeof(float) * siz, cudaMemcpyHostToDevice));

		dim3 blk(64, 16, 1);
		dim3 gid(
			(ZN + blk.x - 1) / blk.x,
			(XN + blk.y - 1) / blk.y,
			(YN + blk.z - 1) / blk.z);

		naive_copyToTwoVolumes << <gid, blk >> > (
			thrust::raw_pointer_cast(&vol[0]),
			thrust::raw_pointer_cast(&ZXY[0]),
			thrust::raw_pointer_cast(&ZYX[0]),
			XN, YN, ZN);

		CUDA_CHECK_RETURN(cudaFree(vol));

		blk.x = 64;
		blk.y = 16;
		blk.z = 1;
		gid.x = (ZN + blk.x - 1) / blk.x;
		gid.y = (YN + blk.y - 1) / blk.y;
		gid.z = 1;

		horizontalIntegral << <gid, blk >> > (
			ZXY, XN + 1, ZN, YN);

		blk.x = 64;
		blk.y = 16;
		blk.z = 1;
		gid.x = (ZN + blk.x - 1) / blk.x;
		gid.y = (XN + blk.y - 1) / blk.y;
		gid.z = 1;

		horizontalIntegral << <gid, blk >> > (
			ZYX, YN + 1, ZN, XN);
	}

	// Main kernel function for projection
	__global__  void DD2_gpu_proj_branchless_sat2d_ker(
		cudaTextureObject_t volTex1,
		cudaTextureObject_t volTex2,
		float* proj,
		float2 s, // source position
		const float2* __restrict__ cossin,
		const float* __restrict__ xds,
		const float* __restrict__ yds,
		const float* __restrict__ bxds,
		const float* __restrict__ byds,
		float2 objCntIdx,
		float dx,
		int XN, int YN, int SLN,
		int DNU, int PN) {

		int slnIdx = threadIdx.x + blockIdx.x * blockDim.x;
		int detIdU = threadIdx.y + blockIdx.y * blockDim.y;
		int angIdx = threadIdx.z + blockIdx.z * blockDim.z;
		if (slnIdx < SLN && detIdU < DNU && angIdx < PN) {

			float2 dir = cossin[angIdx * SLN + slnIdx]; // cossin;

			float2 cursour = make_float2(
				s.x * dir.x - s.y * dir.y,
				s.x * dir.y + s.y * dir.x); // current source position;
			s = dir;

			float2 curDet = make_float2(
				xds[detIdU] * s.x - yds[detIdU] * s.y,
				xds[detIdU] * s.y + yds[detIdU] * s.x);

			float2 curDetL = make_float2(
				bxds[detIdU] * s.x - byds[detIdU] * s.y,
				bxds[detIdU] * s.y + byds[detIdU] * s.x);

			float2 curDetR = make_float2(
				bxds[detIdU + 1] * s.x - byds[detIdU + 1] * s.y,
				bxds[detIdU + 1] * s.y + byds[detIdU + 1] * s.x);

			dir = normalize(curDet - cursour);

			float factL = 0;
			float factR = 0;
			float constVal = 0;
			float obj = 0;
			float realL = 0;
			float realR = 0;
			float intersectLength = 0;

			float invdx = 1.0f / dx;
			//float summ[BLKX];
			float summ;
			if (fabsf(s.x) <= fabsf(s.y)) {

				summ = 0;
				factL = (curDetL.y - cursour.y) / (curDetL.x - cursour.x);
				factR = (curDetR.y - cursour.y) / (curDetR.x - cursour.x);

				constVal = dx / fabsf(dir.x);
#pragma unroll
				for (int ii = 0; ii < XN; ++ii) {
					obj = (ii - objCntIdx.x) * dx;

					realL = (obj - curDetL.x) * factL + curDetL.y;
					realR = (obj - curDetR.x) * factR + curDetR.y;

					intersectLength = realR - realL;
					realL = realL * invdx + objCntIdx.y + 1;
					realR = realR * invdx + objCntIdx.y + 1;

					summ += (tex3D<float>(volTex2, slnIdx + 0.5f, realR, ii + 0.5) - tex3D<float>(volTex2, slnIdx + 0.5, realL, ii + 0.5)) / intersectLength;
				}
				__syncthreads();
				proj[(angIdx * DNU + detIdU) * SLN + slnIdx] = summ * constVal;

			}
			else {
				summ = 0;
				factL = (curDetL.x - cursour.x) / (curDetL.y - cursour.y);
				factR = (curDetR.x - cursour.x) / (curDetR.y - cursour.y);

				constVal = dx / fabsf(dir.y);
#pragma unroll
				for (int ii = 0; ii < YN; ++ii) {
					obj = (ii - objCntIdx.y) * dx;

					realL = (obj - curDetL.y) * factL + curDetL.x;
					realR = (obj - curDetR.y) * factR + curDetR.x;

					intersectLength = realR - realL;
					realL = realL * invdx + objCntIdx.x + 1;
					realR = realR * invdx + objCntIdx.x + 1;

					summ += (tex3D<float>(volTex1, slnIdx + 0.5f, realR, ii + 0.5) - tex3D<float>(volTex1, slnIdx + 0.5, realL, ii + 0.5)) / intersectLength;
				}
				__syncthreads();
				proj[(angIdx * DNU + detIdU) * SLN + slnIdx] = summ * constVal;
				//__syncthreads();
			}
		}
	}

	static void maskingData(const byte* mask, float* vol, const int XN, const int YN, const int SLN) {
		for (int i = 0; i != XN * YN; ++i) {
			byte v = mask[i];
			for (int z = 0; z != SLN; ++z) {
				vol[i * SLN + z] = vol[i * SLN + z] * v;
			}
		}
	}

	void DD2_gpu_proj_branchless_sat2d(float x0, float y0, int DNU,	float* xds, float* yds,
		float imgXCenter, float imgYCenter, float* hangs, int PN, int XN, int YN, int SLN, // SLN is the slice number, it is the same as the rebinned projection slices
		float* vol, float* hprj, float dx, byte* mask, int gpunum) {
		// Masking the data
		maskingData(mask, vol, XN, YN, SLN);
		CUDA_CHECK_RETURN(cudaSetDevice(gpunum));
		CUDA_CHECK_RETURN(cudaDeviceReset());

		float* bxds = new float[DNU + 1];
		float* byds = new float[DNU + 1];

		DD3Boundaries(DNU + 1, xds, bxds);
		DD3Boundaries(DNU + 1, yds, byds);

		float objCntIdxX = (XN - 1.0f) * 0.5f - imgXCenter / dx;
		float objCntIdxY = (YN - 1.0f) * 0.5f - imgYCenter / dx;

		float* SATZXY = nullptr;
		float* SATZYX = nullptr;
		const int nsiz_ZXY = SLN * (XN + 1) * YN; //Only XN or YN dimension changes
		const int nsiz_ZYX = SLN * (YN + 1) * XN;
		CUDA_CHECK_RETURN(cudaMalloc(&SATZXY, sizeof(float) * nsiz_ZXY));
		CUDA_CHECK_RETURN(cudaMalloc(&SATZYX, sizeof(float) * nsiz_ZYX));
		genSAT_for_Volume_MultiSlice(vol, SATZXY, SATZYX, XN, YN, SLN);

		cudaExtent volumeSize1;
		cudaExtent volumeSize2;
		volumeSize1.width = SLN;
		volumeSize1.height = XN + 1;
		volumeSize1.depth = YN;

		volumeSize2.width = SLN;
		volumeSize2.height = YN + 1;
		volumeSize2.depth = XN;

		cudaChannelFormatDesc channelDesc1 = cudaCreateChannelDesc<float>();
		cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<float>();

		cudaArray* d_volumeArray1;
		cudaArray* d_volumeArray2;

		CUDA_CHECK_RETURN(cudaMalloc3DArray(&d_volumeArray1, &channelDesc1, volumeSize1));
		CUDA_CHECK_RETURN(cudaMalloc3DArray(&d_volumeArray2, &channelDesc2, volumeSize2));

		cudaMemcpy3DParms copyParams1 = { 0 };
		copyParams1.srcPtr = make_cudaPitchedPtr((void*)SATZXY,
			volumeSize1.width * sizeof(float),
			volumeSize1.width, volumeSize1.height);
		copyParams1.dstArray = d_volumeArray1;
		copyParams1.extent = volumeSize1;
		copyParams1.kind = cudaMemcpyDeviceToDevice;

		cudaMemcpy3DParms copyParams2 = { 0 };
		copyParams2.srcPtr = make_cudaPitchedPtr((void*)SATZYX,
			volumeSize2.width * sizeof(float),
			volumeSize2.width, volumeSize2.height);
		copyParams2.dstArray = d_volumeArray2;
		copyParams2.extent = volumeSize2;
		copyParams2.kind = cudaMemcpyDeviceToDevice;

		CUDA_CHECK_RETURN(cudaMemcpy3D(&copyParams1));
		CUDA_CHECK_RETURN(cudaMemcpy3D(&copyParams2));

		CUDA_CHECK_RETURN(cudaFree(SATZXY));
		CUDA_CHECK_RETURN(cudaFree(SATZYX));

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

		CUDA_CHECK_RETURN(cudaCreateTextureObject(&texObj1, &resDesc1, &texDesc1, nullptr));
		CUDA_CHECK_RETURN(cudaCreateTextureObject(&texObj2, &resDesc2, &texDesc2, nullptr));

		float* prj = nullptr;
		CUDA_CHECK_RETURN(cudaMalloc(&prj, sizeof(float)* DNU* SLN* PN));
		CUDA_CHECK_RETURN(cudaMemset(prj, 0, sizeof(float)* DNU* SLN* PN));

		float* d_xds = nullptr;
		CUDA_CHECK_RETURN(cudaMalloc(&d_xds, sizeof(float)* DNU));
		CUDA_CHECK_RETURN(cudaMemcpy(d_xds, xds, sizeof(float)* DNU, cudaMemcpyHostToDevice));

		float* d_yds = nullptr;
		CUDA_CHECK_RETURN(cudaMalloc(&d_yds, sizeof(float)* DNU));
		CUDA_CHECK_RETURN(cudaMemcpy(d_yds, yds, sizeof(float)* DNU, cudaMemcpyHostToDevice));

		float* d_bxds = nullptr;
		CUDA_CHECK_RETURN(cudaMalloc(&d_bxds, sizeof(float)* (DNU + 1)));
		CUDA_CHECK_RETURN(cudaMemcpy(d_bxds, bxds, sizeof(float)* (DNU + 1), cudaMemcpyHostToDevice));

		float* d_byds = nullptr;
		CUDA_CHECK_RETURN(cudaMalloc(&d_byds, sizeof(float)* (DNU + 1)));
		CUDA_CHECK_RETURN(cudaMemcpy(d_byds, byds, sizeof(float)* (DNU + 1), cudaMemcpyHostToDevice));

		float* angs = nullptr;
		CUDA_CHECK_RETURN(cudaMalloc(&angs, sizeof(float)* PN * SLN));
		CUDA_CHECK_RETURN(cudaMemcpy(angs, hangs, sizeof(float)* PN* SLN, cudaMemcpyHostToDevice));

		float2* hcossin = new float2[PN * SLN];
		std::transform(hangs, hangs + SLN * PN, hcossin, CTMBIR::Constant_MultiSlice_host(x0, y0));

		float2* cossin = nullptr;
		CUDA_CHECK_RETURN(cudaMalloc(&cossin, sizeof(float2)* PN* SLN));
		CUDA_CHECK_RETURN(cudaMemcpy(cossin, hcossin, sizeof(float2)* PN* SLN, cudaMemcpyHostToDevice));



		dim3 blk(BLKX, BLKY, BLKZ);
		dim3 gid(
			(SLN + blk.x - 1) / blk.x,
			(DNU + blk.y - 1) / blk.y,
			(PN + blk.z - 1) / blk.z);

		DD2_gpu_proj_branchless_sat2d_ker << <gid, blk >> > (texObj1, texObj2,
			prj,
			make_float2(x0, y0),
			cossin,
			d_xds,
			d_yds,
			d_bxds,
			d_byds,
			make_float2(objCntIdxX, objCntIdxY), dx, XN, YN, SLN, DNU, PN);

		CUDA_CHECK_RETURN(cudaMemcpy(hprj, prj, sizeof(float)* DNU* SLN* PN, cudaMemcpyDeviceToHost));

		CUDA_CHECK_RETURN(cudaDestroyTextureObject(texObj1));
		CUDA_CHECK_RETURN(cudaDestroyTextureObject(texObj2));

		CUDA_CHECK_RETURN(cudaFreeArray(d_volumeArray1));
		CUDA_CHECK_RETURN(cudaFreeArray(d_volumeArray2));

		CUDA_CHECK_RETURN(cudaFree(prj));
		CUDA_CHECK_RETURN(cudaFree(angs));
		CUDA_CHECK_RETURN(cudaFree(d_xds));
		CUDA_CHECK_RETURN(cudaFree(d_yds));
		CUDA_CHECK_RETURN(cudaFree(d_bxds));
		CUDA_CHECK_RETURN(cudaFree(d_byds));
		CUDA_CHECK_RETURN(cudaFree(cossin));

		delete[] bxds;
		delete[] byds;
		delete[] hcossin;
	}

}; //End NAMESPACE DD2


void DD2ProjMultiGPUs(
	float x0, float y0,
	int DNU,
	float* xds, float* yds,
	float imgXCenter, float imgYCenter,
	float* hangs, int PN,
	int XN, int YN, int SLN,
	float* hvol, float* hprj,
	float dx, byte * mask, int* startIdx) {
	int nDevices = 0;
	CUDA_CHECK_RETURN(cudaGetDeviceCount(&nDevices));
	const auto processor_count = std::thread::hardware_concurrency();

	// Divide the volume into three parts
	int* SLNn = new int[nDevices];
	float** shvol = new float* [nDevices];
	float** shprj = new float* [nDevices];
	float** shang = new float* [nDevices];
	for (int i = 0; i < nDevices; ++i) {
		if (i == nDevices - 1) {
			SLNn[i] = SLN - startIdx[i];
		}
		else {
			SLNn[i] = startIdx[i + 1] - startIdx[i];
		}
		shvol[i] = new float[XN * YN * SLNn[i]];
		shprj[i] = new float[DNU * SLNn[i] * PN];
		shang[i] = new float[PN * SLNn[i]];
	}

	for (int pIdx = 0; pIdx != PN; ++pIdx) {
		for (int sIdx = 0; sIdx != SLN; ++sIdx) {
			for (int i = 0; i < nDevices; ++i) {
				if (i == nDevices - 1) {
					if (sIdx >= startIdx[i] && sIdx < SLN) {
						shang[i][pIdx * SLNn[i] + (sIdx - startIdx[i])] = hangs[pIdx * SLN + sIdx];
					}
				}
				else {
					if (sIdx >= startIdx[i] && sIdx < startIdx[i + 1]) {
						shang[i][pIdx * SLNn[i] + (sIdx - startIdx[i])] = hangs[pIdx * SLN + sIdx];
					}
				}
			}
		}
	}

	omp_set_num_threads(processor_count);
#pragma omp parallel for
	for (int i = 0; i < XN * YN; ++i) {
		for (int v = 0; v != SLN; ++v) {
			for (int k = 0; k < nDevices; ++k) {
				if (k < nDevices - 1) {
					if (v >= startIdx[k] && v < startIdx[k + 1]) {
						shvol[k][i * SLNn[k] + (v - startIdx[k])] = hvol[i * SLN + v];
					}
				}
				else {
					if (v >= startIdx[k] && v < SLN) {
						shvol[k][i * SLNn[k] + (v - startIdx[k])] = hvol[i * SLN + v];
					}
				}
			}
		}
	}
	omp_set_num_threads(nDevices);
#pragma omp parallel for
	for (int i = 0; i < nDevices; ++i) {
		DD2::DD2_gpu_proj_branchless_sat2d(x0, y0, DNU, xds, yds,
			imgXCenter, imgYCenter, shang[i],
			PN, XN, YN, SLNn[i], shvol[i], shprj[i], dx, mask, i);
	}

	//Gather all
	omp_set_num_threads(processor_count);
#pragma omp parallel for
	for (int i = 0; i < DNU * PN; ++i) {
		for (int v = 0; v != SLN; ++v) {
			float val = 0;
			for (int k = 0; k < nDevices; ++k) {
				if (k < nDevices - 1) {
					if (v >= startIdx[k] && v < startIdx[k + 1]) {
						val = shprj[k][i * SLNn[k] + (v - startIdx[k])];
					}
				}
				else {
					if (v >= startIdx[k] && v < SLN) {
						val = shprj[k][i * SLNn[k] + (v - startIdx[k])];
					}
				}
			}
			hprj[i * SLN + v] = val;
		}
	}

	for (int i = 0; i < nDevices; ++i) {
		delete[] shprj[i];
		delete[] shvol[i];
	}
	delete[] shprj;
	delete[] shvol;
	delete[] SLNn;
}


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    float x0 = *((float*)mxGetData(prhs[0]));
    float y0 = *((float*)mxGetData(prhs[1]));
    int DNU = *((int*)mxGetData(prhs[2]));
    float* xds = (float*)mxGetPr(prhs[3]);
    float* yds = (float*)mxGetPr(prhs[4]);
    float imgXCenter = *((float*)mxGetData(prhs[5]));
    float imgYCenter = *((float*)mxGetData(prhs[6]));
    float* hangs = (float*)mxGetPr(prhs[7]);
    int PN = *((int*)mxGetData(prhs[8]));
    int XN = *((int*)mxGetData(prhs[9]));
    int YN = *((int*)mxGetData(prhs[10]));
    int SLN = *((int*)mxGetData(prhs[11]));
    float* hvol = (float*)mxGetPr(prhs[12]);
    float dx = *((float*)mxGetData(prhs[13]));
    byte* mask = (byte*)mxGetPr(prhs[14]);
    int* startIdx = (int*)mxGetPr(prhs[15]);
       
    const mwSize dims[] = {SLN, DNU, PN};
    plhs[0] = mxCreateNumericArray(3,dims,mxSINGLE_CLASS,mxREAL);
    float* hprj = (float*)mxGetPr(plhs[0]);
    
    DD2ProjMultiGPUs(x0, y0, DNU, xds, yds, imgXCenter, imgYCenter,
        hangs, PN, XN, YN, SLN, hvol, hprj, dx, mask, startIdx);    
}