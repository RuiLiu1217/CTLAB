#include "mex.h"
#include "matrix.h"
#include <fstream>
#include "utilities.cuh"

#define BACK_BLKX 64
#define BACK_BLKY 4
#define BACK_BLKZ 1

enum BackProjectionMethod { _BRANCHLESS, _PSEUDODD, _ZLINEBRANCHLESS, _VOLUMERENDERING };

__global__ void addOneSidedZeroBoarder_multiSlice_Fan(const float* prj_in, float* prj_out, int DNU, int SLN, int PN) {
	int idv = threadIdx.x + blockIdx.x * blockDim.x;
	int idu = threadIdx.y + blockIdx.y * blockDim.y;
	int pn = threadIdx.z + blockIdx.z * blockDim.z;
	if (idu < DNU && idv < SLN && pn < PN) {
		int i = (pn * DNU + idu) * SLN + idv;
		int ni = (pn * (DNU + 1) + (idu + 1)) * SLN + idv;
		prj_out[ni] = prj_in[i];
	}
}

__global__ void heorizontalIntegral_multiSlice_Fan(float* prj, int DNU, int SLN, int PN) {
	int idv = threadIdx.x + blockIdx.x * blockDim.x;
	int pIdx = threadIdx.y + blockIdx.y * blockDim.y;
	if (idv < SLN && pIdx < PN) {
		int headPrt = pIdx * DNU * SLN + idv;
		for (int ii = 1; ii < DNU; ++ii) {
			prj[headPrt + ii * SLN] = prj[headPrt + ii * SLN] + prj[headPrt + (ii - 1) * SLN];
		}
	}
}

float* genSAT_of_Projection_multiSliceFan_dev(
	float* hprj,
	int DNU, int SLN, int PN) {
	const int siz = DNU * SLN * PN;
	const int nsiz = (DNU + 1) * SLN * PN;
	float* prj = nullptr;
	CUDA_CHECK_RETURN(cudaMalloc((void**)(&prj), sizeof(float) * siz));
	CUDA_CHECK_RETURN(cudaMemcpy((void*)prj, (void*)hprj, sizeof(float) * siz, cudaMemcpyHostToDevice));

	dim3 copyBlk(64, 16, 1); //MAY CHANGED
	dim3 copyGid(
		(SLN + copyBlk.x - 1) / copyBlk.x,
		(DNU + copyBlk.y - 1) / copyBlk.y,
		(PN + copyBlk.z - 1) / copyBlk.z);

	float* prjSAT = nullptr;
	CUDA_CHECK_RETURN(cudaMalloc(&prjSAT, sizeof(float) * nsiz));
	addOneSidedZeroBoarder_multiSlice_Fan << <copyGid, copyBlk >> > (
		prj, prjSAT, DNU, SLN, PN);

	const int nDNU = DNU + 1;

	copyBlk.x = 64; // MAY CHANGED
	copyBlk.y = 16;
	copyBlk.z = 1;
	copyGid.x = (SLN + copyBlk.x - 1) / copyBlk.x;
	copyGid.y = (PN + copyBlk.y - 1) / copyBlk.y;
	copyGid.z = 1;

	heorizontalIntegral_multiSlice_Fan << <copyGid, copyBlk >> > (
		prjSAT, nDNU, SLN, PN);
	return prjSAT;
}


void createTextureObject_multiSliceFan(
	cudaTextureObject_t& texObj,
	cudaArray* d_prjArray,
	int Width, int Height, int Depth,
	float* sourceData,
	cudaMemcpyKind memcpyKind,
	cudaTextureAddressMode addressMode,
	cudaTextureFilterMode textureFilterMode,
	cudaTextureReadMode textureReadMode,
	bool isNormalized) {
	cudaExtent prjSize;
	prjSize.width = Width;
	prjSize.height = Height;
	prjSize.depth = Depth;
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

	CUDA_CHECK_RETURN(cudaMalloc3DArray(&d_prjArray, &channelDesc, prjSize));
	cudaMemcpy3DParms copyParams = { 0 };
	copyParams.srcPtr = make_cudaPitchedPtr(
		(void*)sourceData, prjSize.width * sizeof(float),
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

	cudaCreateTextureObject(&texObj, &resDesc, &texDesc, nullptr);
}

void destroyTextureObject_multiSliceFan(cudaTextureObject_t& texObj, cudaArray* d_array) {
	CUDA_CHECK_RETURN(cudaDestroyTextureObject(texObj));
	CUDA_CHECK_RETURN(cudaFreeArray(d_array));
}

__global__ void DD2_gpu_back_ker_multiSlice_Fan(
	cudaTextureObject_t prjTexObj,
	float* vol,
	const byte* __restrict__ msk,
	const float2* __restrict__ cossin,
	float2 s,
	float S2D,
	float2 curvox, // imgCenter index
	float dx, float dbeta, float detCntIdx,
	int2 VN, int SLN, int PN) {

	int3 id;
	id.z = threadIdx.x + __umul24(blockIdx.x, blockDim.x);
	id.x = threadIdx.y + __umul24(blockIdx.y, blockDim.y);
	id.y = threadIdx.z + __umul24(blockIdx.z, blockDim.z);
	if (id.z < SLN && id.x < VN.x && id.y < VN.y) {
		if (msk[id.y * VN.x + id.x] == 0) {
			return;
		}
		curvox = make_float2((id.x - curvox.x) * dx, (id.y - curvox.y) * dx);
		float2 cursour;
		float idxL, idxR;
		float cosVal;
		float summ = 0;

		float2 cossinT;
		float inv_sid = 1.0f / sqrtf(s.x * s.x + s.y * s.y);

		float2 dir;
		float l_square;
		float l;

		float alpha;
		float deltaAlpha;
		//S2D /= ddv;
		dbeta = 1.0 / dbeta;

		float ddv;
		for (int angIdx = 0; angIdx < PN; ++angIdx) {
			cossinT = cossin[angIdx * SLN + id.z];
			cursour = make_float2(
				s.x * cossinT.x - s.y * cossinT.y,
				s.x * cossinT.y + s.y * cossinT.x);

			dir = curvox - cursour;

			l_square = dir.x * dir.x + dir.y * dir.y;

			l = rsqrtf(l_square); // 1 / sqrt(l_square);
			alpha = asinf((cursour.y * dir.x - cursour.x * dir.y) * inv_sid * l);

			if (fabsf(cursour.x) > fabsf(cursour.y)) {
				ddv = dir.x;
			} else {
				ddv = dir.y;
			}

			deltaAlpha = ddv / l_square * dx * 0.5;
			cosVal = dx / ddv * sqrtf(l_square);

			idxL = (alpha - deltaAlpha) * dbeta + detCntIdx + 1.0;
			idxR = (alpha + deltaAlpha) * dbeta + detCntIdx + 1.0;

			summ += (tex3D<float>(prjTexObj, id.z + 0.5, idxR, angIdx + 0.5) -
				tex3D<float>(prjTexObj, id.z + 0.5, idxL, angIdx + 0.5)) * cosVal;
		}
		__syncthreads();
		vol[(id.y * VN.x + id.x) * SLN + id.z] = summ;
	}
}


void DD2_gpu_back(float x0, float y0,
	int DNU,
	float* xds, float* yds,
	float detCntIdx,
	float imgXCenter, float imgYCenter,
	float* hangs, int PN, int XN, int YN, int SLN,
	float* hvol, float* hprj, float dx,
	byte* mask, int gpunum) {
	CUDA_CHECK_RETURN(cudaSetDevice(gpunum));
	CUDA_CHECK_RETURN(cudaDeviceReset());

	float2 objCntIdx = make_float2(
		(XN - 1.0) * 0.5 - imgXCenter / dx,
		(YN - 1.0) * 0.5 - imgYCenter / dx); // set the center of the image
	float2 sour = make_float2(x0, y0);

	byte* msk = nullptr;
	CUDA_CHECK_RETURN(cudaMalloc(&msk, sizeof(byte) * XN * YN));
	CUDA_CHECK_RETURN(cudaMemcpy(msk, mask, sizeof(byte) * XN * YN, cudaMemcpyHostToDevice));
	
	float* vol = nullptr;
	CUDA_CHECK_RETURN(cudaMalloc(&vol, sizeof(float) * XN * YN * SLN));
	CUDA_CHECK_RETURN(cudaMemcpy(vol, hvol, sizeof(float) * XN * YN * SLN, cudaMemcpyHostToDevice));
	
	const float S2D = hypotf(xds[0] - x0, yds[0] - y0);
	
	float2* hcossin = new float2[PN * SLN];
	std::transform(hangs, hangs + PN * SLN, hcossin, CTMBIR::Constant_MultiSlice_host(x0, y0));

	float2* cossin = nullptr;
	CUDA_CHECK_RETURN(cudaMalloc(&cossin, sizeof(float2) * PN * SLN));
	CUDA_CHECK_RETURN(cudaMemcpy(cossin, hcossin, sizeof(float2) * PN * SLN, cudaMemcpyHostToDevice));

	float* angs = nullptr;
	CUDA_CHECK_RETURN(cudaMalloc(&angs, sizeof(float) * PN * SLN));
	CUDA_CHECK_RETURN(cudaMemcpy(angs, hangs, sizeof(float) * PN * SLN, cudaMemcpyHostToDevice));

	//Calculate the corresponding parameters such as
	// return make_float4(detCtrIdxU, detCtrIdxV, dbeta, ddv);
	float detTransverseSize = sqrt(powf(xds[1] - xds[0], 2) + powf(yds[1] - yds[0], 2));
	float dbeta = atanf(detTransverseSize / S2D * 0.5) * 2.0f;

	cudaArray* d_prjArray = nullptr;
	cudaTextureObject_t texObj;

	dim3 blk;
	dim3 gid;

	//Generate the SAT along XY direction;
	float* prjSAT = genSAT_of_Projection_multiSliceFan_dev(hprj, DNU, SLN, PN);

	createTextureObject_multiSliceFan(texObj, d_prjArray, SLN, DNU + 1, PN,
		prjSAT,
		cudaMemcpyDeviceToDevice,
		cudaAddressModeClamp, cudaFilterModeLinear, cudaReadModeElementType, false);
	CUDA_CHECK_RETURN(cudaFree(prjSAT));

	blk.x = BACK_BLKX; // May be changed
	blk.y = BACK_BLKY;
	blk.z = BACK_BLKZ;
	gid.x = (SLN + blk.x - 1) / blk.x;
	gid.y = (XN + blk.y - 1) / blk.y;
	gid.z = (YN + blk.z - 1) / blk.z;

	DD2_gpu_back_ker_multiSlice_Fan << <gid, blk >> > (texObj,
		vol,
		msk,
		cossin,
		make_float2(x0, y0),
		S2D,
		make_float2(objCntIdx.x, objCntIdx.y),
		dx, dbeta, detCntIdx, make_int2(XN, YN), SLN, PN);

	CUDA_CHECK_RETURN(cudaMemcpy(hvol, vol, sizeof(float) * XN * YN * SLN, cudaMemcpyDeviceToHost));

	destroyTextureObject_multiSliceFan(texObj, d_prjArray);

	CUDA_CHECK_RETURN(cudaFree(vol));
	CUDA_CHECK_RETURN(cudaFree(msk));
	CUDA_CHECK_RETURN(cudaFree(angs));
	CUDA_CHECK_RETURN(cudaFree(cossin));
}


void DD2BackMultiGPUs(
	float x0, float y0,
	int DNU,
	float* xds, float* yds,
	float detCntIdx,
	float imgXCenter, float imgYCenter,
	float* hangs, int PN,
	int XN, int YN, int SLN,
	float* hvol, float* hprj,
	float dx,
	byte* mask, int* startIdx) {
	int nDevices;
	CUDA_CHECK_RETURN(cudaGetDeviceCount(&nDevices));
	const auto processorCount = std::thread::hardware_concurrency();

	// Divide the volume into three parts
	int* SLNn = new int[nDevices];
	float** shvol = new float* [nDevices];
	float** shprj = new float* [nDevices];
	float** shang = new float* [nDevices];
	for (int n = 0; n < nDevices; ++n) {
		if (n == nDevices - 1) {
			SLNn[n] = SLN - startIdx[n];
		}
		else {
			SLNn[n] = startIdx[n + 1] - startIdx[n];
		}
		shvol[n] = new float[XN * YN * SLNn[n]];
		shprj[n] = new float[DNU * SLNn[n] * PN];
		shang[n] = new float[PN * SLNn[n]];
	}


	for (int pIdx = 0; pIdx != PN; ++pIdx) {
		for (int sIdx = 0; sIdx != SLN; ++sIdx) {
			for (int n = 0; n < nDevices; ++n) {
				if (((n == nDevices - 1) && (sIdx >= startIdx[n] && sIdx < SLN)) || (sIdx >= startIdx[n] && sIdx < startIdx[n + 1])){
					shang[n][pIdx * SLNn[n] + (sIdx - startIdx[n])] = hangs[pIdx * SLN + sIdx];
				}
			}
		}
	}

	//Split the projection
	omp_set_num_threads(processorCount);
#pragma omp parallel for
	for (int i = 0; i < DNU * PN; ++i) {
		for (int v = 0; v != SLN; ++v) {
			for (int n = 0; n < nDevices; ++n) {
				if (((n == nDevices - 1) && (v >= startIdx[n] && v < SLN)) || (v >= startIdx[n] && v < startIdx[n + 1])) {
					shprj[n][i * SLNn[n] + (v - startIdx[n])] = hprj[i * SLN + v];
				}
			}
		}
	}

	cudaDeviceSynchronize();
	omp_set_num_threads(nDevices);
#pragma omp parallel for
	for (int i = 0; i < nDevices; ++i) {
		DD2_gpu_back(x0, y0, DNU, xds, yds, detCntIdx,
			imgXCenter, imgYCenter, shang[i],
			PN, XN, YN, SLNn[i], shvol[i], shprj[i], dx, mask, i);
	}

	//Gather all
	omp_set_num_threads(processorCount);
#pragma omp parallel for
	for (int i = 0; i < XN * YN; ++i) {
		for (int v = 0; v != SLN; ++v) {
			float val = 0;

			for (int n = 0; n < nDevices; ++n) {
				if (((n == nDevices - 1) && (v >= startIdx[n] && v < SLN)) || (v >= startIdx[n] && v < startIdx[n + 1])){
					val = shvol[n][i * SLNn[n] + (v - startIdx[n])];
				}
			}
			hvol[i * SLN + v] = val;
		}
	}

	for (int n = 0; n < nDevices; ++n) {
		delete[] shprj[n];
		delete[] shvol[n];
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
    float detCntIdx = *((float*)mxGetData(prhs[5]));
    float imgXCenter = *((float*)mxGetData(prhs[6]));
    float imgYCenter = *((float*)mxGetData(prhs[7]));
    float* hangs = (float*)mxGetPr(prhs[8]);
    int PN = *((int*)mxGetData(prhs[9]));
    int XN = *((int*)mxGetData(prhs[10]));
    int YN = *((int*)mxGetData(prhs[11]));
    int SLN = *((int*)mxGetData(prhs[12]));
    float* hprj = (float*)mxGetPr(prhs[13]);
    float dx = *((float*)mxGetData(prhs[14]));
    byte* mask = (byte*)mxGetPr(prhs[15]);
    int* startIdx = (int*)mxGetPr(prhs[16]);
       
    const mwSize dims[] = {SLN, XN, YN};
    plhs[0] = mxCreateNumericArray(3,dims,mxSINGLE_CLASS,mxREAL);
    float* hvol = (float*)mxGetPr(plhs[0]);
    
    DD2BackMultiGPUs(x0, y0, DNU, xds, yds, detCntIdx,
        imgXCenter, imgYCenter, hangs, PN, XN, YN, SLN,
        hvol, hprj, dx, mask, startIdx); 
}