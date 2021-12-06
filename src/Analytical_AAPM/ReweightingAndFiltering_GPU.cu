/**
 * Wake Forest University Health Science
 * University of Massachusetts Lowell
 *
 * Organization:
 * 	Wake Forest University
 *
 * 	reWeiAdFiltr.cpp
 * 	Matlab mex routine for the GPU based reweighting
 * 	and filtering for analytical helical CT
 * 	reconstruction
 *
 * 	Author: Rui Liu
 * 	Email: liurui1217@gmail.com
 * 	Date: 2016-08-04
 *
 * 	Version 1.0
 */

#include "mex.h"
#include "matrix.h"

#include <iostream>

#include <vector>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <cufft.h>
// includes, project
#include <cuda_runtime.h>

#include <iostream>
typedef float2 Complex;


struct genKernel
{
	float _dYL;
	genKernel(const float _dyl):_dYL(_dyl){}

	__host__ __device__ float operator()(float n) const
	{
		return -2.0 / (9.869604401089358 * (4.0 * n * n - 1.0)) / _dYL;
	}
};

// Generate the Hilbert Filtering Kernel

__global__ void copyKernel(cufftComplex* output,const float* HS, const int NN)
{
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id < NN)
	{
		output[id].x = HS[id];
		output[id].y = 0;
	}
}

__global__ void copyKernel(Complex* output, const cufftComplex* HS, const int NN)
{
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id < NN)
	{
		output[id].x = HS[id].x;
		output[id].y = HS[id].y;
	}
}
__global__ void copyKernel(float* output, const cufftComplex* HS, const int NN)
{
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id < NN)
	{
		output[id] = HS[id].x;

	}
}
int genHilbertKer(
		thrust::device_vector<cufftComplex>& FtS,
		const int YL,
		const float dYL)
{
	thrust::device_vector<float> n(2 * YL + 1, 0);
	thrust::sequence(n.begin(),n.end(),static_cast<float>(-YL));
	thrust::transform(n.begin(),n.end(),n.begin(),genKernel(dYL));
	const int NN = pow(2.0,ceil(log2(YL * 3.0)));

	FtS.resize(NN);
	thrust::device_vector<float> HS(NN,0);
	thrust::copy(n.begin(),n.begin() + YL,HS.end() - YL);
	thrust::copy(n.begin()+YL,n.end(),HS.begin());

	// To Do the CUFFT
	cufftHandle plan;
	dim3 blk(1024);
	dim3 gid((NN + blk.x - 1) / blk.x);
	copyKernel<<<gid,blk>>>(thrust::raw_pointer_cast(&FtS[0]),thrust::raw_pointer_cast(&HS[0]), NN);
	cufftPlan1d(&plan, NN, CUFFT_C2C,1);
	cufftExecC2C(plan, thrust::raw_pointer_cast(&FtS[0]), thrust::raw_pointer_cast(&FtS[0]), CUFFT_FORWARD);
	cufftDestroy(plan);
	HS.clear();
	n.clear();
	return NN;

}

__global__ void copyExpandedProjectionData(cufftComplex* output, const float* input,
		const int YL, const int ZL, const int ViewN, const int NN)
{
	int curPos = threadIdx.x + blockIdx.x * blockDim.x;
	int curBatch = threadIdx.y + blockIdx.y * blockDim.y;
	if(curPos < YL && curBatch < ZL * ViewN)
	{
		output[curBatch * NN + curPos].x = input[curBatch * YL + curPos];
		output[curBatch * NN + curPos].y = 0;
	}

}

__global__ void multiplyProjectionWithKernel(cufftComplex* proj, const cufftComplex* kernel, const int kernelLength, const int batchSize)
{
	int kerIdx = threadIdx.x + blockIdx.x * blockDim.x;
	int batIdx = threadIdx.y + blockIdx.y * blockDim.y;
	//__shared__ float kk[KSIZE];
	//kk[threadIdx.x] = kernel[kerIdx];
	//__syncthreads();
	if(kerIdx < kernelLength && batIdx < batchSize)
	{
		cufftComplex kk = kernel[kerIdx];
		cufftComplex res;
		res.x = proj[batIdx * kernelLength + kerIdx].x * kk.x - proj[batIdx * kernelLength + kerIdx].y * kk.y;
		res.y = proj[batIdx * kernelLength + kerIdx].x * kk.y + proj[batIdx * kernelLength + kerIdx].y * kk.x;
		proj[batIdx * kernelLength + kerIdx] = res;
		//proj[batIdx * kernelLength + kerIdx].y *= kk;
	}
}

__global__ void cutProjectionData(float* fpwd, cufftComplex* proj, const int YL, const int NN, const int batSize)
{
	int curIdx = threadIdx.x + blockIdx.x * blockDim.x;
	int batIdx = threadIdx.y + blockIdx.y * blockDim.y;
	if(curIdx < YL && batIdx < batSize)
	{
		fpwd[batIdx * YL + curIdx] = proj[batIdx * NN + curIdx].x / NN;
	}
}

// Let the projection data stored in the addressing order:
// 1. detector cell transversal direction (YL)
// 2. vertical direction (ZL)
// 3. view index (ViewN)
void filtering(
		thrust::device_vector<float>& fpwd, // Filtered projection data
		const thrust::device_vector<float>& Proj, // Projection data
		const int YL, const int ZL, const int ViewN, // Size of the projection data
		const float dYL)
{
	thrust::device_vector<cufftComplex> FtS;
	int NN = genHilbertKer(FtS, YL, dYL);
	//Expand the projection data
	thrust::device_vector<cufftComplex> exProj(NN * ZL * ViewN);
	dim3 copyExpBlk(32,32);
	dim3 copyExpGid(
			(YL + copyExpBlk.x - 1) / copyExpBlk.x,
			(ZL * ViewN + copyExpBlk.y - 1) / copyExpBlk.y);
	copyExpandedProjectionData<<<copyExpGid, copyExpBlk>>>(
			thrust::raw_pointer_cast(&exProj[0]),
			thrust::raw_pointer_cast(&Proj[0]),
			YL, ZL, ViewN, NN);
	// Forward Batch FFT
	cufftHandle plan;
	cufftPlan1d(&plan, NN, CUFFT_C2C, ZL * ViewN);
	cufftExecC2C(plan, thrust::raw_pointer_cast(&exProj[0]),
			thrust::raw_pointer_cast(&exProj[0]),CUFFT_FORWARD);
	// Multiply with the kernel
	dim3 multBlk(32,32);
	dim3 multGid(
			(NN + multBlk.x - 1) / multBlk.x,
			(ZL * ViewN + multBlk.y - 1) / multBlk.y);

	multiplyProjectionWithKernel<<<multGid,multBlk>>>(thrust::raw_pointer_cast(&exProj[0]),
			thrust::raw_pointer_cast(&FtS[0]), NN, ZL * ViewN);

	// Back batch FFT
	cufftExecC2C(plan, thrust::raw_pointer_cast(&exProj[0]),
				thrust::raw_pointer_cast(&exProj[0]),CUFFT_INVERSE);

	// Cut the data
	cutProjectionData<<<copyExpGid, copyExpBlk>>>(thrust::raw_pointer_cast(&fpwd[0]),
			thrust::raw_pointer_cast(&exProj[0]),
			YL, NN, ZL * ViewN);

	cufftDestroy(plan);
	FtS.clear();
	exProj.clear();
}

__global__ void preWeighting_ker(float* Proj,
		const int YL,
		const int ZL,
		const int ViewN,
		const float PLC,
		const float ZLC,
		const float dYL,
		const float dZL,
		const float SO)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x;
	int k = threadIdx.y + blockIdx.y * blockDim.y;
	int v = threadIdx.z + blockIdx.z * blockDim.z;
	if(j < YL && k < ZL && v < ViewN)
	{
		const float t = (j - PLC) * dYL;
		const float b = (k - ZLC) * dZL;
		const float wei = SO * SO / sqrtf(SO * SO * (SO * SO + b * b) - b * t * b * t);
		Proj[(v * ZL + k) * YL + j] *= wei;
	}
}

void preWeighting(
		thrust::device_vector<float>& Proj,
		const int YL,
		const int ZL,
		const int ViewN,
		const float PLC,
		const float ZLC,
		const float dYL,
		const float dZL,
		const float SO)
{
	dim3 blk(16,4,4);
	dim3 gid(
			(YL + blk.x - 1) / blk.x,
			(ZL + blk.y - 1) / blk.y,
			(ViewN + blk.z - 1) / blk.z);
	preWeighting_ker<<<gid,blk>>>(
		thrust::raw_pointer_cast(&Proj[0]),
		YL, ZL, ViewN, PLC, ZLC, dYL, dZL, SO);

}



//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void addressOrder(
		T* proj_ZYV,
		const T* proj_YZV,
		const int YL, const int ZL, const int ViewN)
{
	int zIdx = threadIdx.x + blockIdx.x * blockDim.x;
	int yIdx = threadIdx.y + blockIdx.y * blockDim.y;
	int vIdx = threadIdx.z + blockIdx.z * blockDim.z;
	if(zIdx < ZL && yIdx < YL && vIdx < ViewN)
	{
		proj_ZYV[(vIdx * YL + yIdx) * ZL + zIdx] =
				proj_YZV[(vIdx * ZL + zIdx) * YL + yIdx];

	}
}

template<typename T>
__global__ void addressOrder_2(
		T* proj_YZV,
		const T* proj_ZYV,
		const int YL, const int ZL, const int ViewN)
{
	int yIdx = threadIdx.x + blockIdx.x * blockDim.x;
	int zIdx = threadIdx.y + blockIdx.y * blockDim.y;
	int vIdx = threadIdx.z + blockIdx.z * blockDim.z;
	if(zIdx < ZL && yIdx < YL && vIdx < ViewN)
	{
		proj_YZV[(vIdx * ZL + zIdx) * YL + yIdx]
		         = proj_ZYV[(vIdx * YL + yIdx) * ZL + zIdx];

	}
}



extern "C"
void filtering(float* hfpwd,
		const float* hProj,
		const int YL, const int ZL, const int ViewN,
		const float PLC, const float ZLC,
		const float dYL, const float dZL,
		const float SO, const int GPUID)
{
    cudaSetDevice(GPUID);
	thrust::device_vector<float> Proj(hProj, hProj + YL * ZL * ViewN);
	preWeighting(Proj, YL, ZL, ViewN, PLC,
			 ZLC, dYL, dZL, SO);
	thrust::device_vector<float> fpwd(YL * ZL * ViewN, 0);
	filtering(fpwd,Proj, YL, ZL, ViewN, dYL);
	thrust::copy(fpwd.begin(),fpwd.end(),hfpwd);
	
}


// The function routine should be called
// FPWD = reWeihAdFiltr(Proj, PLC, ZLC, dYL, dZL, SO);
void mexFunction(
		int nlhs, mxArray* plhs[],
		int nrhs, const mxArray* prhs[])
{
	// YL, ZL, ViewN should be calculated from Proj
	if(nlhs != 1 || nrhs != 6)
	{
		std::cerr<<"Wrong number of parameters\n";
		std::cerr<<"This function requires 6 input parameter and provides one return value\n";
		exit(-1);
	}
	const mwSize* siz = mxGetDimensions(prhs[0]);
	const int YL = siz[0];
	const int ZL = siz[1];
	const int ViewN = siz[2];
	// Generate a new array
	plhs[0] = mxCreateNumericArray(3,siz,mxSINGLE_CLASS,mxREAL);
	float* fpwd = (float*)mxGetPr(plhs[0]);
	const float* Proj = (float*)mxGetPr(prhs[0]);
	const float PLC = *((float*)mxGetData(prhs[1]));
	const float ZLC = *((float*)mxGetData(prhs[2]));
	const float dYL = *((float*)mxGetData(prhs[3]));
	const float dZL = *((float*)mxGetData(prhs[4]));
	const float SO = *((float*)mxGetData(prhs[5]));
    

	filtering(fpwd, Proj, YL, ZL, ViewN,
			PLC - 1.0f, ZLC - 1.0f, dYL, dZL, SO, 1);
}














