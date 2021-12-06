/////////////////////////////////////////////////////
// Imaging and Informatics Lab
// Department of Electrical and Computer Engineering
// University of Massachusetts Lowell
// Graduate School of Arts and Sciences
// Wake Forest University
// \brief Rebinning the projection data
// \author Rui Liu, Hengyong Yu
// \date Jan. 15, 2016
// \version 1.0
//////////////////////////////////////////////////////
#include <mex.h>
#include <omp.h>
#include <cmath>
#include <assert.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>

// Let's make the order of the projection data ZL, YL, ViewN
extern "C"
void ParallelRebinningCBCurve_CPU(
		float* outputProj,
		const float* Proj,
		const int YL, const int ZL, const int ViewN,
		const float YLC, const float dYA,
		const float DeltaTheta,
		const float PLC, const float DeltaT, const float DeltaFai,
		const float SO)
{
#pragma omp parallel for
	for(int i = 0; i < ViewN; i++)
	{
		float Theta = i * DeltaTheta; // The View for the parallel projection
		for(int j = 0; j != YL; ++j)
		{
			float t = (j - PLC) * DeltaT;
			float Beta = asinf(t / SO);
			float Fai = Theta + Beta;
			float a = atanf(t / sqrtf(SO*SO-t*t));
			float FaiIndex = (Fai / DeltaFai);
			float UIndex = a / dYA + YLC;
			int FI = ceilf(FaiIndex);
			int UI = ceilf(UIndex);
			float coeXB = FI - FaiIndex;
			float coeXU = 1.0f - coeXB;
			float coeYB = UI - UIndex;
			float coeYU = 1.0f - coeYB;

			int IndexXU(0);
			int IndexXB(0);
			int IndexYU(0);
			int IndexYB(0);

			if(FI <= 0)
			{
				IndexXU = 0;
				IndexXB = 0;
			}
			else if(FI > ViewN - 1)
			{
				IndexXU = ViewN - 1;
				IndexXB = ViewN - 1;
			}
			else
			{
				IndexXU = FI;
				IndexXB = FI - 1.0;
			}

			if(UI <= 0)
			{
				IndexYU = 0;
				IndexYB = 0;
			}
			else if(UI > YL - 1)
			{
				IndexYU = YL - 1;
				IndexYB = YL - 1;
			}
			else
			{
				IndexYU = UI;
				IndexYB = UI - 1;
			}

			for(int k = 0; k != ZL; ++k)
			{
				outputProj[(i * YL + j) * ZL + k] =
						coeXB * coeYB * Proj[(IndexXB * YL + IndexYB) * ZL + k] +
						coeXU * coeYB * Proj[(IndexXU * YL + IndexYB) * ZL + k] +
						coeXB * coeYU * Proj[(IndexXB * YL + IndexYU) * ZL + k] +
						coeXU * coeYU * Proj[(IndexXU * YL + IndexYU) * ZL + k];
			}
		}
	}
}


__global__ void ParallelRebinningCBCurve_GPU_ker(
		float* outputProj,
		const float* Proj,
		const int YL, const int ZL, const int ViewN,
		const float YLC, 
        const float dYA, // Detector corresponding 
		const float DeltaTheta, // view step
        const float PLC,
		const float DeltaT, // ideal detector size
        const float DeltaFai,
        const float SO)
{
	int k = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int i = threadIdx.z + blockIdx.z * blockDim.z;
	if(i < ViewN && j < YL && k < ZL)
	{
        float Theta = i * DeltaTheta;
        float t = (j - PLC) * DeltaT;
		float Beta = asinf(t / SO);
		float Fai = Theta + Beta;
		float a = atanf(t / sqrtf(SO*SO-t*t));
		float FaiIndex = (Fai / DeltaFai);
		float UIndex = a / dYA + YLC;
		int FI = ceilf(FaiIndex);
		int UI = ceilf(UIndex);
		float coeXB = FI - FaiIndex;
		float coeXU = 1.0f - coeXB;
		float coeYB = UI - UIndex;
		float coeYU = 1.0f - coeYB;

		int IndexXU(0);
		int IndexXB(0);
		int IndexYU(0);
		int IndexYB(0);

		if(FI <= 0)
		{
			IndexXU = 0;
			IndexXB = 0;
		}
		else if(FI >= ViewN - 1)
		{
			IndexXU = ViewN - 1;
			IndexXB = ViewN - 1;
		}
		else
		{
			IndexXU = FI;
			IndexXB = FI - 1.0;
		}

		if(UI <= 0)
		{
			IndexYU = 0;
			IndexYB = 0;
		}
		else if(UI >= YL - 1)
		{
			IndexYU = YL - 1;
			IndexYB = YL - 1;
		}
		else
		{
			IndexYU = UI;
			IndexYB = UI - 1;
		}
		outputProj[(i * YL + j) * ZL + k] =
			coeXB * coeYB * Proj[(IndexXB * YL + IndexYB) * ZL + k] +
			coeXU * coeYB * Proj[(IndexXU * YL + IndexYB) * ZL + k] +
			coeXB * coeYU * Proj[(IndexXB * YL + IndexYU) * ZL + k] +
			coeXU * coeYU * Proj[(IndexXU * YL + IndexYU) * ZL + k];
	}
}


void ParallelRebinningCBCurve_GPU_fcn(
		float* outputProj,
		const float* Proj,
		const int YL, const int ZL, const int ViewN,
		const float YLC, 
        const float dYA, // Detector corresponding 
		const float DeltaTheta, // view step
        const float PLC,
		const float DeltaT, // ideal detector size
        const float DeltaFai,
        const float SO, int threadx, int thready, int threadz)
{
	dim3 blk(threadx,thready,threadz);
	dim3 gid(
		(ZL + blk.x - 1) / blk.x,
		(YL + blk.y - 1) / blk.y,
		(ViewN + blk.z - 1) / blk.z);
	thrust::device_vector<float> dProj(Proj, Proj + YL * ZL * ViewN);
	thrust::device_vector<float> doutputProj(dProj.size(),0);
	ParallelRebinningCBCurve_GPU_ker<<<gid,blk>>>(
			thrust::raw_pointer_cast(&doutputProj[0]),
			thrust::raw_pointer_cast(&dProj[0]), YL, ZL, ViewN,
			YLC, dYA, DeltaTheta, PLC, DeltaT, DeltaFai, SO);

	thrust::copy(doutputProj.begin(),doutputProj.end(),outputProj);
	dProj.clear();
	doutputProj.clear();
}


extern "C"
void ParallelRebinningCBCurve_GPU(
		float* outputProj,
		const float* Proj,
		const int YL, const int ZL, const int ViewN,
		const float YLC, 
        const float dYA, // Detector corresponding 
		const float DeltaTheta, // view step
        const float PLC,
		const float DeltaT, // ideal detector size
        const float DeltaFai,
        const float SO, int threadx, int thready, int threadz)
{
    cudaDeviceReset();
    cudaSetDevice(0);
    ParallelRebinningCBCurve_GPU_fcn(outputProj, Proj,
		YL, ZL, ViewN, YLC, dYA, DeltaTheta, PLC, DeltaT, 
        DeltaFai, SO, threadx, thready, threadz);
    cudaDeviceReset();
}






void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // The projection is ordered in ZL, YL, ViewN order
    float* Proj = (float*)mxGetPr(prhs[0]);
    int YL = *((int*)(mxGetPr(prhs[1])));
    int ZL = *((int*)(mxGetPr(prhs[2])));
    int ViewN = *((int*)(mxGetPr(prhs[3])));
    float YLC = *((float*)(mxGetPr(prhs[4])));
    float dYA = *((float*)(mxGetPr(prhs[5])));
    float DeltaTheta = *((float*)(mxGetPr(prhs[6])));
    float PLC = *((float*)(mxGetPr(prhs[7])));
    float DeltaT = *((float*)(mxGetPr(prhs[8])));
    float DeltaFai = *((float*)(mxGetPr(prhs[9])));
    float SO = *((float*)(mxGetPr(prhs[10])));
    int useGPU = *((int*)(mxGetPr(prhs[11])));
    int threadx = *((int*)(mxGetPr(prhs[12])));
    int thready = *((int*)(mxGetPr(prhs[13])));
    int threadz = *((int*)(mxGetPr(prhs[14])));
    
    const mwSize dims[]={static_cast<mwSize>(ZL),static_cast<mwSize>(YL),static_cast<mwSize>(ViewN)};
    plhs[0] = mxCreateNumericArray(3,dims,mxSINGLE_CLASS,mxREAL);
    float* outputProj = (float*)mxGetPr(plhs[0]);
    if(useGPU == 1)
    {
          ParallelRebinningCBCurve_GPU(outputProj, Proj,
              YL, ZL, ViewN, YLC - 1.0f, dYA, DeltaTheta,
              PLC - 1.0f, DeltaT, DeltaFai, SO, threadx, thready, threadz);
    }
    else
    {
          ParallelRebinningCBCurve_CPU(outputProj, Proj,
             YL, ZL, ViewN, YLC - 1.0f, dYA, DeltaTheta,
             PLC - 1.0f, DeltaT, DeltaFai, SO);
    }
}
