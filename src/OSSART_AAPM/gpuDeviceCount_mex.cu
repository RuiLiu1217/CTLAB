#include "mex.h"
#include "matrix.h"

#include <stdio.h>
#include <cuda_runtime.h>

void getGPUNumber(int& nDevices) {
    cudaGetDeviceCount(&nDevices);
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    const mwSize dims[] = {1};
    plhs[0] = mxCreateNumericArray(1,dims,mxINT32_CLASS,mxREAL);
    int* num = (int*)mxGetPr(plhs[0]);
    getGPUNumber(*num);
}