#ifndef CUDA_CHECK_RETURNER_H_
#define CUDA_CHECK_RETURNER_H_

// Usually, for safety, we will use checkCudaErrors or CUDA_CHECK_RETURN before "cuda" API functions
#if DEBUG
#define CUDA_CHECK_RETURN(value) { cudaError_t _m_cudaStat = value; if(_m_cudaStat != cudaSuccess){fprintf(stderr, "Error %s at line %d in file %s\n", cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__); exit(1);}}
#else
#define CUDA_CHECK_RETURN(value) {value;}
#endif


#endif