#include "TextureObjectProvider.h"
#include "utilities.hpp"

void destroyTextureObject(cudaTextureObject_t& texObj, cudaArray* d_array) {
	CUDA_CHECK_RETURN(cudaDestroyTextureObject(texObj));
	CUDA_CHECK_RETURN(cudaFreeArray(d_array));
}