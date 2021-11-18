#ifndef TEXTURE_OBJECT_PROVIDER_H_
#define TEXTURE_OBJECT_PROVIDER_H_
#include <cuda_runtime.h>
#include <cuda.h>
//Create texture object and corresponding cudaArray function
template<typename T>
void createTextureObject(
	cudaTextureObject_t& texObj, //return: texture object pointing to the cudaArray
	cudaArray* d_prjArray, // return: cudaArray storing the data
	int Width, int Height, int Depth, // data size
	T* sourceData, // where is the data
	cudaMemcpyKind memcpyKind, // data from host or memory
	cudaTextureAddressMode addressMode, // how to address the texture (clamp, border ...)
	cudaTextureFilterMode textureFilterMode, // usually linear filtering (double --> int2 use pointer not linear interpolation)
	cudaTextureReadMode textureReadMode, // usually use element wise reading mode.
	bool isNormalized) { // usually false
	cudaExtent prjSize;
	prjSize.width = Width;
	prjSize.height = Height;
	prjSize.depth = Depth;
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<T>();
	CUDA_CHECK_RETURN(cudaMalloc3DArray(&d_prjArray, &channelDesc, prjSize));
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
	CUDA_CHECK_RETURN(cudaCreateTextureObject(&texObj, &resDesc, &texDesc, nullptr));
}

// Destroy a GPU array and corresponding TextureObject
void destroyTextureObject(cudaTextureObject_t& texObj, cudaArray* d_array) {
	CUDA_CHECK_RETURN(cudaDestroyTextureObject(texObj));
	CUDA_CHECK_RETURN(cudaFreeArray(d_array));
}
#endif