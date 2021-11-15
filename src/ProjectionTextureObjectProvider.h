#ifndef _PROJECTION_TEXTURE_OBJECT_PROVIDER_H_
#define _PROJECTION_TEXTURE_OBJECT_PROVIDER_H_
#include <cuda_runtime.h>
class ProjectionTextureObjectProvider {
private:
	cudaTextureObject_t texObj;
	cudaArray* d_prjArray;
public:
	cudaTextureObject_t getTexObj() const;
	cudaArray* getArray() const;
	void destory();
    template<typename T>
	ProjectionTextureObjectProvider(
		int Width, int Height, int Depth, // data size
		T* sourceData, // where is the data
		cudaMemcpyKind memcpyKind, // data from host or memory
		cudaTextureAddressMode addressMode, // how to address the texture (clamp, border ...)
		cudaTextureFilterMode textureFilterMode, // usually linear filtering (double --> int2 use pointer not linear interpolation)
		cudaTextureReadMode textureReadMode, // usually use element wise reading mode.
		bool isNormalized) {
        cudaExtent prjSize;
        prjSize.width = Width;
        prjSize.height = Height;
        prjSize.depth = Depth;
        cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<T>();
        cudaError_t cudaError = cudaMalloc3DArray(&d_prjArray, &channelDesc, prjSize);
        if (cudaError != cudaSuccess) {
            // Throw exception that the cuda memory cannot be allocated.
        }
        cudaMemcpy3DParms copyParams = { 0 };
        copyParams.srcPtr = make_cudaPitchedPtr(
            (void*)sourceData, prjSize.width * sizeof(T),
            prjSize.width, prjSize.height);
        copyParams.dstArray = d_prjArray;
        copyParams.extent = prjSize;
        copyParams.kind = memcpyKind;
        cudaError = cudaMemcpy3D(&copyParams);
        if (cudaError != cudaSuccess) {
            // Todo: return false;
        }
        cudaResourceDesc resDesc;
        memset(&resDesc, 0, sizeof(resDesc));
        assert(&resDesc != NULL);

        resDesc.resType = cudaResourceTypeArray;
        resDesc.res.array.array = d_prjArray;
        cudaTextureDesc texDesc;
        memset(&texDesc, 0, sizeof(texDesc));
        assert(&texDesc != NULL);

        texDesc.addressMode[0] = addressMode;
        texDesc.addressMode[1] = addressMode;
        texDesc.addressMode[2] = addressMode;
        texDesc.filterMode = textureFilterMode;
        texDesc.readMode = textureReadMode;

        texDesc.normalizedCoords = isNormalized;
        cudaError = cudaCreateTextureObject(&texObj, &resDesc, &texDesc, nullptr);
        if (cudaError != cudaSuccess) {
            // Todo: return false;
        }
    }
    ~ProjectionTextureObjectProvider();
};

#endif