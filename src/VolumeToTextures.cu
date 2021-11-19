#include "VolumeToTextures.h"
#include <thrust/device_vector.h>
#include <thrust/pair.h>
#include "cudaCheckReturner.h"

template<typename Ta, typename Tb>
static __global__ void naive_copyToTwoVolumes(Ta* in_ZXY,
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
static __global__ void naive_herizontalIntegral(Ta* in, Tb* out, int N, int ZN) {
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
static __global__ void naive_verticalIntegral(Ta* in, Tb* out, int N, int ZN) {
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
static __global__ void naive_verticalIntegral(double*  in, int2* out, int N, int ZN) {
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
static __global__ void verticalIntegral(T* prj, int ZN, int N) {
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < N) {
		int currentHead = idx * ZN;
		for (int ii = 1; ii < ZN; ++ii) {
			prj[currentHead + ii] += prj[currentHead + ii - 1];
		}
	}
}

template<typename T>
static __global__ void horizontalIntegral(T* prj, int DNU, int DNV, int PN) {
	int idv = threadIdx.x + blockIdx.x * blockDim.x;
	int pIdx = threadIdx.y + blockIdx.y * blockDim.y;
	if (idv < DNV && pIdx < PN)	{
		int headPtr = pIdx * DNU * DNV + idv;
		for (int ii = 1; ii < DNU; ++ii) {
			prj[headPtr + ii * DNV] += prj[headPtr + (ii - 1) * DNV];
		}
	}
}

static thrust::pair<thrust::device_vector<float>, thrust::device_vector<float>> generateSumAreaTableForVolume(thrust::device_vector<float>& vol,
	int XN, int YN, int ZN) {
	const int nZN = ZN + 1;
	const int nXN = XN + 1;
	const int nYN = YN + 1;
	const int nsiz_ZXY = nZN * nXN * YN;
	const int nsiz_ZYX = nZN * nYN * XN;

	thrust::device_vector<float> ZXY(nsiz_ZXY, 0);
	thrust::device_vector<float> ZYX(nsiz_ZYX, 0);
	float* pvol = thrust::raw_pointer_cast(&vol[0]);
	float* pZXY = thrust::raw_pointer_cast(&ZXY[0]);
	float* pZYX = thrust::raw_pointer_cast(&ZYX[0]);

	dim3 blk(64, 16, 1);
	dim3 gid(
		(ZN + blk.x - 1) / blk.x,
		(XN + blk.y - 1) / blk.y,
		(YN + blk.z - 1) / blk.z);

	naive_copyToTwoVolumes << <gid, blk >> > (pvol, pZXY, pZYX, XN, YN, ZN);

	blk.x = 32;
	blk.y = 1;
	blk.z = 1;
	gid.x = (nXN * YN + blk.x - 1) / blk.x;
	gid.y = 1;
	gid.z = 1;
	verticalIntegral << <gid, blk >> > (pZXY, nZN, nXN * YN);

	blk.x = 64;
	blk.y = 16;
	blk.z = 1;
	gid.x = (nZN + blk.x - 1) / blk.x;
	gid.y = (YN + blk.y - 1) / blk.y;
	gid.z = 1;

	horizontalIntegral << <gid, blk >> > (pZXY, nXN, nZN, YN);

	blk.x = 32;  blk.y = 1; blk.z = 1;
	gid.x = (nYN * XN + blk.x - 1) / blk.x;
	gid.y = 1;
	gid.z = 1;
	verticalIntegral << <gid, blk >> > (pZYX, nZN, nXN * YN);

	blk.x = 64;  blk.y = 16;  blk.z = 1;
	gid.x = (nZN + blk.x - 1) / blk.x;
	gid.y = (XN + blk.y - 1) / blk.y;
	gid.z = 1;

	horizontalIntegral << <gid, blk >> > (pZYX,	nYN, nZN, XN);

	return thrust::make_pair(ZXY, ZYX);
}

static thrust::pair<thrust::device_vector<float>, thrust::device_vector<float>> generateSumAreaTableForVolume(float* vol, int XN, int YN, int ZN) {
	const int size = XN * YN * ZN;
	thrust::device_vector<float> dvol(vol, vol + size);
	return generateSumAreaTableForVolume(dvol, XN, YN, ZN);
}

static thrust::pair<cudaArray*, cudaTextureObject_t> generateTexture(thrust::device_vector<float>& SATZXY, int width, int height, int depth, cudaTextureAddressMode textureAddressMode) {
    cudaExtent volumeSize;
    volumeSize.width = width;
    volumeSize.height = height;
    volumeSize.depth = depth;
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

    cudaArray* d_volumeArray;
    CUDA_CHECK_RETURN(cudaMalloc3DArray(&d_volumeArray, &channelDesc, volumeSize));

    cudaMemcpy3DParms copyParams = { 0 };
    copyParams.srcPtr = make_cudaPitchedPtr((void*)
        thrust::raw_pointer_cast(&SATZXY[0]),
        volumeSize.width * sizeof(float),
        volumeSize.width, volumeSize.height);
    copyParams.dstArray = d_volumeArray;
    copyParams.extent = volumeSize;
    copyParams.kind = cudaMemcpyDeviceToDevice;
    CUDA_CHECK_RETURN(cudaMemcpy3D(&copyParams));

    SATZXY.clear();

    cudaResourceDesc resDesc;
    memset(&resDesc, 0, sizeof(resDesc));

    resDesc.resType = cudaResourceTypeArray;
    resDesc.res.array.array = d_volumeArray;

    cudaTextureDesc texDesc;
    memset(&texDesc, 0, sizeof(texDesc));

    texDesc.addressMode[0] = textureAddressMode;
    texDesc.addressMode[1] = textureAddressMode;
    texDesc.addressMode[2] = textureAddressMode;
    texDesc.filterMode = cudaFilterModeLinear;
    texDesc.readMode = cudaReadModeElementType;
    texDesc.normalizedCoords = false;

    cudaTextureObject_t texObj = 0;
    CUDA_CHECK_RETURN(cudaCreateTextureObject(&texObj, &resDesc, &texDesc, nullptr));
    return thrust::make_pair(d_volumeArray, texObj);
}

VolumeToTextures::VolumeToTextures(float* vol, int XN, int YN, int ZN, cudaTextureAddressMode textureAddressMode) {
	
    auto summedAreaTablePair = generateSumAreaTableForVolume(vol, XN, YN, ZN);
	thrust::pair<cudaArray*, cudaTextureObject_t> P1 = generateTexture(summedAreaTablePair.first, ZN + 1, XN + 1, YN, cudaAddressModeClamp);
	thrust::pair<cudaArray*, cudaTextureObject_t> P2 = generateTexture(summedAreaTablePair.second, ZN + 1, YN + 1, XN, cudaAddressModeClamp);
    d_volumeArray1 = P1.first;
    d_volumeArray2 = P2.first;
    texObj1 = P1.second;
    texObj2 = P2.second;
}

VolumeToTextures::~VolumeToTextures() {
    CUDA_CHECK_RETURN(cudaDestroyTextureObject(texObj1));
    CUDA_CHECK_RETURN(cudaFreeArray(d_volumeArray1));	
    CUDA_CHECK_RETURN(cudaDestroyTextureObject(texObj2));
    CUDA_CHECK_RETURN(cudaFreeArray(d_volumeArray2));
}

cudaTextureObject_t VolumeToTextures::getTextureZXY() const {
    return texObj1;
}

cudaTextureObject_t VolumeToTextures::getTextureZYX() const {
    return texObj2;
}


VolumeToOneTexture::VolumeToOneTexture(float* hvol, int XN, int YN, int ZN, cudaTextureAddressMode textureAddressMode) {
	thrust::device_vector<float> vol(hvol, hvol + XN * YN * ZN);
	thrust::pair<cudaArray*, cudaTextureObject_t> arrayTexturePair = generateTexture(vol, ZN, XN, YN, textureAddressMode);
	d_volumeArray = arrayTexturePair.first;
	texObj = arrayTexturePair.second;
}

VolumeToOneTexture::~VolumeToOneTexture() {
	CUDA_CHECK_RETURN(cudaDestroyTextureObject(texObj));
	CUDA_CHECK_RETURN(cudaFreeArray(d_volumeArray));
}

cudaTextureObject_t VolumeToOneTexture::getTexture() const {
	return texObj;
}
