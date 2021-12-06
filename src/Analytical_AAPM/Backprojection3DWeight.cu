//////////////////////////////////
//  Imaging and Informatics Lab
//  Department of Electrical and Computer Engineering
//  University of Massachusetts Lowell
//  Parallel-cone backprojection
//  May 28, 2016
//////////////////////////////////

#include <mex.h>
#include <cmath>
#include <complex>
#include <vector>
#include <cstdlib>
#include <thread>
#include <iostream>
#include <fstream>
#include <thread>
#include <omp.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <cuda_runtime.h>
#include <cassert>
#ifndef nullptr
#define nullptr NULL
#endif

#define TYPE double
const static TYPE pi = 3.14159265358979;
const static float TWOPI = pi * 2.0f;

#if DEBUG
#define CUDA_CHECK_RETURN(value) { cudaError_t _m_cudaStat = value; if(_m_cudaStat != cudaSuccess){fprintf(stderr, "Error %s at line %d in file %s\n", cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__); exit(1);}}
#else
#define CUDA_CHECK_RETURN(value) {value;}
#endif


void Backprojection3DWeightingCPU(float* ctimage, const float* fp, const float* w, const int N_2pi,
		const int XN, const int YN, const int ZN,
		const float XNC, const float YNC, const float ZNC,
		const float dx, const float dy, const float dz,
		const int YL, const int ZL, const int PN,
		const float YLC, const float ZLC, const float dYL, const float dZL,
		const float h, const float BetaS, const float DeltaFai,
		const float N_pi, const float HSCoef, const float SO, const float SD, const float k1) {   
    const float invPI = 1.0f / TWOPI;
    
    int nProcessors=omp_get_max_threads();
    
    omp_set_num_threads(nProcessors);            
#pragma omp parallel for
    for(int yi=0;yi<YN;yi++) {
        const float y = (yi-YNC)*dy;
        for(int xi=0;xi<XN;xi++) {
            const float x  = (xi-XNC)*dx;

            for(int zi = 0; zi < ZN; ++zi) {
                ///compute the projection position for every grid on the image plane
                const float z = (zi - ZNC) * dz;
                const float Beta0 = 2.0f * pi * z / h;
                const int s0 = ceilf((Beta0 - BetaS) / DeltaFai - 0.5f);
                int s1 = s0 - ceilf(N_pi * HSCoef);
                int s2 = s0 + ceilf(N_pi * HSCoef);
                float res = 0; // used to accumulate the results

                if((s1 < PN) && (s2 > 0)) {
                    s1 = (s1 < 0)? 0 : s1;
                    s2 = (s2 > PN - 1)? PN - 1 : s2;
                    for(int ProjInd = s1; ProjInd <= s2; ++ProjInd) {
                        const float View = BetaS + ProjInd * DeltaFai;
                        const int d1   = N_pi-(s0-ProjInd); //d1 = ProjInd;
                        const int d2 = (ProjInd < s0)?(d1 + N_pi) : (d1 - N_pi);
                        const float UU = -x * cosf(View) - y * sinf(View);
                        const float Yr = -x * sinf(View) + y * cosf(View);
                        const float temp1 = sqrtf(SO * SO - Yr * Yr);
                        const float temp2 = (z-h * View * invPI);
                        const float Zr = temp2*(SO*SD)/(UU * temp1+SO * SO - Yr * Yr);
                        const float U1 = Yr/dYL+YLC;
                        const int U  = ceilf(U1);
                        const float V1 = Zr/dZL+ZLC;
                        const int V  = ceilf(V1);
                        const float Dey = U-U1;
                        const float Dez = V-V1;
                        //Linear interploate
                        if ( (U > 0) && (U < YL) && (V > 0) && (V < ZL)) {
                            const float touying =
                                    Dey *          Dez * fp[(ProjInd * ZL + (V-1)) * YL + (U-1)] +
                                    Dey * (1.0f - Dez) * fp[(ProjInd * ZL + (V)) * YL + (U-1)] +
                                    (1.0f - Dey) * Dez * fp[(ProjInd * ZL + (V-1)) * YL + U] +
                                    (1.0f - Dey) * (1.0f - Dez) * fp[(ProjInd * ZL + V) * YL + U];
                            const float weight1 = w[d1];
                            const float weight2 = w[d2];
                            assert(temp1 + UU != 0);
                            const float Gama   = fabsf( temp2 / ( temp1 + UU));
                            assert(temp1 != UU);
                            const float Gama_C = (ProjInd < s0) ? fabsf((temp2 - 0.5f) / (temp1 - UU)) : fabsf((temp2 + 0.5f) / (temp1 - UU));
                            const float m1 = powf(Gama,  k1);    //m1     = std::real(std::pow(Gama,k1*h));
                            const float m2 = powf(Gama_C, k1);  //m2     = std::real(std::pow(Gama_C,k1*h));
                            const float tw = (weight2 * m1 + weight1 * m2);
                            float weight = 0;
                            if(tw != 0) {
                                weight = (weight1*m2) / tw;
                            }                              
                            res += weight*touying*DeltaFai;
                        } // end if linear interpolation
                    } // end for projection
                } // end if range
                ctimage[(yi*XN+xi)*ZN+zi] = res;
            } // end zi
        } // xi
    } // yi
}


__global__ void ParallekFDKHelical3DWeightingKer_texture_grouped(
		float* ctimage,
		cudaTextureObject_t* texObj, const int prjPerGroup,
		float2* __restrict cossinV,
		float* __restrict w,
		const int XN, const int YN, const int ZN, // cannot change
		float x, float y, float z, // used /XNC/ YNC/ ZNC
		float dx, float dy, float dz,
		//const int YL, const int ZL,
		const int PN,
		float YLC, float ZLC, //cannot change
		float invdYL, float invdZL,
		float h_div_twopi,
		float BetaS, const float DeltaFai,
		float N_pi,
		int NHSCOEF, float SOSO, float SOSD, float k1)
{
	int zi = threadIdx.x + blockIdx.x * blockDim.x;
	int xi = threadIdx.y + blockIdx.y * blockDim.y;
	int yi = threadIdx.z + blockIdx.z * blockDim.z;
	if(zi < ZN && xi < XN && yi < YN)
	{

		x = (xi - x) * dx; // z,y,z cannot change
		y = (yi - y) * dy;
		z = (zi - z) * dz;

		// Calculate Beta0 according to the position
		int s0 = ceilf((z / h_div_twopi - BetaS) / DeltaFai - 0.5); //Center projection index ;
		int s1 = s0 - NHSCOEF;
		int s2 = s0 + NHSCOEF;
		dx = 0;

		s1 = (s1 < 0)?0:s1;
		s2 = (s2 > PN-1)?PN-1:s2;
		BetaS = BetaS + s1 * DeltaFai;
		for(; s1 <= s2; ++s1)
		{
			int d1 = N_pi - s0 + s1; 
			int d2 = (s1 < s0)?(d1 + N_pi) : (d1 - N_pi);
			float weight1 = w[d1];
			float weight2 = w[d2];

			dy = (z - BetaS * h_div_twopi);

			float UU = -x * cossinV[s1].x - y * cossinV[s1].y; // x,y NO!
			float Yr = -x * cossinV[s1].y + y * cossinV[s1].x;
			dz = sqrtf(SOSO - Yr * Yr);

			float Zr = dy * (SOSD) / ((UU + dz) * dz);

			float U1 = Yr * invdYL + YLC;
			float V1 = Zr * invdZL + ZLC;
			d1 = s1 / prjPerGroup;
			d2 = s1 - prjPerGroup * d1;
			float touying = tex3D<float>(texObj[d1],U1+0.5, V1+0.5, d2+0.5);
			assert(dz + UU != 0);
			float m1 = powf(fabsf(dy / (dz + UU)), k1) * weight2;
			assert(dz != UU);
			float m2 = powf((s1 < s0) ? fabsf((dy - 0.5) / (dz - UU)) : fabsf((dy + 0.5) / (dz - UU)), k1) * weight1;
			assert(m1 + m2 != 0);
			float weight = m2 / (m1 + m2);

			dx += weight * touying;
			BetaS += DeltaFai;
		}
		ctimage[(yi*XN+xi)*ZN+zi] = dx * DeltaFai;
	}

}

void createTextureObject(
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
		(void*) sourceData, prjSize.width * sizeof(float),
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



void destroyTextureObject(cudaTextureObject_t& texObj, cudaArray* d_array) {
    CUDA_CHECK_RETURN(cudaDestroyTextureObject(texObj));
	CUDA_CHECK_RETURN(cudaFreeArray(d_array));
}


//// use texture to accelerate the backprojection
//void ParallekFDKHelical3DWeighting_GPU_texture_grouped(
//		float* hctimage,
//		const float* hfp,
//		const float* hw, const int N_2pi,
//		const int XN, const int YN, const int ZN,
//		const float XNC, const float YNC, const float ZNC,
//		const float dx, const float dy, const float dz,
//		const int YL, const int ZL, const int PN,
//		const float YLC, const float ZLC,
//		const float dYL, const float dZL,
//		const float h, const float BetaS, const float DeltaFai,
//		const float N_pi, const float HSCoef, const float SO, const float SD, const float k1, int threadx, int thready, int threadz)
//{
//	// Get number of devices
//	int nDevices;
//	cudaGetDeviceCount(&nDevices);
//	std::vector<cudaStream_t> cudaStream(nDevices);
//	std::vector<cudaEvent_t> cudaEvent(nDevices);
//	std::vector<thrust::device_vector<float>> subCTImage(nDevices);
//	std::vector<thrust::device_vector<float>> subW(nDevices);
//	int SubPN = PN / nDevices;
//	if (SubPN * nDevices != PN)
//	{
//		SubPN++;
//	}
//
//	int prjNumPerGroup = 2048; // TODO: further reduce the number of views if the total view is actually very small.
//	int groupNum = SubPN / prjNumPerGroup;
//	int lasGrpNum = PN - prjNumPerGroup * groupNum; // We think that all GPU has the same computing capacity.
//	if (lasGrpNum != 0)
//	{
//		++groupNum;
//	}
//
//	std::vector<thrust::host_vector<thrust::device_vector<float>>> subFp(nDevices);
//	std::vector<std::vector<cudaTextureObject_t>> subTexObj(nDevices);
//	std::vector<std::vector<cudaArray *>> sub_d_prjArray(nDevices);
//
//
//	std::vector<thrust::device_vector<cudaTextureObject_t>> sub_dtexObj(nDevices);// = texObj;
//	std::vector<thrust::host_vector<float2>> sub_cossinV(nDevices); // (PN);
//	std::vector<thrust::device_vector<float2>> sub_dcossinV(nDevices);// = cossinV;
//
//
//	thrust::host_vector<float2> cossinV(SubPN * nDevices);
//	for (int ProjInd = 0; ProjInd != PN; ++ProjInd)
//	{
//		const float View = BetaS + ProjInd * DeltaFai;
//		cossinV[ProjInd] = make_float2(cosf(View), sinf(View));
//	}
//	
//
//
//	for (int deviceIdx = 0; deviceIdx != nDevices; ++deviceIdx)
//	{
//		cudaSetDevice(deviceIdx);
//		cudaStreamCreate(&cudaStream[deviceIdx]);
//		cudaEventCreate(&cudaEvent[deviceIdx]);
//		subCTImage[deviceIdx].resize(XN * YN * ZN);
//		subW[deviceIdx].resize(N_2pi);
//		thrust::copy(hw, hw + N_2pi, subW[deviceIdx].begin());
//		subFp[deviceIdx].resize(groupNum);
//		for (int i = 0; i != groupNum - 1; ++i)
//		{
//			subFp[deviceIdx][i].resize(YL * ZL * prjNumPerGroup);
//			thrust::fill(subFp[deviceIdx][i].begin(), subFp[deviceIdx][i].end(), 0);
//			thrust::copy(hfp + deviceIdx * SubPN * YL * ZL + i * YL * ZL * prjNumPerGroup,
//				hfp + deviceIdx * SubPN * YL * ZL + (i + 1)* YL * ZL * prjNumPerGroup,
//				subFp[deviceIdx][i].begin());
//		}
//		subFp[deviceIdx][groupNum - 1].resize(YL * ZL * prjNumPerGroup);
//		thrust::fill(subFp[deviceIdx][groupNum - 1].begin(), subFp[deviceIdx][groupNum - 1].end(), 0);
//		if (lasGrpNum != 0) // It seems with wrong memory accessing calculation
//		{
//			//thrust::copy(hfp + deviceIdx * YL * ZL * SubPN + (groupNum - 1) * YL * ZL * prjNumPerGroup,
//			//	hfp + deviceIdx * YL * ZL * SubPN + (groupNum - 1) * YL * ZL * prjNumPerGroup
//			//	+ YL * ZL * lasGrpNum, subFp[deviceIdx][groupNum - 1].begin());
//		}
//		else
//		{
//			//thrust::copy(hfp + deviceIdx * YL * ZL * SubPN + (groupNum - 1) * YL * ZL * prjNumPerGroup,
//			//	hfp + deviceIdx * YL * ZL * SubPN + groupNum * YL * ZL * prjNumPerGroup,
//			//	subFp[deviceIdx][groupNum - 1].begin());
//		}
//		subTexObj[deviceIdx].resize(groupNum);
//		sub_d_prjArray[deviceIdx].resize(groupNum);
//
//		cudaTextureAddressMode addressMode = cudaAddressModeBorder;
//		cudaTextureFilterMode textureFilterMode = cudaFilterModeLinear;
//		cudaTextureReadMode textureReadMode = cudaReadModeElementType;
//		for (int i = 0; i != groupNum; ++i)
//		{
//			createTextureObject(subTexObj[deviceIdx][i], sub_d_prjArray[deviceIdx][i],
//				YL, ZL, prjNumPerGroup, thrust::raw_pointer_cast(&subFp[deviceIdx][i][0]),
//				cudaMemcpyDeviceToDevice, addressMode, textureFilterMode, textureReadMode, false);
//		}
//		sub_dtexObj[deviceIdx].resize(groupNum);
//		sub_dtexObj[deviceIdx] = subTexObj[deviceIdx];
//
//		sub_dcossinV[deviceIdx].resize(SubPN);
//		thrust::copy(cossinV.begin() + deviceIdx * SubPN, cossinV.begin() + (deviceIdx + 1) * SubPN, sub_dcossinV[deviceIdx].begin());
//	}
//	cudaDeviceSynchronize();
//
//	dim3 blk(threadx,thready,threadz); // It can be configured as you like!
//	dim3 gid(
//			(ZN + blk.x - 1) / blk.x,
//			(XN + blk.y - 1) / blk.y,
//			(YN + blk.z - 1) / blk.z);
//
//	int NHSCOEF = ceilf(N_pi * HSCoef);
//	float h_div_twopi = h / TWOPI;
//	for (int deviceIdx = 0; deviceIdx != nDevices; ++deviceIdx)
//	{
//		cudaSetDevice(deviceIdx);
//		ParallekFDKHelical3DWeightingKer_texture_grouped << <gid, blk, 0, cudaStream[deviceIdx] >> >(
//			thrust::raw_pointer_cast(&subCTImage[deviceIdx][0]),
//			thrust::raw_pointer_cast(&sub_dtexObj[deviceIdx][0]), prjNumPerGroup,
//			thrust::raw_pointer_cast(&sub_dcossinV[deviceIdx][0]),
//			thrust::raw_pointer_cast(&subW[deviceIdx][0]), XN, YN, ZN, XNC, YNC, ZNC,
//			dx, dy, dz,
//			SubPN, YLC, ZLC, 1.0 / dYL, 1.0 / dZL, h_div_twopi, BetaS, DeltaFai, N_pi,
//			NHSCOEF, SO * SO, SO * SD, k1);
//	}
//
//	cudaDeviceSynchronize();
//	std::vector<thrust::host_vector<float>> subHostCTImage(nDevices);
//	for (int deviceIdx = 0; deviceIdx != nDevices; ++deviceIdx)
//	{
//		cudaSetDevice(deviceIdx);
//		for (int i = 0; i != groupNum; ++i)
//		{
//			destroyTextureObject(subTexObj[deviceIdx][i], sub_d_prjArray[deviceIdx][i]);
//		}
//		subHostCTImage[deviceIdx].resize(XN * YN * ZN);
//		thrust::copy(subCTImage[deviceIdx].begin(), subCTImage[deviceIdx].end(), subHostCTImage[deviceIdx].begin());
//		for (int i = 0; i != XN * YN * ZN; ++i)
//		{
//			hctimage[i] += subHostCTImage[deviceIdx][i];
//		}
//	}
//	cudaDeviceSynchronize();
//	/*
//	thrust::copy(ctimage.begin(),ctimage.end(),hctimage);*/
//}


// use texture to accelerate the backprojection
void ParallekFDKHelical3DWeighting_GPU_texture_grouped_old(
		float* hctimage,
		const float* hfp,
		const float* hw, const int N_2pi,
		const int XN, const int YN, const int ZN,
		const float XNC, const float YNC, const float ZNC,
		const float dx, const float dy, const float dz,
		const int YL, const int ZL, const int PN,
		const float YLC, const float ZLC,
		const float dYL, const float dZL,
		const float h, const float BetaS, const float DeltaFai,
		const float N_pi, const float HSCoef, const float SO, const float SD, const float k1, int threadx, int thready, int threadz)
{

	dim3 blk(threadx,thready,threadz); // It can be configured as you like!
	dim3 gid(
			(ZN + blk.x - 1) / blk.x,
			(XN + blk.y - 1) / blk.y,
			(YN + blk.z - 1) / blk.z);

	thrust::device_vector<float> ctimage(XN * YN * ZN, 0);
	thrust::device_vector<float> w(hw, hw + N_2pi);

	// Divide the projection into several groups
	int prjNumPerGroup = 4096; // We need to change this value according to the requirement, the maximum dim
	int groupNum = PN / prjNumPerGroup; // first N groups
	int lasGrpNum = PN - prjNumPerGroup * groupNum;
	if (lasGrpNum != 0) {
		++groupNum;
	}

	thrust::host_vector<thrust::device_vector<float> > fp(groupNum);
	for(int i = 0; i != groupNum - 1; ++i) {
		fp[i].resize(YL * ZL * prjNumPerGroup);
		thrust::fill(fp[i].begin(),fp[i].end(), 0);
		thrust::copy(hfp + i * YL * ZL * prjNumPerGroup,
				hfp + (i+1) * YL * ZL * prjNumPerGroup, fp[i].begin());
	}
	fp[groupNum - 1].resize(YL * ZL * prjNumPerGroup);
	thrust::fill(fp[groupNum - 1].begin(),fp[groupNum - 1].end(), 0);
	if(lasGrpNum != 0) {
		thrust::copy(hfp + (groupNum - 1) * YL * ZL * prjNumPerGroup,
				hfp + (groupNum - 1) * YL * ZL * prjNumPerGroup
				+ YL * ZL * lasGrpNum, fp[groupNum - 1].begin());
	} else {
		thrust::copy(hfp + (groupNum - 1) * YL * ZL * prjNumPerGroup,
				hfp + groupNum * YL * ZL * prjNumPerGroup,
				fp[groupNum - 1].begin());
	}

	std::vector<cudaTextureObject_t> texObj(groupNum);
	std::vector<cudaArray*> d_prjArray(groupNum);
	cudaTextureAddressMode addressMode = cudaAddressModeBorder;
	cudaTextureFilterMode textureFilterMode = cudaFilterModeLinear;
	cudaTextureReadMode textureReadMode = cudaReadModeElementType;
	for(int i = 0; i != groupNum; ++i) {
		createTextureObject(texObj[i], d_prjArray[i],
			YL, ZL, prjNumPerGroup,
			thrust::raw_pointer_cast(&fp[i][0]),
			cudaMemcpyDeviceToDevice,
			addressMode,
			textureFilterMode,
			textureReadMode,
			false);
	}

	thrust::device_vector<cudaTextureObject_t> dtexObj = texObj;
	thrust::host_vector<float2> cossinV(PN);

	for(int ProjInd = 0; ProjInd != PN; ++ProjInd) {
		const float View = BetaS + ProjInd * DeltaFai;
		cossinV[ProjInd] = make_float2(cosf(View), sinf(View));

	}
	thrust::device_vector<float2> dcossinV = cossinV;

	int NHSCOEF = ceilf(N_pi * HSCoef);
	float h_div_twopi = h / TWOPI;

	ParallekFDKHelical3DWeightingKer_texture_grouped<<<gid,blk>>>(
			thrust::raw_pointer_cast(&ctimage[0]),
			thrust::raw_pointer_cast(&dtexObj[0]), prjNumPerGroup,
			thrust::raw_pointer_cast(&dcossinV[0]),
			thrust::raw_pointer_cast(&w[0]), XN, YN, ZN, XNC, YNC, ZNC,
			dx, dy, dz,
			PN, YLC, ZLC, 1.0 / dYL, 1.0 / dZL, h_div_twopi, BetaS, DeltaFai, N_pi,
			NHSCOEF, SO * SO, SO * SD, k1);

	for(int i = 0; i != groupNum; ++i) {
		destroyTextureObject(texObj[i], d_prjArray[i]);
	}
	thrust::copy(ctimage.begin(),ctimage.end(),hctimage);
}

// use texture to accelerate the backprojection
void ParallekFDKHelical3DWeighting_GPU_texture_grouped_openMP(
		float* hctimage,
		const float* hfp,
		const float* hw, const int N_2pi,
		const int XN, const int YN, const int ZN,
		const float XNC, const float YNC, const float ZNC,
		const float dx, const float dy, const float dz,
		const int YL, const int ZL, const int PN,
		const float YLC, const float ZLC,
		const float dYL, const float dZL,
		const float h, const float BetaS, const float DeltaFai,
		const float N_pi, const float HSCoef, const float SO, const float SD, const float k1, int threadx, int thready, int threadz)
{

	dim3 blk(threadx,thready,threadz); // It can be configured as you like!
	dim3 gid(
			(ZN + blk.x - 1) / blk.x,
			(XN + blk.y - 1) / blk.y,
			(YN + blk.z - 1) / blk.z);

	thrust::device_vector<float> ctimage(XN * YN * ZN, 0);
	thrust::device_vector<float> w(hw, hw + N_2pi);

	// Divide the projection into several groups
	int prjNumPerGroup = 4096; // We need to change this value according to the requirement, the maximum dim
	int groupNum = PN / prjNumPerGroup; // first N groups
	int lasGrpNum = PN - prjNumPerGroup * groupNum;
	if (lasGrpNum != 0)
	{
		++groupNum;
	}
	thrust::host_vector<thrust::device_vector<float> > fp(groupNum);
	for(int i = 0; i != groupNum - 1; ++i)
	{
		fp[i].resize(YL * ZL * prjNumPerGroup);
		thrust::fill(fp[i].begin(),fp[i].end(), 0);
		thrust::copy(hfp + i * YL * ZL * prjNumPerGroup,
				hfp + (i+1) * YL * ZL * prjNumPerGroup, fp[i].begin());
	}
	fp[groupNum - 1].resize(YL * ZL * prjNumPerGroup);
	thrust::fill(fp[groupNum - 1].begin(),fp[groupNum - 1].end(), 0);
	if(lasGrpNum != 0)
	{
		thrust::copy(hfp + (groupNum - 1) * YL * ZL * prjNumPerGroup,
				hfp + (groupNum - 1) * YL * ZL * prjNumPerGroup
				+ YL * ZL * lasGrpNum, fp[groupNum - 1].begin());
	}
	else
	{
		thrust::copy(hfp + (groupNum - 1) * YL * ZL * prjNumPerGroup,
				hfp + groupNum * YL * ZL * prjNumPerGroup,
				fp[groupNum - 1].begin());
	}

	std::vector<cudaTextureObject_t> texObj(groupNum);
	std::vector<cudaArray*> d_prjArray(groupNum);
	cudaTextureAddressMode addressMode = cudaAddressModeBorder;
	cudaTextureFilterMode textureFilterMode = cudaFilterModeLinear;
	cudaTextureReadMode textureReadMode = cudaReadModeElementType;
	for(int i = 0; i != groupNum; ++i)
	{
		createTextureObject(texObj[i], d_prjArray[i],
			YL, ZL, prjNumPerGroup,
			thrust::raw_pointer_cast(&fp[i][0]),
			cudaMemcpyDeviceToDevice,
			addressMode,
			textureFilterMode,
			textureReadMode,
			false);
	}
	thrust::device_vector<cudaTextureObject_t> dtexObj = texObj;
	thrust::host_vector<float2> cossinV(PN);
	for(int ProjInd = 0; ProjInd != PN; ++ProjInd)
	{
		const float View = BetaS + ProjInd * DeltaFai;
		cossinV[ProjInd] = make_float2(cosf(View), sinf(View));

	}
	thrust::device_vector<float2> dcossinV = cossinV;

	int NHSCOEF = ceilf(N_pi * HSCoef);
	float h_div_twopi = h / TWOPI;

	ParallekFDKHelical3DWeightingKer_texture_grouped<<<gid,blk>>>(
			thrust::raw_pointer_cast(&ctimage[0]),
			thrust::raw_pointer_cast(&dtexObj[0]), prjNumPerGroup,
			thrust::raw_pointer_cast(&dcossinV[0]),
			thrust::raw_pointer_cast(&w[0]), XN, YN, ZN, XNC, YNC, ZNC,
			dx, dy, dz,
			PN, YLC, ZLC, 1.0 / dYL, 1.0 / dZL, h_div_twopi, BetaS, DeltaFai, N_pi,
			NHSCOEF, SO * SO, SO * SD, k1);

	for(int i = 0; i != groupNum; ++i)
	{
		destroyTextureObject(texObj[i], d_prjArray[i]);
	}
	thrust::copy(ctimage.begin(),ctimage.end(),hctimage);

}






void Backprojection3DWeightingGPU(
		float* hctimage, const float* hfp, const float* hw, const int N_2pi,
		const int XN, const int YN, const int ZN,
		const float XNC, const float YNC, const float ZNC,
		const float dx, const float dy, const float dz,
		const int YL, const int ZL, const int PN,
		const float YLC, const float ZLC, const float dYL, const float dZL,
		const float h, const float BetaS, const float DeltaFai,
		const float N_pi, const float HSCoef, const float SO, const float SD, const float k1, int threadx, int thready, int threadz)
{
  
  ParallekFDKHelical3DWeighting_GPU_texture_grouped_old(hctimage, hfp, hw, N_2pi,
			XN, YN, ZN, XNC, YNC, ZNC, dx, dy, dz,
			YL, ZL, PN, YLC, ZLC, dYL, dZL, h, BetaS, DeltaFai,
			N_pi, HSCoef, SO, SD, k1, threadx, thready, threadz);

}


void ParallelFDKHelical3DWeighting(
		double* ctimage, // Image to be reconstructed
		const double* fp,
		const double SO,
		const double DO,
		const int YL,
		const int ZL,
		const double DecWidth,
		const double DecHeigh,
		const double YLC,
		const double ZLC,
		const double h,
		const double BetaS,
		const int PN,
		const int N_2pi,
		const double ObjR,
		const int XN,
		const int YN,
		const int ZN,
		const double XNC,
		const double YNC,
		const double ZNC,
        const double dx,
        const double dy,
        const double dz,
		const int delta,
		const double HSCoef,
		const double k1,
        const int useGPU, 
        int threadx,
        int thready,
        int threadz)
{
        const double PI = 3.141592653589793;
        const double N_pi = N_2pi / 2.0;
        const double dYL = DecWidth / YL;
        const double dZL = DecHeigh / ZL;
        const double DeltaFai = 2.0 * PI / N_2pi;
        const double SD = SO + DO;
        

        std::vector<double> w(N_2pi,0);
        const int L = 2.0 * ceil(N_pi * HSCoef) + 1.0;
        const int Shift = N_pi - ceil(N_pi * HSCoef);
    
        ////Producing the weighting function

        for (int k=0;k<L;k++)
        {
            if (0 <= k && k<delta)
                w[k+Shift]= std::pow(cos((pi/2)*(delta-k-0.5)/delta),2);
            else if(L-delta<=k && k < L)
                w[k+Shift]= std::pow(cos((pi/2)*(k-(L-delta)+0.5)/delta),2);
            else
                w[k+Shift] = 1;
        }

        
        if(useGPU==1) {
            // use single floating data to implement the backprojection
             std::vector<float> fw(w.begin(),w.end());
             float* fctimage = new float[XN * YN * ZN];

             float* ffp = new float[YL * ZL * PN];
             std::copy(fp, fp + YL * ZL * PN, ffp);
             Backprojection3DWeightingGPU(fctimage,ffp, &fw[0],N_2pi,
                     XN, YN, ZN, XNC, YNC, ZNC, dx, dy, dz,
                      YL, ZL, PN, YLC, ZLC, dYL, dZL, h, BetaS, DeltaFai,
                      N_pi, HSCoef, SO, SD, k1, threadx, thready, threadz);
             std::copy(fctimage, fctimage + XN * YN * ZN, ctimage);
             delete[] fctimage;
             delete[] ffp;
             
        } else {
            std::vector<float> fctimage(ctimage, ctimage + XN * YN * ZN);
            std::vector<float> ffp(fp, fp + YL * ZL * PN);
            std::vector<float> fw(w.begin(), w.end());

            Backprojection3DWeightingCPU(&fctimage[0], &ffp[0], &fw[0], N_2pi,
                    XN, YN, ZN, XNC, YNC, ZNC, dx, dy, dz,
                     YL, ZL, PN, YLC, ZLC, dYL, dZL, h, BetaS, DeltaFai,
                     N_pi, HSCoef, SO, SD, k1);
            std::copy(fctimage.begin(),fctimage.end(), ctimage);
        }
}


static const TYPE *geom[4];
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    TYPE ObjR,SO,DO,YLC,ZLC,XNC,YNC,ZNC,DecWidth,DecHeigh,h,dYL,dZL,dx,dy,dz,DeltaFai,BetaS,k1,HSCoef;
    int  YL,ZL,XN,YN,ZN,N_2pi,delta,N_pi,PN;
    TYPE *ctimage;
    const TYPE *fp;
    /* Check for proper number of arguments */
    if (nrhs != 7) {
        mexErrMsgTxt("Backward projection needs 7 inputs.");
    }
    if (nlhs != 0) {
        mexErrMsgTxt("Backward projection does not need any output.");
    }
    
    geom[0] = mxGetPr(mxGetFieldByNumber(prhs[0], 0, 0));
    SO            = TYPE(geom[0][0]);
    DO            = TYPE(geom[0][1]);
    YL            = int(TYPE(geom[0][2])+0.5);
    ZL            = int(TYPE(geom[0][3])+0.5);
    DecWidth      = TYPE(geom[0][4]);
    DecHeigh      = TYPE(geom[0][5]);
    YLC           = TYPE(geom[0][6]);
    ZLC           = TYPE(geom[0][7]);
    h             =     TYPE(geom[0][8])*DecHeigh;
        
    geom[1] = mxGetPr(mxGetFieldByNumber(prhs[0], 0, 1));
    BetaS         =     TYPE(geom[1][0]);
    PN            =     int(geom[1][2]);
    N_2pi         = int(geom[1][3]);

    geom[2] = mxGetPr(mxGetFieldByNumber(prhs[0], 0, 2));
    ObjR          =    TYPE(geom[2][0]);
    XN            = int(TYPE(geom[2][1])+0.5);
    YN            = int(TYPE(geom[2][2])+0.5);
    ZN            = int(TYPE(geom[2][3])+0.5);
    XNC           =     TYPE(geom[2][4]);
    YNC           =     TYPE(geom[2][5]);
    ZNC           =     TYPE(geom[2][6]);
    dx            =     TYPE(geom[2][7]);
    dy            =     TYPE(geom[2][8]);
    dz            =     TYPE(geom[2][9]);
    
    geom[3] = mxGetPr(mxGetFieldByNumber(prhs[0], 0, 3));
    delta         = int(TYPE(geom[3][0])+0.5);
    HSCoef        =     TYPE(geom[3][1]);
    k1            =     TYPE(geom[3][2])*TYPE(geom[0][8]);

    fp            =     mxGetPr(prhs[1]);
    ctimage       =     mxGetPr(prhs[2]);
    int useGPU = *((int*)mxGetPr(prhs[3]));
    int threadx   = *((int*)mxGetPr(prhs[4]));
    int thready   = *((int*)mxGetPr(prhs[5]));
    int threadz   = *((int*)mxGetPr(prhs[6]));
    
    N_pi = N_2pi/2.0;
    dYL = DecWidth/YL;
    dZL = DecHeigh/ZL;
    DeltaFai = 2.0*pi/N_2pi;
    
    ParallelFDKHelical3DWeighting(ctimage, fp, SO, DO, YL, ZL,
		DecWidth, DecHeigh, YLC - 1.0f, ZLC - 1.0f, h, BetaS, PN, N_2pi,
		ObjR, XN, YN, ZN, XNC - 1.0f, YNC - 1.0f, ZNC - 1.0f, dx, dy, dz, delta, HSCoef, k1,useGPU, threadx, thready, threadz);

}

