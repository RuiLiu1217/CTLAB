#include "utilities.hpp"
/**
 * The multi-GPUs based 2D multi slices projection and backprojection
 * Author: Rui Liu
 * Date: Sep. 18, 2016
 */
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <vector_types.h>

#include "multiSlices_ker.cuh"
typedef unsigned char byte;

#define BLKX 32
#define BLKY 8
#define BLKZ 1

namespace DD2
{
	struct CosSinFunctor
	{
		__host__ __device__ float2 operator()(float ang)
		{
			return make_float2(cos(ang),sin(ang));
		}
	};

	// Split the projection data
	void splitProjection(
			thrust::host_vector<thrust::host_vector<float> >& subProj,
			thrust::host_vector<thrust::host_vector<float2> >& subCossin,
			float* proj, thrust::host_vector<float2>& cossin,
			const int SLN, const int DNU, const int PN,
			thrust::host_vector<int> sSLN, const int gpuNum)
	{
		int psum = 0;
		for(int i = 0; i != gpuNum; ++i)
		{
			subProj[i].resize(sSLN[i] * DNU * PN);
			subCossin[i].resize(sSLN[i] * PN);
			int curPos = sSLN[i];
			for(int p = 0; p != DNU * PN; ++p)
			{
				for(int s = 0; s != sSLN[i]; ++s)
				{
					int subPos = p * sSLN[i] + s;
					int totPos = p * SLN + (s + psum);
					subProj[i][subPos] = proj[totPos];
				}
			}

			for(int p = 0; p != PN; ++p)
			{
				for(int s = 0; s != sSLN[i]; ++s)
				{
					int subPos = p * sSLN[i] + s;
					int totPos = p * SLN + (s + psum);
					subCossin[i][subPos] = cossin[totPos];
				}
			}

			psum += sSLN[i];
		}
	}


	void combineProjection(
			thrust::host_vector<thrust::host_vector<float> >& subProj,
			float* proj, const int SLN, const int DNU, const int PN,
			std::vector<int>& sSLN, const int gpuNum)
	{
		int psum = 0;
		for(int i = 0; i != gpuNum; ++i)
		{
			int curPos = sSLN[i];
			for(int p = 0; p < DNU * PN; ++p)
			{
				for(int s = 0; s != sSLN[i]; ++s)
				{
					int subPos = p * sSLN[i] + s;
					int totPos = p * SLN + (s + psum);
					proj[totPos] = subProj[i][subPos];
				}
			}
			psum += sSLN[i];
		}
	}


	void combineVolume(
			thrust::host_vector<thrust::host_vector<float> >& subVol,
			float* vol, const int SLN, const int XN, const int YN,
			thrust::host_vector<int>& sSLN, const int gpuNum)
	{
		int psum = 0;
		//omp_set_num_threads();
		for(int i = 0; i < gpuNum; ++i)
		{
			int curPos = sSLN[i];
#pragma omp parallel for
			for(int p = 0; p < XN * YN; ++p)
			{
				for(int s = 0; s != sSLN[i]; ++s)
				{
					int subPos = p * sSLN[i] + s;
					int totPos = p * SLN + (s + psum);
					vol[totPos] = subVol[i][subPos];
				}
			}
			psum += sSLN[i];
		}
	}




	// Split the volume
	void splitVolume(
			std::vector<std::vector<float> >& subVol,
			thrust::host_vector<thrust::host_vector<float2> >& subCossin,
			float* vol,
			thrust::host_vector<float2> cossin,
			const int SLN, const int XN, const int YN, const int PN,
			std::vector<int>& sSLN, const int gpuNum)
	{
		int psum = 0;
		for(int i = 0; i != gpuNum; ++i)
		{
			subVol[i].resize(sSLN[i] * XN * YN);
			int curPos = sSLN[i];
			for(int p = 0; p != XN * YN; ++p)
			{
				for(int s = 0; s != sSLN[i]; ++s)
				{
					int subPos = p * sSLN[i] + s;
					int totPos = p * SLN + (s + psum);
					subVol[i][subPos] = vol[totPos];
				}
			}

			for(int p = 0; p != PN; ++p)
			{
				for(int s = 0; s != sSLN[i]; ++s)
				{
					int subPos = p * sSLN[i] + s;
					int totPos = p * SLN + (s + psum);
					subCossin[i][subPos] = cossin[totPos];
				}
			}

			psum += sSLN[i];
		}
	}



	// Copy the volume from the original to
	template<typename Ta, typename Tb>
	__global__ void naive_copyToTwoVolumes(Ta* in_ZXY,
		Tb* out_ZXY, Tb* out_ZYX,
		int XN, int YN, int ZN)
	{
		int idz = threadIdx.x + blockIdx.x * blockDim.x;
		int idx = threadIdx.y + blockIdx.y * blockDim.y;
		int idy = threadIdx.z + blockIdx.z * blockDim.z;
		if (idx < XN && idy < YN && idz < ZN)
		{
			int i = (idy * XN + idx) * ZN + idz;
			int ni = (idy * (XN + 1) + (idx + 1)) * ZN + idz;
			int nj = (idx * (YN + 1) + (idy + 1)) * ZN + idz;

			out_ZXY[ni] = in_ZXY[i];
			out_ZYX[nj] = in_ZXY[i];
		}
	}


	__global__ void horizontalIntegral(float* prj, int DNU, int DNV, int PN)
	{
		int idv = threadIdx.x + blockIdx.x * blockDim.x;
		int pIdx = threadIdx.y + blockIdx.y * blockDim.y;
		if (idv < DNV && pIdx < PN)
		{
			int headPrt = pIdx * DNU * DNV + idv;
			for (int ii = 1; ii < DNU; ++ii)
			{
				prj[headPrt + ii * DNV] = prj[headPrt + ii * DNV] + prj[headPrt + (ii - 1) * DNV];
			}
		}
	}


	__global__ void addOneSidedZeroBoarder(const float* prj_in, float* prj_out, int DNU, int DNV, int PN)
	{
		int idv = threadIdx.x + blockIdx.x * blockDim.x;
		int idu = threadIdx.y + blockIdx.y * blockDim.y;
		int pn = threadIdx.z + blockIdx.z * blockDim.z;
		if (idu < DNU && idv < DNV && pn < PN)
		{
			int i = (pn * DNU + idu) * DNV + idv;
			int ni = (pn * (DNU + 1) + (idu + 1)) * (DNV + 1) + idv + 1;
			prj_out[ni] = prj_in[i];
		}
	}


	__global__ void verticalIntegral2(float* prj, int ZN, int N)
	{
		int idx = threadIdx.x + blockIdx.x * blockDim.x;
		if (idx < N)
		{
			int currentHead = idx * ZN;
			for (int ii = 1; ii < ZN; ++ii)
			{
				prj[currentHead + ii] = prj[currentHead + ii] + prj[currentHead + ii - 1];
			}
		}
	}


	__global__ void heorizontalIntegral2(float* prj, int DNU, int DNV, int PN)
	{
		int idv = threadIdx.x + blockIdx.x * blockDim.x;
		int pIdx = threadIdx.y + blockIdx.y * blockDim.y;
		if (idv < DNV && pIdx < PN)
		{
			int headPrt = pIdx * DNU * DNV + idv;
			for (int ii = 1; ii < DNU; ++ii)
			{
				prj[headPrt + ii * DNV] = prj[headPrt + ii * DNV] + prj[headPrt + (ii - 1) * DNV];
			}
		}
	}

	__global__ void addOneSidedZeroBoarder_multiSlice_Fan(const float* prj_in, float* prj_out, int DNU, int SLN, int PN)
	{
		int idv = threadIdx.x + blockIdx.x * blockDim.x;
		int idu = threadIdx.y + blockIdx.y * blockDim.y;
		int pn = threadIdx.z + blockIdx.z * blockDim.z;
		if (idu < DNU && idv < SLN && pn < PN)
		{
			int i = (pn * DNU + idu) * SLN + idv;
			int ni = (pn * (DNU + 1) + (idu + 1)) * SLN + idv;
			prj_out[ni] = prj_in[i];
		}
	}


	__global__ void heorizontalIntegral_multiSlice_Fan(float* prj, int DNU, int SLN, int PN)
	{
		int idv = threadIdx.x + blockIdx.x * blockDim.x;
		int pIdx = threadIdx.y + blockIdx.y * blockDim.y;
		if (idv < SLN && pIdx < PN)
		{
			int headPrt = pIdx * DNU * SLN + idv;
			for (int ii = 1; ii < DNU; ++ii)
			{
				prj[headPrt + ii * SLN] = prj[headPrt + ii * SLN] + prj[headPrt + (ii - 1) * SLN];
			}
		}
	}



}



__global__  void MultiSlices_DDPROJ_ker(
	cudaTextureObject_t volTex1,
	cudaTextureObject_t volTex2,
	float* proj,
	float2 s, // source position
	const float2* __restrict__ cossin,
	const float* __restrict__ xds,
	const float* __restrict__ yds,
	const float* __restrict__ bxds,
	const float* __restrict__ byds,
	float2 objCntIdx,
	float dx,
	int XN, int YN, int SLN,
	int DNU, int PN)
{
	int slnIdx = threadIdx.x + blockIdx.x * blockDim.x;
	int detIdU = threadIdx.y + blockIdx.y * blockDim.y;
	int angIdx = threadIdx.z + blockIdx.z * blockDim.z;
	if(slnIdx < SLN && detIdU < DNU && angIdx < PN)
	{
		float2 dir = cossin[angIdx * SLN + slnIdx]; // cossin;

		float2 cursour = make_float2(
			s.x * dir.x - s.y * dir.y,
			s.x * dir.y + s.y * dir.x); // current source position;
		s = dir;

		float2 curDet = make_float2(
				xds[detIdU] * s.x - yds[detIdU] * s.y,
				xds[detIdU] * s.y + yds[detIdU] * s.x);

		float2 curDetL = make_float2(
				bxds[detIdU] * s.x - byds[detIdU] * s.y,
				bxds[detIdU] * s.y + byds[detIdU] * s.x);

		float2 curDetR = make_float2(
				bxds[detIdU+1] * s.x - byds[detIdU+1] * s.y,
				bxds[detIdU+1] * s.y + byds[detIdU+1] * s.x);

		dir = normalize(curDet - cursour);

		float factL = 0;
		float factR = 0;
		float constVal = 0;
		float obj = 0;
		float realL = 0;
		float realR = 0;
		float intersectLength = 0;

		float invdx = 1.0f / dx;
		float summ;
		if(fabsf(s.x) <= fabsf(s.y))
		{

			summ = 0;
			factL = (curDetL.y - cursour.y) / (curDetL.x - cursour.x);
			factR = (curDetR.y - cursour.y) / (curDetR.x - cursour.x);

			constVal = dx / fabsf(dir.x);
#pragma unroll
			for (int ii = 0; ii < XN; ++ii)
			{
				obj = (ii - objCntIdx.x) * dx;

				realL = (obj - curDetL.x) * factL + curDetL.y;
				realR = (obj - curDetR.x) * factR + curDetR.y;

				intersectLength = realR - realL;
				realL = realL * invdx + objCntIdx.y + 1;
				realR = realR * invdx + objCntIdx.y + 1;

				summ += (tex3D<float>(volTex2, slnIdx + 0.5f, realR, ii + 0.5) - tex3D<float>(volTex2, slnIdx + 0.5, realL, ii + 0.5)) / intersectLength;

			}
			__syncthreads();
			proj[(angIdx * DNU + detIdU) * SLN + slnIdx] = summ * constVal;

		}
		else
		{
			summ = 0;
			factL = (curDetL.x - cursour.x) / (curDetL.y - cursour.y);
			factR = (curDetR.x - cursour.x) / (curDetR.y - cursour.y);

			constVal = dx / fabsf(dir.y);
#pragma unroll
			for (int ii = 0; ii < YN; ++ii)
			{
				obj = (ii - objCntIdx.y) * dx;

				realL = (obj - curDetL.y) * factL + curDetL.x;
				realR = (obj - curDetR.y) * factR + curDetR.x;

				intersectLength = realR - realL;
				realL = realL * invdx + objCntIdx.x + 1;
				realR = realR * invdx + objCntIdx.x + 1;

				summ += (tex3D<float>(volTex1, slnIdx + 0.5f, realR, ii + 0.5) - tex3D<float>(volTex1, slnIdx + 0.5, realL, ii + 0.5)) / intersectLength;
			}
			__syncthreads();
			proj[(angIdx * DNU + detIdU) * SLN + slnIdx] = summ * constVal;
		}
	}
}



void MultiSlices_DDPROJ(
		float* hvol, // the pointer to the image
		float* hprj, // the pointer to the projection (SLN, DNU, PN) order
		const float x0, const float y0, //position of the initial source
		float* xds, float* yds, // distribution of the detector cells
		const int DNU, // Number of detector cells
		const int SLN, // Number of slices to be projected or backprojected
		const float imgXCenter, const float imgYCenter, //Center of the image
		const int XN, const int YN, // pixel number of the image
		const float dx, // size of the pixel
		float* h_angs, // view angles SHOULD BE WITH SIZE SLN * PN
		int PN, // # of view angles
		byte* mask,
		int* startidx,
		const int gpuNum)
{
	// Regular the projection
	for(int i = 0; i != XN * YN; ++i)
	{
		byte v = mask[i];
		for(int z = 0; z != SLN; ++z)
		{
			hvol[i * SLN + z] *= v;
		}
	}

	float* bxds = new float[DNU + 1];
	float* byds = new float[DNU + 1];

	DD3Boundaries<float>(DNU + 1, xds, bxds);
	DD3Boundaries<float>(DNU + 1, yds, byds);

	const float objCntIdxX = (XN - 1.0) * 0.5 - imgXCenter / dx;
	const float objCntIdxY = (YN - 1.0) * 0.5 - imgYCenter / dx;

	std::vector<int> startIdx(startidx, startidx + gpuNum);
	std::vector<int> endIdx(gpuNum);
	std::copy(startIdx.begin() + 1,
			startIdx.end(), endIdx.begin());
	endIdx[gpuNum - 1] = SLN;
	std::vector<int> sSLN(gpuNum);

	//Split the volumes
	std::vector<std::vector<float> > subVol(gpuNum);
	thrust::host_vector<thrust::host_vector<float2> > subCossin(gpuNum);
	for(int i = 0; i != gpuNum; ++i)
	{
		//subVol[i].resize(sSLN[i] * XN * YN);
		sSLN[i] = endIdx[i] - startIdx[i];
		subCossin[i].resize(sSLN[i] * PN);
	}
	thrust::host_vector<float2> cossin(PN * SLN);
	thrust::transform(h_angs, h_angs + PN * SLN,
			cossin.begin(),[=](float ang){
		return make_float2(cosf(ang), sinf(ang));
	});

	DD2::splitVolume(subVol, subCossin, hvol, cossin,
			SLN, XN, YN, PN, sSLN, gpuNum);


	// Generate multiple streams
	std::vector<cudaStream_t> stream(gpuNum);

	std::vector<int> siz(gpuNum, 0);
	std::vector<int> nsiz_ZXY(gpuNum, 0);
	std::vector<int> nsiz_ZYX(gpuNum, 0);

	thrust::host_vector<thrust::device_vector<float> > SATZXY(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > SATZYX(gpuNum);

	thrust::host_vector<cudaExtent> volumeSize1(gpuNum);
	thrust::host_vector<cudaExtent> volumeSize2(gpuNum);
	thrust::host_vector<cudaChannelFormatDesc> channelDesc1(gpuNum);
	thrust::host_vector<cudaChannelFormatDesc> channelDesc2(gpuNum);
	thrust::host_vector<cudaArray*> d_volumeArray1(gpuNum);
	thrust::host_vector<cudaArray*> d_volumeArray2(gpuNum);
	thrust::host_vector<cudaMemcpy3DParms> copyParams1(gpuNum);
	thrust::host_vector<cudaMemcpy3DParms> copyParams2(gpuNum);
	thrust::host_vector<cudaResourceDesc> resDesc1(gpuNum);
	thrust::host_vector<cudaResourceDesc> resDesc2(gpuNum);
	thrust::host_vector<cudaTextureDesc> texDesc1(gpuNum);
	thrust::host_vector<cudaTextureDesc> texDesc2(gpuNum);
	thrust::host_vector<cudaTextureObject_t> texObj1(gpuNum);
	thrust::host_vector<cudaTextureObject_t> texObj2(gpuNum);

	thrust::host_vector<thrust::device_vector<float> > d_prj(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > d_xds(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > d_yds(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > d_bxds(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > d_byds(gpuNum);
	thrust::host_vector<thrust::device_vector<float2> > d_cossin(gpuNum);

	thrust::host_vector<thrust::host_vector<float> > h_prj(gpuNum);

	dim3 pblk(BLKX,BLKY,BLKZ);
	thrust::host_vector<dim3> pgid(gpuNum);

	omp_set_num_threads(gpuNum);

#pragma omp parallel for
	for(int i = 0 ; i < gpuNum; ++i)
	{
		cudaSetDevice(i);

		siz[i] = XN * YN * sSLN[i];
		nsiz_ZXY[i] = sSLN[i] * (XN + 1) * YN;
		nsiz_ZYX[i] = sSLN[i] * (YN + 1) * XN;
		SATZXY[i].resize(nsiz_ZXY[i]);
		SATZYX[i].resize(nsiz_ZYX[i]);

		thrust::device_vector<float> vol = subVol[i];

		dim3 blk(64,16,1);
		dim3 gid(
			(sSLN[i] + blk.x - 1) / blk.x,
			(XN      + blk.y - 1) / blk.y,
			(YN      + blk.z - 1) / blk.z);

		DD2::naive_copyToTwoVolumes<<<gid,blk,0,stream[i]>>>(
			(thrust::raw_pointer_cast(&vol[0])),
			(thrust::raw_pointer_cast(&SATZXY[i][0])),
			(thrust::raw_pointer_cast(&SATZYX[i][0])),
			XN,YN,sSLN[i]);
		vol.clear();

		blk.x = 64;
		blk.y = 16;
		blk.z = 1;
		gid.x = (sSLN[i] + blk.x - 1) / blk.x;
		gid.y = (YN + blk.y - 1) / blk.y;
		gid.z = 1;

		DD2::horizontalIntegral<<<gid, blk, 0, stream[i]>>>(
			thrust::raw_pointer_cast(&SATZXY[i][0]),
			XN + 1, sSLN[i], YN);

		blk.x = 64;
		blk.y = 16;
		blk.z = 1;
		gid.x = (sSLN[i] + blk.x - 1) / blk.x;
		gid.y = (XN + blk.y - 1) / blk.y;
		gid.z = 1;

		DD2::horizontalIntegral<<<gid, blk, 0, stream[i]>>>(
			thrust::raw_pointer_cast(&SATZYX[i][0]),
			YN + 1, sSLN[i], XN);

		volumeSize1[i].width = sSLN[i];
		volumeSize1[i].height = XN + 1;
		volumeSize1[i].depth = YN;
		volumeSize2[i].width = sSLN[i];
		volumeSize2[i].height = YN + 1;
		volumeSize2[i].depth = XN;

		channelDesc1[i] = cudaCreateChannelDesc<float>();
		channelDesc2[i] = cudaCreateChannelDesc<float>();

		cudaMalloc3DArray(&d_volumeArray1[i], &channelDesc1[i], volumeSize1[i]);
		cudaMalloc3DArray(&d_volumeArray2[i], &channelDesc2[i], volumeSize2[i]);

		copyParams1[i].srcPtr = make_cudaPitchedPtr((void*)
				thrust::raw_pointer_cast(&SATZXY[i][0]),
				volumeSize1[i].width * sizeof(float),
				volumeSize1[i].width, volumeSize1[i].height);
		copyParams1[i].dstArray = d_volumeArray1[i];
		copyParams1[i].extent = volumeSize1[i];
		copyParams1[i].kind = cudaMemcpyDeviceToDevice;

		copyParams2[i].srcPtr = make_cudaPitchedPtr((void*)
			thrust::raw_pointer_cast(&SATZYX[i][0]),
			volumeSize2[i].width * sizeof(float),
			volumeSize2[i].width, volumeSize2[i].height);
		copyParams2[i].dstArray = d_volumeArray2[i];
		copyParams2[i].extent = volumeSize2[i];
		copyParams2[i].kind = cudaMemcpyDeviceToDevice;

		cudaMemcpy3D(&copyParams1[i]);
		cudaMemcpy3D(&copyParams2[i]);

		SATZXY[i].clear();
		SATZYX[i].clear();

		memset(&resDesc1[i], 0, sizeof(resDesc1[i]));
		memset(&resDesc2[i], 0, sizeof(resDesc2[i]));

		resDesc1[i].resType = cudaResourceTypeArray;
		resDesc2[i].resType = cudaResourceTypeArray;

		resDesc1[i].res.array.array = d_volumeArray1[i];
		resDesc2[i].res.array.array = d_volumeArray2[i];

		memset(&texDesc1[i], 0, sizeof(texDesc1[i]));
		memset(&texDesc2[i], 0, sizeof(texDesc2[i]));

		texDesc1[i].addressMode[0] = cudaAddressModeClamp;
		texDesc1[i].addressMode[1] = cudaAddressModeClamp;
		texDesc1[i].addressMode[2] = cudaAddressModeClamp;

		texDesc2[i].addressMode[0] = cudaAddressModeClamp;
		texDesc2[i].addressMode[1] = cudaAddressModeClamp;
		texDesc2[i].addressMode[2] = cudaAddressModeClamp;

		texDesc1[i].filterMode = cudaFilterModeLinear;
		texDesc2[i].filterMode = cudaFilterModeLinear;

		texDesc1[i].readMode = cudaReadModeElementType;
		texDesc2[i].readMode = cudaReadModeElementType;

		texDesc1[i].normalizedCoords = false;
		texDesc2[i].normalizedCoords = false;

		cudaCreateTextureObject(&texObj1[i], &resDesc1[i], &texDesc1[i], nullptr);
		cudaCreateTextureObject(&texObj2[i], &resDesc2[i], &texDesc2[i], nullptr);

		d_prj[i].resize(sSLN[i] * DNU * PN);
		h_prj[i].resize(sSLN[i] * DNU * PN);


		d_xds[i].resize(DNU);
		thrust::copy(xds,xds+DNU, d_xds[i].begin());

		d_yds[i].resize(DNU);
		thrust::copy(yds,yds+DNU, d_yds[i].begin());

		d_bxds[i].resize(DNU + 1);
		thrust::copy(bxds,bxds+DNU+1, d_bxds[i].begin());

		d_byds[i].resize(DNU + 1);
		thrust::copy(byds,byds+DNU+1, d_byds[i].begin());

		d_cossin[i].resize(sSLN[i] * PN);
		thrust::copy(subCossin[i].begin(), subCossin[i].end(), d_cossin[i].begin());

		pgid[i].x = (sSLN[i] + pblk.x - 1) / pblk.x;
		pgid[i].y = (DNU + pblk.y - 1) / pblk.y;
		pgid[i].z = (PN + pblk.z - 1) / pblk.z;

	}
#pragma omp barrier
#pragma omp parallel for
	for(int i = 0; i < gpuNum; ++i)
	{
		cudaSetDevice(i);
		MultiSlices_DDPROJ_ker<<<pgid[i],pblk, 0, stream[i]>>>
				(texObj1[i],texObj2[i],
				thrust::raw_pointer_cast(&d_prj[i][0]),
				make_float2(x0,y0),
				thrust::raw_pointer_cast(&d_cossin[i][0]),
				thrust::raw_pointer_cast(&d_xds[i][0]),
				thrust::raw_pointer_cast(&d_yds[i][0]),
				thrust::raw_pointer_cast(&d_bxds[i][0]),
				thrust::raw_pointer_cast(&d_byds[i][0]),
				make_float2(objCntIdxX, objCntIdxY),
				dx,XN,YN,sSLN[i],DNU,PN);
		h_prj[i] = d_prj[i];
	}
#pragma omp barrier

	DD2::combineProjection(h_prj,hprj, SLN, DNU, PN, sSLN, gpuNum);

#pragma omp parallel for
	for(int i = 0; i < gpuNum; ++i)
	{
		cudaSetDevice(i);
		subVol[i].clear();
		subCossin[i].clear();
		cudaStreamDestroy(stream[i]);

		cudaFreeArray(d_volumeArray1[i]);
		cudaFreeArray(d_volumeArray2[i]);
		cudaDestroyTextureObject(texObj1[i]);
		cudaDestroyTextureObject(texObj2[i]);

		d_prj[i].clear();
		d_xds[i].clear();
		d_yds[i].clear();
		d_bxds[i].clear();
		d_byds[i].clear();
		d_cossin[i].clear();
		h_prj[i].clear();
	}
#pragma omp barrier
	//Clear the data
	delete[] bxds;
	delete[] byds;
	startIdx.clear();
	endIdx.clear();
	sSLN.clear();
	subVol.clear();
	subCossin.clear();
	cossin.clear();
	stream.clear();
	siz.clear();
	nsiz_ZXY.clear();
	nsiz_ZYX.clear();
	SATZXY.clear();
	SATZYX.clear();
	volumeSize1.clear();
	volumeSize2.clear();
	channelDesc1.clear();
	channelDesc2.clear();
	d_volumeArray1.clear();
	d_volumeArray2.clear();
	copyParams1.clear();
	copyParams2.clear();
	resDesc1.clear();
	resDesc2.clear();
	texDesc1.clear();
	texDesc2.clear();
	texObj1.clear();
	texObj2.clear();

	d_prj.clear();
	d_xds.clear();
	d_yds.clear();
	d_bxds.clear();
	d_byds.clear();
	d_cossin.clear();
	h_prj.clear();
	pgid.clear();
}



__global__ void MultiSlices_DDBACK_ker(
	cudaTextureObject_t prjTexObj,
	float* vol,
	const byte* __restrict__ msk,
	const float2* __restrict__ cossin,
	float2 s,
	float S2D,
	float2 curvox, // imgCenter index
	float dx, float dbeta, float detCntIdx,
	int2 VN, int SLN, int PN)
{
	int3 id;
	id.z = threadIdx.x + __umul24(blockIdx.x, blockDim.x);
	id.x = threadIdx.y + __umul24(blockIdx.y, blockDim.y);
	id.y = threadIdx.z + __umul24(blockIdx.z, blockDim.z);
	if(id.z < SLN && id.x < VN.x && id.y < VN.y)
	{
		if(msk[id.y * VN.x + id.x] != 1)
		{
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
		for(int angIdx = 0; angIdx < PN; ++angIdx)
		{
			cossinT = cossin[angIdx * SLN + id.z];
			cursour = make_float2(
					s.x * cossinT.x - s.y * cossinT.y,
					s.x * cossinT.y + s.y * cossinT.x);

			dir = curvox - cursour;

			l_square = dir.x * dir.x + dir.y * dir.y;

			l = rsqrtf(l_square); // 1 / sqrt(l_square);
			alpha = asinf((cursour.y * dir.x - cursour.x * dir.y) * inv_sid * l);

			if(fabsf(cursour.x) > fabsf(cursour.y))
			{
				ddv = dir.x;
			}
			else
			{
				ddv = dir.y;
			}

			deltaAlpha  = ddv / l_square * dx * 0.5;
			cosVal = dx / ddv * sqrtf(l_square);

			idxL = (alpha - deltaAlpha) * dbeta + detCntIdx + 1.0;
			idxR = (alpha + deltaAlpha) * dbeta + detCntIdx + 1.0;

			summ += (tex3D<float>(prjTexObj,id.z + 0.5, idxR, angIdx + 0.5) -
					 tex3D<float>(prjTexObj,id.z + 0.5, idxL, angIdx + 0.5)) * cosVal;
		}
		__syncthreads();
		vol[(id.y * VN.x + id.x) * SLN + id.z] = summ;

	}
}



void MultiSlices_DDBACK(
		float* hvol, // the pointer to the image
		float* hprj, // the pointer to the projection (SLN, DNU, PN) order
		const float x0, const float y0, //position of the initial source
		float* xds, float* yds, // distribution of the detector cells
		const int DNU, // Number of detector cells
		const int SLN, // Number of slices to be projected or backprojected
		const float imgXCenter, const float imgYCenter, //Center of the image
		const int XN, const int YN, // pixel number of the image
		const float dx, // size of the pixel
		float* h_angs, // view angles SIZE SHOULD BE (SLN x PN)
		int PN, // # of view angles
		byte* mask,
		int* startidx,
		const int gpuNum)
{
	std::vector<int> startIdx(startidx, startidx + gpuNum);
	std::vector<int> endIdx(startIdx.size());
	std::copy(startIdx.begin() + 1, startIdx.end(), endIdx.begin());
	endIdx[gpuNum - 1] = SLN;
	std::vector<int> sSLN(startIdx.size());// = endIdx - startIdx;
	for(int i = 0; i < gpuNum; ++i)
	{
		sSLN[i] = endIdx[i] - startIdx[i];
	}
	startIdx.clear();
	endIdx.clear();


	const float2 objCntIdx(
		make_float2((XN - 1.0) * 0.5 - imgXCenter / dx,
			(YN - 1.0) * 0.5 - imgYCenter / dx));

	const float2 sour(make_float2(x0, y0));

	const float S2D = hypotf(xds[0] - x0, yds[0] - y0);
	const float dbeta = atanf(
			(sqrt(powf(xds[1] - xds[0],2.0) + powf(yds[1] - yds[0],2.0)))
			/ S2D * 0.5f) * 2.0f;

	float* bxds = new float[DNU + 1];
	float* byds = new float[DNU + 1];
	DD3Boundaries(DNU+1, xds, bxds);
	DD3Boundaries(DNU+1, yds, byds);
	//Calculate the most left angle
	const float detCntIdx =  fabsf(atanf(bxds[0] / (y0 - byds[0]))) / dbeta - 0.5f;
	delete[] bxds;
	delete[] byds;
	/////////////////////////////////////////////////////////////////////////////

	thrust::host_vector<float2> h_cossin(SLN * PN);
	thrust::transform(h_angs, h_angs + PN * SLN,
			h_cossin.begin(), [=](float ang)
			{return make_float2(cosf(ang),sinf(ang));});

	thrust::host_vector<thrust::host_vector<float> > subProj(gpuNum);
	thrust::host_vector<thrust::host_vector<float2> > subCossin(gpuNum);

	thrust::host_vector<thrust::device_vector<byte> > d_msk(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > d_vol(gpuNum);
	thrust::host_vector<thrust::host_vector<float> > h_vol(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > d_subProj(gpuNum);
	thrust::host_vector<thrust::device_vector<float2> > d_subCossin(gpuNum);
	thrust::host_vector<cudaArray*> d_prjArray(gpuNum);
	thrust::host_vector<cudaTextureObject_t> texObj(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > d_prjSAT(gpuNum);
	thrust::host_vector<cudaStream_t> stream(gpuNum);

	thrust::host_vector<int> siz(gpuNum);
	thrust::host_vector<int> nsiz(gpuNum);

	thrust::host_vector<cudaExtent> prjSize(gpuNum);
	thrust::host_vector<cudaChannelFormatDesc> channelDesc(gpuNum);

	dim3 copyBlk(64,16,1);
	thrust::host_vector<dim3> copyGid(gpuNum);
#pragma omp parallel for
	for(int i = 0; i < gpuNum; ++i)
	{
		cudaSetDevice(i);
		cudaStreamCreate(&stream[i]);

		subProj[i].resize(sSLN[i] * DNU * PN);
		subCossin[i].resize(sSLN[i] * PN);

		d_msk[i].resize(XN * YN);
		thrust::copy(mask, mask + XN * YN, d_msk[i].begin());

		d_vol[i].resize(sSLN[i] * XN * YN);
		h_vol[i].resize(sSLN[i] * XN * YN);

		d_subProj[i].resize(sSLN[i] * DNU * PN);
		d_subCossin[i].resize(sSLN[i] * PN);

		d_prjSAT[i].resize(sSLN[i] * (DNU + 1) * PN);
	}
	// Split the projection
	DD2::splitProjection(subProj, subCossin, hprj, h_cossin, SLN, DNU,
			PN, sSLN, gpuNum);
	h_cossin.clear();
#pragma omp parallel for
	for(int i = 0; i < gpuNum; ++i)
	{
		cudaSetDevice(i);
		d_subProj[i] = subProj[i];
		subProj[i].clear();
		d_subCossin[i] = subCossin[i];
		subCossin[i].clear();

		siz[i] = DNU * sSLN[i] * PN;
		nsiz[i] = (DNU + 1) * sSLN[i] * PN;

		copyGid[i].x = (sSLN[i] + copyBlk.x - 1) / copyBlk.x;
		copyGid[i].y = (DNU + copyBlk.y - 1) / copyBlk.y;
		copyGid[i].z = (PN + copyBlk.z - 1) / copyBlk.z;

		DD2::addOneSidedZeroBoarder_multiSlice_Fan << <copyGid[i], copyBlk, 0, stream[i] >> >(
			thrust::raw_pointer_cast(&d_subProj[i][0]),
			thrust::raw_pointer_cast(&d_prjSAT[i][0]),
			DNU, sSLN[i], PN);

		copyGid[i].x = (sSLN[i] + copyBlk.x - 1) / copyBlk.x;
		copyGid[i].y = (PN + copyBlk.y - 1) / copyBlk.y;
		copyGid[i].z = 1;

		DD2::heorizontalIntegral_multiSlice_Fan << <copyGid[i], copyBlk, 0, stream[i] >> >(
			thrust::raw_pointer_cast(&d_prjSAT[i][0]), DNU + 1, sSLN[i], PN);
		d_subProj[i].clear();

		/////////////////////////////////////////////////////////////////
		prjSize[i].width = sSLN[i];
		prjSize[i].height=  DNU + 1;
		prjSize[i].depth = PN;

		channelDesc[i] = cudaCreateChannelDesc<float>();
		cudaMalloc3DArray(&d_prjArray[i], &channelDesc[i], prjSize[i]);

		cudaMemcpy3DParms copyParams = { 0 };
		copyParams.srcPtr = make_cudaPitchedPtr(
			(void*) thrust::raw_pointer_cast(&d_prjSAT[i][0]),
			prjSize[i].width * sizeof(float),
			prjSize[i].width, prjSize[i].height);
		copyParams.dstArray = d_prjArray[i];
		copyParams.extent = prjSize[i];
		copyParams.kind = cudaMemcpyDeviceToDevice;
		cudaMemcpy3D(&copyParams);
		cudaResourceDesc resDesc;
		memset(&resDesc, 0, sizeof(resDesc));
		resDesc.resType = cudaResourceTypeArray;
		resDesc.res.array.array = d_prjArray[i];
		cudaTextureDesc texDesc;
		memset(&texDesc, 0, sizeof(texDesc));
		texDesc.addressMode[0] = cudaAddressModeClamp;
		texDesc.addressMode[1] = cudaAddressModeClamp;
		texDesc.addressMode[2] = cudaAddressModeClamp;
		texDesc.filterMode = cudaFilterModeLinear;
		texDesc.readMode = cudaReadModeElementType;
		texDesc.normalizedCoords = false;

		cudaCreateTextureObject(&texObj[i], &resDesc, &texDesc, nullptr);

		d_prjSAT[i].clear();
	}
#pragma omp barrier
	subProj.clear();
	d_prjSAT.clear();
	subCossin.clear();
	d_subProj.clear();

#pragma omp parallel for
	for(int i = 0; i < gpuNum; ++i)
	{
		cudaSetDevice(i);
		dim3 blk(BLKX,BLKY,BLKZ);
		dim3 gid(
				(sSLN[i] + blk.x - 1) / blk.x,
				(XN + blk.y - 1) / blk.y,
				(YN + blk.z - 1) / blk.z);
		MultiSlices_DDBACK_ker<< <gid, blk, 0, stream[i] >> >(texObj[i],
			thrust::raw_pointer_cast(&d_vol[i][0]),
			thrust::raw_pointer_cast(&d_msk[i][0]),
			thrust::raw_pointer_cast(&d_subCossin[i][0]),
			sour, S2D, objCntIdx,
			dx, dbeta, detCntIdx, make_int2(XN, YN), sSLN[i], PN);
		h_vol[i] = d_vol[i];
		d_vol[i].clear();

		d_msk[i].clear();
		d_subCossin[i].clear();

		cudaDestroyTextureObject(texObj[i]);
		cudaFreeArray(d_prjArray[i]);
		cudaStreamDestroy(stream[i]);

	}
#pragma omp barrier

	d_vol.clear();
	d_msk.clear();
	d_subCossin.clear();
	d_prjArray.clear();
	texObj.clear();
	stream.clear();
	siz.clear();
	nsiz.clear();
	prjSize.clear();
	channelDesc.clear();
	copyGid.clear();


	thrust::host_vector<int> sSLNn = sSLN;
	DD2::combineVolume(h_vol, hvol, SLN, XN, YN, sSLNn, gpuNum);
	sSLNn.clear();
	sSLN.clear();

}


















__global__ void MultiSlices_PDPROJ_ker(
		cudaTextureObject_t texObj,
		float* proj,
		float2 s,
		float* __restrict__ xds,
		float* __restrict__ yds,
		float2* __restrict__ cossin, // size should be SLN * PN
		float2 objCntIdx,
		float dx, int SLN, int DNU, int PN, int XN, int YN)
{
	int slnIdx = threadIdx.x + blockIdx.x * blockDim.x;
	int detIdx = threadIdx.y + blockIdx.y * blockDim.y;
	int angIdx = threadIdx.z + blockIdx.z * blockDim.z;
	if(slnIdx < SLN && detIdx < DNU && angIdx < PN)
	{
		float2 cssv = cossin[angIdx * SLN + slnIdx];
		float2 cursour = make_float2(
				s.x * cssv.x - s.y * cssv.y,
				s.x * cssv.y + s.y * cssv.x);

		float summ = xds[detIdx];
		float obj = yds[detIdx];

		float idx = 0;
		float2 curDet = make_float2(
				summ * cssv.x - obj * cssv.y,
				summ * cssv.y + obj * cssv.x);
		float2 dir = normalize(curDet - cursour);


		summ = 0;
		obj = 0;
		if(fabs(cssv.x) <= fabs(cssv.y))
		//if(fabsf(dir.y) <= fabsf(dir.x))
		{
			summ = 0;
#pragma unroll
			for(int ii = 0; ii < XN; ++ii)
			{
				obj = (ii - objCntIdx.x) * dx;
				idx = (obj - curDet.x) / dir.x * dir.y + curDet.y;
				idx = idx / dx + objCntIdx.y + 0.5f;
				summ += tex3D<float>(texObj, slnIdx + 0.5f, ii + 0.5f, idx);
			}
			__syncthreads();
			proj[(angIdx * DNU + detIdx) * SLN + slnIdx] = summ * dx / fabsf(dir.x);
		}
		else
		{

			summ = 0;
#pragma unroll
			for(int jj = 0; jj < YN; ++jj)
			{
				obj = (jj - objCntIdx.y) * dx;
				idx = (obj - curDet.y) / dir.y * dir.x + curDet.x;
				idx = idx / dx + objCntIdx.x + 0.5f;
				summ += tex3D<float>(texObj, slnIdx + 0.5f, idx, jj + 0.5f);
			}
			__syncthreads();
			proj[(angIdx * DNU + detIdx) * SLN + slnIdx] = summ * dx / fabsf(dir.y);
		}
	}
}


void MultiSlices_PDPROJ(
		float* hvol, // the pointer to the image
		float* hprj, // the pointer to the projection (SLN, DNU, PN) order
		const float x0, const float y0, //position of the initial source
		float* xds, float* yds, // distribution of the detector cells
		const int DNU, // Number of detector cells
		const int SLN, // Number of slices to be projected or backprojected
		const float imgXCenter, const float imgYCenter, //Center of the image
		const int XN, const int YN, // pixel number of the image
		const float dx, // size of the pixel
		float* h_angs, // view angles SHOULD BE WITH SIZE SLN * PN
		int PN, // # of view angles
		byte* mask,
		int* startIdx, // This means how many slices will be applied to one GPU
		const int gpuNum)
{
	thrust::host_vector<float> hangs(h_angs,h_angs + PN * SLN);

	// Regular the image volume
	for(int i = 0; i != XN * YN; ++i)
	{
		byte v = mask[i];
		for(int z = 0; z != SLN; ++z)
		{
			hvol[i * SLN + z] *= v;
		}
	}

	const float objCntIdxX = (XN - 1.0) * 0.5 - imgXCenter / dx;
	const float objCntIdxY = (YN - 1.0) * 0.5 - imgYCenter / dx;

	//We do not need the overlapping for projection
	std::vector<int> ObjIdx_Start(gpuNum, -1);
	std::vector<int> ObjIdx_End(gpuNum, -1);

	std::vector<std::vector<float> > subVol(gpuNum);


	std::vector<int> sSLN(gpuNum,0);
	for(int i = 1; i != gpuNum; ++i)
	{
		sSLN[i-1] = startIdx[i] - startIdx[i-1];
	}
	sSLN[gpuNum-1] = SLN - startIdx[gpuNum-1];


	std::vector<cudaStream_t> stream(gpuNum);
	std::vector<cudaExtent> volumeSize(gpuNum);
	thrust::host_vector<thrust::host_vector<float2> > subCossin(gpuNum);
	std::vector<int> siz(gpuNum, 0);

	for(int i = 0; i != gpuNum; ++i)
	{
		cudaSetDevice(i);
		cudaStreamCreate(&stream[i]); // Generate multiple streams
		siz[i] = XN * YN * sSLN[i];
		subCossin[i].resize(sSLN[i] * PN);

	}

//	// precalculate the cossin value
	thrust::host_vector<float2> hcossin(PN * SLN);
	thrust::transform(h_angs, h_angs + PN * SLN,
			hcossin.begin(),[=](float ag){return make_float2(cosf(ag),sinf(ag));});
	// Split the volume
	DD2::splitVolume(subVol, subCossin, hvol, hcossin, SLN, XN, YN, PN, sSLN, gpuNum);

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

	std::vector<cudaArray*> d_volumeArray(gpuNum);

	thrust::host_vector<thrust::device_vector<float> > d_vol(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > d_prj(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > d_xds(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > d_yds(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > d_angs(gpuNum);
	thrust::host_vector<thrust::device_vector<float2> >d_cossin(gpuNum);

	dim3 blk(64,16,1);
	std::vector<dim3> gid(gpuNum);
	std::vector<cudaTextureObject_t> texObj(gpuNum);

	omp_set_num_threads(gpuNum);
#pragma omp parallel for
	for(int i = 0; i < gpuNum; ++i)
	{
		cudaSetDevice(i);
		d_vol[i].resize(siz[i]);
		d_vol[i] = subVol[i];
		subVol[i].clear();

		volumeSize[i].width = sSLN[i];
		volumeSize[i].height= XN;
		volumeSize[i].depth = YN;

		cudaMalloc3DArray(&d_volumeArray[i], &channelDesc, volumeSize[i]);

		cudaMemcpy3DParms copyParams = { 0 };

		copyParams.srcPtr = make_cudaPitchedPtr((void*)
			thrust::raw_pointer_cast(&d_vol[i][0]),
			volumeSize[i].width * sizeof(float),
			volumeSize[i].width, volumeSize[i].height);
		copyParams.dstArray = d_volumeArray[i];
		copyParams.extent = volumeSize[i];
		copyParams.kind = cudaMemcpyDeviceToDevice;

		cudaMemcpy3DAsync(&copyParams, stream[i]);
		d_vol[i].clear();

		cudaResourceDesc resDesc;
		memset(&resDesc, 0, sizeof(resDesc));
		resDesc.resType = cudaResourceTypeArray;
		resDesc.res.array.array = d_volumeArray[i];

		cudaTextureDesc texDesc;
		memset(&texDesc, 0, sizeof(texDesc));
		texDesc.addressMode[0] = cudaAddressModeBorder;
		texDesc.addressMode[1] = cudaAddressModeBorder;
		texDesc.addressMode[2] = cudaAddressModeBorder;
		texDesc.filterMode = cudaFilterModeLinear;
		texDesc.readMode = cudaReadModeElementType;
		texDesc.normalizedCoords = false;
		texObj[i] = 0;
		cudaCreateTextureObject(&texObj[i], &resDesc, &texDesc, nullptr);

		d_prj[i].resize(DNU * PN * sSLN[i]);
		d_xds[i].resize(DNU);
		d_yds[i].resize(DNU);
		thrust::copy(xds, xds + DNU, d_xds[i].begin());
		thrust::copy(yds, yds + DNU, d_yds[i].begin());

		//d_angs[i].resize(PN * SLN);
		//thrust::copy(hangs.begin(), hangs.end(), d_angs[i].begin());

		d_cossin[i].resize(PN * sSLN[i]);
		d_cossin[i] = subCossin[i];
		//thrust::transform(d_angs[i].begin(), d_angs[i].end(),
				//d_cossin[i].begin(), DD2::CosSinFunctor());
		//d_angs[i].clear();

		gid[i].x = (sSLN[i] + blk.x - 1) / blk.x;
		gid[i].y = (DNU + blk.y - 1) / blk.y;
		gid[i].z = (PN + blk.z - 1) / blk.z;

	}



	thrust::host_vector<thrust::host_vector<float> > h_prj(gpuNum);
	// Projection process
	 omp_set_num_threads(gpuNum);
#pragma omp parallel for
	for(int i = 0; i < gpuNum; ++i)
	{
		cudaSetDevice(i);
		MultiSlices_PDPROJ_ker<<<gid[i],blk, 0, stream[i]>>>(texObj[i],
				thrust::raw_pointer_cast(&d_prj[i][0]),
				make_float2(x0, y0),
				thrust::raw_pointer_cast(&d_xds[i][0]),
				thrust::raw_pointer_cast(&d_yds[i][0]),
				thrust::raw_pointer_cast(&d_cossin[i][0]),
				make_float2(objCntIdxX,objCntIdxY), dx,
				sSLN[i], DNU, PN, XN, YN);
		h_prj[i].resize(sSLN[i] * DNU * PN);
		h_prj[i] = d_prj[i];
	}

#pragma omp barrier

	DD2::combineProjection(h_prj, hprj, SLN, DNU, PN, sSLN, gpuNum);

	// Clean the resources
	omp_set_num_threads(gpuNum);
#pragma omp parallel for
	for(int i = 0; i < gpuNum; ++i)
	{
		cudaSetDevice(i);
		subVol[i].clear();
		cudaStreamDestroy(stream[i]);
		cudaFreeArray(d_volumeArray[i]);
		d_vol[i].clear();
		d_prj[i].clear();
		d_xds[i].clear();
		d_yds[i].clear();
		d_angs[i].clear();
		d_cossin[i].clear();
		h_prj[i].clear();
	}
	hangs.clear();
	ObjIdx_Start.clear();
	ObjIdx_End.clear();
	sSLN.clear();
	subVol.clear();
	stream.clear();
	volumeSize.clear();
	d_volumeArray.clear();
	d_vol.clear();
	d_prj.clear();
	d_xds.clear();
	d_yds.clear();
	d_angs.clear();
	d_cossin.clear();
	gid.clear();
	texObj.clear();
	siz.clear();
	h_prj.clear();

}

__global__ void MultiSlices_PDBACK_ker(
		cudaTextureObject_t texObj, // projection texture
		float* vol,
		const byte* __restrict__ msk,
		const float2* __restrict__ cossin, // size should be SLN * PN
		float2 s, // source position
		float S2D,
		float2 objCntIdx,
		float dx,
		float dbeta, /// what is dbeta
		float detCntIdx,
		int SLN, int XN, int YN, int DNU, int PN)
{
	int slnIdx = threadIdx.x + blockIdx.x * blockDim.x;
	int xIdx = threadIdx.y + blockIdx.y * blockDim.y;
	int yIdx = threadIdx.z + blockIdx.z * blockDim.z;
	if(slnIdx < SLN && xIdx < XN && yIdx < YN)
	{
		if(msk[yIdx * XN + xIdx] != 1)
			return;

		float2 curVox =
				make_float2((xIdx - objCntIdx.x) * dx, (yIdx - objCntIdx.y) * dx);

		float2 dir;
		float2 cursour;
		float invsid = rsqrtf(s.x * s.x + s.y * s.y);

		float invl;
		//float idxZ;
		float idxXY;
		float alpha;
		float cosVal;

		float2 cossinT;
		float summ = 0;

		float tempVal;

		dbeta = 1.0 / dbeta;
		for(int angIdx = 0; angIdx != PN; ++angIdx)
		{
			cossinT = cossin[angIdx * SLN + slnIdx];
			cursour = make_float2(
				s.x * cossinT.x - s.y * cossinT.y,
				s.x * cossinT.y + s.y * cossinT.x);

			dir = curVox - cursour;
			tempVal = dir.x * dir.x + dir.y * dir.y;

			invl = rsqrtf(tempVal);

			alpha = asinf((cursour.y * dir.x - cursour.x * dir.y) * invl * invsid);
			if(fabsf(cursour.x) >= fabsf(cursour.y))
			{
				cosVal = fabsf(1.0f / dir.x);
			}
			else
			{
				cosVal = fabsf(1.0f / dir.y);
			}
			cosVal *= (dx * sqrtf(tempVal));
			idxXY = alpha * dbeta + detCntIdx + 0.5;
			summ += tex3D<float>(texObj, slnIdx + 0.5f,
					idxXY, angIdx + 0.5f) * cosVal;
		}
		__syncthreads();
		vol[(yIdx * XN + xIdx) * SLN + slnIdx] = summ;
	}
}


void MultiSlices_PDBACK(
		float* hvol, // the pointer to the image
		float* hprj, // the pointer to the projection (SLN, DNU, PN) order
		const float x0, const float y0, //position of the initial source
		float* xds, float* yds, // distribution of the detector cells
		const int DNU, // Number of detector cells
		const int SLN, // Number of slices to be projected or backprojected
		const float imgXCenter, const float imgYCenter, //Center of the image
		const int XN, const int YN, // pixel number of the image
		const float dx, // size of the pixel
		float* h_angs, // view angles SHOULD BE WITH SIZE SLN*PN
		int PN, // # of view angles
		byte* mask,
		int* startIdx,
		const int gpuNum)
{
	//Set the start and end slices for each GPU
	thrust::host_vector<int> ObjZIdx_Start(startIdx, startIdx + gpuNum);
	thrust::host_vector<int> ObjZIdx_End(ObjZIdx_Start.size());

	std::copy(ObjZIdx_Start.begin() + 1, ObjZIdx_Start.end(), ObjZIdx_End.begin());
	ObjZIdx_End[gpuNum - 1] = SLN;

	float* bxds = new float[DNU + 1];
	float* byds = new float[DNU + 1];

	DD3Boundaries(DNU + 1, xds, bxds);
	DD3Boundaries(DNU + 1, yds, byds);

	float2 dir = normalize(make_float2(-x0, -y0));
	float2 dirL = normalize(make_float2(bxds[0] - x0, byds[0] - y0));
	float2 dirR = normalize(make_float2(bxds[DNU] - x0, byds[DNU] - y0));
	float dbeta = asin(dirL.x * dirR.y - dirL.y * dirR.x) / DNU;
	float minBeta = asin(dir.x * dirL.y - dir.y * dirL.x);
	float detCntIdx = -minBeta / dbeta - 0.5;
	const float S2D = hypotf(xds[0] - x0, yds[0] - y0);

	delete[] bxds;
	delete[] byds;

	thrust::host_vector<int> sSLN = ObjZIdx_End - ObjZIdx_Start;

	const float objCntIdxX = (XN - 1.0f) * 0.5f - imgXCenter / dx;
	const float objCntIdxY = (YN - 1.0f) * 0.5f - imgYCenter / dx;

	thrust::host_vector<float2> sour(gpuNum);
	thrust::host_vector<thrust::device_vector<byte> > msk(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > vol(gpuNum);
	thrust::host_vector<thrust::device_vector<float> > prj(gpuNum);
	thrust::host_vector<thrust::device_vector<float2> > cossin(gpuNum);
	thrust::host_vector<cudaArray*> d_prjArray(gpuNum);
	thrust::host_vector<cudaTextureObject_t> texObj(gpuNum);
	thrust::host_vector<cudaStream_t> stream(gpuNum);

	thrust::host_vector<thrust::host_vector<float> > host_vol(gpuNum);
	dim3 blk(32,16,1);
	thrust::host_vector<dim3> gid(gpuNum);

	// precalculate the cossin value
	thrust::host_vector<float2> hcossin(PN * SLN);
	thrust::transform(h_angs, h_angs + PN * SLN,
			hcossin.begin(),[=](float ag){return make_float2(cosf(ag),sinf(ag));});
	thrust::host_vector<thrust::host_vector<float2> > subCossin(gpuNum);

	//Split the projection data
	thrust::host_vector<thrust::host_vector<float> > sbprj(gpuNum);

	for(int i = 0 ; i != gpuNum; ++i)
	{
		sbprj.resize(sSLN[i] * DNU * PN);
	}

	DD2::splitProjection(sbprj,subCossin, hprj, hcossin, SLN, DNU, PN, sSLN, gpuNum);

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	omp_set_num_threads(gpuNum);
#pragma omp parallel for
	for(int i = 0; i < gpuNum; ++i)
	{
		cudaSetDevice(i);
		cudaStreamCreate(&stream[i]);

		msk[i].resize(XN * YN);
		thrust::copy(mask, mask + XN * YN, msk[i].begin());
		vol[i].resize(sSLN[i] * XN * YN);
		prj[i].resize(sSLN[i] * DNU * PN);
		prj[i] = sbprj[i];

		cudaExtent prjSize;
		prjSize.width = sSLN[i];
		prjSize.height = DNU;
		prjSize.depth = PN;

		cudaMalloc3DArray(&d_prjArray[i], &channelDesc, prjSize);

		cudaMemcpy3DParms copyParams = {0};
		copyParams.srcPtr = make_cudaPitchedPtr(
			(void*) thrust::raw_pointer_cast(&prj[i][0]),
			prjSize.width * sizeof(float),
			prjSize.width, prjSize.height);
		copyParams.dstArray = d_prjArray[i];
		copyParams.extent = prjSize;
		copyParams.kind = cudaMemcpyDeviceToDevice;

		cudaMemcpy3DAsync(&copyParams, stream[i]);

		cudaResourceDesc resDesc;
		memset(&resDesc, 0, sizeof(resDesc));
		resDesc.resType = cudaResourceTypeArray;
		resDesc.res.array.array = d_prjArray[i];
		cudaTextureDesc texDesc;
		memset(&texDesc, 0, sizeof(texDesc));
		texDesc.addressMode[0] = cudaAddressModeBorder;
		texDesc.addressMode[1] = cudaAddressModeBorder;
		texDesc.addressMode[2] = cudaAddressModeBorder;
		texDesc.filterMode = cudaFilterModeLinear;
		texDesc.readMode = cudaReadModeElementType;
		texDesc.normalizedCoords = false;
		cudaCreateTextureObject(&texObj[i], &resDesc, &texDesc, nullptr);
		prj[i].clear();

		cossin[i].resize(PN * sSLN[i]);
		cossin[i] = subCossin[i];

		gid[i].x = (sSLN[i] + blk.x - 1) / blk.x;
		gid[i].y = (DNU + blk.y - 1) / blk.y;
		gid[i].z = (PN + blk.z) / blk.z;
	}

#pragma omp parallel for
	for(int i = 0; i < gpuNum; ++i)
	{
		cudaSetDevice(i);

	    MultiSlices_PDBACK_ker<<<gid[i], blk, 0, stream[i]>>>(
	    	texObj[i],
	    	thrust::raw_pointer_cast(&vol[i][0]),
	    	thrust::raw_pointer_cast(&msk[i][0]),
	    	thrust::raw_pointer_cast(&cossin[i][0]),
	    	make_float2(x0,y0), S2D, make_float2(objCntIdxX,objCntIdxY),
	    	dx, dbeta, detCntIdx, sSLN[i], XN, YN, DNU, PN);
	    host_vol[i].resize(sSLN[i] * XN * YN);
	    host_vol[i] = vol[i];

	}
#pragma omp barrier

	//combine the volume
	DD2::combineVolume(host_vol, hvol, SLN, XN, YN, sSLN, gpuNum);

#pragma omp parallel for
	for(int i = 0; i < gpuNum; ++i)
	{
		cudaSetDevice(i);
		msk[i].clear();
		vol[i].clear();
		prj[i].clear();
		cossin[i].clear();
		cudaDestroyTextureObject(texObj[i]);
		cudaFreeArray(d_prjArray[i]);
		cudaStreamDestroy(stream[i]);
		host_vol[i].clear();
		sbprj[i].clear();
	}

	ObjZIdx_Start.clear();
	ObjZIdx_End.clear();
	sSLN.clear();
	sour.clear();
	msk.clear();
	vol.clear();
	prj.clear();
	cossin.clear();
	d_prjArray.clear();
	texObj.clear();
	stream.clear();
	host_vol.clear();
	gid.clear();
	hcossin.clear();
	//hangs.clear();
	sbprj.clear();
}

extern "C"
void DD2_multiGPU(
		float* hvol, // the pointer to the image
		float* hprj, // the pointer to the projection (SLN, DNU, PN) order
		const int method, // Control to use forward projection or backprojection
		const float x0, const float y0, //position of the initial source
		float* xds, float* yds, // distribution of the detector cells
		const int DNU, // Number of detector cells
		const int SLN, // Number of slices to be projected or backprojected
		const float imgXCenter, const float imgYCenter, //Center of the image
		const int XN, const int YN, // pixel number of the image
		const float dx, // size of the pixel
		float* hangs, // view angles
		int PN, // # of view angles
		byte* mask,
		int* startIdx,
		const int gpuNum)
{
	switch(method)
	{
	case 0: // DD projection
		MultiSlices_DDPROJ(hvol, hprj, x0, y0, xds, yds, DNU, SLN,
			imgXCenter, imgYCenter,XN, YN, dx, hangs, PN, mask,
			startIdx,gpuNum);
		break;
	case 1: // DD backprojection
		MultiSlices_DDBACK(hvol, hprj, x0, y0, xds, yds, DNU, SLN,
			imgXCenter, imgYCenter,XN, YN, dx, hangs, PN, mask,
			startIdx,gpuNum);
		break;
	case 2: // PD projection
		MultiSlices_PDPROJ(hvol, hprj, x0, y0, xds, yds, DNU, SLN,
			imgXCenter, imgYCenter,XN, YN, dx, hangs, PN, mask,
			startIdx,gpuNum);
		break;
	case 3: // PD backprojection
		MultiSlices_PDBACK(hvol, hprj, x0, y0, xds, yds, DNU, SLN,
			imgXCenter, imgYCenter,XN, YN, dx, hangs, PN, mask,
			startIdx,gpuNum);
		break;
	default:
		break;

	}
}
