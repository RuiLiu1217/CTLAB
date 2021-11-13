/*
 * utilities.cuh
 * The utilities functions that used in GPU and CPU
 *  Created on: Oct 2, 2015
 *  Author: Rui Liu
 */

#ifndef UTILITIES_CUH_
#define UTILITIES_CUH_





typedef unsigned char byte;
typedef thrust::device_vector<float> d_vec_t;
typedef thrust::host_vector<float> h_vec_t;


#ifndef nullptr
#define nullptr NULL
#endif




///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
// Get one sub-volume from the whole volume.
// Assume that the volumes are stored in Z, X, Y order
template<typename T>
void getSubVolume(const T* vol,
		const size_t XN, const size_t YN, const size_t ZN,
		const size_t ZIdx_Start, const size_t ZIdx_End, T* subVol)
{
	const size_t SZN = ZIdx_End - ZIdx_Start;
	for (size_t yIdx = 0; yIdx != YN; ++yIdx)
	{
		for (size_t xIdx = 0; xIdx != XN; ++xIdx)
		{
			for (size_t zIdx = ZIdx_Start; zIdx != ZIdx_End; ++zIdx)
			{
				subVol[(yIdx * XN + xIdx) * SZN + (zIdx - ZIdx_Start)] = vol[(yIdx * XN + xIdx) * ZN + zIdx];
			}
		}
	}
}

template<typename T>
void getSubVolume(const T* vol,
		const size_t XYN, const size_t ZN,
		const size_t ZIdx_Start, const size_t ZIdx_End, T* subVol)
{
	const size_t SZN = ZIdx_End - ZIdx_Start;
	for (size_t xyIdx = 0; xyIdx != XYN; ++xyIdx)
	{
		for (size_t zIdx = ZIdx_Start; zIdx != ZIdx_End; ++zIdx)
		{
			subVol[xyIdx * SZN + (zIdx - ZIdx_Start)] = vol[xyIdx * ZN + zIdx];
		}
	}
}
///////////////////////////////////////////////////////////////////////////////////

// For projection, before we divide the volume into serveral sub-volumes, we have
// to calculate the Z index range
template<typename T>
void getVolZIdxPair(const thrust::host_vector<T>& zPos, // Z position of the source.
		//NOTE: We only assume the spiral CT case that zPos is increasing.
		const size_t PrjIdx_Start, const size_t PrjIdx_End,
		const T detCntIdxV, const T detStpZ, const int DNV,
		const T objCntIdxZ,	const T dz, const int ZN, // Size of the volume
		int& ObjIdx_Start, int& ObjIdx_End) // The end is not included
{
	const T lowerPart = (detCntIdxV + 0.5) * detStpZ;
	const T upperPart = DNV * detStpZ - lowerPart;
	const T startPos = zPos[PrjIdx_Start] - lowerPart;
	const T endPos = zPos[PrjIdx_End - 1] + upperPart;

	ObjIdx_Start = floor((startPos / dz) + objCntIdxZ - 1);
	ObjIdx_End = ceil((endPos / dz) + objCntIdxZ + 1) + 1;

	ObjIdx_Start = (ObjIdx_Start < 0) ? 0 : ObjIdx_Start;
	ObjIdx_Start = (ObjIdx_Start > ZN) ? ZN : ObjIdx_Start;

	ObjIdx_End = (ObjIdx_End < 0) ? 0 : ObjIdx_End;
	ObjIdx_End = (ObjIdx_End > ZN) ? ZN : ObjIdx_End;
}

///////////////////////////////////////////////////////////////////////////////////
// For backprojection, after decide the subvolume range, we have to decide the
// projection range to cover the subvolume.
template<typename T>
void getPrjIdxPair(const thrust::host_vector<T>& zPos, // Z Position of the source.
		// NOTE: we assume that it is pre sorted
		const size_t ObjZIdx_Start, const size_t ObjZIdx_End, // sub vol range,
		// NOTE: the objZIdx_End is not included
		const T objCntIdxZ, const T dz, const int ZN,
		const T detCntIdxV, const T detStpZ, const int DNV,
		int& prjIdx_Start, int& prjIdx_End)
{
	const int PN = zPos.size();

	const T lowerPartV = (ObjZIdx_Start - objCntIdxZ - 0.5) * dz;
	const T highrPartV = lowerPartV + (ObjZIdx_End - ObjZIdx_Start) * dz;

	const T lowerPartDet = (detCntIdxV + 0.5) * detStpZ;
	const T upperPartDet = DNV * detStpZ - lowerPartDet;

	//The source position
	const T sourLPos = lowerPartV - upperPartDet;
	const T sourHPos = highrPartV + lowerPartDet;

	prjIdx_Start = thrust::upper_bound(zPos.begin(),zPos.end(),sourLPos) - zPos.begin() - 1;
	prjIdx_End = thrust::upper_bound(zPos.begin(),zPos.end(),sourHPos) - zPos.begin() + 2;
	prjIdx_Start = (prjIdx_Start < 0) ? 0 : prjIdx_Start;
	prjIdx_Start = (prjIdx_Start > PN)? PN: prjIdx_Start;

	prjIdx_End = (prjIdx_End < 0) ? 0 : prjIdx_End;
	prjIdx_End = (prjIdx_End > PN) ? PN : prjIdx_End;
}


////////////////////////////////////////////////////////////////////////////////////
// The volume is also stored in Z, X, Y order
// Not tested yet.
template<typename T>
void combineVolume(
	T* vol, // The volume to be combined
	const int XN, const int YN, const int ZN,
	T** subVol, // All sub volumes
	const int* SZN, // Number of slices for each subVolume
	const int subVolNum) // Number of sub volumes
{
	int kk = 0;
	for (size_t yIdx = 0; yIdx != YN; ++yIdx)
	{
		for (size_t xIdx = 0; xIdx != XN; ++xIdx)
		{
			kk = 0;
			for (size_t volIdx = 0; volIdx != subVolNum; ++volIdx)
			{
				for (size_t zIdx = 0; zIdx != SZN[volIdx]; ++zIdx)
				{
					vol[(yIdx * XN + xIdx) * ZN + kk] = subVol[volIdx][(yIdx * XN + xIdx) * SZN[volIdx] + zIdx];
					kk = kk + 1;
				}
			}
		}
	}
}


template<typename T>
void combineVolume(
	T* vol, // The volume to be combined
	const int XN, const int YN, const int ZN,
	thrust::host_vector<thrust::host_vector<T> >& subVol, // All sub volumes
	const int* SZN, // Number of slices for each subVolume
	const int subVolNum) // Number of sub volumes
{
	int kk = 0;
	for (size_t yIdx = 0; yIdx != YN; ++yIdx)
	{
		for (size_t xIdx = 0; xIdx != XN; ++xIdx)
		{
			kk = 0;
			for (size_t volIdx = 0; volIdx != subVolNum; ++volIdx)
			{
				for (size_t zIdx = 0; zIdx != SZN[volIdx]; ++zIdx)
				{
					vol[(yIdx * XN + xIdx) * ZN + kk] = subVol[volIdx][(yIdx * XN + xIdx) * SZN[volIdx] + zIdx];
					kk = kk + 1;
				}
			}
		}
	}
}



template<typename T>
__host__ __device__ inline T intersectLength(const T& fixedmin, const T& fixedmax, const T& varimin, const T& varimax)
{
	const T left = (fixedmin > varimin) ? fixedmin : varimin;
	const T right = (fixedmax < varimax) ? fixedmax : varimax;
	return abs(right - left) * static_cast<double>(right > left);
}





/// \brief SIDDON kernel function 5
inline __device__  float dev_alphaU_Fun(const float& d, const float& startx, const float& endx)
{
	if (IS_ZERO<float>(startx - endx))
	{
		return 1000.0f;//(d/1e-6);
	}
	return d / fabsf(startx - endx);
}


/// \brief SIDDON kernel function 6
template<typename T>
inline __device__  int dev_iu_Fun(const T& start, const T& end)
{
	return (start < end) ? 1 : -1;
}


enum ForwardDDMethod{PROJ_BRANCHLESS=0,PROJ_PSEUDODISTANCE=2};
enum BackwardDDMethod{BACK_BRANCHLESS=0,BACK_PSEUDODISTANCE=2,BACK_ZLINEBRANCHLESS=3};


template<typename T>
thrust::host_vector<T> operator-(
	const thrust::host_vector<T>& a,
	const thrust::host_vector<T>& b)
{
	thrust::host_vector<T> res(a);
	thrust::transform(res.begin(),res.end(),b.begin(),res.begin(),[=](T aa, T bb){return aa - bb;});
	return res;
}


template<typename T>
std::vector<T> operator-(
	const std::vector<T>& a,
	const std::vector<T>& b)
{
	std::vector<T> res(a);
	thrust::transform(res.begin(),res.end(),b.begin(),res.begin(),[=](T aa, T bb){return aa - bb;});
	return res;
}

#endif /* UTILITIES_CUH_ */

