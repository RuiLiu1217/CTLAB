/*
 * utilities.cuh
 * The utilities functions that used in GPU and CPU
 *  Created on: Oct 2, 2015
 *  Author: Rui Liu
 */

#ifndef UTILITIES_CUH_
#define UTILITIES_CUH_

#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <thrust/transform.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <thrust/functional.h>
#include <thrust/sequence.h>
#include <thrust/find.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/functional.h>
#include <thrust/sequence.h>
#include <thrust/binary_search.h>


#include <omp.h>

#ifndef __PI__
#define __PI__
#define PI		3.141592653589793
#define PI_2		1.570796326794897
#define PI_4		0.785398163397448
#define PI_3_4		2.356194490192344
#define PI_5_4		3.926990816987241
#define PI_7_4		5.497787143782138
#define TWOPI       6.283185307179586
#endif

#define FORCEINLINE 1
#if FORCEINLINE
#define INLINE __forceinline__
#else
#define INLINE inline
#endif

//#if DEBUG
//#define CUDA_CHECK_RETURN(value) { cudaError_t _m_cudaStat = value; if(_m_cudaStat != cudaSuccess){fprintf(stderr, "Error %s at line $d in file %s\n", cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__); exit(1);}}
//#else
//#define CUDA_CHECK_RETURN(value) {value;}
//#endif

#if DEBUG
#define CUDA_CHECK_RETURN(value) {											\
	cudaError_t _m_cudaStat = value;										\
	if (_m_cudaStat != cudaSuccess) {										\
		fprintf(stderr, "Error %s at line %d in file %s\n",					\
				cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);		\
		exit(1);															\
	} }
// Same function as CUDA_CHECK_RETURN
#define CUDA_SAFE_CALL(call) do{ cudaError_t err = call; if (cudaSuccess != err) {  fprintf (stderr, "Cuda error in file '%s' in line %i : %s.", __FILE__, __LINE__, cudaGetErrorString(err) );  exit(EXIT_FAILURE);  } } while (0)
#else
#define CUDA_CHECK_RETURN(value) {value;}
#define CUDA_SAFE_CALL(value) {value;}
#endif



typedef unsigned char byte;
typedef thrust::device_vector<float> d_vec_t;
typedef thrust::host_vector<float> h_vec_t;


#ifndef nullptr
#define nullptr NULL
#endif

INLINE __host__ __device__ const float2 operator/(const float2& a, float b)
{
	return make_float2(a.x / b, a.y / b);
}

INLINE __host__ __device__ const float3 operator+(const float3& a, const float3& b)
{
	return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

INLINE __host__ __device__ const float3 operator-(const float3& a, const float3& b)
{
	return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

INLINE __host__ __device__ const double3 operator-(const double3& a, const double3& b)
{
	return make_double3(a.x - b.x, a.y - b.y, a.z - b.z);
}
INLINE __host__ __device__ const float2 operator-(const float2& a, const float2& b)
{
	return make_float2(a.x - b.x, a.y - b.y);
}

INLINE __host__ __device__ const float3 operator*(const float3& a, const float3& b)
{
	return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}

INLINE __host__ __device__ const float3 operator*(const float3& a, float b)
{
	return make_float3(a.x * b, a.y * b, a.z * b);
}

INLINE __host__ __device__ const float3 operator/(const float3& a, const float3& b)
{
	return make_float3(a.x / b.x, a.y / b.y, a.z / b.z);
}

INLINE __host__ __device__ const float3 operator/(const float3& a, float b)
{
	return make_float3(a.x / b, a.y / b, a.z / b);
}


INLINE __host__ __device__ const double3 operator/(const double3& a, double b)
{
	return make_double3(a.x / b, a.y / b, a.z / b);
}


INLINE __host__ __device__ const float3 operator-(const int3& a, const float3& b)
{
	return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

INLINE __host__ __device__ float length(const float2& a)
{
	return sqrtf(a.x * a.x + a.y * a.y);
}

INLINE __host__ __device__ float length(const float3& a)
{
	return sqrtf(a.x * a.x + a.y * a.y + a.z * a.z);
}

INLINE __host__ __device__ double length(const double3& a)
{
	return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}


INLINE __host__ __device__ const float2 normalize(const float2& a)
{
	return a / length(a);
}

INLINE __host__ __device__ const float3 normalize(const float3& a)
{
	return a / length(a);
}

INLINE __host__ __device__ const double3 normalize(const double3& a)
{
	return a / length(a);
}

INLINE __host__ __device__ float fminf(const float2& a)
{
	return fminf(a.x, a.y);
}

INLINE __host__ __device__ float fminf(const float3& a)
{
	return fminf(a.x, fminf(a.y, a.z));
}

INLINE __host__ __device__ float fmaxf(const float2& a)
{
	return fmaxf(a.x, a.y);
}

INLINE __host__ __device__ float fmaxf(const float3& a)
{
	return fmaxf(a.x, fmaxf(a.y, a.z));
}

INLINE __host__ __device__ const float3 fminf(const float3& a, const float3& b)
{
	return make_float3(fminf(a.x, b.x), fminf(a.y, b.y), fminf(a.z, b.z));
}

INLINE __host__ __device__ const float3 fmaxf(const float3& a, const float3& b)
{
	return make_float3(fmaxf(a.x, b.x), fmaxf(a.y, b.y), fmaxf(a.z, b.z));
}

INLINE __host__ __device__ const float2 fminf(const float2& a, const float2& b)
{
	return make_float2(fminf(a.x, b.x), fminf(a.y, b.y));
}

INLINE __host__ __device__ const float2 fmaxf(const float2& a, const float2& b)
{
	return make_float2(fmaxf(a.x, b.x), fmaxf(a.y, b.y));
}

template<typename T>
INLINE __host__ __device__ T lerp(T v0, T v1, T t)
{
	return fma(t, v1, fma(-t, v0, v0));
}

template<typename T>
INLINE __host__ __device__ T bilierp(T v0, T v1, T v2, T v3, T t1, T t2)
{
	T vv0 = fma(t1, v1, fma(-t1, v0, v0));
	T vv1 = fma(t1, v3, fma(-t1, v2, v2));
	return fma(t2, vv1, fma(-t2, vv0, vv0));
}


INLINE __device__ double bilerp(int2 v0, int2 v1, int2 v2, int2 v3, float t1, float t2)
{
	double v0_ = __hiloint2double(v0.y, v0.x);
	double v1_ = __hiloint2double(v1.y, v1.x);
	double v2_ = __hiloint2double(v2.y, v2.x);
	double v3_ = __hiloint2double(v3.y, v3.x);

	double vv0 = v0_ * (1.0 - t1) + v1_ * t1;
	double vv1 = v2_ * (1.0 - t1) + v3_ * t1;
	return vv0 * (1 - t2) + vv1 * t2;
}


INLINE __device__ double bilerp(int2 v0, int2 v1, int2 v2, int2 v3, double t1, double t2)
{
	double v0_ = __hiloint2double(v0.y, v0.x);
	double v1_ = __hiloint2double(v1.y, v1.x);
	double v2_ = __hiloint2double(v2.y, v2.x);
	double v3_ = __hiloint2double(v3.y, v3.x);

	double vv0 = v0_ * (1.0 - t1) + v1_ * t1;
	double vv1 = v2_ * (1.0 - t1) + v3_ * t1;
	return vv0 * (1 - t2) + vv1 * t2;
}



INLINE __host__ __device__ bool intersectBox(
	const float3& sour,
	const float3& dir,
	const float3& boxmin,
	const float3& boxmax,
	float* tnear, float* tfar)
{
	const float3 invR = make_float3(1.0 / dir.x, 1.0 / dir.y, 1.0 / dir.z);
	const float3 tbot = invR * (boxmin - sour);
	const float3 ttop = invR * (boxmax - sour);

	const float3 tmin = fminf(ttop, tbot);
	const float3 tmax = fmaxf(ttop, tbot);

	const float largest_tmin = fmaxf(tmin);
	const float smallest_tmax = fminf(tmax);
	*tnear = largest_tmin;
	*tfar = smallest_tmax;
	return smallest_tmax > largest_tmin;
}

template<typename T>
INLINE __host__ __device__ T regularizeAngle(T curang)
{
	T c = curang;
	while (c >= TWOPI){ c -= TWOPI; }
	while (c < 0){ c += TWOPI; }
	return c;
}


INLINE __host__ __device__ void invRotVox(
	const float3& curVox,
	float3& virVox,
	const float2& cossinT,
	const float zP)
{
	virVox.x = curVox.x * cossinT.x + curVox.y * cossinT.y;
	virVox.y =-curVox.x * cossinT.y + curVox.y * cossinT.x;
	virVox.z = curVox.z - zP;
}

INLINE __device__ float3 invRot(
	const float3 inV,
	const float2 cossin,
	const float zP)
{
	float3 outV;
	outV.x = inV.x * cossin.x + inV.y * cossin.y;
	outV.x =-inV.x * cossin.y + inV.y * cossin.x;
	outV.z = inV.z - zP;
	return outV;
}


namespace CTMBIR
{
	struct ConstantForBackProjection3
	{
		float x0;
		float y0;
		float z0;

		typedef thrust::tuple<float, float> InTuple;
		typedef thrust::tuple<float3, float3, float2> OutTuple;
		ConstantForBackProjection3(
			const float _x0, const float _y0, const float _z0) :x0(_x0),
			y0(_y0), z0(_z0){}

		__device__ OutTuple operator()(const InTuple& tp)
		{
			float curang = regularizeAngle(thrust::get<0>(tp));
			float zP = thrust::get<1>(tp);
			float cosT = cosf(curang);
			float sinT = sinf(curang);

			float3 cursour = make_float3(
				x0 * cosT - y0 * sinT,
				x0 * sinT + y0 * cosT,
				z0 + zP);

			float2 dirsour = normalize(make_float2(-cursour.x, -cursour.y));
			return thrust::make_tuple(make_float3(cosT, sinT, zP), cursour, dirsour);

		}

	};

	struct ConstantForBackProjection4{

		float x0;
		float y0;
		float z0;

		typedef thrust::tuple<float, float> InTuple;
		ConstantForBackProjection4(const float _x0, const float _y0, const float _z0)
			: x0(_x0), y0(_y0), z0(_z0){}

		__device__ float3 operator()(const InTuple& tp)
		{
			float curang = regularizeAngle(thrust::get<0>(tp));
			float zP = thrust::get<1>(tp);
			float cosT = cosf(curang);
			float sinT = sinf(curang);
			return make_float3(cosT, sinT, zP);
		}

	};

	template<typename T>
	struct ConstantForBackProjection{

		T x0;
		T y0;
		T z0;

		typedef thrust::tuple<T, T> InTuple;
		ConstantForBackProjection(const T _x0, const T _y0, const T _z0)
			: x0(_x0), y0(_y0), z0(_z0){}

		__device__ float3 operator()(const InTuple& tp)
		{
			T curang = regularizeAngle(thrust::get<0>(tp));
			T zP = thrust::get<1>(tp);
			T cosT = cosf(curang);
			T sinT = sinf(curang);
			return make_float3(cosT, sinT, zP);
		}
	};


	template<>
	struct ConstantForBackProjection<double>{

		double x0;
		double y0;
		double z0;

		typedef thrust::tuple<double, double> InTuple;
		ConstantForBackProjection(const double _x0, const double _y0, const double _z0)
			: x0(_x0), y0(_y0), z0(_z0){}

		__device__ double3 operator()(const InTuple& tp)
		{
			double curang = regularizeAngle(thrust::get<0>(tp));
			double zP = thrust::get<1>(tp);
			double cosT = cos(curang);
			double sinT = sin(curang);
			return make_double3(cosT, sinT, zP);
		}
	};

}


template<typename T>
void DD3Boundaries(int nrBoundaries, T*pCenters, T *pBoundaries)
{
	int i;
	if (nrBoundaries >= 3)
	{
		*pBoundaries++ = 1.5 * *pCenters - 0.5 * *(pCenters + 1);
		for (i = 1; i <= (nrBoundaries - 2); i++)
		{
			*pBoundaries++ = 0.5 * *pCenters + 0.5 * *(pCenters + 1);
			pCenters++;
		}
		*pBoundaries = 1.5 * *pCenters - 0.5 * *(pCenters - 1);
	}
	else
	{
		*pBoundaries = *pCenters - 0.5;
		*(pBoundaries + 1) = *pCenters + 0.5;
	}

}




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
__host__ __device__ inline T intersectLength(const T& fixedmin, const T& fixedmax, const T& varimin, const T& varimax)
{
	const T left = (fixedmin > varimin) ? fixedmin : varimin;
	const T right = (fixedmax < varimax) ? fixedmax : varimax;
	return abs(right - left) * static_cast<double>(right > left);
}





template<typename T>
inline __host__ __device__ bool IS_ZERO(const T& x)
{
	return (abs(x) < 1.0E-9);
}

/// \brief SIDDON kernel function 1
inline __device__ float dev_pFun(const float& alpha, const float& pstart, const float&pend)
{
	return pstart + alpha * (pend - pstart);
}

/// \brief SIDDON kernel function 2
inline __device__ float dev_alpha_IFun(const float& b, const float& d, const float& pstart, const float& pend, const unsigned int& i)
{
	if (!IS_ZERO<float>(pend - pstart))
	{
		return ((b + (float) i*d) - pstart) / (pend - pstart);
	}
	else return 1000;//((b + i*d)-pstart)/(1e-6);
}



/// \brief SIDDON kernel function 3
inline __device__ float dev_varphiFun(const float& alpha, const float& b, const float& d, const float& pstart, const float& pend)
{
	return (dev_pFun(alpha, pstart, pend) - b) / d;
}


/// \brief SIDDON kernel function 4
inline	__device__ void dev_minmaxIdxFun(
	const float& pstart, const float& pend,
	const float& b, const float& d,
	const float& alphaMIN, const float& alphaMAX,
	const float& alphaPmin, const float& alphaPmax,
	const unsigned int& Nplane, int* imin, int* imax)
{
	if (pstart < pend)
	{
		if (IS_ZERO<float>(alphaMIN - alphaPmin))
		{
			*imin = 1;
		}
		else
		{
			*imin = static_cast<int>(ceil(dev_varphiFun(alphaMIN, b, d, pstart, pend)));
		}
		if (IS_ZERO<float>(alphaMAX - alphaPmax))
		{
			*imax = Nplane - 1;
		}
		else
		{
			*imax = static_cast<int>(floor(dev_varphiFun(alphaMAX, b, d, pstart, pend)));
		}
	}
	else
	{
		if (IS_ZERO<float>(alphaMIN - alphaPmin))
		{
			*imax = Nplane - 2;
		}
		else
		{
			*imax = static_cast<int>(floor(dev_varphiFun(alphaMIN, b, d, pstart, pend)));
		}
		if (IS_ZERO<float>(alphaMAX - alphaPmax))
		{
			*imin = 0;
		}
		else
		{
			*imin = static_cast<int>(ceil(dev_varphiFun(alphaMAX, b, d, pstart, pend)));
		}
	}
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



/// \brief SIDDON line integral function in 3D
inline __device__ float calSiddonOneRayKer(
	const float startX, const float startY, const float startZ,
	const float endX, const float endY, const float endZ,
	const float __MINOL__, const float __MINOW__, const float __MINOH__,
	const float __OBJSTPL__, const float __OBJSTPW__, const float __OBJSTPH__,
	const unsigned int __OBJLR__, const unsigned int __OBJWR__,const unsigned int __OBJHR__,
	cudaTextureObject_t volTex, float* totWeight)
{
	float dirX = endX - startX;
	float dirY = endY - startY;
	float dirZ = endZ - startZ;

	const float dconv = sqrtf(dirX * dirX + dirY * dirY + dirZ * dirZ);
	int imin(0), imax(0), jmin(0), jmax(0), kmin(0), kmax(0);

	float alphaxmin = 0.0f;
	float alphaxmax = 0.0f;
	float alphaymin = 0.0f;
	float alphaymax = 0.0f;
	float alphazmin = 0.0f;
	float alphazmax = 0.0f;

	dirX = dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, 0);
	dirY = dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, __OBJLR__);
	if (dirX < dirY)
	{
		alphaxmin = dirX;
		alphaxmax = dirY;
	}
	else
	{
		alphaxmin = dirY;
		alphaxmax = dirX;
	}
	dirX = dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, 0);
	dirY = dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, __OBJWR__);
	if (dirX < dirY)
	{
		alphaymin = dirX;
		alphaymax = dirY;
	}
	else
	{
		alphaymin = dirY;
		alphaymax = dirX;
	}
	dirX = dev_alpha_IFun(__MINOH__, __OBJSTPH__, startZ, endZ, 0);
	dirY = dev_alpha_IFun(__MINOH__, __OBJSTPH__, startZ, endZ, __OBJHR__);
	if (dirX < dirY)
	{
		alphazmin = dirX;
		alphazmax = dirY;
	}
	else
	{
		alphazmin = dirY;
		alphazmax = dirX;
	}


	const float alphaMIN = max(alphaxmin, max(alphaymin, alphazmin));
	const float alphaMAX = min(alphaxmax, min(alphaymax, alphazmax));
	dev_minmaxIdxFun(startX, endX, __MINOL__, __OBJSTPL__, alphaMIN, alphaMAX, alphaxmin, alphaxmax, __OBJLR__ + 1, &imin, &imax);
	dev_minmaxIdxFun(startY, endY, __MINOW__, __OBJSTPW__, alphaMIN, alphaMAX, alphaymin, alphaymax, __OBJWR__ + 1, &jmin, &jmax);
	dev_minmaxIdxFun(startZ, endZ, __MINOH__, __OBJSTPH__, alphaMIN, alphaMAX, alphazmin, alphazmax, __OBJHR__ + 1, &kmin, &kmax);

	float alphaX = (startX < endX) ? dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, imin) : dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, imax);
	float alphaY = (startY < endY) ? dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, jmin) : dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, jmax);
	float alphaZ = (startZ < endZ) ? dev_alpha_IFun(__MINOH__, __OBJSTPH__, startZ, endZ, kmin) : dev_alpha_IFun(__MINOH__, __OBJSTPH__, startZ, endZ, kmax);

	int Np = static_cast<int>(abs(static_cast<float>(imax - imin)) + abs(static_cast<float>(jmax - jmin)) + abs(static_cast<float>(kmax - kmin)) + 3.0f);
	const float alphaxu = dev_alphaU_Fun(__OBJSTPL__, startX, endX);
	const float alphayu = dev_alphaU_Fun(__OBJSTPW__, startY, endY);
	const float alphazu = dev_alphaU_Fun(__OBJSTPH__, startZ, endZ);

	float alphaC = alphaMIN;
	//const float minApa = MY_MIN<float>(alphaX,alphaY,alphaZ);

	int i = int(dev_varphiFun(alphaMIN*1.00003f, __MINOL__, __OBJSTPL__, startX, endX));
	int j = int(dev_varphiFun(alphaMIN*1.00003f, __MINOW__, __OBJSTPW__, startY, endY));
	int k = int(dev_varphiFun(alphaMIN*1.00003f, __MINOH__, __OBJSTPH__, startZ, endZ));
	//int i = floor(dev_varphiFun((alphaMIN+minApa) * 0.5f, __MINOL__, __OBJSTPL__, startX, endX));
	//int j = floor(dev_varphiFun((alphaMIN+minApa) * 0.5f, __MINOW__, __OBJSTPW__, startY, endY));
	//int k = floor(dev_varphiFun((alphaMIN+minApa) * 0.5f, __MINOH__, __OBJSTPH__, startZ, endZ));

	const int iuu = (startX < endX) ? 1 : -1;
	const int juu = (startY < endY) ? 1 : -1;
	const int kuu = (startZ < endZ) ? 1 : -1;

	float d12(0);
	float weight(0);


	for (int repIdx(0); repIdx < Np; ++repIdx)
	{
		if (i < 0 || i > static_cast<int>(__OBJLR__ - 1) || j < 0 || j > static_cast<int>(__OBJWR__ - 1) || k < 0 || k > static_cast<int>(__OBJHR__ - 1))
		{
			break;
		}
		if (alphaX <= alphaY&&alphaX <= alphaZ)
		{
			weight = (alphaX - alphaC) * dconv;
			if (weight > 0)
			{
				d12 = d12 + weight * tex3D<float>(volTex, k,i,j);
				(*totWeight) += weight;
			}

			i = i + iuu;
			alphaC = alphaX;
			alphaX = alphaX + alphaxu;
			continue;
		}
		else if (alphaY <= alphaX&&alphaY <= alphaZ)
		{
			weight = (alphaY - alphaC) * dconv;
			if (weight > 0)
			{
				d12 = d12 + weight * tex3D<float>(volTex, k,i,j);
				(*totWeight) += weight;
			}
			j = j + juu;
			alphaC = alphaY;
			alphaY = alphaY + alphayu;

			continue;
		}
		else
		{
			weight = (alphaZ - alphaC) * dconv;
			if (weight > 0)
			{
				d12 = d12 + weight * tex3D<float>(volTex, k,i,j);
				(*totWeight) += weight;
			}
			k = k + kuu;
			alphaC = alphaZ;
			alphaZ = alphaZ + alphazu;
			continue;
		}
	}

	return d12;
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

