/*
 * utilities.cuh
 *
 *  Created on: Sep 23, 2020
 *      Author: Rui Liu
 */

#ifndef UTILITIES_CUH_
#define UTILITIES_CUH_

#include "device_launch_parameters.h"
#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <thread>
#include <vector>

#include <thrust/binary_search.h>
#include <thrust/copy.h>
#include <thrust/find.h>
#include <thrust/functional.h>
#include <thrust/host_vector.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/tuple.h>

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef __PI__
#define __PI__
#define PI			3.14159265358979323846
#define PI_2		1.57079632679489661923
#define PI_4		0.785398163397448309616
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
        
#define DEBUG 1

#if DEBUG
#define CUDA_CHECK_RETURN(value) {											\
	cudaError_t _m_cudaStat = value;										\
	if (_m_cudaStat != cudaSuccess) {										\
		fprintf(stderr, "Error %s at line %d in file %s\n",					\
				cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);		\
		exit(1);															\
	} }
 // Same function as CUDA_CHECK_RETURN
#else
#define CUDA_CHECK_RETURN(value) {value;}
#endif

typedef unsigned char byte;

INLINE __host__ __device__ const float2 operator/(const float2& a, float b) { return make_float2(a.x / b, a.y / b); }
INLINE __host__ __device__ const float3 operator+(const float3& a, const float3& b) { return make_float3(a.x + b.x, a.y + b.y, a.z + b.z); }
INLINE __host__ __device__ const float3 operator-(const float3& a, const float3& b) { return make_float3(a.x - b.x, a.y - b.y, a.z - b.z); }
INLINE __host__ __device__ const double3 operator-(const double3& a, const double3& b) { return make_double3(a.x - b.x, a.y - b.y, a.z - b.z); }
INLINE __host__ __device__ const float2 operator-(const float2& a, const float2& b) { return make_float2(a.x - b.x, a.y - b.y); }
INLINE __host__ __device__ const float3 operator*(const float3& a, const float3& b) { return make_float3(a.x * b.x, a.y * b.y, a.z * b.z); }
INLINE __host__ __device__ const float3 operator*(const float3& a, float b) { return make_float3(a.x * b, a.y * b, a.z * b); }
INLINE __host__ __device__ const float3 operator/(const float3& a, const float3& b) { return make_float3(a.x / b.x, a.y / b.y, a.z / b.z); }
INLINE __host__ __device__ const float3 operator/(const float3& a, float b) { return make_float3(a.x / b, a.y / b, a.z / b); }
INLINE __host__ __device__ const double3 operator/(const double3& a, double b) { return make_double3(a.x / b, a.y / b, a.z / b); }
INLINE __host__ __device__ const float3 operator-(const int3& a, const float3& b) {	return make_float3(a.x - b.x, a.y - b.y, a.z - b.z); }
INLINE __host__ __device__ float length(const float2& a) { return sqrtf(a.x * a.x + a.y * a.y); }
INLINE __host__ __device__ float length(const float3& a) { return sqrtf(a.x * a.x + a.y * a.y + a.z * a.z); }
INLINE __host__ __device__ double length(const double3& a) { return sqrt(a.x * a.x + a.y * a.y + a.z * a.z); }
INLINE __host__ __device__ const float2 normalize(const float2& a) { return a / length(a); }
INLINE __host__ __device__ const float3 normalize(const float3& a) { return a / length(a); }
INLINE __host__ __device__ const double3 normalize(const double3& a) { return a / length(a); }
INLINE __host__ __device__ float fminf(const float2& a) { return fminf(a.x, a.y); }
INLINE __host__ __device__ float fminf(const float3& a) { return fminf(a.x, fminf(a.y, a.z)); }
INLINE __host__ __device__ float fmaxf(const float2& a) { return fmaxf(a.x, a.y); }
INLINE __host__ __device__ float fmaxf(const float3& a) { return fmaxf(a.x, fmaxf(a.y, a.z)); }
INLINE __host__ __device__ const float3 fminf(const float3& a, const float3& b) { return make_float3(fminf(a.x, b.x), fminf(a.y, b.y), fminf(a.z, b.z)); }
INLINE __host__ __device__ const float3 fmaxf(const float3& a, const float3& b) { return make_float3(fmaxf(a.x, b.x), fmaxf(a.y, b.y), fmaxf(a.z, b.z)); }
INLINE __host__ __device__ const float2 fminf(const float2& a, const float2& b) { return make_float2(fminf(a.x, b.x), fminf(a.y, b.y)); }
INLINE __host__ __device__ const float2 fmaxf(const float2& a, const float2& b) { return make_float2(fmaxf(a.x, b.x), fmaxf(a.y, b.y)); }

template<typename T> INLINE __host__ __device__ T lerp(T v0, T v1, T t) { return fma(t, v1, fma(-t, v0, v0)); }
template<typename T> INLINE __host__ __device__ T bilierp(T v0, T v1, T v2, T v3, T t1, T t2) {
	T vv0 = fma(t1, v1, fma(-t1, v0, v0));
	T vv1 = fma(t1, v3, fma(-t1, v2, v2));
	return fma(t2, vv1, fma(-t2, vv0, vv0));
}

INLINE __device__ double bilerp(int2 v0, int2 v1, int2 v2, int2 v3, float t1, float t2) {
	const double v0_ = __hiloint2double(v0.y, v0.x);
	const double v1_ = __hiloint2double(v1.y, v1.x);
	const double v2_ = __hiloint2double(v2.y, v2.x);
	const double v3_ = __hiloint2double(v3.y, v3.x);
	const double vv0 = v0_ * (1.0 - t1) + v1_ * t1;
	const double vv1 = v2_ * (1.0 - t1) + v3_ * t1;
	return vv0 * (1 - t2) + vv1 * t2;
}

INLINE __device__ double bilerp(int2 v0, int2 v1, int2 v2, int2 v3, double t1, double t2) {
	const double v0_ = __hiloint2double(v0.y, v0.x);
	const double v1_ = __hiloint2double(v1.y, v1.x);
	const double v2_ = __hiloint2double(v2.y, v2.x);
	const double v3_ = __hiloint2double(v3.y, v3.x);

	const double vv0 = v0_ * (1.0 - t1) + v1_ * t1;
	const double vv1 = v2_ * (1.0 - t1) + v3_ * t1;
	return vv0 * (1 - t2) + vv1 * t2;
}

INLINE __host__ __device__ bool intersectBox(const float3& sour, const float3& dir,	const float3& boxmin, const float3& boxmax,	float* tnear, float* tfar) {
	const float3 invR = make_float3(1.0f / dir.x, 1.0f / dir.y, 1.0f / dir.z);
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

template<typename T>INLINE __host__ __device__ T regularizeAngle(T curang) {
	T c = curang;
	while (c >= TWOPI) { c -= TWOPI; }
	while (c < 0) { c += TWOPI; }
	return c;
}


INLINE __host__ __device__ void invRotVox(const float3& curVox, float3& virVox,	const float2& cossinT, const float zP) {
	virVox.x = curVox.x * cossinT.x + curVox.y * cossinT.y;
	virVox.y = -curVox.x * cossinT.y + curVox.y * cossinT.x;
	virVox.z = curVox.z - zP;
}

INLINE __device__ float3 invRot(const float3 inV, const float2 cossin, const float zP) {
	float3 outV;
	outV.x = inV.x * cossin.x + inV.y * cossin.y;
	outV.x = -inV.x * cossin.y + inV.y * cossin.x;
	outV.z = inV.z - zP;
	return outV;
}


namespace CTMBIR {
	struct ConstantForBackProjection3 {
		float x0;
		float y0;
		float z0;

		using InTuple = thrust::tuple<float, float>;
		using OutTuple = thrust::tuple<float3, float3, float2>;
		ConstantForBackProjection3(
			const float _x0, const float _y0, const float _z0) :x0(_x0),
			y0(_y0), z0(_z0) {}

		__device__ OutTuple operator()(const InTuple& tp) {
			const float curang = regularizeAngle(thrust::get<0>(tp));
			const float zP = thrust::get<1>(tp);
			const float cosT = cosf(curang);
			const float sinT = sinf(curang);
			float3 cursour = make_float3(x0 * cosT - y0 * sinT,	x0 * sinT + y0 * cosT, z0 + zP);
			float2 dirsour = normalize(make_float2(-cursour.x, -cursour.y));
			return thrust::make_tuple(make_float3(cosT, sinT, zP), cursour, dirsour);
		}
	};

	struct ConstantForBackProjection4 {
		float x0;
		float y0;
		float z0;

		using InTuple = thrust::tuple<float, float> ;
		ConstantForBackProjection4(const float _x0, const float _y0, const float _z0)
			: x0(_x0), y0(_y0), z0(_z0) {}

		__device__ float3 operator()(const InTuple& tp) {
			const float curang = regularizeAngle(thrust::get<0>(tp));
			const float zP = thrust::get<1>(tp);
			const float cosT = cosf(curang);
			const float sinT = sinf(curang);
			return make_float3(cosT, sinT, zP);
		}
	};

	struct Constant_MultiSlice {
		float x0;
		float y0;
		Constant_MultiSlice(const float _x0, const float _y0) :x0(_x0), y0(_y0) {}

		__device__ float2 operator()(const float& tp) {
			const float curang = regularizeAngle(tp);
			return make_float2(cosf(curang), sinf(curang));
		}
	};

	struct Constant_MultiSlice_host {
		float x0;
		float y0;
		Constant_MultiSlice_host(const float _x0, const float _y0) :x0(_x0), y0(_y0) {}

		float2 operator()(const float& tp) {
			const float curang = regularizeAngle(tp);
			return make_float2(cosf(curang), sinf(curang));
		}
	};

	template<typename T> struct ConstantForBackProjection {
		T x0;
		T y0;
		T z0;

		using InTuple = thrust::tuple<T, T>;
		ConstantForBackProjection(const T _x0, const T _y0, const T _z0)
			: x0(_x0), y0(_y0), z0(_z0) {}

		__device__ float3 operator()(const InTuple& tp) {
			const T curang = regularizeAngle(thrust::get<0>(tp));
			const T zP = thrust::get<1>(tp);
			const T cosT = cosf(curang);
			const T sinT = sinf(curang);
			return make_float3(cosT, sinT, zP);
		}
	};


	template<> struct ConstantForBackProjection<double> {
		double x0;
		double y0;
		double z0;

		using InTuple = thrust::tuple<double, double> ;
		ConstantForBackProjection(const double _x0, const double _y0, const double _z0)
			: x0(_x0), y0(_y0), z0(_z0) {}

		__device__ double3 operator()(const InTuple& tp) {
			const double curang = regularizeAngle(thrust::get<0>(tp));
			const double zP = thrust::get<1>(tp);
			const double cosT = cos(curang);
			const double sinT = sin(curang);
			return make_double3(cosT, sinT, zP);
		}
	};
}


template<typename T>
void DD3Boundaries(int nrBoundaries, T* pCenters, T* pBoundaries) {
	int i = 0;
	if (nrBoundaries >= 3) {
		for (i = 0; i < nrBoundaries; i++) {
			if (i == 0) {
				pBoundaries[0] = 1.5 * pCenters[0] - 0.5 * pCenters[1];
			}
			else if (i == nrBoundaries - 1) {
				pBoundaries[i] = 1.5 * pCenters[i - 1] - 0.5 * pCenters[i - 2];
			}
			else {
				pBoundaries[i] = 0.5 * pCenters[i - 1] + 0.5 * pCenters[i];
			}
		}
	} else {
		*pBoundaries = *pCenters - 0.5;
		*(pBoundaries + 1) = *pCenters + 0.5;
	}
}


///////////////////////////////////////////////////////////////////////////////////
// Get one sub-volume from the whole volume.
// Assume that the volumes are stored in Z, X, Y order
template<typename T>
void getSubVolume(const T* vol,
	const size_t XN, const size_t YN, const size_t ZN,
	const size_t ZIdx_Start, const size_t ZIdx_End, T* subVol) {
	const size_t SZN = ZIdx_End - ZIdx_Start;
	for (size_t yIdx = 0; yIdx != YN; ++yIdx) {
		for (size_t xIdx = 0; xIdx != XN; ++xIdx) {
			for (size_t zIdx = ZIdx_Start; zIdx != ZIdx_End; ++zIdx) {
				subVol[(yIdx * XN + xIdx) * SZN + (zIdx - ZIdx_Start)] = vol[(yIdx * XN + xIdx) * ZN + zIdx];
			}
		}
	}
}

template<typename T>
void getSubVolume(const T* vol,
	const size_t XYN, const size_t ZN,
	const size_t ZIdx_Start, const size_t ZIdx_End, T* subVol) {
	const size_t SZN = ZIdx_End - ZIdx_Start;
	for (size_t xyIdx = 0; xyIdx != XYN; ++xyIdx) {
		for (size_t zIdx = ZIdx_Start; zIdx != ZIdx_End; ++zIdx) {
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
	const T objCntIdxZ, const T dz, const int ZN, // Size of the volume
	int& ObjIdx_Start, int& ObjIdx_End) { // The end is not included

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
	int& prjIdx_Start, int& prjIdx_End) {
	const int PN = zPos.size();

	const T lowerPartV = (ObjZIdx_Start - objCntIdxZ - 0.5) * dz;
	const T highrPartV = lowerPartV + (ObjZIdx_End - ObjZIdx_Start) * dz;

	const T lowerPartDet = (detCntIdxV + 0.5) * detStpZ;
	const T upperPartDet = DNV * detStpZ - lowerPartDet;

	//The source position
	const T sourLPos = lowerPartV - upperPartDet;
	const T sourHPos = highrPartV + lowerPartDet;

	prjIdx_Start = thrust::upper_bound(zPos.begin(), zPos.end(), sourLPos) - zPos.begin() - 1;
	prjIdx_End = thrust::upper_bound(zPos.begin(), zPos.end(), sourHPos) - zPos.begin() + 2;
	prjIdx_Start = (prjIdx_Start < 0) ? 0 : prjIdx_Start;
	prjIdx_Start = (prjIdx_Start > PN) ? PN : prjIdx_Start;

	prjIdx_End = (prjIdx_End < 0) ? 0 : prjIdx_End;
	prjIdx_End = (prjIdx_End > PN) ? PN : prjIdx_End;
}


////////////////////////////////////////////////////////////////////////////////////
// The volume is also stored in Z, X, Y order
// Not tested yet.
template<typename T> void combineVolume(
	T* vol, // The volume to be combined
	const int XN, const int YN, const int ZN,
	T** subVol, // All sub volumes
	const int* SZN, // Number of slices for each subVolume
	const int subVolNum) { // Number of sub volumes
	int kk = 0;
	for (size_t yIdx = 0; yIdx != YN; ++yIdx) {
		for (size_t xIdx = 0; xIdx != XN; ++xIdx) {
			kk = 0;
			for (size_t volIdx = 0; volIdx != subVolNum; ++volIdx) {
				for (size_t zIdx = 0; zIdx != SZN[volIdx]; ++zIdx) {
					vol[(yIdx * XN + xIdx) * ZN + kk] = subVol[volIdx][(yIdx * XN + xIdx) * SZN[volIdx] + zIdx];
					kk = kk + 1;
				}
			}
		}
	}
}

#endif /* UTILITIES_CUH_ */