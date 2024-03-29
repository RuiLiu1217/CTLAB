/*
* COPYRIGHT NOTICE
* All rights reserved
*
* \file projUtilities.h
* \brief The utilities functions used in GPU.
*
* \version 1.0
* \author Rui Liu
* \date Jun. 5, 2017
*
*/
#pragma once
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <omp.h>
#include <string>
#include <vector>
#include <cuda_runtime.h>
#include <cuda.h>

#include "device_launch_parameters.h"

#include <thrust/binary_search.h>
#include <thrust/copy.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>


// TODO: TO A STATIC VARIABLES
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

#ifndef nullptr
#define nullptr NULL
#endif


/// Key words for inlining the codes
#define FORCEINLINE 1
#if FORCEINLINE
#define INLINE __forceinline__
#else
#define INLINE inline
#endif

// Usually, for safety, we will use checkCudaErrors or CUDA_CHECK_RETURN before "cuda" API functions
#if DEBUG
#define CUDA_CHECK_RETURN(value) { cudaError_t _m_cudaStat = value; if(_m_cudaStat != cudaSuccess){fprintf(stderr, "Error %s at line %d in file %s\n", cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__); exit(1);}}
#define checkCudaErrors(value) { cudaError_t _m_cudaStat = value; if(_m_cudaStat != cudaSuccess){fprintf(stderr, "Error %s at line %d in file %s\n", cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__); exit(1);}}
#else
#define CUDA_CHECK_RETURN(value) {value;}
#define CUDA_SAFE_CALL(value) {value;}
#define checkCudaErrors(value) {value;}
#endif

using byte = unsigned char;

/// \brief Define that macro for Thrust compiling successful
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS 1
#endif


/// \brief The function judging whether a variable is ZERO according to the EPSILON value
/// \param x input parameter for judging whether it is zero or not
template<typename T>
inline __host__ __device__ bool IS_ZERO(const T& x)
{
	return (abs(x) < 1.0E-9);
}

/// \brief SIDDON kernel function 1
inline __host__ __device__ float dev_pFun(const float& alpha, const float& pstart, const float&pend)
{
	return pstart + alpha * (pend - pstart);
}

/// \brief SIDDON kernel function 2
inline __host__ __device__ float dev_alpha_IFun(const float& b, const float& d, const float& pstart, const float& pend, const unsigned int& i)
{
	if (!IS_ZERO<float>(pend - pstart))
	{
		return ((b + (float)i*d) - pstart) / (pend - pstart);
	}
	else return 1000;//((b + i*d)-pstart)/(1e-6);
}



/// \brief SIDDON kernel function 3
inline __host__ __device__ float dev_varphiFun(const float& alpha, const float& b, const float& d, const float& pstart, const float& pend)
{
	return (dev_pFun(alpha, pstart, pend) - b) / d;
}


/// \brief SIDDON kernel function 4
inline __host__ __device__ void dev_minmaxIdxFun(
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
inline __host__ __device__  float dev_alphaU_Fun(const float& d, const float& startx, const float& endx)
{
	if (IS_ZERO<float>(startx - endx))
	{
		return 1000.0f;//(d/1e-6);
	}
	return d / fabsf(startx - endx);
}


/// \brief SIDDON kernel function 6
template<typename T>
inline __host__ __device__  int dev_iu_Fun(const T& start, const T& end)
{
	return (start < end) ? 1 : -1;
}



template<typename T>
INLINE __host__ __device__ T safeDivide(const T& a, const T& b)
{
	if (IS_ZERO(b))
	{
		return 0;
	}
	else
	{
		return a / b;
	}
}


//INLINE __host__ __device__ float2 operator*(const float2& a, const float2& b) { return make_float2(a.x * b.x, a.y * b.y); };
INLINE __host__ __device__ float2 operator*(const float a, const float2& b) { return make_float2(a * b.x, a * b.y); };
INLINE __host__ __device__ float2 operator*(const float2& a, const float b) { return make_float2(a.x * b, a.y * b); };
INLINE __host__ __device__ float2 operator/(const float2& a, const float2& b) { return make_float2(safeDivide<float>(a.x, b.x), safeDivide<float>(a.y, b.y)); };
INLINE __host__ __device__ float2 operator/(const float2& a, const float b) { return make_float2(safeDivide<float>(a.x, b), safeDivide<float>(a.y, b)); };

INLINE __host__ __device__ float3 operator+(const float3& a, const float3& b) { return make_float3(a.x + b.x, a.y + b.y, a.z + b.z); };
INLINE __host__ __device__ float3 operator+(const float& a, const float3& b) { return make_float3(a + b.x, a + b.y, a + b.z); };
INLINE __host__ __device__ float3 operator+=(const float3& a, const float3& b) { return make_float3(a.x + b.x, a.y + b.y, a.z + b.z); };
INLINE __host__ __device__ float3 operator-(const float3& a, const float3& b) { return make_float3(a.x - b.x, a.y - b.y, a.z - b.z); }
INLINE __host__ __device__ float3 operator-(const int3& a, const float3& b) { return make_float3(a.x - b.x, a.y - b.y, a.z - b.z); }
INLINE __host__ __device__ float3 operator*(const float3& a, const float& b) { return make_float3(a.x * b, a.y * b, a.z * b); };
INLINE __host__ __device__ float3 operator*(const float& a, const float3& b) { return make_float3(a * b.x, a * b.y, a * b.z); };
INLINE __host__ __device__ float3 operator*(const float3& a, const float3& b) { return make_float3(a.x * b.x, a.y * b.y, a.z * b.z); }
INLINE __host__ __device__ float3 operator/(const float3& a, const float& b) { return make_float3(safeDivide<float>(a.x, b), safeDivide<float>(a.y, b), safeDivide<float>(a.z, b)); };
INLINE __host__ __device__ float3 operator/(const float3& a, const float3& b) { return make_float3(safeDivide<float>(a.x, b.x), safeDivide<float>(a.y, b.y), safeDivide<float>(a.z, b.z)); }

INLINE __host__ __device__ double2 operator+(const double2& a, const double2& b) { return make_double2(a.x + b.x, a.y + b.y); }
INLINE __host__ __device__ double2 operator+(const double2& a, const float2& b) { return make_double2(a.x + b.x, a.y + b.y); }
INLINE __host__ __device__ double2 operator+(const float2& a, const double2& b) { return make_double2(a.x + b.x, a.y + b.y); }
INLINE __host__ __device__ double2 operator-(const double2& a, const double2& b) { return make_double2(a.x - b.x, a.y - b.y); }
INLINE __host__ __device__ double2 operator-(const double2& a, const float2& b) { return make_double2(a.x - b.x, a.y - b.y); }
INLINE __host__ __device__ double2 operator-(const float2& a, const double2& b) { return make_double2(a.x - b.x, a.y - b.y); }
INLINE __host__ __device__ double2 operator*(const double2& a, const double2& b) { return make_double2(a.x * b.x, a.y * b.y); }
INLINE __host__ __device__ double2 operator*(const double2& a, const float2& b) { return make_double2(a.x * b.x, a.y * b.y); }
INLINE __host__ __device__ double2 operator*(const float2& a, const double2& b) { return make_double2(a.x * b.x, a.y * b.y); }
INLINE __host__ __device__ double2 operator/(const double2& a, const double2& b) { return make_double2(safeDivide<double>(a.x, b.x), safeDivide<double>(a.y, b.y)); }
INLINE __host__ __device__ double2 operator/(const double2& a, const float2& b) { return make_double2(safeDivide<double>(a.x, b.x), safeDivide<double>(a.y, b.y)); }
INLINE __host__ __device__ double2 operator/(const float2& a, const double2& b) { return make_double2(safeDivide<double>(a.x, b.x), safeDivide<double>(a.y, b.y)); }
INLINE __host__ __device__ double2 operator/(const double2& a, double b) { return make_double2(a.x / b, a.y / b); }

INLINE __host__ __device__ double3 operator+(const double3& a, const double3& b) { return make_double3(a.x + b.x, a.y + b.y, a.z + b.z); }
INLINE __host__ __device__ double3 operator+(const double3& a, const float3& b) { return make_double3(a.x + b.x, a.y + b.y, a.z + b.z); }
INLINE __host__ __device__ double3 operator+(const float3& a, const double3& b) { return make_double3(a.x + b.x, a.y + b.y, a.z + b.z); }
INLINE __host__ __device__ double3 operator-(const double3& a, const double3& b) { return make_double3(a.x - b.x, a.y - b.y, a.z - b.z); }
INLINE __host__ __device__ double3 operator-(const int3& a, const double3& b) { return make_double3(a.x - b.x, a.y - b.y, a.z - b.z); }
INLINE __host__ __device__ double3 operator-(const double3& a, const float3& b) { return make_double3(a.x - b.x, a.y - b.y, a.z - b.z); }
INLINE __host__ __device__ double3 operator-(const float3& a, const double3& b) { return make_double3(a.x - b.x, a.y - b.y, a.z - b.z); }
INLINE __host__ __device__ double3 operator*(const double3& a, const double3& b) { return make_double3(a.x * b.x, a.y * b.y, a.z * b.z); }
INLINE __host__ __device__ double3 operator*(const double3& a, const float3& b) { return make_double3(a.x * b.x, a.y * b.y, a.z * b.z); }
INLINE __host__ __device__ double3 operator*(const float3& a, const double3& b) { return make_double3(a.x * b.x, a.y * b.y, a.z * b.z); }
INLINE __host__ __device__ double3 operator/(const double3& a, const double3& b) { return make_double3(safeDivide<double>(a.x, b.x), safeDivide<double>(a.y, b.y), safeDivide<double>(a.z, b.z)); }
INLINE __host__ __device__ double3 operator/(const double3& a, const float3& b) { return make_double3(safeDivide<double>(a.x, b.x), safeDivide<double>(a.y, b.y), safeDivide<double>(a.z, b.z)); }
INLINE __host__ __device__ double3 operator/(const float3& a, const double3& b) { return make_double3(safeDivide<double>(a.x, b.x), safeDivide<double>(a.y, b.y), safeDivide<double>(a.z, b.z)); }
INLINE __host__ __device__ double3 operator/(const double3& a, const double& b) { return make_double3(safeDivide<double>(a.x, b), safeDivide<double>(a.y, b), safeDivide<double>(a.z, b)); };

INLINE __host__ __device__ float fminf(const float2& a) { return fminf(a.x, a.y); }
INLINE __host__ __device__ float fminf(const float3& a) { return fminf(a.x, fminf(a.y, a.z)); }
INLINE __host__ __device__ float2 fminf(const float2& a, const float2& b) { return make_float2(fmin(a.x, b.x), fmin(a.y, b.y)); }
INLINE __host__ __device__ float3 fminf(const float3& a, const float3& b) { return make_float3(fminf(a.x, b.x), fminf(a.y, b.y), fminf(a.z, b.z)); }
INLINE __host__ __device__ double2 fminf(const double2& a, const double2& b) { return make_double2(fmin(a.x, b.x), fmin(a.y, b.y)); }
INLINE __host__ __device__ double2 fminf(const float2& a, const double2& b) { return make_double2(fminf(a.x, b.x), fminf(a.y, b.y)); }
INLINE __host__ __device__ double2 fminf(const double2& a, const float2& b) { return make_double2(fminf(a.x, b.x), fminf(a.y, b.y)); }
INLINE __host__ __device__ double3 fminf(const double3& a, const double3& b) { return make_double3(fminf(a.x, b.x), fminf(a.y, b.y), fminf(a.z, b.z)); }
INLINE __host__ __device__ double3 fminf(const double3& a, const float3& b) { return make_double3(fminf(a.x, b.x), fminf(a.y, b.y), fminf(a.z, b.z)); }
INLINE __host__ __device__ double3 fminf(const float3& a, const double3& b) { return make_double3(fminf(a.x, b.x), fminf(a.y, b.y), fminf(a.z, b.z)); }
INLINE __host__ __device__ float fmaxf(const float2& a) { return fmaxf(a.x, a.y); }
INLINE __host__ __device__ float fmaxf(const float3& a) { return fmaxf(a.x, fmaxf(a.y, a.z)); }
INLINE __host__ __device__ float2 fmaxf(const float2& a, const float2& b) { return make_float2(fmax(a.x, b.x), fmax(a.y, b.y)); }
INLINE __host__ __device__ float3 fmaxf(const float3& a, const float3& b) { return make_float3(fmaxf(a.x, b.x), fmaxf(a.y, b.y), fmaxf(a.z, b.z)); }
INLINE __host__ __device__ double2 fmaxf(const double2& a, const double2& b) { return make_double2(fmaxf(a.x, b.x), fmaxf(a.y, b.y)); }
INLINE __host__ __device__ double2 fmaxf(const float2& a, const double2& b) { return make_double2(fmaxf(a.x, b.x), fmaxf(a.y, b.y)); }
INLINE __host__ __device__ double2 fmaxf(const double2& a, const float2& b) { return make_double2(fmaxf(a.x, b.x), fmaxf(a.y, b.y)); }
INLINE __host__ __device__ double3 fmaxf(const double3& a, const double3& b) { return make_double3(fmaxf(a.x, b.x), fmaxf(a.y, b.y), fmaxf(a.z, b.z)); }
INLINE __host__ __device__ double3 fmaxf(const double3& a, const float3& b) { return make_double3(fmaxf(a.x, b.x), fmaxf(a.y, b.y), fmaxf(a.z, b.z)); }
INLINE __host__ __device__ double3 fmaxf(const float3& a, const double3& b) { return make_double3(fmaxf(a.x, b.x), fmaxf(a.y, b.y), fmaxf(a.z, b.z)); }

INLINE __host__ __device__ float length(const float2& a) { return sqrtf(a.x * a.x + a.y * a.y); }
INLINE __host__ __device__ float length(const float3& a) { return sqrtf(a.x * a.x + a.y * a.y + a.z * a.z); }
INLINE __host__ __device__ double length(const double2& a) { return sqrt(a.x * a.x + a.y * a.y); }
INLINE __host__ __device__ double length(const double3& a) { return sqrt(a.x * a.x + a.y * a.y + a.z * a.z); }

//INLINE __host__ __device__ const float2 normalize(const float2& a){	return a / length(a);}
INLINE __host__ __device__ const float2 normalize(float2 a) { return a / length(a); }
INLINE __host__ __device__ const float3 normalize(const float3& a) { return a / length(a); }
INLINE __host__ __device__ const double2 normalize(const double2& a) { return a / length(a); }
INLINE __host__ __device__ const double3 normalize(const double3& a) { return a / length(a); }



template<typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b)
{
	std::vector<T> c(a.size(), 0);
	std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::minus<T>());
	return c;
}


/// \brief Fast linear interpolation
template<typename T>
INLINE __host__ __device__ T lerp(T v0, T v1, T t)
{
	return fma(t, v1, fma(-t, v0, v0));
}

/// \brief Fast bilinear interpolation
template<typename T>
INLINE __host__ __device__ T bilierp(T v0, T v1, T v2, T v3, T t1, T t2)
{
	T vv0 = fma(t1, v1, fma(-t1, v0, v0));
	T vv1 = fma(t1, v3, fma(-t1, v2, v2));
	return fma(t2, vv1, fma(-t2, vv0, vv0));
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
	while (c >= TWOPI) { c -= TWOPI; }
	while (c < 0) { c += TWOPI; }
	return c;
}



INLINE __device__ float3 invRot(
	const float3 inV,
	const float2 cossin,
	const float zP)
{
	float3 outV;
	outV.x = inV.x * cossin.x + inV.y * cossin.y;
	outV.x = -inV.x * cossin.y + inV.y * cossin.x;
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
			y0(_y0), z0(_z0) {}

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
	struct ConstantForBackProjection4 {

		float x0;
		float y0;
		float z0;

		typedef thrust::tuple<float, float> InTuple;
		ConstantForBackProjection4(const float _x0, const float _y0, const float _z0)
			: x0(_x0), y0(_y0), z0(_z0) {}

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
	struct ConstantForBackProjection {

		T x0;
		T y0;
		T z0;

		typedef thrust::tuple<T, T> InTuple;
		ConstantForBackProjection(const T _x0, const T _y0, const T _z0)
			: x0(_x0), y0(_y0), z0(_z0) {}

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
	struct ConstantForBackProjection<double> {

		double x0;
		double y0;
		double z0;

		typedef thrust::tuple<double, double> InTuple;
		ConstantForBackProjection(const double _x0, const double _y0, const double _z0)
			: x0(_x0), y0(_y0), z0(_z0) {}

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



/// \brief SIDDON kernel function 1
inline __host__	__device__ double dev_pFun(const double& alpha, const double& pstart, const double&pend)
{
	return pstart + alpha * (pend - pstart);
}


/// \brief SIDDON kernel function 2
inline __host__	__device__ double dev_alpha_IFun(const double& b, const double& d, const double& pstart, const double& pend, const unsigned int& i)
{
	if (!IS_ZERO<double>(pend - pstart))
	{
		return ((b + (double)i*d) - pstart) / (pend - pstart);
	}
	else return 1000;//((b + i*d)-pstart)/(1e-6);
}

/// \brief SIDDON kernel function 3
inline __host__	__device__ double dev_varphiFun(const double& alpha, const double& b, const double& d, const double& pstart, const double& pend)
{
	return (dev_pFun(alpha, pstart, pend) - b) / d;
}

/// \brief SIDDON kernel function 4
inline	__host__ __device__ void dev_minmaxIdxFun(
	const double& pstart, const double& pend,
	const double& b, const double& d,
	const double& alphaMIN, const double& alphaMAX,
	const double& alphaPmin, const double& alphaPmax,
	const unsigned int& Nplane, int* imin, int* imax)
{
	if (pstart < pend)
	{
		if (IS_ZERO<double>(alphaMIN - alphaPmin))
		{
			*imin = 1;
		}
		else
		{
			*imin = static_cast<int>(ceil(dev_varphiFun(alphaMIN, b, d, pstart, pend)));
		}
		if (IS_ZERO<double>(alphaMAX - alphaPmax))
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
		if (IS_ZERO<double>(alphaMIN - alphaPmin))
		{
			*imax = Nplane - 2;
		}
		else
		{
			*imax = static_cast<int>(floor(dev_varphiFun(alphaMIN, b, d, pstart, pend)));
		}
		if (IS_ZERO<double>(alphaMAX - alphaPmax))
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
inline __host__ __device__  double dev_alphaU_Fun(const double& d, const double& startx, const double& endx)
{
	if (IS_ZERO<double>(startx - endx))
	{
		return 1000.0f;//(d/1e-6);
	}
	return d / fabs(startx - endx);
}


template<typename T>
__host__ __device__ inline T intersectLength(const T& fixedmin, const T& fixedmax, const T& varimin, const T& varimax)
{
	const T left = (fixedmin > varimin) ? fixedmin : varimin;
	const T right = (fixedmax < varimax) ? fixedmax : varimax;
	return abs(right - left) * static_cast<double>(right > left);
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


template<typename T>
void DD3Boundaries(int nrBoundaries, std::vector<T>& Centers, std::vector<T>& Boundaries)
{
	DD3Boundaries<T>(nrBoundaries, &Centers[0], &Boundaries[0]);
}

template<typename T>
void DD3Boundaries(int nrBoundaries, thrust::host_vector<T>& Centers, thrust::host_vector<T>& Boundaries)
{
	DD3Boundaries<T>(nrBoundaries, &Centers[0], &Boundaries[0]);
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
	const T objCntIdxZ, const T dz, const int ZN, // Size of the volume
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

