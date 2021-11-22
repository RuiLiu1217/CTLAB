/*
 * COPYRIGHT NOTICE
 * COPYRIGHT (c) 2015, GEGRC, Wake Forest and UMass Lowell
 * All rights reserved
 *
 * \file utilities.hpp
 * \brief The utilities functions used in GPU.
 *
 * \version 1.0
 * \author Rui Liu
 * \date Nov. 21, 2021
 *
 */
#pragma once
#include <algorithm>
#include <cassert>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <omp.h>
#include <random>
#include <sstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <vector_functions.h>
#include <vector_types.h>

#include <cusparse.h>
#include <cusparse_v2.h>
#include "cublas_v2.h"
#include <cuda_runtime.h>
#include <cuda.h>
#include <device_functions.h>
#include "device_launch_parameters.h"

#include <thrust/binary_search.h>
#include <thrust/copy.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <thrust/fill.h>
#include <thrust/find.h>
#include <thrust/functional.h>
#include <thrust/host_vector.h>
#include <thrust/inner_product.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/reduce.h>
#include <thrust/remove.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/transform.h>
#include <thrust/transform_reduce.h>
#include <thrust/tuple.h>
#include <thrust/device_vector.h>


#include "FanEAGeo.h"
#include "FanEDGeo.h"
#include "ConeEAGeo.h"
#include "ConeEDGeo.h"
#include "Image.h"
#include "Volume.h"

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

typedef unsigned char byte;
typedef unsigned int uint;
typedef const unsigned int cuint;
typedef unsigned char byte;

typedef thrust::device_vector<float> d_vec_t;
typedef thrust::host_vector<float> h_vec_t;


/// \brief Define that macro for Thrust compiling successful
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS 1
#endif

#ifndef _EPSILON
#define _EPSILON (1.0E-9)
#endif

/// \brief Ray class on 2D
struct Ray2D
{
public:
	float2 o;///< origin of the light
	float2 d; ///< direction of the light
};

/// \brief Ray class on 3D
struct Ray
{
public:
	float3 o;///< origin of the light
	float3 d; ///< direction of the light
};

/// \brief Grid describing four corners of the image
struct Grid
{
public:
	float3 P[4];///< Four corners of the image
};



/// \brief The function judging whether a variable is ZERO according to the EPSILON value
/// \param x input parameter for judging whether it is zero or not
template<typename T>
inline __host__ __device__ bool IS_ZERO(const T& x) {
	return (abs(x) < 1.0E-9);
}

/// \brief SIDDON kernel function 1
template<typename T>
inline __host__ __device__ T dev_pFun(const T& alpha, const T& pstart, const T&pend) {
	return pstart + alpha * (pend - pstart);
}

/// \brief SIDDON kernel function 2
template<typename T>
inline __host__ __device__ T dev_alpha_IFun(const T& b, const T& d, const T& pstart, const T& pend, const unsigned int& i) {
	if (!IS_ZERO<T>(pend - pstart)) {
		return ((b + i * d) - pstart) / (pend - pstart);
	}
	else return 1000;//((b + i*d)-pstart)/(1e-6);
}

/// \brief SIDDON kernel function 3
template<typename T>
inline __host__	__device__ T dev_varphiFun(const T& alpha, const T& b, const T& d, const T& pstart, const T& pend) {
	return (dev_pFun(alpha, pstart, pend) - b) / d;
}

/// \brief SIDDON kernel function 5
template<typename T>
inline __host__ __device__  T dev_alphaU_Fun(const T& d, const T& startx, const T& endx) {
	if (IS_ZERO<T>(startx - endx))	{
		return 1000.0f;//(d/1e-6);
	}
	return d / fabs(startx - endx);
}


/// \brief SIDDON kernel function 6
template<typename T>
inline __host__ __device__  int dev_iu_Fun(const T& start, const T& end) {
	return (start < end) ? 1 : -1;
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



template<typename T>
INLINE __host__ __device__ T safeDivide(const T& a, const T& b) {
	if(IS_ZERO(b)) {
		return 0;
	} else {
		return a/b;
	}
}

INLINE __host__ __device__ float2 operator+(const float2& a, const float2& b){return make_float2(a.x + b.x, a.y + b.y);};
//INLINE __host__ __device__ float2 operator-(const float2& a, const float2& b){return make_float2(a.x - b.x, a.y - b.y);};
INLINE __host__ __device__ float2 operator-(float2 a, float2 b) { return make_float2(a.x - b.x, a.y - b.y); };
INLINE __host__ __device__ float2 operator*(const float2& a, const float2& b){return make_float2(a.x * b.x, a.y * b.y);};
INLINE __host__ __device__ float2 operator*(const float a, const float2& b){return make_float2(a * b.x, a * b.y);};
INLINE __host__ __device__ float2 operator*(const float2& a, const float b){return make_float2(a.x * b, a.y * b);};
INLINE __host__ __device__ float2 operator/(const float2& a, const float2& b){return make_float2(safeDivide<float>(a.x,b.x),safeDivide<float>(a.y,b.y));};
INLINE __host__ __device__ float2 operator/(const float2& a, const float b){return make_float2(safeDivide<float>(a.x,b), safeDivide<float>(a.y,b));};

INLINE __host__ __device__ float3 operator+(const float3& a, const float3& b){return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);};
INLINE __host__ __device__ float3 operator+(const float& a, const float3& b){return make_float3(a + b.x, a + b.y, a + b.z);};
INLINE __host__ __device__ float3 operator+=(const float3& a, const float3& b){	return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);};
INLINE __host__ __device__ float3 operator-(const float3& a, const float3& b){	return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);}
INLINE __host__ __device__ float3 operator-(const int3& a, const float3& b){	return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);}
INLINE __host__ __device__ float3 operator*(const float3& a, const float& b){return make_float3(a.x * b, a.y * b, a.z * b);};
INLINE __host__ __device__ float3 operator*(const float& a, const float3& b){return make_float3(a * b.x, a * b.y, a * b.z);};
INLINE __host__ __device__ float3 operator*(const float3& a, const float3& b){	return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);}
INLINE __host__ __device__ float3 operator/(const float3& a, const float& b){return make_float3(safeDivide<float>(a.x,b),safeDivide<float>(a.y,b),safeDivide<float>(a.z,b));};
INLINE __host__ __device__ float3 operator/(const float3& a, const float3& b){	return make_float3(safeDivide<float>(a.x, b.x),safeDivide<float>(a.y, b.y),safeDivide<float>(a.z, b.z));}

INLINE __host__ __device__ double2 operator+(const double2& a, const double2& b){ return make_double2(a.x + b.x, a.y + b.y); }
INLINE __host__ __device__ double2 operator+(const double2& a, const float2& b){ return make_double2(a.x + b.x, a.y + b.y); }
INLINE __host__ __device__ double2 operator+(const float2& a, const double2& b){ return make_double2(a.x + b.x, a.y + b.y); }
INLINE __host__ __device__ double2 operator-(const double2& a, const double2& b){ return make_double2(a.x - b.x, a.y - b.y); }
INLINE __host__ __device__ double2 operator-(const double2& a, const float2& b){ return make_double2(a.x - b.x, a.y - b.y); }
INLINE __host__ __device__ double2 operator-(const float2& a, const double2& b){ return make_double2(a.x - b.x, a.y - b.y); }
INLINE __host__ __device__ double2 operator*(const double2& a, const double2& b){ return make_double2(a.x * b.x, a.y * b.y); }
INLINE __host__ __device__ double2 operator*(const double2& a, const float2& b){ return make_double2(a.x * b.x, a.y * b.y); }
INLINE __host__ __device__ double2 operator*(const float2& a, const double2& b){ return make_double2(a.x * b.x, a.y * b.y); }
INLINE __host__ __device__ double2 operator/(const double2& a, const double2& b){ return make_double2(safeDivide<double>(a.x,b.x),safeDivide<double>(a.y,b.y));}
INLINE __host__ __device__ double2 operator/(const double2& a, const float2& b){return make_double2(safeDivide<double>(a.x,b.x),safeDivide<double>(a.y,b.y));}
INLINE __host__ __device__ double2 operator/(const float2& a, const double2& b){return make_double2(safeDivide<double>(a.x,b.x),safeDivide<double>(a.y,b.y));}
INLINE __host__ __device__ double2 operator/(const double2& a, double b) { return make_double2(a.x / b, a.y / b); }

INLINE __host__ __device__ double3 operator+(const double3& a, const double3& b){ return make_double3(a.x + b.x, a.y + b.y, a.z + b.z); }
INLINE __host__ __device__ double3 operator+(const double3& a, const float3& b){ return make_double3(a.x + b.x, a.y + b.y, a.z + b.z); }
INLINE __host__ __device__ double3 operator+(const float3& a, const double3& b){ return make_double3(a.x + b.x, a.y + b.y, a.z + b.z); }
INLINE __host__ __device__ double3 operator-(const double3& a, const double3& b){	return make_double3(a.x - b.x, a.y - b.y, a.z - b.z);}
INLINE __host__ __device__ double3 operator-(const int3& a, const double3& b) { return make_double3(a.x - b.x, a.y - b.y, a.z - b.z); }
INLINE __host__ __device__ double3 operator-(const double3& a, const float3& b){ return make_double3(a.x - b.x, a.y - b.y, a.z - b.z); }
INLINE __host__ __device__ double3 operator-(const float3& a, const double3& b){ return make_double3(a.x - b.x, a.y - b.y, a.z - b.z); }
INLINE __host__ __device__ double3 operator*(const double3& a, const double3& b){ return make_double3(a.x * b.x, a.y * b.y, a.z * b.z); }
INLINE __host__ __device__ double3 operator*(const double3& a, const float3& b){ return make_double3(a.x * b.x, a.y * b.y, a.z * b.z); }
INLINE __host__ __device__ double3 operator*(const float3& a, const double3& b){ return make_double3(a.x * b.x, a.y * b.y, a.z * b.z); }
INLINE __host__ __device__ double3 operator/(const double3& a, const double3& b){return make_double3(safeDivide<double>(a.x,b.x),safeDivide<double>(a.y,b.y),safeDivide<double>(a.z,b.z));}
INLINE __host__ __device__ double3 operator/(const double3& a, const float3& b){return make_double3(safeDivide<double>(a.x,b.x),safeDivide<double>(a.y,b.y),safeDivide<double>(a.z,b.z));}
INLINE __host__ __device__ double3 operator/(const float3& a, const double3& b){return make_double3(safeDivide<double>(a.x,b.x),safeDivide<double>(a.y,b.y),safeDivide<double>(a.z,b.z));}
INLINE __host__ __device__ double3 operator/(const double3& a, const double& b){return make_double3(safeDivide<double>(a.x,b),safeDivide<double>(a.y,b),safeDivide<double>(a.z,b));};


/// \brief pointwise multiplication
/// \param a input vector a
/// \param b input vector b
/// \return a.*b
template<typename T>
std::vector<T> operator*(const std::vector<T>& a, const std::vector<T>& b)
{
	std::vector<T> c(a.size(), 0);
	std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::multiplies<T>());
	return c;
}



template<typename T>
thrust::host_vector<T> operator-(
	const thrust::host_vector<T>& a,
	const thrust::host_vector<T>& b)
{
	thrust::host_vector<T> res(a);
	thrust::transform(res.begin(), res.end(), b.begin(), res.begin(), [=](T aa, T bb) {return aa - bb; });
	return res;
}



template<typename T>
thrust::host_vector<T> operator+(thrust::host_vector<T>& a, thrust::host_vector<T>& b)
{
	thrust::host_vector<T> c(a.size(), 0);
	thrust::transform(a.begin(), a.end(), b.begin(), c.begin(), thrust::plus<T>());
	return c;
}



template<typename T>
thrust::device_vector<T> operator+(thrust::device_vector<T>& a, thrust::device_vector<T>& b)
{
	thrust::device_vector<T> c(a.size(), 0);
	thrust::transform(a.begin(), a.end(), b.begin(), c.begin(), thrust::plus<T>());
	return c;
}



template<typename T>
thrust::device_vector<T> operator-(thrust::device_vector<T>& a, thrust::device_vector<T>& b)
{
	thrust::device_vector<T> c(a.size(), 0);
	thrust::transform(a.begin(), a.end(), b.begin(), c.begin(), thrust::minus<T>());
	return c;
}


/// \brief pointwise multiplication
/// \param a input vector a
/// \param b input vector b
/// \return a.*b
template<typename T>
thrust::device_vector<T> operator*(thrust::device_vector<T>& a, thrust::device_vector<T>& b)
{
	thrust::device_vector<T> c(a.size(), 0);
	thrust::transform(a.begin(), a.end(), b.begin(), c.begin(), thrust::multiplies<T>());
	return c;
}



/// \brief scale-vector multiplication
/// \param a input vector
/// \param s scale in R
/// \return the scaled vector of a
template<typename T>
thrust::device_vector<T> operator*(thrust::device_vector<T>& a, const T& s)
{
	thrust::device_vector<T> res(a.size(), 0);
	thrust::transform(a.begin(), a.end(), res.begin(), [&](T v) {return s * v; });
	return res;
}

/// \brief operator overload for multiplication with device vector
/// \param a input vector
/// \param s scale in R
/// \return the scaled vector of a
template<typename T>
thrust::device_vector<T> operator*(const T& s, thrust::device_vector<T>& a)
{
	thrust::device_vector<T> res(a.size(), 0);
	thrust::transform(a.begin(), a.end(), res.begin(), [&](T v) {return s * v; });
	return res;
}



/// \brief operator overload for inner product of two vectors
/// \param a one device vector
/// \param b one device vector
/// \return the inner product of two vectors
template<typename T>
T operator&(thrust::device_vector<T>& a, thrust::device_vector<T>& b)
{
	return thrust::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}

/// \brief Operator overload
/// \param a input vector a
/// \param b input vector b
/// \return a ./ b
template<typename T>
thrust::device_vector<T> operator/(thrust::device_vector<T>& a, thrust::device_vector<T>& b)
{
	thrust::device_vector<T> c(a.size(), 0);
	thrust::transform(a.begin(), a.end(), b.begin(), c.begin(), thrust::divides<T>());
	return c;
}


template<typename T>
thrust::device_vector<T>& operator+=(thrust::device_vector<T>& lhs, const thrust::device_vector<T>& rhs)
{
	thrust::transform(lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(), thrust::plus<T>());
	return lhs;
}

template<typename T>
thrust::device_vector<T>& operator-=(thrust::device_vector<T>& lhs, const thrust::device_vector<T>& rhs)
{
	thrust::transform(lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(), thrust::minus<T>());
	return lhs;
}

template<typename T>
thrust::device_vector<T>& operator*=(thrust::device_vector<T>& lhs, const thrust::device_vector<T>& rhs)
{
	thrust::transform(lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(), thrust::multiplies<T>());
	return lhs;
}

template<typename T>
thrust::device_vector<T>& operator/=(thrust::device_vector<T>& lhs, const thrust::device_vector<T>& rhs)
{
	thrust::transform(lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(), thrust::divides<T>());
	return lhs;
}

template<typename T>
thrust::device_vector<T>& operator*=(thrust::device_vector<T>& lhs, const T& rhs)
{
	thrust::transform(lhs.begin(), lhs.end(), lhs.begin(), thrust::multiplies<T>());
	return lhs;
}



INLINE __host__ __device__ float fminf(const float2& a){	return fminf(a.x, a.y);}
INLINE __host__ __device__ float fminf(const float3& a){	return fminf(a.x, fminf(a.y, a.z)); }
INLINE __host__ __device__ float2 fminf(const float2& a, const float2& b){ return make_float2(fmin(a.x, b.x), fmin(a.y, b.y)); }
INLINE __host__ __device__ float3 fminf(const float3& a, const float3& b){	return make_float3(fminf(a.x, b.x), fminf(a.y, b.y), fminf(a.z, b.z));}
INLINE __host__ __device__ double2 fminf(const double2& a, const double2& b){ return make_double2(fmin(a.x, b.x), fmin(a.y, b.y)); }
INLINE __host__ __device__ double2 fminf(const float2& a, const double2& b){ return make_double2(fminf(a.x, b.x), fminf(a.y, b.y)); }
INLINE __host__ __device__ double2 fminf(const double2& a, const float2& b){ return make_double2(fminf(a.x, b.x), fminf(a.y, b.y)); }
INLINE __host__ __device__ double3 fminf(const double3& a, const double3& b){ return make_double3(fminf(a.x, b.x), fminf(a.y, b.y), fminf(a.z, b.z)); }
INLINE __host__ __device__ double3 fminf(const double3& a, const float3& b){ return make_double3(fminf(a.x, b.x), fminf(a.y, b.y), fminf(a.z, b.z)); }
INLINE __host__ __device__ double3 fminf(const float3& a, const double3& b){ return make_double3(fminf(a.x, b.x), fminf(a.y, b.y), fminf(a.z, b.z)); }
INLINE __host__ __device__ float fmaxf(const float2& a){	return fmaxf(a.x, a.y);}
INLINE __host__ __device__ float fmaxf(const float3& a){	return fmaxf(a.x, fmaxf(a.y, a.z));}
INLINE __host__ __device__ float2 fmaxf(const float2& a, const float2& b){ return make_float2(fmax(a.x, b.x), fmax(a.y, b.y)); }
INLINE __host__ __device__ float3 fmaxf(const float3& a, const float3& b){	return make_float3(fmaxf(a.x, b.x), fmaxf(a.y, b.y), fmaxf(a.z, b.z));}
INLINE __host__ __device__ double2 fmaxf(const double2& a, const double2& b){ return make_double2(fmaxf(a.x, b.x), fmaxf(a.y, b.y)); }
INLINE __host__ __device__ double2 fmaxf(const float2& a, const double2& b){ return make_double2(fmaxf(a.x, b.x), fmaxf(a.y, b.y)); }
INLINE __host__ __device__ double2 fmaxf(const double2& a, const float2& b){ return make_double2(fmaxf(a.x, b.x), fmaxf(a.y, b.y)); }
INLINE __host__ __device__ double3 fmaxf(const double3& a, const double3& b){ return make_double3(fmaxf(a.x, b.x), fmaxf(a.y, b.y), fmaxf(a.z, b.z)); }
INLINE __host__ __device__ double3 fmaxf(const double3& a, const float3& b){ return make_double3(fmaxf(a.x, b.x), fmaxf(a.y, b.y), fmaxf(a.z, b.z)); }
INLINE __host__ __device__ double3 fmaxf(const float3& a, const double3& b){ return make_double3(fmaxf(a.x, b.x), fmaxf(a.y, b.y), fmaxf(a.z, b.z)); }

INLINE __host__ __device__ float length(const float2& a){	return sqrtf(a.x * a.x + a.y * a.y);}
INLINE __host__ __device__ float length(const float3& a){	return sqrtf(a.x * a.x + a.y * a.y + a.z * a.z);}
INLINE __host__ __device__ double length(const double2& a) { return sqrt(a.x * a.x + a.y * a.y); }
INLINE __host__ __device__ double length(const double3& a){	return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);}

//INLINE __host__ __device__ const float2 normalize(const float2& a){	return a / length(a);}
INLINE __host__ __device__ const float2 normalize(float2 a) { return a / length(a); }
INLINE __host__ __device__ const float3 normalize(const float3& a){	return a / length(a);}
INLINE __host__ __device__ const double2 normalize(const double2& a) { return a / length(a); }
INLINE __host__ __device__ const double3 normalize(const double3& a){	return a / length(a);}




/// \brief square of input value
/// \param a input value
/// \return a * a
template<typename T>__host__ __device__ inline T MY_SQUARE(const T& a){ return a*a; }

/// \brief max value of input value
/// \param a input value
/// \param b input value
/// \return max(a,b)
template<typename T>__host__ __device__ inline const T &MY_MAX(const T& a, const T& b){ return b > a ? b : a; }

/// \brief max value of input value
/// \param a input value
/// \param b input value
/// \return max(a,b)
template<typename T, typename T1, typename T2>__host__ __device__ inline T MY_MAX(const T1& a, const T2& b){	return static_cast<T>(a) > static_cast<T>(b) ? static_cast<T>(a) : static_cast<T>(b);}

/// \brief max value of input value
/// \param a input value
/// \param b input value
/// \return max(a,b)
inline __host__ __device__ float MY_MAX(const double &a, const float &b){ return b > a ? (b) : float(a); }

/// \brief max value of input value
/// \param a input value
/// \param b input value
/// \return max(a,b)
inline __host__ __device__ float MY_MAX(const float &a, const double &b){ return b > a ? float(b) : (a); }

/// \brief max value of input value
/// \param a input value
/// \param b input value
/// \param c input value
/// \return max(max(a,b),c)
template<typename T>__host__ __device__ inline T MY_MAX(const T& a, const T& b, const T& c){	return MY_MAX<T>(MY_MAX<T>(a, b), c);}

/// \brief min value of input value
/// \param a input value
/// \param b input value
/// \return min(a,b)
template<typename T>__host__ __device__ inline const T &MY_MIN(const T& a, const T& b){ return b < a ? b : a; }

/// \brief min value of input value
/// \param a input value
/// \param b input value
/// \return min(a,b)
template<typename T, typename T1, typename T2>__host__ __device__ inline T MY_MIN(const T1& a, const T2& b){	return static_cast<T>(a) < static_cast<T>(b) ? static_cast<T>(a) : static_cast<T>(b);}

/// \brief min value of input value
/// \param a input value
/// \param b input value
/// \return min(a,b)
inline __host__ __device__ float MY_MIN(const double &a, const float &b){ return b < a ? (b) : float(a); }
/// \brief min value of input value
/// \param a input value
/// \param b input value
/// \return min(a,b)
inline __host__ __device__ float MY_MIN(const float &a, const double &b){ return b < a ? float(b) : (a); }
/// \brief min value of input value
/// \param a input value
/// \param b input value
/// \param c input value
/// \return min(min(a,b),c)
template<typename T>__host__ __device__ inline T MY_MIN(const T& a, const T& b, const T& c){	return MY_MIN<T>(MY_MIN<T>(a, b), c);}
template<typename T, typename T1, typename T2>inline __host__ __device__ T MY_SIGN(const T1& a, const T2& b){ return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a); }

/// \brief exchange value of a and b
/// \param a parameter a
/// \param b parameter b
template<class T> inline __host__ __device__ void MY_SWAP(T &a, T &b) { T dum = a; a = b; b = dum; }

/// \brief return the absolute value of input value
/// \param x input value
/// \return abs(x)
template<typename T>inline __host__ __device__ T MY_ABS(const T& x){	return (x > 0) ? x : (-x);}


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
	while (c >= TWOPI){ c -= TWOPI; }
	while (c < 0){ c += TWOPI; }
	return c;
}

namespace CTMBIR {
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

/// \brief the bilinear interpolation of the data
/// \param x input x coordinate
/// \param y input y coordinate
/// \param x1 small x coordinate
/// \param x2 big x coordinate
/// \param y1 small y coordinate
/// \param y2 big y coordinate
/// \param Qmm P(x1,y1)
/// \param Qmb P(x2,y1)
/// \param Qbm P(x1,y2)
/// \param Qbb P(x2,y2)
template<typename T>
inline __host__ __device__  T _bilinearInterpolation_ker(
	T x, T y,
	T x1, T x2,
	T y1, T y2,
	T Qmm, T Qmb,
	T Qbm, T Qbb) {
	if (IS_ZERO(x1 - x2) && IS_ZERO(y1 - y2))
		return (Qmm + Qmb + Qbm + Qbb) * 0.25f;
	else if (IS_ZERO(x1 - x2))
		return ((y2 - y) * (Qmm + Qbm) + (y - y1) * (Qmb + Qbb)) * 0.5f;
	else if (IS_ZERO(y1 - y2))
		return ((x2 - x) * (Qmm + Qmb) + (x - x1) * (Qbm + Qbb)) * 0.5f;
	else
	{
		const T invs(1.0f / ((x2 - x1)*(y2 - y1)));
		return invs * (Qmm * (x2 - x) * (y2 - y) + Qbm * (x - x1) * (y2 - y)
			+ Qmb * (x2 - x) * (y - y1) + Qbb * (x - x1) * (y - y1));
	}
}


/// \brief the trilinear interpolation of the data
template<typename T>
__host__ __device__ inline T _trilinearInterpolation(
	T x, T y, T z,
	T x1, T x2,
	T y1, T y2,
	T z1, T z2,
	T Qmmm, T Qmmb, T Qmbm, T Qmbb,
	T Qbmm, T Qbmb, T Qbbm, T Qbbb)
{
	T upZ = _bilinearInterpolation_ker<T>(x, y, x1, x2, y1, y2, Qmmm, Qmmb, Qmbm, Qmbb);
	T dwZ = _bilinearInterpolation_ker<T>(x, y, x1, x2, y1, y2, Qbmm, Qbmb, Qbbm, Qbbb);

	if (z1 == z2)
	{
		return (upZ + dwZ) * 0.5f;
	}
	return upZ * (z - z1) + dwZ * (z2 - z);
}



template<typename T>
inline	__host__ __device__ T calSiddonOneRayKer2D(
	const T& startX, const T& startY,
	const T& endX, const T& endY,
	const T& __MINOL__, const T& __MINOW__,
	const T& __OBJSTPL__, const T& __OBJSTPW__,
	const unsigned int& __OBJLR__, const unsigned int& __OBJWR__,
	T* dev_vol, T* totWeight)
{
	T dirX = endX - startX;
	T dirY = endY - startY;
	const T dconv = sqrt(dirX * dirX + dirY * dirY);
	int imin(0), imax(0), jmin(0), jmax(0);

	T alphaxmin = 0.0f;
	T alphaxmax = 0.0f;
	T alphaymin = 0.0f;
	T alphaymax = 0.0f;

	dirX = dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, 0);
	dirY = dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, __OBJLR__);
	if (dirX < dirY) {
		alphaxmin = dirX;
		alphaxmax = dirY;
	} else {
		alphaxmin = dirY;
		alphaxmax = dirX;
	}

	dirX = dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, 0);
	dirY = dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, __OBJWR__);
	if (dirX < dirY) {
		alphaymin = dirX;
		alphaymax = dirY;
	} else {
		alphaymin = dirY;
		alphaymax = dirX;
	}

	const T alphaMIN = MY_MAX<T>(alphaxmin, alphaymin); // (alphaxmin > alphaymin) ? alphaxmin : alphaymin;
	const T alphaMAX = MY_MIN<T>(alphaxmax, alphaymax); // (alphaxmax < alphaymax) ? alphaxmax : alphaymax;
	dev_minmaxIdxFun(startX, endX, __MINOL__, __OBJSTPL__, alphaMIN, alphaMAX, alphaxmin, alphaxmax, __OBJLR__ + 1, &imin, &imax);
	dev_minmaxIdxFun(startY, endY, __MINOW__, __OBJSTPW__, alphaMIN, alphaMAX, alphaymin, alphaymax, __OBJWR__ + 1, &jmin, &jmax);

	T alphaX = (startX < endX) ? dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, imin) : dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, imax);
	T alphaY = (startY < endY) ? dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, jmin) : dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, jmax);

	int Np = static_cast<int>(MY_ABS<float>(static_cast<T>(imax - imin)) + MY_ABS<float>(static_cast<T>(jmax - jmin)) + 3.0f); //fabsf(imax - imin) + fabsf(jmax - jmin) + 3.0f;
	const T alphaxu = dev_alphaU_Fun(__OBJSTPL__, startX, endX);
	const T alphayu = dev_alphaU_Fun(__OBJSTPW__, startY, endY);

	T alphaC = alphaMIN;
	//T minalpha = MY_MIN<T>(alphaX,alphaY); //min(alphaX, alphaY);
	//T talpha = MY_MIN<T>(alphaX,alphaY); //fminf(alphaX,alphaY);
	int i = int(dev_varphiFun(alphaMIN * 1.00003f, __MINOL__, __OBJSTPL__, startX, endX));
	int j = int(dev_varphiFun(alphaMIN * 1.00003f, __MINOW__, __OBJSTPW__, startY, endY));
	//int i = floor(dev_varphiFun((talpha + alphaMIN)*0.5, __MINOL__, __OBJSTPL__, startX, endX));
	//int j = floor(dev_varphiFun((talpha + alphaMIN)*0.5, __MINOW__, __OBJSTPW__, startY, endY));

	const int iuu = (startX < endX) ? 1 : -1;
	const int juu = (startY < endY) ? 1 : -1;

	T d12(0);
	T weight(0);


	for (int repIdx(0); repIdx < Np; ++repIdx) {
		if (i < 0 || i >= static_cast<int>(__OBJLR__) || j < 0 || j >= static_cast<int>(__OBJWR__)) {
			break;
		}
		if (alphaX <= alphaY) {
			weight = (alphaX - alphaC) * dconv;
			d12 += weight * dev_vol[j * __OBJLR__ + i];
			i += iuu;
			(*totWeight) += weight;
			alphaC = alphaX;
			alphaX += alphaxu;
			continue;
		} else {
			weight = (alphaY - alphaC) * dconv;
			d12 += weight * dev_vol[j * __OBJLR__ + i];
			j += juu;
			(*totWeight) += weight;
			alphaC = alphaY;
			alphaY += alphayu;
			continue;
		}
	}
	return d12;
}

template<typename T>
void preCalculateSiddon(const T& startX, const T& startY, const T& startZ,
	const T& endX, const T& endY, const T& endZ,
	const T& __MINOL__, const T& __MINOW__, const T& __MINOH__,
	const T& __OBJSTPL__, const T& __OBJSTPW__, const T& __OBJSTPH__,
	const unsigned int& __OBJLR__, const unsigned int& __OBJWR__, const unsigned int& __OBJHR__,
	T* dev_vol, T* totWeight) {
	T dirX = endX - startX;
	T dirY = endY - startY;
	T dirZ = endZ - startZ;

	const T dconv = sqrtf(dirX * dirX + dirY * dirY + dirZ * dirZ);
	int imin(0), imax(0), jmin(0), jmax(0), kmin(0), kmax(0);

	thrust::pair<T, T> alpX = calculateAlpha(__MINOL__, __OBJSTPL__, startX, endX, __OBJLR__);
	thrust::pair<T, T> alpY = calculateAlpha(__MINOW__, __OBJSTPW__, startY, endY, __OBJWR__);
	thrust::pair<T, T> alpZ = calculateAlpha(__MINOH__, __OBJSTPH__, startZ, endZ, __OBJHR__);
	T alphaxmin = alpX.first;
	T alphaxmax = alpX.second;
	T alphaymin = alpY.first;
	T alphaymax = alpY.second;
	T alphazmin = alpZ.first;
	T alphazmax = alpZ.second;


	const T alphaMIN = MY_MAX<T>(alphaxmin, alphaymin, alphazmin);
	const T alphaMAX = MY_MIN<T>(alphaxmax, alphaymax, alphazmax);
	dev_minmaxIdxFun(startX, endX, __MINOL__, __OBJSTPL__, alphaMIN, alphaMAX, alphaxmin, alphaxmax, __OBJLR__ + 1, &imin, &imax);
	dev_minmaxIdxFun(startY, endY, __MINOW__, __OBJSTPW__, alphaMIN, alphaMAX, alphaymin, alphaymax, __OBJWR__ + 1, &jmin, &jmax);
	dev_minmaxIdxFun(startZ, endZ, __MINOH__, __OBJSTPH__, alphaMIN, alphaMAX, alphazmin, alphazmax, __OBJHR__ + 1, &kmin, &kmax);

	T alphaX = (startX < endX) ? dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, imin) : dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, imax);
	T alphaY = (startY < endY) ? dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, jmin) : dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, jmax);
	T alphaZ = (startZ < endZ) ? dev_alpha_IFun(__MINOH__, __OBJSTPH__, startZ, endZ, kmin) : dev_alpha_IFun(__MINOH__, __OBJSTPH__, startZ, endZ, kmax);

	int Np = static_cast<int>(MY_ABS<T>(static_cast<T>(imax - imin)) + MY_ABS<T>(static_cast<T>(jmax - jmin)) + MY_ABS<T>(static_cast<T>(kmax - kmin)) + 3.0f);
	const T alphaxu = dev_alphaU_Fun(__OBJSTPL__, startX, endX);
	const T alphayu = dev_alphaU_Fun(__OBJSTPW__, startY, endY);
	const T alphazu = dev_alphaU_Fun(__OBJSTPH__, startZ, endZ);

	T alphaC = alphaMIN;
	//const T minApa = MY_MIN<T>(alphaX,alphaY,alphaZ);

	int i = int(dev_varphiFun(alphaMIN * 1.00003f, __MINOL__, __OBJSTPL__, startX, endX));
	int j = int(dev_varphiFun(alphaMIN * 1.00003f, __MINOW__, __OBJSTPW__, startY, endY));
	int k = int(dev_varphiFun(alphaMIN * 1.00003f, __MINOH__, __OBJSTPH__, startZ, endZ));
	//int i = floor(dev_varphiFun((alphaMIN+minApa) * 0.5f, __MINOL__, __OBJSTPL__, startX, endX));
	//int j = floor(dev_varphiFun((alphaMIN+minApa) * 0.5f, __MINOW__, __OBJSTPW__, startY, endY));
	//int k = floor(dev_varphiFun((alphaMIN+minApa) * 0.5f, __MINOH__, __OBJSTPH__, startZ, endZ));

	const int iuu = dev_iu_Fun(startX, endX);
	const int juu = dev_iu_Fun(startY, endY);
	const int kuu = dev_iu_Fun(startZ, endZ);
}

template<typename T>
__device__ thrust::pair<T, T> calculateAlpha(T __MINO__, T __OBJSTP__, T start, T end, const unsigned int __OBJ__) {
	T dirX = dev_alpha_IFun(__MINO__, __OBJSTP__, start, end, 0);
	T dirY = dev_alpha_IFun(__MINO__, __OBJSTP__, start, end, __OBJ__);
	if (dirX < dirY) {
		return { dirX, dirY };
	} else {
		return { dirY, dirX };
	}
}

template<typename T>
inline __device__ T calSiddonOneRayKer(const T& startX, const T& startY, const T& startZ,
	const T& endX, const T& endY, const T& endZ,
	const T& __MINOL__, const T& __MINOW__, const T& __MINOH__,
	const T& __OBJSTPL__, const T& __OBJSTPW__, const T& __OBJSTPH__,
	const unsigned int& __OBJLR__, const unsigned int& __OBJWR__, const unsigned int& __OBJHR__,
	T* dev_vol, T* totWeight)
{

	T dirX = endX - startX;
	T dirY = endY - startY;
	T dirZ = endZ - startZ;

	const T dconv = sqrtf(dirX * dirX + dirY * dirY + dirZ * dirZ);
	int imin(0), imax(0), jmin(0), jmax(0), kmin(0), kmax(0);

	thrust::pair<T, T> alpX = calculateAlpha(__MINOL__, __OBJSTPL__, startX, endX, __OBJLR__);
	thrust::pair<T, T> alpY = calculateAlpha(__MINOW__, __OBJSTPW__, startY, endY, __OBJWR__);
	thrust::pair<T, T> alpZ = calculateAlpha(__MINOH__, __OBJSTPH__, startZ, endZ, __OBJHR__);
	T alphaxmin = alpX.first;
	T alphaxmax = alpX.second;
	T alphaymin = alpY.first;
	T alphaymax = alpY.second;
	T alphazmin = alpZ.first;
	T alphazmax = alpZ.second;


	const T alphaMIN = MY_MAX<T>(alphaxmin, alphaymin, alphazmin);
	const T alphaMAX = MY_MIN<T>(alphaxmax, alphaymax, alphazmax);
	dev_minmaxIdxFun(startX, endX, __MINOL__, __OBJSTPL__, alphaMIN, alphaMAX, alphaxmin, alphaxmax, __OBJLR__ + 1, &imin, &imax);
	dev_minmaxIdxFun(startY, endY, __MINOW__, __OBJSTPW__, alphaMIN, alphaMAX, alphaymin, alphaymax, __OBJWR__ + 1, &jmin, &jmax);
	dev_minmaxIdxFun(startZ, endZ, __MINOH__, __OBJSTPH__, alphaMIN, alphaMAX, alphazmin, alphazmax, __OBJHR__ + 1, &kmin, &kmax);

	T alphaX = (startX < endX) ? dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, imin) : dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, imax);
	T alphaY = (startY < endY) ? dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, jmin) : dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, jmax);
	T alphaZ = (startZ < endZ) ? dev_alpha_IFun(__MINOH__, __OBJSTPH__, startZ, endZ, kmin) : dev_alpha_IFun(__MINOH__, __OBJSTPH__, startZ, endZ, kmax);

	int Np = static_cast<int>(MY_ABS<T>(static_cast<T>(imax - imin)) + MY_ABS<T>(static_cast<T>(jmax - jmin)) + MY_ABS<T>(static_cast<T>(kmax - kmin)) + 3.0f);
	const T alphaxu = dev_alphaU_Fun(__OBJSTPL__, startX, endX);
	const T alphayu = dev_alphaU_Fun(__OBJSTPW__, startY, endY);
	const T alphazu = dev_alphaU_Fun(__OBJSTPH__, startZ, endZ);

	T alphaC = alphaMIN;
	//const T minApa = MY_MIN<T>(alphaX,alphaY,alphaZ);

	int i = int(dev_varphiFun(alphaMIN * 1.00003f, __MINOL__, __OBJSTPL__, startX, endX));
	int j = int(dev_varphiFun(alphaMIN * 1.00003f, __MINOW__, __OBJSTPW__, startY, endY));
	int k = int(dev_varphiFun(alphaMIN * 1.00003f, __MINOH__, __OBJSTPH__, startZ, endZ));
	//int i = floor(dev_varphiFun((alphaMIN+minApa) * 0.5f, __MINOL__, __OBJSTPL__, startX, endX));
	//int j = floor(dev_varphiFun((alphaMIN+minApa) * 0.5f, __MINOW__, __OBJSTPW__, startY, endY));
	//int k = floor(dev_varphiFun((alphaMIN+minApa) * 0.5f, __MINOH__, __OBJSTPH__, startZ, endZ));

	const int iuu = dev_iu_Fun(startX, endX);
	const int juu = dev_iu_Fun(startY, endY);
	const int kuu = dev_iu_Fun(startZ, endZ);

	T d12(0);
	T weight(0);


	for (int repIdx(0); repIdx < Np; ++repIdx) {
		if (i < 0 || i > static_cast<int>(__OBJLR__ - 1) || j < 0 || j > static_cast<int>(__OBJWR__ - 1) || k < 0 || k > static_cast<int>(__OBJHR__ - 1)) {
			break;
		}

		if (alphaX <= alphaY && alphaX <= alphaZ) {
			weight = (alphaX - alphaC) * dconv;
			if (weight > 0) {
				d12 = d12 + weight * dev_vol[(k * __OBJWR__ + j) * __OBJLR__ + i];
				(*totWeight) += weight;
			}

			i = i + iuu;
			alphaC = alphaX;
			alphaX = alphaX + alphaxu;
		} else if (alphaY <= alphaX && alphaY <= alphaZ) {
			weight = (alphaY - alphaC) * dconv;
			if (weight > 0) {
				d12 = d12 + weight * dev_vol[(k * __OBJWR__ + j) * __OBJLR__ + i];
				(*totWeight) += weight;
			}

			j = j + juu;
			alphaC = alphaY;
			alphaY = alphaY + alphayu;
		} else {
			weight = (alphaZ - alphaC) * dconv;
			if (weight > 0) {
				d12 = d12 + weight * dev_vol[(k * __OBJWR__ + j) * __OBJLR__ + i];
				(*totWeight) += weight;
			}

			k = k + kuu;
			alphaC = alphaZ;
			alphaZ = alphaZ + alphazu;
		}
	}

	return d12;
}






/// \brief intersection length of the rectangle box
/// \param r ray
/// \param boxmin min boundary of the box
/// \param boxmax max boundary of the obx
/// \param tnear pointer to the near intersection parameter
/// \param tfar pointer to the far intersection parameter
/// \return intersection or not
inline	__host__ __device__ bool intersectBox(Ray2D r, float2 boxmin, float2 boxmax, float *tnear, float *tfar)
{
	// compute intersection of ray with all six bbox planes
	float2 invR = make_float2(1.0f, 1.0f) / r.d;
	float2 tbot = invR * (boxmin - r.o);
	float2 ttop = invR * (boxmax - r.o);

	// re-order intersections to find smallest and largest on each axis
	float2 tmin = fminf(ttop, tbot);
	float2 tmax = fmaxf(ttop, tbot);

	// find the largest tmin and the smallest tmax
	float largest_tmin = fmaxf(tmin.x, tmin.y);
	float smallest_tmax = fminf(tmax.x, tmax.y);

	*tnear = largest_tmin;
	*tfar = smallest_tmax;

	return smallest_tmax > largest_tmin;
}

inline	__host__ __device__ bool intersectBox(Ray2D r, double2 boxmin, double2 boxmax, double *tnear, double *tfar)
{
	// compute intersection of ray with all six bbox planes
	double2 invR = make_double2(1.0f, 1.0f) / r.d;
	double2 tbot = invR * (boxmin - r.o);
	double2 ttop = invR * (boxmax - r.o);

	// re-order intersections to find smallest and largest on each axis
	double2 tmin = fminf(ttop, tbot);
	double2 tmax = fmaxf(ttop, tbot);

	// find the largest tmin and the smallest tmax
	double largest_tmin = fmax(tmin.x, tmin.y);
	double smallest_tmax = fmin(tmax.x, tmax.y);

	*tnear = largest_tmin;
	*tfar = smallest_tmax;

	return smallest_tmax > largest_tmin;
}


/// \brief intersection length of the 3D rectangle box
/// \param r ray
/// \param boxmin min boundary of the box
/// \param boxmax max boundary of the obx
/// \param tnear pointer to the near intersection parameter
/// \param tfar pointer to the far intersection parameter
/// \return intersection or not
inline __host__ __device__ int intersectBox(Ray r, float3 boxmin, float3 boxmax, float *tnear, float *tfar)
{
	// compute intersection of ray with all six bbox planes
	float3 invR = make_float3(1.0f,1.0f,1.0f) / r.d;
	float3 tbot = invR * (boxmin - r.o);
	float3 ttop = invR * (boxmax - r.o);

	// re-order intersections to find smallest and largest on each axis
	float3 tmin = fminf(ttop, tbot);
	float3 tmax = fmaxf(ttop, tbot);

	// find the largest tmin and the smallest tmax
	float largest_tmin = fmaxf(fmaxf(tmin.x, tmin.y), fmaxf(tmin.x, tmin.z));
	float smallest_tmax = fminf(fminf(tmax.x, tmax.y), fminf(tmax.x, tmax.z));

	*tnear = largest_tmin;
	*tfar = smallest_tmax;

	return smallest_tmax > largest_tmin;
}

inline __host__ __device__ int intersectBox(Ray r, double3 boxmin, double3 boxmax, double *tnear, double *tfar)
{
	// compute intersection of ray with all six bbox planes
	double3 invR = make_double3(1.0f, 1.0f, 1.0f) / r.d;
	double3 tbot = invR * (boxmin - r.o);
	double3 ttop = invR * (boxmax - r.o);

	// re-order intersections to find smallest and largest on each axis
	double3 tmin = fminf(ttop, tbot);
	double3 tmax = fmaxf(ttop, tbot);

	// find the largest tmin and the smallest tmax
	double largest_tmin = fmax(fmax(tmin.x, tmin.y), fmax(tmin.x, tmin.z));
	double smallest_tmax = fmin(fmin(tmax.x, tmax.y), fmin(tmax.x, tmax.z));

	*tnear = largest_tmin;
	*tfar = smallest_tmax;

	return smallest_tmax > largest_tmin;
}

template<typename T>
__host__ __device__ inline T intersectLength(const T& fixedmin, const T& fixedmax, const T& varimin, const T& varimax)
{
	const T left = (fixedmin > varimin) ? fixedmin : varimin;
	const T right = (fixedmax < varimax) ? fixedmax : varimax;
	return abs(right - left) * (right > left);
}

/// \brief the core function for generating the projection matrix
/// \param startX the start point x of the ray
/// \param startY the start point y of the ray
/// \param endX the end point x of the ray
/// \param endY the end point y of the ray
/// \param bx the minimum x coordinate of the object
/// \param by the minimum y coordinate of the object
/// \param dx the obj x direction step size
/// \param dy the obj y direction step size
/// \param objResoLen the object length resolution
/// \param objResoWid the object width resolution
/// \param rowidx row index
/// \param rowIdx row index in the matrix
/// \param colIdx col index in the matrix
/// \param wgt weighting coefficient of the matrix
template<typename T>
inline void pushMatrix(T startX, T startY, T endX, T endY,
	T bx, T by,T dx, T dy, const unsigned int objResoLen, const unsigned int  objResoWid,
	int& rowidx,
	std::vector<int>& rowIdx,
	std::vector<int>& colIdx,
	std::vector<T>& wgt) {

	const T dirX(endX - startX);
	const T dirY(endY - startY);
	const T lengthSq = dirX * dirX + dirY * dirY;
	const T dconv = sqrt(lengthSq);
	int imin, imax, jmin, jmax;
	const T alphaxmin = MY_MIN<T>(dev_alpha_IFun(bx, dx, startX, endX, 0), dev_alpha_IFun(bx, dx, startX, endX, objResoLen));
	const T alphaxmax = MY_MAX<T>(dev_alpha_IFun(bx, dx, startX, endX, 0), dev_alpha_IFun(bx, dx, startX, endX, objResoLen));
	const T alphaymin = MY_MIN<T>(dev_alpha_IFun(by, dy, startY, endY, 0), dev_alpha_IFun(by, dy, startY, endY, objResoWid));
	const T alphaymax = MY_MAX<T>(dev_alpha_IFun(by, dy, startY, endY, 0), dev_alpha_IFun(by, dy, startY, endY, objResoWid));

	const T alphaMIN = MY_MAX<T>(alphaxmin, alphaymin);
	const T alphaMAX = MY_MIN<T>(alphaxmax, alphaymax);
	dev_minmaxIdxFun(startX, endX, bx, dx, alphaMIN, alphaMAX, alphaxmin, alphaxmax, objResoLen + 1, &imin, &imax);
	dev_minmaxIdxFun(startY, endY, by, dy, alphaMIN, alphaMAX, alphaymin, alphaymax, objResoWid + 1, &jmin, &jmax);

	T alphaX = (startX < endX) ? dev_alpha_IFun(bx, dx, startX, endX, imin) : dev_alpha_IFun(bx, dx, startX, endX, imax);
	T alphaY = (startY < endY) ? dev_alpha_IFun(by, dy, startY, endY, jmin) : dev_alpha_IFun(by, dy, startY, endY, jmax);

	int Np = static_cast<int>(MY_ABS<T>(imax - imin + 1.0f) + MY_ABS<T>(jmax - jmin + 1.0f) + 4.0f);
	const T alphaxu = dev_alphaU_Fun(dx, startX, endX);
	const T alphayu = dev_alphaU_Fun(dy, startY, endY);

	T alphaC = alphaMIN;
	T minApa = MY_MIN<T>(alphaX, alphaY);
	/*int i = static_cast<int>(dev_varphiFun(alphaMIN* 1.00003f, bx, dx, startX, endX));
	int j = static_cast<int>(dev_varphiFun(alphaMIN* 1.00003f, by, dy, startY, endY));
	*/
	int i = static_cast<int>(dev_varphiFun((alphaMIN + minApa) * 0.5f, bx, dx, startX, endX));
	int j = static_cast<int>(dev_varphiFun((alphaMIN + minApa) * 0.5f, by, dy, startY, endY));

	const int iuu = dev_iu_Fun(startX, endX);
	const int juu = dev_iu_Fun(startY, endY);

	T d12(0.0f);
	T weight(0.0f);
	unsigned int repIdx(0);
	unsigned int colidx(0);
	while (repIdx != Np) {
		if (i < 0 || i >= static_cast<int>(objResoLen) || j < 0 || j >= static_cast<int>(objResoWid)) {
			break;
		}

		colidx = j * objResoLen + i;
		weight = (((alphaX <= alphaY) ? alphaX : alphaY) - alphaC) * dconv;
		wgt.push_back(weight);
		rowIdx.push_back(rowidx);
		colIdx.push_back(colidx);
		d12 += weight;

		if (alphaX <= alphaY) {
			i += iuu;
			alphaC = alphaX;
			alphaX += alphaxu;
		} else {
			j += juu;
			alphaC = alphaY;
			alphaY += alphayu;
		}

		++repIdx;
	}
}

// Calculate the intersection point coordinates
template<typename T>
inline __host__ __device__ void calSV(T initX, T initY,const T& cosT,const T& sinT, T* sour, T* SVA, const unsigned int detIdx, T plusOne) {
	T curX = initX * cosT - initY * sinT;
	T curY = initX * sinT + initY * cosT;
	T dx = curX - sour[0];
	T dy = curY - sour[1];
	T legth = sqrt(dx * dx + dy * dy);
	SVA[0] = dx / legth;
	SVA[1] = dy / legth;
	SVA[2] = static_cast<T>(detIdx) + plusOne;
}


template <typename T>
inline __host__ __device__ void calSVASVB(T* SVA, T* SVB, T* sour, const T& cosT, const T& sinT, const FanEAGeo& FanGeo, const Image& Img, cuint detIdx) {
	T pangle = (detIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
	T initY = -cos(pangle) * FanGeo.m_S2D + FanGeo.m_S2O;
	T initX = sin(pangle) * FanGeo.m_S2D;

	calSV<T>(initX, initY, cosT, sinT, sour, SVA, detIdx, 0.0);

	pangle = pangle + FanGeo.m_DetStp;
	initY = -cos(pangle) * FanGeo.m_S2D + FanGeo.m_S2O;
	initX = sin(pangle) * FanGeo.m_S2D;

	calSV<T>(initX, initY, cosT, sinT, sour, SVB, detIdx, 1.0);
}


template<typename T>
__host__ __device__ void calSVASVB(T* SVA, T* SVB, T* sour, const T& cosT, const T& sinT, const FanEDGeo& FanGeo, const Image& Img, cuint detIdx)
{
	T initX = (detIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
	T initY = -FanGeo.m_O2D;

	calSV<T>(initX, initY, cosT, sinT, sour, SVA, detIdx, 0.0);

	initX = (detIdx - FanGeo.m_DetCntIdx + 1) * FanGeo.m_DetStp;
	initY = -FanGeo.m_O2D;

	calSV<T>(initX, initY, cosT, sinT, sour, SVB, detIdx, 1.0);
}

template<typename T>
__host__ __device__ void dev_swap(T& a, T& b) {
	T c = a;
	a = b;
	b = c;
}

template<typename T>
__host__ __device__ inline void SortProjection(T(&Grid)[4][3])
{
	int i, j;
	T td;
	for (i = 0; i < 3; i++) {
		for (j = i + 1; j < 4; j++) {
			if (Grid[j][2] < Grid[i][2]) {
				dev_swap(Grid[i][0], Grid[j][0]);
				dev_swap(Grid[i][1], Grid[j][1]);
				dev_swap(Grid[i][2], Grid[j][2]);
			}
		}
	}
}



template<typename T>
inline __host__ __device__ T ComputeCoefficient(const T(&Grid)[4][3], const T(&SVA)[3],	const T(&SVB)[3], const T(&SPoint)[2], const T area)
{
	T coef = 0;
	T x0(0), y0(0), a(0), b(0), t(0);
	int AI, BI;
	if (SVA[2] > Grid[3][2]) {
		return 0;
	}
	if (SVB[2] < Grid[0][2]) {
		return 0;
	}

	if (SVA[2] < Grid[0][2])
		AI = 0;
	else if (SVA[2] < Grid[1][2])
		AI = 1;
	else if (SVA[2] < Grid[2][2])
		AI = 2;
	else AI = 3;

	if (SVB[2] < Grid[1][2])
		BI = 1;
	else if (SVB[2] < Grid[2][2])
		BI = 2;
	else if (SVB[2] < Grid[3][2])
		BI = 3;
	else BI = 4;

	switch (AI)	{
	case 0:
	{
		switch (BI)	{
			case 1:// case [0,1]
			{
				x0 = Grid[0][0] - SPoint[0];
				y0 = Grid[0][1] - SPoint[1];
				if (abs(SVB[0] * SVB[1]) > 0)
				{
					t = x0*SVB[1] - y0*SVB[0];
					a = t / SVB[0];
					b = t / SVB[1];
					coef = 0.5*abs(a*b);
				}
				break;
			}
			case 2: // case [0,2]
			{
				x0 = abs(Grid[0][0] - Grid[1][0]);
				y0 = abs(Grid[0][1] - Grid[1][1]);
				if (x0 > y0) // line is on the x-direction
				{
					a = abs((Grid[0][0] - SPoint[0])*SVB[1] / SVB[0] - (Grid[0][1] - SPoint[1]));
					b = abs((Grid[1][0] - SPoint[0])*SVB[1] / SVB[0] - (Grid[1][1] - SPoint[1]));
					coef = (a + b)*x0*0.5;
				}
				else
				{
					a = abs((Grid[0][0] - SPoint[0]) - (Grid[0][1] - SPoint[1])*SVB[0] / SVB[1]);
					b = abs((Grid[1][0] - SPoint[0]) - (Grid[1][1] - SPoint[1])*SVB[0] / SVB[1]);
					coef = (a + b)*y0*0.5;
				}
				break;
			}
			case 3://case [0,3]
			{
				x0 = Grid[3][0] - SPoint[0];
				y0 = Grid[3][1] - SPoint[1];
				if (abs(SVB[0] * SVB[1]) > 0)
				{
					t = x0*SVB[1] - y0*SVB[0];
					a = t / SVB[0];
					b = t / SVB[1];
					coef = 0.5*abs(a*b);
					coef = area - coef;
				}
				else
					coef = area;
				break;
			}
			case 4: // case [0,4]
			{
				coef = area;
				break;
			}
			default: break;
		}
		break;
	}//end case 0 of AI
	case 1:
	{
		switch (BI)
		{
			case 1://case [1,1]
			{
				x0 = Grid[0][0] - SPoint[0];
				y0 = Grid[0][1] - SPoint[1];
				t = x0*SVB[1] - y0*SVB[0];
				if (abs(SVB[0] * SVB[1]) > 0)
				{
					a = t / SVB[0];
					b = t / SVB[1];
					coef = 0.5*abs(a*b);
				}
				t = x0*SVA[1] - y0*SVA[0];
				if (abs(SVA[0] * SVA[1]) > 0)
				{
					a = t / SVA[0];
					b = t / SVA[1];
					coef = abs(coef - 0.5*abs(a*b));
				}
				break;
			}
			case 2://case [1,2]
			{
				x0 = abs(Grid[0][0] - Grid[1][0]);
				y0 = abs(Grid[0][1] - Grid[1][1]);
				if (x0 > y0) // line is on the x-dirction
				{
					a = abs((Grid[0][0] - SPoint[0])*SVB[1] / SVB[0] - (Grid[0][1] - SPoint[1]));
					b = abs((Grid[1][0] - SPoint[0])*SVB[1] / SVB[0] - (Grid[1][1] - SPoint[1]));
					coef = (a + b)*x0*0.5;
				}
				else
				{
					a = abs((Grid[0][0] - SPoint[0]) - (Grid[0][1] - SPoint[1])*SVB[0] / SVB[1]);
					b = abs((Grid[1][0] - SPoint[0]) - (Grid[1][1] - SPoint[1])*SVB[0] / SVB[1]);
					coef = (a + b)*y0*0.5;
				}
				x0 = Grid[0][0] - SPoint[0];
				y0 = Grid[0][1] - SPoint[1];
				if (abs(SVA[0] * SVA[1]) > 0)
				{
					t = x0*SVA[1] - y0*SVA[0];
					a = t / SVA[0];
					b = t / SVA[1];
					coef = abs(0.5*abs(a*b) - coef);
				}
				break;
			}
			case 3://case [1,3]
			{
				x0 = Grid[0][0] - SPoint[0];
				y0 = Grid[0][1] - SPoint[1];
				if (abs(SVA[0] * SVA[1]) > 0)
				{
					t = x0*SVA[1] - y0*SVA[0];
					a = t / SVA[0];
					b = t / SVA[1];
					coef = area - 0.5*abs(a*b);
				}
				else
					coef = area;
				x0 = Grid[3][0] - SPoint[0];
				y0 = Grid[3][1] - SPoint[1];
				if (abs(SVB[0] * SVB[1]) > 0)
				{
					t = x0*SVB[1] - y0*SVB[0];
					a = t / SVB[0];
					b = t / SVB[1];
					coef = coef - 0.5*abs(a*b);
				}
				break;
			}
			case 4://case [1,4]
			{
				x0 = Grid[0][0] - SPoint[0];
				y0 = Grid[0][1] - SPoint[1];
				if (abs(SVA[0] * SVA[1]) > 0)
				{
					t = x0*SVA[1] - y0*SVA[0];
					a = t / SVA[0];
					b = t / SVA[1];
					coef = 0.5*abs(a*b);
					coef = area - coef;
				}
				else
					coef = area;
				break;
			}
			default: break;
		}
		break;
	}//end case 1 of AI
	case 2:
	{
		switch (BI) {
		case 2:
		{
			x0 = abs(Grid[0][0] - Grid[1][0]);
			y0 = abs(Grid[0][1] - Grid[1][1]);
			if (x0 > y0) // line is on the x-dirction
			{
				a = abs((Grid[0][0] - SPoint[0])*SVB[1] / SVB[0] - (Grid[0][1] - SPoint[1]));
				b = abs((Grid[1][0] - SPoint[0])*SVB[1] / SVB[0] - (Grid[1][1] - SPoint[1]));
				coef = (a + b)*x0*0.5;
				a = abs((Grid[0][0] - SPoint[0])*SVA[1] / SVA[0] - (Grid[0][1] - SPoint[1]));
				b = abs((Grid[1][0] - SPoint[0])*SVA[1] / SVA[0] - (Grid[1][1] - SPoint[1]));
				coef = abs(coef - (a + b)*x0*0.5);
			}
			else
			{
				a = abs((Grid[0][0] - SPoint[0]) - (Grid[0][1] - SPoint[1])*SVB[0] / SVB[1]);
				b = abs((Grid[1][0] - SPoint[0]) - (Grid[1][1] - SPoint[1])*SVB[0] / SVB[1]);
				coef = (a + b)*y0*0.5;
				a = abs((Grid[0][0] - SPoint[0]) - (Grid[0][1] - SPoint[1])*SVA[0] / SVA[1]);
				b = abs((Grid[1][0] - SPoint[0]) - (Grid[1][1] - SPoint[1])*SVA[0] / SVA[1]);
				coef = abs(coef - (a + b)*y0*0.5);
			}
			break;
		}
		case 3:
		{
			x0 = abs(Grid[2][0] - Grid[3][0]);
			y0 = abs(Grid[2][1] - Grid[3][1]);
			if (x0 > y0) // line is on the x-dirction
			{
				a = abs((Grid[2][0] - SPoint[0])*SVA[1] / SVA[0] - (Grid[2][1] - SPoint[1]));
				b = abs((Grid[3][0] - SPoint[0])*SVA[1] / SVA[0] - (Grid[3][1] - SPoint[1]));
				coef = (a + b)*x0*0.5;
			}
			else
			{
				a = abs((Grid[2][0] - SPoint[0]) - (Grid[2][1] - SPoint[1])*SVA[0] / SVA[1]);
				b = abs((Grid[3][0] - SPoint[0]) - (Grid[3][1] - SPoint[1])*SVA[0] / SVA[1]);
				coef = (a + b)*y0*0.5;
			}
			x0 = Grid[3][0] - SPoint[0];
			y0 = Grid[3][1] - SPoint[1];
			if (abs(SVB[0] * SVB[1]) > 0)
			{
				t = x0*SVB[1] - y0*SVB[0];
				a = t / SVB[0];
				b = t / SVB[1];
				coef = abs(0.5*abs(a*b) - coef);
			}
			break;
		}
		case 4:
		{
			x0 = abs(Grid[2][0] - Grid[3][0]);
			y0 = abs(Grid[2][1] - Grid[3][1]);
			if (x0 > y0) // line is on the x-dirction
			{
				a = abs((Grid[2][0] - SPoint[0])*SVA[1] / SVA[0] - (Grid[2][1] - SPoint[1]));
				b = abs((Grid[3][0] - SPoint[0])*SVA[1] / SVA[0] - (Grid[3][1] - SPoint[1]));
				coef = (a + b)*x0*0.5;
			}
			else // line is on the y-direction
			{
				a = abs((Grid[2][0] - SPoint[0]) - (Grid[2][1] - SPoint[1])*SVA[0] / SVA[1]);
				b = abs((Grid[3][0] - SPoint[0]) - (Grid[3][1] - SPoint[1])*SVA[0] / SVA[1]);
				coef = (a + b)*y0*0.5;
			}
			break;
		}
		default: break;
		}
		break;
	}//end case 2 of AI
	case 3:
	{
		switch (BI) {
			case 3:
			{
				x0 = Grid[3][0] - SPoint[0];
				y0 = Grid[3][1] - SPoint[1];
				if (abs(SVB[0] * SVB[1]) > 0)
				{
					t = x0*SVB[1] - y0*SVB[0];
					a = t / SVB[0];
					b = t / SVB[1];
					coef = 0.5*abs(a*b);
				}
				if (abs(SVA[0] * SVA[1]) > 0)
				{
					t = x0*SVA[1] - y0*SVA[0];
					a = t / SVA[0];
					b = t / SVA[1];
					coef = abs(coef - 0.5*abs(a*b));
				}
				break;
			}
			case 4:
			{
				x0 = Grid[3][0] - SPoint[0];
				y0 = Grid[3][1] - SPoint[1];
				if (abs(SVA[0] * (SVA[1])) > 0)
				{
					t = x0*SVA[1] - y0*SVA[0];
					a = t / SVA[0];
					b = t / SVA[1];
					coef = 0.5*abs(a*b);
					//mexPrintf("x0=%f,y0=%f,SVA[0]=%f,SVA[1]=%f,t=%f,a=%f,b=%f;coef=%f\n",x0,y0,SVA[0],SVA[1],t,a,b,coef);
				}
				break;
			}
			default: break;
		}
		break;
	}//end case 3 of AI
	}//end of switch AI
	return coef;
}




// Generate projection matrix
// Generate The projection matrix in AIM model
template<typename T>
void genProj_AIM(
		std::vector<int>& rowIdx,
		std::vector<int>& colIdx,
		std::vector<T>& coeffs,
		int angIdx, const T ang, const FanEDGeo FanGeo, const Image Img)
{
	int pubIdx = 0;
	const T ScanR = FanGeo.m_S2O;

	const int XN = Img.m_Reso.x;
	const int YN = Img.m_Reso.y;
	const int DN = FanGeo.m_DetN;

	T* PosAry = new T[(XN + 1) * (YN + 1)];
	for (pubIdx = 0; pubIdx != (XN + 1) * (YN + 1); ++pubIdx)
	{
		PosAry[pubIdx] = 0;
	}
	const T dx = Img.m_Step.x;
	const T dy = Img.m_Step.y;
	const T area = dx * dy;
	const T DLen = FanGeo.m_DetSize;
	const T dd = DLen / DN;
	const T xctr = XN * 0.5;
	const T yctr = YN * 0.5;
	const T dctr = FanGeo.m_DetCntIdx;
	const T cosAng = cos(ang);
	const T sinAng = sin(ang);

	T Ew[2]; T Ed[2];
	Ew[0] = -cosAng;
	Ew[1] = -sinAng;
	Ed[0] = -sinAng;
	Ed[1] = cosAng;
	T SPoint[2];
	SPoint[0] = ScanR * cosAng;
	SPoint[1] = ScanR * sinAng;

	int xi(0), yi(0);
	T xcor(0), ycor(0), dcor(0);
	T Grid[4][3];

	for (yi = 0; yi <= YN; ++yi)
	{
		ycor = (yi - yctr) * dy - SPoint[1];
		for (xi = 0; xi <= XN; ++xi)
		{
			xcor = (xi - xctr) * dx - SPoint[0];
			dcor = ScanR * (xcor* Ed[0] + ycor * Ed[1]) / (xcor * Ew[0] + ycor * Ew[1]);
			dcor = dcor / dd;
			PosAry[yi* (XN + 1) + xi] = dcor + dctr;
		}
	}

	int posim(0);
	T pvalue(0);
	T pcfun(0);
	T pdist(0);
	int MinBV(0);
	int MaxBV(0);
	T temp(0);
	T SVA[3], SVB[3];
	T pangle(0);
	int di(0);
	T coef(0);
	for (yi = 0; yi < YN; ++yi) {
		for (xi = 0; xi < XN; ++xi)	{
			//Fetch the four points of the pixel and their projection positions
			Grid[0][0] = (xi - xctr)*dx;
			Grid[0][1] = (yi - yctr)*dy;
			Grid[0][2] = PosAry[yi*(XN + 1) + xi];

			Grid[1][0] = (xi - xctr + 1)*dx;
			Grid[1][1] = (yi - yctr)*dy;
			Grid[1][2] = PosAry[yi*(XN + 1) + xi + 1];

			Grid[2][0] = (xi - xctr + 1)*dx;
			Grid[2][1] = (yi - yctr + 1)*dy;
			Grid[2][2] = PosAry[(yi + 1)*(XN + 1) + xi + 1];

			Grid[3][0] = (xi - xctr)*dx;
			Grid[3][1] = (yi - yctr + 1)*dy;
			Grid[3][2] = PosAry[(yi + 1)*(XN + 1) + xi];

			SortProjection<T>(Grid);//Sort the projection psotion

			posim = yi*XN + xi;
			pvalue = 0.0;
			pdist = hypot((xi + 0.5 - xctr)*dx - SPoint[0], (yi + 0.5 - yctr) * dy - SPoint[1]);

			MinBV = int(Grid[0][2] + 10) - 10;
			MaxBV = int(Grid[3][2] + 10) - 9;
			if (MinBV < 0)   MinBV = 0;
			if (MaxBV > DN)  MaxBV = DN;
			// Compute the directions of the two lines for the projections
			for (di = MinBV; di < MaxBV; di++)
			{
				temp = (di - dctr) * dd;
				SVA[0] = -temp * sinAng - SPoint[0];
				SVA[1] = temp * cosAng - SPoint[1];
				SVA[2] = di;
				temp = hypot(SVA[0], SVA[1]); // sqrt(pow(SVA[0], 2) + pow(SVA[1], 2));
				SVA[0] = SVA[0] / temp;
				SVA[1] = SVA[1] / temp;

				temp = (di - dctr + 1)*dd;
				SVB[0] = -temp * sinAng - SPoint[0];
				SVB[1] = temp * cosAng - SPoint[1];
				SVB[2] = di + 1;
				temp = hypot(SVB[0], SVB[1]); // sqrt(pow(SVB[0], 2) + pow(SVB[1], 2));
				SVB[0] = SVB[0] / temp;
				SVB[1] = SVB[1] / temp;

				pangle = SVA[0] * SVB[0] + SVA[1] * SVB[1];
				pangle = sqrt(1 - pow(pangle, 2));
				//compute the weighting coefficient for a special projection data
				coef = ComputeCoefficient<T>(Grid, SVA, SVB, SPoint, area);
				coef = coef / (pdist*pangle);
				//Main code to forward projection and back projection
				rowIdx.push_back(angIdx * DN + di);
				colIdx.push_back(posim);
				coeffs.push_back(coef);
			}
		}
	}
	delete [] PosAry;
}

template<typename T>
void genProj_AIM(std::vector<int>& rowIdx, std::vector<int>& colIdx, std::vector<T>& coeffs, int angIdx, const T ang, const FanEAGeo FanGeo, const Image Img)
{
	const T ScanR = FanGeo.m_S2O;
	const T ObjR = Img.m_Size.x;

	const int DN = FanGeo.m_DetN;
	const int XN = Img.m_Reso.x;
	const int YN = Img.m_Reso.y;
	const T dx = Img.m_Step.x;
	const T dy = Img.m_Step.y;
	const T area = dx * dy;
	const T xctr = XN * 0.5;
	const T yctr = YN * 0.5;
	const T dctr = FanGeo.m_DetCntIdx;
	const T DtBeta = FanGeo.m_DetStp;
	T coef;
	T* PosAry = new T[(XN + 1) * (YN + 1)];
	T* PBeta = new T[DN];
	for (int i = 0; i != DN; ++i)
	{
		PBeta[i] = (i - dctr + 0.5) * DtBeta;
	}
	T Ew[2], Ed[2];
	const T cosAng = cos(ang);
	const T sinAng = sin(ang);
	Ew[0] = -cosAng;
	Ew[1] = -sinAng;
	Ed[0] = -sinAng;
	Ed[1] = cosAng;
	T SPoint[2];
	SPoint[0] = ScanR * cosAng;
	SPoint[1] = ScanR * sinAng;
	int xi(0), yi(0);

	T xcor, ycor, dcor;
	T pdist;
	int di;
	for (yi = 0; yi <= YN; ++yi)
	{
		ycor = (yi - yctr)*dy - SPoint[1];
		for (xi = 0; xi <= XN; ++xi)
		{
			xcor = (xi - xctr)*dx - SPoint[0];
			dcor = (xcor*Ed[0] + ycor*Ed[1]) / (xcor*Ew[0] + ycor*Ew[1]);
			dcor = (atan(dcor) - PBeta[0]) / DtBeta;
			PosAry[yi*(XN + 1) + xi] = dcor + 0.5;
		}
	}

	T Grid[4][3];
	int posim(0);
	int MinBV(0);
	int MaxBV(0);
	T pangle, temp, SVA[3], SVB[3];
	for (yi = 0; yi < YN; yi++)
	{
		for (xi = 0; xi < XN; xi++)
		{
			//Fetch the four points of the pixel and their projection positions
			Grid[0][0] = (xi - xctr)*dx;
			Grid[0][1] = (yi - yctr)*dy;
			Grid[0][2] = PosAry[yi*(XN + 1) + xi];
			Grid[1][0] = (xi - xctr + 1)*dx;
			Grid[1][1] = (yi - yctr)*dy;
			Grid[1][2] = PosAry[yi*(XN + 1) + xi + 1];
			Grid[2][0] = (xi - xctr + 1)*dx;
			Grid[2][1] = (yi - yctr + 1)*dy;
			Grid[2][2] = PosAry[(yi + 1)*(XN + 1) + xi + 1];
			Grid[3][0] = (xi - xctr)*dx;
			Grid[3][1] = (yi - yctr + 1)*dy;
			Grid[3][2] = PosAry[(yi + 1)*(XN + 1) + xi];
			SortProjection<T>(Grid);//Sort the projection psotion

			posim = yi*XN + xi;

			//pvalue = Img[posim];
			pdist = hypot((xi + 0.5 - xctr)*dx - SPoint[0], (yi + 0.5 - yctr)*dy - SPoint[1]); // sqrt(pow(, 2) + pow(, 2));

			//Computer the weighting coefficient for every projection position
			MinBV = int(Grid[0][2] + 10) - 10;
			MaxBV = int(Grid[3][2] + 10) - 9;
			if (MinBV < 0)   MinBV = 0;
			if (MaxBV > DN)  MaxBV = DN;
			//double total =0;
			for (di = MinBV; di < MaxBV; di++)
			{
				// Compute the directions of the two lines for the projections
				pangle = PBeta[di] - 0.5 * DtBeta;
				temp = ang + PI - pangle;
				SVA[0] = cos(temp);
				SVA[1] = sin(temp);
				SVA[2] = di;
				// mexPrintf("di=%d,VA0=%10.8f,VA1=%10.8f,angle=%10.8f,Beta=%10.8f\n",di,SVA[0],SVA[1],temp,pangle);
				pangle = PBeta[di] + 0.5 * DtBeta;
				temp = ang + PI - pangle;
				SVB[0] = cos(temp);
				SVB[1] = sin(temp);
				SVB[2] = di + 1;

				//compute the weighting coefficient for a special projection data
				coef = ComputeCoefficient(Grid, SVA, SVB, SPoint, area);
				coef = coef / (pdist*abs(DtBeta));
				rowIdx.push_back(angIdx * DN + di);
				colIdx.push_back(posim);
				coeffs.push_back(coef);

			}
		}
	}

	delete [] PosAry;
	delete [] PBeta;
}


template<typename T, typename FanGeometry>
void genProjectionMatrix_AIM_template_Impl(const FanGeometry FanGeo, const Image Img)
{
	unsigned int angIdx = 0;
	std::vector < int> rowIdx;
	std::vector < int > colIdx;
	std::vector < T > coeffs;
	T ang(0);
	for (angIdx = 0; angIdx < FanGeo.m_ViwN; ++angIdx)
	{
		ang = FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp;
		genProj_AIM(rowIdx, colIdx, coeffs, angIdx, ang, FanGeo, Img);
		std::cout << angIdx << std::endl;
	}
	int nonZ = rowIdx.size();
	std::stringstream ss;
	ss << nonZ;
	std::string f1;
	std::string f2;
	std::string f3;
	if (sizeof(T) == 4)
	{
		f1 = "prjAIM" + ss.str() + "f.row";
		f2 = "prjAIM" + ss.str() + "f.col";
		f3 = "prjAIM" + ss.str() + "f.cof";
	}
	else if (sizeof(T) == 8)
	{
		f1 = "prjAIM" + ss.str() + "d.row";
		f2 = "prjAIM" + ss.str() + "d.col";
		f3 = "prjAIM" + ss.str() + "d.cof";
	}

	std::ofstream rowFile(f1.c_str(), std::ios::binary);
	std::ofstream colFile(f2.c_str(), std::ios::binary);
	std::ofstream coeFile(f3.c_str(), std::ios::binary);
	rowFile.write((char*)&(rowIdx[0]), sizeof(int) * rowIdx.size());
	rowFile.close();
	colFile.write((char*)&(colIdx[0]), sizeof(int) * colIdx.size());
	colFile.close();
	coeFile.write((char*)&(coeffs[0]), sizeof(T) * coeffs.size());
	coeFile.close();

	rowFile.close();
	colFile.close();
	coeFile.close();

}

template<typename T>
void genProjectionMatrix_AIM_template(const FanEDGeo FanGeo, const Image Img) {
	genProjectionMatrix_AIM_template_Impl<T, FanEDGeo>(FanGeo, Img);
}

template<typename T>
void genProjectionMatrix_AIM_template(const FanEAGeo FanGeo, const Image Img) {
	genProjectionMatrix_AIM_template_Impl<T, FanEAGeo>(FanGeo, Img);
}


template<typename T>
static void DD3Boundaries(int nrBoundaries, T*pCenters, T *pBoundaries)
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
std::vector<T> getDD3Boundaries(int nCenters, T* pCenters) {
	std::vector<T> res(nCenters + 1, 0);
	T* pBoundaries = &(res[0]);
	int nrBoundaries = nCenters + 1;
	DD3Boundaries(nrBoundaries, pCenters, pBoundaries);
	return res;
}

template<typename T>
thrust::device_vector<T> getDD3Boundaries(T* pCenters, int nCenters) {
	thrust::device_vector<T> res{ getDD3Boundaries(nCenters, pCenters) };
	return res;
}

template<typename T>
std::vector<T> getDD3Boundaries(const std::vector<T>& centers) {
	typename std::vector<T>::iterator pCenters = centers.begin();
	std::vector<T> boundaries(centers.size() + 1, 0);
	typename std::vector<T>::iterator pBoundaries = boundaries.begin();
	int i;
	if (boundaries.size() >= 3) {
		*pBoundaries++ = 1.5 * *pCenters - 0.5 * *(pCenters + 1);
		for (i = 1; i <= (centers.size() - 1); i++) {
			*pBoundaries++ = 0.5 * *pCenters + 0.5 * *(pCenters + 1);
			pCenters++;
		}
		*pBoundaries = 1.5 * *pCenters - 0.5 * *(pCenters - 1);
	}
	else {
		*pBoundaries = *pCenters - 0.5;
		*(pBoundaries + 1) = *pCenters + 0.5;
	}
	return boundaries;
}

template<typename T>
std::vector<T>& operator+(std::vector<T>& v, const T& ms)
{
	std::transform(v.begin(), v.end(), v.begin(), [=](const T& vv){ return vv + ms; });
	return v;
}
template<typename T>
std::vector<T>& operator-(std::vector<T>& v, const T& ms)
{
	std::transform(v.begin(), v.end(), v.begin(), [=](const T& vv){ return vv - ms; });
	return v;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
// Get one sub-volume from the whole volume.
// Assume that the volumes are stored in Z, X, Y order
template<typename T>
void getSubVolume(const T* vol,
	const size_t XYN, const size_t ZN,
	const size_t ZIdx_Start, const size_t ZIdx_End, T* subVol)
{
	const size_t SZN = ZIdx_End - ZIdx_Start;
	for (size_t xyIdx = 0; xyIdx != XYN; ++xyIdx) {
		for (size_t zIdx = ZIdx_Start; zIdx != ZIdx_End; ++zIdx) {
			subVol[xyIdx * SZN + (zIdx - ZIdx_Start)] = vol[xyIdx * ZN + zIdx];
		}
	}
}
template<typename T>
void getSubVolume(const T* vol,
	const size_t XN, const size_t YN, const size_t ZN,
	const size_t ZIdx_Start, const size_t ZIdx_End, T* subVol) {
	getSubVolume(vol, XN * YN, ZN, ZIdx_Start, ZIdx_End);
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