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
 * \date Sep. 1, 2015
 *
 */
#pragma once
#include <algorithm>
#include <cassert>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cusparse.h>
#include <cusparse_v2.h>
#include "cublas_v2.h"

#include <cuda_runtime.h>
#include "device_launch_parameters.h"
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

#include <vector>
#include <vector_functions.h>
#include <vector_types.h>


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
#define checkCudaErrors(value) {value;}
#endif

typedef unsigned int uint;
typedef const unsigned int cuint;
typedef unsigned char byte;
typedef thrust::device_vector<float> d_vec_t;
typedef thrust::host_vector<float> h_vec_t;


/// \brief Define that macro for Thrust compiling successful
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS 1
#endif


#ifndef PI
#define PI (3.14159265358979323846264) 
#endif

#ifndef _PI 
#define _PI (3.14159265358979323846264)
#endif

#ifndef _PI_2
#define _PI_2 (_PI * 0.5)
#endif

#ifndef _TWOPI
#define _TWOPI (_PI * 2.0)
#endif

#ifndef _PI_4
#define _PI_4 (_PI * 0.25)
#endif

#ifndef _7PI_4
#define _7PI_4 (_PI * 1.75)
#endif

#ifndef _3PI_4
#define _3PI_4 (_PI * 0.75)
#endif


#ifndef _5PI_4
#define _5PI_4 (_PI * 1.25)
#endif

#ifndef _EPSILON
#define _EPSILON (1.0E-9)
#endif


//Calculate Two Arrays inner product
///This function is not used for Users but for calling by system; 
///IT IS A META-PROGRAM
///Author: Rui LIU
///Date: 2013-2-18
///Version: 1.0
template<int Dims,typename T>
class DotProduct
{
public:
	static T result(T* a,T* b)
	{
		return (*a)*(*b) + DotProduct<Dims-1,T>::result(a+1,b+1);
	}
};
template<typename T>
class DotProduct<1,T>
{
public:
	static T result(T* a,T* b)
	{
		return (*a)*(*b);
	}
};

//////////////////////////////////////////////////////////////////////////
/// Call the dot production meta function
template<int Dims,typename T>
inline T innerProd_META(T* a, T* b)
{
	return DotProduct<Dims,T>::result(a,b);
}



//Add two vector together
//This function is not used for Users but for calling by system;
//It is a META-program
//Author: Rui LIU
//Date: 2013-8-9
//Version: 1.0
template<int Dims, typename T>
class AddVector
{
public:
	static void result(T* a, T* b, T* c)
	{
		(*c) = (*a) + (*b);
		AddVector<Dims-1,T>::result(a+1,b+1);
	}
};
template<typename T>
class AddVector<1,T>
{
public:
	static T result(T* a,T* b,T* c)
	{
		(*c) = (*a) + (*b);
	}
};

template<int Dims,typename T>
inline void add_META(T* a, T* b, T* c)
{
	return AddVector<Dims,T>::result(a,b,c);
}


//Sub vector a from b together
//This function is not used for Users but for calling by system;
//It is a META-program
//Author: Rui LIU
//Date: 2013-8-9
//Version: 1.0
template<int Dims, typename T>
class SubVector
{
public:
	static void result(T* a, T* b, T* c)
	{
		(*c) = (*a) - (*b);
		SubVector<Dims-1,T>::result(a+1,b+1);
	}
};
template<typename T>
class SubVector<1,T>
{
public:
	static T result(T* a,T* b, T* c)
	{
		(*c) = (*a) - (*b);
	}
};

template<int Dims,typename T>
inline void sub_META(T* a, T* b, T* c)
{
	return SubVector<Dims,T>::result(a,b,c);
}



///This metaprogram is used to calculate the p power for each element and then add them together
//It is a META-PROGRAM
//Author: Rui Liu
//Date: 2013-08-13
//Version: 1.0
template<int Dims, typename T>
class powerVector
{
public:
	static T result(T* vec, T p)
	{
		return std::pow((*vec),p) + powerVector<Dims-1,T>::result(vec+1,p);
	}
};
template<typename T>
class powerVector<1,T>
{
public:
	static T result(T *vec, T p)
	{
		return std::pow((*vec),p);
	}
};



///This meta function solute the problem that a matrix multiply with
/// a vector
///Author: Rui Liu
/// Use the meta programming method to generate the 
/// cosine and sine series given a series of theta
template<int Dims, typename T>
class CosSinTheta
{
public:
	static void result(T* cosTheta, T* sinTheta, T* theta)
	{
		(*cosTheta) = cos((*theta));
		(*sinTheta) = sin((*theta));
		CosSinTheta<Dims-1,T>::result(cosTheta+1,sinTheta+1,theta+1);
	}
};
template<typename T>
class CosSinTheta<1,T>
{
public:
	static T result(T* cosTheta,T* sinTheta,T* theta)
	{
		(*cosTheta) = cos((*theta));
		(*sinTheta) = sin((*theta));
	}
};


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
inline __host__ __device__ bool IS_ZERO(const T& x)
{
	return (x < 1.0E-9) && (x > -1.0E-9);
}


/* Define the Data Type of the files in enum type*/
typedef enum{
	RUI_UCHAR = 0X00,
	RUI_CHAR = 0X02,
	RUI_UINT = 0X04,
	RUI_INT = 0X08,
	RUI_FLOAT = 0X10,
	RUI_DOUBLE = 0X20,
	RUI_LONGDOUBLE = 0X40,
	RUI_OTHER = 0X80,
}RUI_DATA_TYPE;

template<typename T>
INLINE __host__ __device__ T safeDivide(const T& a, const T& b)
{
	if(IS_ZERO(b))
	{
		return 0;
	}
	else
	{
		return a/b;
	}
}

INLINE __host__ __device__ float2 operator+(const float2& a, const float2& b){return make_float2(a.x + b.x, a.y + b.y);};
INLINE __host__ __device__ float2 operator-(const float2& a, const float2& b){return make_float2(a.x - b.x, a.y - b.y);};
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
INLINE __host__ __device__ double3 operator+(const double3& a, const double3& b){ return make_double3(a.x + b.x, a.y + b.y, a.z + b.z); }
INLINE __host__ __device__ double3 operator+(const double3& a, const float3& b){ return make_double3(a.x + b.x, a.y + b.y, a.z + b.z); }
INLINE __host__ __device__ double3 operator+(const float3& a, const double3& b){ return make_double3(a.x + b.x, a.y + b.y, a.z + b.z); }
INLINE __host__ __device__ double3 operator-(const double3& a, const double3& b){	return make_double3(a.x - b.x, a.y - b.y, a.z - b.z);}
INLINE __host__ __device__ double3 operator-(const double3& a, const float3& b){ return make_double3(a.x - b.x, a.y - b.y, a.z - b.z); }
INLINE __host__ __device__ double3 operator-(const float3& a, const double3& b){ return make_double3(a.x - b.x, a.y - b.y, a.z - b.z); }
INLINE __host__ __device__ double3 operator*(const double3& a, const double3& b){ return make_double3(a.x * b.x, a.y * b.y, a.z * b.z); }
INLINE __host__ __device__ double3 operator*(const double3& a, const float3& b){ return make_double3(a.x * b.x, a.y * b.y, a.z * b.z); }
INLINE __host__ __device__ double3 operator*(const float3& a, const double3& b){ return make_double3(a.x * b.x, a.y * b.y, a.z * b.z); }
INLINE __host__ __device__ double3 operator/(const double3& a, const double3& b){return make_double3(safeDivide<double>(a.x,b.x),safeDivide<double>(a.y,b.y),safeDivide<double>(a.z,b.z));}
INLINE __host__ __device__ double3 operator/(const double3& a, const float3& b){return make_double3(safeDivide<double>(a.x,b.x),safeDivide<double>(a.y,b.y),safeDivide<double>(a.z,b.z));}
INLINE __host__ __device__ double3 operator/(const float3& a, const double3& b){return make_double3(safeDivide<double>(a.x,b.x),safeDivide<double>(a.y,b.y),safeDivide<double>(a.z,b.z));}
INLINE __host__ __device__ double3 operator/(const double3& a, const double& b){return make_double3(safeDivide<double>(a.x,b),safeDivide<double>(a.y,b),safeDivide<double>(a.z,b));};


INLINE __host__ __device__ float fminf(const float2& a){	return fminf(a.x, a.y);}
INLINE __host__ __device__ float fminf(const float3& a){	return fminf(a.x, fminf(a.y, a.z)); }
INLINE __host__ __device__ float2 fminf(const float2& a, const float2& b){ return make_float2(fmin(a.x, b.x), fmin(a.y, b.y)); }
INLINE __host__ __device__ float3 fminf(const float3& a, const float3& b){	return make_float3(fminf(a.x, b.x), fminf(a.y, b.y), fminf(a.z, b.z));}
INLINE __host__ __device__ double2 fminf(const double2& a, const double2& b){ return make_double2(fmin(a.x, b.x), fmin(a.y, b.y)); }
INLINE __host__ __device__ double2 fminf(const float2& a, const double2& b){ return make_double2(fmin(a.x, b.x), fmin(a.y, b.y)); }
INLINE __host__ __device__ double2 fminf(const double2& a, const float2& b){ return make_double2(fmin(a.x, b.x), fmin(a.y, b.y)); }
INLINE __host__ __device__ double3 fminf(const double3& a, const double3& b){ return make_double3(fmin(a.x, b.x), fmin(a.y, b.y), fmin(a.z, b.z)); }
INLINE __host__ __device__ double3 fminf(const double3& a, const float3& b){ return make_double3(fmin(a.x, b.x), fmin(a.y, b.y), fmin(a.z, b.z)); }
INLINE __host__ __device__ double3 fminf(const float3& a, const double3& b){ return make_double3(fmin(a.x, b.x), fmin(a.y, b.y), fmin(a.z, b.z)); }
INLINE __host__ __device__ float fmaxf(const float2& a){	return fmaxf(a.x, a.y);}
INLINE __host__ __device__ float fmaxf(const float3& a){	return fmaxf(a.x, fmaxf(a.y, a.z));}
INLINE __host__ __device__ float2 fmaxf(const float2& a, const float2& b){ return make_float2(fmax(a.x, b.x), fmax(a.y, b.y)); }
INLINE __host__ __device__ float3 fmaxf(const float3& a, const float3& b){	return make_float3(fmaxf(a.x, b.x), fmaxf(a.y, b.y), fmaxf(a.z, b.z));}
INLINE __host__ __device__ double2 fmaxf(const double2& a, const double2& b){ return make_double2(fmax(a.x, b.x), fmax(a.y, b.y)); }
INLINE __host__ __device__ double2 fmaxf(const float2& a, const double2& b){ return make_double2(fmax(a.x, b.x), fmax(a.y, b.y)); }
INLINE __host__ __device__ double2 fmaxf(const double2& a, const float2& b){ return make_double2(fmax(a.x, b.x), fmax(a.y, b.y)); }
INLINE __host__ __device__ double3 fmaxf(const double3& a, const double3& b){ return make_double3(fmax(a.x, b.x), fmax(a.y, b.y), fmax(a.z, b.z)); }
INLINE __host__ __device__ double3 fmaxf(const double3& a, const float3& b){ return make_double3(fmax(a.x, b.x), fmax(a.y, b.y), fmax(a.z, b.z)); }
INLINE __host__ __device__ double3 fmaxf(const float3& a, const double3& b){ return make_double3(fmax(a.x, b.x), fmax(a.y, b.y), fmax(a.z, b.z)); }

INLINE __host__ __device__ float length(const float2& a){	return sqrtf(a.x * a.x + a.y * a.y);}
INLINE __host__ __device__ float length(const float3& a){	return sqrtf(a.x * a.x + a.y * a.y + a.z * a.z);}
INLINE __host__ __device__ double length(const double3& a){	return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);}
INLINE __host__ __device__ const float2 normalize(const float2& a){	return a / length(a);}
INLINE __host__ __device__ const float3 normalize(const float3& a){	return a / length(a);}
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

///// \brief Fast bilinear interpolation in double format
//INLINE __device__ double bilerp(int2 v0, int2 v1, int2 v2, int2 v3, float t1, float t2)
//{
//	double v0_ = __hiloint2double(v0.y, v0.x);
//	double v1_ = __hiloint2double(v1.y, v1.x);
//	double v2_ = __hiloint2double(v2.y, v2.x);
//	double v3_ = __hiloint2double(v3.y, v3.x);
//
//	double vv0 = v0_ * (1.0 - t1) + v1_ * t1;
//	double vv1 = v2_ * (1.0 - t1) + v3_ * t1;
//	return vv0 * (1 - t2) + vv1 * t2;
//}
//
///// \brief Fast bilinear interpolation in double format
//INLINE __device__ double bilerp(int2 v0, int2 v1, int2 v2, int2 v3, double t1, double t2)
//{
//	double v0_ = __hiloint2double(v0.y, v0.x);
//	double v1_ = __hiloint2double(v1.y, v1.x);
//	double v2_ = __hiloint2double(v2.y, v2.x);
//	double v3_ = __hiloint2double(v3.y, v3.x);
//
//	double vv0 = v0_ * (1.0 - t1) + v1_ * t1;
//	double vv1 = v2_ * (1.0 - t1) + v3_ * t1;
//	return vv0 * (1 - t2) + vv1 * t2;
//}



INLINE __host__ __device__ bool intersectBox(
	const float3& sour,
	const float3& dir,
	const float3& boxmin,
	const float3& boxmax,
	float* tnear, float* tfar)
{
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



template<typename T>
INLINE __host__ __device__ T regularizeAngle(T curang)
{
	T c = curang;
	while (c >= _TWOPI){ c -= _TWOPI; }
	while (c < 0){ c += _TWOPI; }
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


/// (x,y) template
template<typename T>
class T2
{
public:
	T x;
	T y;
public:
	T2() :x(0), y(0){}
	T2(const T& t) :x(t), y(t){}
	T2(const T& xx, const T& yy) :x(xx), y(yy){}
	T2(const float2& o) :x(o.x), y(o.y){}
	T2(const double2& o) :x(o.x), y(o.y){}
	T2& operator=(const T2& o) const { x = o.x; y = o.y; return *this; }
	T2& operator=(const float2& o) const { x = o.x; y = o.y; return *this; }
	T2& operator=(const double2& o) const { x = o.x; y = o.y; return *this; }
	friend std::ostream& operator<<(std::ostream& os, const T2& t)
	{
		os << "(" << t.x << "," << t.y << ")" << std::endl;
		return os;
	}
};


/// (x,y,z) template
template<typename T>
class T3
{
public:
	T x;
	T y;
	T z;
public:
	T3() :x(0), y(0), z(0){}
	T3(const T& t) :x(t), y(t), z(t){}
	T3(const T& xx, const T& yy, const T& zz) :x(xx), y(yy), z(zz){}
	T3(const float3& o) :x(o.x), y(o.y), z(o.z){}
	T3(const double3& o) :x(o.x), y(o.y), z(o.z){}
	T3& operator=(const T3& o) const { x = o.x; y = o.y; z = o.z; return *this; }
	T3& operator=(const float3& o) const { x = o.x; y = o.y; z = o.z; return *this; }
	T3& operator=(const double3& o) const { x = o.x; y = o.y; z = o.z; return *this; }
	friend std::ostream& operator<<(std::ostream& os, const T3& t)
	{
		os << "(" << t.x << "," << t.y << "," << t.z << ")" << std::endl;
		return os;
	}
};


/// (x,y,z,w) template
template<typename T>
class T4
{
public:
	T x;
	T y;
	T z;
	T w;
public:
	T4() :x(0), y(0), z(0), w(0){}
	T4(const T& t) :x(t), y(t), z(t), w(t){}
	T4(const T& xx, const T& yy, const T& zz, const T& ww) :x(xx), y(yy), z(zz), w(ww){}
	T4(const float4& o) :x(o.x), y(o.y), z(o.z), w(o.z){}
	T4(const double4& o) :x(o.x), y(o.y), z(o.z), w(o.w){}
	T4(const float3& o) :x(o.x), y(o.y), z(o.z), w(1.0){}
	T4(const double3& o) :x(o.x), y(o.y), z(o.z), w(1.0){}
	T4& operator=(const T4& o) const { x = o.x; y = o.y; z = o.z; w = o.w; return *this; }
	T4& operator=(const float4& o) const { x = o.x; y = o.y; z = o.z; w = o.w; return *this; }
	T4& operator=(const double4& o) const { x = o.x; y = o.y; z = o.z; w = o.w; return *this; }
	T4& operator=(const float3& o) const { x = o.x; y = o.y; z = o.z; w = 1.0; return *this; }
	T4& operator=(const double3& o) const { x = o.x; y = o.y; z = o.z; w = 1.0; return *this; }
	friend std::ostream& operator<<(std::ostream& os, const T4& t)
	{
		os << "(" << t.x << "," << t.y << "," << t.z << "," << t.w << ")" << std::endl;
		return os;
	}
};





/// \brief Fan Beam Equal Angle Detector based CT system
class FanEAGeo
{
public:
	float m_S2O; ///< source to object distance
	float m_O2D; ///< object to detector distance
	float m_S2D; ///< source to detector distance
	uint m_ViwN; ///< view number
	float m_ViwBeg; ///< Begin viewer number
	float m_ViwStp; ///< View step size

	float m_DetArc; ///< Detector Arc angle
	uint m_DetN; ///< Detector cells number in one row
	float m_DetStp; ///< Detector cells size
	float m_DetCntIdx; ///< The index of the center of the detector
public:
	FanEAGeo(void);
	~FanEAGeo(void){};
	FanEAGeo(const FanEAGeo& rhs);
	/// \brief constructor
	/// \param S2O source to object distance
	/// \param O2D object to detector distance
	/// \param ViwN view number
	/// \param ViwBeg the begin view
	/// \param ViwEnd the End view
	/// \param DetArc the detector Arc
	/// \param DetN the number of detector cells on one row
	FanEAGeo(const float S2O, const float O2D, const unsigned int ViwN,
		const float ViwBeg, const float ViwEnd, const float DetArc,
		const unsigned int DetN);
};

/// \brief Fan Beam configuration with equal distance detector
class FanEDGeo
{
public:
	float m_S2O; ///< Source to Object distance
	float m_O2D; ///< Object to Detector distance
	float m_S2D; ///< Source to Detector distance

	uint m_ViwN; ///< view number
	float m_ViwBeg;///< View begin
	float m_ViwStp; ///< View Step

	float m_DetSize; ///< Detector size;
	uint m_DetN; ///< How many detector cells
	float m_DetStp; ///< detector cells size;
	float m_DetCntIdx; ///< detector center index;
public:
	/// \brief constructor
	FanEDGeo(void);
	/// \brief destructor
	~FanEDGeo(){};
	/// \brief copy constructor
	FanEDGeo(const FanEDGeo& rhs);
	/// \brief Constructor
	/// \param S2O source to object distance
	/// \param O2D object to detector distance
	/// \param ViwN View number
	/// \param ViwBeg The begin of the view
	/// \param ViwEnd The end of the view
	/// \param DetSize Detector size
	/// \param DetN detector cells number on one row
	FanEDGeo(const float S2O, const float O2D, const unsigned int ViwN,
		const float ViwBeg, const float ViwEnd, const float DetSize,
		const unsigned int DetN);

};


/// \brief Cone Beam configuration with equal angle detector
class ConeEAGeo{
public:
	float m_S2O; ///< Source to object distance
	float m_O2D; ///< object to detector distance
	float m_S2D;///< source to detector distance

	uint m_ViwN; ///< View number
	float m_ViwBeg;///< Begin of the view
	float m_ViwStp; ///< View step

	float m_DetArc; ///< Detector Arc size
	uint m_DetN; ///< detector cell number
	float m_DetStp; ///< detector cell size

	float m_DetHei;  ///< detector height
	uint m_DetHN; ///< Height cell number of the detector
	float m_DetHStp;///< detector step size on the height direction

	float2 m_DetCntIdx; ///< detector center index
public:
	/// \brief constructor
	ConeEAGeo(void);
	/// \brief destructor
	~ConeEAGeo(){}
	/// \brief copy constructor
	ConeEAGeo(const ConeEAGeo& rhs);
	/// \brief Constructor
	/// \param S2O source to object distance
	/// \param O2D object to detector distance
	/// \param viwN view number
	/// \param ViwBeg Begin of the view
	/// \param ViwEnd End of the view
	/// \param DetArc Arc size of the detector
	/// \param DetN detector cells on one row
	/// \param DetHei detector height size of the detector
	/// \param DetHN detector cells number on the height direction
	ConeEAGeo(const float S2O, const float O2D, const unsigned int viwN,
		const float ViwBeg, const float ViwEnd, const float DetArc, const unsigned int DetN,
		const float DetHei, const unsigned int DetHN);
};




/// \brief Cone Beam geometry with equal distance detector
class ConeEDGeo
{
public:
	float m_S2O; ///< Source to object distance
	float m_O2D; ///< object to detector distance
	float m_S2D;///< source to detector distance

	int m_ViwN; ///< view number
	float m_ViwBeg; ///< begin view number
	float m_ViwStp; ///< view step size

	float2 m_DetSize; ///< detector size on length and height direction
	int2 m_DetN; ///< detector Number on length and height direction
	float2 m_DetStp; ///< detector step size on length and height direction
	float2 m_DetCntIdx; ///< detector center index
public:
	/// \brief constructor
	ConeEDGeo(void);
	/// \brief destructor
	~ConeEDGeo(){}
	/// \brief copy constructor
	ConeEDGeo(const ConeEDGeo& rhs);
	/// \brief Constructor
	/// \param S2O source to object distance
	/// \param O2D object to detector distance
	/// \param ViwN view number
	/// \param ViwBeg Begin of the view
	/// \param ViwEnd End of the view
	/// \param DetSize Detector size on both direction
	/// \param DetN detector cells on both direction
	ConeEDGeo(const float S2O, const float O2D, const  int ViwN,
		const float ViwBeg, const float ViwEnd, const float2 DetSize,
		const int2 DetN);
};



/// \brief Image configuration class
class Image
{
public:
	uint2 m_Reso; ///< Image resolution
	float2 m_Size;///< Image size
	float2 m_Step; ///< Image Step
	float2 m_Bias; ///< The bias of the image
public:
	/// \brief constructor
	Image(void);
	/// \brief destructor
	~Image(void){};
	/// \brief copy constructor
	Image(const Image& rhs);
	/// \brief constructor
	Image(
		const unsigned int resoL,///< resolution on length direction
		const unsigned int resoW,///< resolution on width direction
		const float sizeL, ///< length size of the image
		const float sizeW,///< width size of the image
		const float BiasL, ///< bias on length direction
		const float BiasW ///<bias on width direction
		);
};

/// \brief Volume configuration class
class Volume
{
public:
	uint3 m_Reso; ///< Volume resolution
	float3 m_Size; ///< Image size
	float3 m_Step; ///< Image step size
	float3 m_Bias; ///< Bias of the image
public:
	/// \brief constructor
	Volume(void);
	/// \brief destructor
	~Volume(void){}
	/// \brief copy constructor
	Volume(const Volume& rhs);
	/// \brief
	Volume(
		const unsigned int resoL, ///< object length resolution
		const unsigned int resoW,///<object width resolution
		const unsigned int resoH, ///< object height resolution
		const float sizeL, ///< object size on length direction
		const float sizeW,///< object size on width direction
		const float sizeH,///< object size on height direction
		const float biasL,///< object bias on length direction
		const float biasW,///< object bias on width direction
		const float biasH///< object bias on height direction
		);
};



/// \brief Calculate the scalar transform of a vector, called by the thrust::transform or std::transform
template<typename T>
struct _scale_functor
{
	T _s;
	_scale_functor(const T& s) :_s(s){}
	__host__ __device__ T operator()(const T&a) const
	{
		return _s * a;
	}
};
/// \brief The functor called by std::transform or thrust::transform, given the parameter scale, that output the x + scale * y, where x, y are two vectors, the result is also a vector
template<typename T>
struct _saxpy_functor
{
	T _scal;
	_saxpy_functor(const T& s) :_scal(s){}
	__host__ __device__ T operator()(const T& x, const T& y) const
	{
		return x + _scal * y;
	}
};


/// \brief The self defined functor called by thrust::transform or std::transform which prevent the denominator is 0;
template<typename T>
struct _divide_functor
{
	__host__ __device__ T operator()(const T& a, const T& b)
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
};



/// \brief The function helps to get the lp norm of a device or host vector
template<typename T>
struct _lp_functor
{
private:
	T _p;
public:
	_lp_functor(const T& p) :_p(p){}
	__host__ __device__ T operator()(const T& x)
	{
		return pow(MY_ABS<T>(x), _p);
	}
};

/// \brief the functor helps to get |a-b| for two input value a and b, it is the first step to calculate the distance of two vectors
template<typename T>
struct _abs_minus_functor
{

public:
	__host__ __device__ T operator()(const T& a, const T& b)
	{
		return MY_ABS<T>(a - b);
	}
};



template<typename T>
struct _Update
{
	typedef thrust::tuple<T, T, T, T> TU;
	T _lamb;
	T _keepRange;
	T _minV;
	T _maxV;
	_Update(const T& lambda, bool keepInRange, const T& minV, const T& maxV) :
		_lamb(lambda), _minV(minV), _maxV(maxV), _keepRange(keepInRange){}
	__host__ __device__ T operator()(TU t)
	{
		T img = thrust::get<0>(t);
		T cor = thrust::get<1>(t);
		T weg = thrust::get<2>(t);
		T msk = thrust::get<3>(t);
		img = (img + _lamb * cor / weg) * msk;
		if (_keepRange)
		{
			if (img < _minV)
			{
				return _minV;
			}
			if (img > _maxV)
			{
				return _maxV;
			}
			return img;
		}
		return img;
	}
};


/// \brief Inverse transform for the soft thresholding filtering called by OptimumMU function
template<typename T>
struct _softThreshold_functor
{
	const T _g;
	_softThreshold_functor(const T& g) :_g(g){}
	__host__ __device__ T operator()(const T& r) const
	{
		if (r < -_g)
		{
			return r + _g;
		}
		else if (r > _g)
		{
			return r - _g;
		}
		else
		{
			return 0;
		}
	}
};


/// \brief updating the data in CG algorithm with functor 1
template<typename T>
struct CG_update1_functor
{
public:
	T _alpha;
	CG_update1_functor(const T& alpha) :_alpha(alpha){}
	__host__ __device__ T operator()(const T& X, const T& D){
		return X + _alpha * D;
	}
};

/// \brief updating the data in CG algorithm with functor 2
template<typename T>
struct CG_update2_functor
{
public:
	T _alpha;
	CG_update2_functor(const T& alpha) :_alpha(alpha){}
	__host__ __device__ T operator()(const T& R, const T& tem2)
	{
		return R - _alpha * tem2;
	}
};

/// \brief updating the data in CG algorithm with functor 3
template<typename T>
struct CG_update3_functor
{
public:
	T _beta;
	CG_update3_functor(const T& beta) :_beta(beta){}
	__host__ __device__ T operator()(const T& R, const T& D)
	{
		return R + _beta * D;
	}
};



/// \brief FISTA accelerate the convergence.
template<typename T>
struct _FISTA_demo11
{
public:
	T t1;
	T t2;
	_FISTA_demo11(const T& _t1, const T& _t2) :t1(_t1), t2(_t2){}
	__host__ __device__ T operator()(const T& x1, const T& x0)
	{
		return x1 + (t1 - 1) / t2 * (x1 - x0);
	}
};



template<typename T>
struct _ValueLimit_functor
{
	T _UpLim;
	T _DownLim;
	_ValueLimit_functor(const T& up, const T& down) :_UpLim(up), _DownLim(down){}
	__host__ __device__ T operator()(const T& V)
	{
		if (V > _UpLim)
		{
			return _UpLim;
		}
		else if (V < _DownLim)
		{
			return _DownLim;
		}
		else
			return V;

	}
};





template<typename T>
struct ossart_update_functor
{
	typedef thrust::tuple<T, T, T, T> TU;
	T _lamb;
	ossart_update_functor(const T& lambda) :_lamb(lambda){}
	__device__ T operator()(TU t)
	{
		T img = thrust::get<0>(t);
		T cor = thrust::get<1>(t);
		T weg = thrust::get<2>(t);
		T msk = thrust::get<3>(t);
		if (!IS_ZERO<T>(weg))
		{
			return (img + _lamb * cor / weg) * msk;
		}
		else
		{
			return 0;
		}
	}
};




template<typename T>
struct _ADD_TUPLE_functor
{
	typedef thrust::tuple<T, T> TP;
	__host__ __device__ TP operator()(const TP& a, const TP& b)
	{
		T agk = thrust::get<0>(a);
		T agk1 = thrust::get<1>(a);
		T bgk = thrust::get<0>(b);
		T bgk1 = thrust::get<1>(b);
		return thrust::make_tuple(agk + bgk, agk1 + bgk1);
	}
};

template<typename T, int cases>
struct _betaupdate_functor
{
	typedef thrust::tuple<T, T, T> ITP;
	typedef thrust::tuple<T, T> OTP;
	__host__ __device__ OTP operator()(const ITP& a)
	{
		T gk = thrust::get<0>(a);
		T gk1 = thrust::get<1>(a);
		T d = thrust::get<2>(a);
		switch (cases)
		{
		case 0:
			return thrust::make_tuple(gk1 * gk1, gk * gk);
		case 1:
			return thrust::make_tuple(gk1 * (gk1 - gk), gk * gk);
		case 2:
			return thrust::make_tuple(gk1 *(gk1 - gk), d*(gk1 - gk));
		case 3:
			return thrust::make_tuple(gk1 * gk1, -d * gk);
		case 4:
			return thrust::make_tuple(gk1 * gk1, d * (gk1 - gk));
		default:
			return thrust::make_tuple(gk1 * gk1, gk * gk);
		}
	}
};



/// Update the image with functor
template<typename T>
struct _DEMO17_update_functor
{
	typedef thrust::tuple<T, T, T> TP;
	T coef;
	T minV;
	T maxV;
	_DEMO17_update_functor(const T& c, const T& miv, const T& mav) :coef(c), minV(miv), maxV(mav){}
	_DEMO17_update_functor(const T& c, const T& miv) :coef(c), minV(miv), maxV(10000.0){}
	_DEMO17_update_functor(const T& c) :coef(c), minV(0), maxV(10000.0){}
	__host__ __device__ T operator()(const TP& fdw)
	{
		T f = thrust::get < 0 >(fdw);
		T d = thrust::get < 1 >(fdw);
		T w = thrust::get < 2 >(fdw);
		T res = 0;
		if (!IS_ZERO(w))
		{
			res = f + coef * d / w;
		}
		else
		{
			res = minV;
		}
		if (res > maxV)
		{
			res = maxV;
		}
		if (res < minV)
		{
			res = minV;
		}
		return res;
	}
};




/// Update the image with functor
template<typename T>
struct _DEMO18_update_functor
{
	typedef thrust::tuple<T, T, T> TP;
	T coef;
	T minV;
	T maxV;
	_DEMO18_update_functor(const T& c, const T& miv, const T& mav) :coef(c), minV(miv), maxV(mav){}
	__host__ __device__ T operator()(const TP& fdw)
	{
		T f = thrust::get < 0 >(fdw);
		T d = thrust::get < 1 >(fdw);
		T w = thrust::get < 2 >(fdw);
		T res = 0;
		if (!IS_ZERO(w))
		{
			res = f + coef * d / w;
		}
		else
		{
			res = minV;
		}
		if (res > maxV)
		{
			res = maxV;
		}
		if (res < minV)
		{
			res = minV;
		}
		return res;
	}
};



template<typename T>
struct _DEMO18v4_GenProjLambda_functor
{
	T _maxY;
	_DEMO18v4_GenProjLambda_functor(const T& maxy) :_maxY(maxy){}
	__host__ __device__ T operator()(const T& yi, const T& weight)
	{
		if (IS_ZERO(weight))
		{
			return 0;
		}
		return yi / (_maxY * weight);
	}
};
template<typename T>
struct _DEMO18v4_GenBackLambda_functor
{
	T _val;
	_DEMO18v4_GenBackLambda_functor(const T& vl) :_val(vl){}
	__host__ __device__ T operator()(const T& v)
	{
		if (IS_ZERO(v))
		{
			return 0;
		}
		return _val / v;
	}
};

template<typename T>
struct _DEMO18v4_updateImgSIR_functor
{
	T _lambda;
	_DEMO18v4_updateImgSIR_functor(const T& lamb) :_lambda(lamb){}
	__host__ __device__ T operator()(const T& f, const T& upf)
	{
		return f - _lambda * upf;
	}
};


template<typename T>
class Comparer_functor
{
private:
	enum CMP{ LESS = -1, EQUI = 0, BIGG = 1 };

	int m_comType;
public:
	Comparer_functor(const int comType)
	{
		m_comType = comType;
	}
	__host__ __device__ bool operator()(T& num1, T& num2) const{
		bool res;
		switch (m_comType)
		{
		case LESS:
			res = (num1 < num2);
			break;
		case EQUI:
			res = (num1 == num2);
			break;
		case BIGG:
			res = (num1 > num2);
			break;
		default:
			res = false;
			break;
		}
		return res;
	}
};





//function nu = WindowMeyer(xi,deg)
//	% WindowMeyer -- auxiliary window function for Meyer wavelets.
//	%  Usage
//	%    nu = WindowMeyer(xi,deg)
//	%  Inputs
//	%    xi     abscissa values for window evaluation
//	%    deg    degree of the polynomial defining Nu on [0,1]
//	%           1 <= deg <= 3
//	%  Outputs
//	%    nu     polynomial of degree 'deg' if x in [0,1]
//	%           1 if x > 1 and 0 if x < 0.
//	%  See Also
//	%    MeyerPartition
template<typename T>
struct _windowMeyer_functor1
{
	T __host__ __device__ operator()(const T& xi)
	{
		T res = xi * xi * (3.0 - 2.0 * xi); //The coefficients are gotten from the Wikipedia with Filter Window in signal processing
		if (xi <= 0)
		{
			res = 0;
		}
		else if (xi >= 1)
		{
			res = 1;
		}

		return res;
	}
};

template<typename T>
struct _windowMeyer_functor2
{
	T __host__ __device__ operator()(const T& xi)
	{
		T res = xi * xi * xi * (10.0 - 15.0 * xi + 6.0 * xi * xi);//The coefficients are gotten from the Wikipedia with Filter Window in signal processing
		if (xi <= 0)
		{
			res = 0;
		}
		else if (xi >= 1)
		{
			res = 1;
		}
		return res;
	}
};
template<typename T>
struct _windowMeyer_functor3
{
	T __host__ __device__ operator()(const T& xi)
	{
		T res = xi * xi * xi * xi * (35.0 - 84.0 * xi + 70.0 * xi * xi - 20.0 * xi * xi * xi);//The coefficients are gotten from the Wikipedia with Filter Window in signal processing
		if (xi <= 0)
		{
			res = 0;
		}
		else if (xi >= 1)
		{
			res = 1;
		}
		return res;
	}
};

template<typename T, int deg>
struct _windowMeyer_functor
{
	T __host__ __device__ operator()(const T& xi)
	{
		T res(0);
		switch (deg)
		{
		case 0:
			res = xi;
			break;
		case 1:
			res = xi * xi * (3.0 - 2.0 * xi);//The coefficients are gotten from the Wikipedia with Filter Window in signal processing
			break;
		case 2:
			res = xi * xi * xi * (10.0 - 15.0 * xi + 6.0 * xi * xi);//The coefficients are gotten from the Wikipedia with Filter Window in signal processing
			break;
		case 3:
			res = xi * xi * xi * xi * (35.0 - 84.0 * xi + 70.0 * xi * xi - 20.0 * xi * xi * xi);//The coefficients are gotten from the Wikipedia with Filter Window in signal processing
			break;
		}

		if (xi <= 0)
		{
			res = 0;
		}
		else if (xi >= 1)
		{
			res = 1;
		}
		return res;

	}
};

template<typename T, int deg>
struct _windowMeyerCos_functor
{
	T __host__ __device__ operator()(const T& xi)
	{
		T res(0);
		switch (deg)
		{
		case 0:
			res = xi;
			break;
		case 1:
			res = xi * xi * (3.0 - 2.0 * xi);//The coefficients are gotten from the Wikipedia with Filter Window in signal processing
			break;
		case 2:
			res = xi * xi * xi * (10.0 - 15.0 * xi + 6.0 * xi * xi);//The coefficients are gotten from the Wikipedia with Filter Window in signal processing
			break;
		case 3:
			res = xi * xi * xi * xi * (35.0 - 84.0 * xi + 70.0 * xi * xi - 20.0 * xi * xi * xi);//The coefficients are gotten from the Wikipedia with Filter Window in signal processing
			break;
		}

		if (xi <= 0)
		{
			res = 0;
		}
		else if (xi >= 1)
		{
			res = 1;
		}
		return cos(res * _PI_2);

	}
};



template<typename T>
void WindowMeyer(thrust::host_vector<T>& nu, // output
	thrust::host_vector<T>& xi, // input
	const int deg)
{
	switch (deg)
	{
	case 0:
		nu = xi;
		break;
	case 1:
		thrust::transform(xi.begin(), xi.end(), nu.begin(), _windowMeyer_functor<T, 1>());
		break;
	case 2:
		thrust::transform(xi.begin(), xi.end(), nu.begin(), _windowMeyer_functor<T, 2>());
		break;
	case 3:
		thrust::transform(xi.begin(), xi.end(), nu.begin(), _windowMeyer_functor<T, 3>());
		break;
	}
}


template<typename T>
void WindowMeyer(thrust::device_vector<T>& nu, thrust::device_vector<T>& xi, const int deg)
{
	switch (deg)
	{
	case 0:
		nu = xi;
		break;
	case 1:
		thrust::transform(xi.begin(), xi.end(), nu.begin(), _windowMeyer_functor<T, 1>());
		break;
	case 2:
		thrust::transform(xi.begin(), xi.end(), nu.begin(), _windowMeyer_functor<T, 2>());
		break;
	case 3:
		thrust::transform(xi.begin(), xi.end(), nu.begin(), _windowMeyer_functor<T, 3>());
		break;
	}
}



/// \brief The beta update formula for Conjugate Gradient method. Please DO NOT modify this function unless you exactly know what you are doing.
template<typename T, int cases = 0>
T calculateBetaForCG(thrust::device_vector<T>& gk, thrust::device_vector<T>& gk1, thrust::device_vector<T>& d)
{
	thrust::tuple<T, T> init = thrust::make_tuple<T>(0.0f, 0.0f);
	thrust::tuple<T, T> res = thrust::transform_reduce(
		thrust::make_zip_iterator(thrust::make_tuple(gk.begin(), gk1.begin(), d.begin())),
		thrust::make_zip_iterator(thrust::make_tuple(gk.end(), gk1.end(), d.end())),
		_betaupdate_functor < T, cases>(), init, _ADD_TUPLE_functor<T>());
	return thrust::get<0>(res) / thrust::get<1>(res);
}



/// \brief pointwise add
/// \param a input vector a
/// \param b input vector b
/// \return a + b
template<typename T>
thrust::device_vector<T> operator+(thrust::device_vector<T>& a, thrust::device_vector<T>& b)
{
	thrust::device_vector<T> c(a.size(), 0);
	thrust::transform(a.begin(), a.end(), b.begin(), c.begin(), thrust::plus<T>());
	return c;
}

/// \brief pointwise add
/// \param a input vector a
/// \param b input vector b
/// \return a + b
template<typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
	std::vector<T> c(a.size(), 0);
	std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::plus<T>());
	return c;
}
/// \brief pointwise add
/// \param a input vector a
/// \param b input vector b
/// \return a + b
template<typename T>
thrust::host_vector<T> operator+(thrust::host_vector<T>& a, thrust::host_vector<T>& b)
{
	thrust::host_vector<T> c(a.size(), 0);
	thrust::transform(a.begin(), a.end(), b.begin(), c.begin(), thrust::plus<T>());
	return c;
}

/// \brief pointwise minus
/// \param a input vector a
/// \param b input vector b
/// \return a - b
template<typename T>
thrust::device_vector<T> operator-(thrust::device_vector<T>& a, thrust::device_vector<T>& b)
{
	thrust::device_vector<T> c(a.size(), 0);
	thrust::transform(a.begin(), a.end(), b.begin(), c.begin(), thrust::minus<T>());
	return c;
}

/// \brief pointwise minus
/// \param a input vector a
/// \param b input vector b
/// \return a - b
template<typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b)
{
	std::vector<T> c(a.size(), 0);
	std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::minus<T>());
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



/// \brief scale-vector multiplication
/// \param a input vector
/// \param s scale in R
/// \return the scaled vector of a
template<typename T>
thrust::device_vector<T> operator*(thrust::device_vector<T>& a, const T& s)
{
	thrust::device_vector<T> res(a.size(), 0);
	thrust::transform(a.begin(), a.end(), res.begin(), _scale_functor<T>(s));
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
	thrust::transform(a.begin(), a.end(), res.begin(), _scale_functor<T>(s));
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
	thrust::transform(a.begin(), a.end(), b.begin(), c.begin(), _divide_functor<T>());
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
	thrust::transform(lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(), _divide_functor<T>());
	return lhs;
}



template<typename T>
thrust::device_vector<T>& operator*=(thrust::device_vector<T>& lhs, const T& rhs)
{
	thrust::transform(lhs.begin(), lhs.end(), lhs.begin(), _scale_functor<T>(rhs));
	return lhs;
}


/// \brief return the p norm of a vector
/// \param v the vector
/// \param p the p norm
/// \return (sum(v^p))^(1/p)
template<typename T>
T operator^(const thrust::device_vector<T>& v, const T& p)
{
	return std::pow(thrust::transform_reduce(v.begin(), v.end(), _lp_functor<T>(p), 0.0, thrust::plus<T>()));
}

/// \brief Operator overload |, Get abs(a - b) vector of two vectors a and b
/// \param lhs left vector
/// \param rhs right vector
/// \return abs(lhs - rhs)
template<typename T>
thrust::device_vector<T> operator|(
	const thrust::device_vector<T>& lhs,
	const thrust::device_vector<T>& rhs)
{
	thrust::device_vector<T> res(lhs.size(), 0);
	thrust::transform(lhs.begin(), lhs.end(), rhs.begin(), res.begin(), _abs_minus_functor<T>());
	return res;
}


/// \brief the rotation of the 2D vector according to cos(T) and sin(T)
/// \param p original vector
/// \param cosT cos(T)
/// \param sinT sin(T)
inline __host__	__device__ float2 rotation(const float2& p, const float& cosT, const float& sinT)
{
	float2 curP;
	curP.x = p.x * cosT - p.y * sinT;
	curP.y = p.x * sinT + p.y * cosT;
	return curP;
}

/// \brief the rotation of the 2D vector according to cos(T) and sin(T)
/// \param p original vector
/// \param cosT cos(T)
/// \param sinT sin(T)
inline __host__	__device__ double2 rotation(const double2& p, const double& cosT, const double& sinT)
{
	double2 curP;
	curP.x = p.x * cosT - p.y * sinT;
	curP.y = p.x * sinT + p.y * cosT;
	return curP;
}

/// \brief the inverse rotation of the 2D vector according to cos(T) and sin(T)
/// \param p original vector
/// \param cosT cos(T)
/// \param sinT sin(T)
inline __host__ __device__ float2 invRotation(const float2& p, const float& cosT, const float& sinT)
{
	float2 curP;
	curP.x = p.x * cosT + p.y * sinT;
	curP.y = -p.x * sinT + p.y * cosT;
	return curP;
}


/// \brief the inverse rotation of the 2D vector according to cos(T) and sin(T)
/// \param p original vector
/// \param cosT cos(T)
/// \param sinT sin(T)
inline __host__ __device__ double2 invRotation(const double2& p, const double& cosT, const double& sinT)
{
	double2 curP;
	curP.x = p.x * cosT + p.y * sinT;
	curP.y = -p.x * sinT + p.y * cosT;
	return curP;
}

/// \brief the rotation of the 3D vector according to cos(T) and sin(T)
/// \param p original vector
/// \param cosT cos(T)
/// \param sinT sin(T)
inline __host__ __device__ float3 rotation(const float3& p, const float& cosT, const float& sinT)
{
	float3 curP;
	curP.x = p.x * cosT - p.y * sinT;
	curP.y = p.x * sinT + p.y * cosT;
	curP.z = p.z;
	return curP;
}



/// \brief the rotation of the 3D vector according to cos(T) and sin(T)
/// \param p original vector
/// \param cosT cos(T)
/// \param sinT sin(T)
inline __host__ __device__ double3 rotation(const double3& p, const double& cosT, const double& sinT)
{
	double3 curP;
	curP.x = p.x * cosT - p.y * sinT;
	curP.y = p.x * sinT + p.y * cosT;
	curP.z = p.z;
	return curP;
}



/// \brief the inverse rotation of the 3D vector according to cos(T) and sin(T)
/// \param p original vector
/// \param cosT cos(T)
/// \param sinT sin(T)
inline __host__ __device__ float3 invRotation(const float3& p, const float& cosT, const float& sinT)
{
	float3 curP;
	curP.x = p.x * cosT + p.y * sinT;
	curP.y = -p.x * sinT + p.y * cosT;
	curP.z = p.z;
	return curP;
}


/// \brief the inverse rotation of the 3D vector according to cos(T) and sin(T)
/// \param p original vector
/// \param cosT cos(T)
/// \param sinT sin(T)
inline __host__ __device__ double3 invRotation(const double3& p, const double& cosT, const double& sinT)
{
	double3 curP;
	curP.x = p.x * cosT + p.y * sinT;
	curP.y = -p.x * sinT + p.y * cosT;
	curP.z = p.z;
	return curP;
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
	T Qmm, T Qmb, //Í¬ÑùµÄX;
	T Qbm, T Qbb) //Í¬ÑùµÄY;
{
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


/// \brief SIDDON kernel function 1
inline __host__	__device__ float dev_pFun(const float& alpha, const float& pstart, const float&pend)
{
	return pstart + alpha * (pend - pstart);
}
/// \brief SIDDON kernel function 1
inline __host__	__device__ double dev_pFun(const double& alpha, const double& pstart, const double&pend)
{
	return pstart + alpha * (pend - pstart);
}


/// \brief SIDDON kernel function 2
inline __host__	__device__ float dev_alpha_IFun(const float& b, const float& d, const float& pstart, const float& pend, const unsigned int& i)
{
	if (!IS_ZERO<float>(pend - pstart))
	{
		return ((b + (float) i*d) - pstart) / (pend - pstart);
	}
	else return 1000;//((b + i*d)-pstart)/(1e-6);
}
/// \brief SIDDON kernel function 2
inline __host__	__device__ double dev_alpha_IFun(const double& b, const double& d, const double& pstart, const double& pend, const unsigned int& i)
{
	if (!IS_ZERO<double>(pend - pstart))
	{
		return ((b + (double) i*d) - pstart) / (pend - pstart);
	}
	else return 1000;//((b + i*d)-pstart)/(1e-6);
}



/// \brief SIDDON kernel function 3
inline __host__	__device__ float dev_varphiFun(const float& alpha, const float& b, const float& d, const float& pstart, const float& pend)
{
	return (dev_pFun(alpha, pstart, pend) - b) / d;
}
/// \brief SIDDON kernel function 3
inline __host__	__device__ double dev_varphiFun(const double& alpha, const double& b, const double& d, const double& pstart, const double& pend)
{
	return (dev_pFun(alpha, pstart, pend) - b) / d;
}



/// \brief SIDDON kernel function 4
inline	__host__ __device__ void dev_minmaxIdxFun(
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
inline __host__ __device__  float dev_alphaU_Fun(const float& d, const float& startx, const float& endx)
{
	if (IS_ZERO<float>(startx - endx))
	{
		return 1000.0f;//(d/1e-6);
	}
	return d / fabsf(startx - endx);
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


/// \brief SIDDON kernel function 6
template<typename T>
inline __host__ __device__  int dev_iu_Fun(const T& start, const T& end)
{
	return (start < end) ? 1 : -1;
}


/// \brief SIDDON line integral function in 2D
inline	__host__ __device__ float calSiddonOneRayKer2D(
	const float& startX, const float& startY,
	const float& endX, const float& endY,
	const float& __MINOL__, const float& __MINOW__,
	const float& __OBJSTPL__, const float& __OBJSTPW__,
	const unsigned int& __OBJLR__, const unsigned int& __OBJWR__,
	float* dev_vol, float* totWeight)
{
	float dirX = endX - startX;
	float dirY = endY - startY;
	const float dconv = sqrtf(dirX * dirX + dirY * dirY);
	int imin(0), imax(0), jmin(0), jmax(0);

	float alphaxmin = 0.0f;
	float alphaxmax = 0.0f;
	float alphaymin = 0.0f;
	float alphaymax = 0.0f;

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

	const float alphaMIN = MY_MAX<float>(alphaxmin, alphaymin); // (alphaxmin > alphaymin) ? alphaxmin : alphaymin;
	const float alphaMAX = MY_MIN<float>(alphaxmax, alphaymax); // (alphaxmax < alphaymax) ? alphaxmax : alphaymax;
	dev_minmaxIdxFun(startX, endX, __MINOL__, __OBJSTPL__, alphaMIN, alphaMAX, alphaxmin, alphaxmax, __OBJLR__ + 1, &imin, &imax);
	dev_minmaxIdxFun(startY, endY, __MINOW__, __OBJSTPW__, alphaMIN, alphaMAX, alphaymin, alphaymax, __OBJWR__ + 1, &jmin, &jmax);

	float alphaX = (startX < endX) ? dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, imin) : dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, imax);
	float alphaY = (startY < endY) ? dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, jmin) : dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, jmax);

	int Np = static_cast<int>(MY_ABS<float>(static_cast<float>(imax - imin)) + MY_ABS<float>(static_cast<float>(jmax - jmin)) + 3.0f); //fabsf(imax - imin) + fabsf(jmax - jmin) + 3.0f;
	const float alphaxu = dev_alphaU_Fun(__OBJSTPL__, startX, endX);
	const float alphayu = dev_alphaU_Fun(__OBJSTPW__, startY, endY);

	float alphaC = alphaMIN;
	//float minalpha = MY_MIN<float>(alphaX,alphaY); //min(alphaX, alphaY);
	//float talpha = MY_MIN<float>(alphaX,alphaY); //fminf(alphaX,alphaY);
	int i = int(dev_varphiFun(alphaMIN*1.00003f, __MINOL__, __OBJSTPL__, startX, endX));
	int j = int(dev_varphiFun(alphaMIN*1.00003f, __MINOW__, __OBJSTPW__, startY, endY));
	//int i = floor(dev_varphiFun((talpha + alphaMIN)*0.5, __MINOL__, __OBJSTPL__, startX, endX));
	//int j = floor(dev_varphiFun((talpha + alphaMIN)*0.5, __MINOW__, __OBJSTPW__, startY, endY));

	const int iuu = (startX < endX) ? 1 : -1;
	const int juu = (startY < endY) ? 1 : -1;

	float d12(0);
	float weight(0);


	for (int repIdx(0); repIdx < Np; ++repIdx)
	{
		if (i < 0 || i >= static_cast<int>(__OBJLR__) || j < 0 || j >= static_cast<int>(__OBJWR__))
		{
			break;
		}
		if (alphaX <= alphaY)
		{
			weight = (alphaX - alphaC) * dconv;
			d12 += weight * dev_vol[j * __OBJLR__ + i];
			i += iuu;
			(*totWeight) += weight;
			alphaC = alphaX;
			alphaX += alphaxu;
			continue;
		}
		else
		{
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


/// \brief SIDDON line integral function in 2D
inline	__host__ __device__ double calSiddonOneRayKer2D(
	const double& startX, const double& startY,
	const double& endX, const double& endY,
	const double& __MINOL__, const double& __MINOW__,
	const double& __OBJSTPL__, const double& __OBJSTPW__,
	const unsigned int& __OBJLR__, const unsigned int& __OBJWR__,
	double* dev_vol, double* totWeight)
{
	double dirX = endX - startX;
	double dirY = endY - startY;
	const double dconv = sqrt(dirX * dirX + dirY * dirY);
	int imin(0), imax(0), jmin(0), jmax(0);

	double alphaxmin = 0.0f;
	double alphaxmax = 0.0f;
	double alphaymin = 0.0f;
	double alphaymax = 0.0f;

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

	const double alphaMIN = MY_MAX<double>(alphaxmin, alphaymin); // (alphaxmin > alphaymin) ? alphaxmin : alphaymin;
	const double alphaMAX = MY_MIN<double>(alphaxmax, alphaymax); // (alphaxmax < alphaymax) ? alphaxmax : alphaymax;
	dev_minmaxIdxFun(startX, endX, __MINOL__, __OBJSTPL__, alphaMIN, alphaMAX, alphaxmin, alphaxmax, __OBJLR__ + 1, &imin, &imax);
	dev_minmaxIdxFun(startY, endY, __MINOW__, __OBJSTPW__, alphaMIN, alphaMAX, alphaymin, alphaymax, __OBJWR__ + 1, &jmin, &jmax);

	double alphaX = (startX < endX) ? dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, imin) : dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, imax);
	double alphaY = (startY < endY) ? dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, jmin) : dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, jmax);

	int Np = static_cast<int>(MY_ABS<double>(static_cast<double>(imax - imin)) + MY_ABS<double>(static_cast<double>(jmax - jmin)) + 3.0f); //fabsf(imax - imin) + fabsf(jmax - jmin) + 3.0f;
	const double alphaxu = dev_alphaU_Fun(__OBJSTPL__, startX, endX);
	const double alphayu = dev_alphaU_Fun(__OBJSTPW__, startY, endY);

	double alphaC = alphaMIN;
	//double minalpha = MY_MIN<double>(alphaX,alphaY); //min(alphaX, alphaY);
	//double talpha = MY_MIN<double>(alphaX,alphaY); //fminf(alphaX,alphaY);
	int i = int(dev_varphiFun(alphaMIN*1.00003f, __MINOL__, __OBJSTPL__, startX, endX));
	int j = int(dev_varphiFun(alphaMIN*1.00003f, __MINOW__, __OBJSTPW__, startY, endY));
	//int i = floor(dev_varphiFun((talpha + alphaMIN)*0.5, __MINOL__, __OBJSTPL__, startX, endX));
	//int j = floor(dev_varphiFun((talpha + alphaMIN)*0.5, __MINOW__, __OBJSTPW__, startY, endY));

	const int iuu = (startX < endX) ? 1 : -1;
	const int juu = (startY < endY) ? 1 : -1;

	double d12(0);
	double weight(0);


	for (int repIdx(0); repIdx < Np; ++repIdx)
	{
		if (i < 0 || i >= static_cast<int>(__OBJLR__) || j < 0 || j >= static_cast<int>(__OBJWR__))
		{
			break;
		}
		if (alphaX <= alphaY)
		{
			weight = (alphaX - alphaC) * dconv;
			d12 += weight * dev_vol[j * __OBJLR__ + i];
			i += iuu;
			(*totWeight) += weight;
			alphaC = alphaX;
			alphaX += alphaxu;
			continue;
		}
		else
		{
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




/// \brief SIDDON line integral function in 3D
inline __host__ __device__ float calSiddonOneRayKer(
	const float& startX, const float& startY, const float& startZ,
	const float& endX, const float& endY, const float& endZ,
	const float& __MINOL__, const float& __MINOW__, const float& __MINOH__,
	const float& __OBJSTPL__, const float& __OBJSTPW__, const float& __OBJSTPH__,
	cuint& __OBJLR__, cuint& __OBJWR__, cuint& __OBJHR__,
	float* dev_vol, float* totWeight)
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


	const float alphaMIN = MY_MAX<float>(alphaxmin, alphaymin, alphazmin);
	const float alphaMAX = MY_MIN<float>(alphaxmax, alphaymax, alphazmax);
	dev_minmaxIdxFun(startX, endX, __MINOL__, __OBJSTPL__, alphaMIN, alphaMAX, alphaxmin, alphaxmax, __OBJLR__ + 1, &imin, &imax);
	dev_minmaxIdxFun(startY, endY, __MINOW__, __OBJSTPW__, alphaMIN, alphaMAX, alphaymin, alphaymax, __OBJWR__ + 1, &jmin, &jmax);
	dev_minmaxIdxFun(startZ, endZ, __MINOH__, __OBJSTPH__, alphaMIN, alphaMAX, alphazmin, alphazmax, __OBJHR__ + 1, &kmin, &kmax);

	float alphaX = (startX < endX) ? dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, imin) : dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, imax);
	float alphaY = (startY < endY) ? dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, jmin) : dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, jmax);
	float alphaZ = (startZ < endZ) ? dev_alpha_IFun(__MINOH__, __OBJSTPH__, startZ, endZ, kmin) : dev_alpha_IFun(__MINOH__, __OBJSTPH__, startZ, endZ, kmax);

	int Np = static_cast<int>(MY_ABS<float>(static_cast<float>(imax - imin)) + MY_ABS<float>(static_cast<float>(jmax - jmin)) + MY_ABS<float>(static_cast<float>(kmax - kmin)) + 3.0f);
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

	const int iuu = dev_iu_Fun(startX, endX);
	const int juu = dev_iu_Fun(startY, endY);
	const int kuu = dev_iu_Fun(startZ, endZ);

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
				d12 = d12 + weight * dev_vol[(k * __OBJWR__ + j) * __OBJLR__ + i];
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
				d12 = d12 + weight * dev_vol[(k * __OBJWR__ + j) * __OBJLR__ + i];
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
				d12 = d12 + weight * dev_vol[(k * __OBJWR__ + j) * __OBJLR__ + i];
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

inline __host__ __device__ float calSiddonOneRayKer_Zfirst(
	const float& startX, const float& startY, const float& startZ,
	const float& endX, const float& endY, const float& endZ,
	const float& __MINOL__, const float& __MINOW__, const float& __MINOH__,
	const float& __OBJSTPL__, const float& __OBJSTPW__, const float& __OBJSTPH__,
	cuint& __OBJLR__, cuint& __OBJWR__, cuint& __OBJHR__,
	float* dev_vol, float* totWeight)
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


	const float alphaMIN = MY_MAX<float>(alphaxmin, alphaymin, alphazmin);
	const float alphaMAX = MY_MIN<float>(alphaxmax, alphaymax, alphazmax);
	dev_minmaxIdxFun(startX, endX, __MINOL__, __OBJSTPL__, alphaMIN, alphaMAX, alphaxmin, alphaxmax, __OBJLR__ + 1, &imin, &imax);
	dev_minmaxIdxFun(startY, endY, __MINOW__, __OBJSTPW__, alphaMIN, alphaMAX, alphaymin, alphaymax, __OBJWR__ + 1, &jmin, &jmax);
	dev_minmaxIdxFun(startZ, endZ, __MINOH__, __OBJSTPH__, alphaMIN, alphaMAX, alphazmin, alphazmax, __OBJHR__ + 1, &kmin, &kmax);

	float alphaX = (startX < endX) ? dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, imin) : dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, imax);
	float alphaY = (startY < endY) ? dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, jmin) : dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, jmax);
	float alphaZ = (startZ < endZ) ? dev_alpha_IFun(__MINOH__, __OBJSTPH__, startZ, endZ, kmin) : dev_alpha_IFun(__MINOH__, __OBJSTPH__, startZ, endZ, kmax);

	int Np = static_cast<int>(MY_ABS<float>(static_cast<float>(imax - imin)) + MY_ABS<float>(static_cast<float>(jmax - jmin)) + MY_ABS<float>(static_cast<float>(kmax - kmin)) + 3.0f);
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

	const int iuu = dev_iu_Fun(startX, endX);
	const int juu = dev_iu_Fun(startY, endY);
	const int kuu = dev_iu_Fun(startZ, endZ);

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
				d12 = d12 + weight * dev_vol[(j * __OBJLR__ + i) * __OBJHR__ + k];
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
				d12 = d12 + weight * dev_vol[(j * __OBJLR__ + i) * __OBJHR__ + k];
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
				d12 = d12 + weight * dev_vol[(j * __OBJLR__ + i) * __OBJHR__ + k];
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


/// \brief SIDDON line integral function in 3D
inline __host__ __device__ double calSiddonOneRayKer(
	const double& startX, const double& startY, const double& startZ,
	const double& endX, const double& endY, const double& endZ,
	const double& __MINOL__, const double& __MINOW__, const double& __MINOH__,
	const double& __OBJSTPL__, const double& __OBJSTPW__, const double& __OBJSTPH__,
	cuint& __OBJLR__, cuint& __OBJWR__, cuint& __OBJHR__,
	double* dev_vol, double* totWeight)
{
	double dirX = endX - startX;
	double dirY = endY - startY;
	double dirZ = endZ - startZ;

	const double dconv = sqrt(dirX * dirX + dirY * dirY + dirZ * dirZ);
	int imin(0), imax(0), jmin(0), jmax(0), kmin(0), kmax(0);

	double alphaxmin = 0.0f;
	double alphaxmax = 0.0f;
	double alphaymin = 0.0f;
	double alphaymax = 0.0f;
	double alphazmin = 0.0f;
	double alphazmax = 0.0f;

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


	const double alphaMIN = MY_MAX<double>(alphaxmin, alphaymin, alphazmin);
	const double alphaMAX = MY_MIN<double>(alphaxmax, alphaymax, alphazmax);
	dev_minmaxIdxFun(startX, endX, __MINOL__, __OBJSTPL__, alphaMIN, alphaMAX, alphaxmin, alphaxmax, __OBJLR__ + 1, &imin, &imax);
	dev_minmaxIdxFun(startY, endY, __MINOW__, __OBJSTPW__, alphaMIN, alphaMAX, alphaymin, alphaymax, __OBJWR__ + 1, &jmin, &jmax);
	dev_minmaxIdxFun(startZ, endZ, __MINOH__, __OBJSTPH__, alphaMIN, alphaMAX, alphazmin, alphazmax, __OBJHR__ + 1, &kmin, &kmax);

	double alphaX = (startX < endX) ? dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, imin) : dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, imax);
	double alphaY = (startY < endY) ? dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, jmin) : dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, jmax);
	double alphaZ = (startZ < endZ) ? dev_alpha_IFun(__MINOH__, __OBJSTPH__, startZ, endZ, kmin) : dev_alpha_IFun(__MINOH__, __OBJSTPH__, startZ, endZ, kmax);

	int Np = static_cast<int>(MY_ABS<double>(static_cast<double>(imax - imin)) + MY_ABS<double>(static_cast<double>(jmax - jmin)) + MY_ABS<double>(static_cast<double>(kmax - kmin)) + 3.0f);
	const double alphaxu = dev_alphaU_Fun(__OBJSTPL__, startX, endX);
	const double alphayu = dev_alphaU_Fun(__OBJSTPW__, startY, endY);
	const double alphazu = dev_alphaU_Fun(__OBJSTPH__, startZ, endZ);

	double alphaC = alphaMIN;
	//const double minApa = MY_MIN<double>(alphaX,alphaY,alphaZ);

	int i = int(dev_varphiFun(alphaMIN*1.00003f, __MINOL__, __OBJSTPL__, startX, endX));
	int j = int(dev_varphiFun(alphaMIN*1.00003f, __MINOW__, __OBJSTPW__, startY, endY));
	int k = int(dev_varphiFun(alphaMIN*1.00003f, __MINOH__, __OBJSTPH__, startZ, endZ));
	//int i = floor(dev_varphiFun((alphaMIN+minApa) * 0.5f, __MINOL__, __OBJSTPL__, startX, endX));
	//int j = floor(dev_varphiFun((alphaMIN+minApa) * 0.5f, __MINOW__, __OBJSTPW__, startY, endY));
	//int k = floor(dev_varphiFun((alphaMIN+minApa) * 0.5f, __MINOH__, __OBJSTPH__, startZ, endZ));

	const int iuu = dev_iu_Fun(startX, endX);
	const int juu = dev_iu_Fun(startY, endY);
	const int kuu = dev_iu_Fun(startZ, endZ);

	double d12(0);
	double weight(0);


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
				d12 = d12 + weight * dev_vol[(k * __OBJWR__ + j) * __OBJLR__ + i];
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
				d12 = d12 + weight * dev_vol[(k * __OBJWR__ + j) * __OBJLR__ + i];
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
				d12 = d12 + weight * dev_vol[(k * __OBJWR__ + j) * __OBJLR__ + i];
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



/// \brief intersection length of the rectangle box
/// \param r ray
/// \param boxmin min boundary of the box
/// \param boxmax max boundary of the obx
/// \param tnear pointer to the near intersection parameter
/// \param tfar pointer to the far intersection parameter
/// \return intersection or not
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


/// \brief intersection length of the 3D rectangle box
/// \param r ray
/// \param boxmin min boundary of the box
/// \param boxmax max boundary of the obx
/// \param tnear pointer to the near intersection parameter
/// \param tfar pointer to the far intersection parameter
/// \return intersection or not
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







//////////////////////////////////////////////////////////////////////////
//! This is the function used for test the number x is the power of 2 or
//! not, this is a fast algorithm!
//! Author: Rui LIU
//! Date: 2013-3-17
//! Version: x.x
//! @param x       input value
//! @return        whether the input parameter is power of 2
//////////////////////////////////////////////////////////////////////////
template<typename T>
inline __host__ __device__ bool isPow2(const T& x)
{
	return ((x&(x - 1)) == 0);
}



/// \brief Define the CT reconstruction utilities functions. This is a fast algorithm for get the number P that P >= x, P = 2^n. If x = 3, P = 4; if x = 16, P = 16. This function only useful for unsigned int datatype.
/// \version 1.0
/// \return The number of power 2 that not smaller than x
/// \param x input value
inline __host__ __device__ unsigned int nextPow2(const unsigned int& x)
{
	unsigned int t(x);
	--t;
	t |= t >> 1;
	t |= t >> 2;
	t |= t >> 4;
	t |= t >> 8;
	t |= t >> 16;
	return ++t;
}
//! This is a fast algorithm for get the number P that
//!                 P >= x, P = 2^n
//! If x = 3, P = 4; if x = 16, P = 16
//! This function only useful for unsigned int datatype
//! @version 1.0
//! @return         the number of power 2 that not smaller than x
inline __host__ __device__ void nextPow2(
	//! input the variable x for calculation
	const unsigned int&x,
	//! output the result to *t
	unsigned int* t)
{
	*t = x;
	--(*t);
	(*t) |= (*t) >> 1;
	(*t) |= (*t) >> 2;
	(*t) |= (*t) >> 4;
	(*t) |= (*t) >> 8;
	(*t) |= (*t) >> 16;
	(*t) = (*t) + 1;
}



inline void stripnl(char *str) {
	while (strlen(str) && ((str[strlen(str) - 1] == 13) ||
		(str[strlen(str) - 1] == 10))) {
		str[strlen(str) - 1] = 0;
	}
}

inline int countNumLines1(char *filename) //from http://ubuntuforums.org/archive/index.php/t-1057474.html
{
	FILE *f;
	char c;
	int lines = 0;

	f = fopen(filename, "r");
	if (f == nullptr)
		return 0;

	while ((c = fgetc(f)) != EOF)
		if (c == '\n')
			lines++;

	fclose(f);
	if (c != '\n')
		lines++;

	return lines;
}
inline int countNumLines2(char *s) //from: http://www.cplusplus.com/doc/tutorial/files/
{
	std::string line;
	std::ifstream myfile(s);
	int lines = 0;
	if (myfile.is_open())
	{
		while (myfile.good())
		{
			getline(myfile, line);
			if (!line.empty())
				lines++;
		}
		myfile.close();
	}
	else
		lines = 0;
	return lines;
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
inline void pushMatrix(
	float startX, float startY,
	float endX, float endY,
	float bx, float by,
	float dx, float dy,
	cuint objResoLen,
	cuint objResoWid,
	int& rowidx,
	std::vector<int>& rowIdx,
	std::vector<int>& colIdx,
	std::vector<float>& wgt)
{
	const float dirX(endX - startX);
	const float dirY(endY - startY);
	const float lengthSq = dirX * dirX + dirY * dirY;
	const float dconv = sqrt(lengthSq);
	int imin, imax, jmin, jmax;
	const float alphaxmin = MY_MIN<float>(dev_alpha_IFun(bx, dx, startX, endX, 0), dev_alpha_IFun(bx, dx, startX, endX, objResoLen));
	const float alphaxmax = MY_MAX<float>(dev_alpha_IFun(bx, dx, startX, endX, 0), dev_alpha_IFun(bx, dx, startX, endX, objResoLen));
	const float alphaymin = MY_MIN<float>(dev_alpha_IFun(by, dy, startY, endY, 0), dev_alpha_IFun(by, dy, startY, endY, objResoWid));
	const float alphaymax = MY_MAX<float>(dev_alpha_IFun(by, dy, startY, endY, 0), dev_alpha_IFun(by, dy, startY, endY, objResoWid));

	const float alphaMIN = MY_MAX<float>(alphaxmin, alphaymin);
	const float alphaMAX = MY_MIN<float>(alphaxmax, alphaymax);
	dev_minmaxIdxFun(startX, endX, bx, dx, alphaMIN, alphaMAX, alphaxmin, alphaxmax, objResoLen + 1, &imin, &imax);
	dev_minmaxIdxFun(startY, endY, by, dy, alphaMIN, alphaMAX, alphaymin, alphaymax, objResoWid + 1, &jmin, &jmax);

	float alphaX = (startX < endX) ? dev_alpha_IFun(bx, dx, startX, endX, imin) : dev_alpha_IFun(bx, dx, startX, endX, imax);
	float alphaY = (startY < endY) ? dev_alpha_IFun(by, dy, startY, endY, jmin) : dev_alpha_IFun(by, dy, startY, endY, jmax);

	int Np = static_cast<int>(MY_ABS<float>(imax - imin + 1.0f) + MY_ABS<float>(jmax - jmin + 1.0f) + 4.0f);
	const float alphaxu = dev_alphaU_Fun(dx, startX, endX);
	const float alphayu = dev_alphaU_Fun(dy, startY, endY);

	float alphaC = alphaMIN;
	float minApa = MY_MIN<float>(alphaX, alphaY);
	/*int i = static_cast<int>(dev_varphiFun(alphaMIN* 1.00003f, bx, dx, startX, endX));
	int j = static_cast<int>(dev_varphiFun(alphaMIN* 1.00003f, by, dy, startY, endY));
	*/
	int i = static_cast<int>(dev_varphiFun((alphaMIN + minApa) * 0.5f, bx, dx, startX, endX));
	int j = static_cast<int>(dev_varphiFun((alphaMIN + minApa) * 0.5f, by, dy, startY, endY));

	const int iuu = dev_iu_Fun(startX, endX);
	const int juu = dev_iu_Fun(startY, endY);

	float d12(0.0f);
	float weight(0.0f);
	unsigned int repIdx(0);
	unsigned int colidx(0);
	while (repIdx != Np)
	{
		if (i < 0 || i >= static_cast<int>(objResoLen) || j < 0 || j >= static_cast<int>(objResoWid))
		{
			break;
		}
		if (alphaX <= alphaY)
		{
			colidx = j * objResoLen + i;
			weight = (alphaX - alphaC) * dconv;
			wgt.push_back(weight);
			rowIdx.push_back(rowidx);
			colIdx.push_back(colidx);
			d12 += weight;
			i += iuu;
			alphaC = alphaX;
			alphaX += alphaxu;
		}
		else
		{
			colidx = j * objResoLen + i;
			weight = (alphaY - alphaC) * dconv;
			wgt.push_back(weight);
			rowIdx.push_back(rowidx);
			colIdx.push_back(colidx);
			d12 += weight;
			j += juu;
			alphaC = alphaY;
			alphaY += alphayu;
		}
		++repIdx;
	}
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
inline void pushMatrix(
	double startX, double startY,
	double endX, double endY,
	double bx, double by,
	double dx, double dy,
	cuint objResoLen,
	cuint objResoWid,
	int& rowidx,
	std::vector<int>& rowIdx,
	std::vector<int>& colIdx,
	std::vector<double>& wgt)
{
	const double dirX(endX - startX);
	const double dirY(endY - startY);
	const double lengthSq = dirX * dirX + dirY * dirY;
	const double dconv = sqrt(lengthSq);
	int imin, imax, jmin, jmax;
	const double alphaxmin = MY_MIN<double>(dev_alpha_IFun(bx, dx, startX, endX, 0), dev_alpha_IFun(bx, dx, startX, endX, objResoLen));
	const double alphaxmax = MY_MAX<double>(dev_alpha_IFun(bx, dx, startX, endX, 0), dev_alpha_IFun(bx, dx, startX, endX, objResoLen));
	const double alphaymin = MY_MIN<double>(dev_alpha_IFun(by, dy, startY, endY, 0), dev_alpha_IFun(by, dy, startY, endY, objResoWid));
	const double alphaymax = MY_MAX<double>(dev_alpha_IFun(by, dy, startY, endY, 0), dev_alpha_IFun(by, dy, startY, endY, objResoWid));

	const double alphaMIN = MY_MAX<double>(alphaxmin, alphaymin);
	const double alphaMAX = MY_MIN<double>(alphaxmax, alphaymax);
	dev_minmaxIdxFun(startX, endX, bx, dx, alphaMIN, alphaMAX, alphaxmin, alphaxmax, objResoLen + 1, &imin, &imax);
	dev_minmaxIdxFun(startY, endY, by, dy, alphaMIN, alphaMAX, alphaymin, alphaymax, objResoWid + 1, &jmin, &jmax);

	double alphaX = (startX < endX) ? dev_alpha_IFun(bx, dx, startX, endX, imin) : dev_alpha_IFun(bx, dx, startX, endX, imax);
	double alphaY = (startY < endY) ? dev_alpha_IFun(by, dy, startY, endY, jmin) : dev_alpha_IFun(by, dy, startY, endY, jmax);

	int Np = static_cast<int>(MY_ABS<double>(imax - imin + 1.0f) + MY_ABS<double>(jmax - jmin + 1.0f) + 4.0f);
	const double alphaxu = dev_alphaU_Fun(dx, startX, endX);
	const double alphayu = dev_alphaU_Fun(dy, startY, endY);

	double alphaC = alphaMIN;
	double minApa = MY_MIN<double>(alphaX, alphaY);
	/*int i = static_cast<int>(dev_varphiFun(alphaMIN* 1.00003f, bx, dx, startX, endX));
	int j = static_cast<int>(dev_varphiFun(alphaMIN* 1.00003f, by, dy, startY, endY));
	*/
	int i = static_cast<int>(dev_varphiFun((alphaMIN + minApa) * 0.5f, bx, dx, startX, endX));
	int j = static_cast<int>(dev_varphiFun((alphaMIN + minApa) * 0.5f, by, dy, startY, endY));

	const int iuu = dev_iu_Fun(startX, endX);
	const int juu = dev_iu_Fun(startY, endY);

	double d12(0.0f);
	double weight(0.0f);
	unsigned int repIdx(0);
	unsigned int colidx(0);
	while (repIdx != Np)
	{
		if (i < 0 || i >= static_cast<int>(objResoLen) || j < 0 || j >= static_cast<int>(objResoWid))
		{
			break;
		}
		if (alphaX <= alphaY)
		{
			colidx = j * objResoLen + i;
			weight = (alphaX - alphaC) * dconv;
			wgt.push_back(weight);
			rowIdx.push_back(rowidx);
			colIdx.push_back(colidx);
			d12 += weight;
			i += iuu;
			alphaC = alphaX;
			alphaX += alphaxu;
		}
		else
		{
			colidx = j * objResoLen + i;
			weight = (alphaY - alphaC) * dconv;
			wgt.push_back(weight);
			rowIdx.push_back(rowidx);
			colIdx.push_back(colidx);
			d12 += weight;
			j += juu;
			alphaC = alphaY;
			alphaY += alphayu;
		}
		++repIdx;
	}
}




template <typename T>
inline __host__ __device__ void calSVASVB(T* SVA, T* SVB, T* sour, const T& cosT, const T& sinT, const FanEAGeo& FanGeo, const Image& Img, cuint detIdx)
{
	T pangle = (detIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
	T initY = -cos(pangle) * FanGeo.m_S2D + FanGeo.m_S2O;
	T initX = sin(pangle) * FanGeo.m_S2D;
	T curX = initX * cosT - initY * sinT;
	T curY = initX * sinT + initY * cosT;
	initX = curX - sour[0];
	initY = curY - sour[1];
	T legth = sqrt(initX * initX + initY * initY);
	SVA[0] = initX / legth;
	SVA[1] = initY / legth;
	SVA[2] = static_cast<T>(detIdx);

	pangle = pangle + FanGeo.m_DetStp;
	initY = -cos(pangle) * FanGeo.m_S2D + FanGeo.m_S2O;
	initX = sin(pangle) * FanGeo.m_S2D;
	curX = initX * cosT - initY * sinT;
	curY = initX * sinT + initY * cosT;
	initX = curX - sour[0];
	initY = curY - sour[1];
	legth = sqrt(initX * initX + initY * initY);
	SVB[0] = initX / legth;
	SVB[1] = initY / legth;
	SVB[2] = static_cast<T>(detIdx) +1;
}


template<typename T>
__host__ __device__ void calSVASVB(T* SVA, T* SVB, T* sour, const T& cosT, const T& sinT, const FanEDGeo& FanGeo, const Image& Img, cuint detIdx)
{
	T initX = (detIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
	T initY = -FanGeo.m_O2D;
	T curX = initX * cosT - initY * sinT;
	T curY = initX * sinT + initY * cosT;
	initX = curX - sour[0];
	initY = curY - sour[1];
	T legth = sqrt(initX * initX + initY * initY);
	SVA[0] = initX / legth;
	SVA[1] = initY / legth;
	SVA[2] = static_cast<T>(detIdx);


	initX = (detIdx - FanGeo.m_DetCntIdx + 1) * FanGeo.m_DetStp;
	initY = -FanGeo.m_O2D;
	curX = initX * cosT - initY * sinT;
	curY = initX * sinT + initY * cosT;
	initX = curX - sour[0];
	initY = curY - sour[1];
	legth = sqrt(initX * initX + initY * initY);
	SVB[0] = initX / legth;
	SVB[1] = initY / legth;
	SVB[2] = static_cast<T>(detIdx + 1.0);
}




template<typename T>
inline __host__ __device__ void SortProj(T(&grid)[4][3])
{
	int i(0), j(0);
	T td(0.0);
	for (i = 0; i != 3; ++i)
	{
		for (j = i + 1; j != 4; ++j)
		{
			if (grid[j][2] < grid[i][2])
			{
				td = grid[i][0];
				grid[i][0] = grid[j][0];
				grid[j][0] = td;

				td = grid[i][1];
				grid[i][1] = grid[j][1];
				grid[j][1] = td;

				td = grid[i][2];
				grid[i][2] = grid[j][2];
				grid[j][2] = td;
			}
		}
	}
}


template<typename T>
__host__ __device__ inline void SortProjection(T(&Grid)[4][3])
{
	int i, j;
	T td;
	//mexPrintf("S0=%d,S1=%d,S2=%d,S3=%d\n",SSort[0],SSort[1],SSort[2],SSort[3]);
	for (i = 0; i < 3; i++)
	{
		for (j = i + 1; j < 4; j++)
		{
			if (Grid[j][2] < Grid[i][2])
			{
				td = Grid[i][0];
				Grid[i][0] = Grid[j][0];
				Grid[j][0] = td;

				td = Grid[i][1];
				Grid[i][1] = Grid[j][1];
				Grid[j][1] = td;

				td = Grid[i][2];
				Grid[i][2] = Grid[j][2];
				Grid[j][2] = td;
			}
		}
	}
}
template<typename T>
inline __host__ __device__ T ComputeCoefficient(
	const T(&Grid)[4][3],
	const T(&SVA)[3],
	const T(&SVB)[3],
	const T(&SPoint)[2],
	const T area)
{
	T coef = 0;
	T x0(0), y0(0), a(0), b(0), t(0);
	int AI, BI;
	if (SVA[2] > Grid[3][2])
	{
		return 0;
	}
	if (SVB[2] < Grid[0][2])
	{
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

	switch (AI)
	{
	case 0:
	{
		switch (BI)
		{
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
		{    x0 = Grid[3][0] - SPoint[0];
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
		switch (BI)
		{
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
		switch (BI)
		{
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
	for (yi = 0; yi < YN; ++yi)
	{
		for (xi = 0; xi < XN; ++xi)
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
				temp = ang + _PI - pangle;
				SVA[0] = cos(temp);
				SVA[1] = sin(temp);
				SVA[2] = di;
				// mexPrintf("di=%d,VA0=%10.8f,VA1=%10.8f,angle=%10.8f,Beta=%10.8f\n",di,SVA[0],SVA[1],temp,pangle);
				pangle = PBeta[di] + 0.5 * DtBeta;
				temp = ang + _PI - pangle;
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


template<typename T>
void genProjectionMatrix_AIM_template(const FanEDGeo FanGeo, const Image Img)
{
	unsigned int angIdx = 0;
	std::vector < int> rowIdx;
	std::vector < int > colIdx;
	std::vector < T > coeffs;
	T ang(0);
	for (angIdx = 0; angIdx < FanGeo.m_ViwN; ++angIdx)
	{
		ang = FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp;
		genProj_AIM<T>(rowIdx, colIdx, coeffs, angIdx, ang, FanGeo, Img);
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
	rowFile.write((char*) &(rowIdx[0]), sizeof(int) * rowIdx.size());
	rowFile.close();
	colFile.write((char*) &(colIdx[0]), sizeof(int) * colIdx.size());
	colFile.close();
	coeFile.write((char*) &(coeffs[0]), sizeof(T) * coeffs.size());
	coeFile.close();

	rowFile.close();
	colFile.close();
	coeFile.close();

}

template<typename T>
void genProjectionMatrix_AIM_template(const FanEAGeo FanGeo, const Image Img)
{
	int angIdx = 0;
	std::vector < int> rowIdx;
	std::vector < int > colIdx;
	std::vector < T > coeffs;
	T ang(0);
	for (angIdx = 0; angIdx < FanGeo.m_ViwN; ++angIdx)
	{
		ang = FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp;
		genProj_AIM<T>(rowIdx, colIdx, coeffs, angIdx, ang, FanGeo, Img);
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
	rowFile.write((char*) &(rowIdx[0]), sizeof(int) * rowIdx.size());
	rowFile.close();
	colFile.write((char*) &(colIdx[0]), sizeof(int) * colIdx.size());
	colFile.close();
	coeFile.write((char*) &(coeffs[0]), sizeof(T) * coeffs.size());
	coeFile.close();

	rowFile.close();
	colFile.close();
	coeFile.close();
}




/// \brief Generate the Hilbert Filter used by analytical reconstruction algorithm
/// \version 1.0
/// \date 2014-02-18
/// \param filt Filter
/// \param sigLeng Signal Length
template<typename T>
void hilbertFilter(
	std::vector<T>& filt,
	const unsigned int& sigLeng)
{
	filt.clear();
	filt.resize(2 * sigLeng - 1);
	filt[0] = 0;
	for (unsigned int i(1); i != sigLeng + 1; ++i)
	{
		filt[i] = (1 - std::pow(-1.0, i)) / (3.14159265358979323846264 * i);
	}
	for (unsigned int i(sigLeng); i != 2 * sigLeng - 1; ++i)
	{
		filt[i] = (1 - std::pow(-1.0, static_cast<T>(i) -2.0 *
			static_cast<T>(sigLeng) +1.0)) /
			(3.14159265358979323846264 * (static_cast<T>(i) -2.0 * static_cast<T>(sigLeng) +1.0));
	}
}




//Calculate the l0 norm of a device thrust vector in GPU. And points out where are these non-zeros are.
template<typename T>
struct CalNorm0_functor
{
	typedef thrust::tuple<T, int> MyTuple;
	__host__ __device__ bool operator()(const MyTuple& a)
	{
		T val = thrust::get<0>(a);
		int idx = thrust::get<1>(a);
		return abs(val) < 1.0E-20;
	}
};

template<typename T>
int calZeroNorm(const thrust::device_vector<T>& v, thrust::device_vector<int>& IDX)
{
	thrust::device_vector<T> vB = v;
	int Size = vB.size();
	IDX.clear();	IDX.resize(Size);
	thrust::sequence(IDX.begin(), IDX.end(), 0);
	auto resultIter = thrust::remove_if(
		thrust::make_zip_iterator(thrust::make_tuple(vB.begin(), IDX.begin())),
		thrust::make_zip_iterator(thrust::make_tuple(vB.end(), IDX.end())),
		CalNorm0_functor<T>());
	// Really erase the elements
	vB.erase(thrust::get<0>(resultIter.get_iterator_tuple()), vB.end());
	IDX.erase(thrust::get<1>(resultIter.get_iterator_tuple()), IDX.end());
	return vB.size();
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
void DD3Boundaries(int nrBoundaries,const T*pCenters, T *pBoundaries)
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
	int i;
	T* pBoundaries = &Boundaries[0];
	T* pCenters = &Centers[0];
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
void DD3Boundaries(int nrBoundaries, thrust::host_vector<T>& Centers, thrust::host_vector<T>& Boundaries)
{
	DD3Boundaries<T>(nrBoundaries, &Centers[0], &Boundaries[0]);
}




template<typename T>
T SIGN(const T& v)
{
	if (v > 0)
	{
		return 1;
	}
	else if (v < 0)
	{
		return -1;
	}
	else
	{
		return 0;
	}
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


template<typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T>& dt)
{
	for (unsigned int i = 0; i != dt.size(); ++i)
	{
		os << dt[i] << " ";
	}
	os << std::endl;
	return os;
}


template<typename T>
const T mean(std::vector<T>& v)
{
	T summ = 0;
	for (unsigned int i = 0; i != v.size(); ++i)
	{
		summ += v[i];
	}
	if (v.size() == 0)
	{
		return 0;
	}
	else
	{
		return (summ / static_cast<T>(v.size()));
	}
}


template<typename T>
const T norm(std::vector<T>& v)
{
	return std::sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.0));
}

template<typename T>
void KroneckerTensorProduct(
	std::vector<T>& res,
	const std::vector<T>& A,
	const int mA, const int nA,
	const std::vector<T>& B,
	const int mB, const int nB)
{
	res.resize(mA * nA * mB * nB, 0);
	int iA(0), jA(0), iB(0), jB(0), idx(0);
	T coefA, coefB;
	for (jA = 0; jA != nA; ++jA)
	{
		for (iA = 0; iA != mA; ++iA)
		{
			coefA = A[jA * mA + iA];
			for (jB = 0; jB != nB; ++jB)
			{
				for (iB = 0; iB != mB; ++iB)
				{
					coefB = B[jB * mB + iB];
					// calculate idx
					idx = (nB * jA + jB) * (mA * mB) + (mB * iA + iB);
					res[idx] = coefA * coefB;
				}
			}
		}
	}
}





// The projection data update for SART/OSSART liked algorithm with the OMP
template<typename T>
void prjWeight(T* prj, // projection data to be updated
		T* realPrj, // The real projection data
		T* rowSum, // The row sum of the projection weighting
		int N) // total elements number
{
#pragma omp parallel for
	for (int i = 0; i < N; ++i)
	{
		if (rowSum[i] > _EPSILON)
		{
			prj[i] = (realPrj[i] - prj[i]) / rowSum[i];
		}
		else
		{
			prj[i] = 0;
		}
	}
}

template<typename T>
void prjWeight(std::vector<T>& prj, std::vector<T>& realPrj, std::vector<T>& rowSum)
{
	struct prjFunctor
	{
		typedef thrust::tuple<T,T,T> InputTuple;
		T operator()(InputTuple& in)
		{
			T prj_ = thrust::get<0>(in);
			T realPrj_ = thrust::get<1>(in);
			T rowSum_ = thrust::get<2>(in);
			if(rowSum_> _EPSILON)
			{
				return (realPrj_ - prj_) / rowSum_;
			}
			else
			{
				return 0;
			}
		}

	};

	thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(prj.begin(),realPrj.begin(),rowSum.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(prj.end(),realPrj.end(),rowSum.end())),
			prj.begin(),prjFunctor());
}

// The backprojection data update for SART/OSSART liked algorithm with the OMP
template<typename T>
void bakWeight(T* vol, // temporary backprojected volume
		T* reconImg, // volume to be updated
		T* colSum, // The col sum of the volume weighting
		int N)  // total elements number
{
#pragma omp parallel for
	for (int i = 0; i < N; ++i)
	{
		if (colSum[i] > _EPSILON)
		{
			vol[i] = vol[i] / colSum[i];
		}
		else
		{
			vol[i] = 0;
		}
		reconImg[i] = reconImg[i] + vol[i];
	}

}




template<typename T>
void bakWeight(std::vector<T>& vol, std::vector<T>& reconImg, std::vector<T>& colSum)
{
	struct bakFunctor
	{
		typedef thrust::tuple<T,T,T> InputTuple;
		T operator()(InputTuple& in)
		{
			T vol_ = thrust::get<0>(in);
			T reconImg_ = thrust::get<1>(in);
			T colSum_ = thrust::get<2>(in);
			if(colSum_> _EPSILON)
			{
				vol_ = vol_ / colSum_;
			}
			else
			{
				vol_ = 0;
			}
			return reconImg_ + vol_;
		}
	};

	thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(vol.begin(),reconImg.begin(),colSum.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(vol.end(),reconImg.end(),colSum.end())),
			reconImg.begin(),bakFunctor());
}



template<typename T>
void FISTA(T* lasImg, T* currentImg, T t1, T t2, int N)
{
	struct FISTA_functor
	{
		T t1;
		T t2;
		FISTA_functor(const T& _t1, const T& _t2) :t1(_t1), t2(_t2){}
		__host__ __device__ T operator()(T curImg, T lasImg)
		{
			return curImg + (t1 - 1.0) / t2 * (curImg - lasImg);
		}

	};
	thrust::transform(currentImg, currentImg + N, lasImg, currentImg, FISTA_functor(t1, t2));
}




template<typename T>
void FISTA(std::vector<T>& lasImg, std::vector<T>&  currentImg, T t1, T t2, int N)
{
	struct FISTA_functor
	{
		T t1;
		T t2;
		FISTA_functor(const T& _t1, const T& _t2) :t1(_t1), t2(_t2){}
		__host__ __device__ T operator()(T curImg, T lasImg)
		{
			return curImg + (t1 - 1.0) / t2 * (curImg - lasImg);
		}

	};
	thrust::transform(currentImg.begin(), currentImg.end(), lasImg.begin(), currentImg, FISTA_functor(t1, t2));
}


template<typename T>
__device__ inline T intersectLength_device(const T& fixedmin, const T& fixedmax, const T& varimin, const T& varimax)
{
	const T left = (fixedmin > varimin) ? fixedmin : varimin;
	const T right = (fixedmax < varimax) ? fixedmax : varimax;
	return fabsf(right - left) * static_cast<T>(right > left);
}



template<typename T>
__host__ __device__ inline T intersectLength(const T& fixedmin, const T& fixedmax, const T& varimin, const T& varimax)
{
	const T left = (fixedmin > varimin) ? fixedmin : varimin;
	const T right = (fixedmax < varimax) ? fixedmax : varimax;
	return abs(right - left) * static_cast<double>(right > left);
}