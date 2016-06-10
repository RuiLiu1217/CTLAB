/*
 * FastIterativeShrinkageAlgorithm.h
 *
 *  Created on: Jun 8, 2016
 *      Author: liurui
 */

#ifndef FASTITERATIVESHRINKAGEALGORITHM_H_
#define FASTITERATIVESHRINKAGEALGORITHM_H_

 // Moving to DD3_GPU_RECON
#include <thrust/transform.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <vector>

template<typename T>
struct FISTA_functor
{
	T t1;
	T t2;
	T minV;
	T maxV;
	FISTA_functor(const T& _t1, const T& _t2, const T& _minV, const T& _maxV) :t1(_t1), t2(_t2), minV(_minV), maxV(_maxV){}
	__host__ __device__ T operator()(T curImg, T lasImg)
	{
		T res = curImg + (t1 - 1.0) / t2 * (curImg - lasImg);
		if (res < minV)
		{
			return minV;
		}
		else if (res > maxV)
		{
			return maxV;
		}
		else
		{
			return res;
		}
	}

};

template<typename T>
void FISTA(T* lasImg, T* currentImg, T t1, T t2, int N, T minV = -2000, T maxV = 6000)
{
	thrust::transform(currentImg, currentImg + N, lasImg, currentImg, FISTA_functor<float>(t1, t2, minV, maxV));
}


template<typename T>
void FISTA(thrust::device_vector<T>& lasImg, thrust::device_vector<T>& currentImg, T t1, T t2, T minV = -2000.0f, T maxV = 6000.0f)
{

	thrust::transform(currentImg.begin(), currentImg.end(), lasImg.begin(), currentImg.begin(), FISTA_functor<T>(t1, t2, minV, maxV));
}

template<typename T>
void FISTA(thrust::host_vector<T>& lasImg, thrust::host_vector<T>& currentImg,T t1, T t2, T minV = -2000.0f, T maxV = 6000.0f)
{
	thrust::transform(currentImg.begin(), currentImg.end(), lasImg.begin(), currentImg.begin(), FISTA_functor<T>(t1, t2, minV, maxV));
}


template<typename T>
void FISTA(std::vector<T>& lasImg, std::vector<T>& currentImg,T t1, T t2, T minV = -2000.0f, T maxV = 6000.0f)
{
	thrust::transform(currentImg.begin(), currentImg.end(), lasImg.begin(), currentImg.begin(), FISTA_functor<T>(t1, t2, minV, maxV));
}



#endif /* FASTITERATIVESHRINKAGEALGORITHM_H_ */
