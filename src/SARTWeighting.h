/*!
 * \file SARTWeighting.h
 *
 * \date Jun 8, 2016
 * \author liurui
 * \version 1.0
 */
#ifndef SARTWEIGHTING_H_
#define SARTWEIGHTING_H_

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/tuple.h>
#include <thrust/iterator/zip_iterator.h>
#include <omp.h>
/// \brief The static class for the SART weighting
template<typename Type>
class SART_Weighting
{
public:
	static void Proj(thrust::host_vector<Type>& prj, const thrust::host_vector<Type>& realPrj, const thrust::host_vector<Type>& rowSum)
	{
		thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(prj.begin(),realPrj.begin(),rowSum.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(prj.end(),realPrj.end(),rowSum.end())),
			prj.begin(),prjWeight_functor<Type>());
	}

	static void Proj(thrust::device_vector<Type>& prj, const thrust::device_vector<Type>& realPrj, const thrust::device_vector<Type>& rowSum)
	{
		thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(prj.begin(),realPrj.begin(),rowSum.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(prj.end(),realPrj.end(),rowSum.end())),
			prj.begin(),prjWeight_functor<Type>());
	}

	static void Proj_ptr(Type* prj, const Type* realPrj, const Type* rowSum,const int N)
	{
		omp_set_num_threads(32);
	#pragma omp parallel for
		for (int i = 0; i < N; ++i)
		{
			if (rowSum[i] > 1.0e-7)
			{
				prj[i] = (realPrj[i] - prj[i]) / rowSum[i];
			}
			else
			{
				prj[i] = 0;
			}
		}
	}


	static void Back(thrust::host_vector<Type>& vol, thrust::host_vector<Type>& reconImg, const thrust::host_vector<Type>& colSum)
	{
		thrust::transform(
				thrust::make_zip_iterator(thrust::make_tuple(vol.begin(),reconImg.begin(),colSum.begin())),
				thrust::make_zip_iterator(thrust::make_tuple(vol.end(),reconImg.end(),colSum.end())),
				reconImg.begin(),bakWeight_functor<Type>());
	}

	static void Back(thrust::device_vector<Type>& vol, thrust::device_vector<Type>& reconImg, const thrust::device_vector<Type>& colSum)
	{
		thrust::transform(
				thrust::make_zip_iterator(thrust::make_tuple(vol.begin(),reconImg.begin(),colSum.begin())),
				thrust::make_zip_iterator(thrust::make_tuple(vol.end(),reconImg.end(),colSum.end())),
				reconImg.begin(),bakWeight_functor<Type>());
	}



	static void Back_ptr(Type* vol,Type* reconImg, const Type* colSum, const int N)
	{
		omp_set_num_threads(32);
		#pragma omp parallel for
			for (int i = 0; i < N; ++i)
			{
				if (colSum[i] > 1.0e-7)
				{
					vol[i] /= colSum[i];
				}
				else
				{
					vol[i] = 0;
				}
				reconImg[i] += vol[i];
			}

	}


private:
	static const float epsilon;// = 1.0e-7;
public:
	template<typename T>
	struct prjWeight_functor
	{
		typedef thrust::tuple<T,T,T> inputTuple;
		__host__ __device__ T operator()(const inputTuple& input)
		{
			T prj = thrust::get<0>(input);
			T realPrj = thrust::get<1>(input);
			T rowSum = thrust::get<2>(input);
			if(rowSum > epsilon)
			{
				return (realPrj - prj) / rowSum;
			}
			else
			{
				return 0;
			}
		}
	};
	template<typename T>
	struct bakWeight_functor
	{
		typedef thrust::tuple<T,T,T> inputType;
		__host__ __device__ T operator()(const inputType& input)
		{
			T vol = thrust::get<0>(input);
			T reconImg = thrust::get<1>(input);
			T colSum = thrust::get<2>(input);
			if(colSum > epsilon)
			{
				vol = vol / colSum;

			}
			else
			{
				vol = 0;
			}
			return reconImg + vol;
		}
	};
};
template<typename Type>const float SART_Weighting<Type>::epsilon = 1.0e-7;
#endif /* SARTWEIGHTING_H_ */
