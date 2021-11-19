#include "utilities.hpp"
#include "threeGPUHouseholder.hpp"
#include "cudaCheckReturner.h"

namespace SVD{
template<typename T>
__host__ __device__ inline T SIGN(const T& y)
{
	if (y > 0)
	{
		return 1.0f;
	}
	else if (y < 0)
	{
		return -1.0f;
	}
	else
	{
		return 0;
	}
}

template<typename T>
__global__ void calculateQB(T* B, T* hvec, T* Bhvec, const int m, const int n)
{
	const int i = threadIdx.x + blockIdx.x * blockDim.x; 
	const int j = threadIdx.y + blockIdx.y * blockDim.y; 
	if (i < m && j < n)
	{
		B[j * m + i] = B[j * m + i] - 2.0 * hvec[i] * Bhvec[j];
	}
}
template<typename T>
__global__ void calculateBP(T* B, T* lvec, T* lvecB, const int m, const int n)
{
	const int i = threadIdx.x + blockIdx.x * blockDim.x;
	const int j = threadIdx.y + blockIdx.y * blockDim.y;
	if (i < m && j < n)
	{
		B[j * m + i] = B[j * m + i] - 2.0 * lvecB[i] * lvec[j];
	}
}

template<typename T>
__global__ void CopyBm(T* B, T* Bm, const int m, const int k)
{
	const int i = threadIdx.x + blockDim.x * blockIdx.x;
	if (i < m)
	{
		if (i >= k)
		{
			Bm[i] = B[k * m + i];
		}
		else
		{
			Bm[i] = 0;
		}
	}
}

template<typename T>
__global__ void CopyBn(T* B, T* Bn, const int m, const int n, const int k)
{
	const int i = threadIdx.x + blockDim.x * blockIdx.x;
	if (i < n)
	{
		if (i > k)
		{
			Bn[i] = B[(k + 1 + i) * m + k];
		}
		else
		{
			Bn[i] = 0;
		}
	}
}

template<typename T>
__global__ void updateU(T* U, T* Uhvec, T* hvec, const int m)
{
	const int i = threadIdx.x + blockIdx.x * blockDim.x;
	const int j = threadIdx.y + blockIdx.y * blockDim.y;
	if (i < m && j < m)
	{
		U[j * m + i] = U[j * m + i] - 2.0 * Uhvec[i] * hvec[j];
	}
}


template<typename T>
__global__ void updateV(T* V, T* Vlvec, T* lvec, const int n)
{
	const int i = threadIdx.x + blockIdx.x * blockDim.x;
	const int j = threadIdx.y + blockIdx.y * blockDim.y;
	if (i < n && j < n)
	{
		V[j * n + i] = V[j * n + i] - 2.0 * Vlvec[i] * lvec[j];
	}
}


void HouseHolder_v2(thrust::device_vector<float>& B, const int m, const int n)
{
	cudaStream_t streams[3];
	float summ(0);	float len(0); float tempV(0); float alpha(0); float gamma(0);
	const float ONE = 1.0;
	const float ZERO = 0.0;
	CUDA_CHECK_RETURN(cudaSetDevice(0));
	CUDA_CHECK_RETURN(cudaStreamCreate(&streams[0]));
	cublasHandle_t cublashandle;
	cublasCreate(&cublashandle);
	thrust::device_vector<float> hvec(m, 0);
	thrust::device_vector<float> lvec(n, 0);



	CUDA_CHECK_RETURN(cudaSetDevice(1)); //Used to update Left part
	CUDA_CHECK_RETURN(cudaStreamCreate(&streams[1]));
	thrust::device_vector<float> hvec1(m, 0);
	thrust::device_vector<float> Uhvec1(m, 0);
	thrust::device_vector<float> U(m*m, 0);
	for (size_t i = 0; i < m; i++)
	{
		U[i * m + i] = 1.0;
	}
	cublasHandle_t cublashandle1;
	cublasCreate(&cublashandle1);



	CUDA_CHECK_RETURN(cudaSetDevice(2)); //Used to update right part
	CUDA_CHECK_RETURN(cudaStreamCreate(&streams[2]));
	thrust::device_vector<float> lvec2(n, 0);
	thrust::device_vector<float> Vlvec2(n, 0);
	thrust::device_vector<float> V(n * n, 0);
	for (size_t i = 0; i != n; ++i)
	{
		V[i * n + i] = 1.0;
	}
	cublasHandle_t cublashandle2;
	cublasCreate(&cublashandle2);

	dim3 blk(16, 16);
	dim3 gid(
		(m + blk.x - 1) / blk.x,
		(n + blk.y - 1) / blk.y);

	dim3 mblk(16, 16);
	dim3 mgid(
		(m + mblk.x - 1) / mblk.x,
		(m + mblk.y - 1) / mblk.y);

	dim3 nblk(16, 16);
	dim3 ngid(
		(n + nblk.x - 1) / nblk.x,
		(n + nblk.y - 1) / nblk.y);


	for (size_t k = 0; k < n; k++)
	{
		std::cout << k << std::endl;
		CUDA_CHECK_RETURN(cudaSetDevice(0));
		//Copy Vector;
		summ = 0;
		for (size_t i = 0; i < m - k; i++)
		{
			summ += B[k * m + k + i] * B[k * m + k + i];
		}
		len = std::sqrt(summ);
		tempV = B[k * m + k];
		alpha = -SIGN<float>(tempV) * len;
		gamma = std::sqrt(0.5 * (alpha * alpha - B[k * m + k] * alpha));
		thrust::fill(hvec.begin(), hvec.end(), 0.0f);
		hvec[k] = (B[k * m + k] - alpha) / (2.0 * gamma);
		for (size_t jj = k + 1; jj < m; jj++){ hvec[jj] = B[k * m + jj] / (2.0 * gamma); }

		// Update U
		CUDA_CHECK_RETURN(cudaSetDevice(1));
		hvec1 = hvec;
		//Calculate U
		cublasSgemv(cublashandle1, CUBLAS_OP_N, m, m, &ONE, thrust::raw_pointer_cast(&U[0]), m, thrust::raw_pointer_cast(&hvec1[0]), 1, &ZERO, thrust::raw_pointer_cast(&Uhvec1[0]), 1);
		updateU<float> << <mgid, mblk, 0, streams[1] >> >(thrust::raw_pointer_cast(&U[0]), thrust::raw_pointer_cast(&Uhvec1[0]), thrust::raw_pointer_cast(&hvec1[0]), m);

		// Calculate A' * hvec;
		CUDA_CHECK_RETURN(cudaSetDevice(0));
		cublasSgemv(cublashandle, CUBLAS_OP_T, m, n, &ONE, thrust::raw_pointer_cast(&B[0]), m, thrust::raw_pointer_cast(&hvec[0]), 1, &ZERO, thrust::raw_pointer_cast(&lvec[0]), 1);
		calculateQB<float> << <gid, blk, 0, streams[0] >> >(thrust::raw_pointer_cast(&B[0]), thrust::raw_pointer_cast(&hvec[0]), thrust::raw_pointer_cast(&lvec[0]), m, n);


		// Update U in another device
		if (k < n - 2)
		{
			CUDA_CHECK_RETURN(cudaSetDevice(0));
			summ = 0;
			for (size_t i = 0; i != (n - 1 - k); ++i){ summ += B[(k + 1 + i) * m + k] * B[(k + 1 + i) * m + k]; }
			len = std::sqrt(summ);
			tempV = B[(k + 1) * m + k];
			alpha = -SIGN(tempV) * len;
			gamma = std::sqrt(1.0 / 2.0 * (alpha * alpha - B[(k + 1) * m + k] * alpha));
			thrust::fill(lvec.begin(), lvec.end(), 0.0f);

			lvec[k + 1] = (B[(k + 1) * m + k] - alpha) / (2.0f * gamma);
			for (size_t jj = k + 2; jj < n; jj++){ lvec[jj] = B[jj * m + k] / (2.0f * gamma); }

			CUDA_CHECK_RETURN(cudaSetDevice(2));
			lvec2 = lvec;
			// Calculate V
			cublasSgemv(cublashandle2, CUBLAS_OP_N, n, n, &ONE, thrust::raw_pointer_cast(&V[0]), n, thrust::raw_pointer_cast(&lvec2[0]), 1, &ZERO, thrust::raw_pointer_cast(&Vlvec2[0]), 1);
			updateV<float> << <ngid, nblk, 0, streams[2] >> >(thrust::raw_pointer_cast(&V[0]), thrust::raw_pointer_cast(&Vlvec2[0]), thrust::raw_pointer_cast(&lvec2[0]), n);

			// Calculate B * P
			CUDA_CHECK_RETURN(cudaSetDevice(0));
			cublasSgemv(cublashandle, CUBLAS_OP_N, m, n, &ONE, thrust::raw_pointer_cast(&B[0]), m, thrust::raw_pointer_cast(&lvec[0]), 1, &ZERO, thrust::raw_pointer_cast(&hvec[0]), 1);
			calculateBP<float> << <gid, blk, 0, streams[0] >> >(thrust::raw_pointer_cast(&B[0]), thrust::raw_pointer_cast(&lvec[0]), thrust::raw_pointer_cast(&hvec[0]), m, n);

		}
		CUDA_CHECK_RETURN(cudaDeviceSynchronize());
	}

	CUDA_CHECK_RETURN(cudaSetDevice(0)); //Used to update Left part
	CUDA_CHECK_RETURN(cudaStreamDestroy(streams[0]));
	hvec.clear();
	lvec.clear();
	cublasDestroy(cublashandle);

	CUDA_CHECK_RETURN(cudaSetDevice(1)); //Used to update U
	CUDA_CHECK_RETURN(cudaStreamDestroy(streams[1]));
	hvec1.clear();
	Uhvec1.clear();
	U.clear();
	cublasDestroy(cublashandle1);

	CUDA_CHECK_RETURN(cudaSetDevice(2)); //Used to update V
	CUDA_CHECK_RETURN(cudaStreamDestroy(streams[2]));
	lvec2.clear();
	Vlvec2.clear();
	V.clear();
	cublasDestroy(cublashandle2);

}
void HouseHolder_v2(thrust::device_vector<double>& B, const int m, const int n)
{
	cudaStream_t streams[3];
	double summ(0);	double len(0); double tempV(0); double alpha(0); double gamma(0);
	const double ONE = 1.0;
	const double ZERO = 0.0;
	CUDA_CHECK_RETURN(cudaSetDevice(0));
	CUDA_CHECK_RETURN(cudaStreamCreate(&streams[0]));
	cublasHandle_t cublashandle;
	cublasCreate(&cublashandle);
	thrust::device_vector<double> hvec(m, 0);
	thrust::device_vector<double> lvec(n, 0);



	CUDA_CHECK_RETURN(cudaSetDevice(1)); //Used to update Left part
	CUDA_CHECK_RETURN(cudaStreamCreate(&streams[1]));
	thrust::device_vector<double> hvec1(m, 0);
	thrust::device_vector<double> Uhvec1(m, 0);
	thrust::device_vector<double> U(m*m, 0);
	for (size_t i = 0; i < m; i++)
	{
		U[i * m + i] = 1.0;
	}
	cublasHandle_t cublashandle1;
	cublasCreate(&cublashandle1);



	CUDA_CHECK_RETURN(cudaSetDevice(2)); //Used to update right part
	CUDA_CHECK_RETURN(cudaStreamCreate(&streams[2]));
	thrust::device_vector<double> lvec2(n, 0);
	thrust::device_vector<double> Vlvec2(n, 0);
	thrust::device_vector<double> V(n * n, 0);
	for (size_t i = 0; i != n; ++i)
	{
		V[i * n + i] = 1.0;
	}
	cublasHandle_t cublashandle2;
	cublasCreate(&cublashandle2);

	dim3 blk(16, 16);
	dim3 gid(
		(m + blk.x - 1) / blk.x,
		(n + blk.y - 1) / blk.y);

	dim3 mblk(16, 16);
	dim3 mgid(
		(m + mblk.x - 1) / mblk.x,
		(m + mblk.y - 1) / mblk.y);

	dim3 nblk(16, 16);
	dim3 ngid(
		(n + nblk.x - 1) / nblk.x,
		(n + nblk.y - 1) / nblk.y);


	for (size_t k = 0; k < n; k++)
	{
		std::cout << k << std::endl;
		CUDA_CHECK_RETURN(cudaSetDevice(0));
		//Copy Vector;
		summ = 0;
		for (size_t i = 0; i < m - k; i++)
		{
			summ += B[k * m + k + i] * B[k * m + k + i];
		}
		len = std::sqrt(summ);
		tempV = B[k * m + k];
		alpha = -SIGN<float>(tempV) * len;
		gamma = std::sqrt(0.5 * (alpha * alpha - B[k * m + k] * alpha));
		thrust::fill(hvec.begin(), hvec.end(), 0.0f);
		hvec[k] = (B[k * m + k] - alpha) / (2.0 * gamma);
		for (size_t jj = k + 1; jj < m; jj++){ hvec[jj] = B[k * m + jj] / (2.0 * gamma); }

		// Update U
		CUDA_CHECK_RETURN(cudaSetDevice(1));
		hvec1 = hvec;
		//Calculate U
		cublasDgemv(cublashandle1, CUBLAS_OP_N, m, m, &ONE, thrust::raw_pointer_cast(&U[0]), m, thrust::raw_pointer_cast(&hvec1[0]), 1, &ZERO, thrust::raw_pointer_cast(&Uhvec1[0]), 1);
		updateU<double> << <mgid, mblk, 0, streams[1] >> >(thrust::raw_pointer_cast(&U[0]), thrust::raw_pointer_cast(&Uhvec1[0]), thrust::raw_pointer_cast(&hvec1[0]), m);

		// Calculate A' * hvec;
		CUDA_CHECK_RETURN(cudaSetDevice(0));
		cublasDgemv(cublashandle, CUBLAS_OP_T, m, n, &ONE, thrust::raw_pointer_cast(&B[0]), m, thrust::raw_pointer_cast(&hvec[0]), 1, &ZERO, thrust::raw_pointer_cast(&lvec[0]), 1);
		calculateQB<double> << <gid, blk, 0, streams[0] >> >(thrust::raw_pointer_cast(&B[0]), thrust::raw_pointer_cast(&hvec[0]), thrust::raw_pointer_cast(&lvec[0]), m, n);


		// Update U in another device
		if (k < n - 2)
		{
			CUDA_CHECK_RETURN(cudaSetDevice(0));
			summ = 0;
			for (size_t i = 0; i != (n - 1 - k); ++i){ summ += B[(k + 1 + i) * m + k] * B[(k + 1 + i) * m + k]; }
			len = std::sqrt(summ);
			tempV = B[(k + 1) * m + k];
			alpha = -SIGN(tempV) * len;
			gamma = std::sqrt(1.0 / 2.0 * (alpha * alpha - B[(k + 1) * m + k] * alpha));
			thrust::fill(lvec.begin(), lvec.end(), 0.0f);

			lvec[k + 1] = (B[(k + 1) * m + k] - alpha) / (2.0f * gamma);
			for (size_t jj = k + 2; jj < n; jj++){ lvec[jj] = B[jj * m + k] / (2.0f * gamma); }

			CUDA_CHECK_RETURN(cudaSetDevice(2));
			lvec2 = lvec;
			// Calculate V
			cublasDgemv(cublashandle2, CUBLAS_OP_N, n, n, &ONE, thrust::raw_pointer_cast(&V[0]), n, thrust::raw_pointer_cast(&lvec2[0]), 1, &ZERO, thrust::raw_pointer_cast(&Vlvec2[0]), 1);
			updateV<double> << <ngid, nblk, 0, streams[2] >> >(thrust::raw_pointer_cast(&V[0]), thrust::raw_pointer_cast(&Vlvec2[0]), thrust::raw_pointer_cast(&lvec2[0]), n);

			// Calculate B * P
			CUDA_CHECK_RETURN(cudaSetDevice(0));
			cublasDgemv(cublashandle, CUBLAS_OP_N, m, n, &ONE, thrust::raw_pointer_cast(&B[0]), m, thrust::raw_pointer_cast(&lvec[0]), 1, &ZERO, thrust::raw_pointer_cast(&hvec[0]), 1);
			calculateBP<double> << <gid, blk, 0, streams[0] >> >(thrust::raw_pointer_cast(&B[0]), thrust::raw_pointer_cast(&lvec[0]), thrust::raw_pointer_cast(&hvec[0]), m, n);

		}
		CUDA_CHECK_RETURN(cudaDeviceSynchronize());
	}

	CUDA_CHECK_RETURN(cudaSetDevice(0)); //Used to update Left part
	CUDA_CHECK_RETURN(cudaStreamDestroy(streams[0]));
	hvec.clear();
	lvec.clear();
	cublasDestroy(cublashandle);

	CUDA_CHECK_RETURN(cudaSetDevice(1)); //Used to update U
	CUDA_CHECK_RETURN(cudaStreamDestroy(streams[1]));
	hvec1.clear();
	Uhvec1.clear();
	U.clear();
	cublasDestroy(cublashandle1);

	CUDA_CHECK_RETURN(cudaSetDevice(2)); //Used to update V
	CUDA_CHECK_RETURN(cudaStreamDestroy(streams[2]));
	lvec2.clear();
	Vlvec2.clear();
	V.clear();
	cublasDestroy(cublashandle2);

}

void threeGPUHouseHolder(const std::string& FileName,const std::string& dataType, const int m, const int n)
{
	std::ifstream fin(FileName.c_str(),std::ios::binary);
	if(!fin.is_open())
	{
		std::cout<<"Cannot open the file for Householder transform\n";
		exit(-1);
	}
	if(dataType == "float")
	{
		std::vector<float> hB(m * n, 0);
		fin.read((char*)(&hB[0]), sizeof(float) * m * n);
		fin.close();
		thrust::device_vector<float> B = hB;
		HouseHolder_v2(B,m,n);
	}
	else if(dataType == "double")
	{
		std::vector<double> hB(m * n, 0);
		fin.read((char*)(&hB[0]), sizeof(double) * m * n);
		fin.close();
		thrust::device_vector<double> B = hB;
		HouseHolder_v2(B,m,n);
	}
	
}

int testHouseHolder()
{
	const int m = 22000;
	const int n = 22000;
	std::vector<float> hA(m * n, 0);
	for (size_t i = 0; i < m * n; i++)
	{
		hA[i] = rand() % (m * n);
	}
	CUDA_CHECK_RETURN(cudaDeviceReset());
	CUDA_CHECK_RETURN(cudaSetDevice(0));
	thrust::device_vector<float> A = hA;
	time_t start = clock();
	HouseHolder_v2(A, m, n);
	double duration = double((double(clock()) - double(start))) / CLOCKS_PER_SEC;
	CUDA_CHECK_RETURN(cudaDeviceReset());
	std::cout << duration << std::endl;
	return 0;

}



}
