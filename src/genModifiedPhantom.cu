#include "CTLAB.h"
#include "cudaCheckReturner.h"
#include <string>
#include <iostream>
#include <fstream>
template<typename T>
inline __host__ __device__ T ComGray(T x, T a, T b, int m, int n)
{
	T f;
	if (x >= 1)
		f = 0;
	else
		f = 1 + a * pow(x, m) + b * pow(x, n);
	return f;

}

//! Generate 3D Phantom
//! version 1.0
template<typename T>
__global__ void genModiPhKer(
	T* ptr_data,
	const T* PhPar,
	const T* AxisSquare,
	const T hfx,
	const T hfy,
	const T hfz,
	const T dx,
	const T dy,
	const T dz,
	const unsigned int lR,
	const unsigned int wR,
	const unsigned int hR)
{
	const unsigned int xIdx = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned int yIdx = blockDim.y * blockIdx.y + threadIdx.y;
	const unsigned int zIdx = blockDim.z * blockIdx.z + threadIdx.z;
	if (xIdx < lR&&yIdx < wR&&zIdx < hR)
	{
		T cx, cy, cz;
		cx = (xIdx - hfx)*dx;
		cy = (yIdx - hfy)*dy;
		cz = (zIdx - hfz)*dz;
		T tpdata = 0;
#pragma unroll
		for (unsigned int ElpIndex = 0; ElpIndex < 10; ElpIndex++)
		{
			T tptheta = (PhPar[ElpIndex * 11 + 6]) * 0.01745329251994329576923690634218;
			T RX1 = cx - PhPar[ElpIndex * 11 + 3];
			T RX2 = cy - PhPar[ElpIndex * 11 + 4];
			T RX3 = cz - PhPar[ElpIndex * 11 + 5];
			T tempvar = RX1;
			RX1 = tempvar * cos(tptheta) + RX2 * sin(tptheta);
			RX2 = -tempvar * sin(tptheta) + RX2 * cos(tptheta);

			T rsquare =
				pow(RX1, 2) / AxisSquare[ElpIndex * 3 + 0] +
				pow(RX2, 2) / AxisSquare[ElpIndex * 3 + 1] +
				pow(RX3, 2) / AxisSquare[ElpIndex * 3 + 2];
			if (rsquare < 1)
			{
				T a = -PhPar[ElpIndex * 11 + 9] / (PhPar[ElpIndex * 11 + 9] - PhPar[ElpIndex * 11 + 8]);
				T b = PhPar[ElpIndex * 11 + 8] / (PhPar[ElpIndex * 11 + 9] - PhPar[ElpIndex * 11 + 8]);
				T gs = pow(static_cast<double>(PhPar[ElpIndex * 11 + 10]), static_cast<double>(2.0));
				T GSC = pow((double) gs, static_cast<double>(PhPar[ElpIndex * 11 + 8]));
				T gc = 1.0 / (1 - GSC);


				tempvar = ComGray<T>(rsquare, a, b, PhPar[ElpIndex * 11 + 8], PhPar[ElpIndex * 11 + 9]);

				tempvar = tempvar - GSC * ComGray<T>(rsquare / gs, a, b, PhPar[ElpIndex * 11 + 8], PhPar[ElpIndex * 11 + 9]);
				tpdata = tpdata + tempvar * gc * PhPar[ElpIndex * 11 + 7];
			}
		}
		ptr_data[(xIdx * wR + yIdx) * hR + zIdx] = tpdata;
	}
}



template<typename T>
void generateModiPhantom_template(
	const std::string& FileName,
	const unsigned int& lenReso = 256,
	const unsigned int& widReso = 256,
	const unsigned int& heiReso = 256)
{
	const T PhPar[110] = {
		0.6900f, 0.900f, 0.900f, 0.0f, 0.0f, 0.0f, 0.0f, 2.0f, 20.0f, 40.0f, 0.98f,
		0.6792f, 0.882f, 0.882f, 0.0f, 0.0f, 0.0f, 0.0f, -0.98f, 20.0f, 40.0f, 0.98f,
		0.4100f, 0.160f, 0.210f, -0.22f, 0.0f, -0.25f, 108.0f, -0.02f, 3.0f, 6.0f, 0.20f,
		0.3100f, 0.110f, 0.220f, 0.22f, 0.0f, -0.25f, 72.0f, -0.02f, 3.0f, 6.0f, 0.20f,
		0.2100f, 0.250f, 0.500f, 0.0f, 0.35f, -0.25f, 0.0f, 0.02f, 3.0f, 6.0f, 0.20f,
		0.0460f, 0.046f, 0.046f, 0.0f, 0.1f, -0.25f, 0.0f, 0.02f, 2.0f, 4.0f, 0.10f,
		0.0460f, 0.023f, 0.020f, -0.08f, -0.65f, -0.25f, 0.0f, 0.01f, 2.0f, 4.0f, 0.10f,
		0.0460f, 0.023f, 0.020f, 0.06f, -0.65f, -0.25f, 90.0f, 0.01f, 2.0f, 4.0f, 0.10f,
		0.0560f, 0.040f, 0.100f, 0.06f, -0.105f, 0.625f, 90.0f, 0.02f, 2.0f, 4.0f, 0.20f,
		0.0560f, 0.056f, 0.100f, 0.0f, 0.100f, 0.625f, 0.0f, -0.02f, 2.0f, 4.0f, 0.20f };

	T *d_PhPar;
	T *d_phantom;
	T *phantom;

	T *axis = new T[30];
	T *d_axis;
	const unsigned int totSize = lenReso * widReso * heiReso;
	try	{
		phantom = new T[totSize];
	} catch (std::bad_alloc) {
		// Todo: try to logging it.
		std::cout << "Your requirement is too large!\n";
		exit(-1);
	}
	memset(phantom, 0, totSize);

	for (unsigned int ElpIndex = 0; ElpIndex != 10; ++ElpIndex)	{
		axis[ElpIndex * 3 + 0] = pow(PhPar[ElpIndex * 11 + 0], 2);
		axis[ElpIndex * 3 + 1] = pow(PhPar[ElpIndex * 11 + 1], 2);
		axis[ElpIndex * 3 + 2] = pow(PhPar[ElpIndex * 11 + 2], 2);
	}

	dim3 blockSize(8, 8, 8);
	dim3 gridSize((lenReso + blockSize.x - 1) / blockSize.x, (widReso + blockSize.y - 1) / blockSize.y, (heiReso + blockSize.z - 1) / blockSize.z);

	
	CUDA_CHECK_RETURN(cudaMalloc((void**) &d_PhPar, sizeof(T) * 110));
	CUDA_CHECK_RETURN(cudaMalloc((void**) &d_phantom, sizeof(T)*totSize));
	CUDA_CHECK_RETURN(cudaMalloc((void**) &d_axis, sizeof(T) * 30));
	CUDA_CHECK_RETURN(cudaMemcpy(d_axis, axis, sizeof(T) * 30, cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(cudaMemset(d_phantom, 0, sizeof(T)*totSize));
	CUDA_CHECK_RETURN(cudaMemcpy(d_PhPar, PhPar, sizeof(T) * 110, cudaMemcpyHostToDevice));
	
	T halfx = (lenReso - 1)*0.5f;
	T halfy = (widReso - 1)*0.5f;
	T halfz = (widReso - 1)*0.5f;
	T deltax = 2.0f / lenReso;
	T deltay = 2.0f / widReso;
	T deltaz = 2.0f / heiReso;

	genModiPhKer<T> << <gridSize, blockSize >> >(d_phantom, d_PhPar, d_axis, halfx, halfy, halfz, deltax, deltay, deltaz, lenReso, widReso, heiReso);

	CUDA_CHECK_RETURN(cudaMemcpy(phantom, d_phantom, sizeof(T)*totSize, cudaMemcpyDeviceToHost));
	std::ofstream fout(FileName.c_str(), std::ios::binary);
	fout.write((char*) phantom, sizeof(T)*lenReso*widReso*heiReso);
	fout.close();
}

void generateModiPhantomFloat(const std::string& FileName, const unsigned int& lenReso, const unsigned int& widReso, const unsigned int& heiReso) {
	generateModiPhantom_template<float>(FileName, lenReso, widReso, heiReso);
}
void generateModiPhantomDouble(const std::string& FileName, const unsigned int& lenReso, const unsigned int& widReso, const unsigned int& heiReso) {
	generateModiPhantom_template<double>(FileName, lenReso, widReso, heiReso);
}

