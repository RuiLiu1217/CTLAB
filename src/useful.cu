/*
 * COPYRIGHT NOTICE
 * COPYRIGHT (c) 2015, Wake Forest and UMass Lowell
 * All rights reserved
 *
 * @file FastMatrixVectorMultiplication.cu
 * @brief The implementation GPU based fast matrix and vector multiplication
 *
 * @version 1.0
 * @author Rui Liu
 * @date May. 1, 2015
 *
 */
#include "useful.hpp"
#include "utilities.hpp"
#include "projbackproj.hpp"

#ifndef MYEPSILON
#define MYEPSILON 1.0E-9
#endif

// y = Ax
// A : m-by-n matrix, x : n elements vector, y : m elements vector
// m and n are arbitrary positive integers.
texture<float4, 2, cudaReadModeElementType> texRefA;
#define bx blockIdx.x
#define tx threadIdx.x
#define ty threadIdx.y
__global__ void FMVM_KER(float* y, cudaArray* A, float* x, int m, int n)
{
	__shared__ float xs[16][16];
	__shared__ float Ps[16][16];
	float4 a;
	float* Psptr = (float*) Ps + (ty << 4) + tx;
	int ay = (bx << 4) + ty;
	float *xptr = x + (ty << 4) + tx;
	float *xsptr = (float *) xs + (tx << 2);
	*Psptr = 0.0f;

	int i;
	for (i = 0; i < (n & ~255); i += 256, xptr += 256)
	{
		xs[ty][tx] = *xptr;
		__syncthreads();
		int ax = tx + (i >> 2);
		a = tex2D(texRefA, ax, ay);
		*Psptr += a.x * *xsptr + a.y * *(xsptr + 1) + a.z * *(xsptr + 2) + a.w * *(xsptr + 3);
		a = tex2D(texRefA, ax + 16, ay);
		*Psptr += a.x * *(xsptr + 64) + a.y * *(xsptr + 65) + a.z * *(xsptr + 66) + a.w * *(xsptr + 67);
		a = tex2D(texRefA, ax + 32, ay);
		*Psptr += a.x * *(xsptr + 128) + a.y * *(xsptr + 129) + a.z * *(xsptr + 130) + a.w * *(xsptr + 131);
		a = tex2D(texRefA, ax + 48, ay);
		*Psptr += a.x * *(xsptr + 192) + a.y * *(xsptr + 193) + a.z * *(xsptr + 194) + a.w * *(xsptr + 195);
		__syncthreads();
	}
	if (i + (ty << 4) + tx < n)
	{
		xs[ty][tx] = *xptr;
	}
	__syncthreads();
	int j;
	for (j = 0; j < ((n - i) >> 6); j++, xsptr += 61)
	{
		a = tex2D(texRefA, tx + (i >> 2) + (j << 4), ay);
		*Psptr += a.x * *xsptr++ + a.y * *xsptr++ + a.z * *xsptr++ + a.w * *xsptr;
	}
	__syncthreads();
	int remain = (n - i) & 63;
	if ((tx << 2) < remain)
	{
		a = tex2D(texRefA, tx + (i >> 2) + (j << 4), ay);
		*Psptr += a.x * *xsptr++;
	}
	if ((tx << 2) + 1 < remain) *Psptr += a.y * *xsptr++;
	if ((tx << 2) + 2 < remain) *Psptr += a.z * *xsptr++;
	if ((tx << 2) + 3 < remain) *Psptr += a.w * *xsptr;
	__syncthreads();

	if (tx < 8) *Psptr += *(Psptr + 8);
	if (tx < 4) *Psptr += *(Psptr + 4);
	if (tx < 2) *Psptr += *(Psptr + 2);
	if (tx < 1) *Psptr += *(Psptr + 1);
	__syncthreads();
	if (ty == 0 && (bx << 4) + tx < m) y[(bx << 4) + tx] = Ps[tx][0];

}


// Fast Matrix-Vector Multiplication
// A store in row
// A = [1 2 3;
//      4 5 6] m = 2, n = 3; A[1] = 2;
void FMVM(float* y, float* A, float* x, int m, int n)
{
	int blkNum((m >> 4) + ((m & 15) ? 1 : 0));
	int height(blkNum << 4);
	int width((n & 255) ? (256 * ((n >> 8) + 1)) : n);
	dim3 threads(16, 16), grid(blkNum, 1);
	cudaArray *d_A;
	float* d_x;
	float *d_y;

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float4>();
	checkCudaErrors(cudaMallocArray(&d_A, &channelDesc, width >> 2, height));
	checkCudaErrors(cudaMemcpy2DToArray(d_A, 0, 0, A, n * sizeof(float), n * sizeof(float), m, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaBindTextureToArray(texRefA, d_A));
	checkCudaErrors(cudaMalloc((void**) &d_x, n * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**) &d_y, m * sizeof(float)));

	checkCudaErrors(cudaMemcpy(d_x, x, n * sizeof(float), cudaMemcpyHostToDevice));
	FMVM_KER << < grid, threads >> >(d_y, d_A, d_x, m, n);
	checkCudaErrors(cudaMemcpy(y, d_y, m * sizeof(float), cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaFree(d_y));
	checkCudaErrors(cudaFree(d_x));
	checkCudaErrors(cudaUnbindTexture(texRefA));
	checkCudaErrors(cudaFreeArray(d_A));

}

void FMVM_device(float* d_y, float* dA, float* d_x, int m, int n)
{
	int blkNum((m >> 4) + ((m & 15) ? 1 : 0));
	int height(blkNum << 4);
	int width((n & 255) ? (256 * ((n >> 8) + 1)) : n);
	dim3 threads(16, 16), grid(blkNum, 1);
	cudaArray *d_A;

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float4>();
	checkCudaErrors(cudaMallocArray(&d_A, &channelDesc, width >> 2, height));
	checkCudaErrors(cudaMemcpy2DToArray(d_A, 0, 0, dA, n * sizeof(float), n * sizeof(float), m, cudaMemcpyDeviceToDevice));
	checkCudaErrors(cudaBindTextureToArray(texRefA, d_A));

	FMVM_KER << < grid, threads >> >(d_y, d_A, d_x, m, n);
	checkCudaErrors(cudaUnbindTexture(texRefA));
	checkCudaErrors(cudaFreeArray(d_A));
}
void testFMVM()
{
	float* y = new float[3];
	float* A = new float[9];
	float* x = new float[3];
	A[0] = 1;
	A[1] = 2;
	A[2] = 3;
	A[3] = 4;
	A[4] = 5;
	A[5] = 6;
	A[6] = 7;
	A[7] = 8;
	A[8] = 9;
	x[0] = 1;
	x[1] = 2;
	x[2] = 3;
	float* d_x, *d_y, *d_A;
	checkCudaErrors(cudaMalloc((void**) &d_x, sizeof(float) * 3));
	checkCudaErrors(cudaMalloc((void**) &d_y, sizeof(float) * 3));
	checkCudaErrors(cudaMalloc((void**) &d_A, sizeof(float) * 9));
	checkCudaErrors(cudaMemcpy(d_x, x, sizeof(float) * 3, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_y, y, sizeof(float) * 3, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_A, A, sizeof(float) * 9, cudaMemcpyHostToDevice));

	FMVM_device(d_y, d_A, d_x, 3, 3);
	checkCudaErrors(cudaMemcpy(y, d_y, sizeof(float) * 3, cudaMemcpyDeviceToHost));
	std::cout << y[0] << " " << y[1] << " " << y[2];
}







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
	try
	{
		phantom = new T[totSize];
	}
	catch (std::bad_alloc& a)
	{
		std::cerr << "Bad_alloc caught: "<<a.what()<<"\n";
		std::cerr << "Your requirement is too large!\n";
		exit(-1);
	}
	memset(phantom, 0, totSize);
	for (unsigned int ElpIndex = 0; ElpIndex != 10; ++ElpIndex)
	{
		axis[ElpIndex * 3 + 0] = pow(PhPar[ElpIndex * 11 + 0], 2);
		axis[ElpIndex * 3 + 1] = pow(PhPar[ElpIndex * 11 + 1], 2);
		axis[ElpIndex * 3 + 2] = pow(PhPar[ElpIndex * 11 + 2], 2);
	}

	dim3 blockSize(8, 8, 8);
	dim3 gridSize((lenReso + blockSize.x - 1) / blockSize.x, (widReso + blockSize.y - 1) / blockSize.y, (heiReso + blockSize.z - 1) / blockSize.z);


	checkCudaErrors(cudaMalloc((void**) &d_PhPar, sizeof(T) * 110));
	checkCudaErrors(cudaMalloc((void**) &d_phantom, sizeof(T) * totSize));
	checkCudaErrors(cudaMalloc((void**) &d_axis, sizeof(T) * 30));
	checkCudaErrors(cudaMemcpy(d_axis, axis, sizeof(T) * 30, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemset(d_phantom, 0, sizeof(T)*totSize));
	checkCudaErrors(cudaMemcpy(d_PhPar, PhPar, sizeof(T) * 110, cudaMemcpyHostToDevice));

	T halfx = (lenReso - 1) * 0.5f;
	T halfy = (widReso - 1) * 0.5f;
	T halfz = (widReso - 1) * 0.5f;
	T deltax = 2.0f / lenReso;
	T deltay = 2.0f / widReso;
	T deltaz = 2.0f / heiReso;

	genModiPhKer<T> << <gridSize, blockSize >> >(d_phantom, d_PhPar, d_axis, halfx, halfy, halfz, deltax, deltay, deltaz, lenReso, widReso, heiReso);

	checkCudaErrors(cudaMemcpy(phantom, d_phantom, sizeof(T)*totSize, cudaMemcpyDeviceToHost));
	std::ofstream fout("out.raw", std::ios::binary);
	fout.write((char*) phantom, sizeof(T)*lenReso*widReso*heiReso);
	fout.close();
}

void generateModiPhantom(const std::string& FileName, const unsigned int& lenReso, const unsigned int& widReso, const unsigned int& heiReso)
{
	generateModiPhantom_template<float>(FileName, lenReso, widReso, heiReso);
}
void generateModiPhantom_d(const std::string& FileName, const unsigned int& lenReso, const unsigned int& widReso, const unsigned int& heiReso)
{
	generateModiPhantom_template<double>(FileName, lenReso, widReso, heiReso);
}



float norm(float* v, cuint len, const float p)
{
	thrust::device_ptr<float> pv(v);
	return powf(thrust::transform_reduce(pv, pv + len, _lp_functor<float>(p), 0.0f, thrust::plus<float>()), 1.0f / p);
}

template<typename T>
T norm(thrust::device_vector<T>& v, const T& p)
{
	return powf(thrust::transform_reduce(v.begin(), v.end(), _lp_functor<T>(p), 0.0f, thrust::plus<T>()), 1.0f / p);
}

template<typename T>T innerProd_template(T* a, T* b, cuint L)
{
	thrust::device_ptr<T> pa(a);
	thrust::device_ptr<T> pb(b);
	return thrust::inner_product(pa, pa + L, pb, 0.0f);
}
float innerProd(float* a, float* b, cuint L)
{
	return innerProd_template<float>(a, b, L);
}
double innerProd(double* a, double* b, cuint L)
{
	return innerProd_template<double>(a, b, L);
}





template<typename T>
__global__ void _GradientOfTV_ker(T* f, T* d, const T coef, cuint L, cuint W)
{
	cuint idx = threadIdx.x + blockDim.x * blockIdx.x;
	cuint idy = threadIdx.y + blockDim.y * blockIdx.y;
	if (idx < L && idy < W)
	{
		cuint curid = idy * L + idx;
		const T fij = f[curid];
		const T fi1j = f[idy * L + (idx + 1) % L];
		const T fi_1j = f[idy * L + (idx + L - 1) % L];
		const T fij1 = f[((idy + 1) % W) * L + idx];
		const T fij_1 = f[((idy + W - 1) % W) * L + idx];
		const T fi_1j1 = f[((idy + 1) % W) * L + (idx + L - 1) % L];
		const T fi1j_1 = f[((idy + W - 1) % W) * L + (idx + 1) % L];
		const T dom1 = 1.0f / (sqrtf((fi1j - fij)*(fi1j - fij) + (fij1 - fij) * (fij1 - fij) + MYEPSILON));
		const T dom2 = 1.0f / (sqrtf((fij - fi_1j)*(fij - fi_1j) + (fi_1j1 - fi_1j)*(fi_1j1 - fi_1j) + MYEPSILON));
		const T dom3 = 1.0f / (sqrtf((fi1j_1 - fij_1) * (fi1j_1 - fij_1) + (fij - fij_1)*(fij - fij_1) + MYEPSILON));
		d[curid] = ((2.0 * fij - fi1j - fij1) * dom1 + (fij - fi_1j) * dom2 + (fij - fij_1)*dom3) * coef;

		return;
	}
}

void GradientOfTV(float* f, float* d, const float coef, cuint L, cuint W, const dim3& blk, const dim3& gid)
{
	_GradientOfTV_ker<float> << <gid, blk >> >(f, d, coef, L, W);
}

void GradientOfTV(double* f, double* d, const double coef, cuint L, cuint W, const dim3& blk, const dim3& gid)
{
	_GradientOfTV_ker<double> << <gid, blk >> >(f, d, coef, L, W);
}


template<typename T>
__global__ void _GradientOfTV_SHARED(T* img, T* des, const T coef,
	const unsigned int imgL, const unsigned int imgW)
{
	int i = threadIdx.x + blockIdx.x * (blockDim.x - 2);
	int j = threadIdx.y + blockIdx.y * (blockDim.y - 2);
	if (i < imgL && j < imgW)
	{
		int idx = j * imgL + i;

		__shared__ T simg[32][32];
		simg[threadIdx.y][threadIdx.x] = img[idx];
		__syncthreads();
		if (threadIdx.x == 31 || threadIdx.y == 31 || threadIdx.x == 0 || threadIdx.y == 0)
			return;


		const unsigned int y = threadIdx.y;
		const unsigned int x = threadIdx.x;
		T val =
			((2.0f * simg[y][x] - simg[y][x - 1] - simg[y - 1][x]) / sqrtf(1e-36 + (simg[y][x] - simg[y][x - 1])*(simg[y][x] - simg[y][x - 1]) + (simg[y][x] - simg[y - 1][x])*(simg[y][x] - simg[y - 1][x])))
			- ((simg[y][x + 1] - simg[y][x]) / sqrtf(1e-36 + (simg[y][x + 1] - simg[y][x]) * (simg[y][x + 1] - simg[y][x]) + (simg[y][x + 1] - simg[y - 1][x + 1]) * (simg[y][x + 1] - simg[y - 1][x + 1])))
			- ((simg[y + 1][x] - simg[y][x]) / sqrtf(1e-36 + (simg[y + 1][x] - simg[y + 1][x - 1]) * (simg[y + 1][x] - simg[y + 1][x - 1]) + (simg[y + 1][x] - simg[y][x])*(simg[y + 1][x] - simg[y][x])));
		if (i == imgL - 1 || j == imgW - 1 || i == 0 || j == 0)
			val = 0;
		des[idx] = (val * coef);
	}
}

//
//dim3 blk(32, 32);
//dim3 gid(
//	(L + 29) / 30,
//	(W + 29) / 30);
void GradientOfTV_SHARED(float* f, float* d, const float coef, cuint L, cuint W, const dim3& blk, const dim3& gid)
{
	_GradientOfTV_SHARED<float> << <gid, blk >> >(f, d, coef, L, W);
}

//dim3 blk(32, 32);
//dim3 gid(
//	(L + 29) / 30,
//	(W + 29) / 30);
void GradientOfTV_SHARED(double* f, double* d, const double coef, cuint L, cuint W, const dim3& blk, const dim3& gid)
{
	_GradientOfTV_SHARED<double> << <gid, blk >> >(f, d, coef, L, W);
}





//ŒÆËãÍŒÏñµÄÌÝ¶È±ä»»ÍŒÏñ;
template<typename T>
__global__ void _DiscreteGradientTrans_ker(T* f, T* d, const T coef, cuint L, cuint W)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	if (idx < L && idy < W)
	{
		cuint curid = idy * L + idx;
		cuint nexxd = idy * L + ((idx + 1) % L);
		cuint nexyd = ((idy + 1) % W) * L + idx;
		const T difx = f[nexxd] - f[curid];
		const T dify = f[nexyd] - f[curid];
		d[curid] = sqrtf(difx * difx + dify * dify) * coef;
		return;
	}
}

template<typename T>
__global__ void _DiscreteGradientTrans_ker(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	cuint idz = threadIdx.z + blockIdx.z * blockDim.z;
	if (idx < L && idy < W && idz < H)
	{
		cuint curid = (idz * W + idy) * L + idx;
		cuint nexxd = (idz * W + idy) * L + (idx + 1) % L;
		cuint nexyd = (idz * W + ((idy + 1) % W)) * L + idx;
		cuint nexzd = (((idz + 1) % H) * W + idy) * L + idx;
		const T difx = f[nexxd] - f[curid];
		const T dify = f[nexyd] - f[curid];
		const T difz = f[nexzd] - f[curid];
		d[curid] = sqrtf(difx * difx + dify * dify + difz * difz) * coef;
		return;
	}
}

void DiscreteGradientTrans(float* f, float* d, const float coef, cuint L, cuint W, const dim3& blk, const dim3& gid)
{
	_DiscreteGradientTrans_ker<float> << <gid, blk >> >(f, d, coef, L, W);
}

void DiscreteGradientTrans(double* f, double* d, const double coef, cuint L, cuint W, const dim3& blk, const dim3& gid)
{
	_DiscreteGradientTrans_ker<double> << <gid, blk >> >(f, d, coef, L, W);
}


void DiscreteGradientTrans(float* f, float* d, const float coef, cuint L, cuint W, cuint H, const dim3& blk, const dim3& gid)
{
	_DiscreteGradientTrans_ker << <gid, blk >> >(f, d, coef, L, W, H);
}
void DiscreteGradientTrans(double* f, double* d, const double coef, cuint L, cuint W, cuint H, const dim3& blk, const dim3& gid)
{
	_DiscreteGradientTrans_ker << <gid, blk >> >(f, d, coef, L, W, H);
}



//¶þÎ¬ÍŒÏñ  Éè32x32µÄblock
template<typename T>
__global__ void _DiscreteGradientTrans_SHARED_Ker(T* img, T* tvimg, T coef, cuint imgL, cuint imgW)
{

	//µ±Ç°µÄÔ­ÊŒÍŒÏñË÷Òý;
	cuint idxX = threadIdx.x + blockIdx.x * (blockDim.x - 1);
	cuint idxY = threadIdx.y + blockIdx.y * (blockDim.y - 1);
	if (idxX < imgL && idxY < imgW)
	{
		cuint idx = idxY * imgL + idxX;
		T res(0);
		__shared__ T simg[32][32];
		simg[threadIdx.y][threadIdx.x] = img[idx];
		__syncthreads();
		if (threadIdx.x == 31 || threadIdx.y == 31)
			return;
		res = sqrtf(
			(simg[threadIdx.y][threadIdx.x] - simg[threadIdx.y][threadIdx.x + 1]) * (simg[threadIdx.y][threadIdx.x] - simg[threadIdx.y][threadIdx.x + 1]) +
			(simg[threadIdx.y][threadIdx.x] - simg[threadIdx.y + 1][threadIdx.x]) * (simg[threadIdx.y][threadIdx.x] - simg[threadIdx.y + 1][threadIdx.x]));
		tvimg[idx] = res * coef;
		if (idxX == imgL - 1 || idxY == imgW - 1 || idxX == 1 || idxY == 1)
			tvimg[idx] = 0;
	}
}

//ÈýÎ¬ÍŒÏñ Éè 8x8x8 block
template<typename T>
__global__ void _DiscreteGradientTrans_SHARED_Ker(T* vol, T* tvvol, T coef, cuint imgL, cuint imgW, cuint imgH)
{
	cuint idxX = threadIdx.x + blockIdx.x * (blockDim.x - 1);
	cuint idxY = threadIdx.y + blockIdx.y * (blockDim.y - 1);
	cuint idxZ = threadIdx.z + blockIdx.z * (blockDim.z - 1);
	if (idxX < imgL && idxY < imgW && idxZ < imgH)
	{
		cuint idx = (idxZ * imgW + idxY) * imgL + idxX;
		T res(0);
		__shared__ T svol[8][8][8];
		svol[threadIdx.z][threadIdx.y][threadIdx.x] = vol[idx];
		__syncthreads();
		if (threadIdx.x == 7 || threadIdx.y == 7 || threadIdx.z == 7)
			return;
		res = sqrtf(
			(svol[threadIdx.z][threadIdx.y][threadIdx.x] - svol[threadIdx.z][threadIdx.y][threadIdx.x + 1]) * (svol[threadIdx.z][threadIdx.y][threadIdx.x] - svol[threadIdx.z][threadIdx.y][threadIdx.x + 1]) +
			(svol[threadIdx.z][threadIdx.y][threadIdx.x] - svol[threadIdx.z][threadIdx.y + 1][threadIdx.x]) * (svol[threadIdx.z][threadIdx.y][threadIdx.x] - svol[threadIdx.z][threadIdx.y + 1][threadIdx.x]) +
			(svol[threadIdx.z][threadIdx.y][threadIdx.x] - svol[threadIdx.z + 1][threadIdx.y][threadIdx.x]) * (svol[threadIdx.z][threadIdx.y][threadIdx.x] - svol[threadIdx.z + 1][threadIdx.y][threadIdx.x]));
		tvvol[idx] = res * coef;
		if (idxX == imgL - 1 || idxY == imgW - 1 || idxZ == imgH - 1 || idxX == 1 || idxY == 1 || idxZ == 1)
			tvvol[idx] = 0;
	}
}

void DiscreteGradientTrans_SHARED(float* f, float* d, const float coef, cuint L, cuint W)
{
	dim3 blk(32, 32);
	//dim3 gid((L + 30) / 31,	(W + 30) / 31);
	dim3 gid((L + 29) / 30, (W + 29) / 30);

	_DiscreteGradientTrans_SHARED_Ker<float> << <gid, blk >> >(f, d, coef, L, W);
}
void DiscreteGradientTrans_SHARED(double* f, double* d, const double coef, cuint L, cuint W)
{
	dim3 blk(32, 32);
	//dim3 gid((L + 30) / 31,	(W + 30) / 31);
	dim3 gid((L + 29) / 30, (W + 29) / 30);

	_DiscreteGradientTrans_SHARED_Ker<double> << <gid, blk >> >(f, d, coef, L, W);
}
void DiscreteGradientTrans_SHARED(float* f, float* d, const float coef, cuint L, cuint W, cuint H)
{
	dim3 blk(8, 8, 8);
	dim3 gid(
		(L + 6) / 7,
		(W + 6) / 7,
		(H + 6) / 7);
	_DiscreteGradientTrans_SHARED_Ker<float> << <gid, blk >> >(f, d, coef, L, W, H);
}
void DiscreteGradientTrans_SHARED(double* f, double* d, const double coef, cuint L, cuint W, cuint H)
{
	dim3 blk(8, 8, 8);
	dim3 gid(
		(L + 6) / 7,
		(W + 6) / 7,
		(H + 6) / 7);
	_DiscreteGradientTrans_SHARED_Ker<double> << <gid, blk >> >(f, d, coef, L, W, H);
}




template<typename T>
__global__ void _invDiscreteGradientTransform_ker(T* f, T* d, T* r, const T omega, cuint L, cuint W)
{
	cuint idx = threadIdx.x + blockDim.x * blockIdx.x;
	cuint idy = threadIdx.y + blockDim.y * blockIdx.y;
	if (idx < L && idy < W)
	{
		const unsigned int curid = idy * L + idx;
		const unsigned int lasxd = idy * L + (idx + L - 1) % L;
		const unsigned int nexxd = idy * L + (idx + 1) % L;
		const unsigned int lasyd = ((idy + W - 1) % W) * L + idx;
		const unsigned int nexyd = ((idy + 1) % W) * L + idx;
		T fa(0), fb(0), fc(0);
		if (d[curid] < omega)
		{
			fa = (2.0 * f[curid] + f[nexxd] + f[nexyd]) * 0.25;
		}
		else
		{
			fa = f[curid] - omega * (2.0 * f[curid] - f[nexxd] - f[nexyd]) / (4 * d[curid]);
		}
		if (d[lasxd] < omega)
		{
			fb = (f[curid] + f[lasxd]) * 0.5;
		}
		else
		{
			fb = f[curid] - omega * (f[curid] - f[lasxd]) *0.5 / d[lasxd];
		}
		if (d[lasyd] < omega)
		{
			fc = (f[curid] + f[lasyd]) * 0.5;
		}
		else
		{
			fc = f[curid] - omega * (f[curid] - f[lasyd]) * 0.5 / d[lasyd];
		}
		r[curid] = 0.25 * (2.0f * fa + fb + fc);
	}
}


template<typename T>
__global__ void _invDiscreteGradientTransform_ker(T* vol, T* d, T* res, const T omega, cuint L, cuint W, cuint H)
{
	const unsigned int idX = threadIdx.x + blockIdx.x * blockDim.x;
	const unsigned int idY = threadIdx.y + blockIdx.y * blockDim.y;
	const unsigned int idZ = threadIdx.z + blockIdx.z * blockDim.z;
	if (idX < L && idY < W && idZ < H)
	{
		cuint curId = (idZ * W + idY) * L + idX;
		cuint nexxd = (idZ * W + idY) * L + ((idX + 1) % L);
		cuint nexyd = (idZ * W + ((idY + 1) % W)) * L + idX;
		cuint nexzd = (((idZ + 1) % H) * W + idY) * L + idX;
		cuint lasxd = (idZ * W + idY) * L + ((idX + L - 1) % L);
		cuint lasyd = (idZ * W + ((idY + W - 1) % W)) * L + idX;
		cuint laszd = (((idZ + H - 1) % H) * W + idY) * L + idX;
		T fa, fb, fc, fd;
		if (d[curId] < omega)
		{
			fa = (3.0f * vol[curId] + vol[nexxd] + vol[nexyd] + vol[nexzd]) * 0.1666666666666667f;
		}
		else
		{
			fa = vol[curId] - omega *(3.0f * vol[curId] - vol[nexxd] - vol[nexyd] - vol[nexzd]) * 0.1666666666666667f / d[curId];
		}

		if (d[lasxd] < omega)
		{
			fb = (vol[curId] + vol[lasxd]) * 0.5f;
		}
		else
		{
			fb = vol[curId] - omega * (vol[curId] - vol[lasxd]) * 0.5f / d[lasxd];
		}

		if (d[lasyd] < omega)
		{
			fc = (vol[curId] + vol[lasyd]) * 0.5f;
		}
		else
		{
			fc = vol[curId] - omega * (vol[curId] - vol[lasyd]) * 0.5f / d[lasyd];
		}

		if (d[laszd] < omega)
		{
			fd = (vol[curId] + vol[laszd]) * 0.5f;
		}
		else
		{
			fd = vol[curId] - omega * (vol[curId] - vol[laszd]) * 0.5f / d[laszd];
		}
		res[curId] = 0.1666666666666667f * (3.0f * fa + fb + fc + fd);
	}
}



void invDiscreteGradientTransform(float* f, float* d, float* r, const float omega, cuint L, cuint W, const dim3& blk, const dim3& gid)
{
	_invDiscreteGradientTransform_ker<float> << <gid, blk >> >(f, d, r, omega, L, W);
}

void invDiscreteGradientTransform(double* f, double* d, double* r, const double omega, cuint L, cuint W, const dim3& blk, const dim3& gid)
{
	_invDiscreteGradientTransform_ker<double> << <gid, blk >> >(f, d, r, omega, L, W);
}



void invDiscreteGradientTransform(float* f, float* d, float* r, const float omega, cuint L, cuint W, cuint H, const dim3& blk, const dim3& gid)
{
	_invDiscreteGradientTransform_ker<float> << <gid, blk >> >(f, d, r, omega, L, W, H);
}

void invDiscreteGradientTransform(double* f, double* d, double* r, const double omega, cuint L, cuint W, cuint H, const dim3& blk, const dim3& gid)
{
	_invDiscreteGradientTransform_ker<double> << <gid, blk >> >(f, d, r, omega, L, W, H);
}


//Find the optimum mu
template<typename T>
T OptimizedW0_ker(thrust::device_ptr<T>& TVImg, const T& ObjTV, const unsigned int length)
{
	thrust::device_vector<T> Rec(TVImg, TVImg + length);
	//thrust::device_ptr<T> pRec(&Rec[0]);
	//thrust::copy(TVImg,TVImg + length, pRec);
	thrust::device_vector<T> TCor = Rec;
	thrust::sort(Rec.begin(), Rec.end());
	T TNorm = thrust::reduce(Rec.begin(), Rec.end());

	//thrust::device_ptr<T> pTCor(&TCor[0]);

	if (TNorm <= ObjTV)
	{
		return -1;
	}
	else
	{
		size_t BMin = 0;
		size_t BMax = Rec.size() - 1;
		size_t WInd(static_cast<size_t>((BMin + BMax) * 0.5));
		while (BMax - BMin > 1)
		{
			T tpw = Rec[WInd];
			thrust::transform(Rec.begin(), Rec.end(), TCor.begin(), _softThreshold_functor<T>(tpw));
			T CNorm = thrust::reduce(TCor.begin(), TCor.end());
			if (CNorm > ObjTV)
			{
				BMin = WInd;
			}
			else
			{
				BMax = WInd;
			}
			WInd = static_cast<size_t>(ceil((BMin + BMax) * 0.5));
		}
		return Rec[WInd];
	}
}


template<typename T>
T OptimizedW0_ker(thrust::device_vector<T>& TVImg, const T& ObjTV)
{
	thrust::device_vector<T> Rec = TVImg;
	//thrust::device_ptr<T> pRec(&Rec[0]);
	//thrust::copy(TVImg,TVImg + length, pRec);
	thrust::device_vector<T> TCor = Rec;
	thrust::sort(Rec.begin(), Rec.end());
	T TNorm = thrust::reduce(Rec.begin(), Rec.end());

	//thrust::device_ptr<T> pTCor(&TCor[0]);

	if (TNorm <= ObjTV)
	{
		return -1;
	}
	else
	{
		size_t BMin = 0;
		size_t BMax = Rec.size() - 1;
		size_t WInd(static_cast<size_t>((BMin + BMax) * 0.5f));
		while (BMax - BMin > 1)
		{
			T tpw = Rec[WInd];
			thrust::transform(Rec.begin(), Rec.end(), TCor.begin(), _softThreshold_functor<T>(tpw));
			T CNorm = thrust::reduce(TCor.begin(), TCor.end());
			if (CNorm > ObjTV)
			{
				BMin = WInd;
			}
			else
			{
				BMax = WInd;
			}
			WInd = static_cast<size_t>(ceil((BMin + BMax) * 0.5));
		}
		return Rec[WInd];
	}
}


float OptimizedW0(thrust::device_ptr<float>& TVImg, const float& ObjTV, const unsigned int length)
{
	return OptimizedW0_ker<float>(TVImg, ObjTV, length);
}

double OptimizedW0(thrust::device_ptr<double>& TVImg, const double& ObjTV, const unsigned int length)
{
	return OptimizedW0_ker<double>(TVImg, ObjTV, length);
}

float OptimizedW0(thrust::device_vector<float>& TVImg, const float& ObjTV)
{
	return OptimizedW0_ker<float>(TVImg, ObjTV);
}

double OptimizedW0(thrust::device_vector<double>& TVImg, const double& ObjTV)
{
	return OptimizedW0_ker<double>(TVImg, ObjTV);
}





void SART(thrust::host_vector<float>& himg,
	thrust::host_vector<float>& hprj,
	thrust::host_vector<float>& himgWeg,
	thrust::host_vector<float>& hmsk,
	const FanEAGeo& FanGeo, const Image& Img,
	cuint iterNum, const float lambda)
{
	cuint imgReso = Img.m_Reso.x * Img.m_Reso.y;
	cuint prjReso = FanGeo.m_DetN * FanGeo.m_ViwN;

	thrust::device_vector<float> dimg(himg);
	thrust::device_vector<float> dprj(hprj);
	thrust::device_vector<float> dimgWeg(himgWeg);
	thrust::device_vector<float> dmsk(hmsk);
	thrust::device_vector<float> dcorImg(imgReso, 0);
	thrust::device_vector<float> dcorPrj(prjReso, 0);
	float* pimg = thrust::raw_pointer_cast(&dimg[0]);
	float* pcorImg = thrust::raw_pointer_cast(&dcorImg[0]);
	float* pcorPrj = thrust::raw_pointer_cast(&dcorPrj[0]);
	float* pprj = thrust::raw_pointer_cast(&dprj[0]);

	unsigned int iters = 0;
	//unsigned int prjIdx = 0;
	dim3 prjBlk(128, 8);
	dim3 prjGid(
		(FanGeo.m_DetN + prjBlk.x - 1) / prjBlk.x,
		(FanGeo.m_ViwN + prjBlk.y - 1) / prjBlk.y);
	dim3 bakBlk(32, 32);
	dim3 bakGid(
		(Img.m_Reso.x + bakBlk.x - 1) / bakBlk.x,
		(Img.m_Reso.y + bakBlk.y - 1) / bakBlk.y);

	for (iters = 0; iters != iterNum; ++iters)
	{
		checkCudaErrors(cudaMemset(pcorImg, 0, sizeof(float) * imgReso));
		//Projection
		proj(pimg, pcorPrj, pprj, FanGeo, Img, prjBlk, prjGid);
		//Backprojection
		bakproj_PIXEL(pcorPrj, pcorImg, FanGeo, Img, bakBlk, bakGid);
		//Update Image
		thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(dimg.begin(), dcorImg.begin(), dimgWeg.begin(), dmsk.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(dimg.end(), dcorImg.end(), dimgWeg.end(), dmsk.end())),
			dimg.begin(), _Update<float>(lambda, true, 0, 1500));
	}
	himg = dimg;
}




void SART(thrust::host_vector<float>& himg,
	thrust::host_vector<float>& hprj,
	thrust::host_vector<float>& himgWeg,
	thrust::host_vector<float>& hmsk,
	const FanEDGeo& FanGeo, const Image& Img,
	cuint iterNum, const float lambda)
{
	cuint imgReso = Img.m_Reso.x * Img.m_Reso.y;
	cuint prjReso = FanGeo.m_DetN * FanGeo.m_ViwN;

	thrust::device_vector<float> dimg(himg);
	thrust::device_vector<float> dprj(hprj);
	thrust::device_vector<float> dimgWeg(himgWeg);
	thrust::device_vector<float> dmsk(hmsk);
	thrust::device_vector<float> dcorImg(imgReso, 0);
	thrust::device_vector<float> dcorPrj(prjReso, 0);
	float* pimg = thrust::raw_pointer_cast(&dimg[0]);
	float* pcorImg = thrust::raw_pointer_cast(&dcorImg[0]);
	float* pcorPrj = thrust::raw_pointer_cast(&dcorPrj[0]);
	float* pprj = thrust::raw_pointer_cast(&dprj[0]);

	unsigned int iters = 0;
	//unsigned int prjIdx = 0;
	dim3 prjBlk(128, 8);
	dim3 prjGid(
		(FanGeo.m_DetN + prjBlk.x - 1) / prjBlk.x,
		(FanGeo.m_ViwN + prjBlk.y - 1) / prjBlk.y);
	dim3 bakBlk(32, 32);
	dim3 bakGid(
		(Img.m_Reso.x + bakBlk.x - 1) / bakBlk.x,
		(Img.m_Reso.y + bakBlk.y - 1) / bakBlk.y);

	for (iters = 0; iters != iterNum; ++iters)
	{
		checkCudaErrors(cudaMemset(pcorImg, 0, sizeof(float) * imgReso));
		//Projection
		proj(pimg, pcorPrj, pprj, FanGeo, Img, prjBlk, prjGid); // ÓÐÎÊÌâ;
		//Backprojection
		bakproj_PIXEL(pcorPrj, pcorImg, FanGeo, Img, bakBlk, bakGid);
		//Update Image
		thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(dimg.begin(), dcorImg.begin(), dimgWeg.begin(), dmsk.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(dimg.end(), dcorImg.end(), dimgWeg.end(), dmsk.end())),
			dimg.begin(), _Update<float>(lambda, true, 0, 1500));
	}
	himg = dimg;
}


void OS_SART(thrust::host_vector<float>& himg,
	thrust::host_vector<float>& hprj,
	thrust::host_vector<float>& hweg,
	thrust::host_vector<float>& hmsk,
	const FanEAGeo& FanGeo,
	const Image& Img,
	cuint& iterNum,
	cuint subSetNum,
	const float updateCoef)
{
	cuint imgReso = Img.m_Reso.x * Img.m_Reso.y;
	//	cuint prjReso = FanGeo.m_ViwN * FanGeo.m_DetN;
	uint numPerSubSet = FanGeo.m_ViwN / subSetNum;
	uint realIters = subSetNum * iterNum;
	uint iters = 0;
	uint curSubIdx = 0;

	thrust::device_vector<float> dimg = himg;
	thrust::device_vector<float> dprj = hprj;
	thrust::device_vector<float> dweg = hweg;
	thrust::device_vector<float> dmsk = hmsk;
	thrust::device_vector<float> dcorImg(himg.size(), 0);
	thrust::device_vector<float> dcorPrj(FanGeo.m_DetN * numPerSubSet, 0);

	float* pimg = thrust::raw_pointer_cast(&dimg[0]);
	float* pprj = thrust::raw_pointer_cast(&dprj[0]);
	float* pweg = thrust::raw_pointer_cast(&dweg[0]);
	float* pmsk = thrust::raw_pointer_cast(&dmsk[0]);
	float* pcorImg = thrust::raw_pointer_cast(&dcorImg[0]);
	float* pcorPrj = thrust::raw_pointer_cast(&dcorPrj[0]);

	dim3 prjBlk(256, 4);
	dim3 prjGid((FanGeo.m_DetN + prjBlk.x - 1) / prjBlk.x, (numPerSubSet + prjBlk.y - 1) / prjBlk.y);

	dim3 bakBlk(32, 32);
	dim3 bakGid((Img.m_Reso.x + bakBlk.x - 1) / bakBlk.x, (Img.m_Reso.y + bakBlk.y - 1) / bakBlk.y);

	std::cout << "Iterations begin ";
	for (iters = 0; iters != realIters; ++iters)
	{
		std::cout << ".";

		checkCudaErrors(cudaMemset(pcorImg, 0, sizeof(float) * imgReso));
		//ÔÚsubsetÖÐÍ¶Ó°;
		proj(pimg, pcorPrj, pprj, FanGeo, Img, numPerSubSet, subSetNum, curSubIdx, prjBlk, prjGid);

		//subsetÖÐ·ŽÍ¶Ó°;
		bakproj_PIXEL(pcorPrj, pcorImg, FanGeo, Img, numPerSubSet, subSetNum, curSubIdx, bakBlk, bakGid);
		//update with zip iterator
		thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(dimg.begin(), dcorImg.begin(), dweg.begin(), dmsk.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(dimg.end(), dcorImg.end(), dweg.end(), dmsk.end())),
			dimg.begin(), _Update<float>(updateCoef, false, 0, 0));

		curSubIdx = (curSubIdx + 1) % subSetNum;
	}

	himg = dimg;

	dimg.clear();
	dprj.clear();
	dweg.clear();
	dmsk.clear();
	dcorImg.clear();
	dcorPrj.clear();
}


void OS_SART(float* himg, float* hprj, float* himgWeg, float* hmask, const FanEAGeo& FanGeo, const Image& Img,
	cuint& iterNum, const float& lambda, const unsigned int subSetNum, const float updateCoefficient)
{
	cuint imgSize = Img.m_Reso.x * Img.m_Reso.y;
	cuint prjSize = FanGeo.m_DetN * FanGeo.m_ViwN;

	unsigned int numPerSubSet = FanGeo.m_ViwN / subSetNum; //Ã¿žö×ÓŒ¯ÖÐÓÐ¶àÉÙžöœÇ¶È;
	unsigned int realIters = subSetNum * iterNum;

	unsigned int curSubIdx = 0;

	float* dimg = nullptr;
	float* dprj = nullptr;
	float* dimgWeg = nullptr;
	float* dcorImg = nullptr;
	float* dcorPrj = nullptr;
	float* dmask = nullptr;

	checkCudaErrors(cudaMalloc((void**) &dimg, sizeof(float) * imgSize));
	checkCudaErrors(cudaMalloc((void**) &dprj, sizeof(float) * prjSize));
	checkCudaErrors(cudaMalloc((void**) &dimgWeg, sizeof(float) * imgSize));
	checkCudaErrors(cudaMalloc((void**) &dcorImg, sizeof(float) * imgSize));
	checkCudaErrors(cudaMalloc((void**) &dcorPrj, sizeof(float) * FanGeo.m_DetN * numPerSubSet));
	checkCudaErrors(cudaMalloc((void**) &dmask, sizeof(float) * imgSize));

	checkCudaErrors(cudaMemset(dimg, 0, sizeof(float) * imgSize));
	checkCudaErrors(cudaMemcpy(dprj, hprj, sizeof(float) * prjSize, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dimgWeg, himgWeg, sizeof(float) * imgSize, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dmask, hmask, sizeof(float) * imgSize, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemset(dcorImg, 0, sizeof(float) * imgSize));
	checkCudaErrors(cudaMemset(dcorPrj, 0, sizeof(float) * FanGeo.m_DetN * numPerSubSet));

	unsigned int iters = 0;


	dim3 prjBlk(128, 8);
	dim3 prjGid((FanGeo.m_DetN + prjBlk.x - 1) / prjBlk.x, (numPerSubSet + prjBlk.y - 1) / prjBlk.y);
	dim3 bakBlk(32, 32);
	dim3 bakGid((Img.m_Reso.x + bakBlk.x - 1) / bakBlk.x, (Img.m_Reso.y + bakBlk.y - 1) / bakBlk.y);
	thrust::device_ptr<float> pcorImg(dcorImg);
	thrust::device_ptr<float> pimgWeg(dimgWeg);
	thrust::device_ptr<float> pimg(dimg);
	thrust::device_ptr<float> pmask(dmask);

	std::cout << "Iterations begin ";
	for (iters = 0; iters != realIters; ++iters)
	{
		std::cout << ".";
		checkCudaErrors(cudaMemset(dcorImg, 0, sizeof(float) * imgSize));
		//ÔÚsubsetÖÐÍ¶Ó°;
		proj(dimg, dcorPrj, dprj, true, FanGeo, Img, numPerSubSet, subSetNum, curSubIdx, prjBlk, prjGid);
		//proj_Ker<<<prjGid,prjBlk>>>(dimg, dcorPrj,dprj, true, FanGeo,  Img,  numPerSubSet, subSetNum, curSubIdx);
		//subsetÖÐ·ŽÍ¶Ó°;
		//bakproj_PIXEL_Ker<<<bakGid,bakBlk>>>(dcorPrj, dcorImg, true, true, FanGeo, Img, numPerSubSet, subSetNum, curSubIdx);
		bakproj_PIXEL(dcorPrj, dcorImg, true, true, FanGeo, Img, numPerSubSet, subSetNum, curSubIdx, bakBlk, bakGid);
		//update with zip iterator
		thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(pimg, pcorImg, pimgWeg, pmask)),
			thrust::make_zip_iterator(thrust::make_tuple(pimg + imgSize, pcorImg + imgSize, pimgWeg + imgSize, pmask + imgSize)),
			pimg, _Update<float>(updateCoefficient, false, 0, 0));

		curSubIdx = (curSubIdx + 1) % subSetNum;
	}

	checkCudaErrors(cudaMemcpy(himg, dimg, sizeof(float) * imgSize, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaFree(dprj));
	checkCudaErrors(cudaFree(dimgWeg));
	checkCudaErrors(cudaFree(dcorPrj));
	checkCudaErrors(cudaFree(dimg));
}




void OS_SART_SD(thrust::host_vector<float>& himg,
	thrust::host_vector<float>& hprj,
	thrust::host_vector<float>& hweg,
	thrust::host_vector<float>& hmsk,
	const FanEAGeo& FanGeo,
	const Image& Img,
	cuint& iterNum,
	cuint subSetNum,
	const float updateCoef,
	const float initAlpha,
	const float apa_s,
	cuint descentTime)
{
	cuint imgReso = Img.m_Reso.x * Img.m_Reso.y;
	//	cuint prjReso = FanGeo.m_ViwN * FanGeo.m_DetN;
	uint numPerSubSet = FanGeo.m_ViwN / subSetNum;
	uint realIters = subSetNum * iterNum;
	uint iters = 0;
	uint curSubIdx = 0;
	uint desIdx = 0;

	float apa = 0.0f;
	thrust::device_vector<float> dimg = himg;
	thrust::device_vector<float> dprj = hprj;
	thrust::device_vector<float> dweg = hweg;
	thrust::device_vector<float> dmsk = hmsk;
	thrust::device_vector<float> dcorImg(himg.size(), 0);
	thrust::device_vector<float> dcorPrj(FanGeo.m_DetN * numPerSubSet, 0);
	thrust::device_vector<float> ddiff(himg.size(), 0);

	float* pimg = thrust::raw_pointer_cast(&dimg[0]);
	float* pprj = thrust::raw_pointer_cast(&dprj[0]);
	float* pweg = thrust::raw_pointer_cast(&dweg[0]);
	float* pmsk = thrust::raw_pointer_cast(&dmsk[0]);
	float* pcorImg = thrust::raw_pointer_cast(&dcorImg[0]);
	float* pcorPrj = thrust::raw_pointer_cast(&dcorPrj[0]);
	float* pdiff = thrust::raw_pointer_cast(&ddiff[0]);

	dim3 prjBlk(256, 4);
	dim3 prjGid((FanGeo.m_DetN + prjBlk.x - 1) / prjBlk.x, (numPerSubSet + prjBlk.y - 1) / prjBlk.y);

	dim3 bakBlk(32, 32);
	dim3 bakGid((Img.m_Reso.x + bakBlk.x - 1) / bakBlk.x, (Img.m_Reso.y + bakBlk.y - 1) / bakBlk.y);

	std::cout << "Iterations begin ";
	for (iters = 0; iters != realIters; ++iters)
	{
		std::cout << ".";

		checkCudaErrors(cudaMemset(pcorImg, 0, sizeof(float) * imgReso));
		//ÔÚsubsetÖÐÍ¶Ó°;
		proj(pimg, pcorPrj, pprj, FanGeo, Img, numPerSubSet, subSetNum, curSubIdx, prjBlk, prjGid);

		//subsetÖÐ·ŽÍ¶Ó°;
		bakproj_PIXEL(pcorPrj, pcorImg, FanGeo, Img, numPerSubSet, subSetNum, curSubIdx, bakBlk, bakGid);
		//update with zip iterator
		thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(dimg.begin(), dcorImg.begin(), dweg.begin(), dmsk.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(dimg.end(), dcorImg.end(), dweg.end(), dmsk.end())),
			dimg.begin(), _Update<float>(updateCoef, false, 0, 0));
		if (iters > subSetNum * 2)
		{
			apa = initAlpha;
			//TV SD minimization
			for (desIdx = 0; desIdx != descentTime; ++desIdx)
			{
				//_GradientOfTV_ker<float><<<bakGid,bakBlk>>>(dimg, ddiff, 1.0, Img.m_Reso.x,Img.m_Reso.y);
				GradientOfTV(pimg, pdiff, 1.0, Img.m_Reso.x, Img.m_Reso.y, bakBlk, bakGid);
				thrust::transform(dimg.begin(), dimg.end(), ddiff.begin(), dimg.begin(), _saxpy_functor<float>(-apa));
				apa *= apa_s;
			}
		}
		curSubIdx = (curSubIdx + 1) % subSetNum;
	}

	himg = dimg;

	dimg.clear();
	dprj.clear();
	dweg.clear();
	dmsk.clear();
	dcorImg.clear();
	dcorPrj.clear();
}


void OS_SART_SD(float* himg, float* hprj, float* himgWeg, float* hmask, const FanEAGeo& FanGeo, const Image& Img,
	cuint& iterNum, const float& lambda, const unsigned int subSetNum, float initAlpha, float apa_s, cuint descentTime, const float updateCoefficient)
{
	cuint imgSize = Img.m_Reso.x * Img.m_Reso.y;
	cuint prjSize = FanGeo.m_DetN * FanGeo.m_ViwN;

	unsigned int numPerSubSet = FanGeo.m_ViwN / subSetNum; //Ã¿žö×ÓŒ¯ÖÐÓÐ¶àÉÙžöœÇ¶È;
	unsigned int realIters = subSetNum * iterNum;

	unsigned int curSubIdx = 0;

	float* dimg = nullptr;
	float* dprj = nullptr;
	float* dimgWeg = nullptr;
	float* dcorImg = nullptr;
	float* dcorPrj = nullptr;
	float* dmask = nullptr;

	float* ddiff = nullptr;

	checkCudaErrors(cudaMalloc((void**) &dimg, sizeof(float) * imgSize));
	checkCudaErrors(cudaMalloc((void**) &dprj, sizeof(float) * prjSize));
	checkCudaErrors(cudaMalloc((void**) &dimgWeg, sizeof(float) * imgSize));
	checkCudaErrors(cudaMalloc((void**) &dcorImg, sizeof(float) * imgSize));
	checkCudaErrors(cudaMalloc((void**) &dcorPrj, sizeof(float) * FanGeo.m_DetN * numPerSubSet));
	checkCudaErrors(cudaMalloc((void**) &dmask, sizeof(float) * imgSize));
	checkCudaErrors(cudaMalloc((void**) &ddiff, sizeof(float) * imgSize));


	checkCudaErrors(cudaMemset(dimg, 0, sizeof(float) * imgSize));
	checkCudaErrors(cudaMemcpy(dprj, hprj, sizeof(float) * prjSize, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dimgWeg, himgWeg, sizeof(float) * imgSize, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dmask, hmask, sizeof(float) * imgSize, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemset(dcorImg, 0, sizeof(float) * imgSize));
	checkCudaErrors(cudaMemset(dcorPrj, 0, sizeof(float) * FanGeo.m_DetN * numPerSubSet));
	checkCudaErrors(cudaMemset(ddiff, 0, sizeof(float) * imgSize));

	unsigned int iters = 0;


	dim3 prjBlk(128, 4);
	dim3 prjGid((FanGeo.m_DetN + prjBlk.x - 1) / prjBlk.x, (numPerSubSet + prjBlk.y - 1) / prjBlk.y);
	dim3 bakBlk(32, 32);
	dim3 bakGid((Img.m_Reso.x + bakBlk.x - 1) / bakBlk.x, (Img.m_Reso.y + bakBlk.y - 1) / bakBlk.y);
	thrust::device_ptr<float> pcorImg(dcorImg);
	thrust::device_ptr<float> pimgWeg(dimgWeg);
	thrust::device_ptr<float> pimg(dimg);
	thrust::device_ptr<float> pmask(dmask);
	thrust::device_ptr<float> pdiff(ddiff);
	//unsigned int descentTime = 6;
	unsigned int desIdx = 0;
	//float initAlpha = 0.0002f;
	float apa = initAlpha;
	//float apa_s = 0.995f;

	std::cout << "Iterations begin ";
	for (iters = 0; iters != realIters; ++iters)
	{
		std::cout << ".";
		checkCudaErrors(cudaMemset(dcorImg, 0, sizeof(float) * imgSize));
		//ÔÚsubsetÖÐÍ¶Ó°;
		//proj_Ker<<<prjGid,prjBlk>>>(dimg, dcorPrj,dprj, true, FanGeo,  Img,  numPerSubSet, subSetNum, curSubIdx);
		proj(dimg, dcorPrj, dprj, true, FanGeo, Img, numPerSubSet, subSetNum, curSubIdx, prjBlk, prjGid);
		//subsetÖÐ·ŽÍ¶Ó°;
		//bakproj_PIXEL_Ker<<<bakGid,bakBlk>>>(dcorPrj, dcorImg, true, true, FanGeo, Img, numPerSubSet, subSetNum, curSubIdx);
		bakproj_PIXEL(dcorPrj, dcorImg, true, true, FanGeo, Img, numPerSubSet, subSetNum, curSubIdx, bakBlk, bakGid);
		//update with zip iterator
		thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(pimg, pcorImg, pimgWeg, pmask)),
			thrust::make_zip_iterator(thrust::make_tuple(pimg + imgSize, pcorImg + imgSize, pimgWeg + imgSize, pmask + imgSize)),
			pimg, _Update<float>(updateCoefficient, false, 0, 0));

		apa = initAlpha;
		//TV SD minimization
		for (desIdx = 0; desIdx != descentTime; ++desIdx)
		{
			//_GradientOfTV_ker<float><<<bakGid,bakBlk>>>(dimg, ddiff, 1.0, Img.m_Reso.x,Img.m_Reso.y);
			GradientOfTV(dimg, ddiff, 1.0, Img.m_Reso.x, Img.m_Reso.y, bakBlk, bakGid);
			thrust::transform(pimg, pimg + imgSize, pdiff, pimg, _saxpy_functor<float>(-apa));
			apa *= apa_s;
		}

		curSubIdx = (curSubIdx + 1) % subSetNum;
	}

	checkCudaErrors(cudaMemcpy(himg, dimg, sizeof(float) * imgSize, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaFree(dprj));
	checkCudaErrors(cudaFree(dimgWeg));
	checkCudaErrors(cudaFree(dcorPrj));
	checkCudaErrors(cudaFree(dimg));
}



//1.079698460635159e+07;
void OS_SART_STF(float* himg, float* hprj, float* himgWeg, float* hmask, const FanEAGeo& FanGeo, const Image& Img,
	cuint& iterNum, const float& lambda, const unsigned int subSetNum, const float objTV, const float updateCoefficient)
{
	cuint imgSize = Img.m_Reso.x * Img.m_Reso.y;
	cuint prjSize = FanGeo.m_DetN * FanGeo.m_ViwN;

	unsigned int numPerSubSet = FanGeo.m_ViwN / subSetNum; //Ã¿žö×ÓŒ¯ÖÐÓÐ¶àÉÙžöœÇ¶È;
	unsigned int realIters = subSetNum * iterNum;

	unsigned int curSubIdx = 0;

	float* dimg = nullptr;
	float* dprj = nullptr;
	float* dimgWeg = nullptr;
	float* dcorImg = nullptr;
	float* dcorPrj = nullptr;
	float* dmask = nullptr;
	float* dres = nullptr;
	float* ddiff = nullptr;

	checkCudaErrors(cudaMalloc((void**) &dimg, sizeof(float) * imgSize));
	checkCudaErrors(cudaMalloc((void**) &dprj, sizeof(float) * prjSize));
	checkCudaErrors(cudaMalloc((void**) &dimgWeg, sizeof(float) * imgSize));
	checkCudaErrors(cudaMalloc((void**) &dcorImg, sizeof(float) * imgSize));
	checkCudaErrors(cudaMalloc((void**) &dcorPrj, sizeof(float) * FanGeo.m_DetN * numPerSubSet));
	checkCudaErrors(cudaMalloc((void**) &dmask, sizeof(float) * imgSize));
	checkCudaErrors(cudaMalloc((void**) &ddiff, sizeof(float) * imgSize));
	checkCudaErrors(cudaMalloc((void**) &dres, sizeof(float) * imgSize));

	checkCudaErrors(cudaMemset(dimg, 0, sizeof(float) * imgSize));
	checkCudaErrors(cudaMemcpy(dprj, hprj, sizeof(float) * prjSize, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dimgWeg, himgWeg, sizeof(float) * imgSize, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dmask, hmask, sizeof(float) * imgSize, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemset(dcorImg, 0, sizeof(float) * imgSize));
	checkCudaErrors(cudaMemset(dcorPrj, 0, sizeof(float) * FanGeo.m_DetN * numPerSubSet));
	checkCudaErrors(cudaMemset(ddiff, 0, sizeof(float) * imgSize));
	checkCudaErrors(cudaMemset(dres, 0, sizeof(float) * imgSize));
	unsigned int iters = 0;


	dim3 prjBlk(128, 4);
	dim3 prjGid((FanGeo.m_DetN + prjBlk.x - 1) / prjBlk.x, (numPerSubSet + prjBlk.y - 1) / prjBlk.y);
	dim3 bakBlk(32, 32);
	dim3 bakGid((Img.m_Reso.x + bakBlk.x - 1) / bakBlk.x, (Img.m_Reso.y + bakBlk.y - 1) / bakBlk.y);
	thrust::device_ptr<float> pcorImg(dcorImg);
	thrust::device_ptr<float> pimgWeg(dimgWeg);
	thrust::device_ptr<float> pimg(dimg);
	thrust::device_ptr<float> pmask(dmask);
	thrust::device_ptr<float> pdiff(ddiff);
	thrust::device_ptr<float> pres(dres);

	float mu = 0.0f;
	std::cout << "Iterations begin ";
	for (iters = 0; iters != realIters; ++iters)
	{
		std::cout << ".";
		checkCudaErrors(cudaMemset(dcorImg, 0, sizeof(float) * imgSize));
		//ÔÚsubsetÖÐÍ¶Ó°;
		proj(dimg, dcorPrj, dprj, true, FanGeo, Img, numPerSubSet, subSetNum, curSubIdx, prjBlk, prjGid);
		//subsetÖÐ·ŽÍ¶Ó°;
		bakproj_PIXEL(dcorPrj, dcorImg, true, true, FanGeo, Img, numPerSubSet, subSetNum, curSubIdx, bakBlk, bakGid);

		//update with zip iterator
		thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(pimg, pcorImg, pimgWeg, pmask)),
			thrust::make_zip_iterator(thrust::make_tuple(pimg + imgSize, pcorImg + imgSize, pimgWeg + imgSize, pmask + imgSize)),
			pimg, _Update<float>(updateCoefficient, false, 0, 0));


		//_DiscreteGradientTrans_ker<float><<<bakGid,bakBlk>>>(dimg,ddiff,1.0f, Img.m_Reso.x, Img.m_Reso.y);
		DiscreteGradientTrans(dimg, ddiff, 1.0f, Img.m_Reso.x, Img.m_Reso.y, bakBlk, bakGid);
		mu = OptimizedW0(pdiff, objTV, Img.m_Reso.x * Img.m_Reso.y);
		if (mu == -1)
		{
			std::cout << "No filtering\n";
			curSubIdx = (curSubIdx + 1) % subSetNum;
			continue;
		}
		else
		{
			std::cout << mu << std::endl;
		}
		//_invDiscreteGradientTransform_ker<float><<<bakGid,bakBlk>>>(dimg, ddiff, dres,mu,Img.m_Reso.x,Img.m_Reso.y);
		invDiscreteGradientTransform(dimg, ddiff, dres, mu, Img.m_Reso.x, Img.m_Reso.y, bakBlk, bakGid);
		thrust::swap(pimg, pres);


		curSubIdx = (curSubIdx + 1) % subSetNum;
	}

	checkCudaErrors(cudaMemcpy(himg, dimg, sizeof(float) * imgSize, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaFree(dprj));
	checkCudaErrors(cudaFree(dimgWeg));
	checkCudaErrors(cudaFree(dcorPrj));
	checkCudaErrors(cudaFree(dimg));
}


//1.079698460635159e+07;
void OS_SART_STF(thrust::host_vector<float>& himg,
	thrust::host_vector<float>& hprj,
	thrust::host_vector<float>& hweg,
	thrust::host_vector<float>& hmsk,
	const FanEAGeo& FanGeo, const Image& Img,
	cuint& iterNum, const unsigned int subSetNum,
	const float objTV, const float updateCoef)
{

	cuint imgReso = Img.m_Reso.x * Img.m_Reso.y;
	//	cuint prjReso = FanGeo.m_ViwN * FanGeo.m_DetN;
	uint numPerSubSet = FanGeo.m_ViwN / subSetNum;
	uint realIters = subSetNum * iterNum;
	uint iters = 0;
	uint curSubIdx = 0;



	float mu = -1.0f;
	thrust::device_vector<float> dimg = himg;
	thrust::device_vector<float> dprj = hprj;
	thrust::device_vector<float> dweg = hweg;
	thrust::device_vector<float> dmsk = hmsk;
	thrust::device_vector<float> dcorImg(himg.size(), 0);
	thrust::device_vector<float> dcorPrj(FanGeo.m_DetN * numPerSubSet, 0);
	thrust::device_vector<float> ddiff(himg.size(), 0);
	thrust::device_vector<float> dres(himg.size(), 0);

	float* pimg = thrust::raw_pointer_cast(&dimg[0]);
	float* pprj = thrust::raw_pointer_cast(&dprj[0]);
	float* pweg = thrust::raw_pointer_cast(&dweg[0]);
	float* pmsk = thrust::raw_pointer_cast(&dmsk[0]);
	float* pcorImg = thrust::raw_pointer_cast(&dcorImg[0]);
	float* pcorPrj = thrust::raw_pointer_cast(&dcorPrj[0]);
	float* pdiff = thrust::raw_pointer_cast(&ddiff[0]);
	float* pres = thrust::raw_pointer_cast(&dres[0]);

	dim3 prjBlk(256, 4);
	dim3 prjGid((FanGeo.m_DetN + prjBlk.x - 1) / prjBlk.x, (numPerSubSet + prjBlk.y - 1) / prjBlk.y);

	dim3 bakBlk(32, 32);
	dim3 bakGid((Img.m_Reso.x + bakBlk.x - 1) / bakBlk.x, (Img.m_Reso.y + bakBlk.y - 1) / bakBlk.y);

	std::cout << "Iterations begin ";
	for (iters = 0; iters != realIters; ++iters)
	{
		std::cout << ".";

		checkCudaErrors(cudaMemset(pcorImg, 0, sizeof(float) * imgReso));
		//ÔÚsubsetÖÐÍ¶Ó°;
		proj(pimg, pcorPrj, pprj, FanGeo, Img, numPerSubSet, subSetNum, curSubIdx, prjBlk, prjGid);

		//subsetÖÐ·ŽÍ¶Ó°;
		bakproj_PIXEL(pcorPrj, pcorImg, FanGeo, Img, numPerSubSet, subSetNum, curSubIdx, bakBlk, bakGid);
		//update with zip iterator
		thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(dimg.begin(), dcorImg.begin(), dweg.begin(), dmsk.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(dimg.end(), dcorImg.end(), dweg.end(), dmsk.end())),
			dimg.begin(), _Update<float>(updateCoef, false, 0, 0));

		//_DiscreteGradientTrans_ker<float><<<bakGid,bakBlk>>>(dimg,ddiff,1.0f, Img.m_Reso.x, Img.m_Reso.y);
		DiscreteGradientTrans(pimg, pdiff, 1.0f, Img.m_Reso.x, Img.m_Reso.y, bakBlk, bakGid);
		mu = OptimizedW0(ddiff, objTV);
		if (mu != -1)
		{
			std::cout << "mu is " << mu << std::endl;
			invDiscreteGradientTransform(pimg, pdiff, pres, mu, Img.m_Reso.x, Img.m_Reso.y, bakBlk, bakGid);
			thrust::swap(pimg, pres);
		}
		curSubIdx = (curSubIdx + 1) % subSetNum;
	}

	himg = dimg;

	dimg.clear();
	dprj.clear();
	dweg.clear();
	dmsk.clear();
	dcorImg.clear();
	dcorPrj.clear();
}



__global__ void updateImgMultiSlices(float* dimg, float* dcorImg, float* dwegImg, float* dmask, float lambda, cuint ImgLR, cuint ImgWR, cuint sliceNum)
{
	cuint l = threadIdx.x + blockIdx.x * blockDim.x;
	cuint w = threadIdx.y + blockIdx.y * blockDim.y;
	cuint s = threadIdx.z + blockIdx.z * blockDim.z;
	if (l < ImgLR && w < ImgWR && s < sliceNum)
	{
		cuint i = (s * ImgWR + w) * ImgLR + l;
		dimg[i] += lambda * dcorImg[i] / dwegImg[w * ImgLR + l] * dmask[w * ImgLR + l];
		if (dimg[i] < 0)
		{
			dimg[i] = 0;
		}
	}
}

void SART_MULTISLICES(
	thrust::host_vector<float>& himg,
	thrust::host_vector<float>& hprj,
	thrust::host_vector<float>& hweg,
	thrust::host_vector<float>& hmsk,
	const FanEAGeo& FanGeo,
	const Image& Img,
	cuint iterNum,
	const float lambda,
	cuint sliceNum)
{
	cuint sliceReso = Img.m_Reso.x * Img.m_Reso.y;
	cuint projReso = FanGeo.m_DetN * FanGeo.m_ViwN;
	cuint totVolReso = sliceReso * sliceNum;
	//cuint totPrjReso = projReso * sliceNum;

	thrust::device_vector<float> dimg(himg);
	thrust::device_vector<float> dprj(hprj);
	thrust::device_vector<float> dweg(hweg);
	thrust::device_vector<float> dmsk(hmsk);
	thrust::device_vector<float> dcorImg(sliceReso * sliceNum, 0);
	thrust::device_vector<float> dcorPrj(projReso * sliceNum, 0);

	float* pcorImg = thrust::raw_pointer_cast(&dcorImg[0]);
	float* pcorPrj = thrust::raw_pointer_cast(&dcorPrj[0]);
	float* pimg = thrust::raw_pointer_cast(&dimg[0]);
	float* pprj = thrust::raw_pointer_cast(&dprj[0]);
	float* pweg = thrust::raw_pointer_cast(&dweg[0]);
	float* pmsk = thrust::raw_pointer_cast(&dmsk[0]);

	//Set the projection, backprojection and volume updating block grid sizes
	dim3 prjBlk(128, 8);
	dim3 prjGid(
		(FanGeo.m_DetN + 127) / 128,
		(FanGeo.m_ViwN + 7) / 8);

	dim3 bakBlk(32, 32);
	dim3 bakGid(
		(Img.m_Reso.x + 31) / 32,
		(Img.m_Reso.y + 31) / 32);

	dim3 updBlk(16, 16, 4);
	dim3 updGid(
		(Img.m_Reso.x + 15) / 16,
		(Img.m_Reso.y + 15) / 16,
		(sliceNum + 3) / 4);

	uint iters = 0;
	while (iters != iterNum)
	{
		checkCudaErrors(cudaMemset(pcorImg, 0, sizeof(float) * totVolReso));

		//Proj for all slices;
		proj(pimg, pcorPrj, pprj, FanGeo, Img, prjBlk, prjGid, sliceNum);
		//Back proj for all slices;
		bakproj_PIXEL(pcorPrj, pcorImg, FanGeo, Img, bakBlk, bakGid, sliceNum);

		//Update the volume in all slices;
		updateImgMultiSlices << <updGid, updBlk >> >(pimg, pcorImg, pweg, pmsk, lambda, Img.m_Reso.x, Img.m_Reso.y, sliceNum);
		++iters;
	}

	himg = dimg;
}




void DEMO1() //This is the demo of OS-SART
{
	FanEAGeo FanGeo;
	Image Img;
	cuint imgReso = Img.m_Reso.x * Img.m_Reso.y;
	cuint prjReso = FanGeo.m_DetN * FanGeo.m_ViwN;
	thrust::host_vector<float> himg(imgReso, 0);
	thrust::host_vector<float> hweg(imgReso, 0);
	thrust::host_vector<float> hprj(prjReso, 0);
	thrust::host_vector<float> hmsk(imgReso, 0);


	std::string FileName = "demo1_patientPrj.prj";
	std::string FouName = "demo1_reconImg512x512.raw";
	std::string WegName = "demo1_bakwegEA.raw";
	std::string MaskName = "demo1_newmask2.img";

	//Read the data;


	std::ifstream fid(FileName.c_str(), std::ios::binary);
	std::ifstream fons(WegName.c_str(), std::ios::binary);
	std::ifstream fmsk(MaskName.c_str(), std::ios::binary);

	if (!(fid.is_open() && fons.is_open() && fmsk.is_open()))
	{
		std::cerr << "Cannot open the file file for demo1 \n";
		exit(-1);
	}

	fid.read((char*) (&hprj[0]), sizeof(float) * prjReso);
	fons.read((char*) (&hweg[0]), sizeof(float) * imgReso);
	fmsk.read((char*) (&hmsk[0]), sizeof(float) * imgReso);



	fid.close(); fons.close(); fmsk.close();

	cuint iterNum = 20;
	const float	lambda = 2.0;
	cuint subSetNum = 220;
	const float updateCoef = lambda * 2200 / subSetNum;
	//Begin OS-SART;
	OS_SART(himg, hprj, hweg, hmsk, FanGeo, Img, iterNum, subSetNum, updateCoef);

	std::ofstream fou(FouName.c_str(), std::ios::binary);
	fou.write((char*) (&himg[0]), sizeof(float) * imgReso);
	fou.close();

}





void DEMO1_1(const FanEAGeo& FanGeo, const Image& Img,
	const std::string& FileName,
	const std::string& WegName,
	const std::string& MaskName,
	const std::string& FouName,
	cuint iterNum,
	const float lambda,
	cuint subSetNum) //This is the demo of OS-SART
{
	cuint imgReso = Img.m_Reso.x * Img.m_Reso.y;
	cuint prjReso = FanGeo.m_DetN * FanGeo.m_ViwN;
	thrust::host_vector<float> himg(imgReso, 0);
	thrust::host_vector<float> hweg(imgReso, 0);
	thrust::host_vector<float> hprj(prjReso, 0);
	thrust::host_vector<float> hmsk(imgReso, 0);

	//Read the data;
	std::ifstream fid(FileName.c_str(), std::ios::binary);
	std::ifstream fons(WegName.c_str(), std::ios::binary);
	std::ifstream fmsk(MaskName.c_str(), std::ios::binary);

	if (!(fid.is_open() && fons.is_open() && fmsk.is_open()))
	{
		std::cerr << "Cannot open the file file for demo1 \n";
		exit(-1);
	}

	fid.read((char*) (&hprj[0]), sizeof(float) * prjReso);
	fons.read((char*) (&hweg[0]), sizeof(float) * imgReso);
	fmsk.read((char*) (&hmsk[0]), sizeof(float) * imgReso);

	fid.close(); fons.close(); fmsk.close();

	//cuint iterNum = 20;
	//const float	lambda = 2.0;
	//cuint subSetNum = 220;
	const float updateCoef = lambda * FanGeo.m_ViwN / subSetNum;
	//Begin OS-SART;
	OS_SART(himg, hprj, hweg, hmsk, FanGeo, Img, iterNum, subSetNum, updateCoef);

	std::ofstream fou(FouName.c_str(), std::ios::binary);
	fou.write((char*) (&himg[0]), sizeof(float) * imgReso);
	fou.close();

}


//This is the demo of OS-SART + SD
void DEMO2_1(const FanEAGeo& FanGeo, const Image& Img,
	const std::string& FileName,
	const std::string& WegName,
	const std::string& MaskName,
	const std::string& FouName,
	cuint iterNum,
	const float lambda,
	cuint subSetNum,
	const float initAlpha,
	const float apa_s,
	cuint descentTime)
{
	//FanEAGeo FanGeo;
	//Image Img;
	cuint imgReso = Img.m_Reso.x * Img.m_Reso.y;
	cuint prjReso = FanGeo.m_DetN * FanGeo.m_ViwN;
	thrust::host_vector<float> himg(imgReso, 0);
	thrust::host_vector<float> hweg(imgReso, 0);
	thrust::host_vector<float> hprj(prjReso, 0);
	thrust::host_vector<float> hmsk(imgReso, 0);

	//std::string FileName = "demo1_patientPrj.prj";
	//std::string FouName = "demo2_reconImg512x512.raw";
	//std::string WegName = "demo1_bakwegEA.raw";
	//std::string MaskName = "demo1_newmask2.img";

	//Read the data;

	std::ifstream fid(FileName.c_str(), std::ios::binary);
	std::ifstream fons(WegName.c_str(), std::ios::binary);
	std::ifstream fmsk(MaskName.c_str(), std::ios::binary);

	if (!(fid.is_open() && fons.is_open() && fmsk.is_open()))
	{
		std::cerr << "Cannot open the file file for demo1 \n";
		exit(-1);
	}

	fid.read((char*) (&hprj[0]), sizeof(float) * prjReso);
	fons.read((char*) (&hweg[0]), sizeof(float) * imgReso);
	fmsk.read((char*) (&hmsk[0]), sizeof(float) * imgReso);

	fid.close(); fons.close(); fmsk.close();

	//cuint iterNum = 20;
	//const float	lambda = 2.0;
	//cuint subSetNum = 220;
	const float updateCoef = lambda * FanGeo.m_ViwN / subSetNum;
	//const float initAlpha = 0.0002;
	//const float apa_s = 0.995;
	//cuint descentTime = 6;
	//Begin OS-SART;
	OS_SART_SD(himg, hprj, hweg, hmsk, FanGeo, Img, iterNum, subSetNum, updateCoef, initAlpha, apa_s, descentTime);

	writeData_host(FouName, himg, imgReso);
}





//This is the demo of OS-SART + SD
void DEMO2()
{
	FanEAGeo FanGeo;
	Image Img;
	cuint imgReso = Img.m_Reso.x * Img.m_Reso.y;
	cuint prjReso = FanGeo.m_DetN * FanGeo.m_ViwN;
	thrust::host_vector<float> himg(imgReso, 0);
	thrust::host_vector<float> hweg(imgReso, 0);
	thrust::host_vector<float> hprj(prjReso, 0);
	thrust::host_vector<float> hmsk(imgReso, 0);

	std::string FileName = "demo1_patientPrj.prj";
	std::string FouName = "demo2_reconImg512x512.raw";
	std::string WegName = "demo1_bakwegEA.raw";
	std::string MaskName = "demo1_newmask2.img";

	//Read the data;

	std::ifstream fid(FileName.c_str(), std::ios::binary);
	std::ifstream fons(WegName.c_str(), std::ios::binary);
	std::ifstream fmsk(MaskName.c_str(), std::ios::binary);

	if (!(fid.is_open() && fons.is_open() && fmsk.is_open()))
	{
		std::cerr << "Cannot open the file file for demo1 \n";
		exit(-1);
	}

	fid.read((char*) (&hprj[0]), sizeof(float) * prjReso);
	fons.read((char*) (&hweg[0]), sizeof(float) * imgReso);
	fmsk.read((char*) (&hmsk[0]), sizeof(float) * imgReso);

	fid.close(); fons.close(); fmsk.close();

	cuint iterNum = 20;
	const float	lambda = 2.0;
	cuint subSetNum = 220;
	const float updateCoef = lambda * 2200 / subSetNum;
	const float initAlpha = 0.0002f;
	const float apa_s = 0.995f;
	cuint descentTime = 6;
	//Begin OS-SART;
	OS_SART_SD(himg, hprj, hweg, hmsk, FanGeo, Img, iterNum, subSetNum, updateCoef, initAlpha, apa_s, descentTime);

	writeData_host(FouName, himg, imgReso);
}





void DEMO3()//This is the demo of OS-SART + STF
{
	FanEAGeo FanGeo;
	Image Img;
	cuint imgReso = Img.m_Reso.x * Img.m_Reso.y;
	cuint prjReso = FanGeo.m_DetN * FanGeo.m_ViwN;
	thrust::host_vector<float> himg(imgReso, 0);
	thrust::host_vector<float> hweg(imgReso, 0);
	thrust::host_vector<float> hprj(prjReso, 0);
	thrust::host_vector<float> hmsk(imgReso, 0);

	std::string FileName = "demo1_patientPrj.prj";
	std::string FouName = "demo3_reconImg512x512.raw";
	std::string WegName = "demo1_bakwegEA.raw";
	std::string MaskName = "demo1_newmask2.img";

	//Read the data;

	std::ifstream fid(FileName.c_str(), std::ios::binary);
	std::ifstream fons(WegName.c_str(), std::ios::binary);
	std::ifstream fmsk(MaskName.c_str(), std::ios::binary);

	if (!(fid.is_open() && fons.is_open() && fmsk.is_open()))
	{
		std::cerr << "Cannot open the file file for demo1 \n";
		exit(-1);
	}

	fid.read((char*) (&hprj[0]), sizeof(float) * prjReso);
	fons.read((char*) (&hweg[0]), sizeof(float) * imgReso);
	fmsk.read((char*) (&hmsk[0]), sizeof(float) * imgReso);

	fid.close(); fons.close(); fmsk.close();

	cuint iterNum = 20;
	const float	lambda = 2.0;
	cuint subSetNum = 220;
	const float updateCoef = lambda * 2200 / subSetNum;
	const float objTV = static_cast<float>(1.079698460635159e+07);
	//Begin OS-SART;

	OS_SART_STF(himg, hprj, hweg, hmsk, FanGeo, Img, iterNum, subSetNum, objTV, updateCoef);
	writeData_host(FouName, himg, imgReso);
}






void DEMO3_1(const FanEAGeo& FanGeo, const Image& Img,
	const std::string& FileName,
	const std::string& WegName,
	const std::string& MaskName,
	const std::string& FouName,
	cuint iterNum,
	const float lambda,
	cuint subSetNum,
	const float objTV)//This is the demo of OS-SART + STF
{
	//FanEAGeo FanGeo;
	//Image Img;
	cuint imgReso = Img.m_Reso.x * Img.m_Reso.y;
	cuint prjReso = FanGeo.m_DetN * FanGeo.m_ViwN;
	thrust::host_vector<float> himg(imgReso, 0);
	thrust::host_vector<float> hweg(imgReso, 0);
	thrust::host_vector<float> hprj(prjReso, 0);
	thrust::host_vector<float> hmsk(imgReso, 0);

	/*std::string FileName = "demo1_patientPrj.prj";
	std::string FouName = "demo3_reconImg512x512.raw";
	std::string WegName = "demo1_bakwegEA.raw";
	std::string MaskName = "demo1_newmask2.img";*/

	//Read the data;

	std::ifstream fid(FileName.c_str(), std::ios::binary);
	std::ifstream fons(WegName.c_str(), std::ios::binary);
	std::ifstream fmsk(MaskName.c_str(), std::ios::binary);

	if (!(fid.is_open() && fons.is_open() && fmsk.is_open()))
	{
		std::cerr << "Cannot open the file file for demo1 \n";
		exit(-1);
	}

	fid.read((char*) (&hprj[0]), sizeof(float) * prjReso);
	fons.read((char*) (&hweg[0]), sizeof(float) * imgReso);
	fmsk.read((char*) (&hmsk[0]), sizeof(float) * imgReso);

	fid.close(); fons.close(); fmsk.close();

	//cuint iterNum = 20;
	//const float	lambda = 2.0;
	//cuint subSetNum = 220;
	const float updateCoef = lambda * FanGeo.m_ViwN / subSetNum;
	//const float objTV = 1.079698460635159e+07;
	//Begin OS-SART;

	OS_SART_STF(himg, hprj, hweg, hmsk, FanGeo, Img, iterNum, subSetNum, objTV, updateCoef);
	writeData_host(FouName, himg, imgReso);
}




void DEMO4()
{
	FanEAGeo FanGeo;
	Image Img;

	FanGeo.m_DetArc = 0.95928517242269f;
	FanGeo.m_DetN = 1024;
	FanGeo.m_DetCntIdx = FanGeo.m_DetN * 0.5f;
	FanGeo.m_DetStp = FanGeo.m_DetArc / FanGeo.m_DetN;
	FanGeo.m_O2D = 4.082259521484375e+02;

	FanGeo.m_S2O = 5.385200195312500e+02;
	FanGeo.m_S2D = FanGeo.m_S2O + FanGeo.m_O2D;
	FanGeo.m_ViwBeg = 0.0f;// (-2.668082275390625e+02 / 180.0* 3.14159265358979323846264);
	FanGeo.m_ViwN = 360;
	FanGeo.m_ViwStp = 3.14159265358979323846264f * 2.0f / FanGeo.m_ViwN;

	Img.m_Bias.x = 0.0f;
	Img.m_Bias.y = 0.0f;  //ÕâžöÆ«ÒÆµÄµ¥Î»ÊÇÕæÊµÎïÀíµ¥Î»;
	Img.m_Reso.x = Img.m_Reso.y = 512;
	Img.m_Size.x = Img.m_Size.y = static_cast<float>(4.484740011196460e+02);
	Img.m_Step.x = Img.m_Size.x / Img.m_Reso.x;
	Img.m_Step.y = Img.m_Size.y / Img.m_Reso.y;

	const unsigned int imgReso = Img.m_Reso.x * Img.m_Reso.y;
	const unsigned int prjReso = FanGeo.m_DetN * FanGeo.m_ViwN;
	thrust::host_vector<float> h_oimg(imgReso, 0);
	thrust::host_vector<float> h_prj(prjReso, 0);
	thrust::host_vector<float> h_reim(imgReso, 0);


	//¶ÁÈ¡Ò»žöÍŒÏñÎÄŒþ;
	std::string imageFile = "demo4_img.raw";
	std::ifstream fin(imageFile.c_str(), std::ios::binary);
	fin.read((char*) (&h_oimg[0]), sizeof(float) * imgReso);
	fin.close();

	thrust::device_vector<float> d_oimg = h_oimg;
	thrust::device_vector<float> d_prj = h_prj;
	thrust::device_vector<float> d_reim = h_reim;
	thrust::device_vector<float> d_msk(imgReso, 1.0);

	float* p_oimg = thrust::raw_pointer_cast(&d_oimg[0]);
	float* p_prj = thrust::raw_pointer_cast(&d_prj[0]);
	float* p_reim = thrust::raw_pointer_cast(&d_reim[0]);
	float* p_msk = thrust::raw_pointer_cast(&d_msk[0]);
	dim3 prjBlk(256, 4);
	dim3 prjGid(
		(FanGeo.m_DetN + prjBlk.x - 1) / prjBlk.x,
		(FanGeo.m_ViwN + prjBlk.y - 1) / prjBlk.y);

	ProjectionThenBackProjection(p_oimg, p_reim, FanGeo, Img);
	/*
	dim3 bakBlk(32,32);
	dim3 bakGid(
	(Img.m_Reso.x + bakBlk.x - 1) / bakBlk.x,
	(Img.m_Reso.y + bakBlk.y - 1) / bakBlk.y);


	proj(p_oimg,p_prj,FanGeo,Img,prjBlk,prjGid);

	bakproj_PIXEL(p_prj,p_reim,p_msk,FanGeo,Img,bakBlk,bakGid);*/

	h_reim = d_reim;
	std::ofstream fou("demo4_reimg.raw", std::ios::binary);
	fou.write((char*) (&h_reim[0]), sizeof(float) * imgReso);
	fou.close();
}

//////////////////////////////////////////////////////////////////////////
// Projection and backprojection process
// The function combines A^{T}\cdot A together for Split-Bregman method
void ProjectionThenBackProjection(float* dimg, float* dres, const FanEAGeo& FanGeo, const Image& Img)
{
	float* dprj;
	checkCudaErrors(cudaMalloc((void**) &dprj, sizeof(float) * FanGeo.m_ViwN * FanGeo.m_DetN));
	dim3 prjB(256, 4);
	dim3 prjG(
		(FanGeo.m_DetN + prjB.x - 1) / prjB.x,
		(FanGeo.m_ViwN + prjB.y - 1) / prjB.y);
	dim3 bakB(32, 32);
	dim3 bakG(
		(Img.m_Reso.x + bakB.x - 1) / bakB.x,
		(Img.m_Reso.y + bakB.y - 1) / bakB.y);
	proj(dimg, dprj, FanGeo, Img, prjB, prjG);
	bakproj_PIXEL(dprj, dres, FanGeo, Img, bakB, bakG);

}

void ProjectionThenBackProjection(float* dimg, float* dres, const FanEDGeo& FanGeo, const Image& Img)
{
	float* dprj;
	checkCudaErrors(cudaMalloc((void**) &dprj, sizeof(float) * FanGeo.m_ViwN * FanGeo.m_DetN));
	dim3 prjB(256, 4);
	dim3 prjG(
		(FanGeo.m_DetN + prjB.x - 1) / prjB.x,
		(FanGeo.m_ViwN + prjB.y - 1) / prjB.y);
	dim3 bakB(32, 32);
	dim3 bakG(
		(Img.m_Reso.x + bakB.x - 1) / bakB.x,
		(Img.m_Reso.y + bakB.y - 1) / bakB.y);
	proj(dimg, dprj, FanGeo, Img, prjB, prjG);
	bakproj_PIXEL(dprj, dres, FanGeo, Img, bakB, bakG);
}







__global__ void demo5_minuss(float* orig, float* moti, const unsigned int length)
{
	unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x;
	if (idx < length)
	{
		moti[idx] = orig[idx] - moti[idx];
	}
}

__global__ void demo5_update(float* origImg, float* updImg, const unsigned int length, const float coeff)
{
	unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x;
	if (idx < length)
	{
		origImg[idx] = origImg[idx] + coeff * updImg[idx];
	}
}


//ÓÃÏßÐÔŽúÊýµÄ·œ·šÀŽ×öÖØ¹¹
void DEMO5()
{
	checkCudaErrors(cudaDeviceReset());
	cusparseHandle_t handle = 0;
	cusparseMatDescr_t descr = 0;
	const int detR = 1024;
	const int viwN = 360;
	const int imgR = 512;
	const int rowNum = detR * viwN;
	const int colNum = imgR * imgR;
	const int nnz = 208932896;

	std::ifstream csrRowFile("demo5_csrRowIdx.idx", std::ios::binary);
	std::ifstream csrColFile("demo5_csrColIdx.idx", std::ios::binary);
	std::ifstream csrValFile("demo5_csrValIdx.val", std::ios::binary);
	std::ifstream cscRowFile("demo5_cscRowIdx.idx", std::ios::binary);
	std::ifstream cscColFile("demo5_cscColIdx.idx", std::ios::binary);
	std::ifstream cscValFile("demo5_cscValIdx.val", std::ios::binary);

	if (!(csrRowFile.is_open() && csrColFile.is_open() && csrValFile.is_open()
		&& cscRowFile.is_open() && cscColFile.is_open() && cscValFile.is_open()))
	{
		std::cout << "Cannot open the files to demo5";
		exit(-1);
	}
	int* csrRowIdx;
	int* csrColIdx;
	int* cscRowIdx;
	int* cscColIdx;
	float* csrValIdx;
	float* cscValIdx;
	checkCudaErrors(cudaMallocManaged((void**) &csrRowIdx, (rowNum + 1) * sizeof(int)));
	checkCudaErrors(cudaMallocManaged((void**) &csrColIdx, nnz * sizeof(int)));
	checkCudaErrors(cudaMallocManaged((void**) &csrValIdx, nnz * sizeof(float)));
	checkCudaErrors(cudaMallocManaged((void**) &cscColIdx, (colNum + 1) * sizeof(int)));
	checkCudaErrors(cudaMallocManaged((void**) &cscRowIdx, nnz * sizeof(int)));
	checkCudaErrors(cudaMallocManaged((void**) &cscValIdx, nnz * sizeof(float)));

	csrRowFile.read((char*) csrRowIdx, sizeof(int) * (rowNum + 1));
	csrColFile.read((char*) csrColIdx, sizeof(int) * nnz);
	csrValFile.read((char*) csrValIdx, sizeof(float) *nnz);

	cscRowFile.read((char*) cscRowIdx, sizeof(int) * nnz);
	cscColFile.read((char*) cscColIdx, sizeof(int) * (colNum + 1));
	cscValFile.read((char*) cscValIdx, sizeof(float) *nnz);

	cscRowFile.close();
	cscColFile.close();
	cscValFile.close();
	cscRowFile.close();
	cscColFile.close();
	cscValFile.close();

	cusparseCreate(&handle);
	cusparseCreateMatDescr(&descr);
	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);


	//Establish the matrix
	checkCudaErrors(cudaDeviceSynchronize());


	float* image;
	float* prj;
	float* tempproj;
	float* tempimg;
	checkCudaErrors(cudaMallocManaged((void**) &image, sizeof(float) * colNum));
	checkCudaErrors(cudaMallocManaged((void**) &prj, sizeof(float) * rowNum));

	checkCudaErrors(cudaMallocManaged((void**) &tempimg, sizeof(float) * colNum));
	checkCudaErrors(cudaMallocManaged((void**) &tempproj, sizeof(float) * rowNum));

	checkCudaErrors(cudaMemsetAsync(image, 0, sizeof(float) * colNum));
	checkCudaErrors(cudaMemsetAsync(tempimg, 0, sizeof(float) * colNum));
	checkCudaErrors(cudaMemsetAsync(tempproj, 0, sizeof(float) * rowNum));
	checkCudaErrors(cudaDeviceSynchronize());

	std::ifstream prjID("demo5_prjj.prj", std::ios::binary);
	prjID.read((char*) prj, sizeof(float) * rowNum);
	prjID.close();


	float alpha = 1.0f;
	float beta = 0.0;

	dim3 updb(1024, 1);
	dim3 updg((rowNum + 1023) / 1024, 1);
	dim3 updb2(1024);
	dim3 updg2((colNum + 1023) / 1024);
	float stepSize = 0.000008f;
	float stepDes = 0.998f;
	const unsigned int iterN = 150;
	for (unsigned int i = 0; i != iterN; ++i)
	{
		//Projection

		cusparseScsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, rowNum, colNum, nnz, &alpha, descr, csrValIdx, csrRowIdx, csrColIdx, image, &beta, tempproj);
		//Differential;
		demo5_minuss << <updg, updb >> >(prj, tempproj, rowNum);

		//·ŽÍ¶Ó°;
		cusparseScsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, colNum, rowNum, nnz, &alpha, descr, cscValIdx, cscColIdx, cscRowIdx, tempproj, &beta, tempimg);


		//žüÐÂ;
		demo5_update << <updg2, updb2 >> >(image, tempimg, colNum, stepSize);

		if (i % 5 == 0)
		{
			stepSize *= stepDes;
		}
		std::cout << "iteration " << i << std::endl;
	}

	checkCudaErrors(cudaDeviceSynchronize());
	std::ofstream fou("demo5_reimg.raw", std::ios::binary);
	fou.write((char*) image, sizeof(float) * colNum);
	fou.close();
	checkCudaErrors(cudaDeviceReset());
}

void DEMO5_1()
{
	checkCudaErrors(cudaDeviceReset());

	cusparseHandle_t handle = 0;
	cusparseMatDescr_t descr = 0;


	/************************************************************************/
	/* Create the sparse matrix from COO format                             */
	/************************************************************************/
	const int detR = 1024;
	const int viwN = 360;
	const int imgR = 512;
	const int rowNum = detR * viwN;
	const int colNum = imgR * imgR;
	const int nnz = 208932896;

	int* cooRowIdx;
	int* cooColIdx;
	float* cooValIdx;

	checkCudaErrors(cudaMallocManaged((void**) &cooRowIdx, sizeof(int) * nnz));
	checkCudaErrors(cudaMallocManaged((void**) &cooColIdx, sizeof(int) * nnz));
	checkCudaErrors(cudaMallocManaged((void**) &cooValIdx, sizeof(float) *nnz));
	//Read the file
	std::string cooRowF = "demo5_rowIdx208932896.max";
	std::string cooColF = "demo5_colIdx208932896.max";
	std::string cooValF = "demo5_weight208932896.max";
	std::ifstream cooRowFile(cooRowF.c_str(), std::ios::binary);
	std::ifstream cooColFile(cooColF.c_str(), std::ios::binary);
	std::ifstream cooValFile(cooValF.c_str(), std::ios::binary);

	if (!(cooRowFile.is_open() && cooColFile.is_open() && cooValFile.is_open()))
	{
		std::cerr << "Cannot open the file\n";
		checkCudaErrors(cudaDeviceReset());
		exit(-1);
	}
	cooRowFile.read((char*) cooRowIdx, sizeof(int) * nnz);
	cooColFile.read((char*) cooColIdx, sizeof(int) * nnz);
	cooValFile.read((char*) cooValIdx, sizeof(float) *nnz);

	cusparseCreate(&handle);
	cusparseCreateMatDescr(&descr);
	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

	int* csrRowIdx;
	int* csrColIdx = cooColIdx;
	int* cscRowIdx;
	int* cscColIdx;
	float* cscValIdx;
	float* csrValIdx = cooValIdx;
	checkCudaErrors(cudaMallocManaged((void**) &csrRowIdx, (rowNum + 1) * sizeof(int)));
	checkCudaErrors(cudaMallocManaged((void**) &cscColIdx, (colNum + 1) * sizeof(int)));
	checkCudaErrors(cudaMallocManaged((void**) &cscRowIdx, nnz * sizeof(int)));
	checkCudaErrors(cudaMallocManaged((void**) &cscValIdx, nnz * sizeof(float)));
	cusparseXcoo2csr(handle, cooRowIdx, nnz, rowNum, csrRowIdx, CUSPARSE_INDEX_BASE_ZERO);
	cudaDeviceSynchronize();
	cusparseScsr2csc(handle, rowNum, colNum, nnz, csrValIdx, csrRowIdx, csrColIdx, cscValIdx, cscRowIdx, cscColIdx, CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO);
	checkCudaErrors(cudaDeviceSynchronize());


	float* image;
	float* proj;
	float* tempproj;
	float* tempimg;
	checkCudaErrors(cudaMallocManaged((void**) &image, sizeof(float) * colNum));
	checkCudaErrors(cudaMallocManaged((void**) &proj, sizeof(float) * rowNum));

	checkCudaErrors(cudaMallocManaged((void**) &tempimg, sizeof(float) * colNum));
	checkCudaErrors(cudaMallocManaged((void**) &tempproj, sizeof(float) * rowNum));

	checkCudaErrors(cudaMemsetAsync(image, 0, sizeof(float) * colNum));
	checkCudaErrors(cudaMemsetAsync(tempimg, 0, sizeof(float) * colNum));
	checkCudaErrors(cudaMemsetAsync(tempproj, 0, sizeof(float) * rowNum));
	checkCudaErrors(cudaDeviceSynchronize());

	std::ifstream prj("demo5_1_prjj.prj", std::ios::binary);
	prj.read((char*) proj, sizeof(float) * rowNum);
	prj.close();


	float alpha = 1.0f;
	float beta = 0.0;

	dim3 updb(1024, 1);
	dim3 updg((rowNum + 1023) / 1024, 1);
	dim3 updb2(1024);
	dim3 updg2((colNum + 1023) / 1024);
	float stepSize = 0.000008f;
	float stepDes = 0.998f;
	const unsigned int iterN = 80;
	for (unsigned int i = 0; i != iterN; ++i)
	{
		//Projection

		cusparseScsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, rowNum, colNum, nnz, &alpha, descr, csrValIdx, csrRowIdx, csrColIdx, image, &beta, tempproj);
		//Differential;
		demo5_minuss << <updg, updb >> >(proj, tempproj, rowNum);

		//·ŽÍ¶Ó°;
		cusparseScsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, colNum, rowNum, nnz, &alpha, descr, cscValIdx, cscColIdx, cscRowIdx, tempproj, &beta, tempimg);


		//žüÐÂ;
		demo5_update << <updg2, updb2 >> >(image, tempimg, colNum, stepSize);

		if (i % 5 == 0)
		{
			stepSize *= stepDes;
		}
		std::cout << "iteration " << i << std::endl;
	}


	checkCudaErrors(cudaDeviceSynchronize());
	std::ofstream fou("demo5_1_reimg.raw", std::ios::binary);
	fou.write((char*) image, sizeof(float) * colNum);
	fou.close();

	checkCudaErrors(cudaDeviceReset());
}





//This is the demo test the boxed back projection function in one angle
void DEMO7()
{
	checkCudaErrors(cudaDeviceReset());
	ConeEDGeo ConeGeo(85.0f, 15.0f, 360, 0.0f, _TWOPI, make_float2(34.0f, 34.0f), make_int2(1024, 1024));
	Volume Vol(512, 512, 512, 20.0f, 20.0f, 20.0f, 0.0f, 0.0f, 0.0f);


	size_t volSize = Vol.m_Reso.x * Vol.m_Reso.y * Vol.m_Reso.z;
	size_t prjSize = ConeGeo.m_DetN.x * ConeGeo.m_DetN.y;

	float* h_prj = new float[prjSize];
	float* h_vol = new float[volSize];

	std::string inFileName = "demo7.prj";
	std::ifstream fid(inFileName.c_str(), std::ios::binary);
	if (!(fid.is_open()))
	{
		std::cout << "Cannot open demo7.prj\n";
		exit(-1);
	}
	fid.read((char*) h_prj, sizeof(float) * prjSize);
	fid.close();

	float* d_prj;
	float* d_vol;
	checkCudaErrors(cudaMalloc((void**) &d_prj, sizeof(float) * prjSize));
	checkCudaErrors(cudaMalloc((void**) &d_vol, sizeof(float) * volSize));


	checkCudaErrors(cudaMemcpy(d_prj, h_prj, sizeof(float) * prjSize, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemset(d_vol, 0, sizeof(float) * volSize));
	dim3 blk(1024);
	dim3 gid((Vol.m_Reso.x * Vol.m_Reso.y * Vol.m_Reso.z + blk.x - 1) / blk.x);

	// demo7_allBack;
	bakproj_BOXED(d_prj, d_vol, 0, ConeGeo, Vol, blk, gid);

	checkCudaErrors(cudaMemcpy(h_vol, d_vol, sizeof(float) * volSize, cudaMemcpyDeviceToHost));
	std::ofstream fou("demo7_res.raw", std::ios::binary);
	fou.write((char*) h_vol, sizeof(float) * volSize);
	fou.close();
	cudaDeviceReset();
}

void DEMO7_1()
{
	checkCudaErrors(cudaDeviceReset());
	ConeEDGeo ConeGeo(85.0f, 15.0f, 360, 0.0f, _TWOPI, make_float2(34.0f, 3.4f), make_int2(512, 512));
	Volume Vol(512, 512, 512, 20.0f, 20.0f, 20.0f, 0.0f, 0.0f, 0.0f);


	size_t volSize = Vol.m_Reso.x * Vol.m_Reso.y * Vol.m_Reso.z;
	size_t prjSize = ConeGeo.m_DetN.x * ConeGeo.m_DetN.y * ConeGeo.m_ViwN;

	float* h_prj = new float[prjSize];
	float* h_vol = new float[volSize];

	std::string inFileName = "demo7_1.prj";
	std::ifstream fid(inFileName.c_str(), std::ios::binary);
	if (!(fid.is_open()))
	{
		std::cout << "Cannot open demo7_1.prj\n";
		exit(-1);
	}
	fid.read((char*) h_prj, sizeof(float) * prjSize);
	fid.close();

	float* d_prj;
	float* d_vol;
	checkCudaErrors(cudaMalloc((void**) &d_prj, sizeof(float) * prjSize));
	checkCudaErrors(cudaMalloc((void**) &d_vol, sizeof(float) * volSize));


	checkCudaErrors(cudaMemcpy(d_prj, h_prj, sizeof(float) * prjSize, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemset(d_vol, 0, sizeof(float) * volSize));
	dim3 blk(1024);
	dim3 gid((Vol.m_Reso.x * Vol.m_Reso.y * Vol.m_Reso.z + blk.x - 1) / blk.x);

	// demo7_allBack;
	bakproj_BOXED(d_prj, d_vol, ConeGeo, Vol, blk, gid);



	checkCudaErrors(cudaMemcpy(h_vol, d_vol, sizeof(float) * volSize, cudaMemcpyDeviceToHost));
	std::ofstream fou("demo7_1_res.raw", std::ios::binary);
	fou.write((char*) h_vol, sizeof(float) * volSize);
	fou.close();
	checkCudaErrors(cudaDeviceReset());
}











cudaArray *demo8_volarray = 0;
cudaArray *demo8_prjarray = 0;
texture<float, 3, cudaReadModeElementType> demo8_volTex;
texture<float, 3, cudaReadModeElementType> demo8_prjTex;
template<typename T>
__global__ void proj_demo8(T* d_output, const ConeEDGeo ConeGeo, const Volume Vol)
{

	cuint x = blockIdx.x * blockDim.x + threadIdx.x;
	cuint y = blockIdx.y * blockDim.y + threadIdx.y;
	cuint angIdx = blockIdx.z * blockDim.z + threadIdx.z;
	if (x >= ConeGeo.m_DetN.x || y >= ConeGeo.m_DetN.y || angIdx >= ConeGeo.m_ViwN) return;


	const unsigned int maxSteps = (Vol.m_Reso.x * 10); //need modified

	const float tstep = (Vol.m_Step.x / 8.0); //need modified
	const float3 boxMin = make_float3(
		-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
		-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
		-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
	const float3 boxMax = make_float3(
		Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
		Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
		Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
	float cosT = cosf(ConeGeo.m_ViwBeg + ConeGeo.m_ViwStp * angIdx);
	float sinT = sinf(ConeGeo.m_ViwBeg + ConeGeo.m_ViwStp * angIdx);
	Ray eyeRay;
	eyeRay.o.x = -ConeGeo.m_S2O * sinT;
	eyeRay.o.y = ConeGeo.m_S2O  * cosT;
	eyeRay.o.z = 0.0;
	eyeRay.d = normalize(rotation(make_float3(
		(x - ConeGeo.m_DetCntIdx.x + 0.5f) * ConeGeo.m_DetStp.x,
		-ConeGeo.m_O2D,
		(y - ConeGeo.m_DetCntIdx.y + 0.5f) * ConeGeo.m_DetStp.y), cosT, sinT) - eyeRay.o);

	float tnear, tfar;
	int hit = intersectBox(eyeRay, boxMin, boxMax, &tnear, &tfar);
	if (!hit) return;

	if (tnear < 0.0f) tnear = 0.0f;

	//march along ray from front to back, accumulating color
	float sum(0.0f);
	float t = tnear;
	float3 pos = eyeRay.o + eyeRay.d * tnear;
	float3 step = eyeRay.d * tstep;
	float sample(0.0f);
	for (int i = 0; i != maxSteps; ++i)
	{
		sample = tex3D(demo8_volTex,
			(pos.x * 0.5f + 0.5f) * Vol.m_Reso.x,
			(pos.y * 0.5f + 0.5f) * Vol.m_Reso.y,
			(pos.z * 0.5f + 0.5f) * Vol.m_Reso.z);
		sum += sample * tstep;
		t += tstep;
		if (t > tfar) break;
		pos += step;
	}
	d_output[(angIdx * ConeGeo.m_DetN.y + y) * ConeGeo.m_DetN.x + x] = tex3D(demo8_prjTex, x, y, angIdx) - sum;
}

void DEMO8()
{

	checkCudaErrors(cudaDeviceReset());
	ConeEDGeo ConeGeo;
	Volume Vol;

	cudaExtent volumeSize = make_cudaExtent(Vol.m_Reso.x, Vol.m_Reso.y, Vol.m_Reso.z);
	size_t size = volumeSize.width * volumeSize.height * volumeSize.depth;

	float* h_volume = new float[size];
	std::string inFileName = "demo8.raw";
	std::ifstream fid(inFileName.c_str(), std::ios::binary);
	if (!fid.is_open())
	{
		std::cout << "Cannot open demo8.prj for projection test with texture\n";
		exit(-1);
	}
	fid.read((char*) h_volume, sizeof(float) *size);
	fid.close();

	cudaExtent projSize = make_cudaExtent(ConeGeo.m_DetN.x, ConeGeo.m_DetN.y, ConeGeo.m_ViwN);
	size_t prjS = projSize.width * projSize.height * projSize.depth;
	float* h_proj = new float[prjS];
	for (unsigned int i = 0; i != prjS; ++i)
	{
		h_proj[i] = 0;
	}


	//create 3D array
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	checkCudaErrors(cudaMalloc3DArray(&demo8_volarray, &channelDesc, volumeSize));

	// copy data to 3D array
	cudaMemcpy3DParms copyParams = { 0 };
	copyParams.srcPtr = make_cudaPitchedPtr(h_volume, volumeSize.width*sizeof(float), volumeSize.width, volumeSize.height);
	copyParams.dstArray = demo8_volarray;
	copyParams.extent = volumeSize;
	copyParams.kind = cudaMemcpyHostToDevice;
	checkCudaErrors(cudaMemcpy3D(&copyParams));

	// set texture parameters
	demo8_volTex.normalized = false;                      // access with normalized texture coordinates
	demo8_volTex.filterMode = cudaFilterModeLinear;      // linear interpolation
	demo8_volTex.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
	demo8_volTex.addressMode[1] = cudaAddressModeClamp;

	// bind array to 3D texture
	checkCudaErrors(cudaBindTextureToArray(demo8_volTex, demo8_volarray, channelDesc));

	//create 3D array
	cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<float>();
	checkCudaErrors(cudaMalloc3DArray(&demo8_prjarray, &channelDesc2, projSize));

	// copy data to 3D array
	cudaMemcpy3DParms copyParams2 = { 0 };
	copyParams2.srcPtr = make_cudaPitchedPtr(h_proj, projSize.width*sizeof(float), projSize.width, projSize.height);
	copyParams2.dstArray = demo8_prjarray;
	copyParams2.extent = projSize;
	copyParams2.kind = cudaMemcpyHostToDevice;
	checkCudaErrors(cudaMemcpy3D(&copyParams2));

	// set texture parameters
	demo8_prjTex.normalized = false;                      // access with normalized texture coordinates
	demo8_prjTex.filterMode = cudaFilterModeLinear;      // linear interpolation
	demo8_prjTex.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
	demo8_prjTex.addressMode[1] = cudaAddressModeClamp;

	// bind array to 3D texture
	checkCudaErrors(cudaBindTextureToArray(demo8_prjTex, demo8_prjarray, channelDesc2));

	//initCuda((void*)h_volume, volumeSize);
	delete [] h_volume;
	delete [] h_proj;

	//bool bTestResult = true;
	dim3 blockSize(16, 16, 4);
	dim3 gridSize(
		(ConeGeo.m_DetN.x + 15) / 16,
		(ConeGeo.m_DetN.y + 15) / 16,
		(ConeGeo.m_ViwN + 3) / 4);

	float* d_output = nullptr;
	//float* d_depth;
	checkCudaErrors(cudaMalloc((void **) &d_output, prjS * sizeof(float)));
	checkCudaErrors(cudaMemset(d_output, 0, prjS * sizeof(float)));

	proj_demo8 << <gridSize, blockSize >> >(d_output, ConeGeo, Vol);
	checkCudaErrors(cudaDeviceSynchronize());
	// Get elapsed time and throughput, then log to sample and master logs

	//getLastCudaError("Error: render_kernel() execution FAILED");
	checkCudaErrors(cudaDeviceSynchronize());

	float *h_output = new float[prjS];
	checkCudaErrors(cudaMemcpy(h_output, d_output, sizeof(float) * prjS, cudaMemcpyDeviceToHost));

	checkCudaErrors(cudaFree(d_output));

	std::ofstream fou("demo8_1024x512x360.prj", std::ios::binary);
	fou.write((char*) h_output, sizeof(float) * prjS);
	fou.close();
	free(h_output);
	std::cout << "Finished Projection test\n";
}




void CG_recon(float* himg, float* hprj, const FanEAGeo& FanGeo, const Image& Img, cuint IterNum)
{
	float* X = nullptr;
	float* realp = nullptr;
	float* R = nullptr;
	float* D = nullptr;
	float* tem = nullptr;
	float* tem2 = nullptr;
	float* nR = nullptr;
	cuint imgReso = Img.m_Reso.x * Img.m_Reso.y;
	cuint prjReso = FanGeo.m_DetN * FanGeo.m_ViwN;

	dim3 prjblk(128, 16);
	dim3 prjgid(
		(FanGeo.m_DetN + prjblk.x - 1) / prjblk.x,
		(FanGeo.m_ViwN + prjblk.y - 1) / prjblk.y);
	dim3 bakblk(32, 32);
	dim3 bakgid(
		(Img.m_Reso.x + bakblk.x - 1) / bakblk.x,
		(Img.m_Reso.y + bakblk.y - 1) / bakblk.y);



	checkCudaErrors(cudaMalloc((void**) &X, sizeof(float) * imgReso));
	checkCudaErrors(cudaMalloc((void**) &realp, sizeof(float) * prjReso));
	checkCudaErrors(cudaMalloc((void**) &R, sizeof(float) * imgReso));
	checkCudaErrors(cudaMalloc((void**) &D, sizeof(float) * imgReso));
	checkCudaErrors(cudaMalloc((void**) &tem, sizeof(float) * prjReso));
	checkCudaErrors(cudaMalloc((void**) &tem2, sizeof(float) * imgReso));
	checkCudaErrors(cudaMalloc((void**) &nR, sizeof(float) * imgReso));

	checkCudaErrors(cudaMemcpy(realp, hprj, sizeof(float) * prjReso, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemset(X, 0, sizeof(float) * imgReso));
	bakproj_PIXEL(realp, R, FanGeo, Img, bakblk, bakgid);
	checkCudaErrors(cudaMemcpy(D, R, sizeof(float) * imgReso, cudaMemcpyDeviceToDevice));
	unsigned int iter(0);
	float alpha(0.0f), beta(0.0f);
	thrust::device_ptr<float> pX(X);
	thrust::device_ptr<float> pD(D);
	thrust::device_ptr<float> pR(R);
	thrust::device_ptr<float> ptem2(tem2);
	thrust::device_ptr<float> pnR(nR);
	float* tempfou = new float[imgReso];
	float norm1(0.0f), norm2(0.0f), norm3(0.0f), norm4(0.0f);
	while (iter != IterNum)
	{
		++iter;
		proj(D, tem, FanGeo, Img, prjblk, prjgid);
		norm1 = norm(R, imgReso, 2.0f);
		norm1 = norm1 * norm1;
		norm2 = norm(tem, imgReso, 2.0f);
		norm2 = norm2 * norm2;
		alpha = norm1 / norm2;
		thrust::transform(pX, pX + imgReso, pD, pX, CG_update1_functor<float>(alpha));
		//bakproj_BOXED(tem, tem2, FanGeo, Img, bakblk, bakgid);
		bakproj_PIXEL(tem, tem2, FanGeo, Img, bakblk, bakgid);
		thrust::transform(pR, pR + imgReso, ptem2, pnR, CG_update2_functor<float>(alpha));
		norm3 = norm(nR, imgReso, 2.0f);
		norm4 = norm(R, imgReso, 2.0f);
		norm3 = norm3 * norm3;
		norm4 = norm4 * norm4;
		beta = norm3 / norm4;
		checkCudaErrors(cudaMemcpy(R, nR, sizeof(float) * imgReso, cudaMemcpyDeviceToDevice));
		thrust::transform(pR, pR + imgReso, pD, pD, CG_update3_functor<float>(beta));
		checkCudaErrors(cudaMemcpy(tempfou, X, sizeof(float) * imgReso, cudaMemcpyDeviceToHost));
		std::stringstream ss;
		ss << iter;
		std::string fil = "demo12_iterRecon" + ss.str() + ".raw";
		std::ofstream fou(fil.c_str(), std::ios::binary);
		fou.write((char*) tempfou, sizeof(float) * imgReso);
		fou.close();
		std::cout << iter << std::endl;
	}
}






void DEMO9_1()
{
	//Reconstruct the Image;
	FanEAGeo FanGeo;
	Image Img;

	FanGeo.m_DetArc = 0.95928517242269f;
	FanGeo.m_DetN = 888;
	FanGeo.m_DetCntIdx = 444.75f;
	FanGeo.m_DetStp = FanGeo.m_DetArc / FanGeo.m_DetN;
	FanGeo.m_O2D = static_cast<float>(4.082259521484375e+02);

	FanGeo.m_S2O = static_cast<float>(5.385200195312500e+02);
	FanGeo.m_S2D = FanGeo.m_S2O + FanGeo.m_O2D;
	FanGeo.m_ViwBeg = 0.0f;// (-2.668082275390625e+02 / 180.0* 3.14159265358979323846264);
	FanGeo.m_ViwN = 2200;
	FanGeo.m_ViwStp = _TWOPI / FanGeo.m_ViwN;

	Img.m_Bias.x = 0.0f;
	Img.m_Bias.y = 0.0f;  //ÕâžöÆ«ÒÆµÄµ¥Î»ÊÇÕæÊµÎïÀíµ¥Î»;
	Img.m_Reso.x = Img.m_Reso.y = 512;
	Img.m_Size.x = Img.m_Size.y = static_cast<float>(4.484740011196460e+02);
	Img.m_Step.x = Img.m_Size.x / Img.m_Reso.x;
	Img.m_Step.y = Img.m_Size.y / Img.m_Reso.y;

	thrust::host_vector<float> himg(Img.m_Reso.x *Img.m_Reso.y, 0);
	thrust::host_vector<float> hprj(FanGeo.m_DetN * FanGeo.m_ViwN, 0);
	thrust::host_vector<float> hweg(Img.m_Reso.x *Img.m_Reso.y, 0);
	thrust::host_vector<float> hmsk(Img.m_Reso.x *Img.m_Reso.y, 1.0f);
	std::string FileName = "demo9.prj";
	std::string FouName = "demo9ReconImg.raw";
	std::string WegName = "demo9_bakwegEA.raw";
	//Read the data;
	std::ifstream fid(FileName.c_str(), std::ios::binary);
	std::ifstream fons(WegName.c_str(), std::ios::binary);
	fid.read((char*) &hprj[0], sizeof(float) *FanGeo.m_DetN * FanGeo.m_ViwN);
	fons.read((char*) &hweg[0], sizeof(float) *Img.m_Reso.x *Img.m_Reso.y);
	fid.close();
	fons.close();
	SART(himg, hprj, hweg, hmsk, FanGeo, Img, 20, 1.0f);

	std::ofstream fou(FouName.c_str(), std::ios::binary);
	fou.write((char*) &himg[0], sizeof(float) * Img.m_Reso.x * Img.m_Reso.y);
	fou.close();
}

void DEMO9_2()
{
	//Reconstruct the Image;
	FanEDGeo FanGeo(85.0f, 15.0f, 2200, 0.0f, 3.1415926f*2.0f, 45.0f, 888);
	Image Img(512, 512, 20.0f, 20.0f, 0.0f, 0.0f);

	thrust::host_vector<float> himg(Img.m_Reso.x *Img.m_Reso.y, 0);
	thrust::host_vector<float> hprj(FanGeo.m_DetN * FanGeo.m_ViwN, 0);
	thrust::host_vector<float> hweg(Img.m_Reso.x *Img.m_Reso.y, 0);
	thrust::host_vector<float> hmsk(Img.m_Reso.x *Img.m_Reso.y, 1.0f);
	std::string FileName = "demo9.prj";
	std::string FouName = "demo9ReconImg.raw";
	std::string WegName = "demo9_bakwegEA.raw";
	//Read the data;
	std::ifstream fid(FileName.c_str(), std::ios::binary);
	std::ifstream fons(WegName.c_str(), std::ios::binary);
	fid.read((char*) &hprj[0], sizeof(float) *FanGeo.m_DetN * FanGeo.m_ViwN);
	fons.read((char*) &hweg[0], sizeof(float) *Img.m_Reso.x *Img.m_Reso.y);
	fid.close();
	fons.close();
	SART(himg, hprj, hweg, hmsk, FanGeo, Img, 20, 1.0f);

	std::ofstream fou(FouName.c_str(), std::ios::binary);
	fou.write((char*) &himg[0], sizeof(float) * Img.m_Reso.x * Img.m_Reso.y);
	fou.close();

}



void DEMO10()
{
	//The configuration of the IMAGING system;
	FanEAGeo FanGeo(
		static_cast<float>(5.385200195312500e+02), // S2O
		static_cast<float>(4.082259521484375e+02), // O2D
		2200,				   // viwNum
		0.0f,				   // Viw Begin
		_TWOPI,   // Viw End
		0.95928517242269f,      // Detector Arc,
		888);    			   // Detector Element Num
	FanGeo.m_DetCntIdx = 444.75f;
	Image Img(512, 512, 500.0f, 500.0f, 0.0f, 0.0f);
	cuint sliceNum = 64; //That we need to reconstruct the Volume with 64 slices
	cuint sliceReso = Img.m_Reso.x * Img.m_Reso.y;
	cuint totalVolumeBytes = sliceReso * sliceNum;
	cuint projReso = FanGeo.m_DetN * FanGeo.m_ViwN;
	cuint totalPrjectBytes = projReso * sliceNum;


	thrust::host_vector<float> himg(totalVolumeBytes, 0);
	thrust::host_vector<float> hweg(sliceReso, 0);
	thrust::host_vector<float> hprj(totalPrjectBytes, 0);
	thrust::host_vector<float> hmsk(sliceReso, 0);

	std::string ProjFileName = "demo10_proj.prj";
	std::string ReconVolName = "demo10_recon.raw";
	std::string WeightName = "demo10_weight.raw";
	std::string MaskName = "demo10_mask.img";

	std::ifstream fmsk(MaskName.c_str(), std::ios::binary);
	std::ifstream fprj(ProjFileName.c_str(), std::ios::binary);
	std::ifstream fweg(WeightName.c_str(), std::ios::binary);
	if (!(fmsk.is_open() && fprj.is_open() && fweg.is_open()))
	{
		std::cerr << "Cannot open the file for demo 10\n";
		exit(-1);
	}
	fmsk.read((char*) &hmsk[0], sizeof(float) * sliceReso);
	fprj.read((char*) &hprj[0], sizeof(float) * totalPrjectBytes);
	fweg.read((char*) &hweg[0], sizeof(float) * sliceReso);
	fmsk.close(); fprj.close(); fweg.close();

	cuint iterNum = 6;
	const float lambda = 1.8f;


	SART_MULTISLICES(himg, hprj, hweg, hmsk, FanGeo, Img, iterNum, lambda, sliceNum);


	std::ofstream fou(ReconVolName.c_str(), std::ios::binary);
	if (!fou.is_open())
	{
		std::cout << "Cannot open the file for write in demo 10\n";
		exit(-1);
	}
	fou.write((char*) &himg[0], sizeof(float) * totalVolumeBytes);
	fou.close();


}




__global__ void updateImg_demo11(float* dimg, float* dcorImg, float* dwegImg, float* dmask, float lambda, cuint ImgLR, cuint ImgWR, cuint sliceNum)
{
	cuint l = threadIdx.x + blockIdx.x * blockDim.x;
	cuint w = threadIdx.y + blockIdx.y * blockDim.y;
	cuint s = threadIdx.z + blockIdx.z * blockDim.z;
	if (l < ImgLR && w < ImgWR && s < sliceNum)
	{
		cuint i = (s * ImgWR + w) * ImgLR + l;
		dimg[i] += lambda * dcorImg[i] / dwegImg[w * ImgLR + l] * dmask[w * ImgLR + l];
		if (dimg[i] < 0)
		{
			dimg[i] = 0;
		}
	}
}
void DEMO11()
{
	FanEAGeo FanGeo(
		static_cast<float>(5.385200195312500e+02),
		static_cast<float>(4.082259521484375e+02),
		2200,
		0.0f,
		_TWOPI,
		0.95928517242269f,
		888);
	FanGeo.m_DetCntIdx = 444.75f;
	Image Img(512, 512, 500.0f, 500.0f, 0.0f, 0.0f);
	cuint sliceNum = 22; //Every GPU reconstruct 22 slices

	cuint subSetNum = 50; //œ«Í¶Ó°·Ö³É50žösubset
	cuint numPerSubSet = FanGeo.m_ViwN / subSetNum; //Ã¿žö×ÓŒ¯ÖÐÓÐ44žöÍ¶Ó°;

	cuint sliceReso = Img.m_Reso.x * Img.m_Reso.y;
	cuint projReso = FanGeo.m_DetN * FanGeo.m_ViwN;
	cuint totalVolBytes = sliceReso * sliceNum;
	cuint totalPrjBytes = projReso * sliceNum;

	float** img = new float*[3];
	float** weg = new float*[3];
	float** prj = new float*[3];
	float** mask = new float*[3];
	for (unsigned int i = 0; i != 3; ++i)
	{
		img[i] = new float[totalVolBytes];
		weg[i] = new float[sliceReso];
		prj[i] = new float[totalPrjBytes];
		mask[i] = new float[sliceReso];
	}
	float* totvol = new float[totalVolBytes * 3];

	std::string FileName1 = "demo11_partprj1.prj";
	std::string FileName2 = "demo11_partprj2.prj";
	std::string FileName3 = "demo11_partprj3.prj";

	std::string FouName = "demo11_reconImg512x512x66.raw";
	std::string WegName = "demo11_bakwegEA.raw";
	//Mask
	std::string maskImg = "demo11_newmask.img";

	//Read data
	std::ifstream fid1(FileName1.c_str(), std::ios::binary);
	std::ifstream fid2(FileName2.c_str(), std::ios::binary);
	std::ifstream fid3(FileName3.c_str(), std::ios::binary);
	std::ifstream fweg(WegName.c_str(), std::ios::binary);
	std::ifstream fmsk(maskImg.c_str(), std::ios::binary);
	if (!(fid1.is_open() && fid2.is_open() && fid3.is_open() && fweg.is_open() && fmsk.is_open()))
	{
		std::cout << "Cannot open all the files for demo11\n";
		exit(-1);
	}
	fid1.read((char*) prj[0], sizeof(float) * totalPrjBytes);
	fid2.read((char*) prj[1], sizeof(float) * totalPrjBytes);
	fid3.read((char*) prj[2], sizeof(float) * totalPrjBytes);

	fweg.read((char*) weg[0], sizeof(float) * sliceReso);
	std::memcpy(weg[1], weg[0], sizeof(float) * sliceReso);
	std::memcpy(weg[2], weg[0], sizeof(float) * sliceReso);

	fmsk.read((char*) mask[0], sizeof(float) * sliceReso);
	std::memcpy(mask[1], mask[0], sizeof(float) * sliceReso);
	std::memcpy(mask[2], mask[0], sizeof(float) * sliceReso);

	fid1.close();
	fid2.close();
	fid3.close();
	fweg.close();
	fmsk.close();

	cudaStream_t streams[3]; // Three streams corresponding to three GPUs

	//Corresponding to three GPUs
	float* dimg[3];
	float* dprj[3];
	float* dimgWeg[3];
	float* dcorImg[3];
	float* dcorPrj[3];
	float* dmask[3];
	float* dlas[3];


	//Begin SART-algorithm
	unsigned int i = 0;

	//Set the block and grid size
	dim3 prjBlk(128, 8);
	dim3 prjGid(
		(FanGeo.m_DetN + 127) / 128, //X dimension for detIdx
		(numPerSubSet + 7) / 8); // Y dimension for viwIdx

	dim3 bakBlk(32, 32);
	dim3 bakGid(
		(Img.m_Reso.x + 31) / 32, //X dimension for imgXidx,
		(Img.m_Reso.y + 31) / 32); // Y dimension for imgYidx
	dim3 updBlk(16, 16, 4);
	dim3 updGid(
		(Img.m_Reso.x + 15) / 16,
		(Img.m_Reso.y + 15) / 16,
		(sliceNum + 3) / 4);

	unsigned int iters = 0;
	const float lambda = 2.0; //updating coefficient
	const unsigned int iterNum = 40; //Iteration number
	float curLambda = lambda;
	float t1 = 1; //We would like to use the FISTA acceleration
	float t2 = 1;

	unsigned int cursubIdx = 0;
	unsigned int realIters = subSetNum * iterNum;

	//allocate the memory on devices
	for (i = 0; i != 3; ++i)
	{
		checkCudaErrors(cudaSetDevice(i));
		checkCudaErrors(cudaStreamCreate(&streams[i]));

		//·ÖÅäÏÔŽæ;
		checkCudaErrors(cudaMalloc((void**) &dimg[i], sizeof(float) * totalVolBytes));
		checkCudaErrors(cudaMalloc((void**) &dprj[i], sizeof(float) * totalPrjBytes));
		checkCudaErrors(cudaMalloc((void**) &dimgWeg[i], sizeof(float) * sliceReso));
		checkCudaErrors(cudaMalloc((void**) &dcorImg[i], sizeof(float) * totalVolBytes));
		checkCudaErrors(cudaMalloc((void**) &dcorPrj[i], sizeof(float) * FanGeo.m_DetN * numPerSubSet * sliceNum));
		checkCudaErrors(cudaMalloc((void**) &dmask[i], sizeof(float) * sliceReso));
		checkCudaErrors(cudaMalloc((void**) &dlas[i], sizeof(float) * totalVolBytes));

		checkCudaErrors(cudaMemset(dimg[i], 0, sizeof(float) * totalVolBytes));
		checkCudaErrors(cudaMemcpy(dprj[i], prj[i], sizeof(float) * totalPrjBytes, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(dimgWeg[i], weg[i], sizeof(float) * sliceReso, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemset(dcorImg[i], 0, sizeof(float) * totalVolBytes));
		checkCudaErrors(cudaMemset(dcorPrj[i], 0, sizeof(float) * FanGeo.m_DetN * numPerSubSet  * sliceNum));
		checkCudaErrors(cudaMemcpy(dmask[i], mask[i], sizeof(float) * sliceReso, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemset(dlas[i], 0, sizeof(float) * totalVolBytes));
	}
	thrust::device_ptr<float> pimg0(dimg[0]);
	thrust::device_ptr<float> pimg1(dimg[1]);
	thrust::device_ptr<float> pimg2(dimg[2]);
	thrust::device_ptr<float> plas0(dlas[0]);
	thrust::device_ptr<float> plas1(dlas[1]);
	thrust::device_ptr<float> plas2(dlas[2]);

	while (iters != realIters)
	{
		curLambda *= 0.995f;
		if (curLambda < 1.0)
		{
			curLambda = 1.0;
		}
		t2 = (1 + sqrtf(1 + 4.0f * t1 * t1)) * 0.5f;

		checkCudaErrors(cudaSetDevice(0));
		checkCudaErrors(cudaMemset(dcorImg[0], 0, sizeof(float) * totalVolBytes));
		checkCudaErrors(cudaMemcpy(dlas[0], dimg[0], sizeof(float) * totalVolBytes, cudaMemcpyDeviceToDevice));
		proj(dimg[0], dcorPrj[0], dprj[0], FanGeo, Img, numPerSubSet, subSetNum, cursubIdx, prjBlk, prjGid, sliceNum, streams[0]);
		bakproj_PIXEL(dcorPrj[0], dcorImg[0], FanGeo, Img, numPerSubSet, subSetNum, cursubIdx, bakBlk, bakGid, sliceNum, streams[0]);
		updateImg_demo11 << <updGid, updBlk, 0, streams[0] >> >(dimg[0], dcorImg[0], dimgWeg[0], dmask[0], curLambda * subSetNum, Img.m_Reso.x, Img.m_Reso.y, sliceNum);
		thrust::transform(pimg0, pimg0 + totalVolBytes, plas0, pimg0, _FISTA_demo11<float>(t1, t2));

		checkCudaErrors(cudaSetDevice(1));
		checkCudaErrors(cudaMemset(dcorImg[1], 0, sizeof(float) *totalVolBytes));
		checkCudaErrors(cudaMemcpy(dlas[1], dimg[1], sizeof(float) * totalVolBytes, cudaMemcpyDeviceToDevice));
		proj(dimg[1], dcorPrj[1], dprj[1], FanGeo, Img, numPerSubSet, subSetNum, cursubIdx, prjBlk, prjGid, sliceNum, streams[1]);
		bakproj_PIXEL(dcorPrj[1], dcorImg[1], FanGeo, Img, numPerSubSet, subSetNum, cursubIdx, bakBlk, bakGid, sliceNum, streams[1]);
		updateImg_demo11 << <updGid, updBlk, 0, streams[1] >> >(dimg[1], dcorImg[1], dimgWeg[1], dmask[1], curLambda * subSetNum, Img.m_Reso.x, Img.m_Reso.y, sliceNum);
		thrust::transform(pimg1, pimg1 + totalVolBytes, plas1, pimg1, _FISTA_demo11<float>(t1, t2));


		checkCudaErrors(cudaSetDevice(2));
		checkCudaErrors(cudaMemset(dcorImg[2], 0, sizeof(float) * totalVolBytes));
		checkCudaErrors(cudaMemcpy(dlas[2], dimg[2], sizeof(float) * totalVolBytes, cudaMemcpyDeviceToDevice));
		proj(dimg[2], dcorPrj[2], dprj[2], FanGeo, Img, numPerSubSet, subSetNum, cursubIdx, prjBlk, prjGid, sliceNum, streams[2]);
		bakproj_PIXEL(dcorPrj[2], dcorImg[2], FanGeo, Img, numPerSubSet, subSetNum, cursubIdx, bakBlk, bakGid, sliceNum, streams[2]);
		updateImg_demo11 << <updGid, updBlk, 0, streams[2] >> >(dimg[2], dcorImg[2], dimgWeg[2], dmask[2], curLambda * subSetNum, Img.m_Reso.x, Img.m_Reso.y, sliceNum);
		thrust::transform(pimg2, pimg2 + totalVolBytes, plas2, pimg2, _FISTA_demo11<float>(t1, t2));



		cudaDeviceSynchronize();
		t1 = t2;
		std::cout << "Iters: " << iters << std::endl;
		cursubIdx = (cursubIdx + 1) % subSetNum;
		++iters;
	}
	checkCudaErrors(cudaSetDevice(0));
	checkCudaErrors(cudaMemcpy(totvol, dimg[0],
		sizeof(float) * totalVolBytes, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaSetDevice(1));
	checkCudaErrors(cudaMemcpy(totvol + totalVolBytes, dimg[1],
		sizeof(float) * totalVolBytes, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaSetDevice(2));
	checkCudaErrors(cudaMemcpy(totvol + totalVolBytes * 2, dimg[2],
		sizeof(float) * totalVolBytes, cudaMemcpyDeviceToHost));

	std::ofstream fou(FouName.c_str(), std::ios::binary);
	if (!fou.is_open())
	{
		std::cerr << "Cannot open the file for demo11\n";
		exit(-1);
	}

	fou.write((char*) totvol, sizeof(float) * totalVolBytes * 3);
	fou.close();
}



void DEMO12()
{
	FanEAGeo FanGeo(
		static_cast<float>(5.385200195312500e+02),
		static_cast<float>(4.082259521484375e+02),
		2200,
		0.0f,
		_TWOPI,
		0.95928517242269f,
		888);
	FanGeo.m_DetCntIdx = 444.75f;
	Image Img(512, 512, 500.0f, 500.0f, 0.0f, 0.0f);

	std::ifstream fprj("demo12_proj.prj", std::ios::binary);
	std::ofstream frec("demo12_recon.raw", std::ios::binary);
	if (!(fprj.is_open() && frec.is_open()))
	{
		std::cerr << "Cannot open the file for demo 12\n";
		exit(-1);
	}

	cuint imgReso = Img.m_Reso.x * Img.m_Reso.y;
	cuint prjReso = FanGeo.m_DetN * FanGeo.m_ViwN;
	cuint iterNum = 40;
	float* himg = new float[imgReso];
	float* hprj = new float[prjReso];

	fprj.read((char*) hprj, sizeof(float) * prjReso);

	CG_recon(himg, hprj, FanGeo, Img, iterNum);

	frec.write((char*) himg, sizeof(float) * imgReso);

	fprj.close(); frec.close();
}




void EM(float* hprj, float* himg, const FanEAGeo& FanGeo, const Image& Img, cuint iterNum)
{
	dim3 prjBlk(128, 8);
	dim3 prjGid(
		(FanGeo.m_DetN + prjBlk.x - 1) / prjBlk.x,
		(FanGeo.m_ViwN + prjBlk.y - 1) / prjBlk.y);
	dim3 bakBlk(32, 32);
	dim3 bakGid(
		(Img.m_Reso.x + bakBlk.x - 1) / bakBlk.x,
		(Img.m_Reso.y + bakBlk.y - 1) / bakBlk.y);


	float* dprj;
	float* dimg;
	float* dta;
	float* dtb;
	float* LAM;
	cuint prjReso = FanGeo.m_DetN * FanGeo.m_ViwN;
	cuint imgReso = Img.m_Reso.x * Img.m_Reso.y;
	checkCudaErrors(cudaMalloc((void**) &dprj, sizeof(float) * prjReso));
	checkCudaErrors(cudaMalloc((void**) &dimg, sizeof(float) * imgReso));
	checkCudaErrors(cudaMalloc((void**) &dta, sizeof(float) * prjReso));
	checkCudaErrors(cudaMalloc((void**) &dtb, sizeof(float) * imgReso));
	checkCudaErrors(cudaMalloc((void**) &LAM, sizeof(float) * imgReso));
	thrust::device_ptr<float> pprj(dprj);
	thrust::device_ptr<float> pta(dta);
	thrust::device_ptr<float> pimg(dimg);
	thrust::device_ptr<float> ptb(dtb);
	thrust::device_ptr<float> pLAM(LAM);

	checkCudaErrors(cudaMemcpy(dprj, hprj, sizeof(float) * prjReso, cudaMemcpyHostToDevice));
	/*checkCudaErrors(cudaMemset(dimg,1,sizeof(float) * imgReso));
	checkCudaErrors(cudaMemset(dta,1,sizeof(float) * prjReso));
	checkCudaErrors(cudaMemset(dtb,0,sizeof(float) * imgReso));
	checkCudaErrors(cudaMemset(LAM,0,sizeof(float) * imgReso));*/

	thrust::fill(pimg, pimg + imgReso, 1.0f);
	thrust::fill(pta, pta + prjReso, 1.0f);
	thrust::fill(ptb, ptb + imgReso, 0.0f);
	thrust::fill(pLAM, pLAM + imgReso, 0.0f);


	bakproj_PIXEL(dta, LAM, FanGeo, Img, bakBlk, bakGid);

#if _DEBUG
	float* hta = new float[prjReso];
#endif
	for (unsigned int iter = 0; iter != iterNum; ++iter)
	{
		//projection
		proj(dimg, dta, FanGeo, Img, prjBlk, prjGid);

		thrust::transform(pprj, pprj + prjReso, pta, pta, _divide_functor<float>());
#if _DEBUG
		cudaMemcpy(hta, dta, sizeof(float) * prjReso, cudaMemcpyDeviceToHost);
		std::stringstream ss;
		ss << iter;
		std::string File = "hta" + ss.str() + ".prj";
		std::ofstream Fou(File.c_str(), std::ios::binary);
		Fou.write((char*) hta, sizeof(float) * prjReso);
		Fou.close();
#endif
		cudaMemset(dtb, 0, sizeof(float) * imgReso);
		bakproj_PIXEL(dta, dtb, FanGeo, Img, bakBlk, bakGid);
		thrust::transform(pimg, pimg + imgReso, ptb, pimg, thrust::multiplies<float>());
		thrust::transform(pimg, pimg + imgReso, pLAM, pimg, _divide_functor<float>());

	}

	checkCudaErrors(cudaMemcpy(himg, dimg, sizeof(float) * imgReso, cudaMemcpyDeviceToHost));
}


void DEMO13()
{
	FanEAGeo FanGeo(
		static_cast<float>(5.385200195312500e+02),
		static_cast<float>(4.082259521484375e+02),
		2200,
		0.0f,
		_TWOPI,
		0.95928517242269f,
		888);
	FanGeo.m_DetCntIdx = 444.75f;
	Image Img(512, 512, 500.0f, 500.0f, 0.0f, 0.0f);

	std::ifstream fprj("demo13_proj.prj", std::ios::binary);
	std::ofstream frec("demo13_recon.raw", std::ios::binary);
	if (!(fprj.is_open() && frec.is_open()))
	{
		std::cerr << "Cannot open the file for demo 13\n";
		exit(-1);
	}
	cuint imgReso = Img.m_Reso.x * Img.m_Reso.y;
	cuint prjReso = FanGeo.m_DetN * FanGeo.m_ViwN;
	cuint iterNum = 40;

	float* hprj = new float[prjReso];
	float* himg = new float[imgReso];

	fprj.read((char*) hprj, sizeof(float) * prjReso);

	EM(hprj, himg, FanGeo, Img, iterNum);

	frec.write((char*) himg, sizeof(float) * imgReso);
	fprj.close();
	frec.close();
}



template<typename T>
void CG_recon_AIM_temp(thrust::host_vector<T>& img, thrust::host_vector<T>& prj, const FanEAGeo& FanGeo, const Image& Img, cuint maxIter, T err, cuint cases)
{
	cuint len = FanGeo.m_DetN * FanGeo.m_ViwN;
	cuint wid = Img.m_Reso.x * Img.m_Reso.y;
	thrust::device_vector<T> realp(prj);
	thrust::device_vector<T> X(wid, 0);
	thrust::device_vector<T> R(wid, 0);
	thrust::device_vector<T> D(wid, 0);
	thrust::device_vector<T> tem(len, 0);
	thrust::device_vector<T> nR(wid, 0);
	thrust::device_vector<T> tempBack(wid, 0);

	thrust::device_vector<T> tempDiff(wid, 0);

	T* prealp = thrust::raw_pointer_cast(&(realp[0]));
	T* pD = thrust::raw_pointer_cast(&(D[0]));
	T* ptem = thrust::raw_pointer_cast(&(tem[0]));
	T* pR = thrust::raw_pointer_cast(&(R[0]));
	T* ptempBack = thrust::raw_pointer_cast(&(tempBack[0]));
	dim3 prjBlk(128, 8);
	dim3 prjGid(
		(FanGeo.m_DetN + prjBlk.x - 1) / prjBlk.x,
		(FanGeo.m_ViwN + prjBlk.y - 1) / prjBlk.y);
	dim3 bakBlk(32, 32);
	dim3 bakGid(
		(Img.m_Reso.x + bakBlk.x - 1) / bakBlk.x,
		(Img.m_Reso.y + bakBlk.y - 1) / bakBlk.y);

	bakproj_AIM(pR, prealp, FanGeo, Img, bakBlk, bakGid);
	D = R;
	unsigned int iter = 0;
	T alpha, beta;
	T temp1, temp2;
	while (norm<T>(R, 2.0) > err && iter < maxIter)
	{
		proj_AIM(ptem, pD, FanGeo, Img, prjBlk, prjGid);
		temp1 = R&R;  // inner product
		temp2 = tem & tem;
		alpha = temp1 / temp2;
		//X = X + (alpha * D);
		thrust::transform(X.begin(), X.end(), D.begin(), X.begin(), _saxpy_functor<T>(alpha));

		bakproj_AIM(ptempBack, ptem, FanGeo, Img, bakBlk, bakGid);
		//nR = R - (alpha *tempBack);
		thrust::transform(R.begin(), R.end(), tempBack.begin(), nR.begin(), _saxpy_functor<T>(-alpha));

		switch (cases)
		{
		case 0:
			beta = calculateBetaForCG<T, 0>(R, nR, D);
			break;
		case 1:
			beta = calculateBetaForCG<T, 1>(R, nR, D);
			break;
		case 2:
			beta = calculateBetaForCG<T, 2>(R, nR, D);
			break;
		case 3:
			beta = calculateBetaForCG<T, 3>(R, nR, D);
			break;
		case 4:
			beta = calculateBetaForCG<T, 4>(R, nR, D);
			break;
		default:
			beta = calculateBetaForCG<T, 0>(R, nR, D);
			break;
		}


		R = nR;
		//D = R + (beta * D);
		thrust::transform(R.begin(), R.end(), D.begin(), D.begin(), _saxpy_functor<T>(beta));
		++iter;
		std::cout << "Iters " << iter << std::endl;
	}
	img = X;
}

void CG_recon_AIM(thrust::host_vector<double>& img, thrust::host_vector<double>& prj, const FanEAGeo& FanGeo, const Image& Img, cuint maxIter, double err, cuint cases)
{
	CG_recon_AIM_temp<double>(img, prj, FanGeo, Img, maxIter, err, cases);
}


void CG_recon_AIM(thrust::host_vector<float>& img, thrust::host_vector<float>& prj, const FanEAGeo& FanGeo, const Image& Img, cuint maxIter, float err, cuint cases)
{
	CG_recon_AIM_temp<float>(img, prj, FanGeo, Img, maxIter, err, cases);
}



void DEMO14()
{
	typedef double DataType;

	FanEAGeo FanGeo(
		5.385200195312500e+02,
		4.082259521484375e+02,
		360, 0.0f, 6.2831853071795864769252, 0.95928517242269f / 2,
		111);
	FanGeo.m_DetCntIdx = FanGeo.m_DetN * 0.5;
	Image Img(128, 128, 300, 300, 0.0f, 0.0f);
	std::string FileName = "demo14";
	if (typeid(DataType) == typeid(double))
	{
		GenerateProjectionMatrixd_AIM(FanGeo, Img, FileName);
	}
	else
	{
		GenerateProjectionMatrixf_AIM(FanGeo, Img, FileName);
	}

}

void DEMO14f(const std::string& FileName, const FanEAGeo& FanGeo, const Image& Img)
{
	GenerateProjectionMatrixf_AIM(FanGeo, Img, FileName);
}
void DEMO14d(const std::string& FileName, const FanEAGeo& FanGeo, const Image& Img)
{
	GenerateProjectionMatrixd_AIM(FanGeo, Img, FileName);
}



void DEMO15()
{

	// Image Resolution is increasing from 128 to 2048, each time the projection process is executed 10 times
#define RANGE 20

	float tims[RANGE]; // Storing the projection time
	FanEAGeo FanGeo;
	unsigned int prjReso = FanGeo.m_DetN* FanGeo.m_ViwN;
	float* hprj = new float[prjReso];
	dim3 prjBlk(128, 8);
	dim3 prjGid(
		(FanGeo.m_DetN + prjBlk.x - 1) / prjBlk.x,
		(FanGeo.m_ViwN + prjBlk.y - 1) / prjBlk.y);

	for (int imgLReso = 128; imgLReso < 128 + RANGE; imgLReso++)
	{
		Image Img(imgLReso, imgLReso, 4.484740011196460e+02, 4.484740011196460e+02, 0.0f, 0.0f);
		unsigned int imgReso = imgLReso * imgLReso;

		//ÉêÇëÖ÷»úÄÚŽæ;
		float* himg = new float[imgReso];

		for (unsigned int i = 0; i < prjReso; i++)
		{
			hprj[i] = 0.0;
		}
		for (unsigned int i = 0; i < imgReso; i++)
		{
			himg[i] = 1.0f;
		}


		float* dimg;
		float* dprj;
		checkCudaErrors(cudaMalloc((void**) &dimg, imgReso));
		checkCudaErrors(cudaMalloc((void**) &dprj, prjReso));
		checkCudaErrors(cudaMemcpy(dimg, himg, sizeof(float) * imgReso, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(dprj, hprj, sizeof(float) * prjReso, cudaMemcpyHostToDevice));

		cudaEvent_t start, stop;
		checkCudaErrors(cudaEventCreate(&start));
		checkCudaErrors(cudaEventCreate(&stop));
		checkCudaErrors(cudaEventRecord(start, 0));
		proj(dimg, dprj, FanGeo, Img, prjBlk, prjGid);
		checkCudaErrors(cudaEventRecord(stop, 0));
		checkCudaErrors(cudaEventSynchronize(stop));
		float elapsedTime;
		checkCudaErrors(cudaEventElapsedTime(&elapsedTime, start, stop));
		tims[imgLReso - 128] = elapsedTime;

		checkCudaErrors(cudaMemcpy(hprj, dprj, sizeof(float) * prjReso, cudaMemcpyDeviceToHost));
		std::stringstream ss;
		ss << imgReso;
		std::string FileName = "From" + ss.str() + "image.raw";
		std::ofstream fid(FileName.c_str(), std::ios::binary);
		fid.write((char*) hprj, sizeof(float) * prjReso);
		fid.close();

		delete []himg;

		checkCudaErrors(cudaFree(dimg));
		checkCudaErrors(cudaFree(dprj));



	}
	std::ofstream fou("SiddonProjTime.raw", std::ios::binary);
	fou.write((char*) tims, sizeof(float) *RANGE);
	fou.close();
	delete []hprj;
}

template<typename T>
__global__ void updateImg(T* dimg, T* dcorImg, T* dwegImg, T* dmask, T lambda, cuint ImgLR, cuint ImgWR, cuint sliceNum)
{
	cuint l = threadIdx.x + blockIdx.x * blockDim.x;
	cuint w = threadIdx.y + blockIdx.y * blockDim.y;
	cuint s = threadIdx.z + blockIdx.z * blockDim.z;
	if (l < ImgLR && w < ImgWR && s < sliceNum)
	{
		cuint i = (s * ImgWR + w) * ImgLR + l;
		dimg[i] += lambda * dcorImg[i] / dwegImg[w * ImgLR + l] * dmask[w * ImgLR + l];
		if (dimg[i] < 0)
		{
			dimg[i] = 0;
		}
	}
}

template<typename T>
struct FISTA_functor
{
public:
	T t1;
	T t2;
	FISTA_functor(const T& _t1, const T& _t2) :t1(_t1), t2(_t2){}
	__device__ T operator()(const T& x1, const T& x0)
	{
		return x1 + (t1 - 1) / t2 * (x1 - x0);
	}
};


void DEMO16()
{
	std::cout << "Begin DEMO16\n";

	typedef double DataType;
	FanEAGeo FanGeo;
	Image Img;
	const unsigned int sliceNum = 22;

	FanGeo.m_DetArc = 0.95928517242269;
	FanGeo.m_DetN = 888;
	FanGeo.m_DetCntIdx = 444.75;
	FanGeo.m_DetStp = FanGeo.m_DetArc / FanGeo.m_DetN;
	FanGeo.m_O2D = 4.082259521484375e+02;

	FanGeo.m_S2O = 5.385200195312500e+02;
	FanGeo.m_S2D = FanGeo.m_S2O + FanGeo.m_O2D;
	FanGeo.m_ViwBeg = 0.0f;// (-2.668082275390625e+02 / 180.0* 3.14159265358979323846264);
	FanGeo.m_ViwN = 2200;
	FanGeo.m_ViwStp = 3.14159265358979323846264*2.0 / FanGeo.m_ViwN;


	Img.m_Bias.x = 0.0f;
	Img.m_Bias.y = 0.0f;  //ÕâžöÆ«ÒÆµÄµ¥Î»ÊÇÕæÊµÎïÀíµ¥Î»;
	Img.m_Reso.x = Img.m_Reso.y = 512;
	Img.m_Size.x = Img.m_Size.y = 500;
	Img.m_Step.x = Img.m_Size.x / Img.m_Reso.x;
	Img.m_Step.y = Img.m_Size.y / Img.m_Reso.y;

	const unsigned int imgReso = Img.m_Reso.x * Img.m_Reso.y;
	const unsigned int prjReso = FanGeo.m_DetN * FanGeo.m_ViwN;
	//DataType* img[3];
	DataType* weg[3];
	DataType* prj[3];
	DataType* mask[3];
	/*img[0] = new DataType[imgReso * sliceNum];
	img[1] = new DataType[imgReso * sliceNum];
	img[2] = new DataType[imgReso * sliceNum];*/
	weg[0] = new DataType[imgReso];
	weg[1] = new DataType[imgReso];
	weg[2] = new DataType[imgReso];
	prj[0] = new DataType[prjReso * sliceNum];
	prj[1] = new DataType[prjReso * sliceNum];
	prj[2] = new DataType[prjReso * sliceNum];
	mask[0] = new DataType[imgReso];
	mask[1] = new DataType[imgReso];
	mask[2] = new DataType[imgReso];
	DataType*totVol = new DataType[imgReso * 66];
	//std::vector<DataType> totVol(imgReso * 66, 0);

	std::string FileName1 = "demo16_partprj1d.prj";
	std::string FileName2 = "demo16_partprj2d.prj";
	std::string FileName3 = "demo16_partprj3d.prj";

	std::string FouName = "demo16_recond.raw";
	std::string WegName = "demo16_bakwegEAd.raw";
	//Mask
	std::string maskImg = "demo16_maskd.img";
	std::cout << "Begin DEMO16  2\n";
	//Read the data;
	std::ifstream fid1(FileName1.c_str(), std::ios::binary);
	std::ifstream fid2(FileName2.c_str(), std::ios::binary);
	std::ifstream fid3(FileName3.c_str(), std::ios::binary);

	std::ifstream fweg(WegName.c_str(), std::ios::binary);
	std::ifstream fmsk(maskImg.c_str(), std::ios::binary);

	fid1.read((char*) prj[0], sizeof(DataType) * prjReso * sliceNum);
	fid2.read((char*) prj[1], sizeof(DataType) * prjReso * sliceNum);
	fid3.read((char*) prj[2], sizeof(DataType) * prjReso * sliceNum);

	fweg.read((char*) weg[0], sizeof(DataType) * imgReso);
	std::memcpy(weg[1], weg[0], sizeof(DataType) * imgReso);
	std::memcpy(weg[2], weg[0], sizeof(DataType) * imgReso);

	fmsk.read((char*) mask[0], sizeof(DataType) * imgReso);
	std::memcpy(mask[1], mask[0], sizeof(DataType) * imgReso);
	std::memcpy(mask[2], mask[0], sizeof(DataType) * imgReso);
	std::cout << "Begin DEMO16 3 \n";
	fid1.close();	fid2.close();	fid3.close();	fweg.close();	fmsk.close();

	DataType* dimg[3];
	DataType* dprj[3];
	DataType* dimgWeg[3];
	DataType* dcorImg[3];
	DataType* dcorPrj[3];
	DataType* dmask[3];
	DataType* dlas[3];


	//Begin OS-SART-algorithm


	const unsigned int subSetNum = 100;
	const unsigned int numPerSubSet = FanGeo.m_ViwN / subSetNum;
	std::cout << "Begin DEMO16 4 \n";

	//Set the block and grid size
	dim3 prjBlk(128, 8);
	dim3 prjGid(
		(FanGeo.m_DetN + 127) / 128, //X dimension for detIdx
		(numPerSubSet + 7) / 8); // Y dimension for viwIdx

	dim3 bakBlk(32, 32);
	dim3 bakGid(
		(Img.m_Reso.x + 31) / 32, //X dimension for imgXidx,
		(Img.m_Reso.y + 31) / 32); // Y dimension for imgYidx
	dim3 updBlk(16, 16, 4);
	dim3 updGid(
		(Img.m_Reso.x + 15) / 16,
		(Img.m_Reso.y + 15) / 16,
		(sliceNum + 3) / 4);
	int i = 0;
	int iters = 0;
	int initNum = 1;
	int realIter = initNum *subSetNum;

	std::cout << "Begin to allocating the memory\n";
	cudaStream_t streams[3];
	for (i = 0; i < 3; i++)
	{
		checkCudaErrors(cudaSetDevice(i));
		checkCudaErrors(cudaStreamCreate(&streams[i]));

		checkCudaErrors(cudaMalloc((void**) &dimg[i], sizeof(DataType) * imgReso * sliceNum));
		checkCudaErrors(cudaMalloc((void**) &dprj[i], sizeof(DataType) * prjReso * sliceNum));
		checkCudaErrors(cudaMalloc((void**) &dimgWeg[i], sizeof(DataType) * imgReso));
		checkCudaErrors(cudaMalloc((void**) &dcorImg[i], sizeof(DataType) * imgReso * sliceNum));
		checkCudaErrors(cudaMalloc((void**) &dcorPrj[i], sizeof(DataType) * prjReso * sliceNum));
		checkCudaErrors(cudaMalloc((void**) &dmask[i], sizeof(DataType) * imgReso));
		checkCudaErrors(cudaMalloc((void**) &dlas[i], sizeof(DataType) * imgReso * sliceNum));

		checkCudaErrors(cudaMemsetAsync(dimg[i], 0, sizeof(DataType) * imgReso * sliceNum, streams[i]));
		checkCudaErrors(cudaMemcpyAsync(dprj[i], prj[i], sizeof(DataType) * prjReso * sliceNum, cudaMemcpyHostToDevice, streams[i]));
		checkCudaErrors(cudaMemcpyAsync(dimgWeg[i], weg[i], sizeof(DataType) * imgReso, cudaMemcpyHostToDevice, streams[i]));
		checkCudaErrors(cudaMemsetAsync(dcorImg[i], 0, sizeof(DataType) * imgReso * sliceNum, streams[i]));
		checkCudaErrors(cudaMemsetAsync(dcorPrj[i], 0, sizeof(DataType) * FanGeo.m_DetN * numPerSubSet  * sliceNum, streams[i]));
		checkCudaErrors(cudaMemcpyAsync(dmask[i], mask[i], sizeof(DataType) * imgReso, cudaMemcpyHostToDevice, streams[i]));
		checkCudaErrors(cudaMemsetAsync(dlas[i], 0, sizeof(DataType) * imgReso * sliceNum, streams[i]));
	}
	std::cout << "Begin DEMO16 5 \n";
	thrust::device_ptr<DataType> pimg0(dimg[0]);
	thrust::device_ptr<DataType> pimg1(dimg[1]);
	thrust::device_ptr<DataType> pimg2(dimg[2]);
	thrust::device_ptr<DataType> plas0(dlas[0]);
	thrust::device_ptr<DataType> plas1(dlas[1]);
	thrust::device_ptr<DataType> plas2(dlas[2]);
	//thrust::device_ptr<DataType> pcorImg0(dcorImg[0]);
	//thrust::device_ptr<DataType> pcorImg1(dcorImg[1]);
	//thrust::device_ptr<DataType> pcorImg2(dcorImg[2]);

	std::cout << "Begin DEMO16 6 \n";
	DataType lambda = 1.5;
	DataType t1 = 1;
	DataType t2 = 1;
	int curSubSetIdx = 0;
	std::cout << "Begin iteration\n";
	while (iters != realIter)
	{
		t2 = (1.0 + sqrt(1.0 + 4.0 * t1 * t1)) * 0.5;
		checkCudaErrors(cudaSetDevice(0));
		checkCudaErrors(cudaMemset(dcorImg[0], 0, sizeof(DataType) *imgReso * sliceNum));
		checkCudaErrors(cudaMemcpy(dlas[0], dimg[0], sizeof(DataType) * imgReso * sliceNum, cudaMemcpyDeviceToDevice));
		//plas0 = pimg0;
		proj_AIM(dimg[0], dcorPrj[0], dprj[0], FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx, prjBlk, prjGid, sliceNum, streams[0]);
		bakproj_AIM(dcorPrj[0], dimg[0], FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx, bakBlk, bakGid, sliceNum, streams[0]);
		updateImg<DataType> << <updGid, updBlk, 0, streams[0] >> >(dimg[0], dcorImg[0], dimgWeg[0], dmask[0], lambda, Img.m_Reso.x, Img.m_Reso.y, sliceNum);
		thrust::transform(pimg0, pimg0 + imgReso * sliceNum, plas0, pimg0, FISTA_functor<DataType>(t1, t2));

		checkCudaErrors(cudaSetDevice(1));
		checkCudaErrors(cudaMemset(dcorImg[1], 0, sizeof(DataType) *imgReso * sliceNum));
		//thrust::fill(pcorImg1, pcorImg1 + imgReso * sliceNum, 0.0);
		checkCudaErrors(cudaMemcpy(dlas[1], dimg[1], sizeof(DataType) * imgReso * sliceNum, cudaMemcpyDeviceToDevice));
		plas1 = pimg1;
		proj_AIM(dimg[1], dcorPrj[1], dprj[1], FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx, prjBlk, prjGid, sliceNum, streams[1]);
		bakproj_AIM(dcorPrj[1], dimg[1], FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx, bakBlk, bakGid, sliceNum, streams[1]);
		updateImg<DataType> << <updGid, updBlk, 0, streams[1] >> >(dimg[1], dcorImg[1], dimgWeg[1], dmask[1], lambda, Img.m_Reso.x, Img.m_Reso.y, sliceNum);
		thrust::transform(pimg1, pimg1 + imgReso * sliceNum, plas1, pimg1, FISTA_functor<DataType>(t1, t2));


		checkCudaErrors(cudaSetDevice(2));
		checkCudaErrors(cudaMemset(dcorImg[2], 0, sizeof(DataType) *imgReso * sliceNum));
		checkCudaErrors(cudaMemcpy(dlas[2], dimg[2], sizeof(DataType) * imgReso * sliceNum, cudaMemcpyDeviceToDevice));
		//plas2 = pimg2;
		proj_AIM(dimg[2], dcorPrj[2], dprj[2], FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx, prjBlk, prjGid, sliceNum, streams[2]);
		bakproj_AIM(dcorPrj[2], dimg[2], FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx, bakBlk, bakGid, sliceNum, streams[2]);
		updateImg<DataType> << <updGid, updBlk, 0, streams[2] >> >(dimg[2], dcorImg[2], dimgWeg[2], dmask[2], lambda, Img.m_Reso.x, Img.m_Reso.y, sliceNum);
		thrust::transform(pimg2, pimg2 + imgReso * sliceNum, plas2, pimg2, FISTA_functor<DataType>(t1, t2));

		checkCudaErrors(cudaDeviceSynchronize());
		t1 = t2;
		std::cout << "Iters: " << iters << std::endl;
		curSubSetIdx = (curSubSetIdx + 1) % subSetNum; // 
		++iters;
	}

	checkCudaErrors(cudaMemcpyAsync(totVol, dimg[0],
		sizeof(DataType) * Img.m_Reso.x * Img.m_Reso.y * sliceNum, cudaMemcpyDeviceToHost, streams[0]));
	checkCudaErrors(cudaMemcpyAsync(totVol + (Img.m_Reso.x * Img.m_Reso.y * sliceNum), dimg[1],
		sizeof(DataType) * Img.m_Reso.x * Img.m_Reso.y * sliceNum, cudaMemcpyDeviceToHost, streams[1]));
	checkCudaErrors(cudaMemcpyAsync(totVol + (Img.m_Reso.x * Img.m_Reso.y * sliceNum) + (Img.m_Reso.x * Img.m_Reso.y * sliceNum), dimg[2],
		sizeof(DataType) * Img.m_Reso.x * Img.m_Reso.y * sliceNum, cudaMemcpyDeviceToHost, streams[2]));

	std::cout << "Begin DEMO16 7 \n";
	std::ofstream fou(FouName.c_str(), std::ios::binary);
	fou.write((char*) (&totVol[0]), sizeof(DataType) * Img.m_Reso.x * Img.m_Reso.y * 66);
	fou.close();

}




void DEMO17()
{
	checkCudaErrors(cudaSetDevice(0));
	/// Statistical Iterative Reconstruction OS-SART ÖØ¹¹;
	ConeEDGeo ConeGeo(167.80, 93.91, 720, 0.0f, 6.2744586609196148290406615555556, make_float2(50.6, 50.6), make_int2(1012, 1012));
	ConeGeo.m_DetCntIdx.x = ConeGeo.m_DetN.x * 0.5 + 63.4;

	Volume Vol(512, 512, 512, 32.4, 32.4, 32.4, 0.0, 0.0, 0.0);

	cuint volReso = Vol.m_Reso.x * Vol.m_Reso.y * Vol.m_Reso.z;
	cuint prjReso = ConeGeo.m_DetN.x* ConeGeo.m_DetN.y * ConeGeo.m_ViwN;

	thrust::host_vector<float> hvol(volReso, 0);
	thrust::host_vector<float> hpho(prjReso, 1.0);
	float* hprj = new float[prjReso];// thrust::host_vector<float> hprj(prjReso, 0);
	//thrust::host_vector<float> hmsk(volReso, 1.0f); //È«1 mask

	std::ifstream fid("DEMO17_projectionData.prj", std::ios::binary);
	if (!fid.is_open())
	{
		std::cerr << "Cannot open the projection file for reconstruction\n";
		exit(-1);
	}
	fid.read((char*) hprj, sizeof(float) * prjReso);

	checkCudaErrors(cudaSetDevice(2));
	float* dprj;
	checkCudaErrors(cudaMalloc((void**) &dprj, sizeof(float) * prjReso));
	checkCudaErrors(cudaMemcpy(dprj, hprj, sizeof(float)*prjReso, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaSetDevice(0));




	//Ò»ŽÎÍ¶Ó°Ò»žöœÇ¶È
	dim3 prjBlk(32, 32);
	dim3 prjGid(
		(ConeGeo.m_DetN.x + prjBlk.x - 1) / prjBlk.x,
		(ConeGeo.m_DetN.y + prjBlk.y - 1) / prjBlk.y);

	dim3 volBlk(8, 8, 8);
	dim3 volGid(
		(Vol.m_Reso.x + volBlk.x - 1) / volBlk.x,
		(Vol.m_Reso.y + volBlk.y - 1) / volBlk.y,
		(Vol.m_Reso.z + volBlk.z - 1) / volBlk.z);

	thrust::host_vector<float>::iterator maxYIdx = thrust::max_element(hpho.begin(), hpho.end());
	const float maxY = (*maxYIdx);
	thrust::device_vector<float> dbakLambda(volReso, 0.0);
	thrust::device_vector<float> dalloneproj(ConeGeo.m_DetN.x* ConeGeo.m_DetN.y, 1.0);
	thrust::device_vector<float> drawproj(ConeGeo.m_DetN.x*ConeGeo.m_DetN.y, 0.0);
	thrust::device_vector<float> dweightedproj(ConeGeo.m_DetN.x * ConeGeo.m_DetN.y, 0.0);
	thrust::device_vector<float> dvolume(volReso, 0);
	thrust::device_vector<float> dupdVol(volReso, 0);


	float* pbakLambda = thrust::raw_pointer_cast(&dbakLambda[0]);
	float* palloneproj = thrust::raw_pointer_cast(&dalloneproj[0]);
	float* prawproj = thrust::raw_pointer_cast(&drawproj[0]);
	float* pvolume = thrust::raw_pointer_cast(&dvolume[0]);
	float* pupdVol = thrust::raw_pointer_cast(&dupdVol[0]);
	float* pweightedproj = thrust::raw_pointer_cast(&dweightedproj[0]);
	const unsigned int projR = ConeGeo.m_DetN.x * ConeGeo.m_DetN.y;
	thrust::host_vector<float> hv = dbakLambda;
	//ŒÆËãÈ«1µÄ·ŽÍ¶Ó°;
	std::ifstream weiFile("DEMO17_weight_.raw", std::ios::binary);
	if (!weiFile.is_open())
	{
		std::cout << "Begin all one backprojction for weighting\n";
		for (unsigned int i = 0; i != ConeGeo.m_ViwN; ++i)
		{
			back_proj_DEMO17(palloneproj, pbakLambda, ConeGeo, Vol, i, volBlk, volGid);
			checkCudaErrors(cudaDeviceSynchronize());
			std::cout << i << std::endl;
		}
		hv = dbakLambda;
		for (unsigned int sliceIdx = 0; sliceIdx != 512; ++sliceIdx)
		{
			for (unsigned int jj = 0; jj != 512; ++jj)
			{
				for (unsigned int ii = 0; ii != 512; ++ii)
				{
					unsigned int i = (sliceIdx * 512 + jj) * 512 + ii;
					if (((ii - 256.0 + 0.5) / 256.0)*((ii - 256.0 + 0.5) / 256.0) + ((jj - 256.0 + 0.5) / 256.0)*((jj - 256.0 + 0.5) / 256.0) >= 1.44)
					{
						hv[i] = 0;
					}
				}
			}
		}
		std::ofstream weight("DEMO17_weight_.raw", std::ios::binary);
		weight.write((char*) (&hv[0]), sizeof(float) * volReso);
		weight.close();
	}
	else
	{
		weiFile.read((char*) (&hv[0]), sizeof(float) * volReso);
		weiFile.close();
	}

	dbakLambda = hv;
	std::cout << "Finished the backprojction\n";
	unsigned int iterIdx = 0;
	unsigned int subsetIdx = 0;
	unsigned int subangIdx = 0;
	unsigned int angIdx = 0;
	cuint iterNum = 1;
	int subSetNum = 20;
	int numPerSubSet = ConeGeo.m_ViwN / subSetNum;
	//int realIter = subSetNum *iterNum;
	for (iterIdx = 0; iterIdx != iterNum; ++iterIdx)
	{
		if (iterIdx == 0)
		{
			subSetNum = 144;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 1)
		{
			subSetNum = 120;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 2)
		{
			subSetNum = 90;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 3)
		{
			subSetNum = 80;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 4)
		{
			subSetNum = 72;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 5)
		{
			subSetNum = 60;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 6)
		{
			subSetNum = 48;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 7)
		{
			subSetNum = 24;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 8)
		{
			subSetNum = 12;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 9)
		{
			subSetNum = 6;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 10)
		{
			subSetNum = 3;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx > 10 && iterIdx < 16)
		{
			subSetNum = 3;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx >= 16)
		{
			subSetNum = 2;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}

		for (subsetIdx = 0; subsetIdx != subSetNum; ++subsetIdx)
		{
			//Set the volume as 0;
			thrust::fill(dupdVol.begin(), dupdVol.end(), 0.0f);
			//Calculate current angle;
			for (subangIdx = 0; subangIdx != numPerSubSet; ++subangIdx)
			{
				angIdx = subsetIdx + subangIdx *subSetNum;
				//cudaMemcpy(prawproj, hprj + angIdx *projR, sizeof(float) * projR, cudaMemcpyHostToDevice);
				checkCudaErrors(cudaMemcpy(prawproj, dprj + angIdx *projR, sizeof(float) * projR, cudaMemcpyDeviceToDevice));
				proj_DEMO17(pvolume, pweightedproj, prawproj, ConeGeo, Vol, angIdx, prjBlk, prjGid);
				//projection;
				//bakprojection
				back_proj_DEMO17(pweightedproj, pupdVol, ConeGeo, Vol, angIdx, volBlk, volGid);
				std::cout << "angIdx: " << angIdx << std::endl;
			}

			//Update the volume
			thrust::transform(
				thrust::make_zip_iterator(
				thrust::make_tuple(dvolume.begin(), dupdVol.begin(), dbakLambda.begin())
				),
				thrust::make_zip_iterator(
				thrust::make_tuple(dvolume.end(), dupdVol.end(), dbakLambda.end())
				),
				dvolume.begin(), _DEMO17_update_functor < float >(subSetNum));

		}
		std::cout << "Iteration: " << iterIdx << std::endl;
	}

	hvol = dvolume;
	std::ofstream fou("DEMO17_reconstruction.raw", std::ios::binary);
	fou.write((char*) (&hvol[0]), sizeof(float) * volReso);
	fou.close();
	std::cout << "Finished\n";
}



void DEMO17(
	const ConeEDGeo& ConeGeo,
	const Volume& Vol,
	const std::string& projectionFileName,
	const std::string& weightingFileName,
	const std::string& reconstFileName,
	const std::vector<int>& subsetNumberSeries,
	const int iterNum)
{
	checkCudaErrors(cudaSetDevice(0));
	/// Statistical Iterative Reconstruction OS-SART ÖØ¹¹;

	cuint volReso = Vol.m_Reso.x * Vol.m_Reso.y * Vol.m_Reso.z;
	cuint prjReso = ConeGeo.m_DetN.x* ConeGeo.m_DetN.y * ConeGeo.m_ViwN;

	thrust::host_vector<float> hvol(volReso, 0);
	thrust::host_vector<float> hpho(prjReso, 1.0);
	float* hprj = new float[prjReso];// thrust::host_vector<float> hprj(prjReso, 0);
	//thrust::host_vector<float> hmsk(volReso, 1.0f); //È«1 mask

	std::ifstream fid(projectionFileName.c_str(), std::ios::binary);
	if (!fid.is_open())
	{
		std::cerr << "Cannot open the projection file for reconstruction\n";
		exit(-1);
	}
	fid.read((char*) hprj, sizeof(float) * prjReso);

	checkCudaErrors(cudaSetDevice(2));
	float* dprj;
	checkCudaErrors(cudaMalloc((void**) &dprj, sizeof(float) * prjReso));
	checkCudaErrors(cudaMemcpy(dprj, hprj, sizeof(float)*prjReso, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaSetDevice(0));




	//Ò»ŽÎÍ¶Ó°Ò»žöœÇ¶È
	dim3 prjBlk(32, 32);
	dim3 prjGid(
		(ConeGeo.m_DetN.x + prjBlk.x - 1) / prjBlk.x,
		(ConeGeo.m_DetN.y + prjBlk.y - 1) / prjBlk.y);

	dim3 volBlk(8, 8, 8);
	dim3 volGid(
		(Vol.m_Reso.x + volBlk.x - 1) / volBlk.x,
		(Vol.m_Reso.y + volBlk.y - 1) / volBlk.y,
		(Vol.m_Reso.z + volBlk.z - 1) / volBlk.z);

	thrust::host_vector<float>::iterator maxYIdx = thrust::max_element(hpho.begin(), hpho.end());
	const float maxY = (*maxYIdx);
	thrust::device_vector<float> dbakLambda(volReso, 0.0);
	thrust::device_vector<float> dalloneproj(ConeGeo.m_DetN.x* ConeGeo.m_DetN.y, 1.0);
	thrust::device_vector<float> drawproj(ConeGeo.m_DetN.x*ConeGeo.m_DetN.y, 0.0);
	thrust::device_vector<float> dweightedproj(ConeGeo.m_DetN.x * ConeGeo.m_DetN.y, 0.0);
	thrust::device_vector<float> dvolume(volReso, 0);
	thrust::device_vector<float> dupdVol(volReso, 0);


	float* pbakLambda = thrust::raw_pointer_cast(&dbakLambda[0]);
	float* palloneproj = thrust::raw_pointer_cast(&dalloneproj[0]);
	float* prawproj = thrust::raw_pointer_cast(&drawproj[0]);
	float* pvolume = thrust::raw_pointer_cast(&dvolume[0]);
	float* pupdVol = thrust::raw_pointer_cast(&dupdVol[0]);
	float* pweightedproj = thrust::raw_pointer_cast(&dweightedproj[0]);
	const unsigned int projR = ConeGeo.m_DetN.x * ConeGeo.m_DetN.y;
	thrust::host_vector<float> hv = dbakLambda;
	//ŒÆËãÈ«1µÄ·ŽÍ¶Ó°;
	std::ifstream weiFile(weightingFileName.c_str(), std::ios::binary);
	const float hfL = Vol.m_Reso.x * 0.5;
	const float hfW = Vol.m_Reso.y * 0.5;

	if (!weiFile.is_open())
	{
		std::cout << "Begin all one backprojction for weighting\n";
		for (unsigned int i = 0; i != ConeGeo.m_ViwN; ++i)
		{
			back_proj_DEMO17(palloneproj, pbakLambda, ConeGeo, Vol, i, volBlk, volGid);
			checkCudaErrors(cudaDeviceSynchronize());
			std::cout << i << std::endl;
		}
		hv = dbakLambda;
		//for (unsigned int sliceIdx = 0; sliceIdx != 512; ++sliceIdx)
		//{
		//	for (unsigned int jj = 0; jj != 512; ++jj)
		//	{
		//		for (unsigned int ii = 0; ii != 512; ++ii)
		//		{
		//			unsigned int i = (sliceIdx * 512 + jj) * 512 + ii;
		//			if (((ii - 256.0 + 0.5) / 256.0)*((ii - 256.0 + 0.5) / 256.0) + ((jj - 256.0 + 0.5) / 256.0)*((jj - 256.0 + 0.5) / 256.0) >= 1.44)
		//			{
		//				hv[i] = 0;
		//			}
		//		}
		//	}
		//}
		for (unsigned int sliceIdx = 0; sliceIdx != Vol.m_Reso.z; ++sliceIdx)
		{
			for (unsigned int jj = 0; jj != Vol.m_Reso.y; ++jj)
			{
				for (unsigned int ii = 0; ii != Vol.m_Reso.x; ++ii)
				{
					unsigned int i = (sliceIdx * Vol.m_Reso.y + jj) * Vol.m_Reso.x + ii;
					if (((ii - hfL + 0.5) / hfL)*((ii - hfL + 0.5) / hfL) + ((jj - hfW + 0.5) / hfL)*((jj - hfW + 0.5) / hfW) >= 1.44)
					{
						hv[i] = 0;
					}
				}
			}
		}
		std::ofstream weight(weightingFileName.c_str(), std::ios::binary);
		weight.write((char*) (&hv[0]), sizeof(float) * volReso);
		weight.close();
	}
	else
	{
		weiFile.read((char*) (&hv[0]), sizeof(float) * volReso);
		weiFile.close();
	}

	dbakLambda = hv;
	std::cout << "Finished the backprojction\n";
	unsigned int iterIdx = 0;
	unsigned int subsetIdx = 0;
	unsigned int subangIdx = 0;
	unsigned int angIdx = 0;
	//cuint iterNum = 1;
	int subSetNum = 20;
	int numPerSubSet = ConeGeo.m_ViwN / subSetNum;
	const int subsetNumberSize = subsetNumberSeries.size();
	//int realIter = subSetNum *iterNum;
	for (iterIdx = 0; iterIdx != iterNum; ++iterIdx)
	{
		if (iterIdx < subsetNumberSize)
		{
			subSetNum = subsetNumberSeries[iterIdx];
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		else
		{
			subSetNum = 1;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}

		for (subsetIdx = 0; subsetIdx != subSetNum; ++subsetIdx)
		{
			//Set the volume as 0;
			thrust::fill(dupdVol.begin(), dupdVol.end(), 0.0f);
			//Calculate current angle;
			for (subangIdx = 0; subangIdx != numPerSubSet; ++subangIdx)
			{
				angIdx = subsetIdx + subangIdx *subSetNum;
				//cudaMemcpy(prawproj, hprj + angIdx *projR, sizeof(float) * projR, cudaMemcpyHostToDevice);
				checkCudaErrors(cudaMemcpy(prawproj, dprj + angIdx *projR, sizeof(float) * projR, cudaMemcpyDeviceToDevice));
				proj_DEMO17(pvolume, pweightedproj, prawproj, ConeGeo, Vol, angIdx, prjBlk, prjGid);
				//projection;
				//bakprojection
				back_proj_DEMO17(pweightedproj, pupdVol, ConeGeo, Vol, angIdx, volBlk, volGid);
				std::cout << "angIdx: " << angIdx << std::endl;
			}

			//Update the volume
			thrust::transform(
				thrust::make_zip_iterator(
				thrust::make_tuple(dvolume.begin(), dupdVol.begin(), dbakLambda.begin())
				),
				thrust::make_zip_iterator(
				thrust::make_tuple(dvolume.end(), dupdVol.end(), dbakLambda.end())
				),
				dvolume.begin(), _DEMO17_update_functor < float >(subSetNum, 0, 0.6));

		}
		std::cout << "Iteration: " << iterIdx << std::endl;
	}

	hvol = dvolume;
	std::ofstream fou(reconstFileName.c_str(), std::ios::binary);
	fou.write((char*) (&hvol[0]), sizeof(float) * volReso);
	fou.close();
	std::cout << "Finished\n";
}



void DEMO18()
{
	checkCudaErrors(cudaSetDevice(2));
	/// Statistical Iterative Reconstruction OS-SART ÖØ¹¹;
	ConeEDGeo ConeGeo(167.80, 93.91, 400, 0.0f, 3.4906585039886591538473777777778, make_float2(50.6, 50.6), make_int2(1012, 1012));
	ConeGeo.m_DetCntIdx.x = ConeGeo.m_DetN.x * 0.5 + 65.4;

	Volume Vol(512, 512, 512, 32.4, 32.4, 32.4, 0.0, 0.0, 0.0);

	cuint volReso = Vol.m_Reso.x * Vol.m_Reso.y * Vol.m_Reso.z;
	cuint prjReso = ConeGeo.m_DetN.x* ConeGeo.m_DetN.y * ConeGeo.m_ViwN;

	thrust::host_vector<float> hvol(volReso, 0);
	thrust::host_vector<float> hpho(prjReso, 1.0);
	float* hprj = new float[prjReso];// thrust::host_vector<float> hprj(prjReso, 0);
	//thrust::host_vector<float> hmsk(volReso, 1.0f); //È«1 mask

	std::ifstream fid("DEMO18_projectionData.prj", std::ios::binary);
	if (!fid.is_open())
	{
		std::cerr << "Cannot open the projection file for reconstruction\n";
		exit(-1);
	}
	fid.read((char*) hprj, sizeof(float) * prjReso);

	//Ò»ŽÎÍ¶Ó°Ò»žöœÇ¶È;
	dim3 prjBlk(32, 32);
	dim3 prjGid(
		(ConeGeo.m_DetN.x + prjBlk.x - 1) / prjBlk.x,
		(ConeGeo.m_DetN.y + prjBlk.y - 1) / prjBlk.y);

	dim3 volBlk(8, 8, 8);
	dim3 volGid(
		(Vol.m_Reso.x + volBlk.x - 1) / volBlk.x,
		(Vol.m_Reso.y + volBlk.y - 1) / volBlk.y,
		(Vol.m_Reso.z + volBlk.z - 1) / volBlk.z);

	thrust::host_vector<float>::iterator maxYIdx = thrust::max_element(hpho.begin(), hpho.end());
	const float maxY = (*maxYIdx);
	thrust::device_vector<float> dbakLambda(volReso, 0.0);
	thrust::device_vector<float> dalloneproj(ConeGeo.m_DetN.x* ConeGeo.m_DetN.y, 1.0);
	thrust::device_vector<float> drawproj(ConeGeo.m_DetN.x*ConeGeo.m_DetN.y, 0.0);
	thrust::device_vector<float> dweightedproj(ConeGeo.m_DetN.x * ConeGeo.m_DetN.y, 0.0);
	thrust::device_vector<float> dvolume(volReso, 0);
	thrust::device_vector<float> dupdVol(volReso, 0);


	float* pbakLambda = thrust::raw_pointer_cast(&dbakLambda[0]);
	float* palloneproj = thrust::raw_pointer_cast(&dalloneproj[0]);
	float* prawproj = thrust::raw_pointer_cast(&drawproj[0]);
	float* pvolume = thrust::raw_pointer_cast(&dvolume[0]);
	float* pupdVol = thrust::raw_pointer_cast(&dupdVol[0]);
	float* pweightedproj = thrust::raw_pointer_cast(&dweightedproj[0]);
	const unsigned int projR = ConeGeo.m_DetN.x * ConeGeo.m_DetN.y;
	thrust::host_vector<float> hv = dbakLambda;
	const float hfL = Vol.m_Reso.x * 0.5;
	const float hfW = Vol.m_Reso.y * 0.5;
	//ŒÆËãÈ«1µÄ·ŽÍ¶Ó°;
	std::ifstream weiFile("DEMO18_weight.raw", std::ios::binary);
	if (!weiFile.is_open())
	{
		std::cout << "Begin all one backprojction for weighting\n";
		for (unsigned int i = 0; i != ConeGeo.m_ViwN; ++i)
		{
			back_proj_DEMO17(palloneproj, pbakLambda, ConeGeo, Vol, i, volBlk, volGid);
			checkCudaErrors(cudaDeviceSynchronize());
			std::cout << i << std::endl;
		}
		hv = dbakLambda;
		for (unsigned int sliceIdx = 0; sliceIdx != Vol.m_Reso.z; ++sliceIdx)
		{
			for (unsigned int jj = 0; jj != Vol.m_Reso.y; ++jj)
			{
				for (unsigned int ii = 0; ii != Vol.m_Reso.x; ++ii)
				{
					unsigned int i = (sliceIdx * Vol.m_Reso.y + jj) * Vol.m_Reso.x + ii;
					if (((ii - hfL + 0.5) / hfL)*((ii - hfL + 0.5) / hfL) + ((jj - hfW + 0.5) / hfL)*((jj - hfW + 0.5) / hfW) >= 1.44)
					{
						hv[i] = 0;
					}
				}
			}
		}
		std::ofstream weight("DEMO18_weight.raw", std::ios::binary);
		weight.write((char*) (&hv[0]), sizeof(float) * volReso);
		weight.close();
	}
	else
	{
		weiFile.read((char*) (&hv[0]), sizeof(float) * volReso);
		weiFile.close();
	}

	dbakLambda = hv;
	std::cout << "Finished the backprojction\n";
	unsigned int iterIdx = 0;
	unsigned int subsetIdx = 0;
	unsigned int subangIdx = 0;
	unsigned int angIdx = 0;
	cuint iterNum = 1; // We use 20 iterations.
	int subSetNum = 20;
	int numPerSubSet = ConeGeo.m_ViwN / subSetNum;

	for (iterIdx = 0; iterIdx != iterNum; ++iterIdx)
	{
		if (iterIdx == 0)
		{
			subSetNum = 100;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 1)
		{
			subSetNum = 80;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 2)
		{
			subSetNum = 50;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 3)
		{
			subSetNum = 40;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 4)
		{
			subSetNum = 40;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 5)
		{
			subSetNum = 20;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 6)
		{
			subSetNum = 20;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 7)
		{
			subSetNum = 20;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 8)
		{
			subSetNum = 10;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 9)
		{
			subSetNum = 10;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 10)
		{
			subSetNum = 10;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx > 10 && iterIdx < 16)
		{
			subSetNum = 5;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx >= 16)
		{
			subSetNum = 2;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}

		for (subsetIdx = 0; subsetIdx != subSetNum; ++subsetIdx)
		{
			//Set the volume as 0;
			thrust::fill(dupdVol.begin(), dupdVol.end(), 0.0f);
			//Calculate current angle;
			for (subangIdx = 0; subangIdx != numPerSubSet; ++subangIdx)
			{
				angIdx = subsetIdx + subangIdx *subSetNum;
				checkCudaErrors(cudaMemcpy(prawproj, hprj + angIdx *projR, sizeof(float) * projR, cudaMemcpyHostToDevice));
				//cudaMemcpy(prawproj, dprj + angIdx *projR, sizeof(float) * projR, cudaMemcpyDeviceToDevice);
				proj_DEMO17(pvolume, pweightedproj, prawproj, ConeGeo, Vol, angIdx, prjBlk, prjGid);
				//bakprojection
				back_proj_DEMO17(pweightedproj, pupdVol, ConeGeo, Vol, angIdx, volBlk, volGid);
				std::cout << "angIdx: " << angIdx << std::endl;
			}

			//Update the volume
			thrust::transform(
				thrust::make_zip_iterator(
				thrust::make_tuple(dvolume.begin(), dupdVol.begin(), dbakLambda.begin())
				),
				thrust::make_zip_iterator(
				thrust::make_tuple(dvolume.end(), dupdVol.end(), dbakLambda.end())
				),
				dvolume.begin(), _DEMO18_update_functor< float >(subSetNum, 0, 0.6));

		}
		std::cout << "Iteration: " << iterIdx << std::endl;
	}

	hvol = dvolume;
	std::ofstream fou("DEMO18_reconstruction.raw", std::ios::binary);
	fou.write((char*) (&hvol[0]), sizeof(float) * volReso);
	fou.close();
	std::cout << "Finished\n";
}



void DEMO18(
	const ConeEDGeo& ConeGeo,
	const Volume& Vol,
	const std::string& projectionFileName,
	const std::string& weightingFileName,
	const std::string& reconstrctFileName,
	const std::vector < int >& subSetNumSeries,
	const int iterNum)
{
	checkCudaErrors(cudaSetDevice(0));
	/// Statistical Iterative Reconstruction OS-SART ÖØ¹¹;

	cuint volReso = Vol.m_Reso.x * Vol.m_Reso.y * Vol.m_Reso.z;
	cuint prjReso = ConeGeo.m_DetN.x* ConeGeo.m_DetN.y * ConeGeo.m_ViwN;

	thrust::host_vector<float> hvol(volReso, 0);
	thrust::host_vector<float> hpho(prjReso, 1.0);
	float* hprj = new float[prjReso];// thrust::host_vector<float> hprj(prjReso, 0);
	//thrust::host_vector<float> hmsk(volReso, 1.0f); //È«1 mask

	std::ifstream fid(projectionFileName.c_str(), std::ios::binary);
	if (!fid.is_open())
	{
		std::cerr << "Cannot open the projection file for reconstruction\n";
		exit(-1);
	}
	fid.read((char*) hprj, sizeof(float) * prjReso);

	//Ò»ŽÎÍ¶Ó°Ò»žöœÇ¶È;
	dim3 prjBlk(32, 32);
	dim3 prjGid(
		(ConeGeo.m_DetN.x + prjBlk.x - 1) / prjBlk.x,
		(ConeGeo.m_DetN.y + prjBlk.y - 1) / prjBlk.y);

	dim3 volBlk(8, 8, 8);
	dim3 volGid(
		(Vol.m_Reso.x + volBlk.x - 1) / volBlk.x,
		(Vol.m_Reso.y + volBlk.y - 1) / volBlk.y,
		(Vol.m_Reso.z + volBlk.z - 1) / volBlk.z);

	thrust::host_vector<float>::iterator maxYIdx = thrust::max_element(hpho.begin(), hpho.end());
	const float maxY = (*maxYIdx);
	thrust::device_vector<float> dbakLambda(volReso, 0.0);
	thrust::device_vector<float> dalloneproj(ConeGeo.m_DetN.x* ConeGeo.m_DetN.y, 1.0);
	thrust::device_vector<float> drawproj(ConeGeo.m_DetN.x*ConeGeo.m_DetN.y, 0.0);
	thrust::device_vector<float> dweightedproj(ConeGeo.m_DetN.x * ConeGeo.m_DetN.y, 0.0);
	thrust::device_vector<float> dvolume(volReso, 0);
	thrust::device_vector<float> dupdVol(volReso, 0);


	float* pbakLambda = thrust::raw_pointer_cast(&dbakLambda[0]);
	float* palloneproj = thrust::raw_pointer_cast(&dalloneproj[0]);
	float* prawproj = thrust::raw_pointer_cast(&drawproj[0]);
	float* pvolume = thrust::raw_pointer_cast(&dvolume[0]);
	float* pupdVol = thrust::raw_pointer_cast(&dupdVol[0]);
	float* pweightedproj = thrust::raw_pointer_cast(&dweightedproj[0]);
	const unsigned int projR = ConeGeo.m_DetN.x * ConeGeo.m_DetN.y;
	thrust::host_vector<float> hv = dbakLambda;
	const float hfL = Vol.m_Reso.x * 0.5;
	const float hfW = Vol.m_Reso.y * 0.5;
	//ŒÆËãÈ«1µÄ·ŽÍ¶Ó°;
	std::ifstream weiFile(weightingFileName.c_str(), std::ios::binary);
	if (!weiFile.is_open())
	{
		std::cout << "Begin all one backprojction for weighting\n";
		for (unsigned int i = 0; i != ConeGeo.m_ViwN; ++i)
		{
			back_proj_DEMO17(palloneproj, pbakLambda, ConeGeo, Vol, i, volBlk, volGid);
			checkCudaErrors(cudaDeviceSynchronize());
			std::cout << i << std::endl;
		}
		hv = dbakLambda;
		for (unsigned int sliceIdx = 0; sliceIdx != Vol.m_Reso.z; ++sliceIdx)
		{
			for (unsigned int jj = 0; jj != Vol.m_Reso.y; ++jj)
			{
				for (unsigned int ii = 0; ii != Vol.m_Reso.x; ++ii)
				{
					unsigned int i = (sliceIdx * Vol.m_Reso.y + jj) * Vol.m_Reso.x + ii;
					if (((ii - hfL + 0.5) / hfL)*((ii - hfL + 0.5) / hfL) + ((jj - hfW + 0.5) / hfL)*((jj - hfW + 0.5) / hfW) >= 1.44)
					{
						hv[i] = 0;
					}
				}
			}
		}
		std::ofstream weight(reconstrctFileName.c_str(), std::ios::binary);
		weight.write((char*) (&hv[0]), sizeof(float) * volReso);
		weight.close();
	}
	else
	{
		weiFile.read((char*) (&hv[0]), sizeof(float) * volReso);
		weiFile.close();
	}

	dbakLambda = hv;
	std::cout << "Finished the backprojction\n";
	unsigned int iterIdx = 0;
	unsigned int subsetIdx = 0;
	unsigned int subangIdx = 0;
	unsigned int angIdx = 0;

	int subSetNum = 20;
	int numPerSubSet = ConeGeo.m_ViwN / subSetNum;
	int subsetNumSeriesSize = subSetNumSeries.size();

	for (iterIdx = 0; iterIdx != iterNum; ++iterIdx)
	{
		if (iterIdx <= subsetNumSeriesSize)
		{
			subSetNum = subSetNumSeries[iterIdx];
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		else
		{
			subSetNum = 1;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}

		for (subsetIdx = 0; subsetIdx != subSetNum; ++subsetIdx)
		{
			//Set the volume as 0;
			thrust::fill(dupdVol.begin(), dupdVol.end(), 0.0f);
			//Calculate current angle;
			for (subangIdx = 0; subangIdx != numPerSubSet; ++subangIdx)
			{
				angIdx = subsetIdx + subangIdx *subSetNum;
				checkCudaErrors(cudaMemcpy(prawproj, hprj + angIdx *projR, sizeof(float) * projR, cudaMemcpyHostToDevice));
				//cudaMemcpy(prawproj, dprj + angIdx *projR, sizeof(float) * projR, cudaMemcpyDeviceToDevice);
				proj_DEMO17(pvolume, pweightedproj, prawproj, ConeGeo, Vol, angIdx, prjBlk, prjGid);
				//bakprojection
				back_proj_DEMO17(pweightedproj, pupdVol, ConeGeo, Vol, angIdx, volBlk, volGid);
				std::cout << "angIdx: " << angIdx << std::endl;
			}

			//Update the volume
			thrust::transform(
				thrust::make_zip_iterator(
				thrust::make_tuple(dvolume.begin(), dupdVol.begin(), dbakLambda.begin())
				),
				thrust::make_zip_iterator(
				thrust::make_tuple(dvolume.end(), dupdVol.end(), dbakLambda.end())
				),
				dvolume.begin(), _DEMO18_update_functor< float >(subSetNum, 0, 0.6));

		}
		std::cout << "Iteration: " << iterIdx << std::endl;
	}

	hvol = dvolume;
	std::ofstream fou(reconstrctFileName.c_str(), std::ios::binary);
	fou.write((char*) (&hvol[0]), sizeof(float) * volReso);
	fou.close();
	std::cout << "Finished\n";
}


void DEMO18v2()
{

	checkCudaErrors(cudaSetDevice(0));
	/// Statistical Iterative Reconstruction OS-SART ÖØ¹¹;
	ConeEDGeo ConeGeo(167.80, 93.91, 400, 0.0f, 3.4906585039886591538473777777778, make_float2(50.6, 50.6), make_int2(2024, 2024));
	ConeGeo.m_DetCntIdx.x = ConeGeo.m_DetN.x * 0.5 + 65.4 * 2;

	Volume Vol(512, 512, 512, 32.4, 32.4, 32.4, 0.0, 0.0, 0.0);

	cuint volReso = Vol.m_Reso.x * Vol.m_Reso.y * Vol.m_Reso.z;
	cuint prjReso = ConeGeo.m_DetN.x* ConeGeo.m_DetN.y * ConeGeo.m_ViwN;

	thrust::host_vector<float> hvol(volReso, 0);
	thrust::host_vector<float> hpho(prjReso, 1.0);
	float* hprj = new float[prjReso];// thrust::host_vector<float> hprj(prjReso, 0);
	//thrust::host_vector<float> hmsk(volReso, 1.0f); //È«1 mask

	std::ifstream fid("DEMO18v2_projectionData.prj", std::ios::binary);
	if (!fid.is_open())
	{
		std::cerr << "Cannot open the projection file for reconstruction\n";
		exit(-1);
	}
	fid.read((char*) hprj, sizeof(float) * prjReso);

	//cudaSetDevice(2);
	//float* dprj;
	//cudaMalloc((void**) &dprj, sizeof(float) * prjReso);
	//cudaMemcpy(dprj, hprj, sizeof(float)*prjReso, cudaMemcpyHostToDevice);
	//cudaSetDevice(0);


	//Ò»ŽÎÍ¶Ó°Ò»žöœÇ¶È;
	dim3 prjBlk(32, 32);
	dim3 prjGid(
		(ConeGeo.m_DetN.x + prjBlk.x - 1) / prjBlk.x,
		(ConeGeo.m_DetN.y + prjBlk.y - 1) / prjBlk.y);

	dim3 volBlk(8, 8, 8);
	dim3 volGid(
		(Vol.m_Reso.x + volBlk.x - 1) / volBlk.x,
		(Vol.m_Reso.y + volBlk.y - 1) / volBlk.y,
		(Vol.m_Reso.z + volBlk.z - 1) / volBlk.z);

	thrust::host_vector<float>::iterator maxYIdx = thrust::max_element(hpho.begin(), hpho.end());
	const float maxY = (*maxYIdx);
	thrust::device_vector<float> dbakLambda(volReso, 0.0);
	thrust::device_vector<float> dalloneproj(ConeGeo.m_DetN.x* ConeGeo.m_DetN.y, 1.0);
	thrust::device_vector<float> drawproj(ConeGeo.m_DetN.x*ConeGeo.m_DetN.y, 0.0);
	thrust::device_vector<float> dweightedproj(ConeGeo.m_DetN.x * ConeGeo.m_DetN.y, 0.0);
	thrust::device_vector<float> dvolume(volReso, 0);
	thrust::device_vector<float> dupdVol(volReso, 0);


	float* pbakLambda = thrust::raw_pointer_cast(&dbakLambda[0]);
	float* palloneproj = thrust::raw_pointer_cast(&dalloneproj[0]);
	float* prawproj = thrust::raw_pointer_cast(&drawproj[0]);
	float* pvolume = thrust::raw_pointer_cast(&dvolume[0]);
	float* pupdVol = thrust::raw_pointer_cast(&dupdVol[0]);
	float* pweightedproj = thrust::raw_pointer_cast(&dweightedproj[0]);
	const unsigned int projR = ConeGeo.m_DetN.x * ConeGeo.m_DetN.y;
	thrust::host_vector<float> hv = dbakLambda;
	const float hfL = Vol.m_Reso.x * 0.5;
	const float hfW = Vol.m_Reso.y * 0.5;
	//ŒÆËãÈ«1µÄ·ŽÍ¶Ó°;
	std::ifstream weiFile("DEMO18v2_weight.raw", std::ios::binary);
	if (!weiFile.is_open())
	{
		std::cout << "Begin all one backprojction for weighting\n";
		for (unsigned int i = 0; i != ConeGeo.m_ViwN; ++i)
		{
			back_proj_DEMO17(palloneproj, pbakLambda, ConeGeo, Vol, i, volBlk, volGid);
			checkCudaErrors(cudaDeviceSynchronize());
			std::cout << i << std::endl;
		}
		hv = dbakLambda;
		for (unsigned int sliceIdx = 0; sliceIdx != 512; ++sliceIdx)
		{
			for (unsigned int jj = 0; jj != Vol.m_Reso.y; ++jj)
			{
				for (unsigned int ii = 0; ii != Vol.m_Reso.x; ++ii)
				{
					unsigned int i = (sliceIdx * Vol.m_Reso.y + jj) * Vol.m_Reso.x + ii;
					if (((ii - hfL + 0.5) / hfL)*((ii - hfL + 0.5) / hfL) + ((jj - hfW + 0.5) / hfL)*((jj - hfW + 0.5) / hfW) >= 1.44)
					{
						hv[i] = 0;
					}
				}
			}
		}
		std::ofstream weight("DEMO18v2_weight.raw", std::ios::binary);
		weight.write((char*) (&hv[0]), sizeof(float) * volReso);
		weight.close();
	}
	else
	{
		weiFile.read((char*) (&hv[0]), sizeof(float) * volReso);
		weiFile.close();
	}

	dbakLambda = hv;
	std::cout << "Finished the backprojction\n";
	unsigned int iterIdx = 0;
	unsigned int subsetIdx = 0;
	unsigned int subangIdx = 0;
	unsigned int angIdx = 0;
	cuint iterNum = 1; // We use 20 iterations.
	int subSetNum = 20;
	int numPerSubSet = ConeGeo.m_ViwN / subSetNum;

	for (iterIdx = 0; iterIdx != iterNum; ++iterIdx)
	{
		if (iterIdx == 0)
		{
			subSetNum = 100;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 1)
		{
			subSetNum = 80;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 2)
		{
			subSetNum = 50;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 3)
		{
			subSetNum = 40;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 4)
		{
			subSetNum = 40;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 5)
		{
			subSetNum = 20;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 6)
		{
			subSetNum = 20;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 7)
		{
			subSetNum = 20;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 8)
		{
			subSetNum = 10;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 9)
		{
			subSetNum = 10;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx == 10)
		{
			subSetNum = 10;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx > 10 && iterIdx < 16)
		{
			subSetNum = 5;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}
		if (iterIdx >= 16)
		{
			subSetNum = 2;
			numPerSubSet = ConeGeo.m_ViwN / subSetNum;
		}

		for (subsetIdx = 0; subsetIdx != subSetNum; ++subsetIdx)
		{
			//Set the volume as 0;
			thrust::fill(dupdVol.begin(), dupdVol.end(), 0.0f);
			//Calculate current angle;
			for (subangIdx = 0; subangIdx != numPerSubSet; ++subangIdx)
			{
				angIdx = subsetIdx + subangIdx *subSetNum;
				checkCudaErrors(cudaMemcpy(prawproj, hprj + angIdx *projR, sizeof(float) * projR, cudaMemcpyHostToDevice));
				//cudaMemcpy(prawproj, dprj + angIdx *projR, sizeof(float) * projR, cudaMemcpyDeviceToDevice);
				proj_DEMO17(pvolume, pweightedproj, prawproj, ConeGeo, Vol, angIdx, prjBlk, prjGid);
				//bakprojection
				back_proj_DEMO17(pweightedproj, pupdVol, ConeGeo, Vol, angIdx, volBlk, volGid);
				std::cout << "angIdx: " << angIdx << std::endl;
			}

			//Update the volume
			thrust::transform(
				thrust::make_zip_iterator(
				thrust::make_tuple(dvolume.begin(), dupdVol.begin(), dbakLambda.begin())
				),
				thrust::make_zip_iterator(
				thrust::make_tuple(dvolume.end(), dupdVol.end(), dbakLambda.end())
				),
				dvolume.begin(), _DEMO18_update_functor< float >(subSetNum, 0, 1.2));

		}
		std::cout << "Iteration: " << iterIdx << std::endl;
	}

	hvol = dvolume;
	std::ofstream fou("DEMO18v2_reconstruction.raw", std::ios::binary);
	fou.write((char*) (&hvol[0]), sizeof(float) * volReso);
	fou.close();
	std::cout << "Finished\n";
}

void DEMO18v3()
{

	//cudaSetDevice(2);
	///// Statistical Iterative Reconstruction OS-SART ÖØ¹¹;
	//ConeEDGeo ConeGeo(141.5261, 185.0329 - 141.5261, 400, 0.0f, 3.4906585039886591538473777777778, make_float2(60.0, 60.0), make_int2(1200, 1200));
	//ConeGeo.m_DetCntIdx.x = ConeGeo.m_DetN.x * 0.5 - 1.8766;
	//ConeGeo.m_DetCntIdx.y = ConeGeo.m_DetN.y * 0.5 - 0.1951;
	//const unsigned int projR = ConeGeo.m_DetN.x * ConeGeo.m_DetN.y;
	//Volume Vol(512, 512, 512, 39.168, 39.168, 39.168, 0.0, 0.0, 0.0);

	//cuint volReso = Vol.m_Reso.x * Vol.m_Reso.y * Vol.m_Reso.z;
	//const unsigned int prjReso = projR * ConeGeo.m_ViwN;

	//thrust::host_vector<float> hvol(volReso, 0);
	//thrust::host_vector<float> hpho(prjReso, 1.0);
	//float* hprj = new float[prjReso];// thrust::host_vector<float> hprj(prjReso, 0);
	////thrust::host_vector<float> hmsk(volReso, 1.0f); //È«1 mask

	//std::ifstream fid("DEMO18v3_projectionData.prj", std::ios::binary);
	//if (!fid.is_open())
	//{
	//	std::cerr << "Cannot open the projection file for reconstruction\n";
	//	exit(-1);
	//}
	//fid.read((char*)hprj, sizeof(float) * prjReso);

	////cudaSetDevice(2);
	////float* dprj;
	////cudaMalloc((void**) &dprj, sizeof(float) * prjReso);
	////cudaMemcpy(dprj, hprj, sizeof(float)*prjReso, cudaMemcpyHostToDevice);
	////cudaSetDevice(0);


	////Ò»ŽÎÍ¶Ó°Ò»žöœÇ¶È;
	//dim3 prjBlk(32, 32);
	//dim3 prjGid(
	//	(ConeGeo.m_DetN.x + prjBlk.x - 1) / prjBlk.x,
	//	(ConeGeo.m_DetN.y + prjBlk.y - 1) / prjBlk.y);

	//dim3 volBlk(8, 8, 8);
	//dim3 volGid(
	//	(Vol.m_Reso.x + volBlk.x - 1) / volBlk.x,
	//	(Vol.m_Reso.y + volBlk.y - 1) / volBlk.y,
	//	(Vol.m_Reso.z + volBlk.z - 1) / volBlk.z);

	//thrust::host_vector<float>::iterator maxYIdx = thrust::max_element(hpho.begin(), hpho.end());
	//const float maxY = (*maxYIdx);
	//thrust::device_vector<float> dbakLambda(volReso, 0.0);
	//thrust::device_vector<float> dalloneproj(projR, 1.0);
	//thrust::device_vector<float> drawproj(projR, 0.0);
	//thrust::device_vector<float> dweightedproj(projR, 0.0);
	//thrust::device_vector<float> dvolume(volReso, 0);
	//thrust::device_vector<float> dupdVol(volReso, 0);


	//float* pbakLambda = thrust::raw_pointer_cast(&dbakLambda[0]);
	//float* palloneproj = thrust::raw_pointer_cast(&dalloneproj[0]);
	//float* prawproj = thrust::raw_pointer_cast(&drawproj[0]);
	//float* pvolume = thrust::raw_pointer_cast(&dvolume[0]);
	//float* pupdVol = thrust::raw_pointer_cast(&dupdVol[0]);
	//float* pweightedproj = thrust::raw_pointer_cast(&dweightedproj[0]);

	//thrust::host_vector<float> hv = dbakLambda;
	//const float hfL = Vol.m_Reso.x * 0.5;
	//const float hfW = Vol.m_Reso.y * 0.5;
	////ŒÆËãÈ«1µÄ·ŽÍ¶Ó°;
	//std::ifstream weiFile("DEMO18v3_weight.raw", std::ios::binary);
	//if (!weiFile.is_open())
	//{
	//	std::cout << "Begin all one backprojction for weighting\n";
	//	for (unsigned int i = 0; i != ConeGeo.m_ViwN; ++i)
	//	{
	//		//back_proj_DEMO18v3(palloneproj, pbakLambda, ConeGeo, Vol, i, volBlk, volGid);
	//		cudaDeviceSynchronize();
	//		std::cout << i << std::endl;
	//	}
	//	hv = dbakLambda;
	//	for (unsigned int sliceIdx = 0; sliceIdx != Vol.m_Reso.z; ++sliceIdx)
	//	{
	//		for (unsigned int jj = 0; jj != Vol.m_Reso.y; ++jj)
	//		{
	//			for (unsigned int ii = 0; ii != Vol.m_Reso.x; ++ii)
	//			{
	//				unsigned int i = (sliceIdx * Vol.m_Reso.y + jj) * Vol.m_Reso.x + ii;
	//				if (((ii - hfL + 0.5) / hfL)*((ii - hfL + 0.5) / hfL) + ((jj - hfW + 0.5) / hfL)*((jj - hfW + 0.5) / hfW) >= 1.44)
	//				{
	//					hv[i] = 0;
	//				}
	//			}
	//		}
	//	}
	//	std::ofstream weight("DEMO18v3_weight.raw", std::ios::binary);
	//	weight.write((char*)(&hv[0]), sizeof(float) * volReso);
	//	weight.close();
	//}
	//else
	//{
	//	weiFile.read((char*)(&hv[0]), sizeof(float) * volReso);
	//	weiFile.close();
	//}

	//dbakLambda = hv;
	//std::cout << "Finished the backprojction\n";
	//unsigned int iterIdx = 0;
	//unsigned int subsetIdx = 0;
	//unsigned int subangIdx = 0;
	//unsigned int angIdx = 0;
	//cuint iterNum = 1; // We use 20 iterations.
	//int subSetNum = 20;
	//int numPerSubSet = ConeGeo.m_ViwN / subSetNum;

	//for (iterIdx = 0; iterIdx != iterNum; ++iterIdx)
	//{
	//	if (iterIdx == 0)
	//	{
	//		subSetNum = 100;
	//		numPerSubSet = ConeGeo.m_ViwN / subSetNum;
	//	}
	//	if (iterIdx == 1)
	//	{
	//		subSetNum = 80;
	//		numPerSubSet = ConeGeo.m_ViwN / subSetNum;
	//	}
	//	if (iterIdx == 2)
	//	{
	//		subSetNum = 50;
	//		numPerSubSet = ConeGeo.m_ViwN / subSetNum;
	//	}
	//	if (iterIdx == 3)
	//	{
	//		subSetNum = 40;
	//		numPerSubSet = ConeGeo.m_ViwN / subSetNum;
	//	}
	//	if (iterIdx == 4)
	//	{
	//		subSetNum = 40;
	//		numPerSubSet = ConeGeo.m_ViwN / subSetNum;
	//	}
	//	if (iterIdx == 5)
	//	{
	//		subSetNum = 20;
	//		numPerSubSet = ConeGeo.m_ViwN / subSetNum;
	//	}
	//	if (iterIdx == 6)
	//	{
	//		subSetNum = 20;
	//		numPerSubSet = ConeGeo.m_ViwN / subSetNum;
	//	}
	//	if (iterIdx == 7)
	//	{
	//		subSetNum = 20;
	//		numPerSubSet = ConeGeo.m_ViwN / subSetNum;
	//	}
	//	if (iterIdx == 8)
	//	{
	//		subSetNum = 10;
	//		numPerSubSet = ConeGeo.m_ViwN / subSetNum;
	//	}
	//	if (iterIdx == 9)
	//	{
	//		subSetNum = 10;
	//		numPerSubSet = ConeGeo.m_ViwN / subSetNum;
	//	}
	//	if (iterIdx == 10)
	//	{
	//		subSetNum = 10;
	//		numPerSubSet = ConeGeo.m_ViwN / subSetNum;
	//	}
	//	if (iterIdx > 10 && iterIdx < 16)
	//	{
	//		subSetNum = 5;
	//		numPerSubSet = ConeGeo.m_ViwN / subSetNum;
	//	}
	//	if (iterIdx >= 16)
	//	{
	//		subSetNum = 2;
	//		numPerSubSet = ConeGeo.m_ViwN / subSetNum;
	//	}

	//	for (subsetIdx = 0; subsetIdx != subSetNum; ++subsetIdx)
	//	{
	//		//Set the volume as 0;
	//		thrust::fill(dupdVol.begin(), dupdVol.end(), 0.0f);
	//		//Calculate current angle;
	//		for (subangIdx = 0; subangIdx != numPerSubSet; ++subangIdx)
	//		{
	//			angIdx = subsetIdx + subangIdx *subSetNum;
	//			cudaMemcpy(prawproj, hprj + angIdx *projR, sizeof(float) * projR, cudaMemcpyHostToDevice);
	//			//cudaMemcpy(prawproj, dprj + angIdx *projR, sizeof(float) * projR, cudaMemcpyDeviceToDevice);
	//			proj_DEMO17(pvolume, pweightedproj, prawproj, ConeGeo, Vol, angIdx, prjBlk, prjGid);
	//			//bakprojection
	//			//back_proj_DEMO17(pweightedproj, pupdVol, ConeGeo, Vol, angIdx, volBlk, volGid);
	//			back_proj_DEMO18v3(pweightedproj, pupdVol, ConeGeo, Vol, angIdx, volBlk, volGid);
	//			std::cout << "angIdx: " << angIdx << std::endl;
	//		}

	//		//Update the volume
	//		thrust::transform(
	//			thrust::make_zip_iterator(
	//			thrust::make_tuple(dvolume.begin(), dupdVol.begin(), dbakLambda.begin())
	//			),
	//			thrust::make_zip_iterator(
	//			thrust::make_tuple(dvolume.end(), dupdVol.end(), dbakLambda.end())
	//			),
	//			dvolume.begin(), _DEMO18_update_functor< float >(subSetNum, 0, 1.2));

	//	}
	//	std::cout << "Iteration: " << iterIdx << std::endl;
	//}

	//hvol = dvolume;
	//std::ofstream fou("DEMO18v3_reconstruction.raw", std::ios::binary);
	//fou.write((char*)(&hvol[0]), sizeof(float) * volReso);
	//fou.close();
	//std::cout << "Finished\n";
}

// It is right
void DEMO18v4_2D()
{
	FanEDGeo FanGeo(141.5261, 185.0329 - 141.5261, 400, 0.0, 0.5 * 400 / 180 * 3.14159265358979323846264, 60.0, 1200);
	Image Img(512, 512, 39.168, 39.168, 0.0, 0.0);
	cuint imgReso = Img.m_Reso.x * Img.m_Reso.y;
	cuint prjReso = FanGeo.m_DetN * FanGeo.m_ViwN;

	thrust::host_vector<float> himg(imgReso, 0);
	thrust::host_vector<float> hpho(prjReso, 0); // Represent Y0
	thrust::host_vector<float> hprj(prjReso, 0); // Represent P
	thrust::host_vector<float> hmsk(imgReso, 1.0f);

	//¶ÁÈ¡ÎÄŒþ;
	std::ifstream fid("DEMO18v42D_mouseMidProj.prj", std::ios::binary);
	std::ifstream fi2("DEMO18v42D_mouseMidPhot.pho", std::ios::binary);
	if (!(fid.is_open() && fi2.is_open()))
	{
		std::cerr << "Cannot open the projection file and the photon weighting file\n";
		exit(-1);
	}
	fid.read((char*) (&hprj[0]), sizeof(float) * prjReso);
	fi2.read((char*) (&hpho[0]), sizeof(float) * prjReso);

	for (unsigned int i = 0; i != imgReso; ++i)
	{
		hmsk[i] = 1.0f;
	}

	cuint iterNum = 300;

	cuint subSetNum = 1;
	cuint numPerSubSet = FanGeo.m_ViwN / subSetNum;
	unsigned int realIters = subSetNum * iterNum;


	dim3 prjBlk(256, 4);
	dim3 prjGid(
		(FanGeo.m_DetN + prjBlk.x - 1) / prjBlk.x,
		(numPerSubSet + prjBlk.y - 1) / prjBlk.y);
	dim3 bakBlk(32, 32);
	dim3 bakGid(
		(Img.m_Reso.x + bakBlk.x - 1) / bakBlk.x,
		(Img.m_Reso.y + bakBlk.y - 1) / bakBlk.y);

	thrust::host_vector<float>::iterator maxYidx = thrust::max_element(hpho.begin(), hpho.end());
	std::cout << "The maximum photon number for this projection is " << *maxYidx << std::endl;
	const float maxY = (*maxYidx);

	thrust::device_vector<float> dimg = himg; //Image to be reconstructed;
	thrust::device_vector<float> dprj = hprj; //Projection data;
	thrust::device_vector<float> dpho = hpho; //Photon number 
	thrust::device_vector<float> dmsk = hmsk; //Mask Image;
	thrust::device_vector<float> dprojLambda(prjReso, 0); // Lambda for projection weighting
	thrust::device_vector<float> dbackLambda(imgReso, 0); // Lambda for backprojection weighting
	thrust::device_vector<float> dprjdiff(prjReso, 0);
	thrust::device_vector<float> dupdateImg(imgReso, 0);
	thrust::device_vector<float> doneproj(prjReso, 0);
	thrust::fill(doneproj.begin(), doneproj.end(), 1.0f);


	float* pmsk = thrust::raw_pointer_cast(&dmsk[0]);
	float* pprojLambda = thrust::raw_pointer_cast(&dprojLambda[0]);
	float* pbackLambda = thrust::raw_pointer_cast(&dbackLambda[0]);
	float* pimg = thrust::raw_pointer_cast(&dimg[0]);
	float* pprj = thrust::raw_pointer_cast(&dprj[0]);
	float* pprjdiff = thrust::raw_pointer_cast(&dprjdiff[0]);
	float* pupdateImg = thrust::raw_pointer_cast(&dupdateImg[0]);
	float* ppho = thrust::raw_pointer_cast(&dpho[0]);
	float* poneproj = thrust::raw_pointer_cast(&doneproj[0]);

	//Generate the Lambda matrix for projection weighting
	proj_DEMO18v4_2D(pmsk, pprojLambda, FanGeo, Img, prjBlk, prjGid);
	thrust::transform(dpho.begin(), dpho.end(), dprojLambda.begin(), dprojLambda.begin(), _DEMO18v4_GenProjLambda_functor<float>(*maxYidx));

	thrust::fill(dprjdiff.begin(), dprjdiff.end(), 1.0f);
	bakproj_PIXEL_DEMO18v4_2D(pprjdiff, pbackLambda, pmsk, FanGeo, Img, bakBlk, bakGid);
	thrust::transform(dbackLambda.begin(), dbackLambda.end(), dbackLambda.begin(), _DEMO18v4_GenBackLambda_functor<float>(1.0f));

	const float lambdas = 1.5f;
	for (unsigned int i = 0; i != realIters; ++i)
	{
		//cudaThreadSynchronize();
		proj_DEMO18v4_2D(pimg, pprjdiff, FanGeo, Img, prjBlk, prjGid);
		dprjdiff = dprjdiff - dprj;
		thrust::transform(dprojLambda.begin(), dprojLambda.end(), dprjdiff.begin(), dprjdiff.begin(), thrust::multiplies<float>());
		thrust::fill(dupdateImg.begin(), dupdateImg.end(), 0.0f);
		bakproj_PIXEL_DEMO18v4_2D(pprjdiff, pupdateImg, pmsk, FanGeo, Img, bakBlk, bakGid);
		thrust::transform(dbackLambda.begin(), dbackLambda.end(), dupdateImg.begin(), dupdateImg.begin(), thrust::multiplies<float>());
		thrust::transform(dimg.begin(), dimg.end(), dupdateImg.begin(), dimg.begin(), _DEMO18v4_updateImgSIR_functor<float>(lambdas));
	}

	himg = dimg;
	std::ofstream fou("DEMO18v42D_SIRreconImg300.raw", std::ios::binary);
	if (!fou.is_open())
	{
		std::cout << "Cannot write the file\n";
		exit(-1);
	}
	fou.write((char*) (&himg[0]), sizeof(float) * imgReso);
	fou.close();
}