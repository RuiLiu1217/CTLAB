#include <iostream>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "FastMatrixVectorMultiplication.hpp"
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
	cudaMallocArray(&d_A, &channelDesc, width >> 2, height);
	cudaMemcpy2DToArray(d_A, 0, 0, A, n * sizeof(float), n * sizeof(float), m, cudaMemcpyHostToDevice);
	cudaBindTextureToArray(texRefA, d_A);
	cudaMalloc((void**) &d_x, n * sizeof(float));
	cudaMalloc((void**) &d_y, m * sizeof(float));

	cudaMemcpy(d_x, x, n * sizeof(float), cudaMemcpyHostToDevice);
	FMVM_KER << < grid, threads >> >(d_y, d_A, d_x, m, n);
	cudaMemcpy(y, d_y, m * sizeof(float), cudaMemcpyDeviceToHost);
	cudaFree(d_y);
	cudaFree(d_x);
	cudaUnbindTexture(texRefA);
	cudaFreeArray(d_A);

}

void FMVM_device(float* d_y, float* dA, float* d_x, int m, int n)
{
	int blkNum((m >> 4) + ((m & 15) ? 1 : 0));
	int height(blkNum << 4);
	int width((n & 255) ? (256 * ((n >> 8) + 1)) : n);
	dim3 threads(16, 16), grid(blkNum, 1);
	cudaArray *d_A;

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float4>();
	cudaMallocArray(&d_A, &channelDesc, width >> 2, height);
	cudaMemcpy2DToArray(d_A, 0, 0, dA, n * sizeof(float), n * sizeof(float), m, cudaMemcpyDeviceToDevice);
	cudaBindTextureToArray(texRefA, d_A);

	FMVM_KER << < grid, threads >> >(d_y, d_A, d_x, m, n);
	cudaUnbindTexture(texRefA);
	cudaFreeArray(d_A);
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
	cudaMalloc((void**) &d_x, sizeof(float) * 3);
	cudaMalloc((void**) &d_y, sizeof(float) * 3);
	cudaMalloc((void**) &d_A, sizeof(float) * 9);
	cudaMemcpy(d_x, x, sizeof(float) * 3, cudaMemcpyHostToDevice);
	cudaMemcpy(d_y, y, sizeof(float) * 3, cudaMemcpyHostToDevice);
	cudaMemcpy(d_A, A, sizeof(float) * 9, cudaMemcpyHostToDevice);

	FMVM_device(d_y, d_A, d_x, 3, 3);
	cudaMemcpy(y, d_y, sizeof(float) * 3, cudaMemcpyDeviceToHost);
	std::cout << y[0] << " " << y[1] << " " << y[2];
}