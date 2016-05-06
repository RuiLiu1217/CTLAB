#include "differentialFuncs.hpp"

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS 1
#endif

#include <thrust/device_vector.h>
#include <thrust/extrema.h>
//#include <thrust/pair.h>
#ifndef EPSILON
#define EPSILON 1.0E-9
#endif

typedef const unsigned int cuint;
template<typename T>
__global__ void _Dx_ker(T* f, T* d, const T coef, cuint L, cuint W)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	if (idx < L && idy < W)
	{
		cuint curid = idy * L + idx;
		cuint lasid = idy * L + ((idx + L - 1) % L);
		d[curid] = (f[curid] - f[lasid]) * coef;
		return;
	}
}

template<typename T>
__global__ void _D2x_ker(T* f, T* d, const T coef, cuint L, cuint W)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	if (idx < L && idy < W)
	{
		cuint curid = idy * L + idx;
		cuint lasid = idy * L + ((idx + L - 1) % L);
		cuint las2id = idy * L + ((idx + L - 2) % L);
		d[curid] = (f[curid] + f[las2id] - f[lasid] * 2.0) * coef;
		return;
	}
}


template<typename T>
__global__ void _Dx_ker(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	cuint idz = threadIdx.z + blockIdx.z * blockDim.z;
	if (idx < L && idy < W && idz < H)
	{
		cuint curid = (idz * W + idy) * L + idx;
		cuint lasid = (idz * W + idy) * L + ((idx + L - 1) % L);
		d[curid] = (f[curid] - f[lasid]) * coef;
		return;
	}
}


template<typename T>
__global__ void _D2x_ker(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	cuint idz = threadIdx.z + blockIdx.z * blockDim.z;
	if (idx < L && idy < W && idz < H)
	{
		cuint curid = (idz * W + idy) * L + idx;
		cuint lasid = (idz * W + idy) * L + ((idx + L - 1) % L);
		cuint las2id = (idz * W + idy) * L + ((idx + L - 2) % L);
		d[curid] = (f[curid] + f[las2id] - f[lasid] * 2.0) * coef;
		return;
	}
}

template<typename T>
__global__ void _Dxt_ker(T* f, T* d, const T coef, cuint L, cuint W)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	if (idx < L && idy < W)
	{
		cuint curid = idy * L + idx;
		cuint nexid = idy * L + ((idx + 1) % L);
		d[curid] = (f[curid] - f[nexid]) * coef;
		return;
	}
}

template<typename T>
__global__ void _D2xt_ker(T* f, T* d, const T coef, cuint L, cuint W)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	if (idx < L && idy < W)
	{
		cuint curid = idy * L + idx;
		cuint nexid = idy * L + ((idx + 1) % L);
		cuint nex2id = idy * L + ((idx + 2) % L);
		d[curid] = (f[curid] + f[nex2id] - f[nexid] * 2.0) * coef;
		return;
	}
}

template<typename T>
__global__ void _Dxt_ker(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	cuint idz = threadIdx.z + blockIdx.z * blockDim.z;
	if (idx < L && idy < W && idz < H)
	{
		cuint curid = (idz * W + idy) * L + idx;
		cuint nexid = (idz * W + idy) * L + ((idx + 1) % L);
		d[curid] = (f[curid] - f[nexid]) * coef;
		return;
	}
}

template<typename T>
__global__ void _D2xt_ker(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	cuint idz = threadIdx.z + blockIdx.z * blockDim.z;
	if (idx < L && idy < W && idz < H)
	{
		cuint curid = (idz * W + idy) * L + idx;
		cuint nexid = (idz * W + idy) * L + ((idx + 1) % L);
		cuint nex2id = (idz * W + idy) * L + ((idx + 2) % L);
		d[curid] = (f[curid] + f[nex2id] - f[nexid] * 2.0) * coef;
		return;
	}
}


template<typename T>
__global__ void _Dy_ker(T* f, T* d, const T coef, cuint L, cuint W)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	if (idx < L && idy < W)
	{
		cuint curid = idy * L + idx;
		cuint lasid = ((idy + W - 1) % W) * L + idx;
		d[curid] = (f[curid] - f[lasid]) * coef;
		return;
	}
}

template<typename T>
__global__ void _D2y_ker(T* f, T* d, const T coef, cuint L, cuint W)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	if (idx < L && idy < W)
	{
		cuint curid = idy * L + idx;
		cuint lasid = ((idy + W - 1) % W) * L + idx;
		cuint las2id = ((idy + W - 2) % W) * L + idx;
		d[curid] = (f[curid] + f[las2id] - f[lasid] * 2.0) * coef;
		return;
	}
}

template<typename T>
__global__ void _Dy_ker(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	cuint idz = threadIdx.z + blockIdx.z * blockDim.z;
	if (idx < L && idy < W && idz < H)
	{
		cuint curid = (idz * W + idy) * L + idx;
		cuint lasid = (idz * W + ((idy + W - 1) % W)) * L + idx;
		d[curid] = (f[curid] - f[lasid]) * coef;
		return;
	}
}


template<typename T>
__global__ void _D2y_ker(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	cuint idz = threadIdx.z + blockIdx.z * blockDim.z;
	if (idx < L && idy < W && idz < H)
	{
		cuint curid = (idz * W + idy) * L + idx;
		cuint lasid = (idz * W + ((idy + W - 1) % W)) * L + idx;
		cuint las2id = (idz * W + ((idy + W - 2) % W)) * L + idx;
		d[curid] = (f[curid] + f[las2id] - f[lasid] * 2.0) * coef;
		return;
	}
}


template<typename T>
__global__ void _Dyt_ker(T* f, T* d, const T coef, cuint L, cuint W)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	if (idx < L && idy < W)
	{
		cuint curid = idy * L + idx;
		cuint nexid = ((idy + 1) % W) * L + idx;
		d[curid] = (f[curid] - f[nexid]) * coef;
		return;
	}
}



template<typename T>
__global__ void _D2yt_ker(T* f, T* d, const T coef, cuint L, cuint W)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	if (idx < L && idy < W)
	{
		cuint curid = idy * L + idx;
		cuint nexid = ((idy + 1) % W) * L + idx;
		cuint nex2id = ((idy + 2) % W) * L + idx;
		d[curid] = (f[curid] + f[nex2id] - f[nexid] * 2.0) * coef;
		return;
	}
}




template<typename T>
__global__ void _Dyt_ker(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	cuint idz = threadIdx.z + blockIdx.z * blockDim.z;
	if (idx < L && idy < W && idz < H)
	{
		cuint curid = (idz * W + idy) * L + idx;
		cuint nexid = (idz * W + ((idy + 1) % W)) * L + idx;
		d[curid] = (f[curid] - f[nexid]) * coef;
		return;
	}
}




template<typename T>
__global__ void _D2yt_ker(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	cuint idz = threadIdx.z + blockIdx.z * blockDim.z;
	if (idx < L && idy < W && idz < H)
	{
		cuint curid = (idz * W + idy) * L + idx;
		cuint nexid = (idz * W + ((idy + 1) % W)) * L + idx;
		cuint nex2id = (idz * W + ((idy + 2) % W)) * L + idx;
		d[curid] = (f[curid] + f[nex2id] - f[nexid] * 2.0) * coef;
		return;
	}
}



template<typename T>
__global__ void _Dz_ker(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	cuint idz = threadIdx.z + blockIdx.z * blockDim.z;
	if (idx < L && idy < W && idz < H)
	{
		cuint curid = (idz * W + idy) * L + idx;
		cuint lasid = (((idz + H - 1) % H) * W + idy) * L + idx;
		d[curid] = (f[curid] - f[lasid]) * coef;
		return;
	}
}



template<typename T>
__global__ void _D2z_ker(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	cuint idz = threadIdx.z + blockIdx.z * blockDim.z;
	if (idx < L && idy < W && idz < H)
	{
		cuint curid = (idz * W + idy) * L + idx;
		cuint lasid = (((idz + H - 1) % H) * W + idy) * L + idx;
		cuint las2id = (((idz + H - 2) % H) * W + idy) * L + idx;
		d[curid] = (f[curid] + f[las2id] - f[lasid] * 2.0) * coef;
		return;
	}
}


template<typename T>
__global__ void _Dzt_ker(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	cuint idz = threadIdx.z + blockIdx.z * blockDim.z;
	if (idx < L && idy < W && idz < H)
	{
		cuint curid = (idz * W + idy) * L + idx;
		cuint nexid = (((idz + 1) % H) * W + idy) * L + idx;
		d[curid] = (f[curid] - f[nexid]) * coef;
		return;
	}
}

template<typename T>
__global__ void _D2zt_ker(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	cuint idz = threadIdx.z + blockIdx.z * blockDim.z;
	if (idx < L && idy < W && idz < H)
	{
		cuint curid = (idz * W + idy) * L + idx;
		cuint nexid = (((idz + 1) % H) * W + idy) * L + idx;
		cuint nex2id = (((idz + 2) % H) * W + idy) * L + idx;
		d[curid] = (f[curid] + f[nex2id] - f[nexid] * 2.0) * coef;
		return;
	}
}


template<typename T>
__global__ void _Laplacian_ker(T* f, T* l, const T coef, cuint L, cuint W)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	if (idx < L && idy < W)
	{
		cuint curid = idy * L + idx;
		cuint lasxd = idy * L + ((idx + L - 1) % L);
		cuint nexxd = idy * L + ((idx + 1) % L);
		cuint lasyd = ((idy + W - 1) % W) * L + idx;
		cuint nexyd = ((idy + 1) % W) * L + idx;
		l[curid] = (4.0 * f[curid] - f[lasxd] - f[lasyd] - f[nexxd] - f[nexyd]) * coef;
	}
}



template<typename T>
__global__ void _Laplacian_ker(T* f, T* l, const T coef, cuint L, cuint W, cuint H)
{
	cuint idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuint idy = threadIdx.y + blockIdx.y * blockDim.y;
	cuint idz = threadIdx.z + blockIdx.z * blockDim.z;
	if (idx < L && idy < W && idz < H)
	{
		cuint curid = (idz * W + idy) * L + idx;
		cuint lasxd = (idz * W + idy) * L + ((idx + L - 1) % L);
		cuint nexxd = (idz * W + idy) * L + ((idx + 1) % L);
		cuint lasyd = (idz * W + ((idy + W - 1) % W)) * L + idx;
		cuint nexyd = (idz * W + ((idy + 1) % W)) * L + idx;
		cuint laszd = (((idz + H - 1) % H) * W + idy) * L + idx;
		cuint nexzd = (((idz + 1) % H) * W + idy) * L + idx;
		l[curid] = (6.0 * f[curid] - f[lasxd] - f[lasyd] - f[nexxd] - f[nexyd] - f[laszd] - f[nexzd]) * coef;
	}
}


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
		d[curid] = sqrt(difx * difx + dify * dify) * coef;
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
		cuint nexxd = (idz * W + idy) * L + ((idx + 1) % L);
		cuint nexyd = (idz * W + ((idy + 1) % W)) * L + idx;
		cuint nexzd = (((idz + 1) % H) * W + idy) * L + idx;
		const T difx = f[nexxd] - f[curid];
		const T dify = f[nexyd] - f[curid];
		const T difz = f[nexzd] - f[curid];
		d[curid] = sqrt(difx * difx + dify * dify + difz * difz) * coef;
		return;
	}
}


template<typename T>
__global__ void _DiscreteGradientTrans_Shared_ker(T* img, T* tvimg, T coef, cuint imgL, cuint imgW)
{
	cuint idxX = threadIdx.x + blockIdx.x * (blockDim.x - 1);
	cuint idxY = threadIdx.y + blockIdx.y * (blockDim.y - 1);
	if (idxX < imgL && idxY < imgW)
	{
		cuint idx = idxY * imgL + idxX;
		T res(0);
		__shared__ T svol[32][32];
		svol[threadIdx.y][threadIdx.x] = img[idx];
		__syncthreads();
		if (threadIdx.x == 31 || threadIdx.y == 31 || threadIdx.z == 31)
			return;
		res = sqrt(
			(svol[threadIdx.y][threadIdx.x] - svol[threadIdx.y][threadIdx.x + 1]) * (svol[threadIdx.y][threadIdx.x] - svol[threadIdx.y][threadIdx.x + 1]) +
			(svol[threadIdx.y][threadIdx.x] - svol[threadIdx.y + 1][threadIdx.x]) * (svol[threadIdx.y][threadIdx.x] - svol[threadIdx.y + 1][threadIdx.x]));
		tvimg[idx] = res * coef;
		if (idxX == imgL - 1 || idxY == imgW - 1 || idxX == 0 || idxY == 0)
		{
			tvimg[idx] = 0;
		}
	}
}
template<typename T>
void DiscreteGradientTrans_Shared(T* vol, T* tvvol, T coef, cuint imgL, cuint imgW)
{
	dim3 blk(32, 32);
	dim3 gid((imgL + 30) / 31, (imgW + 30) / 31);
	_DiscreteGradientTrans_Shared_ker<T> << <gid, blk >> >(vol, tvvol, coef, imgL, imgW);
}
void DiscreteGradientTrans_Shared_GPU(float* vol, float* tvvol, float coef, cuint imgL, cuint imgW)
{
	DiscreteGradientTrans_Shared<float>(vol, tvvol, coef, imgL, imgW);
}
void DiscreteGradientTrans_Shared_GPU(double* vol, double* tvvol, double coef, cuint imgL, cuint imgW)
{
	DiscreteGradientTrans_Shared<double>(vol, tvvol, coef, imgL, imgW);
}


template<typename T>
__global__ void _DiscreteGradientTrans_Shared_ker(T* vol, T* tvvol, T coef, cuint imgL, cuint imgW, cuint imgH)
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
		res = sqrt(
			(svol[threadIdx.z][threadIdx.y][threadIdx.x] - svol[threadIdx.z][threadIdx.y][threadIdx.x + 1]) * (svol[threadIdx.z][threadIdx.y][threadIdx.x] - svol[threadIdx.z][threadIdx.y][threadIdx.x + 1]) +
			(svol[threadIdx.z][threadIdx.y][threadIdx.x] - svol[threadIdx.z][threadIdx.y + 1][threadIdx.x]) * (svol[threadIdx.z][threadIdx.y][threadIdx.x] - svol[threadIdx.z][threadIdx.y + 1][threadIdx.x]) +
			(svol[threadIdx.z][threadIdx.y][threadIdx.x] - svol[threadIdx.z + 1][threadIdx.y][threadIdx.x]) * (svol[threadIdx.z][threadIdx.y][threadIdx.x] - svol[threadIdx.z + 1][threadIdx.y][threadIdx.x]));
		tvvol[idx] = res * coef;
		if (idxX == imgL - 1 || idxY == imgW - 1 || idxZ == imgH - 1
			|| idxX == 0 || idxY == 0 || idxZ == 0)
		{
			tvvol[idx] = 0;
		}
	}
}
template<typename T>
void DiscreteGradientTrans_Shared(T* vol, T* tvvol, T coef, cuint imgL, cuint imgW, cuint imgH)
{
	dim3 blk(8, 8, 8);
	dim3 gid(
		(imgL + 6) / 7,
		(imgW + 6) / 7,
		(imgH + 6) / 7);
	_DiscreteGradientTrans_Shared_ker<T> << <gid, blk >> >(vol, tvvol, coef, imgL, imgW, imgH);
}
void DiscreteGradientTrans_Shared_GPU(float* vol, float* tvvol, float coef, cuint imgL, cuint imgW, cuint imgH)
{
	DiscreteGradientTrans_Shared<float>(vol, tvvol, coef, imgL, imgW, imgH);
}
void DiscreteGradientTrans_Shared_GPU(double* vol, double* tvvol, double coef, cuint imgL, cuint imgW, cuint imgH)
{
	DiscreteGradientTrans_Shared<double>(vol, tvvol, coef, imgL, imgW, imgH);
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
		const T dom1 = 1.0 / (sqrt((fi1j - fij)*(fi1j - fij) + (fij1 - fij) * (fij1 - fij) + EPSILON));
		const T dom2 = 1.0 / (sqrt((fij - fi_1j)*(fij - fi_1j) + (fi_1j1 - fi_1j)*(fi_1j1 - fi_1j) + EPSILON));
		const T dom3 = 1.0 / (sqrt((fi1j_1 - fij_1) * (fi1j_1 - fij_1) + (fij - fij_1)*(fij - fij_1) + EPSILON));
		d[curid] = ((2.0 * fij - fi1j - fij1) * dom1 + (fij - fi_1j) * dom2 + (fij - fij_1) * dom3) * coef;
		return;
	}
}




template<typename T>
__global__ void _GradientOfTV_ker(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	cuint idx = threadIdx.x + blockDim.x * blockIdx.x;
	cuint idy = threadIdx.y + blockDim.y * blockIdx.y;
	cuint idz = threadIdx.z + blockDim.z * blockIdx.z;
	if (idx < L && idy < W && idz < H)
	{
		cuint curId = (idz * W + idy) * L + idx;
		const T f_ijk = f[(idz * W + idy) * L + idx];
		const T f_i1jk = f[(idz * W + idy) * L + ((idx + 1) % L)];
		const T f_ij1k = f[(idz * W + ((idy + 1) % W)) * L + idx];
		const T f_ijk1 = f[(((idz + 1) % H) * W + idy) * L + idx];
		const T f_i_1jk = f[(idz * W + idy) * L + ((idx + L - 1) % L)];
		const T f_i_1j1k = f[(idz * W + ((idy + 1) % W)) * L + ((idx + L - 1) % L)];
		const T f_i_1jk1 = f[(((idz + 1) % H) * W + idy) * L + ((idx + L - 1) % L)];
		const T f_ij_1k = f[(idz * W + ((idy + W - 1) % W)) * L + idx];
		const T f_i1j_1k = f[(idz * W + ((idy + W - 1) % W)) * L + ((idx + 1) % L)];
		const T f_ij_1k1 = f[(((idz + 1) % H) * W + ((idy + W - 1) % W)) * L + idx];
		const T f_ijk_1 = f[(((idz + H - 1) % H) * W + idy) * L + idx];
		const T f_i1jk_1 = f[(((idz + H - 1) % H) * W + idy) * L + ((idx + 1) % L)];
		const T f_ij1k_1 = f[(((idz + H - 1) % H) * W + ((idy + 1) % W)) * L + idx];

		const T dom1 = 1.0 / sqrt((f_i1jk - f_ijk)*(f_i1jk - f_ijk) + (f_ij1k - f_ijk)*(f_ij1k - f_ijk) + (f_ijk1 - f_ijk)*(f_ijk1 - f_ijk) + EPSILON);
		const T dom2 = 1.0 / sqrt((f_ijk - f_i_1jk)*(f_ijk - f_i_1jk) + (f_i_1j1k - f_i_1jk)*(f_i_1j1k - f_i_1jk) + (f_i_1jk1 - f_i_1jk)*(f_i_1jk1 - f_i_1jk) + EPSILON);
		const T dom3 = 1.0 / sqrt((f_i1j_1k - f_ij_1k)*(f_i1j_1k - f_ij_1k) + (f_ijk - f_ij_1k)*(f_ijk - f_ij_1k) + (f_ij_1k1 - f_ij_1k)*(f_ij_1k1 - f_ij_1k) + EPSILON);
		const T dom4 = 1.0 / sqrt((f_i1jk_1 - f_ijk_1)*(f_i1jk_1 - f_ijk_1) + (f_ij1k_1 - f_ijk_1)*(f_ij1k_1 - f_ijk_1) + (f_ijk - f_ijk_1)*(f_ijk - f_ijk_1) + EPSILON);

		d[curId] = ((3.0 * f_ijk - f_i1jk - f_ij1k - f_ijk1) * dom1 + (f_ijk - f_i_1jk) * dom2 + (f_ijk - f_ij_1k) * dom3 + (f_ijk - f_ijk_1) * dom4) * coef;
		return;
	}
}



//��άͼ���� 32x32, ���õ�ֻ��30x30
template<typename T>
__global__ void descentDir_SHARED_KER(T* img, T* des,
	cuint imgL, cuint imgW)
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

		cuint y = threadIdx.y;
		cuint x = threadIdx.x;
		T val =
			((2.0 * simg[y][x] - simg[y][x - 1] - simg[y - 1][x]) / sqrt(1.0e-36 + (simg[y][x] - simg[y][x - 1])*(simg[y][x] - simg[y][x - 1]) + (simg[y][x] - simg[y - 1][x])*(simg[y][x] - simg[y - 1][x])))
			- ((simg[y][x + 1] - simg[y][x]) / sqrt(1.0e-36 + (simg[y][x + 1] - simg[y][x]) * (simg[y][x + 1] - simg[y][x]) + (simg[y][x + 1] - simg[y - 1][x + 1]) * (simg[y][x + 1] - simg[y - 1][x + 1])))
			- ((simg[y + 1][x] - simg[y][x]) / sqrt(1.0e-36 + (simg[y + 1][x] - simg[y + 1][x - 1]) * (simg[y + 1][x] - simg[y + 1][x - 1]) + (simg[y + 1][x] - simg[y][x])*(simg[y + 1][x] - simg[y][x])));
		if (i == imgL - 1 || j == imgW - 1 || i == 0 || j == 0)
			val = 0;
		des[idx] = val;
	}
}
template<typename T>
void descentDir_SHARED_template(T* dimg, T* ddes, cuint imgL, cuint imgW)
{
	dim3 blk(32, 32);
	dim3 gid((imgL + 29) / 30, (imgW + 29) / 30);
	descentDir_SHARED_KER<T> << <gid, blk >> >(dimg, ddes, imgL, imgW);
}
void descentDir_SHARED(float* img, float* des, cuint imgL, cuint imgW)
{
	descentDir_SHARED_template<float>(img, des, imgL, imgW);
}
void descentDir_SHARED(double* img, double* des, cuint imgL, cuint imgW)
{
	descentDir_SHARED_template<double>(img, des, imgL, imgW);
}


//��άͼ���� 8x8x8,�����õ�ֻ��6x6x6����ı߱��������;
template<typename T>
__global__ void descentDir_SHARED_KER(T* img, T* des,
	cuint imgL, cuint imgW, cuint imgH)
{
	int i = threadIdx.x + blockIdx.x * (blockDim.x - 2);
	int j = threadIdx.y + blockIdx.y * (blockDim.y - 2);
	int k = threadIdx.z + blockIdx.z * (blockDim.z - 2);
	if (i < imgL && j < imgW && k < imgH)
	{
		int idx = (k * imgW + j) * imgL + i;
		__shared__ T f[8][8][8];
		f[threadIdx.z][threadIdx.y][threadIdx.x] = img[idx];
		__syncthreads();
		if (threadIdx.x == 7 || threadIdx.x == 0 || threadIdx.y == 7 || threadIdx.y == 0 || threadIdx.z == 7 || threadIdx.z == 0)
			return;
		cuint z = threadIdx.z;
		cuint y = threadIdx.y;
		cuint x = threadIdx.x;

		const T inv1 = 1.0 / sqrt((f[z][y][x] - f[z][y][x + 1]) * (f[z][y][x] - f[z][y][x + 1]) + (f[z][y][x] - f[z][y + 1][x]) * (f[z][y][x] - f[z][y + 1][x]) + (f[z][y][x] - f[z + 1][y][x]) * (f[z][y][x] - f[z + 1][y][x]) + 1e-36);
		const T inv2 = 1.0 / sqrt((f[z][y][x - 1] - f[z][y][x]) * (f[z][y][x - 1] - f[z][y][x]) + (f[z][y][x - 1] - f[z][y + 1][x - 1]) * (f[z][y][x - 1] - f[z][y + 1][x - 1]) + (f[z][y][x - 1] - f[z + 1][y][x - 1]) * (f[z][y][x - 1] - f[z + 1][y][x - 1]) + 1e-36);
		const T inv3 = 1.0 / sqrt((f[z][y - 1][x] - f[z][y - 1][x + 1]) * (f[z][y - 1][x] - f[z][y - 1][x + 1]) + (f[z][y - 1][x] - f[z][y][x]) * (f[z][y - 1][x] - f[z][y][x]) + (f[z][y - 1][x] - f[z + 1][y - 1][x]) * (f[z][y - 1][x] - f[z + 1][y - 1][x]) + 1e-36);
		const T inv4 = 1.0 / sqrt((f[z - 1][y][x] - f[z - 1][y][x + 1]) * (f[z - 1][y][x] - f[z - 1][y][x + 1]) + (f[z - 1][y][x] - f[z - 1][y + 1][x]) * (f[z - 1][y][x] - f[z - 1][y + 1][x]) + (f[z - 1][y][x] - f[z][y][x]) * (f[z - 1][y][x] - f[z][y][x]) + 1e-36);

		T val = ((f[z][y][x] - f[z][y][x + 1]) + (f[z][y][x] - f[z][y + 1][x]) + (f[z][y][x] - f[z + 1][y][x])) * inv1 - (f[z][y][x - 1] - f[z][y][x]) * inv2
			- (f[z][y - 1][x] - f[z][y][x]) * inv3 - (f[z - 1][y][x] - f[z][y][x]) * inv4;

		if (i == imgL - 1 || j == imgW - 1 || k == imgH - 1 || i == 0 || j == 0 || k == 0)
			val = 0;
		des[idx] = val;
	}
}
template<typename T>
void descentDir_SHARED_template(T* dimg, T* ddes, cuint imgL, cuint imgW, cuint imgH)
{
	dim3 blk(8, 8, 8);
	dim3 gid((imgL + 5) / 6, (imgW + 5) / 6, (imgH + 5) / 6);
	descentDir_SHARED_KER<T> << <gid, blk >> >(dimg, ddes, imgL, imgW, imgH);
}
void descentDir_SHARED(float* img, float* des, cuint imgL, cuint imgW, cuint imgH)
{
	descentDir_SHARED_template<float>(img, des, imgL, imgW, imgH);
}
void descentDir_SHARED(double* img, double* des, cuint imgL, cuint imgW, cuint imgH)
{
	descentDir_SHARED_template<double>(img, des, imgL, imgW, imgH);
}


void Dx_GPU(float* f, float* d, const float& coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_Dx_ker<float> << <gid, blk >> >(f, d, coef, L, W);
}

void Dx_GPU(double* f, double* d, const double& coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_Dx_ker<double> << <gid, blk >> >(f, d, coef, L, W);
}


void D2x_GPU(float* f, float* d, const float& coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_D2x_ker<float> << <gid, blk >> >(f, d, coef, L, W);
}

void D2x_GPU(double* f, double* d, const double& coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_D2x_ker<double> << <gid, blk >> >(f, d, coef, L, W);
}


void Dx_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_Dx_ker<float> << <gid, blk >> >(f, d, coef, L, W, H);
}

void Dx_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_Dx_ker<double> << <gid, blk >> >(f, d, coef, L, W, H);
}


void D2x_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_D2x_ker<float> << <gid, blk >> >(f, d, coef, L, W, H);
}

void D2x_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_D2x_ker<double> << <gid, blk >> >(f, d, coef, L, W, H);
}

void Dxt_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_Dxt_ker<float> << <gid, blk >> >(f, d, coef, L, W);
}

void Dxt_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_Dxt_ker<double> << <gid, blk >> >(f, d, coef, L, W);
}

void D2xt_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_D2xt_ker<float> << <gid, blk >> >(f, d, coef, L, W);
}

void D2xt_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_D2xt_ker<double> << <gid, blk >> >(f, d, coef, L, W);
}

void Dxt_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_Dxt_ker<float> << <gid, blk >> >(f, d, coef, L, W, H);
}

void Dxt_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_Dxt_ker<double> << <gid, blk >> >(f, d, coef, L, W, H);
}

void D2xt_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_D2xt_ker<float> << <gid, blk >> >(f, d, coef, L, W, H);
}

void D2xt_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_D2xt_ker<double> << <gid, blk >> >(f, d, coef, L, W, H);
}

void Dy_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_Dy_ker<float> << <gid, blk >> >(f, d, coef, L, W);
}
void Dy_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_Dy_ker<double> << <gid, blk >> >(f, d, coef, L, W);
}

void D2y_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_D2y_ker<float> << <gid, blk >> >(f, d, coef, L, W);
}

void D2y_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_D2y_ker<double> << <gid, blk >> >(f, d, coef, L, W);
}


void Dy_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_Dy_ker<float> << <gid, blk >> >(f, d, coef, L, W, H);
}

void Dy_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_Dy_ker<double> << <gid, blk >> >(f, d, coef, L, W, H);
}


void D2y_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_D2y_ker<float> << <gid, blk >> >(f, d, coef, L, W, H);
}

void D2y_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_D2y_ker<double> << <gid, blk >> >(f, d, coef, L, W, H);
}

void Dyt_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_Dyt_ker<float> << <gid, blk >> >(f, d, coef, L, W);
}

void Dyt_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_Dyt_ker<double> << <gid, blk >> >(f, d, coef, L, W);
}


void Dyt_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_Dyt_ker<float> << <gid, blk >> >(f, d, coef, L, W, H);
}

void Dyt_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_Dyt_ker<double> << <gid, blk >> >(f, d, coef, L, W, H);
}

void D2yt_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_D2yt_ker<float> << <gid, blk >> >(f, d, coef, L, W);
}

void D2yt_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_D2yt_ker<double> << <gid, blk >> >(f, d, coef, L, W);
}


void D2yt_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_D2yt_ker<float> << <gid, blk >> >(f, d, coef, L, W, H);
}

void D2yt_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_D2yt_ker<double> << <gid, blk >> >(f, d, coef, L, W, H);
}


void Dz_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_Dz_ker<float> << <gid, blk >> >(f, d, coef, L, W, H);
}
void Dz_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_Dz_ker<double> << <gid, blk >> >(f, d, coef, L, W, H);
}
void D2z_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_D2z_ker<float> << <gid, blk >> >(f, d, coef, L, W, H);
}
void D2z_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_D2z_ker<double> << <gid, blk >> >(f, d, coef, L, W, H);
}

void Dzt_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_Dzt_ker<float> << <gid, blk >> >(f, d, coef, L, W, H);
}
void Dzt_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_Dzt_ker<double> << <gid, blk >> >(f, d, coef, L, W, H);
}
void D2zt_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_D2zt_ker<float> << <gid, blk >> >(f, d, coef, L, W, H);
}
void D2zt_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_D2zt_ker<double> << <gid, blk >> >(f, d, coef, L, W, H);
}


void Laplacian_GPU(float* f, float* l, const float coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_Laplacian_ker<float> << <gid, blk >> >(f, l, coef, L, W);
}

void Laplacian_GPU(double* f, double* l, const double coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_Laplacian_ker<double> << <gid, blk >> >(f, l, coef, L, W);
}

void Laplacian_GPU(float* f, float* l, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_Laplacian_ker<float> << <gid, blk >> >(f, l, coef, L, W, H);
}

void Laplacian_GPU(double* f, double* l, const double coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_Laplacian_ker<double> << <gid, blk >> >(f, l, coef, L, W, H);
}

void DiscreteGradientTrans_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_DiscreteGradientTrans_ker<float> << <gid, blk >> >(f, d, coef, L, W);
}

void DiscreteGradientTrans_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_DiscreteGradientTrans_ker<double> << <gid, blk >> >(f, d, coef, L, W);
}

void DiscreteGradientTrans_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_DiscreteGradientTrans_ker<float> << <gid, blk >> >(f, d, coef, L, W, H);
}

void DiscreteGradientTrans_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_DiscreteGradientTrans_ker<double> << <gid, blk >> >(f, d, coef, L, W, H);
}





void GradientOfTV_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_GradientOfTV_ker<float> << <gid, blk >> >(f, d, coef, L, W);
}

void GradientOfTV_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_GradientOfTV_ker<double> << <gid, blk >> >(f, d, coef, L, W);
}



void GradientOfTV_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_GradientOfTV_ker<float> << <gid, blk >> >(f, d, coef, L, W, H);
}

void GradientOfTV_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid)
{
	_GradientOfTV_ker<double> << <gid, blk >> >(f, d, coef, L, W, H);
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
__global__ void _invDiscreteGradientTransform_ker(T* f, T* d, T* r, const T omega, cuint L, cuint W, cuint H)
{
	cuint idx = threadIdx.x + blockDim.x * blockIdx.x;
	cuint idy = threadIdx.y + blockDim.y * blockIdx.y;
	cuint idz = threadIdx.z + blockDim.z * blockIdx.z;

	if (idx < L && idy < W && idz < H)
	{
		cuint curId = (idz * W + idy) * L + idx;
		cuint nexxd = (idz * W + idy) * L + ((idx + 1) % L);
		cuint lasxd = (idz * W + idy) * L + ((idx + L - 1) % L);
		cuint nexyd = (idz * W + ((idy + 1) % W)) * L + idx;
		cuint lasyd = (idz * W + ((idy + W - 1) % W)) * L + idx;
		cuint nexzd = (((idz + 1) % H) * W + idy) * L + idx;
		cuint laszd = (((idz + H - 1) % H) * W + idy) * L + idx;

		T fa(0), fb(0), fc(0), fd(0);
		if (d[curId] < omega)
		{
			fa = (3.0 * f[curId] + f[nexxd] + f[nexyd] + f[nexzd]) * 0.1666666666666667;
		}
		else
		{
			fa = f[curId] - omega *(3 * f[curId] - f[nexxd] - f[nexyd] - f[nexzd]) * 0.1666666666666667 / d[curId];
		}
		if (d[lasxd] < omega)
		{
			fb = (f[curId] + f[lasxd]) * 0.5;
		}
		else
		{
			fb = f[curId] - omega * (f[curId] - f[lasxd]) * 0.5 / d[lasxd];
		}
		if (d[lasyd] < omega)
		{
			fc = (f[curId] + f[lasyd]) * 0.5;
		}
		else
		{
			fc = f[curId] - omega * (f[curId] - f[lasyd]) * 0.5 / d[lasyd];
		}
		if (d[laszd] < omega)
		{
			fd = (f[curId] + f[laszd]) * 0.5;
		}
		else
		{
			fd = f[curId] - omega * (f[curId] - f[laszd]) * 0.5 / d[laszd];
		}
		r[curId] = 0.1666666666666667 * (3.0 * fa + fb + fc + fd);
	}
}




void invDiscreteGradientTransform_GPU(float* f, float* d, float* r, const float omega, cuint L, cuint W, const dim3& blk, const dim3& gid)
{
	_invDiscreteGradientTransform_ker<float> << <gid, blk >> >(f, d, r, omega, L, W);
}
void invDiscreteGradientTransform_GPU(double* f, double* d, double* r, const double omega, cuint L, cuint W, const dim3& blk, const dim3& gid)
{
	_invDiscreteGradientTransform_ker<double> << <gid, blk >> >(f, d, r, omega, L, W);
}
void invDiscreteGradientTransform_GPU(float* f, float* d, float* r, const float omega, cuint L, cuint W, cuint H, const dim3& blk, const dim3& gid)
{
	_invDiscreteGradientTransform_ker<float> << <gid, blk >> >(f, d, r, omega, L, W, H);
}
void invDiscreteGradientTransform_GPU(double* f, double* d, double* r, const double omega, cuint L, cuint W, cuint H, const dim3& blk, const dim3& gid)
{
	_invDiscreteGradientTransform_ker<double> << <gid, blk >> >(f, d, r, omega, L, W, H);
}


template<typename T>
__global__ void _updateImg_forDEMO2_ker(
	T* dimg, T* dcorImg, T* dwegImg, T* dmask,
	T lambda, cuint ImgLR, cuint ImgWR, cuint sliceNum)
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

void updateImg_forDEMO2_GPU(float* dimg, float* dcorImg, float* dwegImg, float* dmask, float lambda, cuint ImgLR, cuint ImgWR, cuint sliceNum, const dim3& blk, const dim3& gid, const cudaStream_t& streams)
{
	_updateImg_forDEMO2_ker<float> << <gid, blk, 0, streams >> >(dimg, dcorImg, dwegImg, dmask, lambda, ImgLR, ImgWR, sliceNum);
}


template<typename T>
void minmaxValue_template(thrust::device_vector<T>& v, T& minV, T& maxV)
{
	auto ext = thrust::minmax_element(v.begin(), v.end());
	minV = *(ext.first);
	maxV = *(ext.second);
}


void minmaxValue(thrust::device_vector<float>& vec, float& minV, float& maxV)
{
	minmaxValue_template<float>(vec, minV, maxV);
}
void minmaxValue(thrust::device_vector<double>& vec, double& minV, double& maxV)
{
	minmaxValue_template<double>(vec, minV, maxV);
}


template<typename T>
void minmaxValue_template(thrust::device_ptr<T>& vec, cuint length, T& minV, T& maxV)
{
	auto ext = thrust::minmax_element(vec, vec + length);
	minV = *(ext.first);
	maxV = *(ext.second);
}
void minmaxValue(thrust::device_ptr<float>& vec, cuint length, float& minV, float& maxV)
{
	minmaxValue_template<float>(vec, length, minV, maxV);
}

void minmaxValue(thrust::device_ptr<double>& vec, cuint length, double& minV, double& maxV)
{
	minmaxValue_template<double>(vec, length, minV, maxV);
}



template<typename T>
__global__ void _gen01Mask_ker(T* dmsk, const T ratio, cuint L, cuint W)
{
	cuint idx = threadIdx.x + blockDim.x * blockIdx.x;
	cuint idy = threadIdx.y + blockDim.y * blockIdx.y;
	if (idx < L && idy < W)
	{
		T hL = static_cast<T>(L) * 0.5f;
		T hW = static_cast<T>(W) * 0.5f;

		if (
			((static_cast<T>(idx) -hL + 0.5f) / hL)*((static_cast<T>(idx) -hL + 0.5f) / hL) +
			((static_cast<T>(idy) -hW + 0.5f) / hW)*((static_cast<T>(idy) -hW + 0.5f) / hW) < (ratio * ratio)
			)
		{
			dmsk[idy * L + idx] = 1;
		}
		else
		{
			dmsk[idy * L + idx] = 0;
		}
	}
}

void gen01Mask(thrust::device_vector<float>& dmsk, const float& ratio, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	float* pmsk = thrust::raw_pointer_cast(dmsk.data());
	_gen01Mask_ker<float> << <gid, blk >> >(pmsk, ratio, L, W);
}
void gen01Mask(thrust::device_vector<double>& dmsk, const double& ratio, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	double* pmsk = thrust::raw_pointer_cast(dmsk.data());
	_gen01Mask_ker<double> << <gid, blk >> >(pmsk, ratio, L, W);
}

void gen01Mask(float* dmsk, const float& ratio, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_gen01Mask_ker<float> << <gid, blk >> >(dmsk, ratio, L, W);
}
void gen01Mask(double* dmsk, const double& ratio, cuint& L, cuint& W, const dim3& blk, const dim3& gid)
{
	_gen01Mask_ker<double> << <gid, blk >> >(dmsk, ratio, L, W);
}



