/*
 * utilities.cuh
 * The utilities functions that used in GPU and CPU
 *  Created on: Oct 2, 2015
 *  Author: Rui Liu
 */
#include "utilities.hpp"
#ifndef UTILITIES_CUH_
#define UTILITIES_CUH_
#include <thrust/device_vector.h>


/// \brief SIDDON line integral function in 3D
inline __device__ float calSiddonOneRayKer(
	const float startX, const float startY, const float startZ,
	const float endX, const float endY, const float endZ,
	const float __MINOL__, const float __MINOW__, const float __MINOH__,
	const float __OBJSTPL__, const float __OBJSTPW__, const float __OBJSTPH__,
	const unsigned int __OBJLR__, const unsigned int __OBJWR__, const unsigned int __OBJHR__,
	cudaTextureObject_t volTex, float* totWeight)
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


	const float alphaMIN = fmaxf(alphaxmin, fmaxf(alphaymin, alphazmin));
	const float alphaMAX = fminf(alphaxmax, fminf(alphaymax, alphazmax));
	dev_minmaxIdxFun(startX, endX, __MINOL__, __OBJSTPL__, alphaMIN, alphaMAX, alphaxmin, alphaxmax, __OBJLR__ + 1, &imin, &imax);
	dev_minmaxIdxFun(startY, endY, __MINOW__, __OBJSTPW__, alphaMIN, alphaMAX, alphaymin, alphaymax, __OBJWR__ + 1, &jmin, &jmax);
	dev_minmaxIdxFun(startZ, endZ, __MINOH__, __OBJSTPH__, alphaMIN, alphaMAX, alphazmin, alphazmax, __OBJHR__ + 1, &kmin, &kmax);

	float alphaX = (startX < endX) ? dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, imin) : dev_alpha_IFun(__MINOL__, __OBJSTPL__, startX, endX, imax);
	float alphaY = (startY < endY) ? dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, jmin) : dev_alpha_IFun(__MINOW__, __OBJSTPW__, startY, endY, jmax);
	float alphaZ = (startZ < endZ) ? dev_alpha_IFun(__MINOH__, __OBJSTPH__, startZ, endZ, kmin) : dev_alpha_IFun(__MINOH__, __OBJSTPH__, startZ, endZ, kmax);

	int Np = static_cast<int>(abs(static_cast<float>(imax - imin)) + abs(static_cast<float>(jmax - jmin)) + abs(static_cast<float>(kmax - kmin)) + 3.0f);
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

	const int iuu = (startX < endX) ? 1 : -1;
	const int juu = (startY < endY) ? 1 : -1;
	const int kuu = (startZ < endZ) ? 1 : -1;

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
				d12 = d12 + weight * tex3D<float>(volTex, k, i, j);
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
				d12 = d12 + weight * tex3D<float>(volTex, k, i, j);
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
				d12 = d12 + weight * tex3D<float>(volTex, k, i, j);
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


enum ForwardDDMethod{PROJ_BRANCHLESS=0,PROJ_PSEUDODISTANCE=2};
enum BackwardDDMethod{BACK_BRANCHLESS=0,BACK_PSEUDODISTANCE=2,BACK_ZLINEBRANCHLESS=3};

template<typename T>
T* getRawPtr(thrust::device_vector<T>& vec)
{
	return thrust::raw_pointer_cast(&vec[0]);
}

#endif /* UTILITIES_CUH_ */

