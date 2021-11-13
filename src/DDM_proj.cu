/*
 * COPYRIGHT NOTICE
 * COPYRIGHT (c) 2015, Wake Forest and UMass Lowell
 * All rights reserved
 *
 * @file DDM_proj.cu
 * @brief The GPU based DD projection in brute force method
 *
 * @version 1.0
 * @author Rui Liu
 * @date May. 1, 2015
 *
 */


#include "cuda_runtime.h"
#include "DDM_proj.h"


#ifndef PI
#define PI 3.14159265358979323846264
#endif

template<typename T>
__device__ inline T intersectLength_device(const T& fixedmin, const T& fixedmax, const T& varimin, const T& varimax)
{
	const T left = (fixedmin > varimin) ? fixedmin : varimin;
	const T right = (fixedmax < varimax) ? fixedmax : varimax;
	return abs(right - left) * static_cast<double>(right > left);
}



template<typename T>
__global__ void DDM_ED_proj_ker(T* proj, const T* img, const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T detSize, const T detCntIdx,
	const int XN, const int YN, const int DN, const int PN, const T dd, const T dx, const T dy,
	const T* angs)
{
	const int detIdx = threadIdx.x + blockIdx.x * blockDim.x;
	const int angIdx = threadIdx.y + blockIdx.y * blockDim.y;
	if (detIdx < DN && angIdx < PN)
	{
		T curang = angs[angIdx];
		T minP = cos(curang);
		T maxP = sin(curang);
		T cursourX = S2O * minP;
		T cursourY = S2O * maxP;
		T summ = 0;

		T curDetXLeft = -O2D * minP + (detIdx - detCntIdx - 0.5) * dd * maxP; //当前det左边X坐标
		T curDetYLeft = -O2D * maxP - (detIdx - detCntIdx - 0.5) * dd * minP; //当前det左边Y坐标
		T curDetXRight = -O2D * minP + (detIdx - detCntIdx + 0.5) * dd * maxP; //当前det右边X坐标
		T curDetYRight = -O2D * maxP - (detIdx - detCntIdx + 0.5) * dd * minP; //当前det右边Y坐标

		T dirX = -O2D * minP + (detIdx - detCntIdx) * dd * maxP - cursourX;
		T dirY = -O2D * maxP - (detIdx - detCntIdx) * dd * minP - cursourY;
		T obj = hypot(dirX, dirY);
		dirX /= obj;
		dirY /= obj;
		T detPosLeft, detPosRight;
		T temp;
		int ii, jj;

		int minIdx, maxIdx;

		if ((curang > PI * 0.25 && curang <= PI * 0.75) || (curang >= PI * 1.25 && curang < PI * 1.75))
		{

			curang = abs(dirY); //当前光线和Y轴夹角余弦

			detPosLeft = (0 - cursourY) / (curDetYLeft - cursourY) * (curDetXLeft - cursourX) + cursourX; //det左边界X轴上的坐标;
			detPosRight = (0 - cursourY) / (curDetYRight - cursourY) * (curDetXRight - cursourX) + cursourX;//det右边界在x轴上的坐标;
			if (detPosLeft > detPosRight)
			{
				temp = detPosLeft;
				detPosLeft = detPosRight;
				detPosRight = temp;
			}

			for (jj = 0; jj < YN; jj++)
			{
				obj = (jj - YN / 2.0 + 0.5) * dy;
				minP = (obj - cursourY) / (curDetYLeft - cursourY) * (curDetXLeft - cursourX) + cursourX;
				maxP = (obj - cursourY) / (curDetYRight - cursourY) *  (curDetXRight - cursourX) + cursourX;
				if (minP > maxP)
				{
					temp = minP;
					minP = maxP;
					maxP = temp;

				}

				minIdx = floor(minP / dx + XN / 2.0);
				maxIdx = ceil(maxP / dx + XN / 2.0);

				if (maxIdx <= 0)
				{
					continue;
				}
				else if (minIdx > XN)
				{
					continue;
				}

				if (minIdx < 0)
				{
					minIdx = 0;
				}
				if (maxIdx > XN)
				{
					maxIdx = XN;
				}
				minP = (-cursourY) / (obj - cursourY) * ((minIdx - XN / 2.0) * dx - cursourX) + cursourX;
				for (ii = minIdx; ii < maxIdx; ++ii)
				{

					maxP = (-cursourY) / (obj - cursourY) * ((ii + 1 - XN / 2.0) * dx - cursourX) + cursourX;
					summ += img[jj * XN + ii] * intersectLength_device<double>(detPosLeft, detPosRight, minP, maxP);
					minP = maxP;
				}
			}
			proj[angIdx * DN + detIdx] = summ / (curang * (detPosRight - detPosLeft)) * dy;
			return;
		}
		else
		{

			curang = abs(dirX); //与Case1区别;
			detPosLeft = cursourX / (cursourX - curDetXLeft) * (curDetYLeft - cursourY) + cursourY; //det左边界X轴上的坐标;
			detPosRight = cursourX / (cursourX - curDetXRight) * (curDetYRight - cursourY) + cursourY;//det右边界在x轴上的坐标;

			if (detPosLeft > detPosRight)
			{
				temp = detPosLeft;
				detPosLeft = detPosRight;
				detPosRight = temp;
			}

			for (ii = 0; ii < XN; ++ii)
			{

				obj = (ii - YN / 2.0 + 0.5) * dy;
				minP = (obj - cursourX) / (curDetXLeft - cursourX) * (curDetYLeft - cursourY) + cursourY;
				maxP = (obj - cursourX) / (curDetXRight - cursourX) *  (curDetYRight - cursourY) + cursourY;
				if (minP > maxP)
				{
					temp = minP;
					minP = maxP;
					maxP = temp;
				}

				minIdx = floor(minP / dy + YN / 2.0);
				maxIdx = ceil(maxP / dy + YN / 2.0);

				if (maxIdx <= 0)
				{
					continue;
				}
				else if (minIdx > XN)
				{
					continue;
				}

				if (minIdx < 0)
				{
					minIdx = 0;
				}
				if (maxIdx > YN)
				{
					maxIdx = YN;
				}


				minP = (-cursourX) / (obj - cursourX) * ((minIdx - YN / 2.0) * dy - cursourY) + cursourY;
				for (jj = minIdx; jj < maxIdx; ++jj)
				{
					maxP = (-cursourX) / (obj - cursourX) * ((jj + 1 - YN / 2.0) * dy - cursourY) + cursourY;
					summ += img[jj * XN + ii] * intersectLength_device<double>(detPosLeft, detPosRight, minP, maxP);
					minP = maxP;
				}
			}

			proj[angIdx * DN + detIdx] = summ / (curang * (detPosRight - detPosLeft)) * dx;
		}

	}
}


template<typename T>
void DDM_ED_proj_GPU_template(T* proj, const T* img,
	const T S2O, const T O2D, const T objSizeX, const T objSizeY,
	const T detSize, const T detCntIdx,
	const int XN, const int YN, const int DN, const int PN, const T dd, const T dx, const T dy,
	const T* angs, const dim3 blk, const dim3 gid)
{
	DDM_ED_proj_ker<T> << <gid, blk >> >(proj, img, S2O, O2D, objSizeX, objSizeY,
		detSize, detCntIdx, XN, YN, DN, PN, dd, dx, dy, angs);
}

void DDM_ED_proj_GPU(double* proj, const double* img,
	const double S2O, const double O2D, const double objSizeX, const double objSizeY,
	const double detSize, const double detCntIdx,
	const int XN, const int YN, const int DN, const int PN, const double dd, const double dx, const double dy,
	const double* angs, const dim3 blk, const dim3 gid)
{
	DDM_ED_proj_GPU_template<double>(proj, img, S2O, O2D, objSizeX, objSizeY,
		detSize, detCntIdx, XN, YN, DN, PN, dd, dx, dy, angs, blk, gid);
}
void DDM_ED_proj_GPU(float* proj, const float* img,
	const float S2O, const float O2D, const float objSizeX, const float objSizeY,
	const float detSize, const float detCntIdx,
	const int XN, const int YN, const int DN, const int PN, const float dd, const float dx, const float dy,
	const float* angs, const dim3 blk, const dim3 gid)
{
	DDM_ED_proj_GPU_template<float>(proj, img, S2O, O2D, objSizeX, objSizeY,
		detSize, detCntIdx, XN, YN, DN, PN, dd, dx, dy, angs, blk, gid);
}


void DDM_ED_proj_GPU(thrust::device_vector<float>& proj, const thrust::device_vector<float>& img,
	const float S2O, const float O2D, const float objSizeX, const float objSizeY,
	const float detSize, const float detCntIdx,
	const int XN, const int YN, const int DN, const int PN, const float dd, const float dx, const float dy,
	const thrust::device_vector<float>& angs, const dim3 blk, const dim3 gid)
{
	DDM_ED_proj_GPU(thrust::raw_pointer_cast(&proj[0]), thrust::raw_pointer_cast(&img[0]),
		S2O, O2D, objSizeX, objSizeY, detSize, detCntIdx,
		XN, YN, DN, PN, dd, dx, dy, thrust::raw_pointer_cast(&angs[0]), blk, gid);
}
void DDM_ED_proj_GPU(thrust::device_vector<double>& proj, const thrust::device_vector<double>& img,
	const double S2O, const double O2D, const double objSizeX, const double objSizeY,
	const double detSize, const double detCntIdx,
	const int XN, const int YN, const int DN, const int PN, const double dd, const double dx, const double dy,
	const thrust::device_vector<double>& angs, const dim3 blk, const dim3 gid)
{
	DDM_ED_proj_GPU(thrust::raw_pointer_cast(&proj[0]), thrust::raw_pointer_cast(&img[0]),
		S2O, O2D, objSizeX, objSizeY, detSize, detCntIdx,
		XN, YN, DN, PN, dd, dx, dy, thrust::raw_pointer_cast(&angs[0]), blk, gid);
}




template<typename T>
__global__ void DDM3D_ED_proj_GPU_template(T* proj, const T* vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detSizeU, const T detSizeV,
	const T detCntIdxU, const T detCntIdxV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const T ddu, const T ddv, const T dx, const T dy, const T dz,
	const T* angs)
{
	const int detIdU = threadIdx.x + blockIdx.x * blockDim.x;
	const int detIdV = threadIdx.y + blockIdx.y * blockDim.y;
	const int angIdx = threadIdx.z + blockIdx.z * blockDim.z;
	if (detIdU < DNU && detIdV < DNV && angIdx < PN)
	{
		T curang = angs[angIdx];
		T cosT = cos(curang);
		T sinT = sin(curang);

		T cursourx = S2O * cosT;
		T cursoury = S2O * sinT;
		T cursourz = 0;
		T temp;
		T summ = 0;
		T initDetX = -O2D;


		T initDetY = (detIdU - detCntIdxU) * ddu;
		T initDetZ = (detIdV - detCntIdxV) * ddv;

		T initDetLY = (detIdU - detCntIdxU - 0.5) * ddu;
		//initDetLZ = (detIdV - detCntIdxV) * ddv;

		T initDetRY = (detIdU - detCntIdxU + 0.5) * ddu;
		//initDetRZ = (detIdV - detCntIdxV) * ddv;

		T initDetDY = (detIdU - detCntIdxU) * ddu;
		T initDetDZ = (detIdV - detCntIdxV - 0.5) * ddv;

		T initDetUY = (detIdU - detCntIdxU) * ddu;
		T initDetUZ = (detIdV - detCntIdxV + 0.5) * ddv;

		T curDetLX = initDetX * cosT - initDetLY * sinT;
		T curDetLY = initDetX * sinT + initDetLY * cosT;
		//curDetLZ = initDetLZ;

		T curDetRX = initDetX * cosT - initDetRY * sinT;
		T curDetRY = initDetX * sinT + initDetRY * cosT;
		//curDetRZ = initDetRZ;

		T curDetDX = initDetX * cosT - initDetDY * sinT;
		T curDetDY = initDetX * sinT + initDetDY * cosT;
		T curDetDZ = initDetDZ;

		T curDetUX = initDetX * cosT - initDetUY * sinT;
		T curDetUY = initDetX * sinT + initDetUY * cosT;
		T curDetUZ = initDetUZ;

		T curDetX = initDetX * cosT - initDetY * sinT;
		T curDetY = initDetX * sinT + initDetY * cosT;
		//curDetZ = initDetZ;

		T dirX = curDetX - cursourx;
		T dirY = curDetY - cursoury;
		//dirZ = curDetZ - cursourz;

		if ((curang > PI * 0.25 && curang <= PI * 0.75) || (curang >= PI * 1.25 && curang < PI * 1.75))
		{




			T cosAlpha = abs(dirY / sqrt(dirY * dirY + dirX * dirX));
			//cosGamma = abs(dirY / sqrt(dirY * dirY + dirZ * dirZ));
			T cosGamma = abs(sqrt((S2O + O2D)*(S2O + O2D) - initDetZ*initDetZ) / (S2O + O2D));

			T detPosLX = -cursoury * (curDetLX - cursourx) / (curDetLY - cursoury) + cursourx; //左边点在XOZ平面上的投影;
			T detPosRX = -cursoury * (curDetRX - cursourx) / (curDetRY - cursoury) + cursourx;
			T detPosDZ = -cursoury * (curDetDZ - cursourz) / (curDetDY - cursoury) + cursourz;
			T detPosUZ = -cursoury * (curDetUZ - cursourz) / (curDetUY - cursoury) + cursourz;

			T detprojLength = abs(detPosLX - detPosRX);
			T detprojHeight = abs(detPosUZ - detPosDZ);

			//假设左边的小;
			if (detPosLX > detPosRX)
			{
				temp = detPosLX;
				detPosLX = detPosRX;
				detPosRX = temp;
				//std::swap(detPosLX, detPosRX);
			}
			//假设下边的小;
			if (detPosDZ > detPosUZ)
			{
				temp = detPosDZ;
				detPosDZ = detPosUZ;
				detPosUZ = temp;
				//std::swap(detPosDZ, detPosUZ);
			}

			for (size_t jj = 0; jj < YN; jj++)
			{
				T objY = (jj - YN / 2.0 + 0.5) * dy;
				T temp = (objY - cursoury) / (curDetLY - cursoury);

				T minX = temp * (curDetLX - cursourx) + cursourx;
				T maxX = temp * (curDetRX - cursourx) + cursourx;
				T minZ = temp * (curDetDZ - cursourz) + cursourz;
				T maxZ = temp * (curDetUZ - cursourz) + cursourz;
				if (minX > maxX)
				{
					temp = minX;
					minX = maxX;
					maxX = temp;
					//std::swap(minX, maxX);
				}
				if (minZ > maxZ)
				{
					temp = minZ;
					minZ = maxZ;
					maxZ = temp;

					//std::swap(minZ, maxZ);
				}
				int minXIdx = floor(minX / dx + XN / 2.0) - 2;
				int maxXIdx = ceil(maxX / dx + XN / 2.0) + 2;
				int minZIdx = floor(minZ / dz + ZN / 2.0) - 2;
				int maxZIdx = ceil(maxZ / dz + ZN / 2.0) + 2;
				if (maxXIdx < 0){ continue; }
				if (minXIdx > XN){ continue; }
				if (maxZIdx < 0){ continue; }
				if (minZIdx > ZN){ continue; }
				if (minXIdx < 0){ minXIdx = 0; }
				if (maxXIdx > XN){ maxXIdx = XN; }
				if (minZIdx < 0){ minZIdx = 0; }
				if (maxZIdx > ZN){ maxZIdx = ZN; }


				for (size_t ii = minXIdx; ii < maxXIdx; ii++)
				{
					T curminx = (cursourx - (ii - XN / 2.0) * dx) * cursoury / (objY - cursoury) + cursourx;
					T curmaxx = (cursourx - ((ii + 1) - XN / 2.0) * dx) * cursoury / (objY - cursoury) + cursourx;
					T intersectL = intersectLength_device<double>(detPosLX, detPosRX, curminx, curmaxx);
					if (intersectL > 0)
					{
						for (size_t kk = minZIdx; kk < maxZIdx; kk++)
						{

							T curminz = (cursourz - (kk - ZN / 2.0) * dz) * cursoury / (objY - cursoury) + cursourz;
							T curmaxz = (cursourz - ((kk + 1) - ZN / 2.0) * dz) * cursoury / (objY - cursoury) + cursourz;

							T intersectH = intersectLength_device<double>(detPosDZ, detPosUZ, curminz, curmaxz);
							if (intersectH > 0)
							{
								summ += vol[(kk * YN + jj) * XN + ii] * (intersectL * intersectH) / (detprojLength * detprojHeight * cosAlpha * cosGamma) * dx;
							}
						}
					}
					else
					{
						continue;
					}

				}

			}
			proj[(angIdx * DNV + detIdV) * DNU + detIdU] = summ;
		}
		else
		{
			T cosAlpha = abs(dirX / sqrt(dirY * dirY + dirX * dirX));
			//cosGamma = abs(dirY / sqrt(dirY * dirY + dirZ * dirZ));
			T cosGamma = abs(sqrt((S2O + O2D)*(S2O + O2D) - initDetZ*initDetZ) / (S2O + O2D));

			T detPosLY = -cursourx * (curDetLY - cursoury) / (curDetLX - cursourx) + cursoury; //左边点在XOZ平面上的投影;
			T detPosRY = -cursourx * (curDetRY - cursoury) / (curDetRX - cursourx) + cursoury;
			T detPosDZ = -cursourx * (curDetDZ - cursourz) / (curDetDX - cursourx) + cursourz;
			T detPosUZ = -cursourx * (curDetUZ - cursourz) / (curDetUX - cursourx) + cursourz;

			T detprojLength = abs(detPosLY - detPosRY);
			T detprojHeight = abs(detPosUZ - detPosDZ);

			//假设左边的小;
			if (detPosLY > detPosRY)
			{
				temp = detPosLY;
				detPosLY = detPosRY;
				detPosRY = temp;
				//std::swap(detPosLY, detPosRY);
			}
			//假设下边的小;
			if (detPosDZ > detPosUZ)
			{
				temp = detPosDZ;
				detPosDZ = detPosUZ;
				detPosUZ = temp;
				//std::swap(detPosDZ, detPosUZ);
			}

			for (size_t ii = 0; ii < XN; ii++)
			{
				T objX = (ii - XN / 2.0 + 0.5) * dx;
				T temp = (objX - cursourx) / (curDetLX - cursourx);

				T minY = temp * (curDetLY - cursoury) + cursoury;
				T maxY = temp * (curDetRY - cursoury) + cursoury;
				T minZ = temp * (curDetDZ - cursourz) + cursourz;
				T maxZ = temp * (curDetUZ - cursourz) + cursourz;
				if (minY > maxY)
				{
					temp = minY;
					minY = maxY;
					maxY = temp;
					//std::swap(minY, maxY);
				}
				if (minZ > maxZ)
				{
					temp = minZ;
					minZ = maxZ;
					maxZ = temp;
					//std::swap(minZ, maxZ);
				}
				int minYIdx = floor(minY / dy + YN / 2.0) - 2;
				int maxYIdx = ceil(maxY / dy + YN / 2.0) + 2;
				int minZIdx = floor(minZ / dz + ZN / 2.0) - 2;
				int maxZIdx = ceil(maxZ / dz + ZN / 2.0) + 2;
				if (maxYIdx < 0){ continue; }
				if (minYIdx > XN){ continue; }
				if (maxZIdx < 0){ continue; }
				if (minZIdx > ZN){ continue; }
				if (minYIdx < 0){ minYIdx = 0; }
				if (maxYIdx > XN){ maxYIdx = YN; }
				if (minZIdx < 0){ minZIdx = 0; }
				if (maxZIdx > ZN){ maxZIdx = ZN; }


				for (size_t jj = minYIdx; jj < maxYIdx; jj++)
				{
					T curminy = (cursoury - (jj - YN / 2.0) * dy) * cursourx / (objX - cursourx) + cursoury;
					T curmaxy = (cursoury - ((jj + 1) - YN / 2.0) * dy) * cursourx / (objX - cursourx) + cursoury;
					T intersectL = intersectLength_device<double>(detPosLY, detPosRY, curminy, curmaxy);
					if (intersectL > 0)
					{
						for (size_t kk = minZIdx; kk < maxZIdx; kk++)
						{

							T curminz = (cursourz - (kk - ZN / 2.0) * dz) * cursourx / (objX - cursourx) + cursourz;
							T curmaxz = (cursourz - ((kk + 1) - ZN / 2.0) * dz) * cursourx / (objX - cursourx) + cursourz;

							T intersectH = intersectLength_device<double>(detPosDZ, detPosUZ, curminz, curmaxz);
							if (intersectH > 0)
							{
								summ += vol[(kk * YN + jj) * XN + ii] * (intersectL * intersectH) / (detprojLength * detprojHeight * cosAlpha * cosGamma) * dx;
							}
						}
					}
					else
					{
						continue;
					}

				}

			}
			proj[(angIdx * DNV + detIdV) * DNU + detIdU] = summ;
		}
	}
}


void DDM3D_ED_proj_GPU(double* proj, const double* vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeU, const double detSizeV,
	const double detCntIdxU, const double detCntIdxV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double ddu, const double ddv, const double dx, const double dy, const double dz,
	const double* angs, const dim3 blk, const dim3 gid)
{
	DDM3D_ED_proj_GPU_template<double> << <gid, blk >> >(proj, vol, S2O, O2D,
		objSizeX, objSizeY, objSizeZ, detSizeU, detSizeV,
		detCntIdxU, detCntIdxV, XN, YN, ZN, DNU, DNV, PN,
		ddu, ddv, dx, dy, dz, angs);
}

void DDM3D_ED_proj_GPU(float* proj, const float* vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeU, const float detSizeV,
	const float detCntIdxU, const float detCntIdxV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float ddu, const float ddv, const float dx, const float dy, const float dz,
	const float* angs, const dim3 blk, const dim3 gid)
{
	DDM3D_ED_proj_GPU_template<float> << <gid, blk >> >(proj, vol, S2O, O2D,
		objSizeX, objSizeY, objSizeZ, detSizeU, detSizeV,
		detCntIdxU, detCntIdxV, XN, YN, ZN, DNU, DNV, PN,
		ddu, ddv, dx, dy, dz, angs);
}
void DDM3D_ED_proj_GPU(thrust::device_vector<float>& proj, const thrust::device_vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeU, const float detSizeV,
	const float detCntIdxU, const float detCntIdxV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float ddu, const float ddv, const float dx, const float dy, const float dz,
	const thrust::device_vector<float>& angs, const dim3 blk, const dim3 gid)
{
	DDM3D_ED_proj_GPU(thrust::raw_pointer_cast(&proj[0]), thrust::raw_pointer_cast(&vol[0]),
		S2O, O2D, objSizeX, objSizeY, objSizeZ, detSizeU, detSizeV, detCntIdxU, detCntIdxV,
		XN, YN, ZN, DNU, DNV, PN, ddu, ddv, dx, dy, dz, thrust::raw_pointer_cast(&angs[0]),
		blk, gid);
}
void DDM3D_ED_proj_GPU(thrust::device_vector<double>& proj, const thrust::device_vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeU, const double detSizeV,
	const double detCntIdxU, const double detCntIdxV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double ddu, const double ddv, const double dx, const double dy, const double dz,
	const thrust::device_vector<double>& angs, const dim3 blk, const dim3 gid)
{
	DDM3D_ED_proj_GPU(thrust::raw_pointer_cast(&proj[0]), thrust::raw_pointer_cast(&vol[0]),
		S2O, O2D, objSizeX, objSizeY, objSizeZ, detSizeU, detSizeV, detCntIdxU, detCntIdxV,
		XN, YN, ZN, DNU, DNV, PN, ddu, ddv, dx, dy, dz, thrust::raw_pointer_cast(&angs[0]),
		blk, gid);
}



template<typename T>
__global__ void DDM3D_EA_helical_proj_GPU_template(T* proj, const T* vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detArc, const T detSizeV,
	const T detCntIdU, const T detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const T dbeta, const T ddv, const T dx, const T dy, const T dz,
	const T initZPos, const T pitch, const T* angs)
{
	const int detIdU = threadIdx.x + blockIdx.x * blockDim.x;
	const int detIdV = threadIdx.y + blockIdx.y * blockDim.y;
	const int angIdx = threadIdx.z + blockIdx.z * blockDim.z;
	if (detIdU < DNU && detIdV < DNV && angIdx < PN)
	{
		T curang = angs[angIdx];
		T cosT = cos(curang);
		T sinT = sin(curang);

		T cursourx = S2O * cosT;
		T cursoury = S2O * sinT;
		T cursourz = initZPos + angIdx * pitch; //ÕâµãÓëcone beam ²»Í¬
		T summ = 0;
		T beta = (detIdU - detCntIdU) * dbeta;

		T cosBeta = cos(beta);
		T sinBeta = sin(beta);
		T initDetX = (-O2D) * cosBeta - 0 * sinBeta + S2O;
		T initDetY = (-O2D) * sinBeta + 0 * cosBeta;
		T initDetZ = (detIdV - detCntIdV) * ddv;

		beta = (detIdU - detCntIdU - 0.5) * dbeta;
		cosBeta = cos(beta);
		sinBeta = sin(beta);
		T initDetLX = (-O2D) * cosBeta - 0 * sinBeta + S2O;
		T initDetLY = (-O2D) * sinBeta + 0 * cosBeta;


		beta = (detIdU - detCntIdU + 0.5) * dbeta;
		cosBeta = cos(beta);
		sinBeta = sin(beta);
		T initDetRX = (-O2D) * cosBeta - 0 * sinBeta + S2O;
		T initDetRY = (-O2D) * sinBeta + 0 * cosBeta;

		T initDetDZ = (detIdV - detCntIdV - 0.5) * ddv;
		T initDetUZ = (detIdV - detCntIdV + 0.5) * ddv;

		T curDetLX = initDetLX * cosT - initDetLY * sinT;
		T curDetLY = initDetLX * sinT + initDetLY * cosT;

		T curDetRX = initDetRX * cosT - initDetRY * sinT;
		T curDetRY = initDetRX * sinT + initDetRY * cosT;

		T curDetDX = initDetX * cosT - initDetY * sinT;
		T curDetDY = initDetX * sinT + initDetY * cosT;
		T curDetDZ = initDetDZ + initZPos + angIdx * pitch; //

		T curDetUX = curDetDX;
		T curDetUY = curDetDY;
		T curDetUZ = initDetUZ + initZPos + angIdx * pitch; //

		T curDetX = initDetX * cosT - initDetY * sinT;
		T curDetY = initDetX * sinT + initDetY * cosT;

		T dirX = curDetX - cursourx;
		T dirY = curDetY - cursoury;

		if ((curang > PI * 0.25 && curang <= PI * 0.75) || (curang >= PI * 1.25 && curang < PI * 1.75))
		{
			T cosAlpha = abs(dirY / sqrt(dirY * dirY + dirX * dirX));
			T cosGamma = abs(sqrt((S2O + O2D) * (S2O + O2D) - initDetZ * initDetZ) / (S2O + O2D));


			T detPosLX = -cursoury * (curDetLX - cursourx) / (curDetLY - cursoury) + cursourx; //×ó±ßµãÔÚXOZÆœÃæÉÏµÄÍ¶Ó°;
			T detPosRX = -cursoury * (curDetRX - cursourx) / (curDetRY - cursoury) + cursourx;
			T detPosDZ = -cursoury * (curDetDZ - cursourz) / (curDetDY - cursoury) + cursourz;
			T detPosUZ = -cursoury * (curDetUZ - cursourz) / (curDetUY - cursoury) + cursourz;

			T detprojLength = abs(detPosLX - detPosRX);
			T detprojHeight = abs(detPosUZ - detPosDZ);

			//ŒÙÉè×ó±ßµÄÐ¡;
			if (detPosLX > detPosRX)
			{
				T tt = detPosLX;
				detPosLX = detPosRX;
				detPosRX = tt;
				//std::swap(detPosLX, detPosRX);
			}
			//ŒÙÉèÏÂ±ßµÄÐ¡;
			if (detPosDZ > detPosUZ)
			{
				T tt = detPosDZ;
				detPosDZ = detPosUZ;
				detPosUZ = tt;
				//std::swap(detPosDZ, detPosUZ);
			}


			for (size_t jj = 0; jj < YN; jj++)
			{
				T objY = (jj - YN / 2.0 + 0.5) * dy;
				T temp = (objY - cursoury) / (curDetLY - cursoury);

				T minX = temp * (curDetLX - cursourx) + cursourx;
				T maxX = temp * (curDetRX - cursourx) + cursourx;
				T minZ = temp * (curDetDZ - cursourz) + cursourz;
				T maxZ = temp * (curDetUZ - cursourz) + cursourz;
				if (minX > maxX)
				{
					T tt = minX;
					minX = maxX;
					maxX = tt;
					//std::swap(minX, maxX);
				}
				if (minZ > maxZ)
				{
					T tt = minZ;
					minZ = maxZ;
					maxZ = tt;
					//std::swap(minZ, maxZ);
				}
				int minXIdx = floor(minX / dx + XN / 2.0) - 2;
				int maxXIdx = ceil(maxX / dx + XN / 2.0) + 2;
				int minZIdx = floor(minZ / dz + ZN / 2.0) - 2;
				int maxZIdx = ceil(maxZ / dz + ZN / 2.0) + 2;
				if (maxXIdx < 0){ continue; }
				if (minXIdx > XN){ continue; }
				if (maxZIdx < 0){ continue; }
				if (minZIdx > ZN){ continue; }
				if (minXIdx < 0){ minXIdx = 0; }
				if (maxXIdx > XN){ maxXIdx = XN; }
				if (minZIdx < 0){ minZIdx = 0; }
				if (maxZIdx > ZN){ maxZIdx = ZN; }


				for (size_t ii = minXIdx; ii < maxXIdx; ii++)
				{
					T curminx = (cursourx - (ii - XN / 2.0) * dx) * cursoury / (objY - cursoury) + cursourx;
					T curmaxx = (cursourx - ((ii + 1) - XN / 2.0) * dx) * cursoury / (objY - cursoury) + cursourx;
					T intersectL = intersectLength_device<T>(detPosLX, detPosRX, curminx, curmaxx);
					if (intersectL > 0)
					{
						for (size_t kk = minZIdx; kk < maxZIdx; kk++)
						{

							T curminz = (cursourz - (kk - ZN / 2.0) * dz) * cursoury / (objY - cursoury) + cursourz;
							T curmaxz = (cursourz - ((kk + 1) - ZN / 2.0) * dz) * cursoury / (objY - cursoury) + cursourz;

							T intersectH = intersectLength_device<T>(detPosDZ, detPosUZ, curminz, curmaxz);
							if (intersectH > 0)
							{
								summ += vol[(kk * YN + jj) * XN + ii] * (intersectL * intersectH) / (detprojLength * detprojHeight * cosAlpha * cosGamma) * dx;
							}
						}
					}
					else
					{
						continue;
					}
				}
			}
			proj[(angIdx * DNV + detIdV) * DNU + detIdU] = summ;
		}
		else
		{
			T cosAlpha = abs(dirX / sqrt(dirY * dirY + dirX * dirX));
			T cosGamma = abs(sqrt((S2O + O2D) * (S2O + O2D) - initDetZ * initDetZ) / (S2O + O2D));


			T detPosLY = -cursourx * (curDetLY - cursoury) / (curDetLX - cursourx) + cursoury; //×ó±ßµãÔÚXOZÆœÃæÉÏµÄÍ¶Ó°;
			T detPosRY = -cursourx * (curDetRY - cursoury) / (curDetRX - cursourx) + cursoury;
			T detPosDZ = -cursourx * (curDetDZ - cursourz) / (curDetDX - cursourx) + cursourz;
			T detPosUZ = -cursourx * (curDetUZ - cursourz) / (curDetUX - cursourx) + cursourz;

			T detprojLength = abs(detPosLY - detPosRY);
			T detprojHeight = abs(detPosUZ - detPosDZ);

			//ŒÙÉè×ó±ßµÄÐ¡;
			if (detPosLY > detPosRY)
			{
				T tt = detPosLY;
				detPosLY = detPosRY;
				detPosRY = tt;
				//std::swap(detPosLY, detPosRY);
			}
			//ŒÙÉèÏÂ±ßµÄÐ¡;
			if (detPosDZ > detPosUZ)
			{
				T tt = detPosDZ;
				detPosDZ = detPosUZ;
				detPosUZ = tt;
				//std::swap(detPosDZ, detPosUZ);
			}

			for (size_t ii = 0; ii < XN; ii++)
			{
				T objX = (ii - XN / 2.0 + 0.5) * dx;
				T temp = (objX - cursourx) / (curDetLX - cursourx);

				T minY = temp * (curDetLY - cursoury) + cursoury;
				T maxY = temp * (curDetRY - cursoury) + cursoury;
				T minZ = temp * (curDetDZ - cursourz) + cursourz;
				T maxZ = temp * (curDetUZ - cursourz) + cursourz;
				if (minY > maxY)
				{
					T tt = minY;
					minY = maxY;
					maxY = tt;
					//std::swap(minY, maxY);
				}
				if (minZ > maxZ)
				{
					T tt = minZ;
					minZ = maxZ;
					maxZ = tt;
					//std::swap(minZ, maxZ);
				}
				int minYIdx = floor(minY / dy + YN / 2.0) - 2;
				int maxYIdx = ceil(maxY / dy + YN / 2.0) + 2;
				int minZIdx = floor(minZ / dz + ZN / 2.0) - 2;
				int maxZIdx = ceil(maxZ / dz + ZN / 2.0) + 2;
				if (maxYIdx < 0){ continue; }
				if (minYIdx > XN){ continue; }
				if (maxZIdx < 0){ continue; }
				if (minZIdx > ZN){ continue; }
				if (minYIdx < 0){ minYIdx = 0; }
				if (maxYIdx > XN){ maxYIdx = YN; }
				if (minZIdx < 0){ minZIdx = 0; }
				if (maxZIdx > ZN){ maxZIdx = ZN; }


				for (size_t jj = minYIdx; jj < maxYIdx; jj++)
				{
					T curminy = (cursoury - (jj - YN / 2.0) * dy) * cursourx / (objX - cursourx) + cursoury;
					T curmaxy = (cursoury - ((jj + 1) - YN / 2.0) * dy) * cursourx / (objX - cursourx) + cursoury;
					T intersectL = intersectLength_device<T>(detPosLY, detPosRY, curminy, curmaxy);
					if (intersectL > 0)
					{
						for (size_t kk = minZIdx; kk < maxZIdx; kk++)
						{

							T curminz = (cursourz - (kk - ZN / 2.0) * dz) * cursourx / (objX - cursourx) + cursourz;
							T curmaxz = (cursourz - ((kk + 1) - ZN / 2.0) * dz) * cursourx / (objX - cursourx) + cursourz;

							T intersectH = intersectLength_device<T>(detPosDZ, detPosUZ, curminz, curmaxz);
							if (intersectH > 0)
							{
								summ += vol[(kk * YN + jj) * XN + ii] * (intersectL * intersectH) / (detprojLength * detprojHeight * cosAlpha * cosGamma) * dx;
							}
						}
					}
					else
					{
						continue;
					}

				}

			}
			proj[(angIdx * DNV + detIdV) * DNU + detIdU] = summ;
		}

	}
}


void DDM3D_EA_helical_proj_GPU(float* proj, const float*  vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float dbeta, const float ddv, const float dx, const float dy, const float dz,
	const float initZPos, const float pitch, const float* angs, const dim3 blk, const dim3 gid)
{
	DDM3D_EA_helical_proj_GPU_template<float> << <gid, blk >> >(proj, vol, S2O, O2D, objSizeX,
		objSizeY, objSizeZ, detArc, detSizeH, detCntIdU, detCntIdV,
		XN, YN, ZN, DNU, DNV, PN, dbeta, ddv, dx, dy, dz, initZPos, pitch, angs);
}

void DDM3D_EA_helical_proj_GPU(double* proj, const double*  vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double dbeta, const double ddv, const double dx, const double dy, const double dz,
	const double initZPos, const double pitch, const double* angs, const dim3 blk, const dim3 gid)
{
	DDM3D_EA_helical_proj_GPU_template<double> << <gid, blk >> >(proj, vol, S2O, O2D, objSizeX,
		objSizeY, objSizeZ, detArc, detSizeH, detCntIdU, detCntIdV,
		XN, YN, ZN, DNU, DNV, PN, dbeta, ddv, dx, dy, dz, initZPos, pitch, angs);
}






template<typename T>
__global__ void DDM3D_EA_helical_proj_GPU_template(T* proj, const T* vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detArc, const T detSizeV,
	const T detCntIdU, const T detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const T dbeta, const T ddv, const T dx, const T dy, const T dz,
	const T* zShifts, const T* angs)
{
	const int detIdU = threadIdx.x + blockIdx.x * blockDim.x;
	const int detIdV = threadIdx.y + blockIdx.y * blockDim.y;
	const int angIdx = threadIdx.z + blockIdx.z * blockDim.z;
	if (detIdU < DNU && detIdV < DNV && angIdx < PN)
	{
		T curang = angs[angIdx];
		T cosT = cos(curang);
		T sinT = sin(curang);

		T cursourx = S2O * cosT;
		T cursoury = S2O * sinT;
		T cursourz = zShifts[angIdx]; //ÕâµãÓëcone beam ²»Í¬
		T summ = 0;
		T beta = (detIdU - detCntIdU) * dbeta;

		T cosBeta = cos(beta);
		T sinBeta = sin(beta);
		T initDetX = (-O2D) * cosBeta - 0 * sinBeta + S2O;
		T initDetY = (-O2D) * sinBeta + 0 * cosBeta;
		T initDetZ = (detIdV - detCntIdV) * ddv;

		beta = (detIdU - detCntIdU - 0.5) * dbeta;
		cosBeta = cos(beta);
		sinBeta = sin(beta);
		T initDetLX = (-O2D) * cosBeta - 0 * sinBeta + S2O;
		T initDetLY = (-O2D) * sinBeta + 0 * cosBeta;


		beta = (detIdU - detCntIdU + 0.5) * dbeta;
		cosBeta = cos(beta);
		sinBeta = sin(beta);
		T initDetRX = (-O2D) * cosBeta - 0 * sinBeta + S2O;
		T initDetRY = (-O2D) * sinBeta + 0 * cosBeta;

		T initDetDZ = (detIdV - detCntIdV - 0.5) * ddv;
		T initDetUZ = (detIdV - detCntIdV + 0.5) * ddv;

		T curDetLX = initDetLX * cosT - initDetLY * sinT;
		T curDetLY = initDetLX * sinT + initDetLY * cosT;

		T curDetRX = initDetRX * cosT - initDetRY * sinT;
		T curDetRY = initDetRX * sinT + initDetRY * cosT;

		T curDetDX = initDetX * cosT - initDetY * sinT;
		T curDetDY = initDetX * sinT + initDetY * cosT;
		T curDetDZ = initDetDZ + zShifts[angIdx]; //

		T curDetUX = curDetDX;
		T curDetUY = curDetDY;
		T curDetUZ = initDetUZ + zShifts[angIdx]; //

		T curDetX = initDetX * cosT - initDetY * sinT;
		T curDetY = initDetX * sinT + initDetY * cosT;

		T dirX = curDetX - cursourx;
		T dirY = curDetY - cursoury;

		if ((curang > PI * 0.25 && curang <= PI * 0.75) || (curang >= PI * 1.25 && curang < PI * 1.75))
		{
			T cosAlpha = abs(dirY / sqrt(dirY * dirY + dirX * dirX));
			T cosGamma = abs(sqrt((S2O + O2D) * (S2O + O2D) - initDetZ * initDetZ) / (S2O + O2D));


			T detPosLX = -cursoury * (curDetLX - cursourx) / (curDetLY - cursoury) + cursourx; //×ó±ßµãÔÚXOZÆœÃæÉÏµÄÍ¶Ó°;
			T detPosRX = -cursoury * (curDetRX - cursourx) / (curDetRY - cursoury) + cursourx;
			T detPosDZ = -cursoury * (curDetDZ - cursourz) / (curDetDY - cursoury) + cursourz;
			T detPosUZ = -cursoury * (curDetUZ - cursourz) / (curDetUY - cursoury) + cursourz;

			T detprojLength = abs(detPosLX - detPosRX);
			T detprojHeight = abs(detPosUZ - detPosDZ);

			//ŒÙÉè×ó±ßµÄÐ¡;
			if (detPosLX > detPosRX)
			{
				T tt = detPosLX;
				detPosLX = detPosRX;
				detPosRX = tt;
				//std::swap(detPosLX, detPosRX);
			}
			//ŒÙÉèÏÂ±ßµÄÐ¡;
			if (detPosDZ > detPosUZ)
			{
				T tt = detPosDZ;
				detPosDZ = detPosUZ;
				detPosUZ = tt;
				//std::swap(detPosDZ, detPosUZ);
			}


			for (size_t jj = 0; jj < YN; jj++)
			{
				T objY = (jj - YN / 2.0 + 0.5) * dy;
				T temp = (objY - cursoury) / (curDetLY - cursoury);

				T minX = temp * (curDetLX - cursourx) + cursourx;
				T maxX = temp * (curDetRX - cursourx) + cursourx;
				T minZ = temp * (curDetDZ - cursourz) + cursourz;
				T maxZ = temp * (curDetUZ - cursourz) + cursourz;
				if (minX > maxX)
				{
					T tt = minX;
					minX = maxX;
					maxX = tt;
					//std::swap(minX, maxX);
				}
				if (minZ > maxZ)
				{
					T tt = minZ;
					minZ = maxZ;
					maxZ = tt;
					//std::swap(minZ, maxZ);
				}
				int minXIdx = floor(minX / dx + XN / 2.0) - 2;
				int maxXIdx = ceil(maxX / dx + XN / 2.0) + 2;
				int minZIdx = floor(minZ / dz + ZN / 2.0) - 2;
				int maxZIdx = ceil(maxZ / dz + ZN / 2.0) + 2;
				if (maxXIdx < 0){ continue; }
				if (minXIdx > XN){ continue; }
				if (maxZIdx < 0){ continue; }
				if (minZIdx > ZN){ continue; }
				if (minXIdx < 0){ minXIdx = 0; }
				if (maxXIdx > XN){ maxXIdx = XN; }
				if (minZIdx < 0){ minZIdx = 0; }
				if (maxZIdx > ZN){ maxZIdx = ZN; }


				for (size_t ii = minXIdx; ii < maxXIdx; ii++)
				{
					T curminx = (cursourx - (ii - XN / 2.0) * dx) * cursoury / (objY - cursoury) + cursourx;
					T curmaxx = (cursourx - ((ii + 1) - XN / 2.0) * dx) * cursoury / (objY - cursoury) + cursourx;
					T intersectL = intersectLength_device<T>(detPosLX, detPosRX, curminx, curmaxx);
					if (intersectL > 0)
					{
						for (size_t kk = minZIdx; kk < maxZIdx; kk++)
						{

							T curminz = (cursourz - (kk - ZN / 2.0) * dz) * cursoury / (objY - cursoury) + cursourz;
							T curmaxz = (cursourz - ((kk + 1) - ZN / 2.0) * dz) * cursoury / (objY - cursoury) + cursourz;

							T intersectH = intersectLength_device<T>(detPosDZ, detPosUZ, curminz, curmaxz);
							if (intersectH > 0)
							{
								summ += vol[(kk * YN + jj) * XN + ii] * (intersectL * intersectH) / (detprojLength * detprojHeight * cosAlpha * cosGamma) * dx;
							}
						}
					}
					else
					{
						continue;
					}
				}
			}
			proj[(angIdx * DNV + detIdV) * DNU + detIdU] = summ;
		}
		else
		{
			T cosAlpha = abs(dirX / sqrt(dirY * dirY + dirX * dirX));
			T cosGamma = abs(sqrt((S2O + O2D) * (S2O + O2D) - initDetZ * initDetZ) / (S2O + O2D));


			T detPosLY = -cursourx * (curDetLY - cursoury) / (curDetLX - cursourx) + cursoury; //×ó±ßµãÔÚXOZÆœÃæÉÏµÄÍ¶Ó°;
			T detPosRY = -cursourx * (curDetRY - cursoury) / (curDetRX - cursourx) + cursoury;
			T detPosDZ = -cursourx * (curDetDZ - cursourz) / (curDetDX - cursourx) + cursourz;
			T detPosUZ = -cursourx * (curDetUZ - cursourz) / (curDetUX - cursourx) + cursourz;

			T detprojLength = abs(detPosLY - detPosRY);
			T detprojHeight = abs(detPosUZ - detPosDZ);

			//ŒÙÉè×ó±ßµÄÐ¡;
			if (detPosLY > detPosRY)
			{
				T tt = detPosLY;
				detPosLY = detPosRY;
				detPosRY = tt;
				//std::swap(detPosLY, detPosRY);
			}
			//ŒÙÉèÏÂ±ßµÄÐ¡;
			if (detPosDZ > detPosUZ)
			{
				T tt = detPosDZ;
				detPosDZ = detPosUZ;
				detPosUZ = tt;
				//std::swap(detPosDZ, detPosUZ);
			}

			for (size_t ii = 0; ii < XN; ii++)
			{
				T objX = (ii - XN / 2.0 + 0.5) * dx;
				T temp = (objX - cursourx) / (curDetLX - cursourx);

				T minY = temp * (curDetLY - cursoury) + cursoury;
				T maxY = temp * (curDetRY - cursoury) + cursoury;
				T minZ = temp * (curDetDZ - cursourz) + cursourz;
				T maxZ = temp * (curDetUZ - cursourz) + cursourz;
				if (minY > maxY)
				{
					T tt = minY;
					minY = maxY;
					maxY = tt;
					//std::swap(minY, maxY);
				}
				if (minZ > maxZ)
				{
					T tt = minZ;
					minZ = maxZ;
					maxZ = tt;
					//std::swap(minZ, maxZ);
				}
				int minYIdx = floor(minY / dy + YN / 2.0) - 2;
				int maxYIdx = ceil(maxY / dy + YN / 2.0) + 2;
				int minZIdx = floor(minZ / dz + ZN / 2.0) - 2;
				int maxZIdx = ceil(maxZ / dz + ZN / 2.0) + 2;
				if (maxYIdx < 0){ continue; }
				if (minYIdx > XN){ continue; }
				if (maxZIdx < 0){ continue; }
				if (minZIdx > ZN){ continue; }
				if (minYIdx < 0){ minYIdx = 0; }
				if (maxYIdx > XN){ maxYIdx = YN; }
				if (minZIdx < 0){ minZIdx = 0; }
				if (maxZIdx > ZN){ maxZIdx = ZN; }


				for (size_t jj = minYIdx; jj < maxYIdx; jj++)
				{
					T curminy = (cursoury - (jj - YN / 2.0) * dy) * cursourx / (objX - cursourx) + cursoury;
					T curmaxy = (cursoury - ((jj + 1) - YN / 2.0) * dy) * cursourx / (objX - cursourx) + cursoury;
					T intersectL = intersectLength_device<T>(detPosLY, detPosRY, curminy, curmaxy);
					if (intersectL > 0)
					{
						for (size_t kk = minZIdx; kk < maxZIdx; kk++)
						{

							T curminz = (cursourz - (kk - ZN / 2.0) * dz) * cursourx / (objX - cursourx) + cursourz;
							T curmaxz = (cursourz - ((kk + 1) - ZN / 2.0) * dz) * cursourx / (objX - cursourx) + cursourz;

							T intersectH = intersectLength_device<T>(detPosDZ, detPosUZ, curminz, curmaxz);
							if (intersectH > 0)
							{
								summ += vol[(kk * YN + jj) * XN + ii] * (intersectL * intersectH) / (detprojLength * detprojHeight * cosAlpha * cosGamma) * dx;
							}
						}
					}
					else
					{
						continue;
					}

				}

			}
			proj[(angIdx * DNV + detIdV) * DNU + detIdU] = summ;
		}

	}
}

void DDM3D_EA_helical_proj_GPU(float* proj, const float*  vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float dbeta, const float ddv, const float dx, const float dy, const float dz,
	const float* zShifts, const float* angs, const dim3 blk, const dim3 gid)
{
	DDM3D_EA_helical_proj_GPU_template<float> << <gid, blk >> >(proj, vol,
		S2O, O2D, objSizeX, objSizeY, objSizeZ, detArc, detSizeH,
		detCntIdU, detCntIdV, XN, YN, ZN, DNU, DNV, PN,
		dbeta, ddv, dx, dy, dz, zShifts, angs);
}
void DDM3D_EA_helical_proj_GPU(double* proj, const double*  vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double dbeta, const double ddv, const double dx, const double dy, const double dz,
	const double* zShifts, const double* angs, const dim3 blk, const dim3 gid)
{
	DDM3D_EA_helical_proj_GPU_template<double> << <gid, blk >> >(proj, vol,
		S2O, O2D, objSizeX, objSizeY, objSizeZ, detArc, detSizeH,
		detCntIdU, detCntIdV, XN, YN, ZN, DNU, DNV, PN,
		dbeta, ddv, dx, dy, dz, zShifts, angs);
}
void DDM3D_EA_helical_proj_GPU(thrust::device_vector<float>& proj, const thrust::device_vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float dbeta, const float ddv, const float dx, const float dy, const float dz,
	const float initZPos, const float pitch, const thrust::device_vector<float>& angs, const dim3 blk, const dim3 gid)
{
	DDM3D_EA_helical_proj_GPU_template<float> << <gid, blk >> >(thrust::raw_pointer_cast(&proj[0]), thrust::raw_pointer_cast(&vol[0]),
		S2O, O2D, objSizeX, objSizeY, objSizeZ, detArc, detSizeH, detCntIdU, detCntIdV,
		XN, YN, ZN, DNU, DNV, PN, dbeta, ddv, dx, dy, dz, initZPos, pitch, thrust::raw_pointer_cast(&angs[0]));
}
void DDM3D_EA_helical_proj_GPU(thrust::device_vector<double>& proj, const thrust::device_vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double dbeta, const double ddv, const double dx, const double dy, const double dz,
	const double initZPos, const double pitch, const thrust::device_vector<double>& angs, const dim3 blk, const dim3 gid)
{
	DDM3D_EA_helical_proj_GPU_template<double> << <gid, blk >> >(thrust::raw_pointer_cast(&proj[0]), thrust::raw_pointer_cast(&vol[0]),
		S2O, O2D, objSizeX, objSizeY, objSizeZ, detArc, detSizeH, detCntIdU, detCntIdV,
		XN, YN, ZN, DNU, DNV, PN, dbeta, ddv, dx, dy, dz, initZPos, pitch, thrust::raw_pointer_cast(&angs[0]));
}


void DDM3D_EA_helical_proj_GPU(thrust::device_vector<float>& proj, const thrust::device_vector<float>&  vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float dbeta, const float ddv, const float dx, const float dy, const float dz,
	const thrust::device_vector<float>& zShifts, const thrust::device_vector<float>& angs, const dim3 blk, const dim3 gid)
{
	DDM3D_EA_helical_proj_GPU_template<float> << <gid, blk >> >(thrust::raw_pointer_cast(&proj[0]), thrust::raw_pointer_cast(&vol[0]),
		S2O, O2D, objSizeX, objSizeY, objSizeZ, detArc, detSizeH,
		detCntIdU, detCntIdV, XN, YN, ZN, DNU, DNV, PN,
		dbeta, ddv, dx, dy, dz, thrust::raw_pointer_cast(&zShifts[0]), thrust::raw_pointer_cast(&angs[0]));
}
void DDM3D_EA_helical_proj_GPU(thrust::device_vector<double>& proj, const thrust::device_vector<double>&  vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double dbeta, const double ddv, const double dx, const double dy, const double dz,
	const thrust::device_vector<double>& zShifts, const thrust::device_vector<double>& angs, const dim3 blk, const dim3 gid)
{
	DDM3D_EA_helical_proj_GPU_template<double> << <gid, blk >> >(thrust::raw_pointer_cast(&proj[0]), thrust::raw_pointer_cast(&vol[0]),
		S2O, O2D, objSizeX, objSizeY, objSizeZ, detArc, detSizeH,
		detCntIdU, detCntIdV, XN, YN, ZN, DNU, DNV, PN,
		dbeta, ddv, dx, dy, dz, thrust::raw_pointer_cast(&zShifts[0]), thrust::raw_pointer_cast(&angs[0]));
}

