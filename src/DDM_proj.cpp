/*
 * COPYRIGHT NOTICE
 * COPYRIGHT (c) 2015, Wake Forest and UMass Lowell
 * All rights reserved
 *
 * @file DDM_proj.cpp
 * @brief The CPU based DD projection in conventional method
 *
 * @version 1.0
 * @author Rui Liu
 * @date May. 1, 2015
 *
 */


#include "DDM_proj.h"
#include "utilities.hpp"


template<typename T>
inline T intersectLength(const T& fixedmin, const T& fixedmax, const T& varimin, const T& varimax)
{
	const T left = (fixedmin > varimin) ? fixedmin : varimin;
	const T right = (fixedmax < varimax) ? fixedmax : varimax;
	return abs(right - left) * static_cast<double>(right > left);
}


template<typename T>
void findMinMax_X(T& minX, T& maxX, const T x1, const T x2, const T x3, const T x4)
{
	minX = x1;
	if (x2 < minX)
	{
		minX = x2;
	}
	if (x3 < minX)
	{
		minX = x3;
	}
	if (x4 < minX)
	{
		minX = x4;
	}

	maxX = x1;
	if (x2 > maxX)
	{
		maxX = x2;
	}
	if (x3 > maxX)
	{
		maxX = x3;
	}
	if (x4 > maxX)
	{
		maxX = x4;
	}
}

//The ray propagate along Y axis
template<typename T>
void projCase4(std::vector<T>& proj, const std::vector<T>& img, const T cursourX, const T cursourY, const T S2O, const T O2D, const T objSizeX, const T objSizeY, const T detSize, const T detCntIdx, const int XN, const int YN, const int DN, const int PN,
	const T dx, const T dy, const T dd, const T curAng, const T cosT, const T sinT, const int angIdx)
{
	T summ = 0;
	T initX, initY, initYLeft, initYRight, curDetXLeft, curDetXRight, curDetYLeft, curDetYRight;
	T curDetX, curDetY, dirX, dirY, legth, cosAng, detPosLeft, detPosRight;
	T detprojLength;
	T objX;

	T minY, maxY;
	int minYIdx, maxYIdx, detIdx, ii, jj;
	T curminy, curmaxy;
	initX = -O2D;

	initYLeft = -(0 - detCntIdx - 0.5) * dd; //初始det左边Y坐标
	curDetXLeft = initX * cosT - initYLeft * sinT; //当前det左边X坐标
	curDetYLeft = initX * sinT + initYLeft * cosT; //当前det左边Y坐标

	for (detIdx = 0; detIdx != DN; ++detIdx)
	{
		summ = 0;

		initYRight = -(detIdx - detCntIdx - 0.5 + 1) * dd; //初始det右边Y坐标;
		initY = -(detIdx - detCntIdx) * dd; //初始det中心Y坐标;
		curDetXRight = initX * cosT - initYRight * sinT; //当前det右边X坐标
		curDetYRight = initX * sinT + initYRight * cosT; //当前det右边Y坐标
		curDetX = initX * cosT - initY * sinT; //当前det中心X坐标
		curDetY = initX * sinT + initY * cosT; //当前det中心Y坐标
		dirX = curDetX - cursourX;
		dirY = curDetY - cursourY;
		legth = hypot(dirX, dirY);
		dirX /= legth;
		dirY /= legth;
		cosAng = abs(dirX); //与Case1区别;

		//这里是从X坐标算Y坐标;

		detPosLeft = (0 - cursourX) / (curDetXLeft - cursourX) * (curDetYLeft - cursourY) + cursourY; //det左边界X轴上的坐标;
		detPosRight = (0 - cursourX) / (curDetXRight - cursourX) * (curDetYRight - cursourY) + cursourY;//det右边界在x轴上的坐标;
		detprojLength = abs(detPosRight - detPosLeft);

		//沿X方向扫描;
		for (ii = 0; ii < XN; ++ii)
		{

			objX = (ii - YN / 2.0 + 0.5) * dy;
			minY = (objX - cursourX) / (curDetXLeft - cursourX) * (curDetYLeft - cursourY) + cursourY;
			maxY = (objX - cursourX) / (curDetXRight - cursourX) *  (curDetYRight - cursourY) + cursourY;
			if (minY > maxY)
			{
				std::swap(minY, maxY);
			}

			minYIdx = floor(minY / dy + YN / 2.0);
			maxYIdx = ceil(maxY / dy + YN / 2.0);

			if (maxYIdx <= 0)
			{
				continue;
			}
			else if (minYIdx > XN)
			{
				continue;
			}

			if (minYIdx < 0)
			{
				minYIdx = 0;
			}
			if (maxYIdx > YN)
			{
				maxYIdx = YN;
			}


			curminy = (-cursourX) / (objX - cursourX) * ((0 - YN / 2.0) * dy - cursourY) + cursourY;
			for (jj = minYIdx; jj < maxYIdx; ++jj)
			{
				curmaxy = (-cursourX) / (objX - cursourX) * ((jj + 1 - YN / 2.0) * dy - cursourY) + cursourY;
				if (detPosLeft > detPosRight)
				{
					std::swap(detPosLeft, detPosRight);
				}
				summ += img[jj * XN + ii] * intersectLength<double>(detPosLeft, detPosRight, curminy, curmaxy) * dx;
				curminy = curmaxy;
			}
		}

		proj[angIdx * DN + detIdx] = summ / (cosAng * detprojLength);
		initYLeft = initYRight;
		curDetXLeft = curDetXRight;
		curDetYLeft = curDetYRight;
	}
}

//沿X轴传播;
template<typename T>
void projCase1(std::vector<T>& proj, const std::vector<T>& img, const T cursourX, const T cursourY, const T S2O, const T O2D, const T objSizeX, const T objSizeY, const T detSize, const T detCntIdx,
	const int XN, const int YN, const int DN, const int PN, const T dx, const T dy, const T dd, const T curAng, const T cosT, const T sinT, const int angIdx)
{
	int detIdx = 0;
	int ii = 0, jj = 0;
	T initX, initYLeft, initYRight;
	T curDetXLeft, curDetYLeft;
	T curDetXRight, curDetYRight;
	T minX, maxX;
	int minXIdx, maxXIdx;
	T detPosLeft;
	T detPosRight;
	T initY;
	T curDetX, curDetY;
	T dirX, dirY;
	T legth;
	T cosAng;
	T detprojLength = 0;
	T objY;
	T summ = 0;
	T curminx, curmaxx;

	initX = -O2D;
	initYLeft = -(0 - detCntIdx - 0.5) * dd; //初始det左边Y坐标
	curDetXLeft = initX * cosT - initYLeft * sinT; //当前det左边X坐标
	curDetYLeft = initX * sinT + initYLeft * cosT; //当前det左边Y坐标
	for (detIdx = 0; detIdx != DN; ++detIdx)
	{
		summ = 0;
		initYRight = initYLeft - dd;// -(detIdx - detCntIdx - 0.5 + 1) * dd; //初始det右边Y坐标;
		initY = -(detIdx - detCntIdx) * dd; //初始det中心Y坐标;
		curDetXRight = initX * cosT - initYRight * sinT; //当前det右边X坐标
		curDetYRight = initX * sinT + initYRight * cosT; //当前det右边Y坐标
		curDetX = initX * cosT - initY * sinT; //当前det中心X坐标
		curDetY = initX * sinT + initY * cosT; //当前det中心Y坐标
		dirX = curDetX - cursourX;
		dirY = curDetY - cursourY;
		legth = hypot(dirX, dirY);
		dirX /= legth;
		dirY /= legth;
		cosAng = abs(dirY); //当前光线和Y轴夹角余弦

		detPosLeft = (0 - cursourY) / (curDetYLeft - cursourY) * (curDetXLeft - cursourX) + cursourX; //det左边界X轴上的坐标;
		detPosRight = (0 - cursourY) / (curDetYRight - cursourY) * (curDetXRight - cursourX) + cursourX;//det右边界在x轴上的坐标;
		detprojLength = abs(detPosRight - detPosLeft);

		for (jj = 0; jj < YN; jj++)
		{
			objY = (jj - YN / 2.0 + 0.5) * dy;
			minX = (objY - cursourY) / (curDetYLeft - cursourY) * (curDetXLeft - cursourX) + cursourX;
			maxX = (objY - cursourY) / (curDetYRight - cursourY) *  (curDetXRight - cursourX) + cursourX;
			if (minX > maxX)
			{
				std::swap(minX, maxX);
			}

			minXIdx = floor(minX / dx + XN / 2.0);
			maxXIdx = ceil(maxX / dx + XN / 2.0);

			if (maxXIdx <= 0)
			{
				continue;
			}
			else if (minXIdx > XN)
			{
				continue;
			}

			if (minXIdx < 0)
			{
				minXIdx = 0;
			}
			if (maxXIdx > XN)
			{
				maxXIdx = XN;
			}
			curminx = (-cursourY) / (objY - cursourY) * ((minXIdx - XN / 2.0) * dx - cursourX) + cursourX;
			for (ii = minXIdx; ii < maxXIdx; ++ii)
			{

				curmaxx = (-cursourY) / (objY - cursourY) * ((ii + 1 - XN / 2.0) * dx - cursourX) + cursourX;
				if (detPosLeft > detPosRight)
				{
					std::swap(detPosLeft, detPosRight);
				}
				summ += img[jj * XN + ii] * intersectLength<double>(detPosLeft, detPosRight, curminx, curmaxx) * dy;
				curminx = curmaxx;
			}
		}
		proj[angIdx * DN + detIdx] = summ / (cosAng * detprojLength);
		initYLeft = initYRight;
		curDetXLeft = curDetXRight;
		curDetYLeft = curDetYRight;
	}
}


template<typename T>
void DDM_ED_proj_template(std::vector<T>& proj, const std::vector<T>& img,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY,
	const T detSize,
	const T detCntIdx,
	const int XN, const int YN, const int DN, const int PN,
	const std::vector<T>& angs)
{
	//这个函数以x轴正方向为0度方向
	int angIdx(0), detIdx(0);
	T curang(0), cosT(0), sinT(0);
	T cursourx(0), cursoury(0);
	const T dd = detSize / static_cast<T>(DN);
	const T dx = objSizeX / static_cast<T>(XN);
	const T dy = objSizeY / static_cast<T>(YN);

#pragma omp parallel for private(curang, cosT, sinT, cursourx, cursoury)
	for (int angIdx = 0; angIdx < PN; ++angIdx)
	{
		curang = angs[angIdx];// angBeg + angIdx * angStp;
		cosT = std::cos(curang);
		sinT = std::sin(curang);
		cursourx = S2O * cosT; //以X轴为准;
		cursoury = S2O * sinT;

		if ((curang > CONSTVAL<T>::_PI_4 && curang <= CONSTVAL<T>::_3PI_4) || (curang >= CONSTVAL<T>::_5PI_4 && curang < CONSTVAL<T>::_7PI_4)) //按照角度来计算;
		{
			//the ray propagate along Y axis
			projCase1<T>(proj, img, cursourx, cursoury, S2O, O2D,
				objSizeX, objSizeY, detSize, detCntIdx,
				XN, YN, DN, PN, dx, dy, dd, curang, cosT, sinT, angIdx);
		}
		else
		{
			//the ray propagate along X axis
			projCase4<T>(proj, img, cursourx, cursoury, S2O, O2D,
				objSizeX, objSizeY, detSize, detCntIdx,
				XN, YN, DN, PN, dx, dy, dd, curang, cosT, sinT, angIdx);
		}
	}
}


void DDM_ED_proj(std::vector<double>& proj, const std::vector<double>& img,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY,
	const double detSize,
	const double detCntIdx,
	const int XN, const int YN, const int DN, const int PN,
	const std::vector<double>& angs)
{
	DDM_ED_proj_template<double>(proj, img, S2O, O2D, objSizeX, objSizeY, detSize, detCntIdx, XN, YN, DN, PN,
		angs);
}


void DDM_ED_proj(std::vector<float>& proj, const std::vector<float>& img,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY,
	const float detSize,
	const float detCntIdx,
	const int XN, const int YN, const int DN, const int PN,
	const std::vector<float>& angs)
{
	DDM_ED_proj_template<float>(proj, img, S2O, O2D, objSizeX, objSizeY, detSize, detCntIdx, XN, YN, DN, PN,
		angs);
}





//一个角度的投影;
template<typename T>
void DDM3D_ED_proj_template(std::vector<T>& proj, const std::vector<T>& vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detSizeU, const T detSizeV,
	const T detCntIdxU, const T detCntIdxV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<T>& angs)
{
	const T dx = objSizeX / XN;
	const T dy = objSizeY / YN;
	const T dz = objSizeZ / ZN;
	const T ddu = detSizeU / DNU;
	const T ddv = detSizeV / DNV;

	int angIdx;
#pragma omp parallel for
	for (angIdx = 0; angIdx < PN; angIdx++)
	{
		T curang = angs[angIdx];
		T cosT = cos(curang);
		T sinT = sin(curang);

		T cursourx = S2O * cosT;
		T cursoury = S2O * sinT;
		T cursourz = 0;

		if ((curang > CONSTVAL<T>::_PI_4 && curang <= CONSTVAL<T>::_3PI_4) || (curang >= CONSTVAL<T>::_5PI_4 && curang < CONSTVAL<T>::_7PI_4))
		{
			for (int detIdV = 0; detIdV < DNV; ++detIdV)
			{
				for (int detIdU = 0; detIdU < DNU; ++detIdU)
				{
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
						std::swap(detPosLX, detPosRX);
					}
					//假设下边的小;
					if (detPosDZ > detPosUZ)
					{
						std::swap(detPosDZ, detPosUZ);
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
							std::swap(minX, maxX);
						}
						if (minZ > maxZ)
						{
							std::swap(minZ, maxZ);
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
							T intersectL = intersectLength<double>(detPosLX, detPosRX, curminx, curmaxx);
							if (intersectL > 0)
							{
								for (size_t kk = minZIdx; kk < maxZIdx; kk++)
								{

									T curminz = (cursourz - (kk - ZN / 2.0) * dz) * cursoury / (objY - cursoury) + cursourz;
									T curmaxz = (cursourz - ((kk + 1) - ZN / 2.0) * dz) * cursoury / (objY - cursoury) + cursourz;

									T intersectH = intersectLength<double>(detPosDZ, detPosUZ, curminz, curmaxz);
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
		else
		{
			for (int detIdV = 0; detIdV < DNV; ++detIdV)
			{
				for (int detIdU = 0; detIdU < DNU; ++detIdU)
				{
					T summ = 0;
					T initDetX = -O2D;
					T initDetY = (detIdU - detCntIdxU) * ddu;
					T initDetZ = (detIdV - detCntIdxV) * ddv;

					T initDetLY = (detIdU - detCntIdxU - 0.5) * ddu;
					T initDetLZ = (detIdV - detCntIdxV) * ddv;

					T initDetRY = (detIdU - detCntIdxU + 0.5) * ddu;
					T initDetRZ = (detIdV - detCntIdxV) * ddv;

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
						std::swap(detPosLY, detPosRY);
					}
					//假设下边的小;
					if (detPosDZ > detPosUZ)
					{
						std::swap(detPosDZ, detPosUZ);
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
							std::swap(minY, maxY);
						}
						if (minZ > maxZ)
						{
							std::swap(minZ, maxZ);
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
							T intersectL = intersectLength<double>(detPosLY, detPosRY, curminy, curmaxy);
							if (intersectL > 0)
							{
								for (size_t kk = minZIdx; kk < maxZIdx; kk++)
								{

									T curminz = (cursourz - (kk - ZN / 2.0) * dz) * cursourx / (objX - cursourx) + cursourz;
									T curmaxz = (cursourz - ((kk + 1) - ZN / 2.0) * dz) * cursourx / (objX - cursourx) + cursourz;

									T intersectH = intersectLength<double>(detPosDZ, detPosUZ, curminz, curmaxz);
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

	}

}

void DDM3D_ED_proj(std::vector<double>& proj, const std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeU, const double detSizeV,
	const double detCntIdxU, const double detCntIdxV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<double>& angs)
{
	DDM3D_ED_proj_template(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ,
		detSizeU, detSizeV, detCntIdxU, detCntIdxV, XN, YN, ZN, DNU, DNV, PN, angs);
}

void DDM3D_ED_proj(std::vector<float>& proj, const std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeU, const float detSizeV,
	const float detCntIdxU, const float detCntIdxV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<float>& angs)
{
	DDM3D_ED_proj_template(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ,
		detSizeU, detSizeV, detCntIdxU, detCntIdxV, XN, YN, ZN, DNU, DNV, PN, angs);
}



template<typename T>
void DDM_ED_bproj_template(const std::vector<T>& proj, std::vector<T>& img,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY,
	const T detSize,
	const T detCntIdx,
	const int XN, const int YN, const int DN, const int PN,
	const std::vector<T>& angs)
{
	//这个函数以x轴正方向为0度方向
	//int pixelIdx(0), angIdx;
	//T curang(0), cosT(0), sinT(0);
	T cursourx(0), cursoury(0);
	const T dd = detSize / static_cast<T>(DN);
	const T dx = objSizeX / static_cast<T>(XN);
	const T dy = objSizeY / static_cast<T>(YN);
	const int OBJN = XN * YN;
	const T hfXN = XN * 0.5;
	const T hfYN = YN * 0.5;
	const T S2D = S2O + O2D;
	//int ii, jj;
	//T lefPx, lefPy, rghPx, rghPy;
	//T initObjX1, initObjY1, initObjX2, initObjY2;
	//T objYdetPosMin, objYdetPosMax, objYaxisPosMin, objYaxisPosMax;

	//int minDetIdx, maxDetIdx;
	//T summ = 0;
	//T objYaxisLength;
	//T maxDetPos, minDetPos;
	//T cosAng;
	//T ll;

#pragma omp parallel for//private(initObjX1, initObjX2, initObjY1, initObjY2,objYdetPosMin,objYdetPosMax,minDetIdx, maxDetIdx,objYaxisPosMin, objYaxisPosMax,objYaxisLength,maxDetPos,minDetPos,ll)
	for (int pixelIdx = 0; pixelIdx < OBJN; ++pixelIdx)
	{
		//summ = 0;
		int jj = pixelIdx / XN;
		int ii = pixelIdx - jj * XN;
		////计算当前像素四个中点的坐标;
		T summ = 0;

		for (int angIdx = 0; angIdx != angs.size(); ++angIdx)
		{
			T curang = angs[angIdx];
			T cosT = cos(curang);
			T sinT = sin(curang);
			T lefPx(0); //= (ii - hfXN) * dx;
			T lefPy(0); //= (jj - hfYN + 0.5) * dy;
			T rghPx(0); //= (ii - hfXN + 1.0) * dx;
			T rghPy(0);// = (jj - hfYN + 0.5) * dy;

			if ((curang >  CONSTVAL<T>::_PI_4 && curang <= CONSTVAL<T>::_3PI_4) || (curang >= CONSTVAL<T>::_5PI_4 && curang < CONSTVAL<T>::_7PI_4))
			{
				lefPx = (ii - hfXN) * dx;
				lefPy = (jj - hfYN + 0.5) * dy;
				rghPx = (ii - hfXN + 1.0) * dx;
				rghPy = (jj - hfYN + 0.5) * dy;
			}
			else
			{
				lefPx = (ii - hfXN + 0.5) * dx;
				lefPy = (jj - hfYN + 1.0) * dy;
				rghPx = (ii - hfXN + 0.5) * dx;
				rghPy = (jj - hfYN) * dy;
			}

			T initObjX1 = lefPx * cosT + lefPy * sinT;
			T initObjY1 = -lefPx * sinT + lefPy * cosT;
			T initObjX2 = rghPx * cosT + rghPy * sinT;
			T initObjY2 = -rghPx * sinT + rghPy * cosT;

			T objYdetPosMin = initObjY1 * S2D / (S2O - initObjX1);
			T objYdetPosMax = initObjY2 * S2D / (S2O - initObjX2);

			if (objYdetPosMax < objYdetPosMin)
			{
				std::swap(objYdetPosMin, objYdetPosMax);
			}
			int minDetIdx = floor(objYdetPosMax / (-dd) + detCntIdx);
			int maxDetIdx = ceil(objYdetPosMin / (-dd) + detCntIdx);

			if (minDetIdx > DN)
			{
				continue;
			}
			if (maxDetIdx < 0)
			{
				continue;
			}

			T objYaxisPosMin = initObjX1 * initObjY1 / (S2O - initObjX1) + initObjY1; //pixel端点在Y轴上的投影;
			T objYaxisPosMax = initObjX2 * initObjY2 / (S2O - initObjX2) + initObjY2;
			if (objYaxisPosMax < objYaxisPosMin)
			{
				std::swap(objYaxisPosMin, objYaxisPosMax);
			}
			T objYaxisLength = abs(objYaxisPosMax - objYaxisPosMin);

			if (minDetIdx < 0)
			{
				minDetIdx = 0;
			}
			if (maxDetIdx >= DN)
			{
				maxDetIdx = DN;
			}

			for (int detIdx = minDetIdx; detIdx < maxDetIdx; ++detIdx)
			{
				T maxDetPos = (-(detIdx - detCntIdx) * dd) * S2O / S2D;
				T minDetPos = (-(detIdx + 1.0 - detCntIdx) * dd) * S2O / S2D;
				if (maxDetPos < minDetPos)
				{
					std::swap(minDetPos, maxDetPos);
				}
				T ll = hypot(S2O, (-(detIdx + 0.5 - detCntIdx) * dd) * S2O / S2D);

				T cosAng = abs(S2O / ll);
				summ += proj[angIdx * DN + detIdx] * intersectLength<T>(objYaxisPosMin, objYaxisPosMax, minDetPos, maxDetPos) / (objYaxisLength * cosAng);
			}

			if ((curang > CONSTVAL<T>::_PI_4 && curang <= CONSTVAL<T>::_3PI_4) || (curang >= CONSTVAL<T>::_5PI_4 && curang < CONSTVAL<T>::_7PI_4))
			{
				summ *= dy;
			}
			else
			{
				summ *= dx;
			}


		}
		img[pixelIdx] = summ;
	}
}


void DDM_ED_bproj(const std::vector<double>& proj, std::vector<double>& img,
	const double S2O, const double O2D, const double objSizeX, const double objSizeY,
	const double detSize, const double detCntIdx, const int XN, const int YN, const int DN, const int PN,
	const std::vector<double>& angs)
{
	DDM_ED_bproj_template<double>(proj, img, S2O, O2D, objSizeX, objSizeY,
		detSize, detCntIdx, XN, YN, DN, PN, angs);
}


void DDM_ED_bproj(const std::vector<float>& proj, std::vector<float>& img,
	const float S2O, const float O2D, const float objSizeX, const float objSizeY,
	const float detSize, const float detCntIdx, const int XN, const int YN, const int DN, const int PN,
	const std::vector<float>& angs)
{
	DDM_ED_bproj_template<float>(proj, img, S2O, O2D, objSizeX, objSizeY,
		detSize, detCntIdx, XN, YN, DN, PN, angs);
}





template<typename T>
void DDM3D_bproj_template(const std::vector<T>& proj, std::vector<T>& vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detSizeU, const T detSizeV,
	const T detCntIdU, const T detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<T>& angs)
{
	//以X轴正方向为0度方向;
	const T ddu = detSizeU / static_cast<T>(DNU);
	const T ddv = detSizeV / static_cast<T>(DNV);
	const T dx = objSizeX / static_cast<T>(XN);
	const T dy = objSizeY / static_cast<T>(YN);
	const T dz = objSizeZ / static_cast<T>(ZN);
	const T hfXN = XN * 0.5;
	const T hfYN = YN * 0.5;
	const T hfZN = ZN * 0.5;
	const int VOLN = XN * YN * ZN;
	const T S2D = S2O + O2D;
#pragma omp parallel for
	for (int voxelIdx = 0; voxelIdx < VOLN; ++voxelIdx)
	{
		int kk = voxelIdx / (XN * YN);
		int jj = (voxelIdx - kk * XN * YN) / XN;
		int ii = voxelIdx - kk * XN * YN - jj * XN;
		T summ = 0;
		for (int angIdx = 0; angIdx != angs.size(); ++angIdx)
		{
			T curang = angs[angIdx];
			T cosT = cos(curang);
			T sinT = sin(curang);
			T Px = (ii - hfXN + 0.5) * dx;
			T Py = (jj - hfYN + 0.5) * dy;
			T Pz = (kk - hfZN + 0.5) * dz;

			T lefPx(0);
			T lefPy(0);
			//T lefPz(0);

			T rghPx(0);
			T rghPy(0);
			//T rghPz(0);

			T uppPx(0);
			T uppPy(0);
			T uppPz(0);

			T dowPx(0);
			T dowPy(0);
			T dowPz(0);

			if ((curang > CONSTVAL<T>::_PI_4 && curang <= CONSTVAL<T>::_3PI_4) || (curang >= CONSTVAL<T>::_5PI_4 && curang < CONSTVAL<T>::_7PI_4))
			{
				lefPx = Px - 0.5 * dx;
				lefPy = Py;
				//lefPz = Pz;

				rghPx = Px + 0.5 * dx;
				rghPy = Py;
				//rghPz = Pz;

				uppPx = Px;
				uppPy = Py;
				uppPz = Pz + 0.5 * dz;

				dowPx = Px;
				dowPy = Py;
				dowPz = Pz - 0.5 * dz;

			}
			else
			{
				lefPx = Px;
				lefPy = Py - 0.5 * dy;
				//lefPz = Pz;

				rghPx = Px;
				rghPy = Py + 0.5 * dy;
				//rghPz = Pz;

				uppPx = Px;
				uppPy = Py;
				uppPz = Pz + 0.5 * dz;

				dowPx = Px;
				dowPy = Py;
				dowPz = Pz - 0.5 * dz;
			}

			T initObjlefx = lefPx * cosT + lefPy * sinT;
			T initObjlefy = -lefPx * sinT + lefPy * cosT;
			//T initObjlefz = lefPz;

			T initObjrghx = rghPx * cosT + rghPy * sinT;
			T initObjrghy = -rghPx * sinT + rghPy * cosT;
			//T initObjrghz = rghPz;

			T initObjuppx = uppPx * cosT + uppPy * sinT;
			//T initObjuppy = -uppPx * sinT + uppPy * cosT;
			T initObjuppz = uppPz;

			T initObjdowx = dowPx* cosT + dowPy * sinT;
			//T initObjdowy = -dowPx * sinT + dowPy * cosT;
			T initObjdowz = dowPz;

			T objYdetPosUMin = initObjlefy * S2D / (S2O - initObjlefx);
			T objYdetPosUMax = initObjrghy * S2D / (S2O - initObjrghx);
			T objYdetPosVMin = initObjdowz * S2D / (S2O - initObjdowx);
			T objYdetPosVMax = initObjuppz * S2D / (S2O - initObjuppx);

			if (objYdetPosUMin > objYdetPosUMax)
			{
				std::swap(objYdetPosUMin, objYdetPosUMax);
			}
			if (objYdetPosVMin > objYdetPosVMax)
			{
				std::swap(objYdetPosVMin, objYdetPosVMax);
			}
			int minDetUIdx = floor(objYdetPosUMin / ddu + detCntIdU) - 1;
			int maxDetUIdx = ceil(objYdetPosUMax / ddu + detCntIdU) + 1;
			int minDetVIdx = floor(objYdetPosVMin / ddv + detCntIdV) - 1;
			int maxDetVIdx = ceil(objYdetPosVMax / ddv + detCntIdV) + 1;

			if (minDetUIdx > DNU)
			{
				continue;
			}
			if (maxDetUIdx < 0)
			{
				continue;
			}
			if (minDetVIdx > DNV)
			{
				continue;
			}
			if (maxDetVIdx < 0)
			{
				continue;
			}

			T objYOZLength = objYdetPosUMax - objYdetPosUMin;
			T objYOZHeight = objYdetPosVMax - objYdetPosVMin;


			if (minDetUIdx < 0)
			{
				minDetUIdx = 0;
			}
			if (maxDetUIdx > DNU)
			{
				maxDetUIdx = DNU;
			}
			if (minDetVIdx < 0)
			{
				minDetVIdx = 0;
			}
			if (maxDetVIdx > DNV)
			{
				maxDetVIdx = DNV;
			}


			for (int detIdU = minDetUIdx; detIdU < maxDetUIdx; ++detIdU)
			{
				T minDetUPos = (detIdU - detCntIdU - 0.5) * ddu;// *S2O / S2D;
				T maxDetUPos = (detIdU - detCntIdU + 0.5) * ddu;// *S2O / S2D;

				T ll = intersectLength<T>(objYdetPosUMin, objYdetPosUMax, minDetUPos, maxDetUPos);
				if (ll > 0)
				{
					for (int detIdV = minDetVIdx; detIdV < maxDetVIdx; ++detIdV)
					{
						T minDetVPos = (detIdV - detCntIdV - 0.5) * ddv;// *S2O / S2D;
						T maxDetVPos = (detIdV - detCntIdV + 0.5) * ddv;// *S2O / S2D;

						T DU = (detIdU - detCntIdU) * ddu;
						T DV = (detIdV - detCntIdV) * ddv;
						T cosAlphacosGamma = S2D / sqrt(DU*DU + DV*DV + S2D*S2D);
						T mm = intersectLength<T>(objYdetPosVMin, objYdetPosVMax, minDetVPos, maxDetVPos);
						if (mm > 0)
						{
							summ += (proj[(angIdx * DNV + detIdV) * DNU + detIdU] * ll * mm / (objYOZLength * objYOZHeight * cosAlphacosGamma) * dx);
						}
						else
						{
							summ += 0;
						}

					}
				}

			}


		}
		vol[voxelIdx] = summ;
	}
}

void DDM3D_bproj(const std::vector<double>& proj, std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeU, const double detSizeV,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<double>& angs)
{
	DDM3D_bproj_template<double>(proj, vol, S2O, O2D,
		objSizeX, objSizeY, objSizeZ,
		detSizeU, detSizeV, detCntIdU, detCntIdV,
		XN, YN, ZN, DNU, DNV, PN, angs);
}


void DDM3D_bproj(const std::vector<float>& proj, std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeU, const float detSizeV,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<float>& angs)
{
	DDM3D_bproj_template<float>(proj, vol, S2O, O2D,
		objSizeX, objSizeY, objSizeZ,
		detSizeU, detSizeV, detCntIdU, detCntIdV,
		XN, YN, ZN, DNU, DNV, PN, angs);
}





template<typename T>
void pushCase4(std::vector<int>& rowIdx,
	std::vector<int>& colIdx,
	std::vector<T>& weight,
	const T cursourX, const T cursourY, const T S2O, const T O2D, const T objSizeX, const T objSizeY, const T detSize, const T detCntIdx, const int XN, const int YN, const int DN, const int PN,
	const T dx, const T dy, const T dd, const T curAng, const T cosT, const T sinT, const int angIdx)
{
	T summ = 0;
	T initX, initY, initYLeft, initYRight, curDetXLeft, curDetXRight, curDetYLeft, curDetYRight;
	T curDetX, curDetY, dirX, dirY, legth, cosAng, detPosLeft, detPosRight;
	T detprojLength;
	T objX;

	T minY, maxY;
	int minYIdx, maxYIdx, detIdx, ii, jj;
	T curminy, curmaxy;
	initX = -O2D;
	int ridx;

	T w;
	initYLeft = -(0 - detCntIdx - 0.5) * dd; 
	curDetXLeft = initX * cosT - initYLeft * sinT; 
	curDetYLeft = initX * sinT + initYLeft * cosT; 

	for (detIdx = 0; detIdx != DN; ++detIdx)
	{
		ridx = angIdx * DN + detIdx;

		summ = 0;

		initYRight = -(detIdx - detCntIdx - 0.5 + 1) * dd; 
		initY = -(detIdx - detCntIdx) * dd; 
		curDetXRight = initX * cosT - initYRight * sinT; 
		curDetYRight = initX * sinT + initYRight * cosT; 
		curDetX = initX * cosT - initY * sinT; 
		curDetY = initX * sinT + initY * cosT; 
		dirX = curDetX - cursourX;
		dirY = curDetY - cursourY;
		legth = hypot(dirX, dirY);
		dirX /= legth;
		dirY /= legth;
		cosAng = abs(dirX); 

		detPosLeft = (0 - cursourX) / (curDetXLeft - cursourX) * (curDetYLeft - cursourY) + cursourY; 
		detPosRight = (0 - cursourX) / (curDetXRight - cursourX) * (curDetYRight - cursourY) + cursourY;
		detprojLength = abs(detPosRight - detPosLeft);

		for (ii = 0; ii < XN; ++ii)
		{

			objX = (ii - YN / 2.0 + 0.5) * dy;
			minY = (objX - cursourX) / (curDetXLeft - cursourX) * (curDetYLeft - cursourY) + cursourY;
			maxY = (objX - cursourX) / (curDetXRight - cursourX) *  (curDetYRight - cursourY) + cursourY;
			if (minY > maxY)
			{
				std::swap(minY, maxY);
			}

			minYIdx = floor(minY / dy + YN / 2.0);
			maxYIdx = ceil(maxY / dy + YN / 2.0);

			if (maxYIdx <= 0)
			{
				continue;
			}
			else if (minYIdx > XN)
			{
				continue;
			}

			if (minYIdx < 0)
			{
				minYIdx = 0;
			}
			if (maxYIdx > YN)
			{
				maxYIdx = YN;
			}


			curminy = (-cursourX) / (objX - cursourX) * ((0 - YN / 2.0) * dy - cursourY) + cursourY;
			for (jj = minYIdx; jj < maxYIdx; ++jj)
			{
				curmaxy = (-cursourX) / (objX - cursourX) * ((jj + 1 - YN / 2.0) * dy - cursourY) + cursourY;
				if (detPosLeft > detPosRight)
				{
					std::swap(detPosLeft, detPosRight);
				}

				w = intersectLength<double>(detPosLeft, detPosRight, curminy, curmaxy);
				if (w > 0)
				{
					rowIdx.push_back(ridx);
					colIdx.push_back(jj * XN + ii);
					weight.push_back(w * dx / (cosAng * detprojLength));
				}
				//summ += img[jj * XN + ii] * * dx;
				curminy = curmaxy;
			}
		}

		//proj[angIdx * DN + detIdx] = summ / (cosAng * detprojLength);
		initYLeft = initYRight;
		curDetXLeft = curDetXRight;
		curDetYLeft = curDetYRight;
	}
}


template<typename T>
void pushCase1(std::vector<int>& rowIdx,
	std::vector<int>& colIdx,
	std::vector<T>& weight,
	const T cursourX, const T cursourY, const T S2O, const T O2D, const T objSizeX, const T objSizeY, const T detSize, const T detCntIdx, const int XN, const int YN, const int DN, const int PN,
	const T dx, const T dy, const T dd, const T curAng, const T cosT, const T sinT, const int angIdx)
{
	int detIdx = 0;
	int ii = 0, jj = 0;
	T initX, initYLeft, initYRight;
	T curDetXLeft, curDetYLeft;
	T curDetXRight, curDetYRight;
	T minX, maxX;
	int minXIdx, maxXIdx;
	T detPosLeft;
	T detPosRight;
	T initY;
	T curDetX, curDetY;
	T dirX, dirY;
	T legth;
	T cosAng;
	T detprojLength = 0;
	T objY;
	T summ = 0;
	T curminx, curmaxx;
	int ridx;
	T w;
	initX = -O2D;
	initYLeft = -(0 - detCntIdx - 0.5) * dd;
	curDetXLeft = initX * cosT - initYLeft * sinT;
	curDetYLeft = initX * sinT + initYLeft * cosT;
	for (detIdx = 0; detIdx != DN; ++detIdx)
	{
		ridx = angIdx * DN + detIdx;
		summ = 0;
		initYRight = initYLeft - dd;// -(detIdx - detCntIdx - 0.5 + 1) * dd;
		initY = -(detIdx - detCntIdx) * dd;
		curDetXRight = initX * cosT - initYRight * sinT;
		curDetYRight = initX * sinT + initYRight * cosT;
		curDetX = initX * cosT - initY * sinT;
		curDetY = initX * sinT + initY * cosT;
		dirX = curDetX - cursourX;
		dirY = curDetY - cursourY;
		legth = hypot(dirX, dirY);
		dirX /= legth;
		dirY /= legth;
		cosAng = abs(dirY);

		detPosLeft = (0 - cursourY) / (curDetYLeft - cursourY) * (curDetXLeft - cursourX) + cursourX; //det左边界X轴上的坐标;
		detPosRight = (0 - cursourY) / (curDetYRight - cursourY) * (curDetXRight - cursourX) + cursourX;//det右边界在x轴上的坐标;
		detprojLength = abs(detPosRight - detPosLeft);

		for (jj = 0; jj < YN; jj++)
		{
			objY = (jj - YN / 2.0 + 0.5) * dy;
			minX = (objY - cursourY) / (curDetYLeft - cursourY) * (curDetXLeft - cursourX) + cursourX;
			maxX = (objY - cursourY) / (curDetYRight - cursourY) *  (curDetXRight - cursourX) + cursourX;
			if (minX > maxX)
			{
				std::swap(minX, maxX);
			}

			minXIdx = floor(minX / dx + XN / 2.0);
			maxXIdx = ceil(maxX / dx + XN / 2.0);

			if (maxXIdx <= 0)
			{
				continue;
			}
			else if (minXIdx > XN)
			{
				continue;
			}

			if (minXIdx < 0)
			{
				minXIdx = 0;
			}
			if (maxXIdx > XN)
			{
				maxXIdx = XN;
			}
			curminx = (-cursourY) / (objY - cursourY) * ((minXIdx - XN / 2.0) * dx - cursourX) + cursourX;
			for (ii = minXIdx; ii < maxXIdx; ++ii)
			{

				curmaxx = (-cursourY) / (objY - cursourY) * ((ii + 1 - XN / 2.0) * dx - cursourX) + cursourX;
				if (detPosLeft > detPosRight)
				{
					std::swap(detPosLeft, detPosRight);
				}

				w = intersectLength<double>(detPosLeft, detPosRight, curminx, curmaxx);
				if (w > 0)
				{
					rowIdx.push_back(ridx);
					colIdx.push_back(jj * XN + ii);
					weight.push_back(w * dx / (cosAng * detprojLength));
				}

				//summ += img[jj * XN + ii] * intersectLength<double>(detPosLeft, detPosRight, curminx, curmaxx) * dy;
				curminx = curmaxx;
			}
		}
		//proj[angIdx * DN + detIdx] = summ / (cosAng * detprojLength);
		initYLeft = initYRight;
		curDetXLeft = curDetXRight;
		curDetYLeft = curDetYRight;
	}
}



template<typename T>
void genMatrix_DDM_ED_template(
	std::vector<int>& rowIdx,
	std::vector<int>& colIdx,
	std::vector<T>& weight,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY,
	const T detSize,
	const T detCntIdx,
	const int XN, const int YN, const int DN, const int PN,
	const std::vector<T>& angs)
{
	int angIdx(0), detIdx(0);
	T curang(0), cosT(0), sinT(0);
	T cursourx(0), cursoury(0);
	const T dd = detSize / static_cast<T>(DN);
	const T dx = objSizeX / static_cast<T>(XN);
	const T dy = objSizeY / static_cast<T>(YN);
	for (angIdx = 0; angIdx < PN; ++angIdx)
	{
		curang = angs[angIdx];// angBeg + angIdx * angStp;
		cosT = cos(curang);
		sinT = sin(curang);
		cursourx = S2O * cosT;
		cursoury = S2O * sinT;

		if ((curang > CONSTVAL<T>::_PI_4 && curang <= CONSTVAL<T>::_3PI_4) || (curang >= CONSTVAL<T>::_5PI_4 && curang < CONSTVAL<T>::_7PI_4))
		{
			pushCase1<T>(rowIdx, colIdx, weight, cursourx, cursoury, S2O, O2D,
				objSizeX, objSizeY, detSize, detCntIdx,
				XN, YN, DN, PN, dx, dy, dd, curang, cosT, sinT, angIdx);
		}
		else
		{
			pushCase4<T>(rowIdx, colIdx, weight, cursourx, cursoury, S2O, O2D,
				objSizeX, objSizeY, detSize, detCntIdx,
				XN, YN, DN, PN, dx, dy, dd, curang, cosT, sinT, angIdx);
		}
	}


}



void genMatrix_DDM_ED(
	std::vector<int>& rowIdx,
	std::vector<int>& colIdx,
	std::vector<double>& weight,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY,
	const double detSize,
	const double detCntIdx,
	const int XN, const int YN, const int DN, const int PN,
	const std::vector<double>& angs)
{
	genMatrix_DDM_ED_template<double>(rowIdx, colIdx, weight, S2O, O2D, objSizeX, objSizeY, detSize, detCntIdx, XN, YN, DN, PN, angs);
}



void genMatrix_DDM_ED(
	std::vector<int>& rowIdx,
	std::vector<int>& colIdx,
	std::vector<float>& weight,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY,
	const float detSize,
	const float detCntIdx,
	const int XN, const int YN, const int DN, const int PN,
	const std::vector<float>& angs)
{
	genMatrix_DDM_ED_template<float>(rowIdx, colIdx, weight, S2O, O2D, objSizeX, objSizeY, detSize, detCntIdx, XN, YN, DN, PN, angs);
}


