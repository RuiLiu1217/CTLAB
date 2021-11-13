/*
 * COPYRIGHT NOTICE
 * COPYRIGHT (c) 2015, Wake Forest and UMass Lowell
 * All rights reserved
 *
 * @file DDM_proj.cpp
 * @brief The CPU based DD projection in brute force method
 *
 * @version 1.0
 * @author Rui Liu
 * @date May. 1, 2015
 *
 */


#include "DDM_proj.h"
#include "utilities.hpp"


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

		if ((curang > _PI_4 && curang <= _3PI_4) || (curang >= _5PI_4 && curang < _7PI_4)) //按照角度来计算;
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





//DDM projection from one angle
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

		if ((curang > _PI_4 && curang <= _3PI_4) || (curang >= _5PI_4 && curang < _7PI_4))
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
void projCase1_EA(std::vector<T>& proj, const std::vector<T>& img, const T cursourX, const T cursourY,
	const T S2O, const T O2D, const T objSizeX, const T objSizeY, const T detArc, const T detCntIdx,
	const int XN, const int YN, const int DN, const int PN, const T dx, const T dy, const T dd,
	const T curAng, const T cosT, const T sinT, const int angIdx)
{
	const T S2D = S2O + O2D;

	T beta;// = (0 - detCntIdx) * dd;
	T cosBeta;// = cos(beta);
	T sinBeta;// = sin(beta);

	T initX;// = (-O2D) * cosBeta - 0 * sinBeta + S2O;
	T initY;// = (-O2D) * sinBeta + 0 * cosBeta;

	T curDetX;
	T curDetY;

	T initXLeft;// = (-O2D) * cosBeta - 0 * sinBeta + S2O;
	T initYLeft;// = (-O2D) * sinBeta + 0 * cosBeta;

	T initXRight(0);
	T initYRight(0);

	T curDetXLeft;// = initXLeft * cosT - initYLeft * sinT;
	T curDetYLeft;// = initXLeft * sinT + initYLeft * cosT;

	T curDetXRight(0);
	T curDetYRight(0);

	T dirX, dirY;

	T legth;
	T detPosLeft;
	T detPosRight;
	T detprojLength;
	T summ;
	T cosAng;
	T objY, minX, maxX;
	int minXIdx, maxXIdx;
	T curminx, curmaxx;
	int detIdx, ii, jj;

	for (detIdx = 0; detIdx != DN; ++detIdx)
	{
		summ = 0;
		beta = (detIdx - detCntIdx) * dd;
		cosBeta = cos(beta);
		sinBeta = sin(beta);
		initX = (-O2D) * cosBeta - 0 * sinBeta + S2O;
		initY = (-O2D) * sinBeta + 0 * cosBeta;
		curDetX = initX * cosT - initY * sinT;
		curDetY = initX * sinT + initY * cosT;
		//×ó±ß;
		beta = (detIdx - detCntIdx - 0.5) * dd;
		cosBeta = cos(beta);
		sinBeta = sin(beta);
		initXLeft = (-O2D) * cosBeta - 0 * sinBeta + S2O;
		initYLeft = (-O2D) * sinBeta + 0 * cosBeta;
		curDetXLeft = initXLeft * cosT - initYLeft * sinT;
		curDetYLeft = initXLeft * sinT + initYLeft * cosT;

		//ÓÒ±ß;
		beta = (detIdx - detCntIdx + 0.5) * dd;
		cosBeta = cos(beta);
		sinBeta = sin(beta);
		initXRight = (-O2D) * cosBeta - 0 * sinBeta + S2O;
		initYRight = (-O2D) * sinBeta + 0 * cosBeta;
		curDetXRight = initXRight * cosT - initYRight * sinT;
		curDetYRight = initXRight * sinT + initYRight * cosT;

		dirX = curDetX - cursourX;
		dirY = curDetY - cursourY;
		legth = hypot(dirX, dirY);
		dirX /= legth;
		dirY /= legth;

		cosAng = abs(dirY); //µ±Ç°¹âÏßºÍYÖáŒÐœÇÓàÏÒ;

		detPosLeft = (0 - cursourY) / (curDetYLeft - cursourY) * (curDetXLeft - cursourX) + cursourX;
		detPosRight = (0 - cursourY) / (curDetYRight - cursourY) *(curDetXRight - cursourX) + cursourX;
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
	}

}



template<typename T>
void projCase4_EA(std::vector<T>& proj, const std::vector<T>& img, const T cursourX, const T cursourY,
	const T S2O, const T O2D, const T objSizeX, const T objSizeY, const T detArc, const T detCntIdx,
	const int XN, const int YN, const int DN, const int PN, const T dx, const T dy, const T dd,
	const T curAng, const T cosT, const T sinT, const int angIdx)
{
	const T S2D = S2O + O2D;

	T beta;// = (0 - detCntIdx) * dd;
	T cosBeta;// = cos(beta);
	T sinBeta;// = sin(beta);

	T initX;// = (-O2D) * cosBeta - 0 * sinBeta + S2O;
	T initY;// = (-O2D) * sinBeta + 0 * cosBeta;

	T curDetX;
	T curDetY;

	T initXLeft;// = (-O2D) * cosBeta - 0 * sinBeta + S2O;
	T initYLeft;// = (-O2D) * sinBeta + 0 * cosBeta;

	T initXRight(0);
	T initYRight(0);

	T curDetXLeft;// = initXLeft * cosT - initYLeft * sinT;
	T curDetYLeft;// = initXLeft * sinT + initYLeft * cosT;

	T curDetXRight(0);
	T curDetYRight(0);

	T dirX, dirY;

	T legth;
	T detPosLeft;
	T detPosRight;
	T detprojLength;
	T summ;
	T cosAng;
	T objX, minY, maxY;
	int minYIdx, maxYIdx;
	T curminy, curmaxy;
	int detIdx, ii, jj;

	for (detIdx = 0; detIdx != DN; ++detIdx)
	{
		summ = 0;
		beta = (detIdx - detCntIdx) * dd;
		cosBeta = cos(beta);
		sinBeta = sin(beta);
		initX = (-O2D) * cosBeta - 0 * sinBeta + S2O;
		initY = (-O2D) * sinBeta + 0 * cosBeta;
		curDetX = initX * cosT - initY * sinT;
		curDetY = initX * sinT + initY * cosT;
		//×ó±ß;
		beta = (detIdx - detCntIdx - 0.5) * dd;
		cosBeta = cos(beta);
		sinBeta = sin(beta);
		initXLeft = (-O2D) * cosBeta - 0 * sinBeta + S2O;
		initYLeft = (-O2D) * sinBeta + 0 * cosBeta;
		curDetXLeft = initXLeft * cosT - initYLeft * sinT;
		curDetYLeft = initXLeft * sinT + initYLeft * cosT;

		//ÓÒ±ß;
		beta = (detIdx - detCntIdx + 0.5) * dd;
		cosBeta = cos(beta);
		sinBeta = sin(beta);
		initXRight = (-O2D) * cosBeta - 0 * sinBeta + S2O;
		initYRight = (-O2D) * sinBeta + 0 * cosBeta;
		curDetXRight = initXRight * cosT - initYRight * sinT;
		curDetYRight = initXRight * sinT + initYRight * cosT;

		dirX = curDetX - cursourX;
		dirY = curDetY - cursourY;
		legth = hypot(dirX, dirY);
		dirX /= legth;
		dirY /= legth;

		cosAng = abs(dirX); //ÓëCase1Çø±ð;
		//ÕâÀïÊÇŽÓX×ø±êËãY×ø±ê;
		detPosLeft = (0 - cursourX) / (curDetXLeft - cursourX) * (curDetYLeft - cursourY) + cursourY; //det×ó±ßœçXÖáÉÏµÄ×ø±ê;
		detPosRight = (0 - cursourX) / (curDetXRight - cursourX) * (curDetYRight - cursourY) + cursourY;//detÓÒ±ßœçÔÚxÖáÉÏµÄ×ø±ê;
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
				summ += img[jj * XN + ii] * intersectLength<double>(detPosLeft, detPosRight, curminy, curmaxy) * dx;
				curminy = curmaxy;
			}
		}
		proj[angIdx * DN + detIdx] = summ / (cosAng * detprojLength);
	}


}


template<typename T>
void DDM_EA_proj_template(std::vector<T>& proj, const std::vector<T>& img,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY,
	const T detArc,
	const T detCntIdx,
	const int XN, const int YN, const int DN, const int PN,
	const std::vector<T>& angs)
{
	int angIdx(0);
	T curang(0), cosT(0), sinT(0);
	T cursourx(0), cursoury(0);

	const T dd = detArc / static_cast<T>(DN);
	const T dx = objSizeX / static_cast<T>(XN);
	const T dy = objSizeY / static_cast<T>(YN);
#pragma omp parallel for private(curang, cosT, sinT, cursourx, cursoury)
	for (angIdx = 0; angIdx < PN; angIdx++)
	{
		curang = angs[angIdx];
		cosT = cos(curang);
		sinT = sin(curang);
		cursourx = S2O * cosT; //ÒÔXÖáÎª×Œ;
		cursoury = S2O * sinT;

		if ((curang > PI * 0.25 && curang <= PI * 0.75) || (curang >= PI * 1.25 && curang < PI * 1.75))
		{
			projCase1_EA(proj, img, cursourx, cursoury, S2O, O2D,
				objSizeX, objSizeY, detArc, detCntIdx, XN, YN, DN, PN, dx, dy, dd,
				curang, cosT, sinT, angIdx);
		}
		else
		{
			projCase4_EA(proj, img, cursourx, cursoury, S2O, O2D,
				objSizeX, objSizeY, detArc, detCntIdx, XN, YN, DN, PN, dx, dy, dd,
				curang, cosT, sinT, angIdx);
		}

	}
}


void DDM_EA_proj(std::vector<float>& proj, const std::vector<float>& img,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY,
	const float detArc,
	const float detCntIdx,
	const int XN, const int YN, const int DN, const int PN,
	const std::vector<float>& angs)
{
	DDM_EA_proj_template<float>(proj, img, S2O, O2D, objSizeX, objSizeY, detArc,
		detCntIdx, XN, YN, DN, PN, angs);
}
void DDM_EA_proj(std::vector<double>& proj, const std::vector<double>& img,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY,
	const double detArc,
	const double detCntIdx,
	const int XN, const int YN, const int DN, const int PN,
	const std::vector<double>& angs)
{
	DDM_EA_proj_template<double>(proj, img, S2O, O2D, objSizeX, objSizeY, detArc,
		detCntIdx, XN, YN, DN, PN, angs);
}







template<typename T>
void DDM3D_EA_proj_template(std::vector<T>& proj, const std::vector<T>& vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detArc, const T detSizeH,
	const T detCntIdU, const T detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<T>& angs)
{
	const T dx = objSizeX / static_cast<T>(XN);
	const T dy = objSizeY / static_cast<T>(YN);
	const T dz = objSizeZ / static_cast<T>(ZN);
	const T dbeta = detArc / static_cast<T>(DNU);
	const T ddv = detSizeH / static_cast<T>(DNV);


	int angIdx(0);
	//#pragma omp parallel for
	for (angIdx = 0; angIdx < PN; ++angIdx)
	{
		T curang = angs[angIdx];
		T cosT = cos(curang);
		T sinT = sin(curang);

		T cursourx = S2O * cosT;
		T cursoury = S2O * sinT;
		T cursourz = 0;

		if ((curang > PI * 0.25 && curang <= PI * 0.75) || (curang >= PI * 1.25 && curang < PI * 1.75))
		{
			for (int detIdV = 0; detIdV < DNV; ++detIdV)
			{
				for (int detIdU = 0; detIdU < DNU; ++detIdU)
				{
					T summ = 0;
					T beta = (detIdU - detCntIdU) * dbeta; //ÖÐŒä;
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
					T curDetDZ = initDetDZ;

					T curDetUX = curDetDX;
					T curDetUY = curDetDY;
					T curDetUZ = initDetUZ;

					T curDetX = initDetX * cosT - initDetY * sinT;
					T curDetY = initDetX * sinT + initDetY * cosT;

					T dirX = curDetX - cursourx;
					T dirY = curDetY - cursoury;

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
						std::swap(detPosLX, detPosRX);
					}
					//ŒÙÉèÏÂ±ßµÄÐ¡;
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

		} // End if
		else
		{
			for (int detIdV = 0; detIdV < DNV; ++detIdV)
			{
				for (int detIdU = 0; detIdU < DNU; ++detIdU)
				{
					T summ = 0;
					T beta = (detIdU - detCntIdU) * dbeta; //ÖÐŒä;
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
					T curDetDZ = initDetDZ;

					T curDetUX = curDetDX;
					T curDetUY = curDetDY;
					T curDetUZ = initDetUZ;

					T curDetX = initDetX * cosT - initDetY * sinT;
					T curDetY = initDetX * sinT + initDetY * cosT;

					T dirX = curDetX - cursourx;
					T dirY = curDetY - cursoury;

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
						std::swap(detPosLY, detPosRY);
					}
					//ŒÙÉèÏÂ±ßµÄÐ¡;
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
		} // End else
	} // End for
}



void DDM3D_EA_proj(std::vector<float>& proj, const std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<float>& angs)
{
	DDM3D_EA_proj_template<float>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ,
		detArc, detSizeH, detCntIdU, detCntIdV, XN, YN, ZN, DNU, DNV, PN, angs);
}


void DDM3D_EA_proj(std::vector<double>& proj, const std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<double>& angs)
{
	DDM3D_EA_proj_template<double>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ,
		detArc, detSizeH, detCntIdU, detCntIdV, XN, YN, ZN, DNU, DNV, PN, angs);
}





template<typename T>
void DDM3D_EA_helical_proj_template(std::vector<T>& proj, const std::vector<T>& vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detArc, const T detSizeH,
	const T detCntIdU, const T detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const T initZPos, const T pitch, const std::vector<T>& angs)
{
	const T dx = objSizeX / static_cast<T>(XN);
	const T dy = objSizeY / static_cast<T>(YN);
	const T dz = objSizeZ / static_cast<T>(ZN);
	const T dbeta = detArc / static_cast<T>(DNU);
	const T ddv = detSizeH / static_cast<T>(DNV);

	int angIdx(0);
#pragma omp parallel for
	for (angIdx = 0; angIdx < PN; ++angIdx)
	{
		T curang = angs[angIdx];
		T cosT = cos(curang);
		T sinT = sin(curang);

		T cursourx = S2O * cosT;
		T cursoury = S2O * sinT;
		T cursourz = initZPos + angIdx * pitch;

		if ((curang > PI * 0.25 && curang <= PI * 0.75) || (curang >= PI * 1.25 && curang < PI * 1.75))
		{
			for (int detIdV = 0; detIdV < DNV; ++detIdV)
			{
				for (int detIdU = 0; detIdU < DNU; ++detIdU)
				{
					T summ = 0;
					T beta = (detIdU - detCntIdU) * dbeta; //ÖÐŒä;
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
					T curDetDZ = initDetDZ + initZPos + angIdx * pitch;

					T curDetUX = curDetDX;
					T curDetUY = curDetDY;
					T curDetUZ = initDetUZ + initZPos + angIdx * pitch;

					T curDetX = initDetX * cosT - initDetY * sinT;
					T curDetY = initDetX * sinT + initDetY * cosT;

					T dirX = curDetX - cursourx;
					T dirY = curDetY - cursoury;

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
						std::swap(detPosLX, detPosRX);
					}
					//ŒÙÉèÏÂ±ßµÄÐ¡;
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
					T beta = (detIdU - detCntIdU) * dbeta; //ÖÐŒä;
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
					T curDetDZ = initDetDZ + initZPos + angIdx * pitch;//

					T curDetUX = curDetDX;
					T curDetUY = curDetDY;
					T curDetUZ = initDetUZ + initZPos + angIdx * pitch;

					T curDetX = initDetX * cosT - initDetY * sinT;
					T curDetY = initDetX * sinT + initDetY * cosT;

					T dirX = curDetX - cursourx;
					T dirY = curDetY - cursoury;

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
						std::swap(detPosLY, detPosRY);
					}
					//ŒÙÉèÏÂ±ßµÄÐ¡;
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





void DDM3D_EA_helical_proj(std::vector<float>& proj, const std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float initZPos, const float pitch, const std::vector<float>& angs)
{
	DDM3D_EA_helical_proj_template<float>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ,
		detArc, detSizeH, detCntIdU, detCntIdV,
		XN, YN, ZN, DNU, DNV, PN, initZPos, pitch, angs);
}


void DDM3D_EA_helical_proj(std::vector<double>& proj, const std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double initZPos, const double pitch, const std::vector<double>& angs)
{
	DDM3D_EA_helical_proj_template<double>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ,
		detArc, detSizeH, detCntIdU, detCntIdV,
		XN, YN, ZN, DNU, DNV, PN, initZPos, pitch, angs);
}





