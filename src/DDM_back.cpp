#include "DDM_back.h"

#include <omp.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "utilities.hpp"


template<typename T>
void DDM_ED_bproj_template(const std::vector<T>& proj, std::vector<T>& img,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY,
	const T detSize,
	const T detCntIdx,
	const int XN, const int YN, const int DN, const int PN,
	const std::vector<T>& angs)
{

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
		////ŒÆËãµ±Ç°ÏñËØËÄžöÖÐµãµÄ×ø±ê;
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

			if ((curang > PI * 0.25 && curang <= PI * 0.75) || (curang >= PI * 1.25 && curang < PI * 1.75)) //°ŽÕÕœÇ¶ÈÀŽŒÆËã;
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

			T objYaxisPosMin = initObjX1 * initObjY1 / (S2O - initObjX1) + initObjY1; //pixel¶ËµãÔÚYÖáÉÏµÄÍ¶Ó°;
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
				summ += proj[angIdx * DN + detIdx] * intersectLength<T>(objYaxisPosMin, objYaxisPosMax, minDetPos, maxDetPos) / (objYaxisLength * cosAng) * dx;
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
void DDM3D_ED_bproj_template(const std::vector<T>& proj, std::vector<T>& vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detSizeU, const T detSizeV,
	const T detCntIdU, const T detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<T>& angs)
{
	
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

			if ((curang > PI * 0.25 && curang <= PI * 0.75) || (curang >= PI * 1.25 && curang < PI * 1.75)) //°ŽÕÕœÇ¶ÈÀŽŒÆËã;
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

void DDM3D_ED_bproj(const std::vector<double>& proj, std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeU, const double detSizeV,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<double>& angs)
{
	DDM3D_ED_bproj_template<double>(proj, vol, S2O, O2D,
		objSizeX, objSizeY, objSizeZ,
		detSizeU, detSizeV, detCntIdU, detCntIdV,
		XN, YN, ZN, DNU, DNV, PN, angs);
}


void DDM3D_ED_bproj(const std::vector<float>& proj, std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeU, const float detSizeV,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<float>& angs)
{
	DDM3D_ED_bproj_template<float>(proj, vol, S2O, O2D,
		objSizeX, objSizeY, objSizeZ,
		detSizeU, detSizeV, detCntIdU, detCntIdV,
		XN, YN, ZN, DNU, DNV, PN, angs);
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
void DDM_EA_bproj_template(const std::vector<T>& proj, std::vector<T>& img,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY,
	const T detArc,
	const T detCntIdx,
	const int XN, const int YN, const int DN, const int PN,
	const std::vector<T>& angs)
{
	T cursourx(0), cursoury(0);
	const T dd = detArc / static_cast<T>(DN);
	const T dx = objSizeX / static_cast<T>(XN);
	const T dy = objSizeY / static_cast<T>(YN);

	const int OBJN = XN * YN;
	const T hfXN = XN * 0.5;
	const T hfYN = YN * 0.5;
	const T S2D = S2O + O2D;
#pragma omp parallel for
	for (int pixelIdx = 0; pixelIdx < OBJN; ++pixelIdx)
	{
		int jj = pixelIdx / XN;
		int ii = pixelIdx - jj * XN;

		//ŒÆËãµ±Ç°ÏñËØËÄžöÖÐµãµÄ×ø±ê;
		T summ = 0;
		for (int angIdx = 0; angIdx != PN; ++angIdx)
		{
			T curang = angs[angIdx];
			T cosT = cos(curang);
			T sinT = sin(curang);
			T lefPx(0);
			T lefPy(0);
			T rghPx(0);
			T rghPy(0);

			if ((curang > PI * 0.25 && curang <= PI * 0.75) || (curang >= PI * 1.25 && curang < PI * 1.75)) //°ŽÕÕœÇ¶ÈÀŽŒÆËã;
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

			T betaMin = initObjY1 / (S2O - initObjX1);
			T betaMax = initObjY2 / (S2O - initObjX2);
			if (betaMax < betaMin)
			{
				std::swap(betaMax, betaMin);
			}
			int minDetIdx = floor(atan(betaMax) / (-dd) + detCntIdx);
			int maxDetIdx = ceil(atan(betaMin) / (-dd) + detCntIdx);

			if (minDetIdx > DN)
			{
				continue;
			}
			if (maxDetIdx < 0)
			{
				continue;
			}

			T objYaxisPosMin = initObjX1 * initObjY1 / (S2O - initObjX1) + initObjY1; //pixel¶ËµãÔÚYÖáÉÏµÄÍ¶Ó°;
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
				T maxDetPos = -tan((detIdx - detCntIdx) * dd) * S2O;
				T minDetPos = -tan((detIdx - detCntIdx + 1.0) * dd) * S2O;
				if (maxDetPos < minDetPos)
				{
					std::swap(minDetPos, maxDetPos);
				}

				T cosAng = abs(cos((detIdx - detCntIdx + 0.5) * dd));
				summ += proj[angIdx * DN + detIdx] * intersectLength<T>(objYaxisPosMin, objYaxisPosMax, minDetPos, maxDetPos) / (objYaxisLength * cosAng) * dx;
			}
		}
		img[pixelIdx] += summ;
	}
}

void DDM_EA_bproj(const std::vector<double>& proj, std::vector<double>& img,
	const double S2O, const double O2D, const double objSizeX, const double objSizeY,
	const double detArc, const double detCntIdx, const int XN, const int YN, const int DN, const int PN,
	const std::vector<double>& angs)
{

	DDM_EA_bproj_template<double>(proj, img, S2O, O2D, objSizeX, objSizeY, detArc, detCntIdx,
		XN, YN, DN, PN, angs);
}
void DDM_EA_bproj(const std::vector<float>& proj, std::vector<float>& img,
	const float S2O, const float O2D, const float objSizeX, const float objSizeY,
	const float detArc, const float detCntIdx, const int XN, const int YN, const int DN, const int PN,
	const std::vector<float>& angs)
{
	DDM_EA_bproj_template<float>(proj, img, S2O, O2D, objSizeX, objSizeY, detArc, detCntIdx,
		XN, YN, DN, PN, angs);
}



template<typename T>
void DDM3D_EA_bproj_template(const std::vector<T>& proj, std::vector<T>& vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detArc, const T detSizeH,
	const T detCntIdU, const T detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<T>& angs)
{
	const T dbeta = detArc / static_cast<T>(DNU);
	const T ddv = detSizeH / static_cast<T>(DNV);
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

			T rghPx(0);
			T rghPy(0);

			T uppPx(0);
			T uppPy(0);
			T uppPz(0);

			T dowPx(0);
			T dowPy(0);
			T dowPz(0);

			if ((curang > PI * 0.25 && curang <= PI * 0.75) || (curang >= PI * 1.25 && curang < PI * 1.75))
			{
				lefPx = Px - 0.5 * dx;
				lefPy = Py;

				rghPx = Px + 0.5 * dx;
				rghPy = Py;

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

				rghPx = Px;
				rghPy = Py + 0.5 * dy;

				uppPx = Px;
				uppPy = Py;
				uppPz = Pz + 0.5 * dz;

				dowPx = Px;
				dowPy = Py;
				dowPz = Pz - 0.5 * dz;

			}

			T initObjLefx = lefPx * cosT + lefPy * sinT;
			T initObjlefy = -lefPx * sinT + lefPy * cosT;

			T initObjrghx = rghPx * cosT + rghPy * sinT;
			T initObjrghy = -rghPx * sinT + rghPy * cosT;

			T initObjuppx = uppPx * cosT + uppPy * sinT;

			T initObjuppz = uppPz;

			T initObjdowx = dowPx * cosT + dowPy * sinT;

			T initObjdowz = dowPz;

			T objYdetPosUMinbeta = atan(initObjLefx / (S2O - initObjlefy)); //Ò²ÐíÐèÒªÐÞžÄ;
			T objYdetPosUMaxbeta = atan(initObjrghx / (S2O - initObjrghy));
			if (objYdetPosUMinbeta > objYdetPosUMaxbeta)
			{
				std::swap(objYdetPosUMinbeta, objYdetPosUMaxbeta);
			}

			T objYdetPosVMin = initObjdowz * S2D / (S2O - initObjdowx);
			T objYdetPosVMax = initObjuppz * S2D / (S2O - initObjuppx);
			if (objYdetPosVMin > objYdetPosVMax)
			{
				std::swap(objYdetPosVMin, objYdetPosVMax);
			}

			int minDetUIdx = floor(objYdetPosUMinbeta / dbeta + detCntIdU) - 1;
			int maxDetUIdx = ceil(objYdetPosUMaxbeta / dbeta + detCntIdU) + 1;

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

			T objYOZLengthbeta = abs(objYdetPosUMaxbeta - objYdetPosUMinbeta);//žÄ;
			T objYOZHeight = abs(objYdetPosVMax - objYdetPosVMin);

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
				T minDetUPosbeta = (detIdU - detCntIdU - 0.5) * dbeta;
				T maxDetUPosbeta = (detIdU - detCntIdU + 0.5) * dbeta;

				T ll = intersectLength<T>(objYdetPosUMinbeta, objYdetPosUMaxbeta, minDetUPosbeta, maxDetUPosbeta);
				if (ll > 0)
				{
					for (int detIdV = minDetVIdx; detIdV < maxDetVIdx; ++detIdV)
					{
						T minDetVPos = (detIdV - detCntIdV - 0.5) * ddv;
						T maxDetVPos = (detIdV - detCntIdV + 0.5) * ddv;

						T DU = abs(sin((detIdU - detCntIdU) * dbeta) * S2D);
						T DV = (detIdV - detCntIdV) * ddv;
						T cosAlphacosGamma = S2D / sqrt(DU * DU + DV * DV + S2D * S2D);
						T mm = intersectLength<T>(objYdetPosVMin, objYdetPosVMax, minDetVPos, maxDetVPos);
						if (mm > 0)
						{
							summ += (proj[(angIdx * DNV + detIdV) * DNU + detIdU] * ll * mm / (objYOZLengthbeta * objYOZHeight * cosAlphacosGamma) * dx);

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


void DDM3D_EA_bproj(const std::vector<float>& proj, std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<float>& angs)
{
	DDM3D_EA_bproj_template<float>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ,
		detArc, detSizeH, detCntIdU, detCntIdV, XN, YN, ZN, DNU, DNV, PN, angs);
}

void DDM3D_EA_bproj(const std::vector<double>& proj, std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<double>& angs)
{
	DDM3D_EA_bproj_template<double>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ,
		detArc, detSizeH, detCntIdU, detCntIdV, XN, YN, ZN, DNU, DNV, PN, angs);
}




template<typename T>
void DDM3D_EA_helical_bproj_template(const std::vector<T>& proj, std::vector<T>& vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detArc, const T detSizeH,
	const T detCntIdU, const T detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const T initZPos, const T pitch,
	const std::vector<T>& angs)
{
	const T dbeta = detArc / static_cast<T>(DNU);
	const T ddv = detSizeH / static_cast<T>(DNV);
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

			T rghPx(0);
			T rghPy(0);

			T uppPx(0);
			T uppPy(0);
			T uppPz(0);

			T dowPx(0);
			T dowPy(0);
			T dowPz(0);

			if ((curang > PI * 0.25 && curang <= PI * 0.75) || (curang >= PI * 1.25 && curang < PI * 1.75))
			{
				lefPx = Px - 0.5 * dx;
				lefPy = Py;

				rghPx = Px + 0.5 * dx;
				rghPy = Py;

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

				rghPx = Px;
				rghPy = Py + 0.5 * dy;

				uppPx = Px;
				uppPy = Py;
				uppPz = Pz + 0.5 * dz;

				dowPx = Px;
				dowPy = Py;
				dowPz = Pz - 0.5 * dz;

			}

			T initObjLefx = lefPx * cosT + lefPy * sinT;
			T initObjlefy = -lefPx * sinT + lefPy * cosT;

			T initObjrghx = rghPx * cosT + rghPy * sinT;
			T initObjrghy = -rghPx * sinT + rghPy * cosT;

			T initObjuppx = uppPx * cosT + uppPy * sinT;

			T initObjuppz = uppPz - (initZPos + angIdx * pitch);

			T initObjdowx = dowPx * cosT + dowPy * sinT;

			T initObjdowz = dowPz - (initZPos + angIdx * pitch);

			T objYdetPosUMinbeta = atan(initObjLefx / (S2O - initObjlefy)); //Ò²ÐíÐèÒªÐÞžÄ;
			T objYdetPosUMaxbeta = atan(initObjrghx / (S2O - initObjrghy));
			if (objYdetPosUMinbeta > objYdetPosUMaxbeta)
			{
				std::swap(objYdetPosUMinbeta, objYdetPosUMaxbeta);
			}

			T objYdetPosVMin = initObjdowz * S2D / (S2O - initObjdowx);
			T objYdetPosVMax = initObjuppz * S2D / (S2O - initObjuppx);
			if (objYdetPosVMin > objYdetPosVMax)
			{
				std::swap(objYdetPosVMin, objYdetPosVMax);
			}

			int minDetUIdx = floor(objYdetPosUMinbeta / dbeta + detCntIdU) - 1;
			int maxDetUIdx = ceil(objYdetPosUMaxbeta / dbeta + detCntIdU) + 1;

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

			T objYOZLengthbeta = abs(objYdetPosUMaxbeta - objYdetPosUMinbeta);//žÄ;
			T objYOZHeight = abs(objYdetPosVMax - objYdetPosVMin);

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
				T minDetUPosbeta = (detIdU - detCntIdU - 0.5) * dbeta;
				T maxDetUPosbeta = (detIdU - detCntIdU + 0.5) * dbeta;

				T ll = intersectLength<T>(objYdetPosUMinbeta, objYdetPosUMaxbeta, minDetUPosbeta, maxDetUPosbeta);
				if (ll > 0)
				{
					for (int detIdV = minDetVIdx; detIdV < maxDetVIdx; ++detIdV)
					{
						T minDetVPos = (detIdV - detCntIdV - 0.5) * ddv;
						T maxDetVPos = (detIdV - detCntIdV + 0.5) * ddv;

						T DU = abs(sin((detIdU - detCntIdU) * dbeta) * S2D);
						T DV = (detIdV - detCntIdV) * ddv;
						T cosAlphacosGamma = S2D / sqrt(DU * DU + DV * DV + S2D * S2D);
						T mm = intersectLength<T>(objYdetPosVMin, objYdetPosVMax, minDetVPos, maxDetVPos);
						if (mm > 0)
						{
							summ += (proj[(angIdx * DNV + detIdV) * DNU + detIdU] * ll * mm / (objYOZLengthbeta * objYOZHeight * cosAlphacosGamma) * dx);

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

void DDM3D_EA_helical_bproj(const std::vector<float>& proj, std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float initZPos, const float pitch,
	const std::vector<float>& angs)
{
	DDM3D_EA_helical_bproj_template<float>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ,
		detArc, detSizeH, detCntIdU, detCntIdV, XN, YN, ZN, DNU, DNV, PN,
		initZPos, pitch, angs);
}

void DDM3D_EA_helical_bproj(const std::vector<double>& proj, std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double initZPos, const double pitch,
	const std::vector<double>& angs)
{
	DDM3D_EA_helical_bproj_template<double>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ,
		detArc, detSizeH, detCntIdU, detCntIdV, XN, YN, ZN, DNU, DNV, PN,
		initZPos, pitch, angs);
}

