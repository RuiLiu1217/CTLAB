//#include <cstdlib>
//#include <cmath>
//#include <utility>
//#include <vector_types.h>
//#include <vector_functions.h>
//#include <fstream>
//#include <iostream>
//#include <string>
//#include <sstream>
//#include <vector>
//#include <functional>
//#include <algorithm>
//#include <omp.h>
//#include "ConstantValues.h"
//#include "genSysMatrix.hpp"
//#include "utilities.hpp"
//#include "FanEAGeo.h"
//#include "FanEDGeo.h"
//namespace genSysMatrix
//{
//	static const float EPSILON = 1.0E-9;
//
//	template<typename T>
//	inline bool IS_ZERO(const T& x)
//	{
//		return ((x < EPSILON) && (x > -EPSILON));
//	}
//
//
//
//	float2 rotation(const float2& p, const float& cosT, const float& sinT)
//	{
//		float2 curP;
//		curP.x = p.x * cosT - p.y * sinT;
//		curP.y = p.x * sinT + p.y * cosT;
//		return curP;
//	}
//
//	double2 rotation(const double2& p, const double& cosT, const double& sinT)
//	{
//		double2 curP;
//		curP.x = p.x * cosT - p.y * sinT;
//		curP.y = p.x * sinT + p.y * cosT;
//		return curP;
//	}
//
//	struct Ray2Dd
//	{
//	public:
//		double2 o;
//		double2 d;
//	};
//
//	void pushMatrix(
//		float startX, float startY,
//		float endX, float endY,
//		float bx, float by,
//		float dx, float dy,
//		const int objResoLen,
//		const int objResoWid,
//		int& rowidx,
//		std::vector<int>& rowIdx,
//		std::vector<int>& colIdx,
//		std::vector<float>& wgt)
//	{
//		const float dirX(endX - startX);
//		const float dirY(endY - startY);
//		const float lengthSq = dirX * dirX + dirY * dirY;
//		const float dconv = sqrt(lengthSq);
//		int imin, imax, jmin, jmax;
//		const float alphaxmin = fminf(dev_alpha_IFun(bx, dx, startX, endX, 0), dev_alpha_IFun(bx, dx, startX, endX, objResoLen));
//		const float alphaxmax = fmaxf(dev_alpha_IFun(bx, dx, startX, endX, 0), dev_alpha_IFun(bx, dx, startX, endX, objResoLen));
//		const float alphaymin = fminf(dev_alpha_IFun(by, dy, startY, endY, 0), dev_alpha_IFun(by, dy, startY, endY, objResoWid));
//		const float alphaymax = fmaxf(dev_alpha_IFun(by, dy, startY, endY, 0), dev_alpha_IFun(by, dy, startY, endY, objResoWid));
//
//		const float alphaMIN = fmaxf(alphaxmin, alphaymin);
//		const float alphaMAX = fminf(alphaxmax, alphaymax);
//		dev_minmaxIdxFun(startX, endX, bx, dx, alphaMIN, alphaMAX, alphaxmin, alphaxmax, objResoLen + 1, &imin, &imax);
//		dev_minmaxIdxFun(startY, endY, by, dy, alphaMIN, alphaMAX, alphaymin, alphaymax, objResoWid + 1, &jmin, &jmax);
//
//		float alphaX = (startX < endX) ? dev_alpha_IFun(bx, dx, startX, endX, imin) : dev_alpha_IFun(bx, dx, startX, endX, imax);
//		float alphaY = (startY < endY) ? dev_alpha_IFun(by, dy, startY, endY, jmin) : dev_alpha_IFun(by, dy, startY, endY, jmax);
//
//		int Np = static_cast<int>(fabsf(imax - imin + 1.0f) + fabsf(jmax - jmin + 1.0f) + 4.0f);
//		const float alphaxu = dev_alphaU_Fun(dx, startX, endX);
//		const float alphayu = dev_alphaU_Fun(dy, startY, endY);
//
//		float alphaC = alphaMIN;
//
//		int i = static_cast<int>(dev_varphiFun(alphaMIN* 1.00003f, bx, dx, startX, endX));
//		int j = static_cast<int>(dev_varphiFun(alphaMIN* 1.00003f, by, dy, startY, endY));
//
//		const int iuu = dev_iu_Fun(startX, endX);
//		const int juu = dev_iu_Fun(startY, endY);
//
//		float d12(0.0f);
//		float weight(0.0f);
//		unsigned int repIdx(0);
//		unsigned int colidx(0);
//		while (repIdx != Np)
//		{
//			if (i < 0 || i >= objResoLen || j < 0 || j >= objResoWid)
//			{
//				break;
//			}
//			if (alphaX <= alphaY)
//			{
//				colidx = j * objResoLen + i;
//				weight = (alphaX - alphaC) * dconv;
//
//				wgt.push_back(weight);
//				rowIdx.push_back(rowidx);
//				colIdx.push_back(colidx);
//
//				d12 += weight;
//				i += iuu;
//				alphaC = alphaX;
//				alphaX += alphaxu;
//			}
//			else
//			{
//				colidx = j * objResoLen + i;
//				weight = (alphaY - alphaC) * dconv;
//
//				wgt.push_back(weight);
//				rowIdx.push_back(rowidx);
//				colIdx.push_back(colidx);
//
//				d12 += weight;
//				j += juu;
//
//				alphaC = alphaY;
//				alphaY += alphayu;
//			}
//			++repIdx;
//		}
//	}
//
//
//
//
//	void pushMatrix(
//		double startX, double startY,
//		double endX, double endY,
//		double bx, double by,
//		double dx, double dy,
//		const int objResoLen,
//		const int objResoWid,
//		int& rowidx,
//		std::vector<int>& rowIdx,
//		std::vector<int>& colIdx,
//		std::vector<double>& wgt)
//	{
//		const double dirX(endX - startX);
//		const double dirY(endY - startY);
//		const double lengthSq = dirX * dirX + dirY * dirY;
//		const double dconv = sqrt(lengthSq);
//		int imin, imax, jmin, jmax;
//		const double alphaxmin = fminf(dev_alpha_IFun(bx, dx, startX, endX, 0), dev_alpha_IFun(bx, dx, startX, endX, objResoLen));
//		const double alphaxmax = fmaxf(dev_alpha_IFun(bx, dx, startX, endX, 0), dev_alpha_IFun(bx, dx, startX, endX, objResoLen));
//		const double alphaymin = fminf(dev_alpha_IFun(by, dy, startY, endY, 0), dev_alpha_IFun(by, dy, startY, endY, objResoWid));
//		const double alphaymax = fmaxf(dev_alpha_IFun(by, dy, startY, endY, 0), dev_alpha_IFun(by, dy, startY, endY, objResoWid));
//
//		const double alphaMIN = fmaxf(alphaxmin, alphaymin);
//		const double alphaMAX = fminf(alphaxmax, alphaymax);
//		dev_minmaxIdxFun(startX, endX, bx, dx, alphaMIN, alphaMAX, alphaxmin, alphaxmax, objResoLen + 1, &imin, &imax);
//		dev_minmaxIdxFun(startY, endY, by, dy, alphaMIN, alphaMAX, alphaymin, alphaymax, objResoWid + 1, &jmin, &jmax);
//
//		double alphaX = (startX < endX) ? dev_alpha_IFun(bx, dx, startX, endX, imin) : dev_alpha_IFun(bx, dx, startX, endX, imax);
//		double alphaY = (startY < endY) ? dev_alpha_IFun(by, dy, startY, endY, jmin) : dev_alpha_IFun(by, dy, startY, endY, jmax);
//
//		int Np = static_cast<int>(abs(imax - imin + 1.0f) + abs(jmax - jmin + 1.0f) + 4.0f);
//		const double alphaxu = dev_alphaU_Fun(dx, startX, endX);
//		const double alphayu = dev_alphaU_Fun(dy, startY, endY);
//
//		double alphaC = alphaMIN;
//		double talpha = std::min(alphaX, alphaY);
//
//		int i = static_cast<int>(dev_varphiFun((talpha + alphaMIN) * 0.5, bx, dx, startX, endX));
//		int j = static_cast<int>(dev_varphiFun((talpha + alphaMIN) * 0.5, by, dy, startY, endY));
//
//		//int i = static_cast<int>(dev_varphiFun(alphaMIN* 1.00003f, bx, dx, startX, endX));
//		//int j = static_cast<int>(dev_varphiFun(alphaMIN* 1.00003f, by, dy, startY, endY));
//
//		const int iuu = dev_iu_Fun(startX, endX);
//		const int juu = dev_iu_Fun(startY, endY);
//
//		double d12(0.0f);
//		double weight(0.0f);
//		unsigned int repIdx(0);
//		unsigned int colidx(0);
//		while (repIdx != Np)
//		{
//			if (i < 0 || i >= objResoLen || j < 0 || j >= objResoWid)
//			{
//				break;
//			}
//			if (alphaX <= alphaY)
//			{
//				colidx = j * objResoLen + i;
//				weight = (alphaX - alphaC) * dconv;
//
//				wgt.push_back(weight);
//				rowIdx.push_back(rowidx);
//				colIdx.push_back(colidx);
//
//				d12 += weight;
//				i += iuu;
//				alphaC = alphaX;
//				alphaX += alphaxu;
//			}
//			else
//			{
//				colidx = j * objResoLen + i;
//				weight = (alphaY - alphaC) * dconv;
//
//				wgt.push_back(weight);
//				rowIdx.push_back(rowidx);
//				colIdx.push_back(colidx);
//
//				d12 += weight;
//				j += juu;
//
//				alphaC = alphaY;
//				alphaY += alphayu;
//			}
//			++repIdx;
//		}
//	}
//
//
//
//
//
//
//	void genProj_SIDDON(
//		std::vector<int>& rowIdx,
//		std::vector<int>& colIdx,
//		std::vector<float>& weight,
//		const FanEAGeo& FanGeo,
//		const Image& Img,
//		const unsigned int& sampNum)
//	{
//		float2 MINO = make_float2(
//			-Img.m_Size.x / 2.0f + Img.m_Bias.x,
//			-Img.m_Size.y / 2.0f + Img.m_Bias.y);
//
//		float curAng = 0;
//		float cosT = 0;
//		float sinT = 0;
//		Ray2D ray;
//
//		//unsigned int detId;
//		float ang(0); //bias angle from the center of the Fan Beams
//
//		float2 curDetPos; //the current detector element position;
//		//float totWeight(0);
//
//		float smallDetStep = FanGeo.m_DetStp / sampNum; //\CF²\C9\D1\F9\BA\F3\B5\C4detector\B2\BD\B3\A4;
//		float cntSmallDet = sampNum * 0.5f;
//		float realAng = 0;
//		unsigned int angIdx = 0;
//		unsigned int detIdx = 0;
//		unsigned int subDetIdx = 0;
//
//		int rowidx = 0;
//		for (angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
//		{
//			//Current rotation angle;
//			curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//			for (detIdx = 0; detIdx != FanGeo.m_DetN; ++detIdx)
//			{
//
//				rowidx = angIdx * FanGeo.m_DetN + detIdx;
//
//				//Calculate current Angle
//				ang = ((float) detIdx - FanGeo.m_DetCntIdx + 0.5f) * FanGeo.m_DetStp;
//
//				for (subDetIdx = 0; subDetIdx != sampNum; ++subDetIdx)
//				{
//					//correct ang
//					realAng = ang + (static_cast<float>(subDetIdx) -cntSmallDet + 0.5f) *smallDetStep;
//					// current detector element position;
//					curDetPos = rotation(make_float2(sinf(realAng) * FanGeo.m_S2D, -cosf(realAng) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//
//					// X-ray direction
//					ray.d = normalize(curDetPos - ray.o);
//
//					pushMatrix(ray.o.x, ray.o.y,
//						curDetPos.x, curDetPos.y,
//						MINO.x, MINO.y,
//						Img.m_Step.x, Img.m_Step.y,
//						Img.m_Reso.x, Img.m_Reso.y,
//						rowidx, rowIdx, colIdx, weight);
//				}
//			}
//		}
//
//	}
//
//	void genProj_SIDDON(
//		std::vector<int>& rowIdx,
//		std::vector<int>& colIdx,
//		std::vector<double>& weight,
//		const std::vector<double> angs,
//		const double S2O,
//		const double O2D,
//		const double detCntIdx,
//		const double detStp,
//		const double dx,
//		const double dy,
//		const int prjNum,
//		const int detNum,
//		const int XN,
//		const int YN)
//	{
//		int ridx = 0, cidx = 0;
//		bool hit = false;
//		double tnear, tfar;
//		double curAng, cosT, sinT, curSourX, curSourY;
//		double curDetX, curDetY;
//		double boxminx, boxminy, boxmaxx, boxmaxy;
//		const double minX = -static_cast<double>(XN) / 2.0 * dx;
//		const double minY = -static_cast<double>(YN) / 2.0 * dy;
//		double detang;
//		const double S2D = S2O + O2D;
//		double XX, YY;
//		for (unsigned int prjIdx = 0; prjIdx != prjNum; ++prjIdx)
//		{
//			curAng = angs[prjIdx];
//			cosT = cos(curAng);
//			sinT = sin(curAng);
//			curSourX = -sinT * S2O;
//			curSourY = cosT * S2O;
//
//			for (unsigned int detIdx = 0; detIdx != detNum; ++detIdx)
//			{
//
//				//Calculate current detector position
//				//(detIdx - detCntIdx) * detStp
//				detang = ((double) detIdx - detCntIdx) * detStp;
//
//				// current detector element position;
//				XX = sin(detang) * S2D;
//				YY = -cos(detang) * S2D + S2O;
//				curDetX = XX * cosT - YY * sinT;
//				curDetY = XX * sinT + YY * cosT;
//
//				ridx = prjIdx * detNum + detIdx;
//				for (unsigned int yIdx = 0; yIdx != YN; ++yIdx)
//				{
//					for (unsigned int xIdx = 0; xIdx != XN; ++xIdx)
//					{
//						boxminx = xIdx * dx + minX;
//						boxminy = yIdx * dy * minY;
//						boxmaxx = boxminx + dx;
//						boxmaxy = boxminy + dy;
//
//						cidx = yIdx * XN + xIdx;
//
//						tnear = 0;
//						tfar = 0;
//						hit = intersectBox(curSourX, curSourY, curDetX, curDetY,
//							boxminx, boxminy, boxmaxx, boxmaxy, &tnear, &tfar);
//						if (hit)
//						{
//							rowIdx.push_back(ridx);
//							colIdx.push_back(cidx);
//							weight.push_back(tfar - tnear);
//						}
//					}
//				}
//			}
//			std::cout << "current projection angle index = " << prjIdx << "\n";
//		}
//	}
//
//	// generate the COO
//	void pushValuesCOO(
//		const int sarIdx,
//		const int endIdx,
//		const double S2O,
//		const double S2D,
//		const int DN,
//		const int XN,
//		const int YN,
//		const double detCntIdx,
//		const double objCntIdxy,
//		const double objCntIdxx,
//		const double ObjStpx,
//		const double ObjStpy,
//		const double detStp,
//		const std::vector<double>& angs,
//		std::vector<int>& rowIdx,
//		std::vector<int>& colIdx,
//		std::vector<double>& weight)
//	{
//		double cosT, sinT, sourx, soury, beta, sinBeta, cosBeta,
//			initDetX, initDetY, curDetX, curDetY, dirx, diry,
//			rleng, boxminx, boxminy, boxmaxx, boxmaxy, tnear, tfar;
//		bool intersected;
//
//		for (int angIdx = sarIdx; angIdx != endIdx; ++angIdx)
//		{
//			cosT = cos(angs[angIdx]);
//			sinT = sin(angs[angIdx]);
//			sourx = -S2O * sinT;
//			soury = S2O * cosT;
//
//			//std::cout << angIdx << std::endl;
//			for (unsigned int detIdx = 0; detIdx != DN; ++detIdx)
//			{
//				beta = (detIdx - detCntIdx) * detStp;
//				sinBeta = sin(beta);
//				cosBeta = cos(beta);
//				initDetX = S2D * sinBeta;
//				initDetY = -S2D * cosBeta + S2O;
//				curDetX = initDetX * cosT - initDetY * sinT;
//				curDetY = initDetX * sinT + initDetY * cosT;
//				dirx = curDetX - sourx;
//				diry = curDetY - soury;
//				rleng = sqrt(dirx * dirx + diry * diry);
//				dirx /= rleng;
//				diry /= rleng;
//
//				for (unsigned int yIdx = 0; yIdx != YN; ++yIdx)
//				{
//					boxminy = (yIdx - objCntIdxy - 0.5) * ObjStpy;
//					boxmaxy = (yIdx - objCntIdxy + 0.5) * ObjStpy;
//					for (unsigned int xIdx = 0; xIdx != XN; ++xIdx)
//					{
//						boxminx = (xIdx - objCntIdxx - 0.5) * ObjStpx;
//						boxmaxx = (xIdx - objCntIdxx + 0.5) * ObjStpx;
//
//						intersected = intersectBox<double>(sourx, soury,
//							dirx, diry, boxminx, boxminy,
//							boxmaxx, boxmaxy, &tnear, &tfar);
//						if (intersected)
//						{
//							rowIdx.push_back(angIdx * DN + detIdx);
//							colIdx.push_back(yIdx * XN + xIdx);
//							weight.push_back(tfar - tnear);
//						}
//
//					}
//				}
//			}
//		}
//	}
//
//
//	void pushValuesCOO_upSample4(
//		const int sarIdx,
//		const int endIdx,
//		const double S2O,
//		const double S2D,
//		const int DN,
//		const int XN,
//		const int YN,
//		const double detCntIdx,
//		const double objCntIdxy,
//		const double objCntIdxx,
//		const double ObjStpx,
//		const double ObjStpy,
//		const double detStp,
//		const std::vector<double>& angs,
//		std::vector<int>& rowIdx,
//		std::vector<int>& colIdx,
//		std::vector<double>& weight)
//	{
//		double cosT, sinT, sourx, soury, beta, sinBeta, cosBeta,
//			initDetX, initDetY, curDetX, curDetY, dirx, diry,
//			rleng, boxminx, boxminy, boxmaxx, boxmaxy, tnear, tfar;
//		bool intersected;
//		double updetStp = detStp / 8.0;
//		double wweight = 0;
//
//		for (unsigned int yIdx = 0; yIdx != YN; ++yIdx)
//		{
//			boxminy = (yIdx - objCntIdxy - 0.5) * ObjStpy;
//			boxmaxy = (yIdx - objCntIdxy + 0.5) * ObjStpy;
//			for (unsigned int xIdx = 0; xIdx != XN; ++xIdx)
//			{
//				boxminx = (xIdx - objCntIdxx - 0.5) * ObjStpx;
//				boxmaxx = (xIdx - objCntIdxx + 0.5) * ObjStpx;
//				for (int angIdx = sarIdx; angIdx != endIdx; ++angIdx) // Each angle;
//				{
//					cosT = cos(angs[angIdx]);
//					sinT = sin(angs[angIdx]);
//					sourx = -S2O * sinT;
//					soury = S2O * cosT;
//					for (unsigned int detIdx = 0; detIdx != DN; ++detIdx) // Each detector cell;
//					{
//						//intersected = false;
//						wweight = 0;
//						for (unsigned int subdetIdx = 0; subdetIdx != 8; ++subdetIdx)
//						{
//							beta = (detIdx - detCntIdx) * detStp + (subdetIdx - 3.5) * updetStp;
//							sinBeta = sin(beta);
//							cosBeta = cos(beta);
//							initDetX = S2D * sinBeta;
//							initDetY = -S2D * cosBeta + S2O;
//							curDetX = initDetX * cosT - initDetY * sinT;
//							curDetY = initDetX * sinT + initDetY * cosT;
//							dirx = curDetX - sourx;
//							diry = curDetY - soury;
//							rleng = sqrt(dirx * dirx + diry * diry);
//							dirx /= rleng;
//							diry /= rleng;
//							intersected = intersectBox<double>(sourx, soury,
//								dirx, diry, boxminx, boxminy,
//								boxmaxx, boxmaxy, &tnear, &tfar);
//							if (intersected)
//							{
//								wweight += (tfar - tnear);
//							}
//						}
//						if (wweight != 0)
//						{
//							rowIdx.push_back(angIdx * DN + detIdx);
//							colIdx.push_back(yIdx * XN + xIdx);
//							weight.push_back(wweight * 0.125);
//						}
//					}
//				}
//			}
//		}
//	}
//
//
//
//
//
//
//
//
//
//
//
//	void genProj_SIDDON(std::vector<int>& rowIdx,
//		std::vector<int>& colIdx,
//		std::vector<double>& weight,
//		const double S2O, const double O2D,
//		const double detCntIdx, const double detStp,
//		const int DN,
//		const double ObjStpx, const double ObjStpy,
//		const int XN, const int YN,
//		const double objCntIdxx,
//		const double objCntIdxy,
//		const std::vector<double> angs)
//	{
//		const int PN = angs.size();
//		double cosT, sinT;
//		double sourx, soury;
//		double beta = 0;
//		double cosBeta, sinBeta;
//		double initDetX, initDetY;
//		const double S2D = S2O + O2D;
//		double curDetX, curDetY;
//		double boxminy, boxmaxy, boxminx, boxmaxx;
//		double dirx, diry, rleng;
//		bool intersected;
//		double tnear, tfar;
//
//		for (int angIdx = 0; angIdx != PN; ++angIdx)
//		{
//			cosT = cos(angs[angIdx]);
//			sinT = sin(angs[angIdx]);
//			sourx = -S2O * sinT;
//			soury = S2O * cosT;
//
//			std::cout << angIdx << std::endl;
//			for (unsigned int detIdx = 0; detIdx != DN; ++detIdx)
//			{
//				beta = (detIdx - detCntIdx) * detStp;
//				sinBeta = sin(beta);
//				cosBeta = cos(beta);
//				initDetX = S2D * sinBeta;
//				initDetY = -S2D * cosBeta + S2O;
//				curDetX = initDetX * cosT - initDetY * sinT;
//				curDetY = initDetX * sinT + initDetY * cosT;
//				dirx = curDetX - sourx;
//				diry = curDetY - soury;
//				rleng = sqrt(dirx * dirx + diry * diry);
//				dirx /= rleng;
//				diry /= rleng;
//
//				for (unsigned int yIdx = 0; yIdx != YN; ++yIdx)
//				{
//					boxminy = (yIdx - objCntIdxy - 0.5) * ObjStpy;
//					boxmaxy = (yIdx - objCntIdxy + 0.5) * ObjStpy;
//					for (unsigned int xIdx = 0; xIdx != XN; ++xIdx)
//					{
//						boxminx = (xIdx - objCntIdxx - 0.5) * ObjStpx;
//						boxmaxx = (xIdx - objCntIdxx + 0.5) * ObjStpx;
//
//						intersected = intersectBox<double>(sourx, soury,
//							dirx, diry, boxminx, boxminy,
//							boxmaxx, boxmaxy, &tnear, &tfar);
//						if (intersected)
//						{
//							rowIdx.push_back(angIdx * DN + detIdx);
//							colIdx.push_back(yIdx * XN + xIdx);
//							weight.push_back(tfar - tnear);
//						}
//
//					}
//				}
//			}
//		}
//
//
//
//	}
//
//
//
//
//	void genProj_SIDDON_openMP(
//		std::vector<int>& rowIdx,
//		std::vector<int>& colIdx,
//		std::vector<double>& weight,
//		const double S2O, const double O2D,
//		const double detCntIdx, const double detStp,
//		const int DN,
//		const double ObjStpx, const double ObjStpy,
//		const int XN, const int YN,
//		const double objCntIdxx,
//		const double objCntIdxy,
//		const std::vector<double> angs)
//	{
//		const int PN = angs.size();
//		//double cosT, sinT;
//		//double sourx, soury;
//		//double beta = 0;
//		//double cosBeta, sinBeta;
//		//double initDetX, initDetY;
//		const double S2D = S2O + O2D;
//		//double curDetX, curDetY;
//		//double boxminy, boxmaxy, boxminx, boxmaxx;
//		//double dirx, diry, rleng;
//		//bool intersected;
//		//double tnear, tfar;
//
//
//		int threadsNum = omp_get_num_procs();
//		int vPN_perRang = PN / threadsNum;
//		int *rang = new int[threadsNum + 1];
//		rang[0] = 0;
//		for (int i = 1; i <= threadsNum; ++i)
//		{
//			rang[i] = vPN_perRang * i;
//		}
//		rang[threadsNum] = PN;
//
//		std::vector<int> *srowIdx = new std::vector<int>[threadsNum];
//		std::vector<int> *scolIdx = new std::vector<int>[threadsNum];
//		std::vector<double> *sweight = new std::vector<double>[threadsNum];
//
//
//		//#pragma omp parallel for
//		//	for (int threadidx = 0; threadidx < threadsNum; ++threadidx)
//		//	{
//		//		pushValuesCOO(rang[threadidx], rang[threadidx + 1], S2O, S2D,
//		//			DN, XN, YN, detCntIdx, objCntIdxy, objCntIdxx, ObjStpx, ObjStpy,
//		//			detStp, angs, srowIdx[threadidx], scolIdx[threadidx], sweight[threadidx]);
//		//	}
//
//#pragma omp parallel for
//		for (int threadidx = 0; threadidx < threadsNum; ++threadidx)
//		{
//			pushValuesCOO_upSample4(rang[threadidx], rang[threadidx + 1], S2O, S2D,
//				DN, XN, YN, detCntIdx, objCntIdxy, objCntIdxx, ObjStpx, ObjStpy,
//				detStp, angs, srowIdx[threadidx], scolIdx[threadidx], sweight[threadidx]);
//		}
//
//
//		//Connect them together
//		rowIdx.clear();
//		colIdx.clear();
//		weight.clear();
//
//		for (int kk = 0; kk != threadsNum; ++kk)
//		{
//			rowIdx.insert(rowIdx.end(), srowIdx[kk].begin(), srowIdx[kk].end());
//			colIdx.insert(colIdx.end(), scolIdx[kk].begin(), scolIdx[kk].end());
//			weight.insert(weight.end(), sweight[kk].begin(), sweight[kk].end());
//		}
//	}
//
//
//	void genProj_SIDDON(
//		std::vector<int>& rowIdx,
//		std::vector<int>& colIdx,
//		std::vector<double>& weight,
//		const FanEAGeo& FanGeo,
//		const Image& Img,
//		const unsigned int& sampNum)
//	{
//		double2 MINO = make_double2(
//			-Img.m_Size.x / 2.0f + Img.m_Bias.x,
//			-Img.m_Size.y / 2.0f + Img.m_Bias.y);
//
//		double curAng = 0;
//		double cosT = 0;
//		double sinT = 0;
//		Ray2Dd ray;
//
//		//unsigned int detId;
//		double ang(0); //bias angle from the center of the Fan Beams
//
//		double2 curDetPos; //the current detector element position;
//		//double totWeight(0);
//
//		double smallDetStep = FanGeo.m_DetStp / sampNum; //\CF²\C9\D1\F9\BA\F3\B5\C4detector\B2\BD\B3\A4;
//		double cntSmallDet = sampNum * 0.5f;
//		double realAng = 0;
//		unsigned int angIdx = 0;
//		unsigned int detIdx = 0;
//		unsigned int subDetIdx = 0;
//
//		int rowidx = 0;
//		for (angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
//		{
//			//Current rotation angle;
//			curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//			cosT = cos(curAng);
//			sinT = sin(curAng);
//			ray.o = rotation(make_double2(0, FanGeo.m_S2O), cosT, sinT);
//
//			for (detIdx = 0; detIdx != FanGeo.m_DetN; ++detIdx)
//			{
//
//				rowidx = angIdx * FanGeo.m_DetN + detIdx;
//
//				//Calculate current Angle
//				ang = ((double) detIdx - FanGeo.m_DetCntIdx + 0.5f) * FanGeo.m_DetStp;
//
//				for (subDetIdx = 0; subDetIdx != sampNum; ++subDetIdx)
//				{
//					//correct ang
//					realAng = ang + (static_cast<double>(subDetIdx) -cntSmallDet + 0.5f) *smallDetStep;
//					// current detector element position;
//					curDetPos = rotation(make_double2(sin(realAng) * FanGeo.m_S2D, -cos(realAng) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//
//					// X-ray direction
//					ray.d = normalize(curDetPos - ray.o);
//
//					pushMatrix(ray.o.x, ray.o.y,
//						curDetPos.x, curDetPos.y,
//						MINO.x, MINO.y,
//						Img.m_Step.x, Img.m_Step.y,
//						Img.m_Reso.x, Img.m_Reso.y,
//						rowidx, rowIdx, colIdx, weight);
//				}
//			}
//		}
//
//	}
//
//
//
//	
//	void genMatrix_EA_AIM()
//	{
//		FanEAGeo FanGeo;
//		Image Img;
//		std::ifstream fin("configurationFile.txt", std::ios::in);
//		if (!fin.is_open())
//		{
//			std::cout << "Cannot Find configuration file\n";
//			std::cout << "Would you like to use the default configuration?[y/n]";
//			char yes;
//			std::cin >> yes;
//			if (yes == 'y' || yes == 'Y')
//			{
//				FanGeo.m_S2O = 5.385200195312500e+02;
//				FanGeo.m_O2D = 4.082259521484375e+02;
//				FanGeo.m_ViwN = 360;
//				FanGeo.m_ViwBeg = 0;
//				FanGeo.m_ViwStp = 3.14159265358979323846264*2.0 / 359.0;
//				FanGeo.m_DetArc = 0.95928517242269f;
//				FanGeo.m_DetCntIdx = 222;
//				FanGeo.m_DetN = 444;
//				FanGeo.m_DetStp = 0.95928517242269f / 444;
//				FanGeo.m_S2D = 946.7459716796875;
//				Img.m_Size.x = 4.484740011196460e+02;
//				Img.m_Size.y = 4.484740011196460e+02;
//				Img.m_Reso.x = 256;
//				Img.m_Reso.y = 256;
//				Img.m_Bias.x = 0;
//				Img.m_Bias.y = 0;
//				Img.m_Step.x = Img.m_Size.x / Img.m_Reso.x;
//				Img.m_Step.y = Img.m_Size.y / Img.m_Reso.y;
//				int angIdx;
//				std::vector<int> rowIdx;
//				std::vector<int> colIdx;
//				std::vector<double> coeffs;
//
//				std::stringstream nonZ;
//
//				for (angIdx = 0; angIdx < FanGeo.m_ViwN; angIdx++)
//				{
//					double ang = FanGeo.m_ViwBeg + angIdx* FanGeo.m_ViwStp;
//
//					genProj_AIM<double>(rowIdx, colIdx, coeffs, angIdx, ang, FanGeo, Img);
//					std::cout << angIdx << std::endl;
//				}
//				nonZ << rowIdx.size();
//				std::string F1 = "prjAIM" + nonZ.str() + ".row";
//				std::string F2 = "prjAIM" + nonZ.str() + ".col";
//				std::string F3 = "prjAIM" + nonZ.str() + ".cof";
//				std::ofstream rowFile(F1.c_str(), std::ios::binary);
//				std::ofstream colFile(F2.c_str(), std::ios::binary);
//				std::ofstream coeFile(F3.c_str(), std::ios::binary);
//
//				rowFile.write((char*) &(rowIdx[0]), sizeof(int) * rowIdx.size());
//				rowFile.close();
//				colFile.write((char*) &(colIdx[0]), sizeof(int) * colIdx.size());
//				colFile.close();
//				coeFile.write((char*) &(coeffs[0]), sizeof(double) * coeffs.size());
//				coeFile.close();
//				std::cout << coeffs.size() << std::endl;
//			}
//			else
//			{
//				std::cout << "Exit without generating\n";
//				exit(-1);
//			}
//		}
//		std::string s;
//		double val[14];
//		int i = 0;
//		while (fin >> s)
//		{
//			if ((i % 2))
//			{
//				std::stringstream ss;
//				ss << s;
//				ss >> val[i / 2];
//				std::cout << val[i / 2] << std::endl;
//			}
//			++i;
//		}
//		FanGeo.m_S2O = val[0];
//		FanGeo.m_O2D = val[1];
//		FanGeo.m_ViwN = val[2];
//		FanGeo.m_ViwBeg = val[3];
//		FanGeo.m_ViwStp = val[4];
//		FanGeo.m_DetArc = val[5];
//		FanGeo.m_DetN = val[6];
//		FanGeo.m_DetCntIdx = val[7];
//
//		FanGeo.m_DetStp = val[5] / val[6];
//		FanGeo.m_S2D = val[0] + val[1];
//		Img.m_Size.x = val[8];
//		Img.m_Size.y = val[9];
//		Img.m_Reso.x = val[10];
//		Img.m_Reso.y = val[11];
//		Img.m_Bias.x = val[12];
//		Img.m_Bias.y = val[13];
//		Img.m_Step.x = Img.m_Size.x / Img.m_Reso.x;
//		Img.m_Step.y = Img.m_Size.y / Img.m_Reso.y;
//
//		int angIdx;
//		std::vector<int> rowIdx;
//		std::vector <int> colIdx;
//		std::vector <double> coeffs;
//
//		std::stringstream nonZ;
//
//		for (angIdx = 0; angIdx < FanGeo.m_ViwN; angIdx++)
//		{
//			double ang = FanGeo.m_ViwBeg + angIdx* FanGeo.m_ViwStp;
//
//			genProj_AIM<double>(rowIdx, colIdx, coeffs, angIdx, ang, FanGeo, Img);
//			std::cout << angIdx << std::endl;
//		}
//		nonZ << rowIdx.size();
//		std::string F1 = "prjAIM" + nonZ.str() + ".row";
//		std::string F2 = "prjAIM" + nonZ.str() + ".col";
//		std::string F3 = "prjAIM" + nonZ.str() + ".cof";
//		std::ofstream rowFile(F1.c_str(), std::ios::binary);
//		std::ofstream colFile(F2.c_str(), std::ios::binary);
//		std::ofstream coeFile(F3.c_str(), std::ios::binary);
//
//		rowFile.write((char*) &(rowIdx[0]), sizeof(int) * rowIdx.size());
//		rowFile.close();
//		colFile.write((char*) &(colIdx[0]), sizeof(int) * colIdx.size());
//		colFile.close();
//		coeFile.write((char*) &(coeffs[0]), sizeof(double) * coeffs.size());
//		coeFile.close();
//		std::cout << coeffs.size() << std::endl;
//	}
//
//	void genMatrix_EA_SIDDON()
//	{
//		FanEAGeo FanGeo;
//		Image Img;
//		//¶ÁÈ¡ÅäÖÃÎÄŒþ;
//		std::ifstream fin("configurationFile.txt", std::ios::in);
//		if (!fin.is_open())
//		{
//			std::cout << "Cannot Find configuration file\n";
//			std::cout << "Would you like to use the default configuration?[y/n]";
//			char yes;
//			std::cin >> yes;
//			if (yes == 'y' || yes == 'Y')
//			{
//				FanGeo.m_S2O = 5.385200195312500e+02;
//				FanGeo.m_O2D = 4.082259521484375e+02;
//				FanGeo.m_ViwN = 360;
//				FanGeo.m_ViwBeg = 0;
//				FanGeo.m_ViwStp = 3.14159265358979323846264*2.0 / 359.0;
//				FanGeo.m_DetArc = 0.95928517242269f;
//				FanGeo.m_DetCntIdx = 222;
//				FanGeo.m_DetN = 444;
//				FanGeo.m_DetStp = 0.95928517242269f / 444;
//				FanGeo.m_S2D = 946.7459716796875;
//				Img.m_Size.x = 4.484740011196460e+02;
//				Img.m_Size.y = 4.484740011196460e+02;
//				Img.m_Reso.x = 256;
//				Img.m_Reso.y = 256;
//				Img.m_Bias.x = 0;
//				Img.m_Bias.y = 0;
//				Img.m_Step.x = Img.m_Size.x / Img.m_Reso.x;
//				Img.m_Step.y = Img.m_Size.y / Img.m_Reso.y;
//				//int angIdx;
//				std::vector<int> rowIdx;
//				std::vector <int> colIdx;
//				std::vector <double> coeffs;
//
//				std::stringstream nonZ;
//
//				genProj_SIDDON(rowIdx, colIdx, coeffs, FanGeo, Img, 1);
//
//				nonZ << rowIdx.size();
//				std::string F1 = "prjSIDDON" + nonZ.str() + ".row";
//				std::string F2 = "prjSIDDON" + nonZ.str() + ".col";
//				std::string F3 = "prjSIDDON" + nonZ.str() + ".cof";
//				std::ofstream rowFile(F1.c_str(), std::ios::binary);
//				std::ofstream colFile(F2.c_str(), std::ios::binary);
//				std::ofstream coeFile(F3.c_str(), std::ios::binary);
//
//				rowFile.write((char*) &(rowIdx[0]), sizeof(int) * rowIdx.size());
//				rowFile.close();
//				colFile.write((char*) &(colIdx[0]), sizeof(int) * colIdx.size());
//				colFile.close();
//				coeFile.write((char*) &(coeffs[0]), sizeof(double) * coeffs.size());
//				coeFile.close();
//				std::cout << coeffs.size() << std::endl;
//			}
//			else
//			{
//				std::cout << "Exit without generating\n";
//				exit(-1);
//			}
//		}
//		std::string s;
//		double val[14];
//		int i = 0;
//		while (fin >> s)
//		{
//			if ((i % 2))
//			{
//				std::stringstream ss;
//				ss << s;
//				ss >> val[i / 2];
//				std::cout << val[i / 2] << std::endl;
//			}
//			++i;
//		}
//		FanGeo.m_S2O = val[0];
//		FanGeo.m_O2D = val[1];
//		FanGeo.m_ViwN = val[2];
//		FanGeo.m_ViwBeg = val[3];
//		FanGeo.m_ViwStp = val[4];
//		FanGeo.m_DetArc = val[5];
//		FanGeo.m_DetN = val[6];
//		FanGeo.m_DetCntIdx = val[7];
//
//		FanGeo.m_DetStp = val[5] / val[6];
//		FanGeo.m_S2D = val[0] + val[1];
//		Img.m_Size.x = val[8];
//		Img.m_Size.y = val[9];
//		Img.m_Reso.x = val[10];
//		Img.m_Reso.y = val[11];
//		Img.m_Bias.x = val[12];
//		Img.m_Bias.y = val[13];
//		Img.m_Step.x = Img.m_Size.x / Img.m_Reso.x;
//		Img.m_Step.y = Img.m_Size.y / Img.m_Reso.y;
//
//		//int angIdx;
//		std::vector<int> rowIdx;
//		std::vector <int> colIdx;
//		std::vector <double> coeffs;
//
//		std::stringstream nonZ;
//		double S2O = val[0];
//		double O2D = val[1];
//		double detCntIdx = val[7];
//		double detStp = val[5] / val[6];
//		int DN = val[6];
//		double ObjStpx = val[8] / val[10];
//		double ObjStpy = val[9] / val[11];
//		int XN = val[10];
//		int YN = val[11];
//		double objCntIdxx = (XN - 1.0) * 0.5;
//		double objCntIdxy = (YN - 1.0) * 0.5;
//		int PN = val[2];
//		std::vector<double> angs(PN, 0);
//		for (unsigned int i = 0; i != PN; ++i)
//		{
//			angs[i] = val[3] + i * val[4];
//		}
//		//genProj_SIDDON(rowIdx, colIdx, coeffs, FanGeo, Img, 1);
//		genProj_SIDDON_openMP(rowIdx, colIdx, coeffs, S2O, O2D, detCntIdx, detStp, DN, ObjStpx, ObjStpy,
//			XN, YN, objCntIdxx, objCntIdxy, angs);
//
//		nonZ << rowIdx.size();
//		std::string F1 = "prjSIDDON" + nonZ.str() + ".row";
//		std::string F2 = "prjSIDDON" + nonZ.str() + ".col";
//		std::string F3 = "prjSIDDON" + nonZ.str() + ".cof";
//		std::ofstream rowFile(F1.c_str(), std::ios::binary);
//		std::ofstream colFile(F2.c_str(), std::ios::binary);
//		std::ofstream coeFile(F3.c_str(), std::ios::binary);
//
//		rowFile.write((char*) &(rowIdx[0]), sizeof(int) * rowIdx.size());
//		rowFile.close();
//		colFile.write((char*) &(colIdx[0]), sizeof(int) * colIdx.size());
//		colFile.close();
//		coeFile.write((char*) &(coeffs[0]), sizeof(double) * coeffs.size());
//		coeFile.close();
//		std::cout << coeffs.size() << std::endl;
//
//	}
//
//
//
//
//
//	template<typename T>
//	void genProj_AIM_template(
//		std::vector<int>& rowIdx,
//		std::vector<int>& colIdx,
//		std::vector<T>& weight,
//		const T S2O, const T O2D,
//		const T detCntIdx,
//		const T detStp, const int DN,
//		const T ObjStpx, const T ObjStpy,
//		const int XN, const int YN,
//		const T objCntIdxx, const T objCntIdxy,
//		const std::vector<T>& angs)
//	{
//		const int PN = angs.size();
//		T Grid[4][3];
//		T SVA[3], SVB[3];
//		const T area = ObjStpx * ObjStpy;
//		T SPoint[2];
//
//		for (unsigned int angIdx = 0; angIdx != PN; ++angIdx)
//		{
//			for (unsigned int yIdx = 0; yIdx != YN; ++yIdx)
//			{
//				for (unsigned int xIdx = 0; xIdx != XN; ++xIdx)
//				{
//					//Fetch the four points of the pixel and their projection positions
//					Grid[0][0] = (xIdx - objCntIdxx - 0.5) * ObjStpx;
//					Grid[0][1] = (yIdx - objCntIdxy - 0.5) * ObjStpy;
//					//Grid[0][2] = PosAry[yi*(XN + 1) + xi];
//					Grid[1][0] = (xIdx - objCntIdxx + 0.5) * ObjStpx;
//					Grid[1][1] = (yIdx - objCntIdxy - 0.5) * ObjStpy;
//					//Grid[1][2] = PosAry[yi*(XN + 1) + xi + 1];
//					Grid[2][0] = (xIdx - objCntIdxx + 0.5) * ObjStpx;
//					Grid[2][1] = (yIdx - objCntIdxy + 0.5) * ObjStpy;
//					//Grid[2][2] = PosAry[(yi + 1)*(XN + 1) + xi + 1];
//					Grid[3][0] = (xIdx - objCntIdxx - 0.5) * ObjStpx;
//					Grid[3][1] = (yIdx - objCntIdxy + 0.5) * ObjStpy;
//					//Grid[3][2] = PosAry[(yi + 1)*(XN + 1) + xi];
//					SortProjection<T>(Grid);//Sort the projection psotion
//
//					for (unsigned int detIdx = 0; detIdx != DN; ++detIdx)
//					{
//
//					}
//				}
//			}
//		}
//	}
//
//
//	
//	void genMatrix_ED_AIM()
//	{
//		FanEDGeo FanGeo;
//		Image Img;
//		//读取配置文件;
//		std::ifstream fin("configurationFile.txt", std::ios::in);
//
//		if (!fin.is_open())
//		{
//			std::cout << "Cannot Find configuration file\n";
//			std::cout << "Would you like to use the default configuration?[y/n]";
//			char yes;
//			std::cin >> yes;
//			if (yes == 'y' || yes == 'Y')
//			{
//				FanGeo.m_S2O = 5.385200195312500e+02;
//				FanGeo.m_O2D = 4.082259521484375e+02;
//				FanGeo.m_ViwN = 360;
//				FanGeo.m_ViwBeg = 0;
//				FanGeo.m_ViwStp = 3.14159265358979323846264*2.0 / 359.0;
//				FanGeo.m_DetSize = 500.0;
//				FanGeo.m_DetN = 888;
//				FanGeo.m_DetStp = FanGeo.m_DetSize / FanGeo.m_DetN;
//				FanGeo.m_S2D = FanGeo.m_S2O + FanGeo.m_O2D;
//				FanGeo.m_DetCntIdx = 444;
//
//				Img.m_Size.x = 4.484740011196460e+02;
//				Img.m_Size.y = 4.484740011196460e+02;
//				Img.m_Reso.x = 256;
//				Img.m_Reso.y = 256;
//				Img.m_Bias.x = 0;
//				Img.m_Bias.y = 0;
//				Img.m_Step.x = Img.m_Size.x / Img.m_Reso.x;
//				Img.m_Step.y = Img.m_Size.y / Img.m_Reso.y;
//				int angIdx;
//				std::vector<int> rowIdx;
//				std::vector <int> colIdx;
//				std::vector <double> coeffs;
//
//				std::stringstream nonZ;
//
//				for (angIdx = 0; angIdx < FanGeo.m_ViwN; angIdx++)
//				{
//					double ang = FanGeo.m_ViwBeg + angIdx* FanGeo.m_ViwStp;
//
//					genProj_AIM<double>(rowIdx, colIdx, coeffs, angIdx, ang, FanGeo, Img);
//					std::cout << angIdx << std::endl;
//				}
//				nonZ << rowIdx.size();
//				std::string F1 = "prjAIM" + nonZ.str() + ".row";
//				std::string F2 = "prjAIM" + nonZ.str() + ".col";
//				std::string F3 = "prjAIM" + nonZ.str() + ".cof";
//				std::ofstream rowFile(F1.c_str(), std::ios::binary);
//				std::ofstream colFile(F2.c_str(), std::ios::binary);
//				std::ofstream coeFile(F3.c_str(), std::ios::binary);
//
//				rowFile.write((char*) &(rowIdx[0]), sizeof(int) * rowIdx.size());
//				rowFile.close();
//				colFile.write((char*) &(colIdx[0]), sizeof(int) * colIdx.size());
//				colFile.close();
//				coeFile.write((char*) &(coeffs[0]), sizeof(double) * coeffs.size());
//				coeFile.close();
//				std::cout << coeffs.size() << std::endl;
//			}
//			else
//			{
//				std::cout << "Exit without generating\n";
//				exit(-1);
//			}
//		}
//		std::string s;
//		double val[14];
//		int i = 0;
//		while (fin >> s)
//		{
//			if ((i % 2))
//			{
//				std::stringstream ss;
//				ss << s;
//				ss >> val[i / 2];
//				std::cout << val[i / 2] << std::endl;
//			}
//			++i;
//		}
//
//		FanGeo.m_S2O = val[0];
//		FanGeo.m_O2D = val[1];
//		FanGeo.m_ViwN = val[2];
//		FanGeo.m_ViwBeg = val[3];
//		FanGeo.m_ViwStp = val[4];
//		FanGeo.m_DetSize = val[5];
//		FanGeo.m_DetN = val[6];
//		FanGeo.m_DetCntIdx = val[7];
//
//		FanGeo.m_DetStp = val[5] / val[6];
//		FanGeo.m_S2D = val[0] + val[1];
//		Img.m_Size.x = val[8];
//		Img.m_Size.y = val[9];
//		Img.m_Reso.x = val[10];
//		Img.m_Reso.y = val[11];
//		Img.m_Bias.x = val[12];
//		Img.m_Bias.y = val[13];
//		Img.m_Step.x = Img.m_Size.x / Img.m_Reso.x;
//		Img.m_Step.y = Img.m_Size.y / Img.m_Reso.y;
//
//		int angIdx;
//		std::vector<int> rowIdx;
//		std::vector <int> colIdx;
//		std::vector <double> coeffs;
//
//		std::stringstream nonZ;
//
//		for (angIdx = 0; angIdx < FanGeo.m_ViwN; angIdx++)
//		{
//			double ang = FanGeo.m_ViwBeg + angIdx* FanGeo.m_ViwStp;
//
//			genProj_AIM<double>(rowIdx, colIdx, coeffs, angIdx, ang, FanGeo, Img);
//			std::cout << angIdx << std::endl;
//		}
//
//		nonZ << rowIdx.size();
//		std::string F1 = "prjAIM" + nonZ.str() + ".row";
//		std::string F2 = "prjAIM" + nonZ.str() + ".col";
//		std::string F3 = "prjAIM" + nonZ.str() + ".cof";
//		std::ofstream rowFile(F1.c_str(), std::ios::binary);
//		std::ofstream colFile(F2.c_str(), std::ios::binary);
//		std::ofstream coeFile(F3.c_str(), std::ios::binary);
//
//		rowFile.write((char*) &(rowIdx[0]), sizeof(int) * rowIdx.size());
//		rowFile.close();
//		colFile.write((char*) &(colIdx[0]), sizeof(int) * colIdx.size());
//		colFile.close();
//		coeFile.write((char*) &(coeffs[0]), sizeof(double) * coeffs.size());
//		coeFile.close();
//		std::cout << coeffs.size() << std::endl;
//	}
//
//
//
//
//
//
//	void genMatrix_DDM_ED(
//		std::vector<int>& rowIdx,
//		std::vector<int>& colIdx,
//		std::vector<double>& weight,
//		const double S2O, const double O2D,
//		const double objSizeX, const double objSizeY,
//		const double detSize,
//		const double detCntIdx,
//		const int XN, const int YN, const int DN, const int PN,
//		const std::vector<double>& angs);
//
//
//
//
//	template<typename T>
//	inline T intersectLength(const T& fixedmin, const T& fixedmax, const T& varimin, const T& varimax)
//	{
//		const T left = (fixedmin > varimin) ? fixedmin : varimin;
//		const T right = (fixedmax < varimax) ? fixedmax : varimax;
//		return abs(right - left) * static_cast<double>(right > left);
//	}
//
//
//
//
//	//沿Y轴传播;
//	template<typename T>
//	void pushCase4(std::vector<int>& rowIdx,
//		std::vector<int>& colIdx,
//		std::vector<T>& weight,
//		const T cursourX, const T cursourY, const T S2O, const T O2D, const T objSizeX, const T objSizeY, const T detSize, const T detCntIdx, const int XN, const int YN, const int DN, const int PN,
//		const T dx, const T dy, const T dd, const T curAng, const T cosT, const T sinT, const int angIdx)
//	{
//		T summ = 0;
//		T initX, initY, initYLeft, initYRight, curDetXLeft, curDetXRight, curDetYLeft, curDetYRight;
//		T curDetX, curDetY, dirX, dirY, legth, cosAng, detPosLeft, detPosRight;
//		T detprojLength;
//		T objX;
//
//		T minY, maxY;
//		int minYIdx, maxYIdx, detIdx, ii, jj;
//		T curminy, curmaxy;
//		initX = -O2D;
//		int ridx;
//
//		T w;
//		initYLeft = -(0 - detCntIdx - 0.5) * dd; //初始det左边Y坐标
//		curDetXLeft = initX * cosT - initYLeft * sinT; //当前det左边X坐标
//		curDetYLeft = initX * sinT + initYLeft * cosT; //当前det左边Y坐标
//
//		for (detIdx = 0; detIdx != DN; ++detIdx)
//		{
//			ridx = angIdx * DN + detIdx;
//
//			summ = 0;
//
//			initYRight = -(detIdx - detCntIdx - 0.5 + 1) * dd; //初始det右边Y坐标;
//			initY = -(detIdx - detCntIdx) * dd; //初始det中心Y坐标;
//			curDetXRight = initX * cosT - initYRight * sinT; //当前det右边X坐标
//			curDetYRight = initX * sinT + initYRight * cosT; //当前det右边Y坐标
//			curDetX = initX * cosT - initY * sinT; //当前det中心X坐标
//			curDetY = initX * sinT + initY * cosT; //当前det中心Y坐标
//			dirX = curDetX - cursourX;
//			dirY = curDetY - cursourY;
//			legth = hypot(dirX, dirY);
//			dirX /= legth;
//			dirY /= legth;
//			cosAng = abs(dirX); //与Case1区别;
//
//			//这里是从X坐标算Y坐标;
//
//			detPosLeft = (0 - cursourX) / (curDetXLeft - cursourX) * (curDetYLeft - cursourY) + cursourY; //det左边界X轴上的坐标;
//			detPosRight = (0 - cursourX) / (curDetXRight - cursourX) * (curDetYRight - cursourY) + cursourY;//det右边界在x轴上的坐标;
//			detprojLength = abs(detPosRight - detPosLeft);
//
//			//沿X方向扫描;
//			for (ii = 0; ii < XN; ++ii)
//			{
//
//				objX = (ii - YN / 2.0 + 0.5) * dy;
//				minY = (objX - cursourX) / (curDetXLeft - cursourX) * (curDetYLeft - cursourY) + cursourY;
//				maxY = (objX - cursourX) / (curDetXRight - cursourX) *  (curDetYRight - cursourY) + cursourY;
//				if (minY > maxY)
//				{
//					std::swap(minY, maxY);
//				}
//
//				minYIdx = floor(minY / dy + YN / 2.0);
//				maxYIdx = ceil(maxY / dy + YN / 2.0);
//
//				if (maxYIdx <= 0)
//				{
//					continue;
//				}
//				else if (minYIdx > XN)
//				{
//					continue;
//				}
//
//				if (minYIdx < 0)
//				{
//					minYIdx = 0;
//				}
//				if (maxYIdx > YN)
//				{
//					maxYIdx = YN;
//				}
//
//
//				curminy = (-cursourX) / (objX - cursourX) * ((0 - YN / 2.0) * dy - cursourY) + cursourY;
//				for (jj = minYIdx; jj < maxYIdx; ++jj)
//				{
//					curmaxy = (-cursourX) / (objX - cursourX) * ((jj + 1 - YN / 2.0) * dy - cursourY) + cursourY;
//					if (detPosLeft > detPosRight)
//					{
//						std::swap(detPosLeft, detPosRight);
//					}
//
//					w = intersectLength<double>(detPosLeft, detPosRight, curminy, curmaxy);
//					if (w > 0)
//					{
//						rowIdx.push_back(ridx);
//						colIdx.push_back(jj * XN + ii);
//						weight.push_back(w * dx / (cosAng * detprojLength));
//					}
//					//summ += img[jj * XN + ii] * * dx;
//					curminy = curmaxy;
//				}
//			}
//
//			//proj[angIdx * DN + detIdx] = summ / (cosAng * detprojLength);
//			initYLeft = initYRight;
//			curDetXLeft = curDetXRight;
//			curDetYLeft = curDetYRight;
//		}
//	}
//
//	//沿X轴传播;
//	template<typename T>
//	void pushCase1(std::vector<int>& rowIdx,
//		std::vector<int>& colIdx,
//		std::vector<T>& weight,
//		const T cursourX, const T cursourY, const T S2O, const T O2D, const T objSizeX, const T objSizeY, const T detSize, const T detCntIdx, const int XN, const int YN, const int DN, const int PN,
//		const T dx, const T dy, const T dd, const T curAng, const T cosT, const T sinT, const int angIdx)
//	{
//		int detIdx = 0;
//		int ii = 0, jj = 0;
//		T initX, initYLeft, initYRight;
//		T curDetXLeft, curDetYLeft;
//		T curDetXRight, curDetYRight;
//		T minX, maxX;
//		int minXIdx, maxXIdx;
//		T detPosLeft;
//		T detPosRight;
//		T initY;
//		T curDetX, curDetY;
//		T dirX, dirY;
//		T legth;
//		T cosAng;
//		T detprojLength = 0;
//		T objY;
//		T summ = 0;
//		T curminx, curmaxx;
//		int ridx;
//		T w;
//		initX = -O2D;
//		initYLeft = -(0 - detCntIdx - 0.5) * dd; //初始det左边Y坐标
//		curDetXLeft = initX * cosT - initYLeft * sinT; //当前det左边X坐标
//		curDetYLeft = initX * sinT + initYLeft * cosT; //当前det左边Y坐标
//		for (detIdx = 0; detIdx != DN; ++detIdx)
//		{
//			ridx = angIdx * DN + detIdx;
//			summ = 0;
//			initYRight = initYLeft - dd;// -(detIdx - detCntIdx - 0.5 + 1) * dd; //初始det右边Y坐标;
//			initY = -(detIdx - detCntIdx) * dd; //初始det中心Y坐标;
//			curDetXRight = initX * cosT - initYRight * sinT; //当前det右边X坐标
//			curDetYRight = initX * sinT + initYRight * cosT; //当前det右边Y坐标
//			curDetX = initX * cosT - initY * sinT; //当前det中心X坐标
//			curDetY = initX * sinT + initY * cosT; //当前det中心Y坐标
//			dirX = curDetX - cursourX;
//			dirY = curDetY - cursourY;
//			legth = hypot(dirX, dirY);
//			dirX /= legth;
//			dirY /= legth;
//			cosAng = abs(dirY); //当前光线和Y轴夹角余弦
//
//			detPosLeft = (0 - cursourY) / (curDetYLeft - cursourY) * (curDetXLeft - cursourX) + cursourX; //det左边界X轴上的坐标;
//			detPosRight = (0 - cursourY) / (curDetYRight - cursourY) * (curDetXRight - cursourX) + cursourX;//det右边界在x轴上的坐标;
//			detprojLength = abs(detPosRight - detPosLeft);
//
//			for (jj = 0; jj < YN; jj++)
//			{
//				objY = (jj - YN / 2.0 + 0.5) * dy;
//				minX = (objY - cursourY) / (curDetYLeft - cursourY) * (curDetXLeft - cursourX) + cursourX;
//				maxX = (objY - cursourY) / (curDetYRight - cursourY) *  (curDetXRight - cursourX) + cursourX;
//				if (minX > maxX)
//				{
//					std::swap(minX, maxX);
//				}
//
//				minXIdx = floor(minX / dx + XN / 2.0);
//				maxXIdx = ceil(maxX / dx + XN / 2.0);
//
//				if (maxXIdx <= 0)
//				{
//					continue;
//				}
//				else if (minXIdx > XN)
//				{
//					continue;
//				}
//
//				if (minXIdx < 0)
//				{
//					minXIdx = 0;
//				}
//				if (maxXIdx > XN)
//				{
//					maxXIdx = XN;
//				}
//				curminx = (-cursourY) / (objY - cursourY) * ((minXIdx - XN / 2.0) * dx - cursourX) + cursourX;
//				for (ii = minXIdx; ii < maxXIdx; ++ii)
//				{
//
//					curmaxx = (-cursourY) / (objY - cursourY) * ((ii + 1 - XN / 2.0) * dx - cursourX) + cursourX;
//					if (detPosLeft > detPosRight)
//					{
//						std::swap(detPosLeft, detPosRight);
//					}
//
//					w = intersectLength<double>(detPosLeft, detPosRight, curminx, curmaxx);
//					if (w > 0)
//					{
//						rowIdx.push_back(ridx);
//						colIdx.push_back(jj * XN + ii);
//						weight.push_back(w * dx / (cosAng * detprojLength));
//					}
//
//					//summ += img[jj * XN + ii] * intersectLength<double>(detPosLeft, detPosRight, curminx, curmaxx) * dy;
//					curminx = curmaxx;
//				}
//			}
//			//proj[angIdx * DN + detIdx] = summ / (cosAng * detprojLength);
//			initYLeft = initYRight;
//			curDetXLeft = curDetXRight;
//			curDetYLeft = curDetYRight;
//		}
//	}
//
//
//
//	template<typename T>
//	void genMatrix_DDM_ED_template(
//		std::vector<int>& rowIdx,
//		std::vector<int>& colIdx,
//		std::vector<T>& weight,
//		const T S2O, const T O2D,
//		const T objSizeX, const T objSizeY,
//		const T detSize,
//		const T detCntIdx,
//		const int XN, const int YN, const int DN, const int PN,
//		const std::vector<T>& angs)
//	{
//		
//		int angIdx(0), detIdx(0);
//		T curang(0), cosT(0), sinT(0);
//		T cursourx(0), cursoury(0);
//		const T dd = detSize / static_cast<T>(DN);
//		const T dx = objSizeX / static_cast<T>(XN);
//		const T dy = objSizeY / static_cast<T>(YN);
//		for (angIdx = 0; angIdx < PN; ++angIdx)
//		{
//			curang = angs[angIdx];// angBeg + angIdx * angStp;
//			cosT = cos(curang);
//			sinT = sin(curang);
//			cursourx = S2O * cosT; //以X轴为准;
//			cursoury = S2O * sinT;
//
//			if (curang > PI * 0.25 && curang <= PI * 0.75 || curang >= PI * 1.25 && curang < PI * 1.75)
//			{
//				
//				pushCase1<T>(rowIdx, colIdx, weight, cursourx, cursoury, S2O, O2D,
//					objSizeX, objSizeY, detSize, detCntIdx,
//					XN, YN, DN, PN, dx, dy, dd, curang, cosT, sinT, angIdx);
//			}
//			else
//			{
//				
//				pushCase4<T>(rowIdx, colIdx, weight, cursourx, cursoury, S2O, O2D,
//					objSizeX, objSizeY, detSize, detCntIdx,
//					XN, YN, DN, PN, dx, dy, dd, curang, cosT, sinT, angIdx);
//			}
//		}
//
//
//	}
//
//
//
//
//	void genMatrix_DDM_ED(
//		std::vector<int>& rowIdx,
//		std::vector<int>& colIdx,
//		std::vector<double>& weight,
//		const double S2O, const double O2D,
//		const double objSizeX, const double objSizeY,
//		const double detSize,
//		const double detCntIdx,
//		const int XN, const int YN, const int DN, const int PN,
//		const std::vector<double>& angs)
//	{
//		genMatrix_DDM_ED_template<double>(rowIdx, colIdx, weight, S2O, O2D, objSizeX, objSizeY, detSize, detCntIdx, XN, YN, DN, PN, angs);
//	}
//
//
//	void genMatrix_ED_DDM()
//	{
//		std::ifstream fin("configurationFile.txt", std::ios::in);
//		std::vector<int> rowIdx;
//		std::vector<int> colIdx;
//		std::vector<double> weight;
//		std::vector<double> angs;
//		double S2O;
//		double O2D;
//		double detCntIdx;
//		double detStp;
//		int DN;
//		double ObjStpx, ObjStpy;
//		int XN, YN;
//		double objCntIdxx, objCntIdxy;
//
//		if (!fin.is_open())
//		{
//			std::cout << "Cannot Find configuration file\n";
//			std::cout << "Would you like to use the default configuration?[y/n]";
//			char yes;
//			std::cin >> yes;
//			if (yes == 'y' || yes == 'Y')
//			{
//				S2O = 53.85200195312500;
//				O2D = 40.82259521484375;
//				DN = 888;
//				detCntIdx = DN / 2.0 - 0.5;
//				detStp = 45.0 / DN;
//				XN = 256;
//				YN = 256;
//				ObjStpx = 20.0 / XN;
//				ObjStpy = 20.0 / YN;
//				objCntIdxx = XN / 2.0 - 0.5;
//				objCntIdxy = YN / 2.0 - 0.5;
//				angs.resize(360);
//				for (unsigned int i = 0; i != 360; ++i)
//				{
//					angs[i] = i * 0.01745329251994329576923688888889;
//				}
//
//				int angIdx;
//				std::stringstream nonZ;
//
//				//genProj_SIDDONV2(rowIdx, colIdx, weight, S2O, O2D, detCntIdx, detStp, DN, ObjStpx, ObjStpy, XN, YN, objCntIdxx, objCntIdxy, angs);
//				genMatrix_DDM_ED(rowIdx, colIdx, weight, S2O, O2D, 20.0, 20.0, 45.0,
//					444, 256, 256, 888, 360, angs);
//
//				nonZ << rowIdx.size();
//				std::string F1 = "prjDDM" + nonZ.str() + ".row";
//				std::string F2 = "prjDDM" + nonZ.str() + ".col";
//				std::string F3 = "prjDDM" + nonZ.str() + ".cof";
//				std::ofstream rowFile(F1.c_str(), std::ios::binary);
//				std::ofstream colFile(F2.c_str(), std::ios::binary);
//				std::ofstream coeFile(F3.c_str(), std::ios::binary);
//
//				rowFile.write((char*) &(rowIdx[0]), sizeof(int) * rowIdx.size());
//				rowFile.close();
//				colFile.write((char*) &(colIdx[0]), sizeof(int) * colIdx.size());
//				colFile.close();
//				coeFile.write((char*) &(weight[0]), sizeof(double) * weight.size());
//				coeFile.close();
//				std::cout << weight.size() << std::endl;
//				return;
//			}
//			else
//			{
//				std::cout << "Exit without generating\n";
//				exit(-1);
//			}
//		}
//		std::string s;
//		double val[14];
//		int i = 0;
//		while (fin >> s)
//		{
//			if ((i % 2))
//			{
//				std::stringstream ss;
//				ss << s;
//				ss >> val[i / 2];
//				std::cout << val[i / 2] << std::endl;
//			}
//			++i;
//		}
//
//
//
//
//
//		S2O = val[0];
//		O2D = val[1];
//		double ViwN = val[2];
//		double ViwBeg = val[3];
//		double ViwStp = val[4];
//		double DetSize = val[5];
//
//		DN = val[6];
//		detCntIdx = val[7];
//		detStp = val[5] / val[6];
//		double S2D = val[0] + val[1];
//
//		double Sizex = val[8];
//		double Sizey = val[9];
//
//		XN = val[10];
//		YN = val[11];
//
//		double Biasx = val[12];
//		double Biasy = val[13];
//
//		ObjStpx = Sizex / XN;
//		ObjStpy = Sizey / YN;
//
//		angs.resize(ViwN);
//
//		for (unsigned int i = 0; i != ViwN; ++i)
//		{
//			angs[i] = ViwBeg + i * ViwStp;
//		}
//
//		objCntIdxx = XN / 2.0 - 0.5;
//		objCntIdxy = YN / 2.0 - 0.5;
//
//		//int angIdx;
//
//		std::stringstream nonZ;
//
//
//		//genProj_SIDDONV2(rowIdx, colIdx, weight, S2O, O2D, detCntIdx, detStp, DN, ObjStpx, ObjStpy, XN, YN, objCntIdxx, objCntIdxy, angs);
//		genMatrix_DDM_ED(rowIdx, colIdx, weight, S2O, O2D, Sizex, Sizey, DetSize, detCntIdx, XN, YN, DN, ViwN, angs);
//
//		//genProj_SIDDON(rowIdx, colIdx, coeffs, FanGeo, Img);
//		nonZ << rowIdx.size();
//		std::string F1 = "prjDDM" + nonZ.str() + ".row";
//		std::string F2 = "prjDDM" + nonZ.str() + ".col";
//		std::string F3 = "prjDDM" + nonZ.str() + ".cof";
//		std::ofstream rowFile(F1.c_str(), std::ios::binary);
//		std::ofstream colFile(F2.c_str(), std::ios::binary);
//		std::ofstream coeFile(F3.c_str(), std::ios::binary);
//
//		rowFile.write((char*) &(rowIdx[0]), sizeof(int) * rowIdx.size());
//		rowFile.close();
//		colFile.write((char*) &(colIdx[0]), sizeof(int) * colIdx.size());
//		colFile.close();
//		coeFile.write((char*) &(weight[0]), sizeof(double) * weight.size());
//		coeFile.close();
//		std::cout << weight.size() << std::endl;
//	}
//
//
//
//
//	void genProj_SIDDONV2(std::vector<int>& rowIdx, std::vector<int>& colIdx, std::vector<double>& weight,
//		const double S2O, const double O2D,
//		const double detCntIdx,
//		const double detStp,
//		const int DN,
//		const double ObjStpx, const double ObjStpy,
//		const int XN, const int YN,
//		const double objCntIdxx, const double objCntIdxy,
//		const std::vector<double> angs)
//	{
//		const unsigned int PN = angs.size();
//		double cosT, sinT;
//		double sourx, soury;
//		double initDetx, initDety;
//		double detx, dety;
//		double dirx, diry, rleng;
//		double boxminx, boxminy, boxmaxx, boxmaxy;
//		bool intersected;
//		double tnear, tfar, ww;
//		for (unsigned int angIdx = 0; angIdx != PN; ++angIdx)
//		{
//			cosT = cos(angs[angIdx]);
//			sinT = sin(angs[angIdx]);
//			sourx = -S2O * sinT;
//			soury = S2O * cosT;
//			for (unsigned int detIdx = 0; detIdx != DN; ++detIdx)
//			{
//				initDetx = (detIdx - detCntIdx) * detStp;
//				initDety = -O2D;
//				detx = initDetx * cosT - initDety * sinT;
//				dety = initDetx * sinT + initDety * cosT;
//				dirx = detx - sourx;
//				diry = dety - soury;
//				rleng = sqrt(dirx * dirx + diry * diry);
//				dirx /= rleng;
//				diry /= rleng;
//
//				for (unsigned int yIdx = 0; yIdx != YN; ++yIdx)
//				{
//					boxminy = (yIdx - objCntIdxy - 0.5) * ObjStpy;
//					boxmaxy = (yIdx - objCntIdxy + 0.5) * ObjStpy;
//					for (unsigned int xIdx = 0; xIdx != XN; ++xIdx)
//					{
//						boxminx = (xIdx - objCntIdxx - 0.5) * ObjStpx;
//						boxmaxx = (xIdx - objCntIdxx + 0.5) * ObjStpx;
//
//						intersected = intersectBox<double>(sourx, soury,
//							dirx, diry, boxminx, boxminy,
//							boxmaxx, boxmaxy, &tnear, &tfar);
//						if (intersected)
//						{
//							rowIdx.push_back(angIdx * DN + detIdx);
//							colIdx.push_back(yIdx * XN + xIdx);
//							weight.push_back(tfar - tnear);
//						}
//
//					}
//				}
//
//
//			}
//			std::cout << angIdx << std::endl;
//		}
//	}
//
//
//
//	void genMatrix_ED_SIDDON()
//	{
//		std::ifstream fin("configurationFile.txt", std::ios::in);
//		std::vector<int> rowIdx;
//		std::vector<int> colIdx;
//		std::vector<double> weight;
//		std::vector<double> angs;
//		double S2O;
//		double O2D;
//		double detCntIdx;
//		double detStp;
//		int DN;
//		double ObjStpx, ObjStpy;
//		int XN, YN;
//		double objCntIdxx, objCntIdxy;
//
//		if (!fin.is_open())
//		{
//			std::cout << "Cannot Find configuration file\n";
//			std::cout << "Would you like to use the default configuration?[y/n]";
//			char yes;
//			std::cin >> yes;
//			if (yes == 'y' || yes == 'Y')
//			{
//				S2O = 53.85200195312500;
//				O2D = 40.82259521484375;
//				DN = 888;
//				detCntIdx = DN / 2.0 - 0.5;
//				detStp = 45.0 / DN;
//				XN = 256;
//				YN = 256;
//				ObjStpx = 20.0 / XN;
//				ObjStpy = 20.0 / YN;
//				objCntIdxx = XN / 2.0 - 0.5;
//				objCntIdxy = YN / 2.0 - 0.5;
//				angs.resize(360);
//				for (unsigned int i = 0; i != 360; ++i)
//				{
//					angs[i] = i * 0.01745329251994329576923688888889;
//				}
//
//				int angIdx;
//				std::stringstream nonZ;
//
//				genProj_SIDDONV2(rowIdx, colIdx, weight, S2O, O2D, detCntIdx, detStp, DN, ObjStpx, ObjStpy, XN, YN, objCntIdxx, objCntIdxy, angs);
//				nonZ << rowIdx.size();
//				std::string F1 = "prjSIDDON" + nonZ.str() + ".row";
//				std::string F2 = "prjSIDDON" + nonZ.str() + ".col";
//				std::string F3 = "prjSIDDON" + nonZ.str() + ".cof";
//				std::ofstream rowFile(F1.c_str(), std::ios::binary);
//				std::ofstream colFile(F2.c_str(), std::ios::binary);
//				std::ofstream coeFile(F3.c_str(), std::ios::binary);
//
//				rowFile.write((char*) &(rowIdx[0]), sizeof(int) * rowIdx.size());
//				rowFile.close();
//				colFile.write((char*) &(colIdx[0]), sizeof(int) * colIdx.size());
//				colFile.close();
//				coeFile.write((char*) &(weight[0]), sizeof(double) * weight.size());
//				coeFile.close();
//				std::cout << weight.size() << std::endl;
//				return ;
//			}
//			else
//			{
//				std::cout << "Exit without generating\n";
//				exit(-1);
//			}
//		}
//		std::string s;
//		double val[14];
//		int i = 0;
//		while (fin >> s)
//		{
//			if ((i % 2))
//			{
//				std::stringstream ss;
//				ss << s;
//				ss >> val[i / 2];
//				std::cout << val[i / 2] << std::endl;
//			}
//			++i;
//		}
//
//
//
//
//
//		S2O = val[0];
//		O2D = val[1];
//		double ViwN = val[2];
//		double ViwBeg = val[3];
//		double ViwStp = val[4];
//		double DetSize = val[5];
//
//		DN = val[6];
//		detCntIdx = val[7];
//		detStp = val[5] / val[6];
//		double S2D = val[0] + val[1];
//
//		double Sizex = val[8];
//		double Sizey = val[9];
//
//		XN = val[10];
//		YN = val[11];
//
//		double Biasx = val[12];
//		double Biasy = val[13];
//
//		ObjStpx = Sizex / XN;
//		ObjStpy = Sizey / YN;
//
//		angs.resize(ViwN);
//
//		for (unsigned int i = 0; i != ViwN; ++i)
//		{
//			angs[i] = ViwBeg + i * ViwStp;
//		}
//
//		objCntIdxx = XN / 2.0 - 0.5;
//		objCntIdxy = YN / 2.0 - 0.5;
//
//		//int angIdx;
//
//		std::stringstream nonZ;
//		genProj_SIDDONV2(rowIdx, colIdx, weight, S2O, O2D, detCntIdx, detStp, DN, ObjStpx, ObjStpy, XN, YN, objCntIdxx, objCntIdxy, angs);
//
//		//genProj_SIDDON(rowIdx, colIdx, coeffs, FanGeo, Img);
//		nonZ << rowIdx.size();
//		std::string F1 = "prjSIDDON" + nonZ.str() + ".row";
//		std::string F2 = "prjSIDDON" + nonZ.str() + ".col";
//		std::string F3 = "prjSIDDON" + nonZ.str() + ".cof";
//		std::ofstream rowFile(F1.c_str(), std::ios::binary);
//		std::ofstream colFile(F2.c_str(), std::ios::binary);
//		std::ofstream coeFile(F3.c_str(), std::ios::binary);
//
//		rowFile.write((char*) &(rowIdx[0]), sizeof(int) * rowIdx.size());
//		rowFile.close();
//		colFile.write((char*) &(colIdx[0]), sizeof(int) * colIdx.size());
//		colFile.close();
//		coeFile.write((char*) &(weight[0]), sizeof(double) * weight.size());
//		coeFile.close();
//		std::cout << weight.size() << std::endl;
//	}
//
//};
//
//
//
//
//
//
//template<typename T>
//void pushCase4(std::vector<int>& rowIdx,
//	std::vector<int>& colIdx,
//	std::vector<T>& weight,
//	const T cursourX, const T cursourY, const T S2O, const T O2D, const T objSizeX, const T objSizeY, const T detSize, const T detCntIdx, const int XN, const int YN, const int DN, const int PN,
//	const T dx, const T dy, const T dd, const T curAng, const T cosT, const T sinT, const int angIdx)
//{
//	T summ = 0;
//	T initX, initY, initYLeft, initYRight, curDetXLeft, curDetXRight, curDetYLeft, curDetYRight;
//	T curDetX, curDetY, dirX, dirY, legth, cosAng, detPosLeft, detPosRight;
//	T detprojLength;
//	T objX;
//
//	T minY, maxY;
//	int minYIdx, maxYIdx, detIdx, ii, jj;
//	T curminy, curmaxy;
//	initX = -O2D;
//	int ridx;
//
//	T w;
//	initYLeft = -(0 - detCntIdx - 0.5) * dd;
//	curDetXLeft = initX * cosT - initYLeft * sinT;
//	curDetYLeft = initX * sinT + initYLeft * cosT;
//
//	for (detIdx = 0; detIdx != DN; ++detIdx)
//	{
//		ridx = angIdx * DN + detIdx;
//
//		summ = 0;
//
//		initYRight = -(detIdx - detCntIdx - 0.5 + 1) * dd;
//		initY = -(detIdx - detCntIdx) * dd;
//		curDetXRight = initX * cosT - initYRight * sinT;
//		curDetYRight = initX * sinT + initYRight * cosT;
//		curDetX = initX * cosT - initY * sinT;
//		curDetY = initX * sinT + initY * cosT;
//		dirX = curDetX - cursourX;
//		dirY = curDetY - cursourY;
//		legth = hypot(dirX, dirY);
//		dirX /= legth;
//		dirY /= legth;
//		cosAng = abs(dirX);
//
//		detPosLeft = (0 - cursourX) / (curDetXLeft - cursourX) * (curDetYLeft - cursourY) + cursourY;
//		detPosRight = (0 - cursourX) / (curDetXRight - cursourX) * (curDetYRight - cursourY) + cursourY;
//		detprojLength = abs(detPosRight - detPosLeft);
//
//		for (ii = 0; ii < XN; ++ii)
//		{
//
//			objX = (ii - YN / 2.0 + 0.5) * dy;
//			minY = (objX - cursourX) / (curDetXLeft - cursourX) * (curDetYLeft - cursourY) + cursourY;
//			maxY = (objX - cursourX) / (curDetXRight - cursourX) *  (curDetYRight - cursourY) + cursourY;
//			if (minY > maxY)
//			{
//				std::swap(minY, maxY);
//			}
//
//			minYIdx = floor(minY / dy + YN / 2.0);
//			maxYIdx = ceil(maxY / dy + YN / 2.0);
//
//			if (maxYIdx <= 0)
//			{
//				continue;
//			}
//			else if (minYIdx > XN)
//			{
//				continue;
//			}
//
//			if (minYIdx < 0)
//			{
//				minYIdx = 0;
//			}
//			if (maxYIdx > YN)
//			{
//				maxYIdx = YN;
//			}
//
//
//			curminy = (-cursourX) / (objX - cursourX) * ((0 - YN / 2.0) * dy - cursourY) + cursourY;
//			for (jj = minYIdx; jj < maxYIdx; ++jj)
//			{
//				curmaxy = (-cursourX) / (objX - cursourX) * ((jj + 1 - YN / 2.0) * dy - cursourY) + cursourY;
//				if (detPosLeft > detPosRight)
//				{
//					std::swap(detPosLeft, detPosRight);
//				}
//
//				w = intersectLength<double>(detPosLeft, detPosRight, curminy, curmaxy);
//				if (w > 0)
//				{
//					rowIdx.push_back(ridx);
//					colIdx.push_back(jj * XN + ii);
//					weight.push_back(w * dx / (cosAng * detprojLength));
//				}
//				//summ += img[jj * XN + ii] * * dx;
//				curminy = curmaxy;
//			}
//		}
//
//		//proj[angIdx * DN + detIdx] = summ / (cosAng * detprojLength);
//		initYLeft = initYRight;
//		curDetXLeft = curDetXRight;
//		curDetYLeft = curDetYRight;
//	}
//}
//
//
//template<typename T>
//void pushCase1(std::vector<int>& rowIdx,
//	std::vector<int>& colIdx,
//	std::vector<T>& weight,
//	const T cursourX, const T cursourY, const T S2O, const T O2D, const T objSizeX, const T objSizeY, const T detSize, const T detCntIdx, const int XN, const int YN, const int DN, const int PN,
//	const T dx, const T dy, const T dd, const T curAng, const T cosT, const T sinT, const int angIdx)
//{
//	int detIdx = 0;
//	int ii = 0, jj = 0;
//	T initX, initYLeft, initYRight;
//	T curDetXLeft, curDetYLeft;
//	T curDetXRight, curDetYRight;
//	T minX, maxX;
//	int minXIdx, maxXIdx;
//	T detPosLeft;
//	T detPosRight;
//	T initY;
//	T curDetX, curDetY;
//	T dirX, dirY;
//	T legth;
//	T cosAng;
//	T detprojLength = 0;
//	T objY;
//	T summ = 0;
//	T curminx, curmaxx;
//	int ridx;
//	T w;
//	initX = -O2D;
//	initYLeft = -(0 - detCntIdx - 0.5) * dd;
//	curDetXLeft = initX * cosT - initYLeft * sinT;
//	curDetYLeft = initX * sinT + initYLeft * cosT;
//	for (detIdx = 0; detIdx != DN; ++detIdx)
//	{
//		ridx = angIdx * DN + detIdx;
//		summ = 0;
//		initYRight = initYLeft - dd;// -(detIdx - detCntIdx - 0.5 + 1) * dd;
//		initY = -(detIdx - detCntIdx) * dd;
//		curDetXRight = initX * cosT - initYRight * sinT;
//		curDetYRight = initX * sinT + initYRight * cosT;
//		curDetX = initX * cosT - initY * sinT;
//		curDetY = initX * sinT + initY * cosT;
//		dirX = curDetX - cursourX;
//		dirY = curDetY - cursourY;
//		legth = hypot(dirX, dirY);
//		dirX /= legth;
//		dirY /= legth;
//		cosAng = abs(dirY);
//
//		detPosLeft = (0 - cursourY) / (curDetYLeft - cursourY) * (curDetXLeft - cursourX) + cursourX; //det左边界X轴上的坐标;
//		detPosRight = (0 - cursourY) / (curDetYRight - cursourY) * (curDetXRight - cursourX) + cursourX;//det右边界在x轴上的坐标;
//		detprojLength = abs(detPosRight - detPosLeft);
//
//		for (jj = 0; jj < YN; jj++)
//		{
//			objY = (jj - YN / 2.0 + 0.5) * dy;
//			minX = (objY - cursourY) / (curDetYLeft - cursourY) * (curDetXLeft - cursourX) + cursourX;
//			maxX = (objY - cursourY) / (curDetYRight - cursourY) *  (curDetXRight - cursourX) + cursourX;
//			if (minX > maxX)
//			{
//				std::swap(minX, maxX);
//			}
//
//			minXIdx = floor(minX / dx + XN / 2.0);
//			maxXIdx = ceil(maxX / dx + XN / 2.0);
//
//			if (maxXIdx <= 0)
//			{
//				continue;
//			}
//			else if (minXIdx > XN)
//			{
//				continue;
//			}
//
//			if (minXIdx < 0)
//			{
//				minXIdx = 0;
//			}
//			if (maxXIdx > XN)
//			{
//				maxXIdx = XN;
//			}
//			curminx = (-cursourY) / (objY - cursourY) * ((minXIdx - XN / 2.0) * dx - cursourX) + cursourX;
//			for (ii = minXIdx; ii < maxXIdx; ++ii)
//			{
//
//				curmaxx = (-cursourY) / (objY - cursourY) * ((ii + 1 - XN / 2.0) * dx - cursourX) + cursourX;
//				if (detPosLeft > detPosRight)
//				{
//					std::swap(detPosLeft, detPosRight);
//				}
//
//				w = intersectLength<double>(detPosLeft, detPosRight, curminx, curmaxx);
//				if (w > 0)
//				{
//					rowIdx.push_back(ridx);
//					colIdx.push_back(jj * XN + ii);
//					weight.push_back(w * dx / (cosAng * detprojLength));
//				}
//
//				//summ += img[jj * XN + ii] * intersectLength<double>(detPosLeft, detPosRight, curminx, curmaxx) * dy;
//				curminx = curmaxx;
//			}
//		}
//		//proj[angIdx * DN + detIdx] = summ / (cosAng * detprojLength);
//		initYLeft = initYRight;
//		curDetXLeft = curDetXRight;
//		curDetYLeft = curDetYRight;
//	}
//}
//
//
//
//template<typename T>
//void genMatrix_DDM_ED_template(
//	std::vector<int>& rowIdx,
//	std::vector<int>& colIdx,
//	std::vector<T>& weight,
//	const T S2O, const T O2D,
//	const T objSizeX, const T objSizeY,
//	const T detSize,
//	const T detCntIdx,
//	const int XN, const int YN, const int DN, const int PN,
//	const std::vector<T>& angs)
//{
//	int angIdx(0), detIdx(0);
//	T curang(0), cosT(0), sinT(0);
//	T cursourx(0), cursoury(0);
//	const T dd = detSize / static_cast<T>(DN);
//	const T dx = objSizeX / static_cast<T>(XN);
//	const T dy = objSizeY / static_cast<T>(YN);
//	for (angIdx = 0; angIdx < PN; ++angIdx)
//	{
//		curang = angs[angIdx];// angBeg + angIdx * angStp;
//		cosT = cos(curang);
//		sinT = sin(curang);
//		cursourx = S2O * cosT;
//		cursoury = S2O * sinT;
//
//		if ((curang > ConstantValues::PI_f * 0.25 && curang <= ConstantValues::PI_f * 0.75) || 
//			(curang >= ConstantValues::PI_f * 1.25 && curang < ConstantValues::PI_f * 1.75))
//		{
//			pushCase1<T>(rowIdx, colIdx, weight, cursourx, cursoury, S2O, O2D,
//				objSizeX, objSizeY, detSize, detCntIdx,
//				XN, YN, DN, PN, dx, dy, dd, curang, cosT, sinT, angIdx);
//		}
//		else
//		{
//			pushCase4<T>(rowIdx, colIdx, weight, cursourx, cursoury, S2O, O2D,
//				objSizeX, objSizeY, detSize, detCntIdx,
//				XN, YN, DN, PN, dx, dy, dd, curang, cosT, sinT, angIdx);
//		}
//	}
//
//
//}
//
//
//
//void genMatrix_DDM_ED(
//	std::vector<int>& rowIdx,
//	std::vector<int>& colIdx,
//	std::vector<double>& weight,
//	const double S2O, const double O2D,
//	const double objSizeX, const double objSizeY,
//	const double detSize,
//	const double detCntIdx,
//	const int XN, const int YN, const int DN, const int PN,
//	const std::vector<double>& angs)
//{
//	genMatrix_DDM_ED_template<double>(rowIdx, colIdx, weight, S2O, O2D, objSizeX, objSizeY, detSize, detCntIdx, XN, YN, DN, PN, angs);
//}
//
//
//
//void genMatrix_DDM_ED(
//	std::vector<int>& rowIdx,
//	std::vector<int>& colIdx,
//	std::vector<float>& weight,
//	const float S2O, const float O2D,
//	const float objSizeX, const float objSizeY,
//	const float detSize,
//	const float detCntIdx,
//	const int XN, const int YN, const int DN, const int PN,
//	const std::vector<float>& angs)
//{
//	genMatrix_DDM_ED_template<float>(rowIdx, colIdx, weight, S2O, O2D, objSizeX, objSizeY, detSize, detCntIdx, XN, YN, DN, PN, angs);
//}
//
//
