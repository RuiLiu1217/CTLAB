
#include <cstdlib>
#include <cmath>
#include <utility>
#include <vector_types.h>
#include <vector_functions.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <functional>
#include <algorithm>
#include <omp.h>

#include "genSysMatrix.hpp"

namespace genSysMatrix
{
	static const float PI = 3.14159265359f;
	static const float EPSILON = 1.0E-9;
	
	/// \brief Fan Beam Equal Angle Detector based CT system
	class FanEAGeo
	{
	public:
		float m_S2O; ///< source to object distance
		float m_O2D; ///< object to detector distance
		float m_S2D; ///< source to detector distance
		int m_ViwN; ///< view number
		float m_ViwBeg; ///< Begin viewer number
		float m_ViwStp; ///< View step size

		float m_DetArc; ///< Detector Arc angle
		int m_DetN; ///< Detector cells number in one row
		float m_DetStp; ///< Detector cells size
		float m_DetCntIdx; ///< The index of the center of the detector
	public:
		FanEAGeo(void);
		~FanEAGeo(void) {};
		FanEAGeo(const FanEAGeo& rhs);
		/// \brief constructor
		/// \param S2O source to object distance
		/// \param O2D object to detector distance
		/// \param ViwN view number
		/// \param ViwBeg the begin view
		/// \param ViwEnd the End view
		/// \param DetArc the detector Arc
		/// \param DetN the number of detector cells on one row
		FanEAGeo(const float S2O, const float O2D, const  int ViwN,
			const float ViwBeg, const float ViwEnd, const float DetArc,
			const  int DetN);
	};


	FanEAGeo::FanEAGeo(void)
	{
		m_DetArc = 0.95928517242269f / 2;
		m_DetN = 444;
		m_DetCntIdx = 222;
		m_DetStp = m_DetArc / m_DetN;
		m_O2D = 4.082259521484375e+02;

		m_S2O = 5.385200195312500e+02;
		m_S2D = m_S2O + m_O2D;
		m_ViwBeg = 0.0f;// (-2.668082275390625e+02 / 180.0* 3.14159265358979323846264);
		m_ViwN = 61;
		m_ViwStp = 3.14159265358979323846264f * 2.0f / m_ViwN;
	}

	FanEAGeo::FanEAGeo(const FanEAGeo& rhs)
	{
		m_DetArc = rhs.m_DetArc;
		m_DetN = rhs.m_DetN;
		m_DetCntIdx = rhs.m_DetCntIdx;
		m_DetStp = rhs.m_DetStp;
		m_O2D = rhs.m_O2D;
		m_S2O = rhs.m_S2O;
		m_S2D = rhs.m_S2D;
		m_ViwBeg = rhs.m_ViwBeg;
		m_ViwN = rhs.m_ViwN;
		m_ViwStp = rhs.m_ViwStp;
	}


	FanEAGeo::FanEAGeo(const float S2O, const float O2D, const  int ViwN,
		const float ViwBeg, const float ViwEnd, const float DetArc,
		const  int DetN) :m_S2O(S2O), m_O2D(O2D), m_S2D(S2O + O2D),
		m_ViwN(ViwN), m_ViwBeg(ViwBeg), m_ViwStp((ViwEnd - ViwBeg) / float(ViwN)),
		m_DetArc(DetArc), m_DetN(DetN), m_DetStp(DetArc / DetN),
		m_DetCntIdx(DetN * 0.5f - 0.5f) {}



	/// \brief Image configuration class
	class Image
	{
	public:
		int2 m_Reso; ///< Image resolution
		float2 m_Size;///< Image size
		float2 m_Step; ///< Image Step
		float2 m_Bias; ///< The bias of the image
	public:
		/// \brief constructor
		Image(void);
		/// \brief destructor
		~Image(void) {};
		/// \brief copy constructor
		Image(const Image& rhs);
		/// \brief constructor
		Image(
			const int resoL,///< resolution on length direction
			const int resoW,///< resolution on width direction
			const float sizeL, ///< length size of the image
			const float sizeW,///< width size of the image
			const float BiasL, ///< bias on length direction
			const float BiasW ///<bias on width direction
			);
	};
	
	Image::Image(void)
	{
		m_Bias.x = 0.0f;
		m_Bias.y = 0.0f;  //ÕâžöÆ«ÒÆµÄµ¥Î»ÊÇÕæÊµÎïÀíµ¥Î»;
		m_Reso.x = m_Reso.y = 512;
		m_Size.x = m_Size.y = 4.484740011196460e+02;
		m_Step.x = m_Size.x / m_Reso.x;
		m_Step.y = m_Size.y / m_Reso.y;
	}

	Image::Image(const Image& rhs)
	{
		m_Bias = rhs.m_Bias;
		m_Reso = rhs.m_Reso;
		m_Size = rhs.m_Size;
		m_Step = rhs.m_Step;
	}



	Image::Image(
		const int resoL,
		const int resoW,
		const float sizeL,
		const float sizeW,
		const float BiasL,
		const float BiasW) :m_Reso(make_int2(resoL, resoW)),
		m_Size(make_float2(sizeL, sizeW)),
		m_Step(make_float2(sizeL / resoL, sizeW / resoW)),
		m_Bias(make_float2(BiasL, BiasW)) {}


	template<typename T>
	inline bool IS_ZERO(const T& x)
	{
		return ((x < EPSILON) && (x > -EPSILON));
	}



	float2 rotation(const float2& p, const float& cosT, const float& sinT)
	{
		float2 curP;
		curP.x = p.x * cosT - p.y * sinT;
		curP.y = p.x * sinT + p.y * cosT;
		return curP;
	}

	double2 rotation(const double2& p, const double& cosT, const double& sinT)
	{
		double2 curP;
		curP.x = p.x * cosT - p.y * sinT;
		curP.y = p.x * sinT + p.y * cosT;
		return curP;
	}

	struct Ray2D
	{
	public:
		float2 o;
		float2 d;
	};

	struct Ray2Dd
	{
	public:
		double2 o;
		double2 d;
	};


	double2 operator-(const double2& a, const double2& b)
	{
		double2 res;
		res.x = a.x - b.x;
		res.y = a.y - b.y;
		return res;
	}

	float2 operator-(const float2& a, const float2& b)
	{
		float2 res;
		res.x = a.x - b.x;
		res.y = a.y - b.y;
		return res;
	}

	double2 normalize(const double2& o)
	{
		double l = hypot(o.x, o.y);
		double2 res;
		res.x = o.x / l;
		res.y = o.y / l;
		return res;
	}


	float2 normalize(const float2& o)
	{
		float l = hypot(o.x, o.y);
		float2 res;
		res.x = o.x / l;
		res.y = o.y / l;
		return res;
	}

	inline float dev_pFun(const float& alpha, const float& pstart, const float&pend)
	{
		return pstart + alpha * (pend - pstart);
	}


	inline double dev_pFun(const double& alpha, const double& pstart, const double&pend)
	{
		return pstart + alpha * (pend - pstart);
	}


	inline float dev_alpha_IFun(const float& b, const float& d, const float& pstart, const float& pend, const unsigned int& i)
	{
		if (!IS_ZERO(pend - pstart))
		{
			return ((b + (float) i*d) - pstart) / (pend - pstart);
		}
		else return 1000;//((b + i*d)-pstart)/(1e-6);
	}

	inline double dev_alpha_IFun(const double& b, const double& d, const double& pstart, const double& pend, const unsigned int& i)
	{
		if (!IS_ZERO(pend - pstart))
		{
			return ((b + (double) i*d) - pstart) / (pend - pstart);
		}
		else return 1000;//((b + i*d)-pstart)/(1e-6);
	}


	inline float dev_varphiFun(const float& alpha, const float& b, const float& d, const float& pstart, const float& pend)
	{
		return (dev_pFun(alpha, pstart, pend) - b) / d;
	}


	inline double dev_varphiFun(const double& alpha, const double& b, const double& d, const double& pstart, const double& pend)
	{
		return (dev_pFun(alpha, pstart, pend) - b) / d;
	}


	inline void dev_minmaxIdxFun(
		const float& pstart, const float& pend,
		const float& b, const float& d,
		const float& alphaMIN, const float& alphaMAX,
		const float& alphaPmin, const float& alphaPmax,
		const unsigned int& Nplane, int* imin, int* imax)
	{
		if (pstart < pend)
		{
			if (IS_ZERO(alphaMIN - alphaPmin))
			{
				*imin = 1;
			}
			else
			{
				*imin = static_cast<int>(ceil(dev_varphiFun(alphaMIN, b, d, pstart, pend)));
			}
			if (IS_ZERO(alphaMAX - alphaPmax))
			{
				*imax = Nplane - 1;
			}
			else
			{
				*imax = static_cast<int>(dev_varphiFun(alphaMAX, b, d, pstart, pend));
			}
		}
		else
		{
			if (IS_ZERO(alphaMIN - alphaPmin))
			{
				*imax = Nplane - 2;
			}
			else
			{
				*imax = static_cast<int>(dev_varphiFun(alphaMIN, b, d, pstart, pend));
			}
			if (IS_ZERO(alphaMAX - alphaPmax))
			{
				*imin = 0;
			}
			else
			{
				*imin = static_cast<int>(ceil(dev_varphiFun(alphaMAX, b, d, pstart, pend)));
			}
		}
	}


	inline void dev_minmaxIdxFun(
		const double& pstart, const double& pend,
		const double& b, const double& d,
		const double& alphaMIN, const double& alphaMAX,
		const double& alphaPmin, const double& alphaPmax,
		const unsigned int& Nplane, int* imin, int* imax)
	{
		if (pstart < pend)
		{
			if (IS_ZERO(alphaMIN - alphaPmin))
			{
				*imin = 1;
			}
			else
			{
				*imin = static_cast<int>(ceil(dev_varphiFun(alphaMIN, b, d, pstart, pend)));
			}
			if (IS_ZERO(alphaMAX - alphaPmax))
			{
				*imax = Nplane - 1;
			}
			else
			{
				*imax = static_cast<int>(dev_varphiFun(alphaMAX, b, d, pstart, pend));
			}
		}
		else
		{
			if (IS_ZERO(alphaMIN - alphaPmin))
			{
				*imax = Nplane - 2;
			}
			else
			{
				*imax = static_cast<int>(dev_varphiFun(alphaMIN, b, d, pstart, pend));
			}
			if (IS_ZERO(alphaMAX - alphaPmax))
			{
				*imin = 0;
			}
			else
			{
				*imin = static_cast<int>(ceil(dev_varphiFun(alphaMAX, b, d, pstart, pend)));
			}
		}
	}



	inline float dev_alphaU_Fun(const float& d, const float& startx, const float& endx)
	{
		if (IS_ZERO(startx - endx))
		{
			return 1000.0f;//(d/1e-6);
		}
		return d / abs(startx - endx);
	}

	inline double dev_alphaU_Fun(const double& d, const double& startx, const double& endx)
	{
		if (IS_ZERO(startx - endx))
		{
			return 1000.0f;//(d/1e-6);
		}
		return d / abs(startx - endx);
	}


	inline  int dev_iu_Fun(const float& start, const float& end)
	{
		return (start < end) ? 1 : -1;
	}


	inline  int dev_iu_Fun(const double& start, const double& end)
	{
		return (start < end) ? 1 : -1;
	}




	void pushMatrix(
		float startX, float startY,
		float endX, float endY,
		float bx, float by,
		float dx, float dy,
		const int objResoLen,
		const int objResoWid,
		int& rowidx,
		std::vector<int>& rowIdx,
		std::vector<int>& colIdx,
		std::vector<float>& wgt)
	{
		const float dirX(endX - startX);
		const float dirY(endY - startY);
		const float lengthSq = dirX * dirX + dirY * dirY;
		const float dconv = sqrt(lengthSq);
		int imin, imax, jmin, jmax;
		const float alphaxmin = fminf(dev_alpha_IFun(bx, dx, startX, endX, 0), dev_alpha_IFun(bx, dx, startX, endX, objResoLen));
		const float alphaxmax = fmaxf(dev_alpha_IFun(bx, dx, startX, endX, 0), dev_alpha_IFun(bx, dx, startX, endX, objResoLen));
		const float alphaymin = fminf(dev_alpha_IFun(by, dy, startY, endY, 0), dev_alpha_IFun(by, dy, startY, endY, objResoWid));
		const float alphaymax = fmaxf(dev_alpha_IFun(by, dy, startY, endY, 0), dev_alpha_IFun(by, dy, startY, endY, objResoWid));

		const float alphaMIN = fmaxf(alphaxmin, alphaymin);
		const float alphaMAX = fminf(alphaxmax, alphaymax);
		dev_minmaxIdxFun(startX, endX, bx, dx, alphaMIN, alphaMAX, alphaxmin, alphaxmax, objResoLen + 1, &imin, &imax);
		dev_minmaxIdxFun(startY, endY, by, dy, alphaMIN, alphaMAX, alphaymin, alphaymax, objResoWid + 1, &jmin, &jmax);

		float alphaX = (startX < endX) ? dev_alpha_IFun(bx, dx, startX, endX, imin) : dev_alpha_IFun(bx, dx, startX, endX, imax);
		float alphaY = (startY < endY) ? dev_alpha_IFun(by, dy, startY, endY, jmin) : dev_alpha_IFun(by, dy, startY, endY, jmax);

		int Np = static_cast<int>(fabsf(imax - imin + 1.0f) + fabsf(jmax - jmin + 1.0f) + 4.0f);
		const float alphaxu = dev_alphaU_Fun(dx, startX, endX);
		const float alphayu = dev_alphaU_Fun(dy, startY, endY);

		float alphaC = alphaMIN;

		int i = static_cast<int>(dev_varphiFun(alphaMIN* 1.00003f, bx, dx, startX, endX));
		int j = static_cast<int>(dev_varphiFun(alphaMIN* 1.00003f, by, dy, startY, endY));

		const int iuu = dev_iu_Fun(startX, endX);
		const int juu = dev_iu_Fun(startY, endY);

		float d12(0.0f);
		float weight(0.0f);
		unsigned int repIdx(0);
		unsigned int colidx(0);
		while (repIdx != Np)
		{
			if (i < 0 || i >= objResoLen || j < 0 || j >= objResoWid)
			{
				break;
			}
			if (alphaX <= alphaY)
			{
				colidx = j * objResoLen + i;
				weight = (alphaX - alphaC) * dconv;

				wgt.push_back(weight);
				rowIdx.push_back(rowidx);
				colIdx.push_back(colidx);

				d12 += weight;
				i += iuu;
				alphaC = alphaX;
				alphaX += alphaxu;
			}
			else
			{
				colidx = j * objResoLen + i;
				weight = (alphaY - alphaC) * dconv;

				wgt.push_back(weight);
				rowIdx.push_back(rowidx);
				colIdx.push_back(colidx);

				d12 += weight;
				j += juu;

				alphaC = alphaY;
				alphaY += alphayu;
			}
			++repIdx;
		}
	}




	void pushMatrix(
		double startX, double startY,
		double endX, double endY,
		double bx, double by,
		double dx, double dy,
		const int objResoLen,
		const int objResoWid,
		int& rowidx,
		std::vector<int>& rowIdx,
		std::vector<int>& colIdx,
		std::vector<double>& wgt)
	{
		const double dirX(endX - startX);
		const double dirY(endY - startY);
		const double lengthSq = dirX * dirX + dirY * dirY;
		const double dconv = sqrt(lengthSq);
		int imin, imax, jmin, jmax;
		const double alphaxmin = fminf(dev_alpha_IFun(bx, dx, startX, endX, 0), dev_alpha_IFun(bx, dx, startX, endX, objResoLen));
		const double alphaxmax = fmaxf(dev_alpha_IFun(bx, dx, startX, endX, 0), dev_alpha_IFun(bx, dx, startX, endX, objResoLen));
		const double alphaymin = fminf(dev_alpha_IFun(by, dy, startY, endY, 0), dev_alpha_IFun(by, dy, startY, endY, objResoWid));
		const double alphaymax = fmaxf(dev_alpha_IFun(by, dy, startY, endY, 0), dev_alpha_IFun(by, dy, startY, endY, objResoWid));

		const double alphaMIN = fmaxf(alphaxmin, alphaymin);
		const double alphaMAX = fminf(alphaxmax, alphaymax);
		dev_minmaxIdxFun(startX, endX, bx, dx, alphaMIN, alphaMAX, alphaxmin, alphaxmax, objResoLen + 1, &imin, &imax);
		dev_minmaxIdxFun(startY, endY, by, dy, alphaMIN, alphaMAX, alphaymin, alphaymax, objResoWid + 1, &jmin, &jmax);

		double alphaX = (startX < endX) ? dev_alpha_IFun(bx, dx, startX, endX, imin) : dev_alpha_IFun(bx, dx, startX, endX, imax);
		double alphaY = (startY < endY) ? dev_alpha_IFun(by, dy, startY, endY, jmin) : dev_alpha_IFun(by, dy, startY, endY, jmax);

		int Np = static_cast<int>(abs(imax - imin + 1.0f) + abs(jmax - jmin + 1.0f) + 4.0f);
		const double alphaxu = dev_alphaU_Fun(dx, startX, endX);
		const double alphayu = dev_alphaU_Fun(dy, startY, endY);

		double alphaC = alphaMIN;
		double talpha = std::min(alphaX, alphaY);

		int i = static_cast<int>(dev_varphiFun((talpha + alphaMIN) * 0.5, bx, dx, startX, endX));
		int j = static_cast<int>(dev_varphiFun((talpha + alphaMIN) * 0.5, by, dy, startY, endY));

		//int i = static_cast<int>(dev_varphiFun(alphaMIN* 1.00003f, bx, dx, startX, endX));
		//int j = static_cast<int>(dev_varphiFun(alphaMIN* 1.00003f, by, dy, startY, endY));

		const int iuu = dev_iu_Fun(startX, endX);
		const int juu = dev_iu_Fun(startY, endY);

		double d12(0.0f);
		double weight(0.0f);
		unsigned int repIdx(0);
		unsigned int colidx(0);
		while (repIdx != Np)
		{
			if (i < 0 || i >= objResoLen || j < 0 || j >= objResoWid)
			{
				break;
			}
			if (alphaX <= alphaY)
			{
				colidx = j * objResoLen + i;
				weight = (alphaX - alphaC) * dconv;

				wgt.push_back(weight);
				rowIdx.push_back(rowidx);
				colIdx.push_back(colidx);

				d12 += weight;
				i += iuu;
				alphaC = alphaX;
				alphaX += alphaxu;
			}
			else
			{
				colidx = j * objResoLen + i;
				weight = (alphaY - alphaC) * dconv;

				wgt.push_back(weight);
				rowIdx.push_back(rowidx);
				colIdx.push_back(colidx);

				d12 += weight;
				j += juu;

				alphaC = alphaY;
				alphaY += alphayu;
			}
			++repIdx;
		}
	}






	void genProj_SIDDON(
		std::vector<int>& rowIdx,
		std::vector<int>& colIdx,
		std::vector<float>& weight,
		const FanEAGeo& FanGeo,
		const Image& Img,
		const unsigned int& sampNum)
	{
		float2 MINO = make_float2(
			-Img.m_Size.x / 2.0f + Img.m_Bias.x,
			-Img.m_Size.y / 2.0f + Img.m_Bias.y);

		float curAng = 0;
		float cosT = 0;
		float sinT = 0;
		Ray2D ray;

		//unsigned int detId;
		float ang(0); //bias angle from the center of the Fan Beams

		float2 curDetPos; //the current detector element position;
		//float totWeight(0);

		float smallDetStep = FanGeo.m_DetStp / sampNum; //\CF²\C9\D1\F9\BA\F3\B5\C4detector\B2\BD\B3\A4;
		float cntSmallDet = sampNum * 0.5f;
		float realAng = 0;
		unsigned int angIdx = 0;
		unsigned int detIdx = 0;
		unsigned int subDetIdx = 0;

		int rowidx = 0;
		for (angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
		{
			//Current rotation angle;
			curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
			cosT = cosf(curAng);
			sinT = sinf(curAng);
			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);

			for (detIdx = 0; detIdx != FanGeo.m_DetN; ++detIdx)
			{

				rowidx = angIdx * FanGeo.m_DetN + detIdx;

				//Calculate current Angle
				ang = ((float) detIdx - FanGeo.m_DetCntIdx + 0.5f) * FanGeo.m_DetStp;

				for (subDetIdx = 0; subDetIdx != sampNum; ++subDetIdx)
				{
					//correct ang
					realAng = ang + (static_cast<float>(subDetIdx) -cntSmallDet + 0.5f) *smallDetStep;
					// current detector element position;
					curDetPos = rotation(make_float2(sinf(realAng) * FanGeo.m_S2D, -cosf(realAng) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);

					// X-ray direction
					ray.d = normalize(curDetPos - ray.o);

					pushMatrix(ray.o.x, ray.o.y,
						curDetPos.x, curDetPos.y,
						MINO.x, MINO.y,
						Img.m_Step.x, Img.m_Step.y,
						Img.m_Reso.x, Img.m_Reso.y,
						rowidx, rowIdx, colIdx, weight);
				}
			}
		}

	}

	template<typename T>
	inline bool intersectBox(
		const T& sourX, const T& sourY, const T& detX, const T& detY,
		const T& boxminx, const T& boxminy, const T& boxmaxx, const T& boxmaxy,
		T *tnear, T *tfar)
	{
		T dist = std::hypot(detX - sourX, detY - sourY);
		T dirX = (detX - sourX) / dist;
		T dirY = (detY - sourY) / dist;

		T invRx = 1.0 / dirX;
		T invRy = 1.0 / dirY;
		T tbotx = invRx * (boxminx - sourX);
		T tboty = invRy * (boxminy - sourY);
		T ttopx = invRx * (boxmaxx - sourX);
		T ttopy = invRy * (boxmaxy - sourY);

		T tminx = std::min(ttopx, tbotx);
		T tminy = std::min(ttopy, tboty);
		T tmaxx = std::max(ttopx, tbotx);
		T tmaxy = std::max(ttopy, tboty);

		T largest_tmin = std::max(tminx, tminy);
		T smallest_tmax = std::min(tmaxx, tmaxy);

		*tnear = largest_tmin;
		*tfar = smallest_tmax;

		return smallest_tmax > largest_tmin;
	}




	void genProj_SIDDON(
		std::vector<int>& rowIdx,
		std::vector<int>& colIdx,
		std::vector<double>& weight,
		const std::vector<double> angs,
		const double S2O,
		const double O2D,
		const double detCntIdx,
		const double detStp,
		const double dx,
		const double dy,
		const int prjNum,
		const int detNum,
		const int XN,
		const int YN)
	{
		int ridx = 0, cidx = 0;
		bool hit = false;
		double tnear, tfar;
		double curAng, cosT, sinT, curSourX, curSourY;
		double curDetX, curDetY;
		double boxminx, boxminy, boxmaxx, boxmaxy;
		const double minX = -static_cast<double>(XN) / 2.0 * dx;
		const double minY = -static_cast<double>(YN) / 2.0 * dy;
		double detang;
		const double S2D = S2O + O2D;
		double XX, YY;
		for (unsigned int prjIdx = 0; prjIdx != prjNum; ++prjIdx)
		{
			curAng = angs[prjIdx];
			cosT = cos(curAng);
			sinT = sin(curAng);
			curSourX = -sinT * S2O;
			curSourY = cosT * S2O;

			for (unsigned int detIdx = 0; detIdx != detNum; ++detIdx)
			{

				//Calculate current detector position
				//(detIdx - detCntIdx) * detStp
				detang = ((double) detIdx - detCntIdx) * detStp;

				// current detector element position;
				XX = sin(detang) * S2D;
				YY = -cos(detang) * S2D + S2O;
				curDetX = XX * cosT - YY * sinT;
				curDetY = XX * sinT + YY * cosT;

				ridx = prjIdx * detNum + detIdx;
				for (unsigned int yIdx = 0; yIdx != YN; ++yIdx)
				{
					for (unsigned int xIdx = 0; xIdx != XN; ++xIdx)
					{
						boxminx = xIdx * dx + minX;
						boxminy = yIdx * dy * minY;
						boxmaxx = boxminx + dx;
						boxmaxy = boxminy + dy;

						cidx = yIdx * XN + xIdx;

						tnear = 0;
						tfar = 0;
						hit = intersectBox(curSourX, curSourY, curDetX, curDetY,
							boxminx, boxminy, boxmaxx, boxmaxy, &tnear, &tfar);
						if (hit)
						{
							rowIdx.push_back(ridx);
							colIdx.push_back(cidx);
							weight.push_back(tfar - tnear);
						}
					}
				}
			}
			std::cout << "current projection angle index = " << prjIdx << "\n";
		}
	}



	template<typename T>
	inline bool intersectBox(const T sourx, const T soury,
		const T dirx, const T diry,
		const T boxminx, const T boxminy,
		const T boxmaxx, const T boxmaxy,
		T& tnear, T& tfar)
	{
		// compute intersection of ray with all six bbox planes
		T invRx = 1.0 / dirx;
		T invRy = 1.0 / diry;

		if (dirx != 0)
		{
			invRx = 1.0 / dirx;
		}
		else
		{
			invRx = 1.0E120;
		}

		if (diry != 0)
		{
			invRy = 1.0 / diry;
		}
		else
		{
			invRy = 1.0E120;
		}
		T tbotx = invRx * (boxminx - sourx);
		T tboty = invRy * (boxminy - soury);
		T ttopx = invRx * (boxmaxx - sourx);
		T ttopy = invRy * (boxmaxy - soury);
		T tminx = (ttopx < tbotx) ? ttopx : tbotx;// min(ttopx, tbotx);
		T tminy = (ttopy < tboty) ? ttopy : tboty;// min(ttopy, tboty);
		T tmaxx = (ttopx >= tbotx) ? ttopx : tbotx;// max(ttopx, tbotx);
		T tmaxy = (ttopy >= tboty) ? ttopy : tboty;// max(ttopy, tboty);

		T largest_tmin = (tminx >= tminy) ? tminx : tminy;// max(tminx, tminy);
		T smalles_tmax = (tmaxx < tmaxy) ? tmaxx : tmaxy;// min(tmaxx, tmaxy);
		tnear = largest_tmin;
		tfar = smalles_tmax;

		return smalles_tmax > largest_tmin;
	}




	// generate the COO
	void pushValuesCOO(
		const int sarIdx,
		const int endIdx,
		const double S2O,
		const double S2D,
		const int DN,
		const int XN,
		const int YN,
		const double detCntIdx,
		const double objCntIdxy,
		const double objCntIdxx,
		const double ObjStpx,
		const double ObjStpy,
		const double detStp,
		const std::vector<double>& angs,
		std::vector<int>& rowIdx,
		std::vector<int>& colIdx,
		std::vector<double>& weight)
	{
		double cosT, sinT, sourx, soury, beta, sinBeta, cosBeta,
			initDetX, initDetY, curDetX, curDetY, dirx, diry,
			rleng, boxminx, boxminy, boxmaxx, boxmaxy, tnear, tfar;
		bool intersected;

		for (int angIdx = sarIdx; angIdx != endIdx; ++angIdx)
		{
			cosT = cos(angs[angIdx]);
			sinT = sin(angs[angIdx]);
			sourx = -S2O * sinT;
			soury = S2O * cosT;

			//std::cout << angIdx << std::endl;
			for (unsigned int detIdx = 0; detIdx != DN; ++detIdx)
			{
				beta = (detIdx - detCntIdx) * detStp;
				sinBeta = sin(beta);
				cosBeta = cos(beta);
				initDetX = S2D * sinBeta;
				initDetY = -S2D * cosBeta + S2O;
				curDetX = initDetX * cosT - initDetY * sinT;
				curDetY = initDetX * sinT + initDetY * cosT;
				dirx = curDetX - sourx;
				diry = curDetY - soury;
				rleng = sqrt(dirx * dirx + diry * diry);
				dirx /= rleng;
				diry /= rleng;

				for (unsigned int yIdx = 0; yIdx != YN; ++yIdx)
				{
					boxminy = (yIdx - objCntIdxy - 0.5) * ObjStpy;
					boxmaxy = (yIdx - objCntIdxy + 0.5) * ObjStpy;
					for (unsigned int xIdx = 0; xIdx != XN; ++xIdx)
					{
						boxminx = (xIdx - objCntIdxx - 0.5) * ObjStpx;
						boxmaxx = (xIdx - objCntIdxx + 0.5) * ObjStpx;

						intersected = intersectBox<double>(sourx, soury,
							dirx, diry, boxminx, boxminy,
							boxmaxx, boxmaxy, tnear, tfar);
						if (intersected)
						{
							rowIdx.push_back(angIdx * DN + detIdx);
							colIdx.push_back(yIdx * XN + xIdx);
							weight.push_back(tfar - tnear);
						}

					}
				}
			}
		}
	}


	void pushValuesCOO_upSample4(
		const int sarIdx,
		const int endIdx,
		const double S2O,
		const double S2D,
		const int DN,
		const int XN,
		const int YN,
		const double detCntIdx,
		const double objCntIdxy,
		const double objCntIdxx,
		const double ObjStpx,
		const double ObjStpy,
		const double detStp,
		const std::vector<double>& angs,
		std::vector<int>& rowIdx,
		std::vector<int>& colIdx,
		std::vector<double>& weight)
	{
		double cosT, sinT, sourx, soury, beta, sinBeta, cosBeta,
			initDetX, initDetY, curDetX, curDetY, dirx, diry,
			rleng, boxminx, boxminy, boxmaxx, boxmaxy, tnear, tfar;
		bool intersected;
		double updetStp = detStp / 8.0;
		double wweight = 0;

		for (unsigned int yIdx = 0; yIdx != YN; ++yIdx)
		{
			boxminy = (yIdx - objCntIdxy - 0.5) * ObjStpy;
			boxmaxy = (yIdx - objCntIdxy + 0.5) * ObjStpy;
			for (unsigned int xIdx = 0; xIdx != XN; ++xIdx)
			{
				boxminx = (xIdx - objCntIdxx - 0.5) * ObjStpx;
				boxmaxx = (xIdx - objCntIdxx + 0.5) * ObjStpx;
				for (int angIdx = sarIdx; angIdx != endIdx; ++angIdx) // Each angle;
				{
					cosT = cos(angs[angIdx]);
					sinT = sin(angs[angIdx]);
					sourx = -S2O * sinT;
					soury = S2O * cosT;
					for (unsigned int detIdx = 0; detIdx != DN; ++detIdx) // Each detector cell;
					{
						//intersected = false;
						wweight = 0;
						for (unsigned int subdetIdx = 0; subdetIdx != 8; ++subdetIdx)
						{
							beta = (detIdx - detCntIdx) * detStp + (subdetIdx - 3.5) * updetStp;
							sinBeta = sin(beta);
							cosBeta = cos(beta);
							initDetX = S2D * sinBeta;
							initDetY = -S2D * cosBeta + S2O;
							curDetX = initDetX * cosT - initDetY * sinT;
							curDetY = initDetX * sinT + initDetY * cosT;
							dirx = curDetX - sourx;
							diry = curDetY - soury;
							rleng = sqrt(dirx * dirx + diry * diry);
							dirx /= rleng;
							diry /= rleng;
							intersected = intersectBox<double>(sourx, soury,
								dirx, diry, boxminx, boxminy,
								boxmaxx, boxmaxy, tnear, tfar);
							if (intersected)
							{
								wweight += (tfar - tnear);
							}
						}
						if (wweight != 0)
						{
							rowIdx.push_back(angIdx * DN + detIdx);
							colIdx.push_back(yIdx * XN + xIdx);
							weight.push_back(wweight * 0.125);
						}
					}
				}
			}
		}
	}











	void genProj_SIDDON(std::vector<int>& rowIdx,
		std::vector<int>& colIdx,
		std::vector<double>& weight,
		const double S2O, const double O2D,
		const double detCntIdx, const double detStp,
		const int DN,
		const double ObjStpx, const double ObjStpy,
		const int XN, const int YN,
		const double objCntIdxx,
		const double objCntIdxy,
		const std::vector<double> angs)
	{
		const int PN = angs.size();
		double cosT, sinT;
		double sourx, soury;
		double beta = 0;
		double cosBeta, sinBeta;
		double initDetX, initDetY;
		const double S2D = S2O + O2D;
		double curDetX, curDetY;
		double boxminy, boxmaxy, boxminx, boxmaxx;
		double dirx, diry, rleng;
		bool intersected;
		double tnear, tfar;

		for (int angIdx = 0; angIdx != PN; ++angIdx)
		{
			cosT = cos(angs[angIdx]);
			sinT = sin(angs[angIdx]);
			sourx = -S2O * sinT;
			soury = S2O * cosT;

			std::cout << angIdx << std::endl;
			for (unsigned int detIdx = 0; detIdx != DN; ++detIdx)
			{
				beta = (detIdx - detCntIdx) * detStp;
				sinBeta = sin(beta);
				cosBeta = cos(beta);
				initDetX = S2D * sinBeta;
				initDetY = -S2D * cosBeta + S2O;
				curDetX = initDetX * cosT - initDetY * sinT;
				curDetY = initDetX * sinT + initDetY * cosT;
				dirx = curDetX - sourx;
				diry = curDetY - soury;
				rleng = sqrt(dirx * dirx + diry * diry);
				dirx /= rleng;
				diry /= rleng;

				for (unsigned int yIdx = 0; yIdx != YN; ++yIdx)
				{
					boxminy = (yIdx - objCntIdxy - 0.5) * ObjStpy;
					boxmaxy = (yIdx - objCntIdxy + 0.5) * ObjStpy;
					for (unsigned int xIdx = 0; xIdx != XN; ++xIdx)
					{
						boxminx = (xIdx - objCntIdxx - 0.5) * ObjStpx;
						boxmaxx = (xIdx - objCntIdxx + 0.5) * ObjStpx;

						intersected = intersectBox<double>(sourx, soury,
							dirx, diry, boxminx, boxminy,
							boxmaxx, boxmaxy, tnear, tfar);
						if (intersected)
						{
							rowIdx.push_back(angIdx * DN + detIdx);
							colIdx.push_back(yIdx * XN + xIdx);
							weight.push_back(tfar - tnear);
						}

					}
				}
			}
		}



	}




	void genProj_SIDDON_openMP(
		std::vector<int>& rowIdx,
		std::vector<int>& colIdx,
		std::vector<double>& weight,
		const double S2O, const double O2D,
		const double detCntIdx, const double detStp,
		const int DN,
		const double ObjStpx, const double ObjStpy,
		const int XN, const int YN,
		const double objCntIdxx,
		const double objCntIdxy,
		const std::vector<double> angs)
	{
		const int PN = angs.size();
		//double cosT, sinT;
		//double sourx, soury;
		//double beta = 0;
		//double cosBeta, sinBeta;
		//double initDetX, initDetY;
		const double S2D = S2O + O2D;
		//double curDetX, curDetY;
		//double boxminy, boxmaxy, boxminx, boxmaxx;
		//double dirx, diry, rleng;
		//bool intersected;
		//double tnear, tfar;


		int threadsNum = omp_get_num_procs();
		int vPN_perRang = PN / threadsNum;
		int *rang = new int[threadsNum + 1];
		rang[0] = 0;
		for (int i = 1; i <= threadsNum; ++i)
		{
			rang[i] = vPN_perRang * i;
		}
		rang[threadsNum] = PN;

		std::vector<int> *srowIdx = new std::vector<int>[threadsNum];
		std::vector<int> *scolIdx = new std::vector<int>[threadsNum];
		std::vector<double> *sweight = new std::vector<double>[threadsNum];


		//#pragma omp parallel for
		//	for (int threadidx = 0; threadidx < threadsNum; ++threadidx)
		//	{
		//		pushValuesCOO(rang[threadidx], rang[threadidx + 1], S2O, S2D,
		//			DN, XN, YN, detCntIdx, objCntIdxy, objCntIdxx, ObjStpx, ObjStpy,
		//			detStp, angs, srowIdx[threadidx], scolIdx[threadidx], sweight[threadidx]);
		//	}

#pragma omp parallel for
		for (int threadidx = 0; threadidx < threadsNum; ++threadidx)
		{
			pushValuesCOO_upSample4(rang[threadidx], rang[threadidx + 1], S2O, S2D,
				DN, XN, YN, detCntIdx, objCntIdxy, objCntIdxx, ObjStpx, ObjStpy,
				detStp, angs, srowIdx[threadidx], scolIdx[threadidx], sweight[threadidx]);
		}


		//Connect them together
		rowIdx.clear();
		colIdx.clear();
		weight.clear();

		for (int kk = 0; kk != threadsNum; ++kk)
		{
			rowIdx.insert(rowIdx.end(), srowIdx[kk].begin(), srowIdx[kk].end());
			colIdx.insert(colIdx.end(), scolIdx[kk].begin(), scolIdx[kk].end());
			weight.insert(weight.end(), sweight[kk].begin(), sweight[kk].end());
		}
	}


	void genProj_SIDDON(
		std::vector<int>& rowIdx,
		std::vector<int>& colIdx,
		std::vector<double>& weight,
		const FanEAGeo& FanGeo,
		const Image& Img,
		const unsigned int& sampNum)
	{
		double2 MINO = make_double2(
			-Img.m_Size.x / 2.0f + Img.m_Bias.x,
			-Img.m_Size.y / 2.0f + Img.m_Bias.y);

		double curAng = 0;
		double cosT = 0;
		double sinT = 0;
		Ray2Dd ray;

		//unsigned int detId;
		double ang(0); //bias angle from the center of the Fan Beams

		double2 curDetPos; //the current detector element position;
		//double totWeight(0);

		double smallDetStep = FanGeo.m_DetStp / sampNum; //\CF²\C9\D1\F9\BA\F3\B5\C4detector\B2\BD\B3\A4;
		double cntSmallDet = sampNum * 0.5f;
		double realAng = 0;
		unsigned int angIdx = 0;
		unsigned int detIdx = 0;
		unsigned int subDetIdx = 0;

		int rowidx = 0;
		for (angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
		{
			//Current rotation angle;
			curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
			cosT = cos(curAng);
			sinT = sin(curAng);
			ray.o = rotation(make_double2(0, FanGeo.m_S2O), cosT, sinT);

			for (detIdx = 0; detIdx != FanGeo.m_DetN; ++detIdx)
			{

				rowidx = angIdx * FanGeo.m_DetN + detIdx;

				//Calculate current Angle
				ang = ((double) detIdx - FanGeo.m_DetCntIdx + 0.5f) * FanGeo.m_DetStp;

				for (subDetIdx = 0; subDetIdx != sampNum; ++subDetIdx)
				{
					//correct ang
					realAng = ang + (static_cast<double>(subDetIdx) -cntSmallDet + 0.5f) *smallDetStep;
					// current detector element position;
					curDetPos = rotation(make_double2(sin(realAng) * FanGeo.m_S2D, -cos(realAng) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);

					// X-ray direction
					ray.d = normalize(curDetPos - ray.o);

					pushMatrix(ray.o.x, ray.o.y,
						curDetPos.x, curDetPos.y,
						MINO.x, MINO.y,
						Img.m_Step.x, Img.m_Step.y,
						Img.m_Reso.x, Img.m_Reso.y,
						rowidx, rowIdx, colIdx, weight);
				}
			}
		}

	}




	// AIM based projection model
	template<typename T>
	inline void SortProjection(T(&Grid)[4][3])
	{
		int i, j;
		T td;
		//mexPrintf("S0=%d,S1=%d,S2=%d,S3=%d\n",SSort[0],SSort[1],SSort[2],SSort[3]);
		for (i = 0; i < 3; i++)
		{
			for (j = i + 1; j < 4; j++)
			{
				if (Grid[j][2] < Grid[i][2])
				{
					td = Grid[i][0];
					Grid[i][0] = Grid[j][0];
					Grid[j][0] = td;

					td = Grid[i][1];
					Grid[i][1] = Grid[j][1];
					Grid[j][1] = td;

					td = Grid[i][2];
					Grid[i][2] = Grid[j][2];
					Grid[j][2] = td;
				}
			}
		}
	}

	template<typename T>
	inline T ComputeCoefficient(const T Grid[4][3], const T SVA[3], const T SVB[3], const T SPoint[2], const T area)
	{
		T coef = 0;
		T x0, y0, a, b, t;
		int AI, BI;

		///compute the value of AI,BI
		if (SVA[2] < Grid[0][2])
			AI = 0;
		else if (SVA[2] < Grid[1][2])
			AI = 1;
		else if (SVA[2] < Grid[2][2])
			AI = 2;
		else AI = 3;

		if (SVB[2] < Grid[1][2])
			BI = 1;
		else if (SVB[2] < Grid[2][2])
			BI = 2;
		else if (SVB[2] < Grid[3][2])
			BI = 3;
		else BI = 4;

		//mexPrintf("AI=%d,PA=%f;BI=%d,PB=%f\n",AI,Grid[AI][2],BI,Grid[BI-1][2]);

		switch (AI)
		{
		case 0:
		{
			switch (BI)
			{
			case 1:// case [0,1]
			{
				x0 = Grid[0][0] - SPoint[0];
				y0 = Grid[0][1] - SPoint[1];
				if (abs(SVB[0] * SVB[1]) > 0)
				{
					t = x0*SVB[1] - y0*SVB[0];
					a = t / SVB[0];
					b = t / SVB[1];
					coef = 0.5 * abs(a*b);
				}
				break;
			}
			case 2: // case [0,2]
			{
				x0 = abs(Grid[0][0] - Grid[1][0]);
				y0 = abs(Grid[0][1] - Grid[1][1]);
				if (x0 > y0) // line is on the x-dirction
				{
					a = abs((Grid[0][0] - SPoint[0])*SVB[1] / SVB[0] - (Grid[0][1] - SPoint[1]));
					b = abs((Grid[1][0] - SPoint[0])*SVB[1] / SVB[0] - (Grid[1][1] - SPoint[1]));
					coef = (a + b) * x0 * 0.5;
				}
				else
				{
					a = abs((Grid[0][0] - SPoint[0]) - (Grid[0][1] - SPoint[1])*SVB[0] / SVB[1]);
					b = abs((Grid[1][0] - SPoint[0]) - (Grid[1][1] - SPoint[1])*SVB[0] / SVB[1]);
					coef = (a + b)*y0*0.5;
				}
				break;
			}
			case 3://case [0,3]
			{
				x0 = Grid[3][0] - SPoint[0];
				y0 = Grid[3][1] - SPoint[1];
				if (abs(SVB[0] * SVB[1]) > 0)
				{
					t = x0*SVB[1] - y0*SVB[0];
					a = t / SVB[0];
					b = t / SVB[1];
					coef = 0.5 * abs(a*b);
					coef = area - coef;
				}
				else
					coef = area;

				break;
			}
			case 4: // case [0,4]
			{
				coef = area;
				break;
			}
			default: break;
			}
			break;
		}//end case 0 of AI
		case 1:
		{
			switch (BI)
			{
			case 1://case [1,1]
			{
				x0 = Grid[0][0] - SPoint[0];
				y0 = Grid[0][1] - SPoint[1];
				t = x0*SVB[1] - y0*SVB[0];
				if (abs(SVB[0] * SVB[1]) > 0)
				{
					a = t / SVB[0];
					b = t / SVB[1];
					coef = 0.5*abs(a*b);
				}
				t = x0*SVA[1] - y0*SVA[0];
				if (abs(SVA[0] * SVA[1]) > 0)
				{
					a = t / SVA[0];
					b = t / SVA[1];
					coef = abs(coef - 0.5*abs(a*b));
				}
				break;
			}
			case 2://case [1,2]
			{
				x0 = abs(Grid[0][0] - Grid[1][0]);
				y0 = abs(Grid[0][1] - Grid[1][1]);
				if (x0 > y0) // line is on the x-dirction
				{
					a = abs((Grid[0][0] - SPoint[0])*SVB[1] / SVB[0] - (Grid[0][1] - SPoint[1]));
					b = abs((Grid[1][0] - SPoint[0])*SVB[1] / SVB[0] - (Grid[1][1] - SPoint[1]));
					coef = (a + b)*x0*0.5;
				}
				else
				{
					a = abs((Grid[0][0] - SPoint[0]) - (Grid[0][1] - SPoint[1])*SVB[0] / SVB[1]);
					b = abs((Grid[1][0] - SPoint[0]) - (Grid[1][1] - SPoint[1])*SVB[0] / SVB[1]);
					coef = (a + b)*y0*0.5;
				}
				x0 = Grid[0][0] - SPoint[0];
				y0 = Grid[0][1] - SPoint[1];
				if (abs(SVA[0] * SVA[1]) > 0)
				{
					t = x0*SVA[1] - y0*SVA[0];
					a = t / SVA[0];
					b = t / SVA[1];
					coef = abs(0.5*abs(a*b) - coef);
				}
				break;
			}
			case 3://case [1,3]
			{
				x0 = Grid[0][0] - SPoint[0];
				y0 = Grid[0][1] - SPoint[1];
				if (abs(SVA[0] * SVA[1]) > 0)
				{
					t = x0*SVA[1] - y0*SVA[0];
					a = t / SVA[0];
					b = t / SVA[1];
					coef = area - 0.5*abs(a*b);
				}
				else
					coef = area;
				x0 = Grid[3][0] - SPoint[0];
				y0 = Grid[3][1] - SPoint[1];
				if (abs(SVB[0] * SVB[1]) > 0)
				{
					t = x0*SVB[1] - y0*SVB[0];
					a = t / SVB[0];
					b = t / SVB[1];
					coef = coef - 0.5*abs(a*b);
				}
				break;
			}
			case 4://case [1,4]
			{
				x0 = Grid[0][0] - SPoint[0];
				y0 = Grid[0][1] - SPoint[1];
				if (abs(SVA[0] * SVA[1]) > 0)
				{
					t = x0*SVA[1] - y0*SVA[0];
					a = t / SVA[0];
					b = t / SVA[1];
					coef = 0.5*abs(a*b);
					coef = area - coef;
				}
				else
					coef = area;
				break;
			}
			default: break;
			}
			break;
		}//end case 1 of AI
		case 2:
		{
			switch (BI)
			{
			case 2:
			{
				x0 = abs(Grid[0][0] - Grid[1][0]);
				y0 = abs(Grid[0][1] - Grid[1][1]);
				if (x0 > y0) // line is on the x-dirction
				{
					a = abs((Grid[0][0] - SPoint[0])*SVB[1] / SVB[0] - (Grid[0][1] - SPoint[1]));
					b = abs((Grid[1][0] - SPoint[0])*SVB[1] / SVB[0] - (Grid[1][1] - SPoint[1]));
					coef = (a + b)*x0*0.5;
					a = abs((Grid[0][0] - SPoint[0])*SVA[1] / SVA[0] - (Grid[0][1] - SPoint[1]));
					b = abs((Grid[1][0] - SPoint[0])*SVA[1] / SVA[0] - (Grid[1][1] - SPoint[1]));
					coef = abs(coef - (a + b)*x0*0.5);
				}
				else
				{
					a = abs((Grid[0][0] - SPoint[0]) - (Grid[0][1] - SPoint[1])*SVB[0] / SVB[1]);
					b = abs((Grid[1][0] - SPoint[0]) - (Grid[1][1] - SPoint[1])*SVB[0] / SVB[1]);
					coef = (a + b)*y0*0.5;
					a = abs((Grid[0][0] - SPoint[0]) - (Grid[0][1] - SPoint[1])*SVA[0] / SVA[1]);
					b = abs((Grid[1][0] - SPoint[0]) - (Grid[1][1] - SPoint[1])*SVA[0] / SVA[1]);
					coef = abs(coef - (a + b)*y0*0.5);
				}
				break;
			}
			case 3:
			{
				x0 = abs(Grid[2][0] - Grid[3][0]);
				y0 = abs(Grid[2][1] - Grid[3][1]);
				if (x0 > y0) // line is on the x-dirction
				{
					a = abs((Grid[2][0] - SPoint[0])*SVA[1] / SVA[0] - (Grid[2][1] - SPoint[1]));
					b = abs((Grid[3][0] - SPoint[0])*SVA[1] / SVA[0] - (Grid[3][1] - SPoint[1]));
					coef = (a + b)*x0*0.5;
				}
				else
				{
					a = abs((Grid[2][0] - SPoint[0]) - (Grid[2][1] - SPoint[1])*SVA[0] / SVA[1]);
					b = abs((Grid[3][0] - SPoint[0]) - (Grid[3][1] - SPoint[1])*SVA[0] / SVA[1]);
					coef = (a + b)*y0*0.5;
				}
				x0 = Grid[3][0] - SPoint[0];
				y0 = Grid[3][1] - SPoint[1];
				if (abs(SVB[0] * SVB[1]) > 0)
				{
					t = x0*SVB[1] - y0*SVB[0];
					a = t / SVB[0];
					b = t / SVB[1];
					coef = abs(0.5 * abs(a*b) - coef);
				}
				break;
			}
			case 4:
			{
				x0 = abs(Grid[2][0] - Grid[3][0]);
				y0 = abs(Grid[2][1] - Grid[3][1]);
				if (x0 > y0) // line is on the x-dirction
				{
					a = abs((Grid[2][0] - SPoint[0])*SVA[1] / SVA[0] - (Grid[2][1] - SPoint[1]));
					b = abs((Grid[3][0] - SPoint[0])*SVA[1] / SVA[0] - (Grid[3][1] - SPoint[1]));
					coef = (a + b)*x0*0.5;
				}
				else // line is on the y-direction
				{
					a = abs((Grid[2][0] - SPoint[0]) - (Grid[2][1] - SPoint[1])*SVA[0] / SVA[1]);
					b = abs((Grid[3][0] - SPoint[0]) - (Grid[3][1] - SPoint[1])*SVA[0] / SVA[1]);
					coef = (a + b)*y0*0.5;
				}
				break;
			}
			default: break;
			}
			break;
		}//end case 2 of AI
		case 3:
		{
			switch (BI)
			{
			case 3:
			{
				x0 = Grid[3][0] - SPoint[0];
				y0 = Grid[3][1] - SPoint[1];
				if (abs(SVB[0] * SVB[1]) > 0)
				{
					t = x0*SVB[1] - y0*SVB[0];
					a = t / SVB[0];
					b = t / SVB[1];
					coef = 0.5 * abs(a*b);
				}
				if (abs(SVA[0] * SVA[1]) > 0)
				{
					t = x0*SVA[1] - y0*SVA[0];
					a = t / SVA[0];
					b = t / SVA[1];
					coef = abs(coef - 0.5 * abs(a*b));
				}
				break;
			}
			case 4:
			{
				x0 = Grid[3][0] - SPoint[0];
				y0 = Grid[3][1] - SPoint[1];
				if (abs(SVA[0] * (SVA[1])) > 0)
				{
					t = x0*SVA[1] - y0*SVA[0];
					a = t / SVA[0];
					b = t / SVA[1];
					coef = 0.5*abs(a*b);

				}
				break;
			}
			default: break;
			}
			break;
		}//end case 3 of AI
		}//end of switch AI
		return coef;
	}



	template<typename T>
	void genProj_AIM(std::vector<int>& rowIdx, std::vector<int>& colIdx, std::vector<T>& coeffs, int angIdx, const T ang, const FanEAGeo FanGeo, const Image Img)
	{
		const T ScanR = FanGeo.m_S2O;
		//const T ObjR = Img.m_Size.x;
		//	const int PN = FanGeo.m_ViwN;
		const int DN = FanGeo.m_DetN;
		const int XN = Img.m_Reso.x;
		const int YN = Img.m_Reso.y;
		const T dx = Img.m_Step.x;
		const T dy = Img.m_Step.y;
		const T area = dx * dy;
		const T xctr = XN * 0.5;
		const T yctr = YN * 0.5;
		const T dctr = FanGeo.m_DetCntIdx;
		const T DtBeta = FanGeo.m_DetStp;
		T coef;
		T* PosAry = new T[(XN + 1) * (YN + 1)];
		T* PBeta = new T[DN];
		for (int i = 0; i != DN; ++i)
		{
			PBeta[i] = (i - dctr + 0.5) * DtBeta;
		}
		T Ew[2], Ed[2];
		const T cosAng = std::cos(ang);
		const T sinAng = std::sin(ang);
		Ew[0] = -cosAng;
		Ew[1] = -sinAng;
		Ed[0] = -sinAng;
		Ed[1] = cosAng;
		T SPoint[2];
		SPoint[0] = ScanR * cosAng;
		SPoint[1] = ScanR * sinAng;
		int xi(0), yi(0);

		T xcor, ycor, dcor;
		T pdist;
		int di;
		for (yi = 0; yi <= YN; ++yi)
		{
			ycor = (yi - yctr)*dy - SPoint[1];
			for (xi = 0; xi <= XN; ++xi)
			{
				xcor = (xi - xctr)*dx - SPoint[0];
				dcor = (xcor*Ed[0] + ycor*Ed[1]) / (xcor*Ew[0] + ycor*Ew[1]);
				dcor = (std::atan(dcor) - PBeta[0]) / DtBeta;
				PosAry[yi*(XN + 1) + xi] = dcor + 0.5;
			}
		}

		T Grid[4][3];
		int posim(0);
		int MinBV(0);
		int MaxBV(0);
		T pangle, temp, SVA[3], SVB[3];
		for (yi = 0; yi < YN; yi++)
		{
			for (xi = 0; xi < XN; xi++)
			{
				//Fetch the four points of the pixel and their projection positions
				Grid[0][0] = (xi - xctr)*dx;
				Grid[0][1] = (yi - yctr)*dy;
				Grid[0][2] = PosAry[yi*(XN + 1) + xi];
				Grid[1][0] = (xi - xctr + 1)*dx;
				Grid[1][1] = (yi - yctr)*dy;
				Grid[1][2] = PosAry[yi*(XN + 1) + xi + 1];
				Grid[2][0] = (xi - xctr + 1)*dx;
				Grid[2][1] = (yi - yctr + 1)*dy;
				Grid[2][2] = PosAry[(yi + 1)*(XN + 1) + xi + 1];
				Grid[3][0] = (xi - xctr)*dx;
				Grid[3][1] = (yi - yctr + 1)*dy;
				Grid[3][2] = PosAry[(yi + 1)*(XN + 1) + xi];
				SortProjection<T>(Grid);//Sort the projection psotion

				posim = yi*XN + xi;

				//pvalue = Img[posim];
				pdist = std::hypot((xi + 0.5 - xctr)*dx - SPoint[0], (yi + 0.5 - yctr)*dy - SPoint[1]); // sqrt(pow(, 2) + pow(, 2));

				//Computer the weighting coefficient for every projection position
				MinBV = int(Grid[0][2] + 10) - 10;
				MaxBV = int(Grid[3][2] + 10) - 9;
				if (MinBV < 0)   MinBV = 0;
				if (MaxBV > DN)  MaxBV = DN;
				//double total =0;
				for (di = MinBV; di < MaxBV; di++)
				{
					// Compute the directions of the two lines for the projections
					pangle = PBeta[di] - 0.5 * DtBeta;
					temp = ang + PI - pangle;
					SVA[0] = std::cos(temp);
					SVA[1] = std::sin(temp);
					SVA[2] = di;
					// mexPrintf("di=%d,VA0=%10.8f,VA1=%10.8f,angle=%10.8f,Beta=%10.8f\n",di,SVA[0],SVA[1],temp,pangle);
					pangle = PBeta[di] + 0.5*DtBeta;
					temp = ang + PI - pangle;
					SVB[0] = std::cos(temp);
					SVB[1] = std::sin(temp);
					SVB[2] = di + 1;

					//compute the weighting coefficient for a special projection data
					coef = ComputeCoefficient<T>(Grid, SVA, SVB, SPoint, area);
					coef = coef / (pdist * std::abs(DtBeta));
					rowIdx.push_back(angIdx * DN + di);
					colIdx.push_back(posim);
					coeffs.push_back(coef);

				}
			}
		}

		delete [] PosAry;
		delete [] PBeta;
	}


	void genMatrix_EA_AIM()
	{
		FanEAGeo FanGeo;
		Image Img;
		std::ifstream fin("configurationFile.txt", std::ios::in);
		if (!fin.is_open())
		{
			std::cout << "Cannot Find configuration file\n";
			std::cout << "Would you like to use the default configuration?[y/n]";
			char yes;
			std::cin >> yes;
			if (yes == 'y' || yes == 'Y')
			{
				FanGeo.m_S2O = 5.385200195312500e+02;
				FanGeo.m_O2D = 4.082259521484375e+02;
				FanGeo.m_ViwN = 360;
				FanGeo.m_ViwBeg = 0;
				FanGeo.m_ViwStp = 3.14159265358979323846264*2.0 / 359.0;
				FanGeo.m_DetArc = 0.95928517242269f;
				FanGeo.m_DetCntIdx = 222;
				FanGeo.m_DetN = 444;
				FanGeo.m_DetStp = 0.95928517242269f / 444;
				FanGeo.m_S2D = 946.7459716796875;
				Img.m_Size.x = 4.484740011196460e+02;
				Img.m_Size.y = 4.484740011196460e+02;
				Img.m_Reso.x = 256;
				Img.m_Reso.y = 256;
				Img.m_Bias.x = 0;
				Img.m_Bias.y = 0;
				Img.m_Step.x = Img.m_Size.x / Img.m_Reso.x;
				Img.m_Step.y = Img.m_Size.y / Img.m_Reso.y;
				int angIdx;
				std::vector<int> rowIdx;
				std::vector<int> colIdx;
				std::vector<double> coeffs;

				std::stringstream nonZ;

				for (angIdx = 0; angIdx < FanGeo.m_ViwN; angIdx++)
				{
					double ang = FanGeo.m_ViwBeg + angIdx* FanGeo.m_ViwStp;

					genProj_AIM<double>(rowIdx, colIdx, coeffs, angIdx, ang, FanGeo, Img);
					std::cout << angIdx << std::endl;
				}
				nonZ << rowIdx.size();
				std::string F1 = "prjAIM" + nonZ.str() + ".row";
				std::string F2 = "prjAIM" + nonZ.str() + ".col";
				std::string F3 = "prjAIM" + nonZ.str() + ".cof";
				std::ofstream rowFile(F1.c_str(), std::ios::binary);
				std::ofstream colFile(F2.c_str(), std::ios::binary);
				std::ofstream coeFile(F3.c_str(), std::ios::binary);

				rowFile.write((char*) &(rowIdx[0]), sizeof(int) * rowIdx.size());
				rowFile.close();
				colFile.write((char*) &(colIdx[0]), sizeof(int) * colIdx.size());
				colFile.close();
				coeFile.write((char*) &(coeffs[0]), sizeof(double) * coeffs.size());
				coeFile.close();
				std::cout << coeffs.size() << std::endl;
			}
			else
			{
				std::cout << "Exit without generating\n";
				exit(-1);
			}
		}
		std::string s;
		double val[14];
		int i = 0;
		while (fin >> s)
		{
			if ((i % 2))
			{
				std::stringstream ss;
				ss << s;
				ss >> val[i / 2];
				std::cout << val[i / 2] << std::endl;
			}
			++i;
		}
		FanGeo.m_S2O = val[0];
		FanGeo.m_O2D = val[1];
		FanGeo.m_ViwN = val[2];
		FanGeo.m_ViwBeg = val[3];
		FanGeo.m_ViwStp = val[4];
		FanGeo.m_DetArc = val[5];
		FanGeo.m_DetN = val[6];
		FanGeo.m_DetCntIdx = val[7];

		FanGeo.m_DetStp = val[5] / val[6];
		FanGeo.m_S2D = val[0] + val[1];
		Img.m_Size.x = val[8];
		Img.m_Size.y = val[9];
		Img.m_Reso.x = val[10];
		Img.m_Reso.y = val[11];
		Img.m_Bias.x = val[12];
		Img.m_Bias.y = val[13];
		Img.m_Step.x = Img.m_Size.x / Img.m_Reso.x;
		Img.m_Step.y = Img.m_Size.y / Img.m_Reso.y;

		int angIdx;
		std::vector<int> rowIdx;
		std::vector <int> colIdx;
		std::vector <double> coeffs;

		std::stringstream nonZ;

		for (angIdx = 0; angIdx < FanGeo.m_ViwN; angIdx++)
		{
			double ang = FanGeo.m_ViwBeg + angIdx* FanGeo.m_ViwStp;

			genProj_AIM<double>(rowIdx, colIdx, coeffs, angIdx, ang, FanGeo, Img);
			std::cout << angIdx << std::endl;
		}
		nonZ << rowIdx.size();
		std::string F1 = "prjAIM" + nonZ.str() + ".row";
		std::string F2 = "prjAIM" + nonZ.str() + ".col";
		std::string F3 = "prjAIM" + nonZ.str() + ".cof";
		std::ofstream rowFile(F1.c_str(), std::ios::binary);
		std::ofstream colFile(F2.c_str(), std::ios::binary);
		std::ofstream coeFile(F3.c_str(), std::ios::binary);

		rowFile.write((char*) &(rowIdx[0]), sizeof(int) * rowIdx.size());
		rowFile.close();
		colFile.write((char*) &(colIdx[0]), sizeof(int) * colIdx.size());
		colFile.close();
		coeFile.write((char*) &(coeffs[0]), sizeof(double) * coeffs.size());
		coeFile.close();
		std::cout << coeffs.size() << std::endl;
	}

	void genMatrix_EA_SIDDON()
	{
		FanEAGeo FanGeo;
		Image Img;
		//¶ÁÈ¡ÅäÖÃÎÄŒþ;
		std::ifstream fin("configurationFile.txt", std::ios::in);
		if (!fin.is_open())
		{
			std::cout << "Cannot Find configuration file\n";
			std::cout << "Would you like to use the default configuration?[y/n]";
			char yes;
			std::cin >> yes;
			if (yes == 'y' || yes == 'Y')
			{
				FanGeo.m_S2O = 5.385200195312500e+02;
				FanGeo.m_O2D = 4.082259521484375e+02;
				FanGeo.m_ViwN = 360;
				FanGeo.m_ViwBeg = 0;
				FanGeo.m_ViwStp = 3.14159265358979323846264*2.0 / 359.0;
				FanGeo.m_DetArc = 0.95928517242269f;
				FanGeo.m_DetCntIdx = 222;
				FanGeo.m_DetN = 444;
				FanGeo.m_DetStp = 0.95928517242269f / 444;
				FanGeo.m_S2D = 946.7459716796875;
				Img.m_Size.x = 4.484740011196460e+02;
				Img.m_Size.y = 4.484740011196460e+02;
				Img.m_Reso.x = 256;
				Img.m_Reso.y = 256;
				Img.m_Bias.x = 0;
				Img.m_Bias.y = 0;
				Img.m_Step.x = Img.m_Size.x / Img.m_Reso.x;
				Img.m_Step.y = Img.m_Size.y / Img.m_Reso.y;
				//int angIdx;
				std::vector<int> rowIdx;
				std::vector <int> colIdx;
				std::vector <double> coeffs;

				std::stringstream nonZ;

				genProj_SIDDON(rowIdx, colIdx, coeffs, FanGeo, Img, 1);

				nonZ << rowIdx.size();
				std::string F1 = "prjSIDDON" + nonZ.str() + ".row";
				std::string F2 = "prjSIDDON" + nonZ.str() + ".col";
				std::string F3 = "prjSIDDON" + nonZ.str() + ".cof";
				std::ofstream rowFile(F1.c_str(), std::ios::binary);
				std::ofstream colFile(F2.c_str(), std::ios::binary);
				std::ofstream coeFile(F3.c_str(), std::ios::binary);

				rowFile.write((char*) &(rowIdx[0]), sizeof(int) * rowIdx.size());
				rowFile.close();
				colFile.write((char*) &(colIdx[0]), sizeof(int) * colIdx.size());
				colFile.close();
				coeFile.write((char*) &(coeffs[0]), sizeof(double) * coeffs.size());
				coeFile.close();
				std::cout << coeffs.size() << std::endl;
			}
			else
			{
				std::cout << "Exit without generating\n";
				exit(-1);
			}
		}
		std::string s;
		double val[14];
		int i = 0;
		while (fin >> s)
		{
			if ((i % 2))
			{
				std::stringstream ss;
				ss << s;
				ss >> val[i / 2];
				std::cout << val[i / 2] << std::endl;
			}
			++i;
		}
		FanGeo.m_S2O = val[0];
		FanGeo.m_O2D = val[1];
		FanGeo.m_ViwN = val[2];
		FanGeo.m_ViwBeg = val[3];
		FanGeo.m_ViwStp = val[4];
		FanGeo.m_DetArc = val[5];
		FanGeo.m_DetN = val[6];
		FanGeo.m_DetCntIdx = val[7];

		FanGeo.m_DetStp = val[5] / val[6];
		FanGeo.m_S2D = val[0] + val[1];
		Img.m_Size.x = val[8];
		Img.m_Size.y = val[9];
		Img.m_Reso.x = val[10];
		Img.m_Reso.y = val[11];
		Img.m_Bias.x = val[12];
		Img.m_Bias.y = val[13];
		Img.m_Step.x = Img.m_Size.x / Img.m_Reso.x;
		Img.m_Step.y = Img.m_Size.y / Img.m_Reso.y;

		//int angIdx;
		std::vector<int> rowIdx;
		std::vector <int> colIdx;
		std::vector <double> coeffs;

		std::stringstream nonZ;
		double S2O = val[0];
		double O2D = val[1];
		double detCntIdx = val[7];
		double detStp = val[5] / val[6];
		int DN = val[6];
		double ObjStpx = val[8] / val[10];
		double ObjStpy = val[9] / val[11];
		int XN = val[10];
		int YN = val[11];
		double objCntIdxx = (XN - 1.0) * 0.5;
		double objCntIdxy = (YN - 1.0) * 0.5;
		int PN = val[2];
		std::vector<double> angs(PN, 0);
		for (unsigned int i = 0; i != PN; ++i)
		{
			angs[i] = val[3] + i * val[4];
		}
		//genProj_SIDDON(rowIdx, colIdx, coeffs, FanGeo, Img, 1);
		genProj_SIDDON_openMP(rowIdx, colIdx, coeffs, S2O, O2D, detCntIdx, detStp, DN, ObjStpx, ObjStpy,
			XN, YN, objCntIdxx, objCntIdxy, angs);

		nonZ << rowIdx.size();
		std::string F1 = "prjSIDDON" + nonZ.str() + ".row";
		std::string F2 = "prjSIDDON" + nonZ.str() + ".col";
		std::string F3 = "prjSIDDON" + nonZ.str() + ".cof";
		std::ofstream rowFile(F1.c_str(), std::ios::binary);
		std::ofstream colFile(F2.c_str(), std::ios::binary);
		std::ofstream coeFile(F3.c_str(), std::ios::binary);

		rowFile.write((char*) &(rowIdx[0]), sizeof(int) * rowIdx.size());
		rowFile.close();
		colFile.write((char*) &(colIdx[0]), sizeof(int) * colIdx.size());
		colFile.close();
		coeFile.write((char*) &(coeffs[0]), sizeof(double) * coeffs.size());
		coeFile.close();
		std::cout << coeffs.size() << std::endl;

	}

};


