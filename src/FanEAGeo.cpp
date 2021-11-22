#include "FanEAGeo.h"

FanEAGeo::FanEAGeo(void)
{
	m_DetArc = 0.95928517242269f;
	m_DetN = 888;
	m_DetCntIdx = 444.75;
	m_DetStp = m_DetArc / m_DetN;
	m_O2D = 4.082259521484375e+02;

	m_S2O = 5.385200195312500e+02;
	m_S2D = m_S2O + m_O2D;
	m_ViwBeg = 0.0f;// (-2.668082275390625e+02 / 180.0* 3.14159265358979323846264);
	m_ViwN = 360;
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


FanEAGeo::FanEAGeo(const float S2O, const float O2D, const unsigned int ViwN,
	const float ViwBeg, const float ViwEnd, const float DetArc,
	const unsigned int DetN) :m_S2O(S2O), m_O2D(O2D), m_S2D(S2O + O2D),
	m_ViwN(ViwN), m_ViwBeg(ViwBeg), m_ViwStp((ViwEnd - ViwBeg) / float(ViwN)),
	m_DetArc(DetArc), m_DetN(DetN), m_DetStp(DetArc / DetN),
	m_DetCntIdx(DetN * 0.5f - 0.5f) {}