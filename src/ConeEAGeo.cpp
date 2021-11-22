#include "ConeEAGeo.h"

ConeEAGeo::ConeEAGeo(void) {

	m_DetCntIdx.x = 444.75;

	m_ViwBeg = 0.0f;// (-2.668082275390625e+02 / 180.0* 3.14159265358979323846264);


	m_S2O = 5.385200195312500e+02;
	m_O2D = 4.082259521484375e+02;
	m_S2D = m_S2O + m_O2D;

	m_ViwN = 2200;
	m_ViwBeg = 0.0f;
	m_ViwStp = 3.14159265358979323846264f * 2.0f / m_ViwN;

	m_DetArc = 0.95928517242269f;
	m_DetN = 888;
	m_DetStp = m_DetArc / m_DetN;

	m_DetHei = 64.0f;
	m_DetHN = 64;
	m_DetHStp = m_DetHei / m_DetHN;

	m_DetCntIdx.x = m_DetN * 0.5f - 0.5f;
	m_DetCntIdx.y = m_DetN * 0.5f - 0.5f;
}

ConeEAGeo::ConeEAGeo(const ConeEAGeo& rhs)
{
	m_DetCntIdx = rhs.m_DetCntIdx;
	m_DetArc = rhs.m_DetArc;
	m_DetHei = rhs.m_DetHei;
	m_DetHN = rhs.m_DetHN;
	m_DetHStp = rhs.m_DetHStp;
	m_DetN = rhs.m_DetN;
	m_DetStp = rhs.m_DetStp;
	m_O2D = rhs.m_O2D;
	m_S2D = rhs.m_S2D;
	m_S2O = rhs.m_S2O;
	m_ViwBeg = rhs.m_ViwBeg;
	m_ViwN = rhs.m_ViwN;
	m_ViwStp = rhs.m_ViwStp;
}

ConeEAGeo::ConeEAGeo(const float S2O, const float O2D, const unsigned int viwN,
	const float ViwBeg, const float ViwEnd, const float DetArc, const unsigned int DetN,
	const float DetHei, const unsigned int DetHN) :m_S2O(S2O), m_O2D(O2D),
	m_S2D(S2O + O2D), m_ViwN(viwN), m_ViwBeg(ViwBeg), m_ViwStp((ViwEnd - ViwBeg) / viwN),
	m_DetArc(DetArc), m_DetN(DetN), m_DetStp(DetArc / DetN), m_DetHei(DetHei),
	m_DetHN(DetHN), m_DetHStp(DetHei / DetHN)
{
	m_DetCntIdx.x = m_DetN * 0.5f - 0.5f;
	m_DetCntIdx.y = m_DetHN * 0.5f - 0.5f;
}