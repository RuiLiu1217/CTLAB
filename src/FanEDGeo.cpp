#include "FanEDGeo.h"


FanEDGeo::FanEDGeo(void)
{
	m_S2O = 85.0f;
	m_O2D = 15.0f;
	m_S2D = 100.0f;

	m_ViwN = 360; // view number
	m_ViwBeg = 0.0f;
	m_ViwStp = 3.14159265358979323846264f * 2.0f / m_ViwN;

	m_DetSize = 30.0f; // Detector size;
	m_DetN = 888; // How many detector cells
	m_DetStp = m_DetSize / m_DetN; // detector cells size;
	m_DetCntIdx = m_DetN * 0.5f - 0.5f; //detector center index;
}

FanEDGeo::FanEDGeo(const FanEDGeo& rhs)
{
	m_S2O = rhs.m_S2O;
	m_O2D = rhs.m_O2D;
	m_S2D = rhs.m_S2D;

	m_ViwN = rhs.m_ViwN;
	m_ViwBeg = rhs.m_ViwBeg;
	m_ViwStp = rhs.m_ViwStp;

	m_DetSize = rhs.m_DetSize;
	m_DetN = rhs.m_DetN;
	m_DetStp = rhs.m_DetStp;
	m_DetCntIdx = rhs.m_DetCntIdx;
}


FanEDGeo::FanEDGeo(const float S2O, const float O2D, const unsigned int ViwN,
	const float ViwBeg, const float ViwEnd, const float DetSize,
	const unsigned int DetN) :m_S2O(S2O), m_O2D(O2D), m_S2D(S2O + O2D),
	m_ViwN(ViwN), m_ViwBeg(ViwBeg), m_ViwStp((ViwEnd - ViwBeg) / ViwN),
	m_DetSize(DetSize), m_DetN(DetN), m_DetStp(DetSize / DetN),
	m_DetCntIdx(DetN * 0.5f - 0.5f) {}
