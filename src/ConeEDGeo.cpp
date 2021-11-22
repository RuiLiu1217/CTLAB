#include "ConeEDGeo.h"


ConeEDGeo::ConeEDGeo()
{
	m_S2O = 85.0f;
	m_O2D = 15.0f;
	m_S2D = m_S2O + m_O2D;
	m_ViwN = 360;
	m_ViwBeg = 0.0f;
	m_ViwStp = 3.14159265358979323846264f * 2.0f / m_ViwN;

	m_DetSize = make_float2(34.0f, 17.0f);
	m_DetN = make_int2(1024, 512);
	m_DetStp = make_float2(m_DetSize.x / m_DetN.x, m_DetSize.y / m_DetN.y);
	m_DetCntIdx = make_float2(m_DetN.x * 0.5f - 0.5f, m_DetN.y * 0.5f - 0.5f);
}

ConeEDGeo::ConeEDGeo(const ConeEDGeo& rhs)
{
	m_DetCntIdx = rhs.m_DetCntIdx;
	m_DetN = rhs.m_DetN;
	m_DetSize = rhs.m_DetSize;
	m_DetStp = rhs.m_DetStp;
	m_O2D = rhs.m_O2D;
	m_S2D = rhs.m_S2D;
	m_S2O = rhs.m_S2O;
	m_ViwBeg = rhs.m_ViwBeg;
	m_ViwN = rhs.m_ViwN;
	m_ViwStp = rhs.m_ViwStp;
}


ConeEDGeo::ConeEDGeo(const float S2O, const float O2D, const  int ViwN,
	const float ViwBeg, const float ViwEnd, const float2 DetSize,
	const int2 DetN) :m_S2O(S2O), m_O2D(O2D), m_S2D(S2O + O2D),
	m_ViwN(ViwN), m_ViwBeg(ViwBeg), m_ViwStp((ViwEnd - ViwBeg) / ViwN),
	m_DetSize(DetSize), m_DetN(DetN)
{
	m_DetStp.x = m_DetSize.x / m_DetN.x;
	m_DetStp.y = m_DetSize.y / m_DetN.y;
	m_DetCntIdx.x = m_DetN.x * 0.5f - 0.5f;
	m_DetCntIdx.y = m_DetN.y * 0.5f - 0.5f;
}
