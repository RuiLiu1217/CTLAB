#ifndef CONE_ED_GEO_H_
#define CONE_ED_GEO_H_
#include <cuda_runtime.h>
/// \brief Cone Beam geometry with equal distance detector
class ConeEDGeo
{
public:
	float m_S2O; ///< Source to object distance
	float m_O2D; ///< object to detector distance
	float m_S2D;///< source to detector distance

	int m_ViwN; ///< view number
	float m_ViwBeg; ///< begin view number
	float m_ViwStp; ///< view step size

	float2 m_DetSize; ///< detector size on length and height direction
	int2 m_DetN; ///< detector Number on length and height direction
	float2 m_DetStp; ///< detector step size on length and height direction
	float2 m_DetCntIdx; ///< detector center index
public:
	/// \brief constructor
	ConeEDGeo(void);
	/// \brief destructor
	~ConeEDGeo(){}
	/// \brief copy constructor
	ConeEDGeo(const ConeEDGeo& rhs);
	/// \brief Constructor
	/// \param S2O source to object distance
	/// \param O2D object to detector distance
	/// \param ViwN view number
	/// \param ViwBeg Begin of the view
	/// \param ViwEnd End of the view
	/// \param DetSize Detector size on both direction
	/// \param DetN detector cells on both direction
	ConeEDGeo(const float S2O, const float O2D, const  int ViwN,
		const float ViwBeg, const float ViwEnd, const float2 DetSize,
		const int2 DetN);
};

#endif