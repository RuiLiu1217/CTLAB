#ifndef CONE_EA_GEO_H_
#define CONE_EA_GEO_H_
#include <cuda_runtime.h>
/// \brief Cone Beam configuration with equal angle detector
class ConeEAGeo{
public:
	float m_S2O; ///< Source to object distance
	float m_O2D; ///< object to detector distance
	float m_S2D;///< source to detector distance

	unsigned int m_ViwN; ///< View number
	float m_ViwBeg;///< Begin of the view
	float m_ViwStp; ///< View step

	float m_DetArc; ///< Detector Arc size
	unsigned int m_DetN; ///< detector cell number
	float m_DetStp; ///< detector cell size

	float m_DetHei;  ///< detector height
	unsigned int m_DetHN; ///< Height cell number of the detector
	float m_DetHStp;///< detector step size on the height direction

	float2 m_DetCntIdx; ///< detector center index
public:
	/// \brief constructor
	ConeEAGeo(void);
	/// \brief destructor
	~ConeEAGeo(){}
	/// \brief copy constructor
	ConeEAGeo(const ConeEAGeo& rhs);
	/// \brief Constructor
	/// \param S2O source to object distance
	/// \param O2D object to detector distance
	/// \param viwN view number
	/// \param ViwBeg Begin of the view
	/// \param ViwEnd End of the view
	/// \param DetArc Arc size of the detector
	/// \param DetN detector cells on one row
	/// \param DetHei detector height size of the detector
	/// \param DetHN detector cells number on the height direction
	ConeEAGeo(const float S2O, const float O2D, const unsigned int viwN,
		const float ViwBeg, const float ViwEnd, const float DetArc, const unsigned int DetN,
		const float DetHei, const unsigned int DetHN);
};

#endif