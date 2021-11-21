#ifndef FAN_EA_GEO_H_
#define FAN_EA_GEO_H_
/// \brief Fan Beam Equal Angle Detector based CT system
class FanEAGeo
{
public:
	float m_S2O; ///< source to object distance
	float m_O2D; ///< object to detector distance
	float m_S2D; ///< source to detector distance
	unsigned int m_ViwN; ///< view number
	float m_ViwBeg; ///< Begin viewer number
	float m_ViwStp; ///< View step size

	float m_DetArc; ///< Detector Arc angle
	unsigned int m_DetN; ///< Detector cells number in one row
	float m_DetStp; ///< Detector cells size
	float m_DetCntIdx; ///< The index of the center of the detector
public:
	FanEAGeo(void);
	~FanEAGeo(void){};
	FanEAGeo(const FanEAGeo& rhs);
	/// \brief constructor
	/// \param S2O source to object distance
	/// \param O2D object to detector distance
	/// \param ViwN view number
	/// \param ViwBeg the begin view
	/// \param ViwEnd the End view
	/// \param DetArc the detector Arc
	/// \param DetN the number of detector cells on one row
	FanEAGeo(const float S2O, const float O2D, const unsigned int ViwN,
		const float ViwBeg, const float ViwEnd, const float DetArc,
		const unsigned int DetN);
};
#endif