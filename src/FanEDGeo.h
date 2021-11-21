#ifndef FAN_ED_GEO_H_
#define FAN_ED_GEO_H_
/// \brief Fan Beam configuration with equal distance detector
class FanEDGeo
{
public:
	float m_S2O; ///< Source to Object distance
	float m_O2D; ///< Object to Detector distance
	float m_S2D; ///< Source to Detector distance

	unsigned int m_ViwN; ///< view number
	float m_ViwBeg;///< View begin
	float m_ViwStp; ///< View Step

	float m_DetSize; ///< Detector size;
	unsigned int m_DetN; ///< How many detector cells
	float m_DetStp; ///< detector cells size;
	float m_DetCntIdx; ///< detector center index;
public:
	/// \brief constructor
	FanEDGeo(void);
	/// \brief destructor
	~FanEDGeo(){};
	/// \brief copy constructor
	FanEDGeo(const FanEDGeo& rhs);
	/// \brief Constructor
	/// \param S2O source to object distance
	/// \param O2D object to detector distance
	/// \param ViwN View number
	/// \param ViwBeg The begin of the view
	/// \param ViwEnd The end of the view
	/// \param DetSize Detector size
	/// \param DetN detector cells number on one row
	FanEDGeo(const float S2O, const float O2D, const unsigned int ViwN,
		const float ViwBeg, const float ViwEnd, const float DetSize,
		const unsigned int DetN);

};
#endif