/*!
 * \file SF_proj.hpp
 * COPYRIGHT NOTICE
 * COPYRIGHT (c) 2015, Wake Forest and UMass Lowell
 * All rights reserved
 *
 * \brief The interface of CPU based SF projection/backprojection in conventional method
 * \version 1.0
 * \author Rui Liu
 * \date May. 1, 2015
 */
#include <vector>
#include <algorithm>
/// \brief Separable Footprint with Trapezoid-Rectangle A1 projection with equal distance detector in single floating datatype
/// \param proj projection dataset
/// \param vol image volume to be projected
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object size along X direction
/// \param objSizeY object size along Y direction
/// \param objSizeZ object size aLONG z direction
/// \param detSizeS detector size along channel direction
/// \param detSizeT detector size along bench moving direction
/// \param XN image pixel number along X direction
/// \param YN image pixel number along Y direction
/// \param ZN image pixel number along Z direction
/// \param DNS detector cell number along channel direction
/// \param DNT detector cell number bench moving direction
/// \param PN number of views
/// \param angs all views
/// \param angIdx current index of the angle
void SFTRA1_3D_ED_proj(std::vector<float>& proj, const std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx);

/// \brief Separable Footprint with Trapezoid-Rectangle A1 projection with equal distance detector in double floating datatype
/// \param proj projection dataset
/// \param vol image volume to be projected
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object size along X direction
/// \param objSizeY object size along Y direction
/// \param objSizeZ object size aLONG z direction
/// \param detSizeS detector size along channel direction
/// \param detSizeT detector size along bench moving direction
/// \param XN image pixel number along X direction
/// \param YN image pixel number along Y direction
/// \param ZN image pixel number along Z direction
/// \param DNS detector cell number along channel direction
/// \param DNT detector cell number bench moving direction
/// \param PN number of views
/// \param angs all views
/// \param angIdx current index of the angle
void SFTRA1_3D_ED_proj(std::vector<double>& proj, const std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx);

/// \brief Separable Footprint with Trapezoid-Rectangle A1 backprojection with equal distance detector in single floating datatype
/// \param proj projection dataset
/// \param vol image volume to be projected
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object size along X direction
/// \param objSizeY object size along Y direction
/// \param objSizeZ object size aLONG z direction
/// \param detSizeS detector size along channel direction
/// \param detSizeT detector size along bench moving direction
/// \param XN image pixel number along X direction
/// \param YN image pixel number along Y direction
/// \param ZN image pixel number along Z direction
/// \param DNS detector cell number along channel direction
/// \param DNT detector cell number bench moving direction
/// \param PN number of views
/// \param angs all views
/// \param angIdx current index of the angle
void SFTRA1_3D_ED_bproj(const std::vector<float>& proj, std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx);

/// \brief Separable Footprint with Trapezoid-Rectangle A1 backprojection with equal distance detector in double floating datatype
/// \param proj projection dataset
/// \param vol image volume to be projected
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object size along X direction
/// \param objSizeY object size along Y direction
/// \param objSizeZ object size aLONG z direction
/// \param detSizeS detector size along channel direction
/// \param detSizeT detector size along bench moving direction
/// \param XN image pixel number along X direction
/// \param YN image pixel number along Y direction
/// \param ZN image pixel number along Z direction
/// \param DNS detector cell number along channel direction
/// \param DNT detector cell number bench moving direction
/// \param PN number of views
/// \param angs all views
/// \param angIdx current index of the angle
void SFTRA1_3D_ED_bproj(const std::vector<double>& proj, std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx);

/// \brief Separable Footprint with Trapezoid-Rectangle A2 projection with equal distance detector in single floating datatype
/// \param proj projection dataset
/// \param vol image volume to be projected
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object size along X direction
/// \param objSizeY object size along Y direction
/// \param objSizeZ object size aLONG z direction
/// \param detSizeS detector size along channel direction
/// \param detSizeT detector size along bench moving direction
/// \param XN image pixel number along X direction
/// \param YN image pixel number along Y direction
/// \param ZN image pixel number along Z direction
/// \param DNS detector cell number along channel direction
/// \param DNT detector cell number bench moving direction
/// \param PN number of views
/// \param angs all views
/// \param angIdx current index of the angle
void SFTRA2_3D_ED_proj(std::vector<float>& proj, const std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx);

/// \brief Separable Footprint with Trapezoid-Rectangle A2 projection with equal distance detector in double floating datatype
/// \param proj projection dataset
/// \param vol image volume to be projected
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object size along X direction
/// \param objSizeY object size along Y direction
/// \param objSizeZ object size aLONG z direction
/// \param detSizeS detector size along channel direction
/// \param detSizeT detector size along bench moving direction
/// \param XN image pixel number along X direction
/// \param YN image pixel number along Y direction
/// \param ZN image pixel number along Z direction
/// \param DNS detector cell number along channel direction
/// \param DNT detector cell number bench moving direction
/// \param PN number of views
/// \param angs all views
/// \param angIdx current index of the angle
void SFTRA2_3D_ED_proj(std::vector<double>& proj, const std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx);
	
/// \brief Separable Footprint with Trapezoid-Rectangle A2 backprojection with equal distance detector in single floating datatype
/// \param proj projection dataset
/// \param vol image volume to be projected
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object size along X direction
/// \param objSizeY object size along Y direction
/// \param objSizeZ object size aLONG z direction
/// \param detSizeS detector size along channel direction
/// \param detSizeT detector size along bench moving direction
/// \param XN image pixel number along X direction
/// \param YN image pixel number along Y direction
/// \param ZN image pixel number along Z direction
/// \param DNS detector cell number along channel direction
/// \param DNT detector cell number bench moving direction
/// \param PN number of views
/// \param angs all views
/// \param angIdx current index of the angle	
	void SFTRA2_3D_ED_bproj(const std::vector<float>& proj, std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx);

/// \brief Separable Footprint with Trapezoid-Rectangle A2 backprojection with equal distance detector in double floating datatype
/// \param proj projection dataset
/// \param vol image volume to be projected
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object size along X direction
/// \param objSizeY object size along Y direction
/// \param objSizeZ object size aLONG z direction
/// \param detSizeS detector size along channel direction
/// \param detSizeT detector size along bench moving direction
/// \param XN image pixel number along X direction
/// \param YN image pixel number along Y direction
/// \param ZN image pixel number along Z direction
/// \param DNS detector cell number along channel direction
/// \param DNT detector cell number bench moving direction
/// \param PN number of views
/// \param angs all views
/// \param angIdx current index of the angle
void SFTRA2_3D_ED_bproj(const std::vector<double>& proj, std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx);
	
/// \brief Separable Footprint with Trapezoid-Trapezoid A1 projection with equal distance detector in single floating datatype
/// \param proj projection dataset
/// \param vol image volume to be projected
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object size along X direction
/// \param objSizeY object size along Y direction
/// \param objSizeZ object size aLONG z direction
/// \param detSizeS detector size along channel direction
/// \param detSizeT detector size along bench moving direction
/// \param XN image pixel number along X direction
/// \param YN image pixel number along Y direction
/// \param ZN image pixel number along Z direction
/// \param DNS detector cell number along channel direction
/// \param DNT detector cell number bench moving direction
/// \param PN number of views
/// \param angs all views
/// \param angIdx current index of the angle	
void SFTTA1_3D_ED_proj(std::vector<float>& proj, const std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx);

/// \brief Separable Footprint with Trapezoid-Trapezoid A1 projection with equal distance detector in double floating datatype
/// \param proj projection dataset
/// \param vol image volume to be projected
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object size along X direction
/// \param objSizeY object size along Y direction
/// \param objSizeZ object size aLONG z direction
/// \param detSizeS detector size along channel direction
/// \param detSizeT detector size along bench moving direction
/// \param XN image pixel number along X direction
/// \param YN image pixel number along Y direction
/// \param ZN image pixel number along Z direction
/// \param DNS detector cell number along channel direction
/// \param DNT detector cell number bench moving direction
/// \param PN number of views
/// \param angs all views
/// \param angIdx current index of the angle
void SFTTA1_3D_ED_proj(std::vector<double>& proj, const std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx);

/// \brief Separable Footprint with Trapezoid-Trapezoid A1 backprojection with equal distance detector in single floating datatype
/// \param proj projection dataset
/// \param vol image volume to be projected
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object size along X direction
/// \param objSizeY object size along Y direction
/// \param objSizeZ object size aLONG z direction
/// \param detSizeS detector size along channel direction
/// \param detSizeT detector size along bench moving direction
/// \param XN image pixel number along X direction
/// \param YN image pixel number along Y direction
/// \param ZN image pixel number along Z direction
/// \param DNS detector cell number along channel direction
/// \param DNT detector cell number bench moving direction
/// \param PN number of views
/// \param angs all views
/// \param angIdx current index of the angle	
void SFTTA1_3D_ED_bproj(const std::vector<float>& proj, std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx);

/// \brief Separable Footprint with Trapezoid-Trapezoid A1 backprojection with equal distance detector in double floating datatype
/// \param proj projection dataset
/// \param vol image volume to be projected
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object size along X direction
/// \param objSizeY object size along Y direction
/// \param objSizeZ object size aLONG z direction
/// \param detSizeS detector size along channel direction
/// \param detSizeT detector size along bench moving direction
/// \param XN image pixel number along X direction
/// \param YN image pixel number along Y direction
/// \param ZN image pixel number along Z direction
/// \param DNS detector cell number along channel direction
/// \param DNT detector cell number bench moving direction
/// \param PN number of views
/// \param angs all views
/// \param angIdx current index of the angle	
void SFTTA1_3D_ED_bproj(const std::vector<double>& proj, std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx);



/// \brief Separable Footprint with Trapezoid-Trapezoid A2 projection with equal distance detector in single floating datatype
/// \param proj projection dataset
/// \param vol image volume to be projected
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object size along X direction
/// \param objSizeY object size along Y direction
/// \param objSizeZ object size aLONG z direction
/// \param detSizeS detector size along channel direction
/// \param detSizeT detector size along bench moving direction
/// \param XN image pixel number along X direction
/// \param YN image pixel number along Y direction
/// \param ZN image pixel number along Z direction
/// \param DNS detector cell number along channel direction
/// \param DNT detector cell number bench moving direction
/// \param PN number of views
/// \param angs all views
/// \param angIdx current index of the angle	
void SFTTA2_3D_ED_proj(std::vector<float>& proj, const std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx);

/// \brief Separable Footprint with Trapezoid-Trapezoid A2 projection with equal distance detector in double floating datatype
/// \param proj projection dataset
/// \param vol image volume to be projected
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object size along X direction
/// \param objSizeY object size along Y direction
/// \param objSizeZ object size aLONG z direction
/// \param detSizeS detector size along channel direction
/// \param detSizeT detector size along bench moving direction
/// \param XN image pixel number along X direction
/// \param YN image pixel number along Y direction
/// \param ZN image pixel number along Z direction
/// \param DNS detector cell number along channel direction
/// \param DNT detector cell number bench moving direction
/// \param PN number of views
/// \param angs all views
/// \param angIdx current index of the angle	
void SFTTA2_3D_ED_proj(std::vector<double>& proj, const std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx);
	
/// \brief Separable Footprint with Trapezoid-Trapezoid A2 backprojection with equal distance detector in float floating datatype
/// \param proj projection dataset
/// \param vol image volume to be projected
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object size along X direction
/// \param objSizeY object size along Y direction
/// \param objSizeZ object size aLONG z direction
/// \param detSizeS detector size along channel direction
/// \param detSizeT detector size along bench moving direction
/// \param XN image pixel number along X direction
/// \param YN image pixel number along Y direction
/// \param ZN image pixel number along Z direction
/// \param DNS detector cell number along channel direction
/// \param DNT detector cell number bench moving direction
/// \param PN number of views
/// \param angs all views
/// \param angIdx current index of the angle	
void SFTTA2_3D_ED_bproj(const std::vector<float>& proj, std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx);

/// \brief Separable Footprint with Trapezoid-Trapezoid A2 backprojection with equal distance detector in double floating datatype
/// \param proj projection dataset
/// \param vol image volume to be projected
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object size along X direction
/// \param objSizeY object size along Y direction
/// \param objSizeZ object size aLONG z direction
/// \param detSizeS detector size along channel direction
/// \param detSizeT detector size along bench moving direction
/// \param XN image pixel number along X direction
/// \param YN image pixel number along Y direction
/// \param ZN image pixel number along Z direction
/// \param DNS detector cell number along channel direction
/// \param DNT detector cell number bench moving direction
/// \param PN number of views
/// \param angs all views
/// \param angIdx current index of the angle	
void SFTTA2_3D_ED_bproj(const std::vector<double>& proj, std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx);
