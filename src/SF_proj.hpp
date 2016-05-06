/*
 * COPYRIGHT NOTICE
 * COPYRIGHT (c) 2015, Wake Forest and UMass Lowell
 * All rights reserved
 *
 * @file SF_proj.hpp
 * @brief The interface of CPU based SF projection/backprojection in conventional method
 *
 * @version 1.0
 * @author Rui Liu
 * @date May. 1, 2015
 *
 */


#include <vector>
#include <algorithm>
void SFTRA1_3D_ED_proj(std::vector<float>& proj, const std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx);
void SFTRA1_3D_ED_proj(std::vector<double>& proj, const std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx);


void SFTRA1_3D_ED_bproj(const std::vector<float>& proj, std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx);
void SFTRA1_3D_ED_bproj(const std::vector<double>& proj, std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx);




void SFTRA2_3D_ED_proj(std::vector<float>& proj, const std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx);
void SFTRA2_3D_ED_proj(std::vector<double>& proj, const std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx);
	
	
	
void SFTRA2_3D_ED_bproj(const std::vector<float>& proj, std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx);
void SFTRA2_3D_ED_bproj(const std::vector<double>& proj, std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx);
	
	
void SFTTA1_3D_ED_proj(std::vector<float>& proj, const std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx);
void SFTTA1_3D_ED_proj(std::vector<double>& proj, const std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx);

	
void SFTTA1_3D_ED_bproj(const std::vector<float>& proj, std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx);
void SFTTA1_3D_ED_bproj(const std::vector<double>& proj, std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx);




void SFTTA2_3D_ED_proj(std::vector<float>& proj, const std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx);
void SFTTA2_3D_ED_proj(std::vector<double>& proj, const std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx);
	
	
	
void SFTTA2_3D_ED_bproj(const std::vector<float>& proj, std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx);
void SFTTA2_3D_ED_bproj(const std::vector<double>& proj, std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx);
