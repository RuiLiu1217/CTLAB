///*! 
// * \file genSysMatrix.hpp
// * \brief Generate the system matrix with different projection models
// *
// * \author Rui Liu
// * \date Jun. 1, 2015
// * \version 1.0
// */
//
///// \brief Generate the system Matrix namespace
//namespace genSysMatrix{
//	/// \brief system matrix generation with Siddon model for equal angular detector
//	void genMatrix_EA_SIDDON();
//	/// \brief system matrix generation with Siddon model for equal distance detector
//	void genMatrix_ED_SIDDON();
//
//	/// \brief system matrix generation with Area integral model for equal angular detector
//	void genMatrix_EA_AIM();
//	/// \brief system matrix generation with Area integral model for equal distance detector
//	void genMatrix_ED_AIM();
//	/// \brief system matrix generation with distance driven model for equal distance detector
//	void genMatrix_ED_DDM();
//}
//
//
///// \brief Generate the 2D system matrix with Equal distance distance in double floating datatype
///// The system matrix only stores the non-zero elements in (rowIdx, colIdx, weight)
///// \param rowIdx row index of the system matrix that is not 0
///// \param colIdx column index of the system matrix that is not 0
///// \param weight non-zero values
///// \param S2O source to iso-center distance
///// \param O2D object to detector distance
///// \param objSizeX size of the object along X direction
///// \param objSizeY size of the object along Y direction
///// \param detSize detector size
///// \param detCntIdx detector center index
///// \param XN image pixel number along X direction
///// \param YN image pixel number along Y direction
///// \param DN detector cell number along channel direction
///// \param PN number of views
///// \param angs all view angles
//void genMatrix_DDM_ED(
//	std::vector<int>& rowIdx,
//	std::vector<int>& colIdx,
//	std::vector<double>& weight,
//	const double S2O, const double O2D,
//	const double objSizeX, const double objSizeY,
//	const double detSize,
//	const double detCntIdx,
//	const int XN, const int YN, const int DN, const int PN,
//	const std::vector<double>& angs);
//
///// \brief Generate the 2D system matrix with Equal distance distance in single floating datatype
///// The system matrix only stores the non-zero elements in (rowIdx, colIdx, weight)
///// \param rowIdx row index of the system matrix that is not 0
///// \param colIdx column index of the system matrix that is not 0
///// \param weight non-zero values
///// \param S2O source to iso-center distance
///// \param O2D object to detector distance
///// \param objSizeX size of the object along X direction
///// \param objSizeY size of the object along Y direction
///// \param detSize detector size
///// \param detCntIdx detector center index
///// \param XN image pixel number along X direction
///// \param YN image pixel number along Y direction
///// \param DN detector cell number along channel direction
///// \param PN number of views
///// \param angs all view angles	
//void genMatrix_DDM_ED(
//	std::vector<int>& rowIdx,
//	std::vector<int>& colIdx,
//	std::vector<float>& weight,
//	const float S2O, const float O2D,
//	const float objSizeX, const float objSizeY,
//	const float detSize,
//	const float detCntIdx,
//	const int XN, const int YN, const int DN, const int PN,
//	const std::vector<float>& angs);