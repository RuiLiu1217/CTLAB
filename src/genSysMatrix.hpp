/*! 
 * \file genSysMatrix.hpp
 * \brief Generate the system matrix with different projection models
 *
 * \author Rui Liu
 * \date Jun. 1, 2015
 * \version 1.0
 */

/// \brief Generate the system Matrix namespace
namespace genSysMatrix{
	/// \brief system matrix generation with Siddon model for equal angular detector
	void genMatrix_EA_SIDDON();
	/// \brief system matrix generation with Siddon model for equal distance detector
	void genMatrix_ED_SIDDON();

	/// \brief system matrix generation with Area integral model for equal angular detector
	void genMatrix_EA_AIM();
	/// \brief system matrix generation with Area integral model for equal distance detector
	void genMatrix_ED_AIM();
	/// \brief system matrix generation with distance driven model for equal distance detector
	void genMatrix_ED_DDM();
}

