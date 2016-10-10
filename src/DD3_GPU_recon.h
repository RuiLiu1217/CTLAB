/*!
 * \file DD3_GPU_recon.h
 * Wake Forest Health Sciences & University of Massachusetts Lowell
 * E-mail: liurui1217@gmail.com
 * \date Jun 8, 2016
 * \author Rui Liu
 * \version 1.0
 */
#ifndef DD3_GPU_RECON_H_
#define DD3_GPU_RECON_H_

#include "utilities.cuh"
/// \brief The updating methods for Conjugate Gradient method
enum CGMethod{CG_FR, ///<
	      CG_PRP,///<
	      CG_CW, ///<
	      CG_DI, ///<
	      CG_DY  ///<
	     };

/// \brief Evaluate the projection/backprojection with all one input
/// \param forwardProjection Forward projection (true) or backprojection (false)
/// \param forwardMethod Forward Projection mode selection
/// \param backwardMethod Backprojection mode selection
void evalProjBack(
		bool forwardProjection, // true: forward projection; false: backprojection
		ForwardDDMethod forwardMethod, // projection model
		BackwardDDMethod backwardMethod); // backprojection model

/// \brief The OS-SART algorithm with one GPU (Titan X)
/// \param reconImg The image volume to be reconstructed
/// \param initImg The initial image volume
/// \param hrpj Projection data
/// \param hangs Projection views
/// \param hzPos Source positions in different views
/// \param mask Image mask along in-plane direction
/// \param volumeName The reconstructed volume name to be stored (Z, X, Y) order in float datatype
/// \param osNum Number of OS
/// \param iterNum Number of total iterations
/// \param outputMedRes Output the intermediate results
/// \param sid Source to iso-center distance
/// \param sdd Source to detector distance
/// \param DNU Number of detector elements along in-plane direction
/// \param DNV Number of detector elements along bench moving direction
/// \param PN Number of total views
/// \param imgXCenter Volume center coordinate X
/// \param imgYCenter Volume center coordinate Y
/// \param imgZCenter Volume center coordinate Z
/// \param XN Number of volume pixel along X
/// \param YN Number of volume pixel along Y
/// \param ZN Number of volume pixel along Z
/// \param dx Image pixel size along in-plane direction
/// \param dz Slice thickness
/// \param col_size Detector elements size along in-plane direction
/// \param row_size Detector elements size along bench moving direction
/// \param col_offset Number (rational) of detector elements that offset from the center of the iso-center along in-plane direction
/// \param row_offset Number (rational) of detector elements that offset from the center of the iso-center along bench moving direction
/// \param forwardMethod Forward projection model
/// \param backwardMethod Back projection model
void OS_SART(thrust::host_vector<float>& reconImg, thrust::host_vector<float>& initImg,  thrust::host_vector<float>& hprj, 
thrust::host_vector<float>& hangs, thrust::host_vector<float>& hzPos, thrust::host_vector<byte>& mask, 
const std::string& VolumeName, const int osNum,  const int iterNum, const bool outputMedRes, const float sid, const float sdd,
const int DNU, const int DNV, const int PN, const float imgXCenter, const float imgYCenter, const float imgZCenter, 
const int XN, const int YN, const int ZN, const float dx, const float dz, const float col_size, const float row_size, 
const float col_offset, const float row_offset,	const ForwardDDMethod forwardMethod, const BackwardDDMethod backwardMethod); // Back projection model

/// \brief The CG algorithm with one GPU (Titan X)
/// \param reconImg The image volume to be reconstructed
/// \param initImg The initial image volume
/// \param hrpj Projection data
/// \param hangs Projection views
/// \param hzPos Source positions in different views
/// \param mask Image mask along in-plane direction
/// \param volumeName The reconstructed volume name to be stored (Z, X, Y) order in float datatype
/// \param osNum Number of OS
/// \param iterNum Number of total iterations
/// \param outputMedRes Output the intermediate results
/// \param sid Source to iso-center distance
/// \param sdd Source to detector distance
/// \param DNU Number of detector elements along in-plane direction
/// \param DNV Number of detector elements along bench moving direction
/// \param PN Number of total views
/// \param imgXCenter Volume center coordinate X
/// \param imgYCenter Volume center coordinate Y
/// \param imgZCenter Volume center coordinate Z
/// \param XN Number of volume pixel along X
/// \param YN Number of volume pixel along Y
/// \param ZN Number of volume pixel along Z
/// \param dx Image pixel size along in-plane direction
/// \param dz Slice thickness
/// \param col_size Detector elements size along in-plane direction
/// \param row_size Detector elements size along bench moving direction
/// \param col_offset Number (rational) of detector elements that offset from the center of the iso-center along in-plane direction
/// \param row_offset Number (rational) of detector elements that offset from the center of the iso-center along bench moving direction
/// \param forwardMethod Forward projection model
/// \param backwardMethod Back projection model
/// \param cgm CG beta updating method
void CG(thrust::host_vector<float>& reconImg, thrust::host_vector<float>& initImg, thrust::host_vector<float>& hprj, 
thrust::host_vector<float>& hangs, thrust::host_vector<float>& hzPos, thrust::host_vector<byte>& mask,
const std::string& VolumeName, const int osNum, const int iterNum, const bool outputMedRes, const float sid, const float sdd, 
const int DNU, const int DNV, const int PN, const float imgXCenter, const float imgYCenter, const float imgZCenter, 
const int XN, const int YN, const int ZN, const float dx, const float dz, const float col_size, const float row_size, 
const float col_offset,	const float row_offset, const ForwardDDMethod forwardMethod, const BackwardDDMethod backwardMethod, const CGMethod cgm); // CG beta updating method
#endif /* DD3_GPU_RECON_H_ */
