/*
 * DD3_GPU_recon.h
 *
 *  Created on: Jun 8, 2016
 *      Author: Rui Liu
 *      E-mail: liurui1217@gmail.com
 */

#ifndef DD3_GPU_RECON_H_
#define DD3_GPU_RECON_H_

#include "utilities.cuh"
enum CGMethod{CG_FR, CG_PRP, CG_CW, CG_DI, CG_DY};

// Evaluate the projection/backprojection with all one input
void evalProjBack(
		bool forwardProjection, // true: forward projection; false: backprojection
		ForwardDDMethod forwardMethod, // projection model
		BackwardDDMethod backwardMethod); // backprojection model

// The OS-SART algorithm with one GPU (Titan X)
void OS_SART(thrust::host_vector<float>& reconImg, // The image volume to be reconstructed
		thrust::host_vector<float>& initImg, // The initial image volume
		thrust::host_vector<float>& hprj, // Projection data
		thrust::host_vector<float>& hangs, // projection views
		thrust::host_vector<float>& hzPos, // source positions in different views
		thrust::host_vector<byte>& mask, // image mask along in-plane direction
		const std::string& VolumeName, // The reconstructed volume name to be stored (Z, X, Y) order in float datatype
		const int osNum, // # of OS
		const int iterNum, // # of total iterations
		const bool outputMedRes, // output the intermediate results
		const float sid, // source to iso-center distance
		const float sdd, // source to detector distance
		const int DNU, // # of detector elements along in-plane direction
		const int DNV, // # of detector elements along bench moving direction
		const int PN, // # total view #
		const float imgXCenter, // volume center coordinate x
		const float imgYCenter, // volume center coordinate y
		const float imgZCenter, // volume center coordinate z
		const int XN, // # of volume pixel along x
		const int YN, // # of volume pixel along y
		const int ZN, // # of volume pixel along z
		const float dx, // image pixel size along in-plane direction
		const float dz, // slice thickness
		const float col_size, // detector elements size along in-plane direction
		const float row_size, // detector elements size along bench moving direction
		const float col_offset, // #(rational) of detector elements that offset from the center of the iso-center along in-plane direction
		const float row_offset, // #(rational) of detector elements that offset from the center of the iso-center along bench moving direction
		const ForwardDDMethod forwardMethod, // Forward projection model
		const BackwardDDMethod backwardMethod); // Back projection model

// The CG algorithm with one GPU (Titan X)
void CG(thrust::host_vector<float>& reconImg, // The image volume to be reconstructed
		thrust::host_vector<float>& initImg, // The initial image volume
		thrust::host_vector<float>& hprj, // Projection data
		thrust::host_vector<float>& hangs, // projection views
		thrust::host_vector<float>& hzPos, // source positions in different views
		thrust::host_vector<byte>& mask, // image mask along in-plane direction
		const std::string& VolumeName, // The reconstructed volume name to be stored (Z, X, Y) order in float datatype
		const int osNum, // # of OS
		const int iterNum, // # of total iterations
		const bool outputMedRes, // output the intermediate results
		const float sid, // source to iso-center distance
		const float sdd, // source to detector distance
		const int DNU, // # of detector elements along in-plane direction
		const int DNV, // # of detector elements along bench moving direction
		const int PN, // # total view #
		const float imgXCenter, // volume center coordinate x
		const float imgYCenter, // volume center coordinate y
		const float imgZCenter, // volume center coordinate z
		const int XN, // # of volume pixel along x
		const int YN, // # of volume pixel along y
		const int ZN, // # of volume pixel along z
		const float dx, // image pixel size along in-plane direction
		const float dz, // slice thickness
		const float col_size, // detector elements size along in-plane direction
		const float row_size, // detector elements size along bench moving direction
		const float col_offset, // #(rational) of detector elements that offset from the center of the iso-center along in-plane direction
		const float row_offset, // #(rational) of detector elements that offset from the center of the iso-center along bench moving direction
		const ForwardDDMethod forwardMethod, // Forward projection model
		const BackwardDDMethod backwardMethod, // Back projection model
		const CGMethod cgm); // CG beta updating method



#endif /* DD3_GPU_RECON_H_ */
