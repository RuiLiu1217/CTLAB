/*
 * COPYRIGHT NOTICE
 * COPYRIGHT (c) 2015, Wake Forest and UMass Lowell
 * All rights reserved
 *
 * @file DDM_proj.hpp
 * @brief The interface of the CPU/GPU based DD projection in conventional method
 *
 * @version 1.0
 * @author Rui Liu
 * @date May. 1, 2015
 * @email liurui1217@gmail.com
 *
 */


#ifndef __DDM_PROJ_H__
#define __DDM_PROJ_H__

#include "cuda_runtime.h"
#include <vector>

void DDM_ED_proj(std::vector<float>& proj, const std::vector<float>& img,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY,
	const float detSize,
	const float detCntIdx,
	const int XN, const int YN, const int DN, const int PN,
	const std::vector<float>& angs);
void DDM_ED_proj(std::vector<double>& proj, const std::vector<double>& img,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY,
	const double detSize,
	const double detCntIdx,
	const int XN, const int YN, const int DN, const int PN,
	const std::vector<double>& angs);

void DDM3D_ED_proj(std::vector<double>& proj, const std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeU, const double detSizeV,
	const double detCntIdxU, const double detCntIdxV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<double>& angs);

void DDM3D_ED_proj(std::vector<float>& proj, const std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeU, const float detSizeV,
	const float detCntIdxU, const float detCntIdxV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<float>& angs);


void DDM_ED_bproj(const std::vector<double>& proj, std::vector<double>& img,
	const double S2O, const double O2D, const double objSizeX, const double objSizeY,
	const double detSize, const double detCntIdx, const int XN, const int YN, const int DN, const int PN,
	const std::vector<double>& angs);
void DDM_ED_bproj(const std::vector<float>& proj, std::vector<float>& img,
	const float S2O, const float O2D, const float objSizeX, const float objSizeY,
	const float detSize, const float detCntIdx, const int XN, const int YN, const int DN, const int PN,
	const std::vector<float>& angs);



void DDM3D_bproj(const std::vector<double>& proj, std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeU, const double detSizeV,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<double>& angs);

void DDM3D_bproj(const std::vector<float>& proj, std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeU, const float detSizeV,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<float>& angs);


void DDM_ED_proj_GPU(double* proj, const double* img,
	const double S2O, const double O2D, const double objSizeX, const double objSizeY,
	const double detSize, const double detCntIdx,
	const int XN, const int YN, const int DN, const int PN, const double dd, const double dx, const double dy,
	const double* angs, const dim3 blk, const dim3 gid);

void DDM_ED_proj_GPU(float* proj, const float* img,
	const float S2O, const float O2D, const float objSizeX, const float objSizeY,
	const float detSize, const float detCntIdx,
	const int XN, const int YN, const int DN, const int PN, const float dd, const float dx, const float dy,
	const float* angs, const dim3 blk, const dim3 gid);


void DDM3D_ED_proj_GPU(double* proj, const double* vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeU, const double detSizeV,
	const double detCntIdxU, const double detCntIdxV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double ddu, const double ddv, const double dx, const double dy, const double dz,
	const double* angs, const dim3 blk, const dim3 gid);

void DDM3D_ED_proj_GPU(float* proj, const float* vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeU, const float detSizeV,
	const float detCntIdxU, const float detCntIdxV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float ddu, const float ddv, const float dx, const float dy, const float dz,
	const float* angs, const dim3 blk, const dim3 gid);



void DDM_ED_bproj_GPU(const double* proj, double* img,
	const double S2O, const double O2D, const double objSizeX, const double objSizeY,
	const double detSize, const double detCntIdx, const int XN, const int YN, const int DN, const int PN,
	const double dd, const double dx, const double dy, const double hfXN, const double hfYN,
	const double S2D,
	const double* angs, const dim3 blk, const dim3 gid);

void DDM_ED_bproj_GPU(const float* proj, float* img,
	const float S2O, const float O2D, const float objSizeX, const float objSizeY,
	const float detSize, const float detCntIdx, const int XN, const int YN, const int DN, const int PN,
	const float dd, const float dx, const float dy, const float hfXN, const float hfYN,
	const float S2D,
	const float* angs, const dim3 blk, const dim3 gid);

void DDM3D_ED_bproj_GPU(const double* proj, double* vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeU, const double detSizeV,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double ddu, const double ddv, const double dx, const double dy, const double dz,
	const double hfXN, const double hfYN, const double hfZN, const double S2D,
	const double* angs, const dim3 blk, const dim3 gid);

void DDM3D_ED_bproj_GPU(const float* proj, float* vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeU, const float detSizeV,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float ddu, const float ddv, const float dx, const float dy, const float dz,
	const float hfXN, const float hfYN, const float hfZN, const float S2D,
	const float* angs, const dim3 blk, const dim3 gid);


void genMatrix_DDM_ED(
	std::vector<int>& rowIdx,
	std::vector<int>& colIdx,
	std::vector<double>& weight,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY,
	const double detSize,
	const double detCntIdx,
	const int XN, const int YN, const int DN, const int PN,
	const std::vector<double>& angs);
	
void genMatrix_DDM_ED(
	std::vector<int>& rowIdx,
	std::vector<int>& colIdx,
	std::vector<float>& weight,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY,
	const float detSize,
	const float detCntIdx,
	const int XN, const int YN, const int DN, const int PN,
	const std::vector<float>& angs);

#endif
