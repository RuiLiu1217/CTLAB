/*
 *
 * \file DDM_proj.hpp
 * \brief The interface of the CPU/GPU based DD projection in conventional method
 * COPYRIGHT NOTICE
 * COPYRIGHT (c) 2015, Wake Forest and UMass Lowell
 * All rights reserved
 *
 * \version 1.0
 * \author Rui Liu
 * \date May. 1, 2015
 * \email liurui1217@gmail.com
 *
 */


#ifndef __DDM_PROJ_H__
#define __DDM_PROJ_H__

#include "cuda_runtime.h"
#include <vector>
#include <thrust/device_vector.h>

/////////////////////////////////////////////////////////////////////////////////////////
// Equal distance detector projection
/////////////////////////////////////////////////////////////////////////////////////////

/// \brief Equal Distance DD Projection in Fan Beam Geometry
/// \param proj projection data
/// \param img image volume
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object Size in X direction
/// \param objSizeY object Size in Y direction
/// \param detSize detector size
/// \param detCntIdx detector center index
/// \param XN number of image pixels in X direction
/// \param YN number of image pixels in Y direction
/// \param DN number of detector cells
/// \param PN number of views (may be eliminated)
/// \param angs projection views
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

/// \brief Equal Distance DD Projection in Cone Beam Geometry
/// \param proj projection data
/// \param vol image volume
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object Size in X direction
/// \param objSizeY object Size in Y direction
/// \param objSizeZ object Size in Z direction
/// \param detSizeU detector size along xy plane
/// \param detSizeV detector size along z direction
/// \param detCntIdU detector center index along xy plane
/// \param detCntIdV detector center index along z plane
/// \param XN number of image pixels in X direction
/// \param YN number of image pixels in Y direction
/// \param ZN number of image pixels in Z direction
/// \param DNU number of detector cells along xy direction
/// \param DNV number of detector cells along z direction
/// \param PN number of views (may be eliminated)
/// \param angs projection views
void DDM3D_ED_proj(std::vector<float>& proj, const std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeU, const float detSizeV,
	const float detCntIdxU, const float detCntIdxV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<float>& angs);
void DDM3D_ED_proj(std::vector<double>& proj, const std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeU, const double detSizeV,
	const double detCntIdxU, const double detCntIdxV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<double>& angs);

/// \brief Equal distance detector DD projection in GPU in Fan Beam Geometry (all data are all in GPU). NOTE: pointers represent in device memory
/// \param proj projection data
/// \param img image volume
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object Size in X direction
/// \param objSizeY object Size in Y direction
/// \param detSize detector size along xy plane
/// \param detCntIdx detector center index along xy plane
/// \param XN number of image pixels in X direction
/// \param YN number of image pixels in Y direction
/// \param DN number of detector cells along xy direction
/// \param PN number of views (may be eliminated)
/// \param dd detector cell size
/// \param dx image pixel size along x
/// \param dy image pixel size along y
/// \param hfXN image center index
/// \param hfYN image center index
/// \param S2D source to detector distance
/// \param angs projection views
/// \param blk threads block size
/// \param gid block grid size
void DDM_ED_proj_GPU(float* proj, const float* img,
	const float S2O, const float O2D, const float objSizeX, const float objSizeY,
	const float detSize, const float detCntIdx,
	const int XN, const int YN, const int DN, const int PN, const float dd, const float dx, const float dy,
	const float* angs, const dim3 blk, const dim3 gid);
void DDM_ED_proj_GPU(double* proj, const double* img,
	const double S2O, const double O2D, const double objSizeX, const double objSizeY,
	const double detSize, const double detCntIdx,
	const int XN, const int YN, const int DN, const int PN, const double dd, const double dx, const double dy,
	const double* angs, const dim3 blk, const dim3 gid);
void DDM_ED_proj_GPU(thrust::device_vector<float>& proj, const thrust::device_vector<float>& img,
	const float S2O, const float O2D, const float objSizeX, const float objSizeY,
	const float detSize, const float detCntIdx,
	const int XN, const int YN, const int DN, const int PN, const float dd, const float dx, const float dy,
	const thrust::device_vector<float>& angs, const dim3 blk, const dim3 gid);
void DDM_ED_proj_GPU(thrust::device_vector<double>& proj, const thrust::device_vector<double>& img,
	const double S2O, const double O2D, const double objSizeX, const double objSizeY,
	const double detSize, const double detCntIdx,
	const int XN, const int YN, const int DN, const int PN, const double dd, const double dx, const double dy,
	const thrust::device_vector<double>& angs, const dim3 blk, const dim3 gid);
	

/// \brief Equal distance DD Projection in Cone Beam Geometry in GPU (all data are all in GPU). NOTE: pointers represent in device memory
/// \param proj projection data
/// \param vol image volume
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object Size in X direction
/// \param objSizeY object Size in Y direction
/// \param objSizeZ object Size in Z direction
/// \param detSizeU detector size along xy plane
/// \param detSizeV detector size along z plane
/// \param detCntIdU detector center index along xy plane
/// \param detCntIdV detector center index along z plane
/// \param XN number of image pixels in X direction
/// \param YN number of image pixels in Y direction
/// \param ZN number of image pixels in Z direction
/// \param DNU number of detector cells along xy direction
/// \param DNV number of detector cells along z direction
/// \param PN number of views (may be eliminated)
/// \param ddu detector cell size along xy direction
/// \param ddv detector cell size along z direction
/// \param dx image pixel size along x
/// \param dy image pixel size along y
/// \param dz image pixel size along z
/// \param angs projection views
/// \param blk threads block size
/// \param gid block grid size
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
void DDM3D_ED_proj_GPU(thrust::device_vector<float>& proj, const thrust::device_vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeU, const float detSizeV,
	const float detCntIdxU, const float detCntIdxV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float ddu, const float ddv, const float dx, const float dy, const float dz,
	const thrust::device_vector<float>& angs, const dim3 blk, const dim3 gid);
void DDM3D_ED_proj_GPU(thrust::device_vector<double>& proj, const thrust::device_vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeU, const double detSizeV,
	const double detCntIdxU, const double detCntIdxV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double ddu, const double ddv, const double dx, const double dy, const double dz,
	const thrust::device_vector<double>& angs, const dim3 blk, const dim3 gid);

/// \brief Equal Angular DD Backprojection in Fan Beam Geometry
/// \param proj projection data
/// \param img image volume
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object Size in X direction
/// \param objSizeY object Size in Y direction
/// \param detArc detector arc
/// \param detCntIdx detector center index
/// \param XN number of image pixels in X direction
/// \param YN number of image pixels in Y direction
/// \param DN number of detector cells
/// \param PN number of views (may be eliminated)
/// \param angs projection views
void DDM_EA_proj(std::vector<float>& proj, const std::vector<float>& img,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY,
	const float detArc,
	const float detCntIdx,
	const int XN, const int YN, const int DN, const int PN,
	const std::vector<float>& angs);
void DDM_EA_proj(std::vector<double>& proj, const std::vector<double>& img,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY,
	const double detArc,
	const double detCntIdx,
	const int XN, const int YN, const int DN, const int PN,
	const std::vector<double>& angs);

/// \brief Equal Angular DD Projection in Cone Beam Geometry
/// \param proj projection data
/// \param vol image volume
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object Size in X direction
/// \param objSizeY object Size in Y direction
/// \param objSizeZ object Size in Z direction
/// \param detArc detector arc along xy plane
/// \param detSizeV detector size along z direction
/// \param detCntIdU detector center index along xy plane
/// \param detCntIdV detector center index along z plane
/// \param XN number of image pixels in X direction
/// \param YN number of image pixels in Y direction
/// \param ZN number of image pixels in Z direction
/// \param DNU number of detector cells along xy direction
/// \param DNV number of detector cells along z direction
/// \param PN number of views (may be eliminated)
/// \param angs projection views
void DDM3D_EA_proj(std::vector<float>& proj, const std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<float>& angs);
void DDM3D_EA_proj(std::vector<double>& proj, const std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<double>& angs);

/// \brief Equal Angular DD Projection in Helical Beam Geometry
/// \param proj projection data
/// \param vol image volume
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object Size in X direction
/// \param objSizeY object Size in Y direction
/// \param objSizeZ object Size in Z direction
/// \param detArc detector arc along xy plane
/// \param detSizeV detector size along z direction
/// \param detCntIdU detector center index along xy plane
/// \param detCntIdV detector center index along z plane
/// \param XN number of image pixels in X direction
/// \param YN number of image pixels in Y direction
/// \param ZN number of image pixels in Z direction
/// \param DNU number of detector cells along xy direction
/// \param DNV number of detector cells along z direction
/// \param PN number of views (may be eliminated)
/// \param initZPos initial source z position
/// \param pitch pitch size
/// \param angs projection views
void DDM3D_EA_helical_proj(std::vector<float>& proj, const std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float initZPos, const float pitch, const std::vector<float>& angs);
void DDM3D_EA_helical_proj(std::vector<double>& proj, const std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double initZPos, const double pitch, const std::vector<double>& angs);

/// \brief Equal angular DD Backprojection in Helical Beam Geometry in GPU (all data are all in GPU). NOTE: pointers represent in device memory
/// \param proj projection data
/// \param vol image volume
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object Size in X direction
/// \param objSizeY object Size in Y direction
/// \param objSizeZ object Size in Z direction
/// \param detArc detector arc along xy plane
/// \param detSizeH detector size along z plane
/// \param detCntIdU detector center index along xy plane
/// \param detCntIdV detector center index along z plane
/// \param XN number of image pixels in X direction
/// \param YN number of image pixels in Y direction
/// \param ZN number of image pixels in Z direction
/// \param DNU number of detector cells along xy direction
/// \param DNV number of detector cells along z direction
/// \param PN number of views (may be eliminated)
/// \param dbeta detector cell size along xy direction
/// \param ddv detector cell size along z direction
/// \param dx image pixel size along x
/// \param dy image pixel size along y
/// \param dz image pixel size along z
/// \param initZPos initial Z position of the source/detector
/// \param pitch pitch size
/// \param angs projection views
/// \param blk threads block size
/// \param gid block grid size
void DDM3D_EA_helical_proj_GPU(float* proj, const float*  vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float dbeta, const float ddv, const float dx, const float dy, const float dz,
	const float initZPos, const float pitch, const float* angs, const dim3 blk, const dim3 gid);
void DDM3D_EA_helical_proj_GPU(double* proj, const double*  vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double dbeta, const double ddv, const double dx, const double dy, const double dz,
	const double initZPos, const double pitch, const double* angs, const dim3 blk, const dim3 gid);

/// \brief Equal angular DD Backprojection in Helical Beam Geometry in GPU (all data are all in GPU). NOTE: pointers represent in device memory
/// \param proj projection data
/// \param vol image volume
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object Size in X direction
/// \param objSizeY object Size in Y direction
/// \param objSizeZ object Size in Z direction
/// \param detArc detector arc along xy plane
/// \param detSizeH detector size along z plane
/// \param detCntIdU detector center index along xy plane
/// \param detCntIdV detector center index along z plane
/// \param XN number of image pixels in X direction
/// \param YN number of image pixels in Y direction
/// \param ZN number of image pixels in Z direction
/// \param DNU number of detector cells along xy direction
/// \param DNV number of detector cells along z direction
/// \param PN number of views (may be eliminated)
/// \param dbeta detector cell size along xy direction
/// \param ddv detector cell size along z direction
/// \param dx image pixel size along x
/// \param dy image pixel size along y
/// \param dz image pixel size along z
/// \param zShifts Z positions of the source/detector
/// \param angs projection views
/// \param blk threads block size
/// \param gid block grid size
void DDM3D_EA_helical_proj_GPU(float* proj, const float*  vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float dbeta, const float ddv, const float dx, const float dy, const float dz,
	const float* zShifts, const float* angs, const dim3 blk, const dim3 gid);
void DDM3D_EA_helical_proj_GPU(double* proj, const double*  vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double dbeta, const double ddv, const double dx, const double dy, const double dz,
	const double* zShifts, const double* angs, const dim3 blk, const dim3 gid);

/// \brief Equal angular DD Backprojection in Helical Beam Geometry in GPU (all data are all in GPU). NOTE: pointers represent in device memory
/// \param proj projection data
/// \param vol image volume
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object Size in X direction
/// \param objSizeY object Size in Y direction
/// \param objSizeZ object Size in Z direction
/// \param detArc detector arc along xy plane
/// \param detSizeH detector size along z plane
/// \param detCntIdU detector center index along xy plane
/// \param detCntIdV detector center index along z plane
/// \param XN number of image pixels in X direction
/// \param YN number of image pixels in Y direction
/// \param ZN number of image pixels in Z direction
/// \param DNU number of detector cells along xy direction
/// \param DNV number of detector cells along z direction
/// \param PN number of views (may be eliminated)
/// \param dbeta detector cell size along xy direc/tion
/// \param ddv detector cell size along z direction
/// \param dx image pixel size along x
/// \param dy image pixel size along y
/// \param dz image pixel size along z
/// \param initZPos initial Z position of the source/detector
/// \param pitch pitch size
/// \param angs projection views
/// \param blk threads block size
/// \param gid block grid size
void DDM3D_EA_helical_proj_GPU(thrust::device_vector<float>& proj, const thrust::device_vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float dbeta, const float ddv, const float dx, const float dy, const float dz,
	const float initZPos, const float pitch, const thrust::device_vector<float>& angs, const dim3 blk, const dim3 gid);
void DDM3D_EA_helical_proj_GPU(thrust::device_vector<double>& proj, const thrust::device_vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double dbeta, const double ddv, const double dx, const double dy, const double dz,
	const double initZPos, const double pitch, const thrust::device_vector<double>& angs, const dim3 blk, const dim3 gid);

/// \brief Equal angular DD Backprojection in Helical Beam Geometry in GPU (all data are all in GPU). NOTE: pointers represent in device memory
/// \param proj projection data
/// \param vol image volume
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object Size in X direction
/// \param objSizeY object Size in Y direction
/// \param objSizeZ object Size in Z direction
/// \param detArc detector arc along xy plane
/// \param detSizeH detector size along z plane
/// \param detCntIdU detector center index along xy plane
/// \param detCntIdV detector center index along z plane
/// \param XN number of image pixels in X direction
/// \param YN number of image pixels in Y direction
/// \param ZN number of image pixels in Z direction
/// \param DNU number of detector cells along xy direction
/// \param DNV number of detector cells along z direction
/// \param PN number of views (may be eliminated)
/// \param dbeta detector cell size along xy direction
/// \param ddv detector cell size along z direction
/// \param dx image pixel size along x
/// \param dy image pixel size along y
/// \param dz image pixel size along z
/// \param zShifts Z positions of the source/detector
/// \param angs projection views
/// \param blk threads block size
/// \param gid block grid size
void DDM3D_EA_helical_proj_GPU(thrust::device_vector<float>& proj, const thrust::device_vector<float>&  vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float dbeta, const float ddv, const float dx, const float dy, const float dz,
	const thrust::device_vector<float>& zShifts, const thrust::device_vector<float>& angs, const dim3 blk, const dim3 gid);
void DDM3D_EA_helical_proj_GPU(thrust::device_vector<double>& proj, const thrust::device_vector<double>&  vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double dbeta, const double ddv, const double dx, const double dy, const double dz,
	const thrust::device_vector<double>& zShifts, const thrust::device_vector<double>& angs, const dim3 blk, const dim3 gid);


#endif
