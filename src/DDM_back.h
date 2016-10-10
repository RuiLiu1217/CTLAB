/*
* \file DDM_back.h
* \brief The head file of brute force backprojection
*
* \date Mar 4, 2016
* \author Rui Liu
* \email liurui1217@gmail.com
* \version 1.0
*/
#ifndef _DDM_BACK_H_
#define _DDM_BACK_H_
#include <vector>
#include <thrust/device_vector.h>
/// \brief Equal Distance DD Backprojection in double floating point precision 2D geometry
/// \param proj pinter to projection data
/// \param img pointer to image volume
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
/// \param dd detector cell size
/// \param dx size of the pixel along x direction
/// \param dy size of the pixel along y direction
/// \param hfXN center index of the image along X cordinate
/// \param hfYN center index of the image along Y cordinate
/// \param angs projection views
/// \param blk Block size configuration for projection in GPU
/// \param gid Grid size configuration for projection in GPU
void DDM_ED_bproj_GPU(const double* proj, double* img, const double S2O, const double O2D, 
const double objSizeX, const double objSizeY, const double detSize, const double detCntIdx, 
const int XN, const int YN, const int DN, const int PN,	const double dd, const double dx, const double dy, 
const double hfXN, const double hfYN, const double S2D, const double* angs, const dim3 blk, const dim3 gid);

/// \brief Equal Distance DD Backprojection in single floating point precision 2D Geometry
/// \param proj pointer to projection data
/// \param img pointer to image volume
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
/// \param dd detector cell size
/// \param dx size of the pixel along x direction
/// \param dy size of the pixel along y direction
/// \param hfXN center index of the image along X cordinate
/// \param hfYN center index of the image along Y cordinate
/// \param angs projection views
/// \param blk Block size configuration for projection in GPU
/// \param gid Grid size configuration for projection in GPU
void DDM_ED_bproj_GPU(const float* proj, float* img, const float S2O, const float O2D, 
const float objSizeX, const float objSizeY, const float detSize, const float detCntIdx, 
const int XN, const int YN, const int DN, const int PN,	const float dd, const float dx, const float dy, 
const float hfXN, const float hfYN, const float S2D, const float* angs, const dim3 blk, const dim3 gid);

/// \brief Equal Distance DD Backprojection in double floating point precision 3D Geometry
/// \param proj pointer to projection data
/// \param img pointer to image volume
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object Size in X direction
/// \param objSizeY object Size in Y direction
/// \param objSizeZ object Size in Z direction
/// \param detSizeU detector size along detector channel direction
/// \param detSizeV detector size along detector bench moving direction
/// \param detCntIdU detector center index along channel direction
/// \param detCntIdV detector center index along bench moving direction
/// \param XN number of image pixels in X direction
/// \param YN number of image pixels in Y direction
/// \param ZN number of image pixels in Y direction
/// \param DNU number of detector cells along channel direction
/// \param DNV number of detector cells along bench moving direction
/// \param PN number of views (may be eliminated)
/// \param ddu detector cell size along channel direction
/// \param ddv detector cell size along bench moving direction
/// \param dx size of the pixel along x direction
/// \param dy size of the pixel along y direction
/// \param dz size of the pixel along z direction
/// \param hfXN center index of the image along X cordinate
/// \param hfYN center index of the image along Y cordinate
/// \param hfZN center index of the image along Z cordinate
/// \param S2D source to detector distance
/// \param angs projection views
/// \param blk Block size configuration for projection in GPU
/// \param gid Grid size configuration for projection in GPU
void DDM3D_ED_bproj_GPU(const double* proj, double* vol, const double S2O, const double O2D,
const double objSizeX, const double objSizeY, const double objSizeZ, const double detSizeU, const double detSizeV,
const double detCntIdU, const double detCntIdV, const int XN, const int YN, const int ZN, const int DNU, const int DNV, 
const int PN, const double ddu, const double ddv, const double dx, const double dy, const double dz,
	const double hfXN, const double hfYN, const double hfZN, const double S2D,
	const double* angs, const dim3 blk, const dim3 gid);

/// \brief Equal Distance DD Backprojection in double floating point precision 3D Geometry
/// \param proj pointer to projection data
/// \param img pointer to image volume
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object Size in X direction
/// \param objSizeY object Size in Y direction
/// \param objSizeZ object Size in Z direction
/// \param detSizeU detector size along detector channel direction
/// \param detSizeV detector size along detector bench moving direction
/// \param detCntIdU detector center index along channel direction
/// \param detCntIdV detector center index along bench moving direction
/// \param XN number of image pixels in X direction
/// \param YN number of image pixels in Y direction
/// \param ZN number of image pixels in Y direction
/// \param DNU number of detector cells along channel direction
/// \param DNV number of detector cells along bench moving direction
/// \param PN number of views (may be eliminated)
/// \param ddu detector cell size along channel direction
/// \param ddv detector cell size along bench moving direction
/// \param dx size of the pixel along x direction
/// \param dy size of the pixel along y direction
/// \param dz size of the pixel along z direction
/// \param hfXN center index of the image along X cordinate
/// \param hfYN center index of the image along Y cordinate
/// \param hfZN center index of the image along Z cordinate
/// \param S2D source to detector distance
/// \param angs projection views
/// \param blk Block size configuration for projection in GPU
/// \param gid Grid size configuration for projection in GPU
void DDM3D_ED_bproj_GPU(const float* proj, float* vol, const float S2O, const float O2D, 
const float objSizeX, const float objSizeY, const float objSizeZ, const float detSizeU, const float detSizeV,
const float detCntIdU, const float detCntIdV, const int XN, const int YN, const int ZN, const int DNU, const int DNV, 
const int PN, const float ddu, const float ddv, const float dx, const float dy, const float dz,
const float hfXN, const float hfYN, const float hfZN, const float S2D,	const float* angs, const dim3 blk, const dim3 gid);

/// \brief Equal Distance DD Backprojection in double floating datatype in Fan Beam Geometry
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
void DDM_ED_bproj(const std::vector<double>& proj, std::vector<double>& img, const double S2O, const double O2D,
const double objSizeX, const double objSizeY, const double detSize, const double detCntIdx, 
const int XN, const int YN, const int DN, const int PN,	const std::vector<double>& angs);

/// \brief Equal Distance DD Backprojection in single floating datatype in Fan Beam Geometry
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
void DDM_ED_bproj(const std::vector<float>& proj, std::vector<float>& img, const float S2O, const float O2D, 
const float objSizeX, const float objSizeY, const float detSize, const float detCntIdx, 
const int XN, const int YN, const int DN, const int PN,	const std::vector<float>& angs);

/// \brief Equal Angular DD Backprojection in double floating Fan Beam Geometry
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
void DDM_EA_bproj(const std::vector<double>& proj, std::vector<double>& img, const double S2O, const double O2D,
const double objSizeX, const double objSizeY, const double detArc, const double detCntIdx, const int XN, const int YN,
const int DN, const int PN, const std::vector<double>& angs);

/// \brief Equal Angular DD Backprojection in double floating Fan Beam Geometry
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
void DDM_EA_bproj(const std::vector<float>& proj, std::vector<float>& img, const float S2O, const float O2D,
const float objSizeX, const float objSizeY, const float detArc, const float detCntIdx, const int XN, const int YN, 
const int DN, const int PN, const std::vector<float>& angs);

/// \brief Equal Distance DD Backprojection in Cone Beam Geometry in double floating datatype
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
void DDM3D_ED_bproj(const std::vector<double>& proj, std::vector<double>& vol, const double S2O, const double O2D,
const double objSizeX, const double objSizeY, const double objSizeZ, const double detSizeU, const double detSizeV,
const double detCntIdU, const double detCntIdV, const int XN, const int YN, const int ZN, const int DNU, const int DNV, 
const int PN, const std::vector<double>& angs);

/// \brief Equal Distance DD Backprojection in Cone Beam Geometry in single floating datatype
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
void DDM3D_ED_bproj(const std::vector<float>& proj, std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeU, const float detSizeV,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const std::vector<float>& angs);

/// \brief Equal Angular DD Backprojection in Cone Beam Geometry in single floating datatype 
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
void DDM3D_EA_bproj(const std::vector<float>& proj, std::vector<float>& vol, const float S2O, const float O2D,
const float objSizeX, const float objSizeY, const float objSizeZ, const float detArc, const float detSizeH,
const float detCntIdU, const float detCntIdV, const int XN, const int YN, const int ZN, 
const int DNU, const int DNV, const int PN, const std::vector<float>& angs);

/// \brief Equal Angular DD Backprojection in Cone Beam Geometry in double floating datatype
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
void DDM3D_EA_bproj(const std::vector<double>& proj, std::vector<double>& vol, const double S2O, const double O2D,
const double objSizeX, const double objSizeY, const double objSizeZ, const double detArc, const double detSizeH,
const double detCntIdU, const double detCntIdV, const int XN, const int YN, const int ZN, 
const int DNU, const int DNV, const int PN, const std::vector<double>& angs);

/// \brief Equal Angular DD Backprojection in Helical Beam Geometry in single floating datatype
/// \param proj projection data
/// \param vol image volume
/// \param S2O source to object distance
/// \param O2D object to detector distance
/// \param objSizeX object Size in X direction
/// \param objSizeY object Size in Y direction
/// \param objSizeZ object Size in Z direction
/// \param detArc detector arc along xy plane
/// \param detSizeH detector size along z direction
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
void DDM3D_EA_helical_bproj(const std::vector<float>& proj, std::vector<float>& vol, const float S2O, const float O2D,
const float objSizeX, const float objSizeY, const float objSizeZ, const float detArc, const float detSizeH,
const float detCntIdU, const float detCntIdV, const int XN, const int YN, const int ZN, const int DNU, const int DNV,
const int PN, const float initZPos, const float pitch, const std::vector<float>& angs);

/// \brief Equal Angular DD Backprojection in Helical Beam Geometry in double floating datatype
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
void DDM3D_EA_helical_bproj(const std::vector<double>& proj, std::vector<double>& vol, const double S2O, const double O2D,
const double objSizeX, const double objSizeY, const double objSizeZ, const double detArc, const double detSizeH,
const double detCntIdU, const double detCntIdV,	const int XN, const int YN, const int ZN, const int DNU, const int DNV, 
const int PN, const double initZPos, const double pitch, const std::vector<double>& angs);

/// \brief Equal distance DD Backprojection in Fan Beam Geometry in GPU (all data are all in GPU)
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
void DDM_ED_bproj_GPU(const double* proj, double* img, const double S2O, const double O2D, 
const double objSizeX, const double objSizeY, const double detSize, const double detCntIdx, 
const int XN, const int YN, const int DN, const int PN, const double dd, const double dx, const double dy, 
const double hfXN, const double hfYN, const double S2D, const double* angs, const dim3 blk, const dim3 gid);

/// \brief Equal distance DD Backprojection in Fan Beam Geometry in GPU (all data are all in GPU)
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
void DDM_ED_bproj_GPU(const float* proj, float* img,
	const float S2O, const float O2D, const float objSizeX, const float objSizeY,
	const float detSize, const float detCntIdx, const int XN, const int YN, const int DN, const int PN,
	const float dd, const float dx, const float dy, const float hfXN, const float hfYN,
	const float S2D,
	const float* angs, const dim3 blk, const dim3 gid);
void DDM_ED_bproj_GPU(const thrust::device_vector<double>& proj, thrust::device_vector<double>& img,
	const double S2O, const double O2D, const double objSizeX, const double objSizeY,
	const double detSize, const double detCntIdx, const int XN, const int YN, const int DN, const int PN,
	const double dd, const double dx, const double dy, const double hfXN, const double hfYN,
	const double S2D,
	const thrust::device_vector<double>& angs, const dim3 blk, const dim3 gid);
void DDM_ED_bproj_GPU(const thrust::device_vector<float>& proj, thrust::device_vector<float>& img,
	const float S2O, const float O2D, const float objSizeX, const float objSizeY,
	const float detSize, const float detCntIdx, const int XN, const int YN, const int DN, const int PN,
	const float dd, const float dx, const float dy, const float hfXN, const float hfYN,
	const float S2D,
	const thrust::device_vector<float>& angs, const dim3 blk, const dim3 gid);



/// \brief Equal distance DD Backprojection in Cone Beam Geometry in GPU (all data are all in GPU)
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
/// \param hfXN image center index along x
/// \param hfYN image center index along y
/// \param hfZN image center index along z
/// \param S2D source to detector distance
/// \param angs projection views
/// \param blk threads block size
/// \param gid block grid size
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
void DDM3D_ED_bproj_GPU(const thrust::device_vector<double>& proj, thrust::device_vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeU, const double detSizeV,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double ddu, const double ddv, const double dx, const double dy, const double dz,
	const double hfXN, const double hfYN, const double hfZN, const double S2D,
	const thrust::device_vector<double>& angs, const dim3 blk, const dim3 gid);
void DDM3D_ED_bproj_GPU(const thrust::device_vector<float>& proj, thrust::device_vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeU, const float detSizeV,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float ddu, const float ddv, const float dx, const float dy, const float dz,
	const float hfXN, const float hfYN, const float hfZN, const float S2D,
	const thrust::device_vector<float>& angs, const dim3 blk, const dim3 gid);

/// \brief Equal distance DD Backprojection in Helical Beam Geometry in GPU (all data are all in GPU)
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
/// \param hfXN image center index along x
/// \param hfYN image center index along y
/// \param hfZN image center index along z
/// \param S2D source to detector distance
/// \param initZPos initial Z position of the source/detector
/// \param pitch pitch size
/// \param angs projection views
/// \param blk threads block size
/// \param gid block grid size
void DDM3D_EA_helical_bproj(const float* proj, float* vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float dbeta, const float ddv, const float dx, const float dy, const float dz,
	const float hfXN, const float hfYN, const float hfZN, const float S2D,
	const float initZPos, const float pitch,
	const float*  angs, const dim3 blk, const dim3 gid);
void DDM3D_EA_helical_bproj(const double* proj, double* vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double dbeta, const double ddv, const double dx, const double dy, const double dz,
	const double hfXN, const double hfYN, const double hfZN, const double S2D,
	const double initZPos, const double pitch,
	const double*  angs, const dim3 blk, const dim3 gid);

/// \brief Equal distance DD Backprojection in Helical Beam Geometry in GPU (all data are all in GPU)
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
/// \param hfXN image center index along x
/// \param hfYN image center index along y
/// \param hfZN image center index along z
/// \param S2D source to detector distance
/// \param zShifts Z positions for each views
/// \param angs projection views
/// \param blk threads block size
/// \param gid block grid size
void DDM3D_EA_helical_bproj(const float* proj, float* vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float dbeta, const float ddv, const float dx, const float dy, const float dz,
	const float hfXN, const float hfYN, const float hfZN, const float S2D,
	const float* zShifts,
	const float*  angs, const dim3 blk, const dim3 gid);
void DDM3D_EA_helical_bproj(const double* proj, double* vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double dbeta, const double ddv, const double dx, const double dy, const double dz,
	const double hfXN, const double hfYN, const double hfZN, const double S2D,
	const double* zShifts,
	const double*  angs, const dim3 blk, const dim3 gid);

/// \brief Equal distance DD Backprojection in Helical Beam Geometry in GPU (all data are all in GPU)
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
/// \param hfXN image center index along x
/// \param hfYN image center index along y
/// \param hfZN image center index along z
/// \param S2D source to detector distance
/// \param initZPos initial Z position of the source/detector
/// \param pitch pitch size
/// \param angs projection views
/// \param blk threads block size
// \param gid block grid size
void DDM3D_EA_helical_bproj(const thrust::device_vector<float>& proj, thrust::device_vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float dbeta, const float ddv, const float dx, const float dy, const float dz,
	const float hfXN, const float hfYN, const float hfZN, const float S2D,
	const float initZPos, const float pitch,
	const thrust::device_vector<float>& angs, const dim3 blk, const dim3 gid);
void DDM3D_EA_helical_bproj(const thrust::device_vector<double>& proj, thrust::device_vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double dbeta, const double ddv, const double dx, const double dy, const double dz,
	const double hfXN, const double hfYN, const double hfZN, const double S2D,
	const double initZPos, const double pitch,
	const thrust::device_vector<double>&, const dim3 blk, const dim3 gid);

/// \brief Equal distance DD Backprojection in Helical Beam Geometry in GPU (all data are all in GPU)
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
/// \param hfXN image center index along x
/// \param hfYN image center index along y
/// \param hfZN image center index along z
/// \param S2D source to detector distance
/// \param zShifts Z positions for each views
/// \param angs projection views
/// \param blk threads block size
/// \param gid block grid size
void DDM3D_EA_helical_bproj(const thrust::device_vector<float>& proj, thrust::device_vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float dbeta, const float ddv, const float dx, const float dy, const float dz,
	const float hfXN, const float hfYN, const float hfZN, const float S2D,
	const thrust::device_vector<float>& zShifts,
	const thrust::device_vector<float>& angs, const dim3 blk, const dim3 gid);
void DDM3D_EA_helical_bproj(const thrust::device_vector<double>& proj, thrust::device_vector<double>&  vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double dbeta, const double ddv, const double dx, const double dy, const double dz,
	const double hfXN, const double hfYN, const double hfZN, const double S2D,
	const thrust::device_vector<double>& zShifts,
	const thrust::device_vector<double>& angs, const dim3 blk, const dim3 gid);

#endif /* _DDM_BACK_H_ */
