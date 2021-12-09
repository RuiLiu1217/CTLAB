/*!
 * \file DD3_GPU_Proj.h
 * \brief The GPU based branchless DD projection routine
 * Wake Forest Health Sciences & University of Massachusetts Lowell
 * E-mail liurui1217@gmail.com
 * \date Oct 2, 2015
 * \author Rui Liu
 * \version 1.0 
 */
#ifndef _DD3_GPU_PROJ_H_
#define _DD3_GPU_PROJ_H_
#include "Image.h"
#include "Projection.h"
#include <thrust/device_vector.h>

using byte = unsigned char;

/// \brief Projection interface in one GPU
/// \param x0 The initial X coordinate of the source position
/// \param y0 The initial Y coordinate of the source position
/// \param z0 The initial Z coordinate of the source position
/// \param DNU The number of detector cells along channel direction
/// \param DNV The number of detector cells along bench moving direction
/// \param xds The pointer to the X coordinates distribution of the detector cell at the initial position (SIZE=DNU)
/// \param yds The pointer to the Y coordinates distribution of the detector cell at the initial position (SIZE=DNU)
/// \param zds The pointer to the Z coordinates distribution of the detector cell at the initial position (SIZE=DNV)
/// \param imgXCenter The X coordinate of the center of the image in world coordinates
/// \param imgYCenter The Y coordinate of the center of the image in world coordinates
/// \param imgZCenter The Z coordinate of the center of the image in world coordinates
/// \param hangs The pointer to all views of the projection
/// \param hzPos The pointer to all z coordinates of the source during sampling
/// \param PN Number of views
/// \param XN Number of pixels along X direction
/// \param YN Number of pixels along Y direction
/// \param ZN Number of pixels along Z direction
/// \param hvol The poniter to the image volume
/// \param hprj The pointer to the projection data
/// \param dx The size of the pixel along xy plane
/// \param dz The size of the pixel along z direction
/// \param mask The pointer to the image mask for reconstruction along xy plane
/// \param gpuNum The ID of the GPU to be used
/// \param prjMode Projection mode 0:branchless DD, 1:volume rendering, 2,3:pseudo branchless DD, 4: siddon's method (DO NOT USE), 5:brute forece dd (DO NOT USE)
extern "C"
void DD3Proj_gpu(float x0, float y0, float z0, int DNU, int DNV, float* xds, float* yds, float* zds,
float imgXCenter, float imgYCenter, float imgZCenter, float* hangs, float* hzPos, int PN, 
int XN, int YN, int ZN, float* hvol, float* hprj, float dx, float dz, byte* mask, int gpunum, int prjMode);

extern "C"
void Proj(Image & image, Projection & projection, byte * mask, int gpunum, int prjMode);

/// \brief Projection interface in three GPUs
/// 
/// This projection routine is old
/// \param x0 The initial X coordinate of the source position
/// \param y0 The initial Y coordinate of the source position
/// \param z0 The initial Z coordinate of the source position
/// \param DNU The number of detector cells along channel direction
/// \param DNV The number of detector cells along bench moving direction
/// \param xds The pointer to the X coordinates distribution of the detector cell at the initial position (SIZE=DNU)
/// \param yds The pointer to the Y coordinates distribution of the detector cell at the initial position (SIZE=DNU)
/// \param zds The pointer to the Z coordinates distribution of the detector cell at the initial position (SIZE=DNV)
/// \param imgXCenter The X coordinate of the center of the image in world coordinates
/// \param imgYCenter The Y coordinate of the center of the image in world coordinates
/// \param imgZCenter The Z coordinate of the center of the image in world coordinates
/// \param hangs The pointer to all views of the projection
/// \param hzPos The pointer to all z coordinates of the source during sampling
/// \param PN Number of views
/// \param XN Number of pixels along X direction
/// \param YN Number of pixels along Y direction
/// \param ZN Number of pixels along Z direction
/// \param hvol The poniter to the image volume
/// \param hprj The pointer to the projection data
/// \param dx The size of the pixel along xy plane
/// \param dz The size of the pixel along z direction
/// \param mask The pointer to the image mask for reconstruction along xy plane
/// \param methodId Projection mode 0:branchless DD, 1:volume rendering, 2,3:pseudo branchless DD, 4: siddon's method (DO NOT USE), 5:brute forece dd (DO NOT USE)
/// \param startPN The start index of the view for each GPU
extern "C"
void DD3ProjHelical_3GPU(float x0, float y0, float z0, int DNU, int DNV, float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter, float* hangs, float* hzPos, int PN, int XN, int YN, int ZN,
	float* hvol, float* hprj, float dx, float dz, byte* mask, int methodId, int (&startPN)[3]);

/// \brief Projection interface in four GPUs
/// 
/// This projection routine is old
/// \param x0 The initial X coordinate of the source position
/// \param y0 The initial Y coordinate of the source position
/// \param z0 The initial Z coordinate of the source position
/// \param DNU The number of detector cells along channel direction
/// \param DNV The number of detector cells along bench moving direction
/// \param xds The pointer to the X coordinates distribution of the detector cell at the initial position (SIZE=DNU)
/// \param yds The pointer to the Y coordinates distribution of the detector cell at the initial position (SIZE=DNU)
/// \param zds The pointer to the Z coordinates distribution of the detector cell at the initial position (SIZE=DNV)
/// \param imgXCenter The X coordinate of the center of the image in world coordinates
/// \param imgYCenter The Y coordinate of the center of the image in world coordinates
/// \param imgZCenter The Z coordinate of the center of the image in world coordinates
/// \param hangs The pointer to all views of the projection
/// \param hzPos The pointer to all z coordinates of the source during sampling
/// \param PN Number of views
/// \param XN Number of pixels along X direction
/// \param YN Number of pixels along Y direction
/// \param ZN Number of pixels along Z direction
/// \param hvol The poniter to the image volume
/// \param hprj The pointer to the projection data
/// \param dx The size of the pixel along xy plane
/// \param dz The size of the pixel along z direction
/// \param mask The pointer to the image mask for reconstruction along xy plane
/// \param methodId Projection mode 0:branchless DD, 1:volume rendering, 2,3:pseudo branchless DD, 4: siddon's method (DO NOT USE), 5:brute forece dd (DO NOT USE)
/// \param startPN The start index of the view for each GPU
extern "C"
void DD3ProjHelical_4GPU(float x0, float y0, float z0, int DNU, int DNV, float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter, float* hangs, float* hzPos, int PN, int XN, int YN, int ZN,
	float* hvol, float* hprj, float dx, float dz, byte* mask, int methodId, int (&startPN)[4]);
	
/// \brief Projection interface in multi-GPUs
/// 
/// This projection routine is new
/// \param x0 The initial X coordinate of the source position
/// \param y0 The initial Y coordinate of the source position
/// \param z0 The initial Z coordinate of the source position
/// \param DNU The number of detector cells along channel direction
/// \param DNV The number of detector cells along bench moving direction
/// \param xds The pointer to the X coordinates distribution of the detector cell at the initial position (SIZE=DNU)
/// \param yds The pointer to the Y coordinates distribution of the detector cell at the initial position (SIZE=DNU)
/// \param zds The pointer to the Z coordinates distribution of the detector cell at the initial position (SIZE=DNV)
/// \param imgXCenter The X coordinate of the center of the image in world coordinates
/// \param imgYCenter The Y coordinate of the center of the image in world coordinates
/// \param imgZCenter The Z coordinate of the center of the image in world coordinates
/// \param hangs The pointer to all views of the projection
/// \param hzPos The pointer to all z coordinates of the source during sampling
/// \param PN Number of views
/// \param XN Number of pixels along X direction
/// \param YN Number of pixels along Y direction
/// \param ZN Number of pixels along Z direction
/// \param hvol The poniter to the image volume
/// \param hprj The pointer to the projection data
/// \param dx The size of the pixel along xy plane
/// \param dz The size of the pixel along z direction
/// \param mask The pointer to the image mask for reconstruction along xy plane
/// \param prjMode Projection mode 0:branchless DD, 1:volume rendering, 2,3:pseudo branchless DD, 4: siddon's method (DO NOT USE), 5:brute forece dd (DO NOT USE)
/// \param startPN The start index of the view for each GPU
/// \param gpuNum Number of GPUs to be applied
extern "C"
void DD3Proj_multiGPU(float x0, float y0, float z0, int DNU, int DNV, float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter, float* hangs, float* hzPos, int PN, int XN, int YN, int ZN,
	float* hvol, float* hprj, float dx, float dz, byte* mask, int prjMode, int* startPN, int gpuNum);

#endif
