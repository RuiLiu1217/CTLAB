/*
 *  DD3_GPU_Back.h
 *  Wake Forest Health Sciences & University of Massachusetts Lowell
 *  Created on: Oct 2, 2015
 *  Author: Rui Liu
 *  Email: liurui1217@gmail.com
 */

#ifndef DD3_GPU_BACK_H_
#define DD3_GPU_BACK_H_

typedef unsigned char byte;


/// \brief Backprojection interface in one GPU
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
/// \param squared A fake parameter
/// \param prjMode The backprojection mode: 0:branchless DD; 1: volume rendering; 2: pseudo DD; 3: z-line based branchless DD
extern "C"
void DD3Back_gpu(float x0, float y0, float z0, int DNU, int DNV, 
float* xds, float* yds, float* zds, float imgXCenter, float imgYCenter, float imgZCenter,
float* hangs, float* hzPos, int PN, int XN, int YN, int ZN, float* hvol, float* hprj,
float dx, float dz, byte* mask, int gpunum, int squared, int prjMode); 

/// \brief Backprojection interface in three GPUs
/// 
/// This backprojection function routine is old.
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
/// \param methodId The backprojection mode: 0:branchless DD; 1: volume rendering; 2: pseudo DD; 3: z-line based branchless DD
/// \param startVOL Start backprojection index of the slices for each GPU
extern "C"
void DD3BackHelical_3GPU(float x0, float y0, float z0, int DNU, int DNV, float* xds, float* yds, float* zds, 
	float imgXCenter, float imgYCenter, float imgZCenter, float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN,	float* hvol, float* hprj, float dx, float dz, byte* mask, int methodId, int (&startVOL)[3]);

/// \brief Backprojection interface in four GPUs
/// 
/// This backprojection function routine is old.
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
/// \param methodId The backprojection mode: 0:branchless DD; 1: volume rendering; 2: pseudo DD; 3: z-line based branchless DD
/// \param startVOL Start backprojection index of the slices for each GPU
extern "C"
void DD3BackHelical_4GPU(float x0, float y0, float z0, int DNU, int DNV, float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter, float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN, float* hvol, float* hprj, float dx, float dz, byte* mask, int methodId, int (&startVOL)[4]);
	
/// \brief Backprojection interface in multiple GPUs
/// 
/// This backprojection function routine is new.
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
/// \param methodId The backprojection mode: 0:branchless DD; 1: volume rendering; 2: pseudo DD; 3: z-line based branchless DD
/// \param startVOL Start backprojection index of the slices for each GPU
/// \param gpuNum Number of GPU to be applied
extern "C"
void DD3Back_multiGPU(float x0, float y0, float z0, int DNU, int DNV, float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter, float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN, float* hvol, float* hprj, float dx, float dz, byte* mask, int methodId, 
        int* startVOL, int gpuNum);

#endif /* DD3_GPU_BACK_H_ */

