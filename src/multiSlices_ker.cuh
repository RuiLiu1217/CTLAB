/*!
 * 
 * \file multiSlices_ker.cuh
 * The multi-slices based projection/backprojection model
 * \date Sep 19, 2016
 * \author liurui
 * \version 1.0
 */

#ifndef MULTISLICES_KER_CUH_
#define MULTISLICES_KER_CUH_

typedef unsigned char byte;

/// \brief The multi slices geometry based projection and backprojection routine with multi-GPUs
/// \param hvol The pointer to the image
/// \param hprj The pointer to the projection data
/// \param method Control to use forward projection or backprojection 0:branchless projection, 1:pseudo projection, 2: branchless backprojection, 3: pseudo backprojection
/// \param x0 initial X coordinate of the source
/// \param y0 initial Y coordinate of the source
/// \param xds pointer to the initial X coordinates of the detector cells
/// \param yds pointer to the initial Y coordinates of the detector cells
/// \param DNU number of detector cells along channel direction
/// \param SLN number of slices to be reconstructed
/// \param imgXCenter center of the image X coordinate
/// \param imgYCenter center of the image Y Coordinate
/// \param XN pixel number of the image along X direction
/// \param YN pixel number of the image along Y direction
/// \param dx size of the pixel
/// \param hangs pointer to all projection views
/// \param PN number of views
/// \param mask pointer to image mask
/// \param startIdx start slice index for each GPU
/// \param gpuNum number of GPUs to be applied
extern "C"
void DD2_multiGPU(float* hvol, float* hprj, const int method, const float x0, const float y0,
	float* xds, float* yds, const int DNU, const int SLN, const float imgXCenter, const float imgYCenter, 
	const int XN, const int YN, const float dx, float* hangs, int PN, byte* mask, int* startIdx, const int gpuNum);
#endif /* MULTISLICES_KER_CUH_ */
