/*
 * multiSlices_ker.cuh
 *
 *  Created on: Sep 19, 2016
 *      Author: liurui
 */

#ifndef MULTISLICES_KER_CUH_
#define MULTISLICES_KER_CUH_

typedef unsigned char byte;

extern "C"
void DD2_multiGPU(
		float* hvol, // the pointer to the image
		float* hprj, // the pointer to the projection (SLN, DNU, PN) order
		const int method, // Control to use forward projection or backprojection
		const float x0, const float y0, //position of the initial source
		float* xds, float* yds, // distribution of the detector cells
		const int DNU, // Number of detector cells
		const int SLN, // Number of slices to be projected or backprojected
		const float imgXCenter, const float imgYCenter, //Center of the image
		const int XN, const int YN, // pixel number of the image
		const float dx, // size of the pixel
		float* hangs, // view angles (the size should be SLN * PN)
		int PN, // # of view angles
		byte* mask,
		int* startIdx,
		const int gpuNum);

#endif /* MULTISLICES_KER_CUH_ */
