/*
 * DD3_GPU_Proj.h
 * The GPU based branchless DD projection routine
 * Created on: Oct 2, 2015
 * Author: Rui Liu
 * Email: liurui1217@gmail.com
 * Wake Forest Health Sciences & University of Massachusetts Lowell
 */

#ifndef _DD3_GPU_PROJ_H_
#define _DD3_GPU_PROJ_H_
#include <thrust/device_vector.h>

typedef unsigned char byte;


/* 
 * The Branchless DD GPU projection interface in one GPU
 */
extern "C"
void DD3Proj_gpu(
float x0, float y0, float z0, // source position
int DNU, int DNV, // # of detector cells along channel and bench moving directions
float* xds, float* yds, float* zds, // the coordinates of the detector cells at initial positions
float imgXCenter, float imgYCenter, float imgZCenter, // image center
float* hangs, float* hzPos,  // view angle and corresponding source position
int PN,  // # of views
int XN, int YN, int ZN, // # of pixels along x,y,z directions
float* hvol, float* hprj, // image volume, projection data
float dx, float dz, // size of pixels along xy direction and z direction
byte* mask,  // image mask along xy plane
int gpunum,  // which gpu is to use
int prjMode); // projection mode 0:branchless DD, 1:volume rendering, 2,3:pseudo branchless DD, 4: siddon's method (DO NOT USE), 5:brute forece dd (DO NOT USE)



//Use the split-collect method to do the projection
extern "C"
void DD3ProjHelical_3GPU(
	float x0, float y0, float z0,
	int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN,
	float* hvol, float* hprj,
	float dx, float dz,
	byte* mask, int methodId, int (&startPN)[3]);

extern "C"
void DD3ProjHelical_4GPU(
	float x0, float y0, float z0,
	int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN,
	float* hvol, float* hprj,
	float dx, float dz,
	byte* mask, int methodId, int (&startPN)[4]);



#endif
