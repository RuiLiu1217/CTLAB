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


// Backprojection interface in one GPU
extern "C"
void DD3Back_gpu(
float x0, float y0, float z0, // initial source position
int DNU, int DNV, // # of detector cells along channel and bench moving direction
float* xds, float* yds, float* zds, // the initial detector cell positions // cylinderal detector assumed
float imgXCenter, float imgYCenter, float imgZCenter, // center of the image
float* hangs, float* hzPos, // views and source z positions in different views
int PN, // #number of projection views
int XN, int YN, int ZN, // pixel number 
float* hvol, float* hprj, // image volume and the projection data
float dx, float dz, // size of the pixel along xy direction and z direction
byte* mask, // image mask along xy plane
int gpunum, // which gpu is to use
int squared, // NOT implemented (FOR Preconditioning conjugate gradient PCG)
int prjMode); // 0:branchless DD; 1: volume rendering; 2: pseudo DD; 3: z-line based branchless DD



extern "C"
void DD3BackHelical_3GPU(
	float x0, float y0, float z0,
	int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN,
	float* hvol, float* hprj,
	float dx, float dz,
	byte* mask, int methodId, int (&startVOL)[3]);


extern "C"
void DD3BackHelical_4GPU(
	float x0, float y0, float z0,
	int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN,
	float* hvol, float* hprj,
	float dx, float dz,
	byte* mask, int methodId, int (&startVOL)[4]);
	
	
extern "C"
void DD3Back_multiGPU(
	float x0, float y0, float z0,
	int DNU, int DNV,
	float* xds, float* yds, float* zds,
	float imgXCenter, float imgYCenter, float imgZCenter,
	float* hangs, float* hzPos, int PN,
	int XN, int YN, int ZN,
	float* hvol, float* hprj,
	float dx, float dz,
	byte* mask, int bakMode, int* startVOL, int gpuNum);

#endif /* DD3_GPU_BACK_H_ */

