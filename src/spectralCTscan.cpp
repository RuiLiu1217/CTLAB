#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include "cublas_v2.h"
#include "spectralCTscan.hpp"
namespace SPECTRALCT{
//Test maanscanning
int test1()
{
	maanscanCylinder();
	return 0;
}

//Test conebeam scan the cylinder in multiEnergies
int test2()
{
	cublasHandle_t handle;
	cublasCreate(&handle);
	const double initx = 0;
	const double inity = 5.0;
	const double initz = -0.5;
	const double initdetcntx = 0;
	const double initdetcnty = -5.0;
	const double initdetcntz = -0.5;
	const double detustp = 0.008;
	const double detvstp = 0.008;
	const int DNi = 512;
	const int DNj = 512;
	const int angN = 360;
	const double angStp = 0.01745329251994329576923688888889 * 2;

	conebeamScanCylinderMultiEnergiesV4(handle, initx, inity, initz, initdetcntx, initdetcnty, initdetcntz,
		detustp, detvstp, DNi, DNj, angN, angStp);
	return 0;
}


//Test conebeam scanning the head phantom in multi energies
int test3()
{
	cublasHandle_t handle;
	cublasCreate(&handle);
	const double initx = 0;
	const double inity = 5.0;
	const double initz = -0.5;
	const double initdetcntx = 0;
	const double initdetcnty = -5.0;
	const double initdetcntz = -0.5;
	const double detustp = 0.008;
	const double detvstp = 0.008;
	const int DNi = 512;
	const int DNj = 512;
	const int angN = 360;
	const double angStp = 0.01745329251994329576923688888889 * 2;
	const double pitch = 2;
	conebeamScanHeadPhantomMultiEnergies(handle, initx, inity, initz, initdetcntx, initdetcnty, initdetcntz, detustp, detvstp, DNi, DNj, angN, angStp);

	
	return 0;
}


//Test sprial beam scanning the head phantom in multi energies
int test4()
{
	cublasHandle_t handle;
	cublasCreate(&handle);
	const double initx = 0;
	const double inity = 5.0;
	const double initz = -0.5;
	const double initdetcntx = 0;
	const double initdetcnty = -5.0;
	const double initdetcntz = -0.5;
	const double detustp = 0.008;
	const double detvstp = 0.008;
	const int DNi = 512;
	const int DNj = 512;
	const int angN = 360;
	const double angStp = 0.01745329251994329576923688888889 * 2;
	const double pitch = 2;
	helicalScanHeadPhantomMultiEnergies(handle, initx, inity, initz, initdetcntx, initdetcnty, initdetcntz, pitch, detustp, detvstp, DNi, DNj, angN, angStp);

	
	return 0;
}

}


