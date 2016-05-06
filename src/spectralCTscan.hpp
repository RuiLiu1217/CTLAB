#include "cublas_v2.h"
namespace SPECTRALCT{
void conebeamScanCylinderMultiEnergiesV4(
cublasHandle_t handle, const double initx, const double inity, const double initz,
const double initdetcntx, const double initdetcnty, const double initdetcntz,
const double detustp, const double detvstp,
const int DNi, const int DNj, const int angN, const double angStp);

void conebeamScanHeadPhantomMultiEnergies(
cublasHandle_t handle, const double initx, const double inity, const double initz,
const double initdetcntx, const double initdetcnty, const double initdetcntz,
const double detustp, const double detvstp,
const int DNi, const int DNj, const int angN, const double angStp);


void helicalScanHeadPhantomMultiEnergies(
cublasHandle_t handle, const double initx, const double inity, const double initz,
const double initdetcntx, const double initdetcnty, const double initdetcntz, const double pitch,
const double detustp, const double detvstp,
const int DNi, const int DNj, const int angN, const double angStp);


void maanscanCylinder();

int test1();
int test2();
int test3();
int test4();

// 这个函数才是正确的，上面的所有函数都有问题;
// 从801个通道合成一个通道，并输出8个通道的值;
// 这个将各个通道的光子大约弄成一样多的方式;
void conebeamScanCylinderMultiEnergiesNew(
	cublasHandle_t handle, const double initx, const double inity, const double initz,
	const double initdetcntx, const double initdetcnty, const double initdetcntz,
	const double detustp, const double detvstp,
	const int DNi, const int DNj, const int angN, const double angStp);

// 这个函数实际上没有数据材料数据，因为都是保存成double的，但是如果有相应数据还是可以用的;
void conebeamScanCylinderMultiEnergiesNew(
	cublasHandle_t handle, const float initx, const float inity, const float initz,
	const float initdetcntx, const float initdetcnty, const float initdetcntz,
	const float detustp, const float detvstp,
	const int DNi, const int DNj, const int angN, const float angStp);
void testNewGeneratePhotonNumProj(); //Conebeam scanning	
	
}
