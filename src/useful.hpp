/*
 * COPYRIGHT NOTICE
 * COPYRIGHT (c) 2015, Wake Forest and UMass Lowell
 * All rights reserved
 *
 * \file useful.hpp
 * \brief The useful functions declaring. Including read and write functions called by users. Also
 * it includes some useful classes that used in reconstruction process. This file mainly includes 
 * some useful functions that can be called by the users. The author strongly suggest the users do 
 * NOT modify any codes but you can add your own functions and classes.  
 *
 * \version 1.0
 * \author Rui Liu
 * \date Aug. 24, 2014
 *
 */
#pragma once

#include "utilities.hpp"

void generateModiPhantom(const std::string& FileName, const unsigned int& lenReso = 256, const unsigned int& widReso = 256, const unsigned int& heiReso = 256);

void generateModiPhantom_d(const std::string& FileName, const unsigned int& lenReso = 256, const unsigned int& widReso = 256, const unsigned int& heiReso = 256);


/// \brief Transform the number to string
/// \param num input the "int" value that will be expressed in string form
/// \return the string form of the number
std::string num2String(int& num);
/// \brief Transform the number to string
/// \param num input the "unsigned int" value that will be expressed in string form
/// \return the string form of the number
std::string num2String(unsigned int& num);
/// \brief Transform the number to string
/// \param num input the "const int" value that will be expressed in string form
/// \return the string form of the number
std::string num2String(const int& num);
/// \brief Transform the number to string
/// \param num input the "const unsigned int" value that will be expressed in string form
/// \return the string form of the number
std::string num2String(const unsigned int& num);

/// \brief With the input index I, we get the output i, j, k
/// \param i input index of the total image
/// \param objL the length direction resolution of the image
/// \param iL output the L direction index
/// \param iW output the W direction index
inline void getIJ(cuint& i, cuint& objL, uint& iL, uint& iW){ iL = i % objL;	iW = i / objL; }
/// \brief With the input index I, we get the output i, j, k
/// \param i input index of the total image
/// \param objL the length direction resolution of the image
/// \param objW the width direction resolution of the image
/// \param iL output the L direction index
/// \param iW output the W direction index
/// \param iH output the H direction index
inline void getIJK(cuint& i, cuint& objL, cuint& objW, uint& iL, uint& iW, uint& iH){ cuint objWL(objL * objW);	iH = i / objWL;	uint iwl(i % objWL);	iW = iwl / objL;	iL = iwl % objL; }
/// \brief With the input the index of the image in two/three directions, we get the linear index of the image
/// \param i input length direction index
/// \param j input width direction index
/// \param L length direction resolution
/// \return linear index of the 2D image
inline unsigned int getIndex(cuint& i, cuint& j, cuint& L){ return j * L + i; }
/// \brief With the input the index of the image in two/three directions, we get the linear index of the image
/// \param i input length direction index
/// \param j input width direction index
/// \param k input height direction index
/// \param L length direction resolution
/// \param W width direction resolution
/// \return linear index of the 3D image
inline unsigned int getIndex(cuint& i, cuint& j, cuint& k, cuint& L, cuint& W){ return (k * W + j) * L + i; }



/// \brief Transform the Cartesian coordinates to Polar Coordinate
/// \param x Cartesian coordinate
/// \param y Cartesian coordinate
/// \param r radius of the vector in Polar coordinate
/// \param theta the angle in Polar coordinate
void toPolarCoords(const float& x, const float& y, float& r, float& theta);
/// \brief Transform the Cartesian coordinates to Polar Coordinate
/// \param x Cartesian coordinate
/// \param y Cartesian coordinate
/// \param z Cartesian coordinate
/// \param r radius of the vector in Polar coordinate
/// \param s Polar angle
/// \param t Elevation angle
void toSphericalCoords(const float& x, const float& y, const float& z, float& r, float& s, float& t);

/// \brief the cross product of two 3D vector
double3 crossProduct(const double3& a, const double3& b);
/// \brief the cross product of two 3D vector
float3 crossProduct(const float3& a, const float3& b);


//
//
///// \brief Generate the projection matrix in COO format with the given geometry configuration
///// \param FanGeo The geometry configuration of the CT system with equal angular detector
///// \param Img The Image parameters correspondingly
//void GenerateProjectionMatrix(const FanEAGeo& FanGeo, const Image& Img);
///// \brief Generate the projection matrix in COO format with the given geometry configuration
///// \param FanGeo The geometry configuration of the CT system with equal distance detector
///// \param Img The Image parameters correspondingly
//void GenerateProjectionMatrix(const FanEDGeo& FanGeo, const Image& Img);
//
//
////void GenerateProjectionMatrix_OPENMP(const FanEAGeo& FanGeo, const Image& Img, std::string FileName);
////void GenerateProjectionMatrix_OPENMP(const FanEDGeo& FanGeo, const Image& Img, std::string FileName);
//
//
///// \brief It combines the projection and backprojection processes in one function with the ray driven projection process and pixel driven backprojection process. This function will be modified in the future. The function combines \f$A^{T}\cdot A\f$ for Split-Bregman method
///// \param dimg the image pointer in device memory
///// \param dres the image be projected and then backprojection result
///// \param FanGeo the Fan Beam geometry with equal angular detector
///// \param Img the Image configuration
//void ProjectionThenBackProjection(float* dimg, float* dres, const FanEAGeo& FanGeo, const Image& Img);
///// \brief It combines the projection and backprojection processes in one function with the ray driven projection process and pixel driven backprojection process. This function will be modified in the future. The function combines \f$A^{T}\cdot A\f$ for Split-Bregman method
///// \param dimg the image pointer in device memory
///// \param dres the image be projected and then backprojection result
///// \param FanGeo the Fan Beam geometry with equal distance detector
///// \param Img the Image configuration
//void ProjectionThenBackProjection(float* dimg, float* dres, const FanEDGeo& FanGeo, const Image& Img);
//


/// \brief Calculate the p norm of a device vector, it cannot be used for host pointer
/// \param v the device pointer pointing to the beginning of an array
/// \param len the length of the vector (how many elements are there)
/// \param p the norm of the vector
/// \return we get the value \f$\left(\sum_{i=1}^{len}|v_{i}|^{p}\right)^{1/p}\f$
float norm(float* v, cuint len, const float p);
/// \brief Calculate the inner-product of two vectors. These two vectors usually are located on device memory.
/// \param a the first vector pointer
/// \param b the second vector pointer
/// \param len the length of the vector (how many elements are there)
/// \return the innerproduct of two vectors \f$<a,b>\f$
float innerProd(float* a, float* b, cuint len);
/// \brief Calculate the inner-product of two vectors. These two vectors usually are located on device memory.
/// \param a the first vector pointer
/// \param b the second vector pointer
/// \param L the length of the vector (how many elements are there)
/// \return the innerproduct of two vectors \f$<a,b>\f$
double innerProd(double* a, double* b, cuint L);

/// \brief The derivation of the TV norm, it is implemented in GPU, more details can be get from http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=6857986
/// \param f the pointer to the image
/// \param d the output image storing the TV derivation
/// \param coef the coefficient for transformed d = d * coef
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param blk the block configuration of the image in CUDA
/// \param gid the grid configuration of the image in CUDA
void GradientOfTV(float* f, float* d, const float coef, cuint L, cuint W, const dim3& blk, const dim3& gid); //use this one;
/// \brief The derivation of the TV norm, it is implemented in GPU, more details can be get from "GPU-Based Acceleration for Interior Tomography"
/// \param f the pointer to the image
/// \param d the output image storing the TV derivation
/// \param coef the coefficient for transformed d = d * coef
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param blk the block configuration of the image in CUDA
/// \param gid the grid configuration of the image in CUDA
void GradientOfTV(double* f, double* d, const double coef, cuint L, cuint W, const dim3& blk, const dim3& gid);
/// \brief The derivation of the TV norm, it is implemented in GPU, it is based on SHARED memory, and it theoretically more effective, more details can be get from "GPU-Based Acceleration for Interior Tomography"
/// \param f the pointer to the image
/// \param d the output image storing the TV derivation
/// \param coef the coefficient control the magnitude of d
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param blk the block configuration of the image in CUDA
/// \param gid the grid configuration of the image in CUDA
void GradientOfTV_SHARED(float* f, float* d, const float coef, cuint L, cuint W, const dim3& blk, const dim3& gid);
/// \brief The derivation of the TV norm, it is implemented in GPU, it is based on SHARED memory, and it theoretically more effective, more details can be get from "GPU-Based Acceleration for Interior Tomography"
/// \param f the pointer to the image
/// \param d the output image storing the TV derivation
/// \param coef the coefficient control the magnitude of d
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param blk the block configuration of the image in CUDA
/// \param gid the grid configuration of the image in CUDA
void GradientOfTV_SHARED(double* f, double* d, const double coef, cuint L, cuint W, const dim3& blk, const dim3& gid);


//void GradientOfTV_SHARED(float* f, float* d, const float coef, cuint L, cuint W, cuint H, const dim3& blk, const dim3& gid);
//void GradientOfTV_SHARED(double* f, double* d, const double coef, cuint L, cuint W, cuint H, const dim3& blk, const dim3& gid);

//void GradientOfTV_CPU(float* f, float* d, const float coef, cuint L, cuint W);
//void GradientOfTV_CPU(double* f, double* d, const double coef, cuint L, cuint W);
//void GradientOfTV_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H);
//void GradientOfTV_CPU(double* f, double* d, const double coef, cuint L, cuint W, cuint H);
//void GradientOfTV_GPU(float* f, float* d, const float coef, cuint& L ,cuint& W, const dim3& blk, const dim3& gid);
//void GradientOfTV_GPU(double* f, double* d, const double coef, cuint& L ,cuint& W, const dim3& blk, const dim3& gid);
//void GradientOfTV_GPU(float* f, float* d, const float coef, cuint& L ,cuint& W, cuint& H, const dim3& blk, const dim3& gid);
//void GradientOfTV_GPU(double* f, double* d, const double coef, cuint& L ,cuint& W, cuint& H, const dim3& blk, const dim3& gid);
//
//


/// \brief Discrete Gradient Transform of 2D image in GPU.
/// \param f image to be transformed
/// \param d the output image
/// \param coef the coefficient controlling the magnitude of d
/// \param L length resolution of image
/// \param W width resolution of image
/// \param blk block size configuration in CUDA
/// \param gid grid size configuration in CUDA
void DiscreteGradientTrans(float* f, float* d, const float coef, cuint L, cuint W, const dim3& blk, const dim3& gid); //Use this one;
/// \brief Discrete Gradient Transform of 2D image in GPU.
/// \param f image to be transformed
/// \param d the output image
/// \param coef the coefficient controlling the magnitude of d
/// \param L length resolution of image
/// \param W width resolution of image
/// \param blk block size configuration in CUDA
/// \param gid grid size configuration in CUDA
void DiscreteGradientTrans(double* f, double* d, const double coef, cuint L, cuint W, const dim3& blk, const dim3& gid);
/// \brief Discrete Gradient Transform of 3D image in GPU.
/// \param f image to be transformed
/// \param d the output image
/// \param coef the coefficient controlling the magnitude of d
/// \param L length resolution of image
/// \param W width resolution of image
/// \param H height resolution of image
/// \param blk block size configuration in CUDA
/// \param gid grid size configuration in CUDA
void DiscreteGradientTrans(float* f, float* d, const float coef, cuint L, cuint W, cuint H, const dim3& blk, const dim3& gid); //Use this one;
/// \brief Discrete Gradient Transform of 3D image in GPU.
/// \param f image to be transformed
/// \param d the output image
/// \param coef the coefficient controlling the magnitude of d
/// \param L length resolution of image
/// \param W width resolution of image
/// \param H height resolution of image
/// \param blk block size configuration in CUDA
/// \param gid grid size configuration in CUDA
void DiscreteGradientTrans(double* f, double* d, const double coef, cuint L, cuint W, cuint H, const dim3& blk, const dim3& gid);
/// \brief Discrete Gradient Transform of 2D image in GPU, it uses SHARED memory, more details are explained in "GPU-Based Acceleration for Interior Tomography"
/// \param f image to be transformed
/// \param d the output image
/// \param coef the coefficient controlling the magnitude of d
/// \param L length resolution of image
/// \param W width resolution of image
void DiscreteGradientTrans_SHARED(float* f, float* d, const float coef, cuint L, cuint W);
/// \brief Discrete Gradient Transform of 2D image in GPU, it uses SHARED memory, more details are explained in "GPU-Based Acceleration for Interior Tomography"
/// \param f image to be transformed
/// \param d the output image
/// \param coef the coefficient controlling the magnitude of d
/// \param L length resolution of image
/// \param W width resolution of image
void DiscreteGradientTrans_SHARED(double* f, double* d, const double coef, cuint L, cuint W);
/// \brief Discrete Gradient Transform of 2D image in GPU, it uses SHARED memory, more details are explained in "GPU-Based Acceleration for Interior Tomography"
/// \param f image to be transformed
/// \param d the output image
/// \param coef the coefficient controlling the magnitude of d
/// \param L length resolution of image
/// \param W width resolution of image
/// \param H height resolution of image
void DiscreteGradientTrans_SHARED(float* f, float* d, const float coef, cuint L, cuint W, cuint H);
/// \brief Discrete Gradient Transform of 2D image in GPU, it uses SHARED memory, more details are explained in "GPU-Based Acceleration for Interior Tomography"
/// \param f image to be transformed
/// \param d the output image
/// \param coef the coefficient controlling the magnitude of d
/// \param L length resolution of image
/// \param W width resolution of image
/// \param H height resolution of image
void DiscreteGradientTrans_SHARED(double* f, double* d, const double coef, cuint L, cuint W, cuint H);

//void DiscreteGradientTrans_CPU(float* f, float* d, const float coef, cuint L, cuint W);
//void DiscreteGradientTrans_CPU(double* f, double* d, const double coef, cuint L, cuint W);
//void DiscreteGradientTrans_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H);
//void DiscreteGradientTrans_CPU(double* f, double* d, const double coef, cuint L, cuint W, cuint H);
//void DiscreteGradientTrans_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid);
//void DiscreteGradientTrans_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid);
//void DiscreteGradientTrans_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid);
//void DiscreteGradientTrans_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid);
////ÒÑÓÐglobalº¯Êý£¬»¹Ã»ÓÐÊµÏÖ
//void DiscreteGradientTrans_Shared_GPU(float* f, float* d, const float coef, cuint& L, cuint& W);
//void DiscreteGradientTrans_Shared_GPU(double* f, double* d, const double coef, cuint& L, cuint& W);
//void DiscreteGradientTrans_Shared_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H);
//void DiscreteGradientTrans_Shared_GPU(double* f, double* d, const double coef, cuint& L, cuint& W, cuint& H);


/// \brief the pseudo inverse transform based on the soft thresholding method in CT reconstruction. 
/// \param f image to be transformed
/// \param d the TV transform of f
/// \param r the output result
/// \param omega the prior knowledge of the image that is TV norm of the ideal image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param blk block Size of the CUDA configuration
/// \param gid grid size of the CUDA configuration
void invDiscreteGradientTransform(float* f, float* d, float* r, const float omega, cuint L, cuint W, const dim3& blk, const dim3& gid); //Use this one;
/// \brief the pseudo inverse transform based on the soft thresholding method in CT reconstruction. 
/// \param f image to be transformed
/// \param d the TV transform of f
/// \param r the output result
/// \param omega the prior knowledge of the image that is TV norm of the ideal image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param blk block Size of the CUDA configuration
/// \param gid grid size of the CUDA configuration
void invDiscreteGradientTransform(double* f, double* d, double* r, const double omega, cuint L, cuint W, const dim3& blk, const dim3& gid);
/// \brief the pseudo inverse transform based on the soft thresholding method in CT reconstruction. 
/// \param f image to be transformed
/// \param d the TV transform of f
/// \param r the output result
/// \param omega the prior knowledge of the image that is TV norm of the ideal image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
/// \param blk block Size of the CUDA configuration
/// \param gid grid size of the CUDA configuration
void invDiscreteGradientTransform(float* f, float* d, float* r, const float omega, cuint L, cuint W, cuint H, const dim3& blk, const dim3& gid);
/// \brief the pseudo inverse transform based on the soft thresholding method in CT reconstruction. 
/// \param f image to be transformed
/// \param d the TV transform of f
/// \param r the output result
/// \param omega the prior knowledge of the image that is TV norm of the ideal image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
/// \param blk block Size of the CUDA configuration
/// \param gid grid size of the CUDA configuration
void invDiscreteGradientTransform(double* f, double* d, double* r, const double omega, cuint L, cuint W, cuint H, const dim3& blk, const dim3& gid);





//
//void invDiscreteGradientTransform_GPU(float* f, float* d, float* r, const float omega, cuint L, cuint W, const dim3& blk, const dim3& gid);
//void invDiscreteGradientTransform_GPU(double* f, double* d, double* r, const double omega, cuint L, cuint W, const dim3& blk, const dim3& gid);
//void invDiscreteGradientTransform_GPU(float* f, float* d, float* r, const float omega, cuint L, cuint W, cuint H, const dim3& blk, const dim3& gid);
//void invDiscreteGradientTransform_GPU(double* f, double* d, double* r, const double omega, cuint L, cuint W, cuint H,const dim3& blk, const dim3& gid);
//void invDiscreteGradientTransform_CPU(float* f, float* d, float* r, const float omega, cuint L, cuint W);
//void invDiscreteGradientTransform_CPU(double* f, double* d, double* r, const double omega, cuint L, cuint W);
//void invDiscreteGradientTransform_CPU(float* f, float* d, float* r, const float omega, cuint L, cuint W, cuint H);
//void invDiscreteGradientTransform_CPU(double* f, double* d, double* r, const double omega, cuint L, cuint W, cuint H);
//


/// \brief Find the optimum thresholding value with the prior knowledge of the image
/// \param TVImg the DGT transformed image
/// \param ObjTV the prior knowledge of the original image or an estimation
/// \param length the length of TVimg
float OptimizedW0(thrust::device_ptr<float>& TVImg, const float& ObjTV, const unsigned int length);
/// \brief Find the optimum thresholding value with the prior knowledge of the image
/// \param TVImg the DGT transformed image
/// \param ObjTV the prior knowledge of the original image or an estimation
/// \param length the length of TVimg
double OptimizedW0(thrust::device_ptr<double>& TVImg, const double& ObjTV, const unsigned int length);
/// \brief Find the optimum thresholding value with the prior knowledge of the image
/// \param TVImg the DGT transformed image
/// \param ObjTV the prior knowledge of the original image or an estimation
float OptimizedW0(thrust::device_vector<float>& TVImg, const float& ObjTV);
/// \brief Find the optimum thresholding value with the prior knowledge of the image
/// \param TVImg the DGT transformed image
/// \param ObjTV the prior knowledge of the original image or an estimation
double OptimizedW0(thrust::device_vector<double>& TVImg, const double& ObjTV);

//
///// \brief OS-SART reconstruction algorithm
///// \param himg the image in host device to be reconstructed
///// \param hprj the raw projection
///// \param himgWeg the back projection weighting
///// \param hmask the image mask usually a circle or the all one matrix
///// \param FanGeo the Equal angular detector
///// \param Img the Image configuration
///// \param iterNum the maximum iteration
///// \param lambda the update coefficient
///// \param subSetNum subset number for the reconstruction
//void OS_SART(float* himg, float* hprj, float* himgWeg, float* hmask, const FanEAGeo& FanGeo, const Image& Img,
//	cuint& iterNum, const float& lambda, const unsigned int subSetNum);
//
///// \brief OS-SART reconstruction algorithm
///// \param himg the image in host device to be reconstructed
///// \param hprj the raw projection
///// \param himgWeg the back projection weighting
///// \param hmsk the image mask usually a circle or the all one matrix
///// \param FanGeo the Equal angular detector
///// \param Img the Image configuration
///// \param iterNum the maximum iteration
///// \param lambda the update coefficient
//void SART(thrust::host_vector<float>& himg,
//	thrust::host_vector<float>& hprj,
//	thrust::host_vector<float>& himgWeg,
//	thrust::host_vector<float>& hmsk,
//	const FanEAGeo& FanGeo, const Image& Img,
//	cuint iterNum, const float lambda);
//
//
//
///// \brief OS-SART reconstruction algorithm with Steepest Descent of TV for CS
///// \param himg the image in host device to be reconstructed
///// \param hprj the projection data in host device for reconstruction
///// \param himgWeg the image weighting with all one back-projections
///// \param hmask the mask of the image for reconstruction, it is usually a circle or a all one matrix
///// \param FanGeo the equal angular detector
///// \param Img the Image configuration
///// \param iterNum the maximum iteration number
///// \param lambda the update coefficient
///// \param subSetNum how many subsets are used for reconstruction
///// \param initAlpha the initial alpha value for descent coefficient
///// \param apa_s each time we need to decrease initAlpha = initAlpha * apa_s
///// \param descentTime how many times for SD of TV
//void OS_SART_SD(float* himg, float* hprj, float* himgWeg, float* hmask, const FanEAGeo& FanGeo, const Image& Img,
//	cuint& iterNum, const float& lambda, const unsigned int subSetNum, float initAlpha, float apa_s, cuint descentTime);
//
//
///// \brief OS-SART reconstruction algorithm with Steepest Descent of TV for CS
///// \param himg the image in host device to be reconstructed
///// \param hprj the projection data in host device for reconstruction
///// \param hweg the image weighting with all one back-projections
///// \param hmsk the mask of the image for reconstruction, it is usually a circle or a all one matrix
///// \param FanGeo the equal angular detector
///// \param Img the Image configuration
///// \param iterNum the maximum iteration number
///// \param subSetNum how many subsets are used for reconstruction
///// \param updateCoef the update coefficient
///// \param initAlpha the initial alpha value for descent coefficient
///// \param apa_s each time we need to decrease initAlpha = initAlpha * apa_s
///// \param descentTime how many times for SD of TV
//void OS_SART_SD(thrust::host_vector<float>& himg,
//	thrust::host_vector<float>& hprj,
//	thrust::host_vector<float>& hweg,
//	thrust::host_vector<float>& hmsk,
//	const FanEAGeo& FanGeo,
//	const Image& Img,
//	cuint& iterNum,
//	cuint subSetNum,
//	const float updateCoef,
//	const float initAlpha,
//	const float apa_s,
//	cuint descentTime);
//
//
//
///// \brief OS-SART reconstruction algorithm with Soft Thresholding Filtering of TV for CS
///// \param himg the image in host device to be reconstructed
///// \param hprj the projection data in host device for reconstruction
///// \param himgWeg the image weighting with all one back-projections
///// \param hmask the mask of the image for reconstruction, it is usually a circle or a all one matrix
///// \param FanGeo the equal angular detector
///// \param Img the Image configuration
///// \param iterNum the maximum iteration number
///// \param lambda the update coefficient
///// \param subSetNum how many subsets are used in reconstruction
///// \param objTV the prior knowledge of the image
//void OS_SART_STF(float* himg, float* hprj, float* himgWeg, float* hmask, const FanEAGeo& FanGeo, const Image& Img,
//	cuint& iterNum, const float& lambda, const unsigned int subSetNum, const float objTV);
//
///// \brief OS-SART reconstruction algorithm with Soft Thresholding Filtering of TV for CS
///// \param himg the image in host device to be reconstructed
///// \param hprj the projection data in host device for reconstruction
///// \param hweg the image weighting with all one back-projections
///// \param hmsk the mask of the image for reconstruction, it is usually a circle or a all one matrix
///// \param FanGeo the equal angular detector
///// \param Img the Image configuration
///// \param iterNum the maximum iteration number
///// \param subSetNum how many subsets are used in reconstruction
///// \param objTV the prior knowledge of the image
///// \param updateCoef the update coefficient
//void OS_SART_STF(thrust::host_vector<float>& himg,
//	thrust::host_vector<float>& hprj,
//	thrust::host_vector<float>& hweg,
//	thrust::host_vector<float>& hmsk,
//	const FanEAGeo& FanGeo, const Image& Img,
//	cuint& iterNum, const unsigned int subSetNum,
//	const float objTV, const float updateCoef);
//
//
///// \brief The initial version of EM algorithm. 
///// \param hprj projection data in host memory for reconstruction
///// \param himg the image to be reconstructed
///// \param FanGeo Fan beam geometrical in equil angle detector
///// \param Img the image configuration
///// \param iterNum the maximum iteration number
//void EM(float* hprj, float* himg, const FanEAGeo& FanGeo, const Image& Img, cuint iterNum);
//
//
///// \brief The Conjugate Gradient method for reconstruction, but it seems that this algorithm doesn't coverge
///// \param himg the host image to be reconstructed
///// \param hprj the projection data for reconstruction
///// \param FanGeo the Fan beam geometry with equil angular detector
///// \param Img the image configuration
///// \param IterNum the maximum iteration number
//void CG_recon(float* himg, float* hprj, const FanEAGeo& FanGeo, const Image& Img, cuint IterNum);
//
///// \brief The conjugate gradient method for reconstruction, This algorithm is based on AIM model and executed in GPU
///// \param img Image to be reconstructed in main memory
///// \param prj Projection data in main memory
///// \param FanGeo FanEAGeo configuration
///// \param Img Image configuration
///// \param maxIter maximum iteration limit
///// \param err stop criteria
///// \param cases which beta updating method will be applied, details can be get from 
//void CG_recon_AIM(thrust::host_vector<float>& img, thrust::host_vector<float>& prj, const FanEAGeo& FanGeo, const Image& Img, cuint maxIter, float err, cuint cases);
///// \brief The conjugate gradient method for reconstruction, This algorithm is based on AIM model and executed in GPU
///// \param img Image to be reconstructed in main memory
///// \param prj Projection data in main memory
///// \param FanGeo FanEAGeo configuration
///// \param Img Image configuration
///// \param maxIter maximum iteration limit
///// \param err stop criteria
///// \param cases which beta updating method will be applied, details can be get from 
//void CG_recon_AIM(thrust::host_vector<double>& img, thrust::host_vector<double>& prj, const FanEAGeo& FanGeo, const Image& Img, cuint maxIter, double err, cuint cases);
//
//
///// \brief Generate the projection matrix with area integral model (AIM), it generates the model very slow, please be patient.
///// \param FanGeo FanEAGeo configuration
///// \param Img Image configuration
///// \param FileName The file name for storing the matrix
///// \param coeff non zero coefficients in the projection matrix
///// \param rowIdx non zero coefficients row Index in the matrix begin with 0
///// \param colIdx non zero coefficients column index in the matrix begin with 0
//void GenerateProjectionMatrixf_AIM(const FanEAGeo& FanGeo, const Image& Img, const std::string& FileName);
//
///// \brief Generate the projection matrix with area integral model (AIM), it generates the model very slow, please be patient.
///// \param FanGeo FanEAGeo configuration
///// \param Img Image configuration
///// \param FileName The file name for storing the matrix
///// \param coeff non zero coefficients in the projection matrix
///// \param rowIdx non zero coefficients row Index in the matrix begin with 0
///// \param colIdx non zero coefficients column index in the matrix begin with 0
//void GenerateProjectionMatrixd_AIM(const FanEAGeo& FanGeo, const Image& Img, const std::string& FileName);
//
///// \brief Generate the projection matrix with area integral model (AIM), it generates the model very slow, please be patient. This method use OpenMP to accelerate your program. The threads number is 30, so the projection number should be better the multiple of 30
///// \param FanGeo FanEAGeo configuration
///// \param Img Image configuration
///// \param FileName The file name for storing the matrix
//void GenerateMatrix_OpenMP_AIM_float(const FanEAGeo& FanGeo, const Image& Img, const std::string& FileName);
///// \brief Generate the projection matrix with area integral model (AIM), it generates the model very slow, please be patient. This method use OpenMP to accelerate your program. The threads number is 30, so the projection number should be better the multiple of 30
///// \param FanGeo FanEAGeo configuration
///// \param Img Image configuration
///// \param FileName The file name for storing the matrix
//void GenerateMatrix_OpenMP_AIM_double(const FanEAGeo& FanGeo, const Image& Img, const std::string& FileName);
//
///// \brief Generate the projection matrix with area integral model (AIM), it is a new version in single float format
///// \param FanGeo The FanEDGeo configuration;
///// \param Img The Image configuration;
//void genProjectionMatrix_AIMf(const FanEDGeo FanGeo, const Image Img);
///// \brief Generate the projection matrix with area integral model (AIM), it is a new version in double float format
///// \param FanGeo The FanEDGeo configuration;
///// \param Img The Image configuration;
//void genProjectionMatrix_AIMd(const FanEDGeo FanGeo, const Image Img);
//
///// \brief Generate the projection matrix with area integral model (AIM), it is a new version in single float format
///// \param FanGeo The FanEAGeo configuration;
///// \param Img The Image configuration;
//void genProjectionMatrix_AIMf(const FanEAGeo FanGeo, const Image Img);
///// \brief Generate the projection matrix with area integral model (AIM), it is a new version in double float format
///// \param FanGeo The FanEAGeo configuration;
///// \param Img The Image configuration;
//void genProjectionMatrix_AIMd(const FanEAGeo FanGeo, const Image Img);
//
//
//
//void genProjectionMatrix_Siddonf(const FanEDGeo FanGeo, const Image Img);
//void genProjectionMatrix_Siddond(const FanEDGeo FanGeo, const Image Img);
//
///// \brief This is the DEMO for Conjugate Gradient reconstruction, but it is wrong!
//void DEMO12();
//
///// \brief The Demo of OS-SART
//void DEMO1();
//
///// \brief The demo of OS-SART
///// \param FanGeo The equal angular detector based CT system geometry
///// \param Img the image configuration
///// \param FileName the projection file for reconstruction
///// \param WegName the all one backprojection file
///// \param MaskName the mask file for reconstruction
///// \param FouName The file that will store the reconstructed image
///// \param iterNum the maximum iterations
///// \param lambda the update coefficient 
///// \param subSetNum the subset number 
//void DEMO1_1(const FanEAGeo& FanGeo, const Image& Img,
//	const std::string& FileName,
//	const std::string& WegName,
//	const std::string& MaskName,
//	const std::string& FouName,
//	cuint iterNum,
//	const float lambda,
//	cuint subSetNum);
//
//
///// \brief This is the demo of OS-SART + SD
//void DEMO2();
//
///// \brief This is the demo of OS-SART + SD
///// \param FanGeo The CT configuration with equal angle detector system
///// \param Img the image configuration
///// \param FileName the projection data for reconstruction
///// \param WegName the weighting of all one back projection
///// \param MaskName the mask image
///// \param FouName the output image
///// \param iterNum the maximum iteration for reconstruction
///// \param lambda updating coefficient
///// \param subSetNum how many subsets are used
///// \param initAlpha the SD step size
///// \param apa_s initAlpha = initAlpha * apa_s for each decreasing step
///// \param descentTime for each iteration, how many descent steps are used
//void DEMO2_1(const FanEAGeo& FanGeo, const Image& Img,
//	const std::string& FileName,
//	const std::string& WegName,
//	const std::string& MaskName,
//	const std::string& FouName,
//	cuint iterNum,
//	const float lambda,
//	cuint subSetNum,
//	const float initAlpha,
//	const float apa_s,
//	cuint descentTime);
//
//
///// \brief the demo of OS-SART + STF
//void DEMO3();
//
//
///// \brief The demo of OS-SART + STF
///// \param FanGeo the fan beam configuration with equil angle detector
///// \param Img the Image configuration
///// \param FileName the projection file
///// \param WegName the all one backprojection for weighting
///// \param MaskName the mask with a circle or all one matrix
///// \param FouName The output image
///// \param iterNum the maximum iterations
///// \param lambda updating coefficient
///// \param subSetNum the subsets number for OS
///// \param objTV the prior knowledge to the image
//void DEMO3_1(const FanEAGeo& FanGeo, const Image& Img,
//	const std::string& FileName,
//	const std::string& WegName,
//	const std::string& MaskName,
//	const std::string& FouName,
//	cuint iterNum,
//	const float lambda,
//	cuint subSetNum,
//	const float objTV);
//
//


/// \brief test once projection and then back projection result
void DEMO4();

/// \brief Linear algebraic reconstruction from CSR and CSC file
void DEMO5();
/// \brief Linear algebraic reconstruction from COO file, this uses the unified memory
void DEMO5_1();

/// \brief Generate the projection from the given image
void DEMO6();

/// \brief Boxed based backprojection in one direction;
void DEMO7();
void DEMO7_1();

/// \brief The texture based projection process, but is has some problems
void DEMO8();


/// \brief The demo of SART
void DEMO9_1(); //OK

void DEMO9_2(); //NOT TESTED YET


//SART DEMO for multi slices reconstructed simultaneously
void DEMO10(); //OK

/// \brief OS-SART DEMO for multi slices reconstructed simultaneously with 3 GPUs, CUDA 6.5
void DEMO11();
/// \brief EM algorithm demo
void DEMO13();


/// \brief Demo for generating the projection matrix in AIM model;
void DEMO14();
///// \brief Demo for generating the projection matrix in AIM model, you can specify the file name and the configuration
///// \param FileName the file name you want to store the matrix in your disk
///// \param FanGeo FanEAGeo configuration
///// \param Img Image configuration
//void DEMO14f(const std::string& FileName, const FanEAGeo& FanGeo, const Image& Img);
///// \brief Demo for generating the projection matrix in AIM model, you can specify the file name and the configuration
///// \param FileName the file name you want to store the matrix in your disk
///// \param FanGeo FanEAGeo configuration
///// \param Img Image configuration
//void DEMO14d(const std::string& FileName, const FanEAGeo& FanGeo, const Image& Img);
//


/// \brief Demo for testing the time for several projection processes
void DEMO15();

/// \brief Multi slices simultaneously recosntructed in three GPUs
void DEMO16();

//
///// \brief Reconstruct the multi-features phantom with two GPUs, one is applied for storing the projection 
///// and the other is applied for reconstruction.
//void DEMO17();
///// \brief The generalized version
//void DEMO17(
//	const ConeEDGeo& ConeGeo,
//	const Volume& Vol,
//	const std::string& projectionFileName,
//	const std::string& weightingFileName,
//	const std::string& reconstFileName,
//	const std::vector<int>& subsetNumberSeries,
//	const int iterNum);
//
///// \brief Reconstruct the mouse species with 400 angles with two GPUs, one is applied for storing the projection
///// and the other is applied for reconstruction.
//void DEMO18();
//
//void DEMO18(
//	const ConeEDGeo& ConeGeo,
//	const Volume& Vol,
//	const std::string& projectionFileName,
//	const std::string& weightingFileName,
//	const std::string& reconstrctFileName,
//	const std::vector < int >& subSetNumSeries,
//	const int iterNum);

/// \brief The same as Reconstruction but with 2024x2024 sized projection data
void DEMO18v2();

/// \brief Reconstruction the mouse species with 400 angles with two GPUs, one is applied for storing the projection
/// and the other is applied for reconstruction. It is the older dataset based reconstruction.
/// ÓÐÎÊÌâ;
void DEMO18v3();

/// \brief Reconstruction the mouse species with 400 angles with two GPUs, one is applied for storing the projection
/// and the other is applied for reconstruction. It is the older dataset based reconstruction. We apply the SIR reconstruction here.
void DEMO18v4_2D();
/// ÎŽÍê³É;
void DEMO18v4_3D();

