/*!
* \file headers.hpp
* \brief This file mainly includes differential of the image in CPU and GPU
* \author Rui Liu
* \version 1.0
* \date 08/24/14
*/
#include "utilities.hpp"
#pragma once
#include <vector>
#include <thrust/device_vector.h>

/// \brief definition of const unsigned int
typedef const unsigned int cuint;


/// \brief X direction differential operation in CPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
void Dx_CPU(float* f, float* d, const float coef, cuint L, cuint W);
/// \brief X direction differential operation in CPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
void Dx_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H);
/// \brief second order of X direction differential operation in CPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
void D2x_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H);
/// \brief second order of X direction differential operation in CPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
void D2x_CPU(float* f, float* d, const float coef, cuint L, cuint W);
/// \brief The transport of X direction differential operation in CPU in matrix conception
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
void Dxt_CPU(float* f, float* d, const float coef, cuint L, cuint W);
/// \brief The transport of X direction differential operation in CPU in matrix conception
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
void Dxt_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H);
/// \brief The transport of X direction second order differential operation in CPU in matrix conception
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
void D2xt_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H);
/// \brief The transport of X direction second order differential operation in CPU in matrix conception
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
void D2xt_CPU(float* f, float* d, const float coef, cuint L, cuint W);
/// \brief Y direction differential operation in CPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
void Dy_CPU(float* f, float* d, const float coef, cuint L, cuint W);
/// \brief Y direction differential operation in CPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
void Dy_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H);
/// \brief Y direction second order differential operation in CPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
void D2y_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H);
/// \brief Y direction second order differential operation in CPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
void D2y_CPU(float* f, float* d, const float coef, cuint L, cuint W);
/// \brief The transport of Y direction differential operation in CPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
void Dyt_CPU(float* f, float* d, const float coef, cuint L, cuint W);
/// \brief The transport of Y direction differential operation in CPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
void Dyt_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H);
/// \brief The transport of second order of Y direction differential operation in CPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
void D2yt_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H);
/// \brief The transport of second order of Y direction differential operation in CPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
void D2yt_CPU(float* f, float* d, const float coef, cuint L, cuint W);
/// \brief Z direction differential operation in CPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
void Dz_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H);
/// \brief The second order of Z direction differential operation in CPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
void D2z_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H);
/// \brief The transport of Z direction differential operation in CPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
void Dzt_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H);
/// \brief The transport of second order Z direction differential operation in CPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
void D2zt_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H);
/// \brief The Laplacian transform of an image in CPU
/// \param f image
/// \param l result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
void Laplacian_CPU(float*f, float* l, const float coef, cuint L, cuint W);
/// \brief The Laplacian transform of an image in CPU
/// \param f image
/// \param l result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
void Laplacian_CPU(float* f, float* l, const float coef, cuint L, cuint W, cuint H);

/// \brief The Discrete Gradient Trans of an image in CPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
void DiscreteGradientTrans_CPU(float* f, float* d, const float coef, cuint L, cuint W);
/// \brief The Discrete Gradient Trans of an image in CPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
void DiscreteGradientTrans_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H);
/// \brief The differential of Discrete Gradient Trans of an image in CPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
void GradientOfTV_CPU(float* f, float* d, const float coef, cuint L, cuint W);
/// \brief The differential of Discrete Gradient Trans of an image in CPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
void GradientOfTV_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H);


/// \brief X direction differential operation in GPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void Dx_GPU(float* f, float* d, const float& coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid);
/// \brief second order of X direction differential operation in GPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void D2x_GPU(float* f, float* d, const float& coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid);
/// \brief X direction differential operation in GPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void Dx_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid);
/// \brief second order of X direction differential operation in GPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void D2x_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid);
/// \brief The transport of X direction differential operation in GPU in matrix conception
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void Dxt_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid);
/// \brief The transport of X direction second order differential operation in GPU in matrix conception
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void D2xt_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid);
/// \brief The transport of X direction differential operation in GPU in matrix conception
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void Dxt_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid);
/// \brief The transport of X direction second order differential operation in GPU in matrix conception
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void D2xt_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid);
/// \brief Y direction differential operation in GPU in matrix conception
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void Dy_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid);
/// \brief Y direction second order differential operation in GPU in matrix conception
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void D2y_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid);
/// \brief Y direction differential operation in GPU in matrix conception
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void Dy_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid);
/// \brief Y direction second order differential operation in GPU in matrix conception
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void D2y_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid);
/// \brief The transport of Y direction differential operation in GPU in matrix conception
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void Dyt_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid);
/// \brief The transport of Y direction differential operation in GPU in matrix conception
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void Dyt_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid);
/// \brief The transport of Y direction second order differential operation in GPU in matrix conception
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void D2yt_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid);
/// \brief The transport of Y direction second order differential operation in GPU in matrix conception
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void D2yt_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid);
/// \brief Z direction differential operation in GPU in matrix conception
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void Dz_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid);
/// \brief The transport of Z direction second order differential operation in GPU in matrix conception
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void D2z_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid);
/// \brief The transport of Z direction differential operation in GPU in matrix conception
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void Dzt_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid);
/// \brief The transport of Z direction second order differential operation in GPU in matrix conception
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void D2zt_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid);
/// \brief The Laplacian transform of image in GPU
/// \param f image
/// \param l result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void Laplacian_GPU(float* f, float* l, const float coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid);
/// \brief The Laplacian transform of image in GPU
/// \param f image
/// \param l result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void Laplacian_GPU(float* f, float* l, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid);
/// \brief The Discrete Gradient Transform of image in GPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void DiscreteGradientTrans_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid);
/// \brief The Discrete Gradient Transform of image in GPU
/// \param f image
/// \param d result
/// \param coef coefficient for updating the image
/// \param L length resolution of the image
/// \param W width resolution of the image
/// \param H height resolution of the image
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void DiscreteGradientTrans_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid);


/// \brief Share memory based discrete gradient transform
/// \param vol image to be transformed
/// \param tvvol discrete gradient transform result
/// \param coef the coefficient multiplying the tvvol
/// \param imgL image length resolution
/// \param imgW image width resolution
void DiscreteGradientTrans_Shared_GPU(float* vol, float* tvvol, float coef, cuint imgL, cuint imgW);
/// \brief Share memory based discrete gradient transform
/// \param vol image to be transformed
/// \param tvvol discrete gradient transform result
/// \param coef the coefficient multiplying the tvvol
/// \param imgL image length resolution
/// \param imgW image width resolution
void DiscreteGradientTrans_Shared_GPU(double* vol, double* tvvol, double coef, cuint imgL, cuint imgW);
/// \brief Share memory based discrete gradient transform
/// \param vol image to be transformed
/// \param tvvol discrete gradient transform result
/// \param coef the coefficient multiplying the tvvol
/// \param imgL image length resolution
/// \param imgW image width resolution
/// \param imgH image height resolution
void DiscreteGradientTrans_Shared_GPU(float* vol, float* tvvol, float coef, cuint imgL, cuint imgW, cuint imgH);
/// \brief Share memory based discrete gradient transform
/// \param vol image to be transformed
/// \param tvvol discrete gradient transform result
/// \param coef the coefficient multiplying the tvvol
/// \param imgL image length resolution
/// \param imgW image width resolution
/// \param imgH image height resolution
void DiscreteGradientTrans_Shared_GPU(double* vol, double* tvvol, double coef, cuint imgL, cuint imgW, cuint imgH);



/// \brief The gradient transform of TV
/// \param f image
/// \param d result
/// \param coef coefficient for image update
/// \param L image length resolution
/// \param W image width resolution
/// \param blk block size for CUDA
/// \param gid block size for CUDA
void GradientOfTV_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, const dim3& blk, const dim3& gid);
/// \brief The gradient transform of TV
/// \param f image
/// \param d result
/// \param coef coefficient for image update
/// \param L image length resolution
/// \param W image width resolution
/// \param H image height resolution
/// \param blk block size for CUDA
/// \param gid block size for CUDA
void GradientOfTV_GPU(float* f, float* d, const float coef, cuint& L, cuint& W, cuint& H, const dim3& blk, const dim3& gid);


/// \brief The inverse gradient transform of TV for STF in GPU
/// \param f image
/// \param d gradient image
/// \param r result
/// \param omega thresholding value
/// \param L image length resolution
/// \param W image width resolution
/// \param blk block size for CUDA
/// \param gid block size for CUDA
void invDiscreteGradientTransform_GPU(float* f, float* d, float* r, const float omega, cuint L, cuint W, const dim3& blk, const dim3& gid);
/// \brief The inverse gradient transform of TV for STF in GPU
/// \param f image
/// \param d gradient image
/// \param r result
/// \param omega thresholding value
/// \param L image length resolution
/// \param W image width resolution
/// \param H image height resolution
/// \param blk block size for CUDA
/// \param gid block size for CUDA
void invDiscreteGradientTransform_GPU(float* f, float* d, float* r, const float omega, cuint L, cuint W, cuint H, const dim3& blk, const dim3& gid);
/// \brief The inverse gradient transform of TV for STF in CPU
/// \param f image
/// \param d gradient image
/// \param r result
/// \param omega thresholding value
/// \param L image length resolution
/// \param W image width resolution
void invDiscreteGradientTransform_CPU(float* f, float* d, float* r, const float omega, cuint L, cuint W);
/// \brief The inverse gradient transform of TV for STF in CPU
/// \param f image
/// \param d gradient image
/// \param r result
/// \param omega thresholding value
/// \param L image length resolution
/// \param W image width resolution
/// \param H image height resolution
void invDiscreteGradientTransform_CPU(float* f, float* d, float* r, const float omega, cuint L, cuint W, cuint H);




/// \brief transform the raw data in float datatype to unsigned char type for volume rendering
/// \param inData input data
/// \param ouData output data
void rawToUnchar(std::vector<float>& inData, std::vector<unsigned char>& ouData);

/// \brief transform the raw data in float datatype to unsigned char type for volume rendering
/// \param inData input data
/// \param len length of the vector
/// \param ouData output data
void rawToUnchar(float* inData, const unsigned int len, std::vector<unsigned char>& ouData);



/// \brief Finding out the minimum and maximum value of a vector 
/// \param vec input vector
/// \param minV minimum value of the vector
/// \param maxV maximum value of the vector
void minmaxValue(thrust::device_vector<float>& vec, float& minV, float& maxV);
/// \brief Finding out the minimum and maximum value of a vector 
/// \param vec input vector
/// \param length length of the vector
/// \param minV minimum value of the vector
/// \param maxV maximum value of the vector
void minmaxValue(thrust::device_ptr<float>& vec, cuint length, float& minV, float& maxV);

/// \brief Generating the circle mask 
/// \param dmsk mask in device memory
/// \param ratio ratio controls the radius of the mask
/// \param L length of the image
/// \param W width of the image
/// \param blk block size in CUDA
/// \param gid grid size in CUDA
void gen01Mask(thrust::device_vector<float>& dmsk, const float& ratio, cuint& L, cuint& W, const dim3& blk, const dim3& gid);
/// \brief Generating the circle mask 
/// \param dmsk mask in device memory
/// \param ratio ratio controls the radius of the mask
/// \param L length of the image
/// \param W width of the image
/// \param blk block size in CUDA
/// \param gid grid size in CUDA
void gen01Mask(float* dmsk, const float& ratio, cuint& L, cuint& W, const dim3& blk, const dim3& gid);

/// \brief The first derivate of TV, get the descent direction with Shared memory in GPU
/// \param img image to be calculated
/// \param des descent direction
/// \param imgL image length resolution
/// \param imgW image width resolution
void descentDir_SHARED(float* img, float* des, cuint imgL, cuint imgW);
/// \brief The first derivate of TV, get the descent direction with Shared memory in GPU
/// \param img image to be calculated
/// \param des descent direction
/// \param imgL image length resolution
/// \param imgW image width resolution
void descentDir_SHARED(double* img, double* des, cuint imgL, cuint imgW);
/// \brief The first derivate of TV, get the descent direction with Shared memory in GPU
/// \param img image to be calculated
/// \param des descent direction
/// \param imgL image length resolution
/// \param imgW image width resolution
/// \param imgH image height resolution
void descentDir_SHARED(float* img, float* des, cuint imgL, cuint imgW, cuint imgH);
/// \brief The first derivate of TV, get the descent direction with Shared memory in GPU
/// \param img image to be calculated
/// \param des descent direction
/// \param imgL image length resolution
/// \param imgW image width resolution
/// \param imgH image height resolution
void descentDir_SHARED(double* img, double* des, cuint imgL, cuint imgW, cuint imgH);
