/*!
 * \file FastMatrixVectorMultiplication.h
 * \brief Fast Matrix-Vector Multiplication in GPU
 * \author Rui Liu
 * \version 1.0
 * \date Jun. 1, 2015
 */

/// \brief Fast Matrix-Vector Multiplication;
/// \param y = Ax, on host memory
/// \param A matrix, on host memory
/// \param x vector to be multiplied, on host memory
/// \param m row number of the matrix
/// \param n col number of the matrix
void FMVM(float* y, float* A, float* x, int m, int n);

/// \brief Fast Matrix-Vector Multiplication;
/// \param d_y = Ax, on device memory
/// \param d_A matrix, on device memory
/// \param d_x vector to be multiplied, on device memory
/// \param m row number of the matrix
/// \param n col number of the matrix
void FMVM_device(float* d_y, float* dA, float* d_x, int m, int n);

/// \brief Test the Fast Matrix Vector Multiplication function
void testFMVM();
