/*
 * Wake Forest University Health Sciences.
 * Recast 8 uint8 number to two float numbers
 * author: Rui Liu
 * date: 2015.09.01
 * version: 1.0
 */

#include "mex.h"
#include "matrix.h"
#include <cstring>
#include <iostream>

union ReCast
{
	float Float;
	unsigned char Byte[4];
};


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{ 
  unsigned char* pt = (unsigned char*)mxGetPr(prhs[0]);
  ReCast Num1;
  Num1.Byte[0] = pt[0];
  Num1.Byte[1] = pt[1];
  Num1.Byte[2] = pt[2];
  Num1.Byte[3] = pt[3];

  ReCast Num2;
  Num2.Byte[0] = pt[4];
  Num2.Byte[1] = pt[5];
  Num2.Byte[2] = pt[6];
  Num2.Byte[3] = pt[7];

  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
  double* pp1 = mxGetPr(plhs[0]);
  double* pp2 = mxGetPr(plhs[1]);

  *pp1 = static_cast<double>(Num1.Float);
  *pp2 = static_cast<double>(Num2.Float);
}

