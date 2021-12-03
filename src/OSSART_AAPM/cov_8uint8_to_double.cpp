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
	double Double;
	unsigned char Byte[8];
};


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{ 
  unsigned char* pt = (unsigned char*)mxGetPr(prhs[0]);
  ReCast Num1;
  Num1.Byte[0] = pt[0];
  Num1.Byte[1] = pt[1];
  Num1.Byte[2] = pt[2];
  Num1.Byte[3] = pt[3];
  Num1.Byte[4] = pt[4];
  Num1.Byte[5] = pt[5];
  Num1.Byte[6] = pt[6];
  Num1.Byte[7] = pt[7];
  
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  double* pp1 = mxGetPr(plhs[0]);
  
  *pp1 = static_cast<double>(Num1.Double);
}

