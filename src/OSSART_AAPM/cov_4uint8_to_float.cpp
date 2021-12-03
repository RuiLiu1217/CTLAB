/*
 * Wake Forest University Health Sciences.
 * Recast 4 uint8 number to a float number 
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
  ReCast Num;
  Num.Byte[0] = pt[0];
  Num.Byte[1] = pt[1];
  Num.Byte[2] = pt[2];
  Num.Byte[3] = pt[3];

  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  double* pp = mxGetPr(plhs[0]);
  *pp = static_cast<double>(Num.Float);

}

