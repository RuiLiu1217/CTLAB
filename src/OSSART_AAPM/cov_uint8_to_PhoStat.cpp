/*
 * Wake Forest University Health Sciences.
 * Recast 4 uint8 number to a serious of float numbers
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
  int* challenNum = (int*)mxGetPr(prhs[1]); 
   
  plhs[0] = mxCreateDoubleMatrix(*challenNum,1,mxREAL);
  double* pp = mxGetPr(plhs[0]);

  for(int i = 0; i != *challenNum; ++i)
  {
     ReCast Num;
     Num.Byte[0] = pt[i * 4 + 0];
     Num.Byte[1] = pt[i * 4 + 1]; 
     Num.Byte[2] = pt[i * 4 + 2];
     Num.Byte[3] = pt[i * 4 + 3];
     pp[i] = static_cast<double>(Num.Float);
  }

}

