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
#include <string>
#include <cstdlib>

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{ 
  char* pt = (char*)mxGetPr(prhs[0]);

  double a = atof(pt);
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  double* pp = mxGetPr(plhs[0]);
  *pp = a;


}

