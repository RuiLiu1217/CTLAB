/*
 * Wake Forest Health Sciences & University of Massachusetts Lowell
 * Organization: 
 *  Wake Forest Health Sciences
 *
 * HelicalToFanFunc_mex.cpp
 * Matlab routine for Helical projection rebinning in CPU
 *
 * author: Rui Liu (Wake Forest Health Sciences)
 * date: 2016.04.21
 * version: 1.0
 */

#include "mex.h"
#include "matrix.h"
#include <cstring>
#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <new>
#include <exception>
#include <omp.h>

namespace AAPM{
	template<typename T>
	const T mod(const T& lambda, const T& regV)
	{
		T v = lambda;
		while(v > regV)
		{
			v -= regV;
		}
		return v;
	}
}

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

extern "C"
void HelicalToFan(
		float* Proj,  					// rebinned projection; in order (Channel, View Index, Slice Index)  TODO: need permute after call by MATLAB
		float* Views, 					// rebinned views (View Index, Slice Index) TODO: need permute after call by MATLAB
		float* proj,  					// raw projection data In order : (Height Index, Channel Index, View Index(Total View))
		float* zPos,  					// sampling position
		const int SLN,                  // slice number
		const float SD, 				// source to detector distance
		const float SO,					// source to iso-center distance
		const float BVAngle,        	// The begin view
		const int DetWidth,         	// number of detector columns
		const int DetHeight,        	// number of detector rows
		const float PerDetW,        	// Detector cell size along channel direction
		const float PerDetH,        	// Detector cell size along bench moving direction
		const int DefTimes,         	// Number of views per rotation
		const float DetCenterW,     	// Detector Center Index
		const float SpiralPitchFactor 	// Pitch defined in SIEMENS
		)
{
	// Calculate necessary parameters
	const float PerDetA = 2.0f * atanf(PerDetW / SD * 0.5f);
	const float delta = 2.0 * M_PI / DefTimes;
	const float DetCenterHreal = (DetHeight - 1.0) * 0.5 / PerDetH;
	const float h = SpiralPitchFactor * PerDetH * DetHeight * SO / SD;
	const float deltaZ = h / DefTimes;

    omp_set_num_threads(32);
#pragma omp parallel for
	for(int iii = 0; iii < SLN; iii++)
	{
		float z = zPos[iii];
        float* PData = new float[DefTimes * DetWidth];
        int* mark = new int[DefTimes * DetWidth];
        float* ViewAngles = new float[DefTimes];
        
		std::fill(PData, PData + DefTimes * DetWidth, 0.0f);
		std::fill(mark, mark + DefTimes * DetWidth, 0.0f);
		int start = floorf((z - h * 0.5f) / deltaZ);
		std::fill(ViewAngles, ViewAngles + DefTimes, 0.0f);
 		for(size_t i = 1; i <= start + DefTimes - 1; ++i)
 		{
 			float* data = proj + (i-1) * DetWidth * DetHeight; // In order (Height Idx, Channel Idx)

 			if(i >= start)
 			{
                float zz = (i - 1) * deltaZ;
 				float lambda = (i - 1) * delta + BVAngle;
 				float lambdaFan = AAPM::mod<float>(lambda, 2 * M_PI);
 				int ibeta = i - start + 1;
 				ViewAngles[ibeta - 1] = lambdaFan;
 
 				for(size_t j = 1; j <= DetWidth; ++j)
 				{
 					float ang = (j - DetCenterW) * PerDetA; //
 					float SM = SO * cosf(ang); //[Noo, 1999, Single-Slice rebinning method]
 					float tanAngV = (z - zz) / SM;
 					float v = SD * tanAngV; // curve detector
 					float dv = (DetCenterHreal + v) / PerDetH + 1.0;
 					int idv = floorf(dv);
 					float t = dv - idv;
 					float cosAngV = SD / sqrtf(SD * SD + v * v);
 					// --------------------- linear interpolation
 					float temp = 0;
 					int tempmark = 0;
 
 					if ((idv >= 1) && (idv <= DetHeight) &&
 						(idv+1 >= 1) && (idv + 1<=DetHeight))
                    {
 						temp = data[(j-1) * DetHeight + (idv - 1)] * (1 - t)
 							 + data[(j-1) * DetHeight + idv] * t; // problem here
 						tempmark = 1;
 					}
                    if(((ibeta - 1) * DetWidth + (j-1)) > DefTimes * DetWidth)
                    {
                        std::cout<<ibeta<<" "<<j<<"\n";
                    }
                    // In order (Channel Index, View Index);
 					PData[(ibeta - 1) * DetWidth + (j - 1)] = cosAngV * temp;
 					mark[(ibeta - 1) * DetWidth + (j - 1)] = tempmark;
 				} // end for j
			}// end if
 		}// end for i

 		int numOfab = 0;
 		for(size_t i = 1; i <= DefTimes; ++i)
 		{
 			for(size_t j = 1; j <= DetWidth; ++j)
 			{
  				if (mark[(i-1) * DetWidth + (j-1)] == 0)
  				{
  					numOfab = numOfab + 1;
  					std::fill(mark + (i-1) * DetWidth, mark + i * DetWidth, 0);
  					std::fill(PData + (i-1) * DetWidth, PData + i * DetWidth, 0);
  					break;
  				}
 			}
 		}
        
 		float ViewAngle = 360 * (DefTimes - numOfab) / DefTimes;
 		std::cout<<"The total rotation angle (degree) = " << iii << std::endl;
 
 
 		for(size_t angIdx = 0; angIdx != DefTimes; ++angIdx)
 		{
 			for(size_t detIdx = 0; detIdx != DetWidth; ++detIdx)
 			{
 				Proj[(iii * DefTimes + angIdx) * DetWidth + detIdx] =
 						PData[angIdx * DetWidth + detIdx];
 			}
 			Views[iii * DefTimes + angIdx] = ViewAngles[angIdx];
 		}
        delete[] PData;
        delete[] mark;
        delete[] ViewAngles;
	}
}




void HelicalToFan_MATLAB(
   mxArray* plhs[], // [view det-row det-col]
   const mxArray* mx_proj,
   const mxArray* mx_zPos,
   const mxArray* mx_SLN,
   const mxArray* mx_SD,
   const mxArray* mx_SO,
   const mxArray* mx_BVAngle,
   const mxArray* mx_DetWidth,
   const mxArray* mx_DetHeight,
   const mxArray* mx_PerDetW,
   const mxArray* mx_PerDetH,
   const mxArray* mx_DefTimes,
   const mxArray* mx_DetCenterW,
   const mxArray* mx_SpiralPitchFactor)
{
    const int SLN = *((int*)mxGetData(mx_SLN));
    const float SD = *((float*)mxGetData(mx_SD));
    const float SO = *((float*)mxGetData(mx_SO));
    const float BVAngle = *((float*)mxGetData(mx_BVAngle));
    const int DetWidth = *((int*)mxGetData(mx_DetWidth));
    const int DetHeight = *((int*)mxGetData(mx_DetHeight));
    const float PerDetW = *((float*)mxGetData(mx_PerDetW));
    const float PerDetH = *((float*)mxGetData(mx_PerDetH));
    const int DefTimes = *((int*)mxGetData(mx_DefTimes));
    const float DetCenterW = *((float*)mxGetData(mx_DetCenterW));
    const float SpiralPitchFactor = *((float*)mxGetData(mx_SpiralPitchFactor));
    
    //Create output array of class mxREAL and mxSINGLE_CLASS
    const mwSize dims_proj[] = {DetWidth, DefTimes, SLN};
    plhs[0] = mxCreateNumericArray(3,dims_proj,mxSINGLE_CLASS,mxREAL);
    const mwSize dims_view[] = {DefTimes, SLN};
    plhs[1] = mxCreateNumericArray(2,dims_view,mxSINGLE_CLASS,mxREAL);
    
    float* Proj = (float*) mxGetPr(plhs[0]);
    float* Views = (float*) mxGetPr(plhs[1]);
    float* proj = (float*) mxGetPr(mx_proj);
    float* zPos = (float*) mxGetPr(mx_zPos);
    
    HelicalToFan(Proj, Views, proj, zPos, SLN, SD, SO, BVAngle,  
		DetWidth, DetHeight, PerDetW, PerDetH, DefTimes, DetCenterW, SpiralPitchFactor);
    
}


void HelicalToFan_mex(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  if(nrhs < 1)
  {
    std::cerr<<"Error: usage\n";
    exit(-1);
  }
  HelicalToFan_MATLAB(plhs,prhs[0],prhs[1],prhs[2],prhs[3],prhs[4],prhs[5],prhs[6],prhs[7],prhs[8],prhs[9],prhs[10],prhs[11],prhs[12]);

}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  HelicalToFan_mex(nlhs, plhs, nrhs, prhs);
}


