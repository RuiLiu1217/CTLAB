/*
**  CT/Micro CT lab
**  Department of Radiology
**  University of Iowa
**  Version of 2004.03.13
*/

#include <mex.h>
#include <math.h>
#define DOUBLE_PRECISION
#define TYPE double
#define PI 3.14159265358979
static const double PhPar[10][8] = {
{0.6900,  0.920,  0.900,   0,     0,      0,     0,     2.0},
{0.6624,  0.874,  0.880,   0,     0,      0,     0,   -0.98},
//{0.4100,  0.500,  0.500,   0,     0,      0,   108,      0.02}};
{0.4100,  0.160,  0.210,  -0.22,  0,    -0.25,  108,  -0.02},
{0.3100,  0.110,  0.220,   0.22,  0,    -0.25,   72,  -0.02},
{0.2100,  0.250,  0.500,   0,     0.35, -0.25,    0,   0.02},
{0.0460,  0.046,  0.046,   0,     0.1,  -0.25,    0,   0.02},
{0.0460,  0.023,  0.020,  -0.08, -0.65, -0.25,    0,   0.01},
{0.0460,  0.023,  0.020,   0.06, -0.65, -0.25,    90,  0.01},
{0.0560,  0.040,  0.100,   0.06, -0.105, 0.625,   90,  0.02},
{0.0560,  0.056,  0.100,   0,     0.100, 0.625,    0, -0.02}};

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    TYPE theta, pou, heli, DeltaU, DeltaV;
	int  YL, ZL;
	TYPE *ProjSlice;

    /* Check for proper number of arguments */
    if (nrhs != 7) {
        mexErrMsgTxt("cone beam projection need five inputs.");
    }
    if (nlhs != 1) {
        mexErrMsgTxt("cone beam projection need one output.");
    }
    /* The input parameters*/
	theta = mxGetScalar(prhs[0]);
    pou   = mxGetScalar(prhs[1]);
	heli  = mxGetScalar(prhs[2]);
	YL    = int(mxGetScalar(prhs[3]));
	ZL    = int(mxGetScalar(prhs[4]));
	DeltaU= mxGetScalar(prhs[5])/YL;
	DeltaV= mxGetScalar(prhs[6])/ZL;
    //mexPrintf("theta=%10.4f,pou=%10.4f,heli=%10.4f,YL=%d,ZL=%d",theta,pou,heli,YL,ZL);

    if (mxGetClassID(prhs[1]) == mxSINGLE_CLASS) 
	{
        mexErrMsgTxt("Single precision is not supported in this version.\n");
		return;
	} 
	else 
	{
        plhs[0] = mxCreateNumericMatrix(YL, ZL, mxDOUBLE_CLASS, mxREAL);
        /* Assign pointers to the various parameters */
        ProjSlice = mxGetPr(plhs[0]);
		//The following is the actural projection part
        //computer some necessary constant
		TYPE z_heigh,sintheta,costheta,YCenter,ZCenter;

		z_heigh  = theta*heli/(2*PI);
		sintheta = sin(theta);
		costheta = cos(theta);
        YCenter  = (YL-1)*0.5;
		ZCenter  = (ZL-1)*0.5;

        int Yindex,Zindex,ElpIndex;
		TYPE PS[3],SPS[3],XS[3],SXS[3],DX[3],tpdata,tempvar,AA,BB,CC;
		TYPE pcos[10],psin[10],AxisSquare[10][3];
		
		XS[0] = pou*costheta;
		XS[1] = pou*sintheta;
		XS[2] = z_heigh;
		for (ElpIndex=0;ElpIndex<10;ElpIndex++)
		{
			tpdata = PhPar[ElpIndex][6]*PI/180;
			pcos[ElpIndex] = cos(tpdata);
			psin[ElpIndex] = sin(tpdata);
			AxisSquare[ElpIndex][0] = pow(PhPar[ElpIndex][0],2);
			AxisSquare[ElpIndex][1] = pow(PhPar[ElpIndex][1],2);
			AxisSquare[ElpIndex][2] = pow(PhPar[ElpIndex][2],2);
		}

		for (Yindex=0;Yindex<YL;Yindex++)
		{
			PS[1] = (Yindex-YCenter)*DeltaU;
			PS[0] =-PS[1]*sintheta;
			PS[1] = PS[1]*costheta;

			for (Zindex=0;Zindex<ZL;Zindex++)
			{
                PS[2] = (Zindex-ZCenter)*DeltaV+z_heigh;
                tpdata = 0;

				for (ElpIndex=0;ElpIndex<10;ElpIndex++)
				{
					SXS[0] = XS[0]-PhPar[ElpIndex][3];
					SPS[0] = PS[0]-PhPar[ElpIndex][3];
					SXS[1] = XS[1]-PhPar[ElpIndex][4];
					SPS[1] = PS[1]-PhPar[ElpIndex][4];
					SXS[2] = XS[2]-PhPar[ElpIndex][5];
					SPS[2] = PS[2]-PhPar[ElpIndex][5];
                    tempvar= SXS[0];
					SXS[0] = tempvar*pcos[ElpIndex]+SXS[1]*psin[ElpIndex];
					SXS[1] =-tempvar*psin[ElpIndex]+SXS[1]*pcos[ElpIndex];
					tempvar= SPS[0];
					SPS[0] = tempvar*pcos[ElpIndex]+SPS[1]*psin[ElpIndex];
					SPS[1] =-tempvar*psin[ElpIndex]+SPS[1]*pcos[ElpIndex];

					DX[0]  = SXS[0]-SPS[0];
					DX[1]  = SXS[1]-SPS[1];
					DX[2]  = SXS[2]-SPS[2];
					tempvar= sqrt(pow(DX[0],2)+pow(DX[1],2)+pow(DX[2],2));
					DX[0]  = DX[0]/tempvar;
					DX[1]  = DX[1]/tempvar;
					DX[2]  = DX[2]/tempvar;
					
                    AA = pow(DX[0],2)/AxisSquare[ElpIndex][0]+
						 pow(DX[1],2)/AxisSquare[ElpIndex][1]+ 
						 pow(DX[2],2)/AxisSquare[ElpIndex][2];
					BB = DX[0]*SPS[0]/AxisSquare[ElpIndex][0]+
						 DX[1]*SPS[1]/AxisSquare[ElpIndex][1]+ 
						 DX[2]*SPS[2]/AxisSquare[ElpIndex][2];
                    CC = pow(SPS[0],2)/AxisSquare[ElpIndex][0]+ 
						 pow(SPS[1],2)/AxisSquare[ElpIndex][1]+ 
						 pow(SPS[2],2)/AxisSquare[ElpIndex][2]-1;
					tempvar=pow(BB,2)-AA*CC;
					if(tempvar>=0)
					{
						tpdata = tpdata+2*sqrt(tempvar)*PhPar[ElpIndex][7]/AA;
					}
				}//for (ElpIndex=0;ElpIndex<10;ElpIndex++)
				ProjSlice[Zindex*YL+Yindex]=tpdata;
			}//for (Zindex=0;Zindex<ZL;Zindex++)
		}//for (Yindex=0;Yindex<YL;Yindex++)        
		//////////////////////////////////////
    }//    if (mxGetClassID(prhs[1]) == mxSINGLE_CLASS) 
}
