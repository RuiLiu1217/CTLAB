#include <mex.h>
#include <math.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <vector>
#include <algorithm>
#include <iostream>

#define TWOPI (6.283185307179586)
#define INV_TWOPI (0.1590250231624044)
#define TYPE double
#define DATATYPE DoubleReconData
#ifndef PI
#define PI 3.14159265358979
#endif
#define DeltaMax 0.0000001;

typedef struct {
    float* image;
    const float* proj;
    double* work_space;
} FloatReconData;

typedef struct {
    double* image;
    const double* proj;
    double* work_space;
} DoubleReconData;

static const double* geom[3];

inline TYPE bilinearInterpolation(TYPE YY, TYPE ZZ, int YLZL, int YL, int ProjIndex, const TYPE* Proj) {
    int YYD = int(YY);
    int YYU = YYD + 1;
    int ZZD = int(ZZ);
    int ZZU = ZZD + 1;

    TYPE alfa = YY - YYD;
    TYPE beta = ZZ - ZZD;
    return (Proj[ProjIndex * YLZL + ZZD * YL + YYD] * (1 - alfa) * (1 - beta) +
        Proj[ProjIndex * YLZL + ZZD * YL + YYU] * alfa * (1 - beta) +
        Proj[ProjIndex * YLZL + ZZU * YL + YYD] * (1 - alfa) * beta +
        Proj[ProjIndex * YLZL + ZZU * YL + YYU] * alfa * beta);
}


/*
 ============================================================================
 Name        : katsevich_backprojection.cu
 Author      : Rui Liu
 Version     :
 Copyright   : Your copyright notice
 Description : Compute sum of reciprocals using STL on CPU and Thrust on GPU
 ============================================================================
 */


 /**
 * This macro checks return value of the CUDA runtime call and exits
 * the application if the call failed.
 */
#if DEBUG
#define CUDA_CHECK_RETURN(value) {											\
	cudaError_t _m_cudaStat = value;										\
	if (_m_cudaStat != cudaSuccess) {										\
		fprintf(stderr, "Error %s at line %d in file %s\n",					\
				cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);		\
		exit(1);															\
	} }
 // Same function as CUDA_CHECK_RETURN
#define CUDA_SAFE_CALL(call) do{ cudaError_t err = call; if (cudaSuccess != err) {  fprintf (stderr, "Cuda error in file '%s' in line %i : %s.", __FILE__, __LINE__, cudaGetErrorString(err) );  exit(EXIT_FAILURE);  } } while (0)
#else
#define CUDA_CHECK_RETURN(value) {value;}
#define CUDA_SAFE_CALL(value) {value;}
#endif

#ifndef nullptr
#define nullptr NULL
#endif

#ifndef EPSILON
#define EPSILON (0.0000001)
#endif




 //Create texture object and corresponding cudaArray function
template<typename T>
void createTextureObject(
    cudaTextureObject_t& texObj, //return: texture object pointing to the cudaArray
    cudaArray* d_prjArray, // return: cudaArray storing the data
    int Width, int Height, int Depth, // data size
    T* sourceData, // where is the data
    cudaMemcpyKind memcpyKind, // data from host or memory
    cudaTextureAddressMode addressMode, // how to address the texture (clamp, border ...)
    cudaTextureFilterMode textureFilterMode, // usually linear filtering (double --> int2 use pointer not linear interpolation)
    cudaTextureReadMode textureReadMode, // usually use element wise reading mode.
    bool isNormalized) // usually false
{
    cudaExtent prjSize;
    prjSize.width = Width;
    prjSize.height = Height;
    prjSize.depth = Depth;
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<T>();
    cudaMalloc3DArray(&d_prjArray, &channelDesc, prjSize);
    cudaMemcpy3DParms copyParams = { 0 };
    copyParams.srcPtr = make_cudaPitchedPtr(
        (void*)sourceData, prjSize.width * sizeof(T),
        prjSize.width, prjSize.height);
    copyParams.dstArray = d_prjArray;
    copyParams.extent = prjSize;
    copyParams.kind = memcpyKind;
    cudaMemcpy3D(&copyParams);
    cudaResourceDesc resDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeArray;
    resDesc.res.array.array = d_prjArray;
    cudaTextureDesc texDesc;
    memset(&texDesc, 0, sizeof(texDesc));
    texDesc.addressMode[0] = addressMode;
    texDesc.addressMode[1] = addressMode;
    texDesc.addressMode[2] = addressMode;
    texDesc.filterMode = textureFilterMode;
    texDesc.readMode = textureReadMode;

    texDesc.normalizedCoords = isNormalized;
    CUDA_SAFE_CALL(cudaCreateTextureObject(&texObj, &resDesc, &texDesc, nullptr));
}


// Destroy a GPU array and corresponding TextureObject
void destroyTextureObject(cudaTextureObject_t& texObj, cudaArray* d_array)
{
    cudaDestroyTextureObject(texObj);
    cudaFreeArray(d_array);
}




template<typename T>
__device__ __host__ inline void PISegment(const T x, const T y, const T z, T& BAngle, T& TAngle)
{
    T delta = 1.0;
    T bmax = z * TWOPI;
    T bmin = bmax - TWOPI;

    T r2 = x * x + y * y;

    T sb = 0;
    T st = 0;

    while ((bmax - bmin > EPSILON) && (delta > EPSILON))
    {
        sb = (bmax + bmin) * 0.5;
        T sinsb = sin(sb);
        T cossb = cos(sb);
        T tempcos = 2.0 * (1.0 - y * sinsb - x * cossb);
        assert(tempcos != 0);
        T t = (1.0 - r2) / tempcos;
        T templan = acos((y * cossb - x * sinsb) / sqrt(tempcos + r2 - 1.0));
        st = 2.0 * templan + sb;
        T zz = (sb * t + (1.0 - t) * st) * INV_TWOPI;
        if (zz < z)
        {
            bmin = sb;
        }
        else
        {
            bmax = sb;
        }
        delta = fabs(zz - z);
    }

    BAngle = sb;
    TAngle = st;
}


// Note: this projection do not consider the edging situation for backprojection
__global__ void backProjectionKer(
    float* Image,
    cudaTextureObject_t ProjTex,
    int RecMX, int RecMY, int RecMZ,
    float ObjRSquare,
    float* __restrict__ xCor, float* __restrict__ yCor, float* __restrict__ zCor,
    float ScanR,
    float StdDis,
    float DeltaU, float DeltaV,
    float HalfY, float HalfZ,
    float HelicP,
    float DeltaL,
    float ProjCtr,
    int ProjBeginIndex,
    int ProjNum,
    int ProjScale)
{
    int Zindex = threadIdx.x + blockIdx.x * blockDim.x;
    int Xindex = threadIdx.y + blockIdx.y * blockDim.y;
    int Yindex = threadIdx.z + blockIdx.z * blockDim.z;
    const size_t idx = (Xindex * RecMY + Yindex) * RecMZ + Zindex;
    if (Xindex < RecMX && Yindex < RecMY && Zindex < RecMZ)
    {
        const float X = xCor[Xindex];
        const float Y = yCor[Yindex];
        const float Z = zCor[Zindex];
        if (std::pow(X, 2.0f) + std::pow(Y, 2.0f) >= ObjRSquare)
            return;

        double BAngle;
        double TAngle;
        PISegment<double>(X / ScanR, Y / ScanR, Z / HelicP, BAngle, TAngle);
        BAngle = BAngle / DeltaL + ProjCtr - ProjBeginIndex;
        TAngle = TAngle / DeltaL + ProjCtr - ProjBeginIndex;

        int Bindex = int(BAngle);
        int Tindex = ceil(TAngle);
        if (Bindex < 0)
        {
            Bindex = 0;
        }
        if (Bindex > ProjNum - 1)
        {
            Bindex = ProjNum - 1;
        }
        if (Tindex < 0)
        {
            Tindex = 0;
        }
        if (Tindex > ProjNum - 1)
        {
            Tindex = ProjNum - 1;
        }
        float tpdata = 0.0f;
        for (int ProjIndex = Bindex; ProjIndex <= Tindex; ProjIndex++)
        {
            float theta = (ProjIndex + ProjBeginIndex - ProjCtr) * DeltaL;
            float cost = cosf(theta);
            float sint = sinf(theta);

            float DPSx = X - ScanR * cost;
            float DPSy = Y - ScanR * sint;
            float DPSz = Z - HelicP * theta * INV_TWOPI;
            float factor = sqrtf(DPSx * DPSx + DPSy * DPSy + DPSz * DPSz);
            float fenmu = -(DPSx * cost + DPSy * sint);
            float YY = DPSy * cost - DPSx * sint;
            YY = YY * StdDis / (fenmu * DeltaU) + HalfY;
            float ZZ = DPSz * StdDis / (fenmu * DeltaV) + HalfZ;
            float temp = tex3D<float>(ProjTex, YY + 0.5f, ZZ + 0.5f, ProjIndex + 0.5f);
            tpdata += temp / factor;
        }
        tpdata = -tpdata / ProjScale;
        Image[idx] = tpdata;

    }
}



void backProjection(
    thrust::host_vector<float>& hImage,
    thrust::host_vector<float>& hProj,
    int RecMX, int RecMY, int RecMZ,
    float ObjRSquare,
    const thrust::host_vector<float>& hxCor,
    const thrust::host_vector<float>& hyCor,
    const thrust::host_vector<float>& hzCor,
    float ScanR,
    float StdDis,
    float DeltaU, float DeltaV,
    float HalfY, float HalfZ,
    int YL, int YLZL,
    float HelicP,
    float DeltaL,
    float ProjCtr,
    float ProjBeginIndex,
    int ProjNum, // number of projections
    int ProjScale,
    int threadidx, int threadidy, int threadidz)
{
    thrust::device_vector<float> Image = hImage;
    thrust::device_vector<float> xCor = hxCor;
    thrust::device_vector<float> yCor = hyCor;
    thrust::device_vector<float> zCor = hzCor;

    dim3 blk(threadidx, threadidy, threadidz);
    dim3 gid(
        (RecMZ + blk.x - 1) / blk.x,
        (RecMX + blk.y - 1) / blk.y,
        (RecMY + blk.z - 1) / blk.z);
    int ZL = YLZL / YL;
    cudaTextureObject_t projTex;
    cudaArray* d_projArray = nullptr;
    createTextureObject<float>(projTex, d_projArray,
        YL, ZL, ProjNum,
        &(hProj[0]),
        cudaMemcpyHostToDevice,
        cudaAddressModeClamp,
        cudaFilterModeLinear,
        cudaReadModeElementType, false);

    backProjectionKer << <gid, blk >> > (
        thrust::raw_pointer_cast(&Image[0]),
        projTex,
        RecMX, RecMY, RecMZ, ObjRSquare,
        thrust::raw_pointer_cast(&xCor[0]),
        thrust::raw_pointer_cast(&yCor[0]),
        thrust::raw_pointer_cast(&zCor[0]),
        ScanR, StdDis, DeltaU, DeltaV, HalfY, HalfZ,
        HelicP, DeltaL, ProjCtr, ProjBeginIndex,
        ProjNum, ProjScale);
    hImage = Image;
    destroyTextureObject(projTex, d_projArray);
}


void backProjection(
    double* hImage, const double* hProj,
    int RecMX, int RecMY, int RecMZ,
    double ObjRSquare,
    double* hxCor, double* hyCor, double* hzCor,
    double ScanR,
    double StdDis,
    double DeltaU, double DeltaV,
    double HalfY, double HalfZ,
    int YL, int YLZL,
    double HelicP, double DeltaL, double ProjCtr,
    int ProjBeginIndex,
    int ProjNum, // number of projections
    int ProjScale,
    int threadidx, int threadidy, int threadidz)
{
    thrust::host_vector<float> Image(RecMX * RecMY * RecMZ, 0.0f);
    thrust::host_vector<float> Proj(ProjNum * YLZL, 0.0f);
    thrust::host_vector<float> xCor(ProjNum, 0.0f);
    thrust::host_vector<float> yCor(ProjNum, 0.0f);
    thrust::host_vector<float> zCor(ProjNum, 0.0f);

    thrust::fill(Image.begin(), Image.end(), 0.0f);

    thrust::copy(hProj, hProj + ProjNum * YLZL, &(Proj[0]));
    thrust::copy(hxCor, hxCor + RecMX, &(xCor[0]));
    thrust::copy(hyCor, hyCor + RecMY, &(yCor[0]));
    thrust::copy(hzCor, hzCor + RecMZ, &(zCor[0]));

    backProjection(Image, Proj, RecMX, RecMY, RecMZ, ObjRSquare,
        xCor, yCor, zCor,
        ScanR, StdDis, DeltaU, DeltaV, HalfY, HalfZ, YL, YLZL,
        HelicP, DeltaL, ProjCtr, ProjBeginIndex, ProjNum,
        ProjScale, threadidx, threadidy, threadidz);
    thrust::copy(&(Image[0]), &(Image[0]) + RecMX * RecMY * RecMZ, hImage);
}


void PISegment(TYPE x, TYPE y, TYPE z, TYPE& BAngle, TYPE& TAngle);

void VoxelDrivenBkProj(DATATYPE* data) {
    /*Compute system parameter*/
    TYPE ScanR, ObjR, StdDis, HelicP, DecWidth, DecHeigh, ProjCtr;
    int  ProjScale, ProjBeginIndex, ProjEndIndex, ProjNum, YL, ZL, RecMX, RecMY, RecMZ;

    //ScanR  StdDis HelicP ProjScale ProjBeginIndex ProjEndIndex ProjCenter
    //DecWidth YL DecHeigh ZL
    //ObjR RecSize RecSize RecSize
    ScanR = geom[0][0];
    ObjR = geom[2][0];
    StdDis = geom[0][1];
    HelicP = geom[0][2];
    DecWidth = geom[1][0];
    DecHeigh = geom[1][2];
    ProjCtr = geom[0][6];
    ProjScale = int(geom[0][3]);
    ProjBeginIndex = int(geom[0][4]);
    ProjEndIndex = int(geom[0][5]);
    ProjNum = ProjEndIndex - ProjBeginIndex + 1;
    YL = int(geom[1][1]);
    ZL = int(geom[1][3]);
    RecMX = int(geom[2][1]);
    RecMY = int(geom[2][2]);
    RecMZ = int(geom[2][3]);

    printf("ProjBeginIndex=%d, ProjEndIndex=%d \n\t", ProjBeginIndex, ProjEndIndex);


    /*Compute some often used constant*/
    TYPE DeltaL, DeltaU, DeltaV, DeltaX, DeltaY, DeltaZ;
    TYPE HalfY, HalfZ, XCenter, YCenter, ZCenter;

    DeltaL = 2 * PI / ProjScale;
    DeltaU = DecWidth / YL;
    DeltaV = DecHeigh / ZL;
    HalfY = (YL - 1) * 0.5;
    HalfZ = (ZL - 1) * 0.5;
    DeltaX = 2 * ObjR / RecMX;
    DeltaY = 2 * ObjR / RecMY;
    DeltaZ = 2 * ObjR / RecMZ;
    XCenter = (RecMX - 1) * 0.5;
    YCenter = (RecMY - 1) * 0.5;
    ZCenter = (RecMZ - 1) * 0.5;

    /*Compute the coordinates of scanning locus and its corresponding local coordinate*/
    int  loop;
    TYPE temp;

    TYPE* VectorS = new TYPE[3 * ProjNum];
    TYPE* VectorE1 = new TYPE[2 * ProjNum];
    TYPE* VectorE2 = new TYPE[2 * ProjNum];
    TYPE* xCor = new TYPE[RecMX];
    TYPE* yCor = new TYPE[RecMY];
    TYPE* zCor = new TYPE[RecMZ];
    TYPE* Coef = new TYPE[ProjNum];

    for (loop = 0; loop < ProjNum; loop++) {
        temp = (loop + ProjBeginIndex - ProjCtr) * DeltaL;
        VectorS[loop * 3] = ScanR * cos(temp);
        VectorS[loop * 3 + 1] = ScanR * sin(temp);
        VectorS[loop * 3 + 2] = HelicP * temp / (2 * PI);
        VectorE1[loop * 2] = -cos(temp);
        VectorE1[loop * 2 + 1] = -sin(temp);
        VectorE2[loop * 2] = -sin(temp);
        VectorE2[loop * 2 + 1] = cos(temp);
    }

    for (loop = 0; loop < RecMX; loop++)
        xCor[loop] = (loop - XCenter) * DeltaX;
    for (loop = 0; loop < RecMY; loop++)
        yCor[loop] = (loop - YCenter) * DeltaY;
    for (loop = 0; loop < RecMZ; loop++)
        zCor[loop] = (loop - ZCenter) * DeltaZ;

    //BackProjection to reconstruct the object
    int  Xindex, Yindex, Zindex, Bindex, Tindex, ProjIndex, LoopBegin, LoopEnd;
    int  YYU, YYD, ZZU, ZZD, YLZL, RecMXY;
    TYPE ObjRSquare, x, y, z, BAngle, TAngle, tpdata;
    TYPE DPSx, DPSy, DPSz, factor, fenmu, YY, ZZ, alfa, beta;
    const TYPE* Proj;
    TYPE* Image;

    YLZL = YL * ZL;
    RecMXY = RecMX * RecMY;
    Image = data->image;
    Proj = data->proj;
    ObjRSquare = std::pow(ObjR, 2);

    bool useGPU = true;
    if (useGPU) {
        int threadidx = 128;
        int threadidy = 4;
        int threadidz = 1;
        backProjection(Image, Proj, RecMX, RecMY, RecMZ,
            ObjRSquare, xCor, yCor, zCor,
            ScanR, StdDis, DeltaU, DeltaV, HalfY, HalfZ,
            YL, YLZL, HelicP, DeltaL, ProjCtr, ProjBeginIndex,
            ProjNum, ProjScale, threadidx, threadidy, threadidz);


    }
    else {

        for (Xindex = 0; Xindex < RecMX; Xindex++) {
            for (Yindex = 0; Yindex < RecMY; Yindex++) {
                if (std::pow(xCor[Xindex], 2) + std::pow(yCor[Yindex], 2) < ObjRSquare) {
                    for (Zindex = 0; Zindex < RecMZ; Zindex++)
                        //+ pow(zCor[Zindex],2) < ObjRSquare)
                    {
                        x = xCor[Xindex] / ScanR;
                        y = yCor[Yindex] / ScanR;
                        z = zCor[Zindex] / HelicP;
                        PISegment(x, y, z, BAngle, TAngle);
                        BAngle = BAngle / DeltaL + ProjCtr - ProjBeginIndex;
                        TAngle = TAngle / DeltaL + ProjCtr - ProjBeginIndex;

                        if (!((BAngle >= ProjNum - 1) || (TAngle <= 0)))
                        {
                            Bindex = int(BAngle + 10) - 10;
                            Tindex = int(TAngle + 10) - 9;
                            for (loop = 0; loop < ProjNum; loop++)
                                Coef[loop] = 1;

                            if (Bindex < 0)
                            {
                                LoopBegin = 1;
                            }
                            else if (Bindex == 0)
                            {
                                temp = 1 - BAngle + Bindex;
                                Coef[1] = 0.5 + (1 - temp * 0.5) * temp;
                                LoopBegin = 1;
                            }
                            else if (Bindex == ProjNum - 2)
                            {
                                temp = 1 - BAngle + Bindex;
                                Coef[ProjNum - 2] = std::pow(temp, 2) * 0.5;
                                LoopBegin = ProjNum - 2;
                            }
                            else
                            {
                                temp = 1 - BAngle + Bindex;
                                Coef[Bindex] = std::pow(temp, 2) * 0.5;
                                Coef[Bindex + 1] = 0.5 + (1 - temp * 0.5) * temp;
                                LoopBegin = Bindex;
                            }
                            ///////////////////////////////////
                            if (Tindex > ProjNum - 1) {
                                LoopEnd = ProjNum - 2;
                            }
                            else if (Tindex == ProjNum - 1)
                            {
                                temp = 1 + TAngle - Tindex;
                                Coef[ProjNum - 2] = 0.5 + (1 - temp * 0.5) * temp;
                                LoopEnd = ProjNum - 2;
                            }
                            else if (Tindex == 1)
                            {
                                LoopEnd = 1;
                                temp = 1 + TAngle - Tindex;
                                Coef[1] = std::pow(temp, 2) * 0.5;
                            }
                            else
                            {
                                LoopEnd = Tindex;
                                temp = 1 + TAngle - Tindex;
                                Coef[Tindex - 1] = 0.5 + (1 - temp * 0.5) * temp;
                                Coef[Tindex] = std::pow(temp, 2) * 0.5;
                            }

                            tpdata = 0;
                            for (ProjIndex = LoopBegin; ProjIndex <= LoopEnd; ProjIndex++)
                            {
                                DPSx = xCor[Xindex] - VectorS[ProjIndex * 3];
                                DPSy = yCor[Yindex] - VectorS[ProjIndex * 3 + 1];
                                DPSz = zCor[Zindex] - VectorS[ProjIndex * 3 + 2];
                                factor = sqrt(DPSx * DPSx + DPSy * DPSy + DPSz * DPSz);
                                fenmu = DPSx * VectorE1[ProjIndex * 2] + DPSy * VectorE1[ProjIndex * 2 + 1];
                                YY = DPSx * VectorE2[ProjIndex * 2] + DPSy * VectorE2[ProjIndex * 2 + 1];
                                YY = YY * StdDis / (fenmu * DeltaU) + HalfY;
                                ZZ = DPSz * StdDis / (fenmu * DeltaV) + HalfZ;
                                YYD = int(YY);
                                YYU = YYD + 1;
                                ZZD = int(ZZ);
                                ZZU = ZZD + 1;
                                alfa = YY - YYD;
                                beta = ZZ - ZZD;
                                temp = Proj[ProjIndex * YLZL + ZZD * YL + YYD] * (1 - alfa) * (1 - beta) +
                                    Proj[ProjIndex * YLZL + ZZD * YL + YYU] * alfa * (1 - beta) +
                                    Proj[ProjIndex * YLZL + ZZU * YL + YYD] * (1 - alfa) * beta +
                                    Proj[ProjIndex * YLZL + ZZU * YL + YYU] * alfa * beta;
                                //YYD = int(YY+0.5);
                                //ZZD = int(ZZ+0.5);
                                //temp = Proj[ProjIndex*YLZL+ZZD*YL+YYD];
                                temp = temp * Coef[ProjIndex];
                                tpdata = tpdata + temp / factor;
                            }
                            tpdata = -tpdata / ProjScale;
                            ProjIndex = (Xindex * RecMY + Yindex) * RecMZ + Zindex;
                            Image[ProjIndex] = tpdata;
                        }//if((BAngle>0)&(TAngle<ProjNum-1))
                    }//if( pow(xCor[Xindex],2) + pow(yCor[Yindex],2)...
                }
            }
        }
    }


    delete[] VectorS;
    delete[] VectorE1;
    delete[] VectorE2;
    delete[] xCor;
    delete[] yCor;
    delete[] zCor;
    delete[] Coef;
}

void PISegment(TYPE x, TYPE y, TYPE z, TYPE& BAngle, TYPE& TAngle) {
    TYPE delta, dm, LanbudaInital, bmin, bmax, r2;
    TYPE sb, st, t, tempcos, templan, zz;

    delta = 1;
    dm = DeltaMax;
    LanbudaInital = z * 2 * PI;
    bmin = LanbudaInital - 2 * PI;
    bmax = LanbudaInital;
    r2 = x * x + y * y;

    while (((bmax - bmin) > dm) && (delta > dm)) {
        sb = (bmax + bmin) * 0.5;
        tempcos = 2 * (1 - y * sin(sb) - x * cos(sb));
        t = (1 - r2) / tempcos;
        templan = acos((y * cos(sb) - x * sin(sb)) / sqrt(tempcos + r2 - 1));
        st = 2 * templan + sb;
        zz = (sb * t + (1 - t) * st) / (2 * PI);
        if (zz < z)
            bmin = sb;
        else
            bmax = sb;
        delta = fabs(zz - z);
    }

    BAngle = sb;
    TAngle = st;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int dims[3];
    /* Check for proper number of arguments */
    if (nrhs != 3) {
        mexErrMsgTxt("bkproj requires two input arguments.");
    }
    if (nlhs != 0) {
        mexErrMsgTxt("bkproj requires one output argument.");
    }
    /* argument about scanning locus: ScanR, StdDis, HelicP, ProjScale, ProjNumber */
	geom[0] = mxGetPr(mxGetFieldByNumber(prhs[0], 0, 0));

    /* argument about detector: DecWidth, YL, DecHeigh, ZL */
    geom[1] = mxGetPr(mxGetFieldByNumber(prhs[0], 0, 1));

	/* argument about object:  ObjR RecMX RecMY RecMZ*/
	geom[2] = mxGetPr(mxGetFieldByNumber(prhs[0], 0, 2));
	
    if (mxGetClassID(prhs[1]) == mxSINGLE_CLASS) 
	{
        FloatReconData data;
        mexErrMsgTxt("Single precision is not supported in this version.\n");
    } 
	else 
	{
        DoubleReconData data;
        /* Assign pointers to the various parameters */
        data.image = mxGetPr(prhs[2]);
        data.proj  = mxGetPr(prhs[1]);
        /* Do the actual computations */
        VoxelDrivenBkProj(&data);
    }
}
