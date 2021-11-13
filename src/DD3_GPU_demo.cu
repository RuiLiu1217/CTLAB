//#include "DD3_GPU_demo.h"
//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"
//
//#include <stdio.h>
//#include "utilities.cuh"
//#include "DD3_GPU_Proj.h"
//#include "DD3_GPU_Back.h"
//#include <sstream>
//#include <string>
//#include <thrust/device_vector.h>
//#include <thrust/inner_product.h>
//#include <thrust/transform.h>
//#include <ctime>
//#include <iostream>
//
//#include <thrust/tuple.h>
//#include <thrust/iterator/zip_iterator.h>
//
//
//template<typename T>
//struct CGop
//{
//	T alpha;
//	CGop(T a):alpha(a){}
//	T operator()(T x, T y)
//	{
//		return x + alpha * y;
//	}
//};
//
//template<typename T>
//struct prjWeight_functor
//{
//	typedef thrust::tuple<T,T,T> inputTuple;
//
//	__host__ __device__ T operator()(const inputTuple& input)
//	{
//		T prj = thrust::get<0>(input);
//		T realPrj = thrust::get<1>(input);
//		T rowSum = thrust::get<2>(input);
//		if(rowSum > 1.0E-7)
//		{
//			return (realPrj - prj) / rowSum;
//		}
//		else
//		{
//			return 0;
//		}
//	}
//};
//
////
//template<typename T>
//void prjWeight_v2(thrust::host_vector<T>& prj,thrust::host_vector<T>& realPrj,thrust::host_vector<T>& rowSum)
//{
//	thrust::transform(
//			thrust::make_zip_iterator(thrust::make_tuple(prj.begin(),realPrj.begin(),rowSum.begin())),
//			thrust::make_zip_iterator(thrust::make_tuple(prj.end(),realPrj.end(),rowSum.end())),
//			prj.begin(),prjWeight_functor<T>());
//}
//
//template<typename T>
//void prjWeight_v2(thrust::device_vector<T>& prj,thrust::device_vector<T>& realPrj,thrust::device_vector<T>& rowSum)
//{
//	thrust::transform(
//			thrust::make_zip_iterator(thrust::make_tuple(prj.begin(),realPrj.begin(),rowSum.begin())),
//			thrust::make_zip_iterator(thrust::make_tuple(prj.end(),realPrj.end(),rowSum.end())),
//			prj.begin(),prjWeight_functor<T>());
//}
//
//
//
//void prjWeight(float* prj, float* realPrj, float* rowSum, int N)
//{
//#pragma omp parallel for
//	for (int i = 0; i < N; ++i)
//	{
//		if (rowSum[i] > 1.0e-7)
//		{
//			prj[i] = (realPrj[i] - prj[i]) / rowSum[i];
//		}
//		else
//		{
//			prj[i] = 0;
//		}
//
//	}
//
//}
//
//template<typename T>
//struct bakWeight_functor
//{
//	typedef thrust::tuple<T,T,T> inputType;
//	__host__ __device__ T operator()(const inputType& input)
//	{
//		T vol = thrust::get<0>(input);
//		T reconImg = thrust::get<1>(input);
//		T colSum = thrust::get<2>(input);
//		if(colSum > 1.0E-7)
//		{
//			vol = vol / colSum;
//
//		}
//		else
//		{
//			vol = 0;
//		}
//		return reconImg + vol;
//	}
//};
//
//template<typename T>
//void bakWeight_v2(thrust::host_vector<T>& vol,thrust::host_vector<T>& reconImg,thrust::host_vector<T>& colSum)
//{
//	thrust::transform(
//			thrust::make_zip_iterator(thrust::make_tuple(vol.begin(),reconImg.begin(),colSum.begin())),
//			thrust::make_zip_iterator(thrust::make_tuple(vol.end(),reconImg.end(),colSum.end())),
//			reconImg.begin(),bakWeight_functor<T>());
//}
//template<typename T>
//void bakWeight_v2(thrust::device_vector<T>& vol,thrust::device_vector<T>& reconImg,thrust::device_vector<T>& colSum)
//{
//	thrust::transform(
//			thrust::make_zip_iterator(thrust::make_tuple(vol.begin(),reconImg.begin(),colSum.begin())),
//			thrust::make_zip_iterator(thrust::make_tuple(vol.end(),reconImg.end(),colSum.end())),
//			reconImg.begin(),bakWeight_functor<T>());
//}
//
//void bakWeight(float* vol, float* reconImg, float* colSum, int N)
//{
//#pragma omp parallel for
//	for (int i = 0; i < N; ++i)
//	{
//		if (colSum[i] > 1.0e-7)
//		{
//			vol[i] = vol[i] / colSum[i];
//		}
//		else
//		{
//			vol[i] = 0;
//		}
//		reconImg[i] = reconImg[i] + vol[i];
//	}
//
//}
//
////Test the correctness of the backprojection in curve detector.
//void testBackprojection()
//{
//
//	float sid = 541.0f;
//	float sdd = 949.0f;
//
//	float x0(0.0f);
//	float y0(sid);
//	float z0(0.0f);
//
//
//	int DNU(888);
//	int DNV(64);
//	int PN(1200);
//	float imgXCenter(0);
//	float imgYCenter(0);
//	float imgZCenter(0);
//
//	int XN(512);
//	int YN(512);
//	int ZN(64);
//
//
//	float dx(500.0f / 512.0f);
//	float dz(0.625);
//
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//	//Generate the positions of the detectors
//	float col_size = 1.0239;
//	float row_size = 1.0963;
//
//	float col_offset = 0;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//	imgXCenter = 0;
//	imgYCenter = 0;
//	imgZCenter = 0;
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(PN);
//		hzPos[ii] = 0;// (ii - PN / 2) * 0.0015;
//	}
//
//
//
//
//	float* hvol = new float[XN * YN * ZN];
//	float* hprj = new float[DNU * DNV * PN];
//
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != XN * YN; ++i)
//	{
//		mask[i] = 1;
//	}
//
//
//	for (int ii = 0; ii != DNU * DNV * PN; ++ii)
//	{
//		hprj[ii] = 1.0;
//	}
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//
//		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs + ii, hzPos, 1,
//			XN, YN, ZN, hvol, hprj + ii * DNU * DNV, dx, dz, mask, 0, 0, 0);
//		std::cout << ".";
//		std::stringstream ss;
//		ss << ii;
//		std::string name = "bak" + ss.str() + ".raw";
//		std::ofstream fou3(name.c_str(), std::ios::binary);
//		fou3.write((char*)hvol, sizeof(float) * XN * YN * ZN);
//		fou3.close();
//
//	}
//}
//
//void testBackprojection2()
//{
//
//	float sid = 541.0f;
//	float sdd = 949.0f;
//
//	float x0(0.0f);
//	float y0(sid);
//	float z0(0.0f);
//
//
//	int DNU(888);
//	int DNV(64);
//	int PN(1200);
//	float imgXCenter(0);
//	float imgYCenter(0);
//	float imgZCenter(0);
//
//	int XN(512);
//	int YN(512);
//	int ZN(64);
//
//
//	float dx(500.0f / 512.0f);
//	float dz(0.625);
//
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//	//Generate the positions of the detectors
//	float col_size = 1.0239;
//	float row_size = 1.0963;
//
//	float col_offset = 0;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//	imgXCenter = 0;
//	imgYCenter = 0;
//	imgZCenter = 0;
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(PN);
//		hzPos[ii] = 0;// (ii - PN / 2) * 0.0015;
//	}
//
//
//
//
//	float* hvol = new float[XN * YN * ZN];
//	float* hprj = new float[DNU * DNV * PN];
//
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != XN * YN; ++i)
//	{
//		mask[i] = 1;
//	}
//
//
//	for (int ii = 0; ii != DNU * DNV * PN; ++ii)
//	{
//		hprj[ii] = 1.0;
//	}
//
//	DD3Back_3gpus(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, hvol, hprj, dx, dz, mask, 0, 0, 0);
//	std::string name = "bak.raw";
//	std::ofstream fou3(name.c_str(), std::ios::binary);
//	fou3.write((char*) hvol, sizeof(float) * XN * YN * ZN);
//	fou3.close();
//}
//
//
//
//
//void testProjection()
//{
//	float sid = 541.0f;
//	float sdd = 949.0f;
//
//	float x0(0.0f);
//	float y0(sid);
//	float z0(0.0f);
//
//
//	int DNU(888);
//	int DNV(64);
//	int PN(984);
//	float imgXCenter(0);
//	float imgYCenter(0);
//	float imgZCenter(0);
//
//	int XN(512);
//	int YN(512);
//	int ZN(64);
//
//
//	float dx(500.0f / 512.0f);
//	float dz(0.625);
//
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//	//Generate the positions of the detectors
//	float col_size = 1.0239;
//	float row_size = 1.0963;
//
//	float col_offset = 0;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//	imgXCenter = 0;
//	imgYCenter = 0;
//	imgZCenter = 0;
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(PN);
//		hzPos[ii] = 0;// (ii - PN / 2) * 0.0015;
//	}
//
//	float* hvol = new float[XN * YN * ZN];
//	float* hprj = new float[DNU * DNV * PN];
//	for (int i = 0; i != XN * YN * ZN; ++i)
//	{
//		hvol[i] = 1.0;
//	}
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != XN * YN; ++i)
//	{
//		mask[i] = 1;
//	}
//
//
//	for (int ii = 0; ii != DNU * DNV * PN; ++ii)
//	{
//		hprj[ii] = 0.0;
//	}
//
//
//	DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter,
//		hangs, hzPos, PN, XN, YN, ZN, hvol, hprj,
//		dx, dz, mask, 0, 0);
//
//	DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter,
//		hangs, hzPos, PN, XN, YN, ZN, hvol, hprj,
//		dx, dz, mask, 0, 1);
//
//	DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter,
//		hangs, hzPos, PN, XN, YN, ZN, hvol, hprj,
//		dx, dz, mask, 0, 3);
//
//	DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter,
//		hangs, hzPos, PN, XN, YN, ZN, hvol, hprj,
//		dx, dz, mask, 0, 4);
////	std::string name = "testProj20160229.prj";
////	std::ofstream fou3(name.c_str(), std::ios::binary);
////	fou3.write((char*) hprj, sizeof(float) * DNU * DNV * PN);
////	fou3.close();
//}
//
//
//
////Test the correctness of panel detector.
//void testPanelDetectorBackprojection()
//{
//
//	float sid = 541.0f;
//	float sdd = 949.0f;
//
//	float x0(0.0f);
//	float y0(sid);
//	float z0(0.0f);
//
//
//	int DNU(888);
//	int DNV(64);
//	int PN(984);
//	float imgXCenter(0);
//	float imgYCenter(0);
//	float imgZCenter(0);
//
//	int XN(512);
//	int YN(512);
//	int ZN(64);
//
//
//	float dx(500.0f / 512.0f);
//	float dz(0.625);
//
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//	//Generate the positions of the detectors
//	float col_size = 1.0239;
//	float row_size = 1.0963;
//
//	//float col_offset = -1.25;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	//float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		xds[ii] = (ii - DNU / 2) * 500.0 / 888.0;
//		yds[ii] = 541.0 - 949.0;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//	imgXCenter = 0;
//	imgYCenter = 0;
//	imgZCenter = 0;
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(PN);
//		hzPos[ii] = 0;// (ii - PN / 2) * 0.0015;
//	}
//
//
//
//
//	float* hvol = new float[XN * YN * ZN];
//	float* hprj = new float[DNU * DNV * PN];
//
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != XN * YN; ++i)
//	{
//		mask[i] = 1;
//	}
//
//
//	for (int ii = 0; ii != DNU * DNV * PN; ++ii)
//	{
//		hprj[ii] = 1.0;
//	}
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//
//
//		DD3_panel_gpu_back_branchless_sat2d( x0, y0, z0,
//			DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs + ii, hzPos, 1,
//			XN, YN, ZN, hvol, hprj + ii * DNU * DNV, dx, dz, mask, 0, 0);
//
//		std::cout << ".";
//		std::stringstream ss;
//		ss << ii;
//		std::string name = "bak" + ss.str() + ".raw";
//		std::ofstream fou3(name.c_str(), std::ios::binary);
//		fou3.write((char*)hvol, sizeof(float) * XN * YN * ZN);
//		fou3.close();
//
//	}
//}
//
//
////Generate the projection with mouse volume
//void genMobyMouseProj()
//{
//
//	float sid = 541.0f;
//	float sdd = 949.0f;
//
//	float x0(0.0f);
//	float y0(sid);
//	float z0(0.0f);
//
//
//	int DNU(888 * 4); // will shrink to 888
//	int DNV(64 * 4); // will shrink to 64
//	int PN(1200);
//	float imgXCenter(0);
//	float imgYCenter(0);
//	float imgZCenter(0);
//
//	int XN(512);
//	int YN(512);
//	int ZN(128); //mouse volume is larger
//
//
//	float dx(500.0f / 512.0f);
//	float dz(0.625);
//
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//	//Generate the positions of the detectors
//	float col_size = 1.0239 / 4;
//	float row_size = 1.0963 / 4;
//
//	float col_offset = 0;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//	imgXCenter = 0;
//	imgYCenter = 0;
//	imgZCenter = 0;
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(PN);
//		hzPos[ii] = 0;// (ii - PN / 2) * 0.0015;
//	}
//
//	float* hvol = new float[XN * YN * ZN];
//	float* hprj = new float[DNU * DNV]; // Each time we only perform on projection
//
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != XN * YN; ++i)
//	{
//		mask[i] = 1;
//	}
//
//
//	for (int ii = 0; ii != DNU * DNV; ++ii)
//	{
//		hprj[ii] = 1.0;
//	}
//
//
//	for (int chIdx = 127; chIdx <= 128; ++chIdx)
//	{
//		std::stringstream ss;
//		ss << chIdx;
//		std::string volumeName = "D:\\MobyMouseVolume\\moby_" + ss.str() + ".raw";
//
//		std::ifstream fin(volumeName.c_str(), std::ios::binary);
//		fin.read((char*)hvol, sizeof(float) * XN * YN * ZN);
//		fin.close();
//
//		for (int ii = 0; ii != PN; ++ii)
//		{
//			memset(hprj, 0, sizeof(float) * DNU * DNV);
//
//			DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//				imgXCenter, imgYCenter, imgZCenter, hangs + ii, hzPos, 1,
//				XN, YN, ZN, hvol, hprj, dx, dz, mask, 0, 0);
//
//			std::stringstream ss2;
//			ss2 << ii;
//			std::string prjName = "D:\\MobyMouseVolume\\Projection\\prj_" + ss.str() + "\\proj_" + ss2.str() + ".prj";
//			std::ofstream fou(prjName.c_str(), std::ios::binary);
//			fou.write((char*)hprj, sizeof(float) * DNU * DNV);
//			fou.close();
//		}
//		std::cout << "channel " << chIdx << "\n";
//	}
//
//}
//
//#include <thrust/transform.h>
//template<typename T>
//struct FISTA_functor
//{
//	T t1;
//	T t2;
//	T minV;
//	T maxV;
//	FISTA_functor(const T& _t1, const T& _t2, const T& _minV, const T& _maxV) :t1(_t1), t2(_t2), minV(_minV), maxV(_maxV){}
//	__host__ __device__ T operator()(T curImg, T lasImg)
//	{
//		T res = curImg + (t1 - 1.0) / t2 * (curImg - lasImg);
//		if (res < minV)
//		{
//			return minV;
//		}
//		else if (res > maxV)
//		{
//			return maxV;
//		}
//		else
//		{
//			return res;
//		}
//	}
//
//};
//
//template<typename T>
//void FISTA(T* lasImg, T* currentImg, T t1, T t2, int N, T minV = -2000, T maxV = 6000)
//{
//	thrust::transform(currentImg, currentImg + N, lasImg, currentImg, FISTA_functor<float>(t1, t2, minV, maxV));
//}
//
//
//template<typename T>
//void FISTA(thrust::device_vector<T>& lasImg, thrust::device_vector<T>& currentImg, T t1, T t2, T minV, T maxV)
//{
//
//	thrust::transform(currentImg.begin(), currentImg.end(), lasImg.begin(), currentImg.begin(), FISTA_functor<T>(t1, t2, minV, maxV));
//}
//
//
//void reconMobyMouse()
//{
//
//	float sid = 541.0f;
//	float sdd = 949.0f;
//
//	float x0(0.0f);
//	float y0(sid);
//	float z0(0.0f);
//
//	int DNU(888);
//	int DNV(64);
//	int PN(1200);
//	float imgXCenter(0);
//	float imgYCenter(0);
//	float imgZCenter(0);
//
//	int XN(512);
//	int YN(512);
//	int ZN(64); //mouse volume is larger
//
//
//	float dx(500.0f / 512.0f);
//	float dz(0.625);
//
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//	//Generate the positions of the detectors
//	float col_size = 1.0239;
//	float row_size = 1.0963;
//
//	float col_offset = 0;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//	imgXCenter = 0;
//	imgYCenter = 0;
//	imgZCenter = 0;
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(PN);
//		hzPos[ii] = 0;// (ii - PN / 2) * 0.0015;
//	}
//
//	float* hvol = new float[XN * YN * ZN];
//	float* hprj = new float[DNU * DNV * PN]; // Each time we only perform on projection
//
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != YN; ++i)
//	{
//		for (int j = 0; j != XN; ++j)
//		{
//			if (sqrtf(powf((i - YN / 2.0 + 0.5) / (YN / 2.0), 2.0) + powf((j - XN / 2.0 + 0.5) / (XN / 2.0), 2.0)) < 0.98)
//			{
//				mask[i * XN + j] = 1;
//			}
//			else
//			{
//				mask[i * XN + j] = 0;
//			}
//		}
//	}
//
//	std::ifstream fin("D:\\MobyMouseVolume\\Projection\\prj_1\\prj.raw", std::ios::binary);
//	fin.read((char*)hprj, sizeof(float) * DNU * DNV * PN);
//	fin.close();
//
//	float* oneVol = new float[XN * YN * ZN];
//	float* onePrj = new float[DNU * DNV * PN];
//	float* rowSum = new float[DNU * DNV * PN];
//	float* colSum = new float[XN * YN * ZN];
//	float* lasImg = new float[XN * YN * ZN];
//
//	for (int ii = 0; ii != XN * YN * ZN; ++ii)
//	{
//		hvol[ii] = 0;
//		oneVol[ii] = 1.0;
//		colSum[ii] = 0;
//	}
//	for (int ii = 0; ii != DNU * DNV * PN; ++ii)
//	{
//		onePrj[ii] = 1.0;
//		rowSum[ii] = 0.0;
//	}
//
//	DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, oneVol, rowSum, dx, dz, mask, 0, 0);
//
//	DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, colSum, onePrj, dx, dz, mask, 0, 0, 0);
//	float t1 = 1.0;
//	float t2 = 1.0;
//	int totN = XN * YN * ZN;
//	for (int ii = 0; ii != 500; ++ii)
//	{
//		memcpy(lasImg, hvol, sizeof(float) * totN);
//
//		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, hvol, onePrj, dx, dz, mask, 0, 0);
//
//		prjWeight(onePrj, hprj, rowSum, DNU*DNV*PN);
//
//		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, oneVol, onePrj, dx, dz, mask, 0, 0, 0);
//
//		bakWeight(oneVol, hvol, colSum, totN);
//
//		t2 = (1.0 + sqrtf(1.0 + 4.0 * t1 * t1)) / 2.0;
//		FISTA(lasImg, hvol, t1, t2, totN);
//		t1 = t2;
//
//		std::cout << ".";
//		if ((ii + 1) % 10 == 0)
//		{
//			std::stringstream ss;
//			ss << ii;
//			std::string name = "rev" + ss.str() + ".raw";
//			std::ofstream fou3(name.c_str(), std::ios::binary);
//			fou3.write((char*)hvol, sizeof(float) * XN * YN * ZN);
//			fou3.close();
//		}
//	}
//}
//void reconSheppLogan()
//{
//	float sid = 541.0f;
//	float sdd = 949.0f;
//
//	float x0(0.0f);
//	float y0(sid);
//	float z0(0.0f);
//
//	int DNU(888);
//	int DNV(64);
//	int PN(984);
//	float imgXCenter(0);
//	float imgYCenter(0);
//	float imgZCenter(0);
//
//	int XN(512);
//	int YN(512);
//	int ZN(64); //mouse volume is larger
//
//
//	float dx(500.0f / 512.0f);
//	float dz(0.625);
//
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//	//Generate the positions of the detectors
//	float col_size = 1.0239;
//	float row_size = 1.0963;
//
//	float col_offset = 0;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//	imgXCenter = 0;
//	imgYCenter = 0;
//	imgZCenter = 0;
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(PN);
//		hzPos[ii] = 0;// (ii - PN / 2) * 0.0015;
//	}
//
//	float* hvol = new float[XN * YN * ZN];
//	float* hprj = new float[DNU * DNV * PN]; // Each time we only perform on projection
//
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != YN; ++i)
//	{
//		for (int j = 0; j != XN; ++j)
//		{
//			if (sqrtf(powf((i - YN / 2.0 + 0.5) / (YN / 2.0), 2.0) + powf((j - XN / 2.0 + 0.5) / (XN / 2.0), 2.0)) < 0.98)
//			{
//				mask[i * XN + j] = 1;
//			}
//			else
//			{
//				mask[i * XN + j] = 0;
//			}
//		}
//	}
//
//	std::ifstream fin("testProj.prj", std::ios::binary);
//	fin.read((char*) hprj, sizeof(float) * DNU * DNV * PN);
//	fin.close();
//
//	float* oneVol = new float[XN * YN * ZN];
//	float* onePrj = new float[DNU * DNV * PN];
//	float* rowSum = new float[DNU * DNV * PN];
//	float* colSum = new float[XN * YN * ZN];
//	float* lasImg = new float[XN * YN * ZN];
//
//	for (int ii = 0; ii != XN * YN * ZN; ++ii)
//	{
//		hvol[ii] = 0;
//		oneVol[ii] = 1.0;
//		colSum[ii] = 0;
//	}
//	for (int ii = 0; ii != DNU * DNV * PN; ++ii)
//	{
//		onePrj[ii] = 1.0;
//		rowSum[ii] = 0.0;
//	}
//
//	DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, oneVol, rowSum, dx, dz, mask, 0, 0);
//
//	DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, colSum, onePrj, dx, dz, mask, 1, 0, 0);
//	float t1 = 1.0;
//	float t2 = 1.0;
//	int totN = XN * YN * ZN;
//	for (int ii = 0; ii != 2; ++ii)
//	{
//		memcpy(lasImg, hvol, sizeof(float) * totN);
//
//		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, hvol, onePrj, dx, dz, mask, 2, 0);
//
//		prjWeight(onePrj, hprj, rowSum, DNU*DNV*PN);
//
//		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, oneVol, onePrj, dx, dz, mask, 3, 0, 0);
//
//		bakWeight(oneVol, hvol, colSum, totN);
//
//		t2 = (1.0 + sqrtf(1.0 + 4.0 * t1 * t1)) / 2.0;
//		FISTA(lasImg, hvol, t1, t2, totN);
//		t1 = t2;
//
//		std::cout << ".";
//		if ((ii + 1) % 10 == 0)
//		{
//			std::stringstream ss;
//			ss << ii;
//			std::string name = "SheppRecon" + ss.str() + ".raw";
//			std::ofstream fou3(name.c_str(), std::ios::binary);
//			fou3.write((char*) hvol, sizeof(float) * XN * YN * ZN);
//			fou3.close();
//		}
//	}
//}
//
//
//
//std::string num2str(int s)
//{
//	std::stringstream ss;
//	ss << s;
//	return ss.str();
//}
//void reconMobyMouse(
//	int chIdx,
//	int iterNum,
//	int storePeriod,
//	int gpuID)
//{
//
//	float sid = 541.0f;
//	float sdd = 949.0f;
//
//	float x0(0.0f);
//	float y0(sid);
//	float z0(0.0f);
//
//	int DNU(888);
//	int DNV(64);
//	int PN(1200);
//	float imgXCenter(0);
//	float imgYCenter(0);
//	float imgZCenter(0);
//
//	int XN(512);
//	int YN(512);
//	int ZN(64); //mouse volume is larger
//
//
//	float dx(500.0f / 512.0f);
//	float dz(0.625);
//
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//	//Generate the positions of the detectors
//	float col_size = 1.0239;
//	float row_size = 1.0963;
//
//	float col_offset = 0;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//	imgXCenter = 0;
//	imgYCenter = 0;
//	imgZCenter = 0;
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(PN);
//		hzPos[ii] = 0;// (ii - PN / 2) * 0.0015;
//	}
//
//	float* hvol = new float[XN * YN * ZN];
//	float* hprj = new float[DNU * DNV * PN]; // Each time we only perform on projection
//
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != YN; ++i)
//	{
//		for (int j = 0; j != XN; ++j)
//		{
//			if (sqrtf(powf((i - YN / 2.0 + 0.5) / (YN / 2.0), 2.0) + powf((j - XN / 2.0 + 0.5) / (XN / 2.0), 2.0)) < 0.98)
//			{
//				mask[i * XN + j] = 1;
//			}
//			else
//			{
//				mask[i * XN + j] = 0;
//			}
//		}
//	}
//
//	std::string prjName = "D:\\MobyMouseVolume\\Projection\\ch" + num2str(chIdx) + "\\noisyProj1000000.raw";
//	std::ifstream fin(prjName.c_str(), std::ios::binary);
//	fin.read((char*) hprj, sizeof(float) * DNU * DNV * PN);
//	fin.close();
//
//	float* oneVol = new float[XN * YN * ZN];
//	float* onePrj = new float[DNU * DNV * PN];
//	float* rowSum = new float[DNU * DNV * PN];
//	float* colSum = new float[XN * YN * ZN];
//	float* lasImg = new float[XN * YN * ZN];
//
//	for (int ii = 0; ii != XN * YN * ZN; ++ii)
//	{
//		hvol[ii] = 0;
//		oneVol[ii] = 1.0;
//		colSum[ii] = 0;
//	}
//	for (int ii = 0; ii != DNU * DNV * PN; ++ii)
//	{
//		onePrj[ii] = 1.0;
//		rowSum[ii] = 0.0;
//	}
//
//	DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, oneVol, rowSum, dx, dz, mask, gpuID, 0);
//
//	DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, colSum, onePrj, dx, dz, mask, gpuID, 0, 0);
//	float t1 = 1.0;
//	float t2 = 1.0;
//	int totN = XN * YN * ZN;
//	for (int ii = 0; ii != iterNum; ++ii)
//	{
//		memcpy(lasImg, hvol, sizeof(float) * totN);
//
//		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, hvol, onePrj, dx, dz, mask, gpuID, 0);
//
//		prjWeight(onePrj, hprj, rowSum, DNU*DNV*PN);
//
//		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, oneVol, onePrj, dx, dz, mask, gpuID, 0, 0);
//
//		bakWeight(oneVol, hvol, colSum, totN);
//
//		t2 = (1.0 + sqrtf(1.0 + 4.0 * t1 * t1)) / 2.0;
//		FISTA(lasImg, hvol, t1, t2, totN);
//		t1 = t2;
//
//		std::cout << ".";
//		if ((ii) % storePeriod == 0)
//		{
//			std::stringstream ss;
//			ss << ii;
//			std::string name = "D:\\MobyMouseVolume\\Projection\\ch" + num2str(chIdx) + "\\3Drecon\\recon_" + num2str(ii) + "from1000000pho.img";
//			std::ofstream fou3(name.c_str(), std::ios::binary);
//			fou3.write((char*) hvol, sizeof(float) * XN * YN * ZN);
//			fou3.close();
//		}
//	}
//}
//void reconFORBILD()
//{
//	// Set the geometry
//	float sid = 541.0f;
//	float sdd = 949.0f;
//
//
//	float x0(0.0f);
//	float y0(sid);
//	float z0(0.0f);
//
//
//	int DNU(888);
//	int DNV(64);
//	int PN(2200);
//	float imgXCenter(0);
//	float imgYCenter(0);
//	float imgZCenter(0);
//
//	int XN(512);
//	int YN(512);
//	int ZN(128);
//	float dx(500.0f / 512.0f);
//	float dz(0.625);
//
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//	//Generate the positions of the detectors
//	float col_size = 1.0239;
//	float row_size = 1.0963;
//
//	float col_offset = -1.25;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//	imgXCenter = 0;
//	imgYCenter = 0;
//	imgZCenter = 0;
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(PN);
//		hzPos[ii] = 0;// (ii - PN / 2) * 0.0015;
//	}
//
//	float* hvol = new float[XN * YN * ZN];
//	for(int i = 0; i != XN * YN * ZN; ++i)
//	{
//		hvol[i] = 1.0;
//	}
//	float* hprj = new float[DNU * DNV * PN];
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != XN * YN; ++i)
//	{
//		mask[i] = 1;
//	}
//
//	std::ifstream fou1("FORBILD_proj.prj", std::ios::binary);
//	fou1.read((char*)hprj, sizeof(float) * DNU * DNV * PN);
//	fou1.close();
//
//	//Cone beam recon
//
//
//	float* oneVol = new float[XN * YN * ZN];
//	float* onePrj = new float[DNU * DNV * PN];
//	float* rowSum = new float[DNU * DNV * PN];
//	float* colSum = new float[XN * YN * ZN];
//
//	for (int ii = 0; ii != XN * YN * ZN; ++ii)
//	{
//		hvol[ii] = 0;
//		oneVol[ii] = 1.0;
//		colSum[ii] = 0;
//	}
//	for (int ii = 0; ii != DNU * DNV * PN; ++ii)
//	{
//		onePrj[ii] = 1.0;
//		rowSum[ii] = 0.0;
//	}
//
//	DD3Proj_3gpus(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, oneVol, rowSum, dx, dz, mask, 0, 0);
//	std::cout<<"Proj sum finished\n";
//	DD3Back_3gpus(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, colSum, onePrj, dx, dz, mask, 0, 0, 0);
//
//	std::ofstream bakSum("BakSum.raw",std::ios::binary);
//	bakSum.write((char*)colSum,sizeof(float) * XN * YN * ZN);
//	bakSum.close();
//	std::cout<<"Backproj sum finished\n";
//
//
//
//	for (int ii = 0; ii != 400; ++ii)
//	{
//		DD3Proj_3gpus(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, hvol, onePrj, dx, dz, mask, 0, 0);
//
//		prjWeight(onePrj, hprj, rowSum, DNU*DNV*PN);
//
//		DD3Back_3gpus(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, oneVol, onePrj, dx, dz, mask, 0, 0, 0);
//
//		bakWeight(oneVol, hvol, colSum, XN * YN * ZN);
//		if((ii + 1)%5  == 0)
//		{
//			std::stringstream ss;
//			ss<<ii;
//			std::string name = "forbildRecon" + ss.str() + ".raw";
//			std::ofstream fou3(name.c_str(), std::ios::binary);
//			fou3.write((char*)hvol, sizeof(float) * XN * YN * ZN);
//			fou3.close();
//
//		}
//		std::cout<<ii<<"\n";
//	}
//
//}
//
//
//template<typename T>
//struct FISTA_20151105
//{
//	T t0;
//	T t1;
//	FISTA_20151105(const T& t_0, const T& t_1):t0(t_0),t1(t_1){}
//	T operator()(const T& lastV, const T& curV)
//	{
//		T res =curV + (t0 - 1.0) / t1 * (curV - lastV);
//		if(res<0)
//			return 0;
//		else
//			return res;
//	}
//};
//
//
//void reconFORBILD20151105()
//{
//	// Set the geometry
//	float sid = 541.0f;
//	float sdd = 949.0f;
//
//
//	float x0(0.0f);
//	float y0(sid);
//	float z0(0.0f);
//
//
//	int DNU(888);
//	int DNV(64);
//	int PN(2200);
//	float imgXCenter(0);
//	float imgYCenter(0);
//	float imgZCenter(0);
//
//	int XN(512);
//	int YN(512);
//	int ZN(128);
//	float dx(500.0f / 512.0f);
//	float dz(0.625);
//
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//	//Generate the positions of the detectors
//	float col_size = 1.0239;
//	float row_size = 1.0963;
//
//	float col_offset = -1.25;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//	imgXCenter = 0;
//	imgYCenter = 0;
//	imgZCenter = 0;
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(PN);
//		hzPos[ii] = 0;// (ii - PN / 2) * 0.0015;
//	}
//
//	float* hvol = new float[XN * YN * ZN];
//	for(int i = 0; i != XN * YN * ZN; ++i)
//	{
//		hvol[i] = 1.0;
//	}
//	float* hprj = new float[DNU * DNV * PN];
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != XN * YN; ++i)
//	{
//		mask[i] = 1;
//	}
//
//	std::ifstream fou1("forbildReal.raw", std::ios::binary);
//	fou1.read((char*)hvol, sizeof(float) * XN * YN * ZN);
//	fou1.close();
//	DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, hvol, hprj, dx, dz, mask, 0, 0);
//	//Cone beam recon
//
//
//	float* oneVol = new float[XN * YN * ZN];
//	float* onePrj = new float[DNU * DNV * PN];
//	float* rowSum = new float[DNU * DNV * PN];
//	float* colSum = new float[XN * YN * ZN];
//
//	for (int ii = 0; ii != XN * YN * ZN; ++ii)
//	{
//		hvol[ii] = 0;
//		oneVol[ii] = 1.0;
//		colSum[ii] = 0;
//	}
//	for (int ii = 0; ii != DNU * DNV * PN; ++ii)
//	{
//		onePrj[ii] = 1.0;
//		rowSum[ii] = 0.0;
//	}
//
//	DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, oneVol, rowSum, dx, dz, mask, 0, 0);
//	std::cout<<"Proj sum finished\n";
//	DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, colSum, onePrj, dx, dz, mask, 0, 0, 0);
//
//	std::ofstream bakSum("BakSum.raw",std::ios::binary);
//	bakSum.write((char*)colSum,sizeof(float) * XN * YN * ZN);
//	bakSum.close();
//	std::cout<<"Backproj sum finished\n";
//
//
//	thrust::host_vector<float> lastVol(hvol,hvol+XN * YN * ZN);
//	float t0 = 1.0f;
//	float t1 = 1.0f;
//	for (int ii = 0; ii != 800; ++ii)
//	{
//		thrust::copy(hvol,hvol+XN*YN*ZN,lastVol.begin());
//
//		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, hvol, onePrj, dx, dz, mask, 0, 0);
//
//		prjWeight(onePrj, hprj, rowSum, DNU*DNV*PN);
//
//		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, oneVol, onePrj, dx, dz, mask, 0, 0, 0);
//
//		bakWeight(oneVol, hvol, colSum, XN * YN * ZN);
//
//		t1 = (1.0 + sqrtf(1.0 + 4.0 * t0 * t0)) * 0.5;
//
//		thrust::transform(lastVol.begin(),lastVol.end(),hvol,hvol,FISTA_20151105<float>(t0,t1));
//		t0 = t1;
//		if(ii % 10  == 0)
//		{
//			std::stringstream ss;
//			ss<<ii;
//			std::string name = "forbildRecon" + ss.str() + ".raw";
//			std::ofstream fou3(name.c_str(), std::ios::binary);
//			fou3.write((char*)hvol, sizeof(float) * XN * YN * ZN);
//			fou3.close();
//
//		}
//		std::cout<<ii<<"\n";
//	}
//
//}
//
//
//
//void projFORBILD()
//{
//	// Set the geometry
//	float sid = 541.0f;
//	float sdd = 949.0f;
//
//
//	float x0(0.0f);
//	float y0(sid);
//	float z0(0.0f);
//
//
//	int DNU(888);
//	int DNV(64);
//	int PN(2200);
//	float imgXCenter(0);
//	float imgYCenter(0);
//	float imgZCenter(0);
//
//	int XN(512);
//	int YN(512);
//	int ZN(512);
//	float dx(500.0f / 512.0f);
//	float dz(0.625);
//
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//	//Generate the positions of the detectors
//	float col_size = 1.0239;
//	float row_size = 1.0963;
//
//	float col_offset = -1.25;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//	imgXCenter = 0;
//	imgYCenter = 0;
//	imgZCenter = 0;
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(PN);
//		hzPos[ii] = 0;// (ii - PN / 2) * 0.0015;
//	}
//
//	float* hvol = new float[XN * YN * ZN];
//	std::ifstream fin("ForbildHead.raw",std::ios::binary);
//	fin.read((char*)hvol, sizeof(float) * XN * YN * ZN);
//	fin.close();
//
//
//	float* hprj = new float[DNU * DNV * PN];
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != XN * YN; ++i)
//	{
//		mask[i] = 1;
//	}
//
//	DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, hvol, hprj, dx, dz, mask, 0, 0);
//
//	std::ofstream fou3("FORBILD_proj.prj", std::ios::binary);
//	fou3.write((char*)hprj, sizeof(float) * DNU * DNV * PN);
//	fou3.close();
//}
//
//
//
//void reconBigPatient()
//{
//
//
//	//
//	//
//	//// Reconstruct the big patient with GPU Branchless DD
//	//	float sid = 541.0f;
//	//	float sdd = 949.0f;
//	//
//	//
//	//	float x0(0.0f);
//	//	float y0(sid);
//	//	float z0(0.0f);
//	//
//	//
//	//	int DNU(888);
//	//	int DNV(64);
//	//	int PN(3210);
//	//	float imgXCenter(0);
//	//	float imgYCenter(0);
//	//	float imgZCenter(0);
//	//
//	//	int XN(512);
//	//	int YN(512);
//	//	int ZN(64 * 5);
//	//
//	//
//	//	float dx(700.0f / 512.0f);
//	//	float dz(0.625);
//	//
//	//
//	//	float* xds = new float[DNU];
//	//	float* yds = new float[DNU];
//	//	float* zds = new float[DNV];
//	//	//Generate the positions of the detectors
//	//	float col_size = 1.0239;
//	//	float row_size = 1.0963;
//	//
//	//	float col_offset = -1.25;
//	//	float row_offset = 0;
//	//
//	//
//	//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	//	float curBeta = 0;
//	//	for (int ii = 0; ii != DNU; ++ii)
//	//	{
//	//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//	//		xds[ii] = sinf(curBeta) * sdd;
//	//		yds[ii] = sid - cosf(curBeta) * sdd;
//	//	}
//	//
//	//	for (int ii = 0; ii != DNV; ++ii)
//	//	{
//	//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	//	}
//	//
//	//	imgXCenter = 0;
//	//	imgYCenter = 0;
//	//	imgZCenter = 0;
//	//
//	//	float* hangs = new float[PN];
//	//	float* hzPos = new float[PN];
//	//
//	//	for (int ii = 0; ii != PN; ++ii)
//	//	{
//	//		hangs[ii] = ii *  0.006385350921930 - 3.295689518577625;
//	//		hzPos[ii] = (ii - PN / 2) * 0.0559;
//	//	}
//	//
//	//
//	//
//	//
//	//	float* hvol = new float[XN * YN * ZN];
//	//	float* hprj = new float[DNU * DNV * PN];
//	//
//	//	std::ifstream fou1("prj3210_v2.raw", std::ios::binary);
//	//	fou1.read((char*)hprj, sizeof(float) * DNU * DNV * PN);
//	//	fou1.close();
//	//	byte* mask = new byte[XN * YN];
//	//	//for (int i = 0; i != YN; ++i)
//	//	//{
//	//	//	for (int j = 0; j != XN; ++j)
//	//	//	{
//	//	//		if (sqrtf(powf((i - YN / 2.0 + 0.5) / (YN / 2.0), 2.0) + powf((j - XN / 2.0 + 0.5) / (XN / 2.0), 2.0)) < 0.98)
//	//	//		{
//	//	//			mask[i * XN + j] = 1;
//	//	//		}
//	//	//		else
//	//	//		{
//	//	//			mask[i * XN + j] = 0;
//	//	//		}
//	//	//	}
//	//	//
//	//	//}
//	//
//	//	for (int i = 0; i != YN; ++i)
//	//	{
//	//		for (int j = 0; j != XN; ++j)
//	//		{
//	//			mask[i * XN + j] = 1;
//	//		}
//	//	}
//	//
//	//	float* oneVol = new float[XN * YN * ZN];
//	//	float* onePrj = new float[DNU * DNV * PN];
//	//	float* rowSum = new float[DNU * DNV * PN];
//	//	float* colSum = new float[XN * YN * ZN];
//	//
//	//	for (int ii = 0; ii != XN * YN * ZN; ++ii)
//	//	{
//	//		hvol[ii] = 0;
//	//		oneVol[ii] = 1.0;
//	//		colSum[ii] = 0;
//	//	}
//	//	for (int ii = 0; ii != DNU * DNV * PN; ++ii)
//	//	{
//	//		onePrj[ii] = 1.0;
//	//		rowSum[ii] = 0.0;
//	//	}
//	//
//	//	DD3Proj_3gpus(x0, y0, z0, DNU, DNV, xds, yds, zds,
//	//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//	//		XN, YN, ZN, oneVol, rowSum, dx, dz, mask, 0, 0);
//	//
//	//	DD3Back_3gpus(x0, y0, z0, DNU, DNV, xds, yds, zds,
//	//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//	//		XN, YN, ZN, colSum, onePrj, dx, dz, mask, 0, 0, 0);
//	//
//	//	for (int ii = 0; ii != 500; ++ii)
//	//	{
//	//		DD3Proj_3gpus(x0, y0, z0, DNU, DNV, xds, yds, zds,
//	//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//	//			XN, YN, ZN, hvol, onePrj, dx, dz, mask, 0, 0);
//	//
//	//		prjWeight(onePrj, hprj, rowSum, DNU*DNV*PN);
//	//
//	//		DD3Back_3gpus(x0, y0, z0, DNU, DNV, xds, yds, zds,
//	//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//	//			XN, YN, ZN, oneVol, onePrj, dx, dz, mask, 0, 0, 0);
//	//
//	//		bakWeight(oneVol, hvol, colSum, XN * YN * ZN);
//	//		std::cout << ".";
//	//		if ((ii + 1) % 10 == 0)
//	//		{
//	//			std::stringstream ss;
//	//			ss << ii;
//	//			std::string name = "rev" + ss.str() + ".raw";
//	//			std::ofstream fou3(name, std::ios::binary);
//	//			fou3.write((char*)hvol, sizeof(float) * XN * YN * ZN);
//	//			fou3.close();
//	//		}
//	//
//	//	}
//	//
//
//		//std::ofstream fou3("rev.raw", std::ios::binary);
//		//fou3.write((char*)hvol, sizeof(float) * XN * YN * ZN);
//		//fou3.close();
//		//return 0;
//}
//
//
////Test the projection time of HD geometry of different methods
//void testProjTime()
//{
//	// Set the geometry
//	float sid = 541.0f;
//	float sdd = 949.0f;
//
//	float x0(0.0f);
//	float y0(sid);
//	float z0(0.0f);
//
//
//	int DNU(888);
//	int DNV(64);
//	int PN(984);
//	float imgXCenter(0);
//	float imgYCenter(0);
//	float imgZCenter(0);
//
//	int XN(512);
//	int YN(512);
//	int ZN(64);
//	float dx(500.0f / 512.0f);
//	float dz(0.625);
//
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//	//Generate the positions of the detectors
//	float col_size = 1.0239;
//	float row_size = 1.0963;
//
//	float col_offset = -1.25;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//	imgXCenter = 0;
//	imgYCenter = 0;
//	imgZCenter = 0;
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(PN);
//		hzPos[ii] = 0;// (ii - PN / 2) * 0.0015;
//	}
//
//	float* hvol = new float[XN * YN * ZN];
//	for(int i = 0 ;i != XN * YN * ZN; ++i)
//	{
//		hvol[i] = 1.0;
//	}
//
//	float* hprj = new float[DNU * DNV * PN];
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != XN * YN; ++i)
//	{
//		mask[i] = 1;
//	}
//
//	//Proj with DD
//	DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, hvol, hprj, dx, dz, mask, 0, 0);
//	std::ofstream fou1("proj_DD.prj",std::ios::binary);
//	fou1.write((char*)hprj,sizeof(float) * DNU * DNV * PN);
//	fou1.close();
//	for(int i = 0; i != DNU * DNV * PN; ++i)
//	{
//		hprj[i] = 0;
//	}
//	//Proj with volume Rendering
//	DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, hvol, hprj, dx, dz, mask, 0, 1);
//	std::ofstream fou2("proj_VR.prj",std::ios::binary);
//	fou2.write((char*)hprj,sizeof(float) * DNU * DNV * PN);
//	fou2.close();
//	for(int i = 0; i != DNU * DNV * PN; ++i)
//	{
//		hprj[i] = 0;
//	}
//	//Proj with pseudo DD
//	DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, hvol, hprj, dx, dz, mask, 0, 3);
//	std::ofstream fou3("proj_PD.prj",std::ios::binary);
//	fou3.write((char*)hprj,sizeof(float) * DNU * DNV * PN);
//	fou3.close();
//	for(int i = 0; i != DNU * DNV * PN; ++i)
//	{
//		hprj[i] = 0;
//	}
//}
//
//void TestProjectionTime()
//{
//	//Set the HD Geometry
//	// Set the geometry
//	float sid = 541.0f;
//	float sdd = 949.0f;
//
//	float x0(0.0f);
//	float y0(sid);
//	float z0(0.0f);
//
//
//	int DNU(888);
//	int DNV(64);
//	int PN(984);
//	float imgXCenter(0);
//	float imgYCenter(0);
//	float imgZCenter(0);
//
//	int XN(512);
//	int YN(512);
//	int ZN(64);
//	float dx(500.0f / 512.0f);
//	float dz(0.625);
//
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//	//Generate the positions of the detectors
//	float col_size = 1.0239;
//	float row_size = 1.0963;
//
//	float col_offset = 0;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//	imgXCenter = 0;
//	imgYCenter = 0;
//	imgZCenter = 0;
//
//	//Test different view number with HD geometry
//
//	float* projTime = new float[440];
//	float* backTime = new float[440];
//	int prjNum = 0;
//
//
//	float* hvol = new float[XN * YN * ZN];
//	for(int i = 0 ;i != XN * YN * ZN; ++i)
//	{
//		hvol[i] = 1.0;
//	}
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != XN * YN; ++i)
//	{
//		mask[i] = 1;
//	}
//
//	for(prjNum = 10; prjNum <= 4400; prjNum+=10)
//	{
//		PN = prjNum;
//		float* hangs = new float[PN];
//		float* hzPos = new float[PN];
//
//		for (int ii = 0; ii != PN; ++ii)
//		{
//			hangs[ii] = ii * TWOPI / static_cast<float>(PN);
//			hzPos[ii] = 0;// (ii - PN / 2) * 0.0015;
//		}
//		float* hprj = new float[DNU * DNV * PN];
//
//		//Proj with DD
//		clock_t start;
//		clock_t stop;
//		start = clock();
//		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, hvol, hprj, dx, dz, mask, 0, 0);
//		stop = clock();
//		projTime[prjNum/10-1] = static_cast<double>(stop - start) / static_cast<double>(CLOCKS_PER_SEC);
//
//		start = clock();
//		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, hvol, hprj, dx, dz, mask, 0, 0, 0);
//		stop = clock();
//		backTime[prjNum / 10 - 1] = static_cast<double>(stop - start) / static_cast<double>(CLOCKS_PER_SEC);
//
//
//		delete[] hprj;
//		delete[] hangs;
//		delete[] hzPos;
//
//		std::cout<<"Projection Time : "<<projTime[prjNum/10-1]<<"\n";
//		std::cout<<"Backprojection Time : "<<backTime[prjNum/10-1]<<"\n";
//
//
//	}
//	delete[] mask;
//	delete[] hvol;
//	std::ofstream f1("prjTime.raw",std::ios::binary);
//	f1.write((char*)projTime,sizeof(float) * 440);
//	f1.close();
//	std::ofstream f2("bakTime.raw",std::ios::binary);
//	f2.write((char*)backTime,sizeof(float) * 440);
//	f2.close();
//
//}
//
//
//
//void TestProjectionTime2()
//{
//	float sid = 541.0f;
//	float sdd = 949.0f;
//
//	float x0(0.0f);
//	float y0(sid);
//	float z0(0.0f);
//
//
//	int DNU(888);
//	int DNV(64);
//	int PN(984);
//	float imgXCenter(0);
//	float imgYCenter(0);
//	float imgZCenter(0);
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//	//Generate the positions of the detectors
//	float col_size = 1.0239;
//	float row_size = 1.0963;
//
//	float col_offset = 0;
//	float row_offset = 0;
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//	imgXCenter = 0;
//	imgYCenter = 0;
//	imgZCenter = 0;
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(PN);
//		hzPos[ii] = 0;// (ii - PN / 2) * 0.0015;
//	}
//
//	float* projTime = new float[200];
//	float* backTime = new float[200];
//	float* hprj = new float[DNU * DNV * PN];
//	// Set the geometry
//	for(int si = 128; si <= 1536; si += 8)
//	{
//		int XN(si);
//		int YN(si);
//		int ZN(64);
//		float dx(500.0f / static_cast<float>(si));
//		float dz(0.625);
//
//
//		float* hvol = new float[XN * YN * ZN];
//		for(int i = 0 ;i != XN * YN * ZN; ++i)
//		{
//			hvol[i] = 1.0;
//		}
//
//
//		byte* mask = new byte[XN * YN];
//		for (int i = 0; i != XN * YN; ++i)
//		{
//			mask[i] = 1;
//		}
//
//		//Proj with DD
//		clock_t start;
//		clock_t stop;
//		start = clock();
//		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, hvol, hprj, dx, dz, mask, 0, 0);
//		stop = clock();
//		projTime[(si -128)/ 8] = static_cast<double>(stop - start) / static_cast<double>(CLOCKS_PER_SEC);
//
//		start = clock();
//		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, hvol, hprj, dx, dz, mask, 0, 0, 0);
//		stop = clock();
//		backTime[(si -128)/ 8] = static_cast<double>(stop - start) / static_cast<double>(CLOCKS_PER_SEC);
//
//		std::cout<<"Projection Time : "<<projTime[(si -128)/ 8]<<"\n";
//		std::cout<<"Backprojection Time : "<<backTime[(si -128)/ 8]<<"\n";
//
//		delete[] hvol;
//		delete[] mask;
//
//
//	}
//
//	std::ofstream f1("prjTime_128to1536.raw",std::ios::binary);
//	f1.write((char*)projTime,sizeof(float) * 200);
//	f1.close();
//	std::ofstream f2("bakTime_128to1536.raw",std::ios::binary);
//	f2.write((char*)backTime,sizeof(float) * 200);
//	f2.close();
//
//	delete[] xds;
//	delete[] yds;
//	delete[] zds;
//	delete[] hangs;
//	delete[] hzPos;
//
//	delete[] hprj;
//
//
//
//}
//
//
//void TestProjectionTime3()
//{
//	// Set the geometry
//	float sid = 541.0f;
//	float sdd = 949.0f;
//
//	float x0(0.0f);
//	float y0(sid);
//	float z0(0.0f);
//
//
//	int DNU(888);
//	int DNV(64);
//	int PN(2200);
//	float imgXCenter(0);
//	float imgYCenter(0);
//	float imgZCenter(0);
//
//	int XN(512);
//	int YN(512);
//	int ZN(64);
//	float dx(500.0f / 512.0f);
//	float dz(0.625);
//
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//	//Generate the positions of the detectors
//	float col_size = 1.0239;
//	float row_size = 1.0963;
//
//	float col_offset = 0;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//	imgXCenter = 0;
//	imgYCenter = 0;
//	imgZCenter = 0;
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(PN);
//		hzPos[ii] = 0;// (ii - PN / 2) * 0.0015;
//	}
//
//	float* hvol = new float[XN * YN * ZN];
//	for(int i = 0 ;i != XN * YN * ZN; ++i)
//	{
//		hvol[i] = 1.0;
//	}
//
//	float* hprj = new float[DNU * DNV * PN];
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != XN * YN; ++i)
//	{
//		mask[i] = 1;
//	}
//
//	//Proj with DD
//	clock_t start = clock();
//	DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, hvol, hprj, dx, dz, mask, 0, 0);
//	clock_t stop = clock();
//	std::cout<<"Projection Time is : " <<static_cast<float>(stop - start) / static_cast<float>(CLOCKS_PER_SEC)<<"\n";
//
//	start = clock();
//	DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, hvol, hprj, dx, dz, mask, 0, 0, 0);
//	stop = clock();
//	std::cout<<"Backprojection Time is : "<<static_cast<float>(stop - start) / static_cast<float>(CLOCKS_PER_SEC)<<"\n";
//
//}
//
//
//
//
//
////Test projection backprojection time with different number of GPUs
//void TestProjectionTimeWithNumOfGPUs(int SIZE1)
//{
//	//Set the HD Geometry
//	// Set the geometry
//	float sid = 541.0f;
//	float sdd = 949.0f;
//
//	float x0(0.0f);
//	float y0(sid);
//	float z0(0.0f);
//
//	int DNU(888);
//	int DNV(64);
//	int PN(984 * 4);
//	float imgXCenter(0);
//	float imgYCenter(0);
//	float imgZCenter(0);
//
//	int XN(512);
//	int YN(512);
//	int ZN(64);
//	float dx(500.0f / 512.0f);
//	float dz(0.625);
//
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//	//Generate the positions of the detectors
//	float col_size = 1.0239 ;
//	float row_size = 1.0963 ;
//	float col_offset = 0;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//	imgXCenter = 0;
//	imgYCenter = 0;
//	imgZCenter = 0;
//
//	//Test different view number with HD geometry
//
//	float* hvol = new float[XN * YN * ZN];
//	for(int i = 0 ;i != XN * YN * ZN; ++i)
//	{
//		hvol[i] = 1.0;
//	}
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != XN; ++i)
//	{
//		for(int j = 0; j != YN; ++j)
//		{
//			if(sqrt(pow((i-XN/2+0.5)/(XN/2),2.0)+ pow((j-YN/2+0.5)/(YN/2),2.0))<0.98)
//			{
//				mask[i * YN + j] = 1;
//			}
//			else
//			{
//				mask[i * YN + j] = 0;
//			}
//
//		}
//
//	}
//
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(PN);
//		hzPos[ii] = 0;// (ii - PN / 2) * 0.0015;
//	}
//	float* hprj = new float[DNU * DNV * PN];
//
//	//Proj with DD
//	clock_t start;
//	clock_t stop;
//	start = clock();
//
//	DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, hvol, hprj, dx, dz, mask, 0, 0, 0);
//	stop = clock();
//	std::cout<<"Projection With Two GPUs spends "<<static_cast<double>(stop - start) / static_cast<double>(CLOCKS_PER_SEC)<<" seconds\n";
//
//}
//
//void TestProjectionTimeWithNumOfGPUs(int SIZE1,int SIZE2)
//{
//	//Set the HD Geometry
//	// Set the geometry
//	float sid = 541.0f;
//	float sdd = 949.0f;
//
//	float x0(0.0f);
//	float y0(sid);
//	float z0(0.0f);
//
//	int DNU(888);
//	int DNV(64);
//	int PN(984 * 4);
//	float imgXCenter(0);
//	float imgYCenter(0);
//	float imgZCenter(0);
//
//	int XN(512);
//	int YN(512);
//	int ZN(64);
//	float dx(500.0f / 512.0f);
//	float dz(0.625);
//
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//	//Generate the positions of the detectors
//	float col_size = 1.0239 ;
//	float row_size = 1.0963 ;
//	float col_offset = 0;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//	imgXCenter = 0;
//	imgYCenter = 0;
//	imgZCenter = 0;
//
//	//Test different view number with HD geometry
//
//
//
//	float* hvol = new float[XN * YN * ZN];
//	for(int i = 0 ;i != XN * YN * ZN; ++i)
//	{
//		hvol[i] = 1.0;
//	}
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != XN; ++i)
//	{
//		for(int j = 0; j != YN; ++j)
//		{
//			if(sqrt(pow((i-XN/2+0.5)/(XN/2),2.0)+ pow((j-YN/2+0.5)/(YN/2),2.0))<0.98)
//			{
//				mask[i * YN + j] = 1;
//			}
//			else
//			{
//				mask[i * YN + j] = 0;
//			}
//
//		}
//
//	}
//
//
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(PN);
//		hzPos[ii] = 0;// (ii - PN / 2) * 0.0015;
//	}
//	float* hprj = new float[DNU * DNV * PN];
//
//	//Proj with DD
//	clock_t start;
//	clock_t stop;
//	start = clock();
//	DD3Back_3gpus_fortable(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, hvol, hprj, dx, dz, mask, 0, 0, 0,SIZE1, SIZE2);
//	stop = clock();
//	std::cout<<"Projection With Two GPUs spends "<<static_cast<double>(stop - start) / static_cast<double>(CLOCKS_PER_SEC)<<" seconds\n";
//
//}
//
//
//
//class ImgRecon
//{
//private:
//	float sid;
//	float sdd;
//
//	float x0; // initial source position
//	float y0;
//	float z0;
//
//	int DNU;
//	int DNV;
//	int PN;
//
//	//Generate the positions of the detectors
//	float col_size; // Unit: mm;
//	float row_size;
//
//	float col_offset; // Unit: one detector cell
//	float row_offset;
//
//	//Information of the volume
//	float imgXCenter; //Center of the image position
//	float imgYCenter;
//	float imgZCenter;
//	float dx;
//	float dz;
//
//	int XN;
//	int YN;
//	int ZN;
//
//
//	thrust::host_vector<float> xds;
//	thrust::host_vector<float> yds;
//	thrust::host_vector<float> zds;
//
//	thrust::host_vector<float> hangs;
//	thrust::host_vector<float> hzpos;
//
//	thrust::host_vector<byte> mask;
//
//	thrust::host_vector<float> hprj;
//	thrust::host_vector<float> hvol;
//
//public:
//	ImgRecon():sid(538.5200193125),sdd(946.745971679688),x0(0),y0(538.5200193125),z0(0),
//	DNU(888),DNV(64),PN(984),col_size(1.02390003204346),row_size(1.0963),
//	col_offset(-1.25),row_offset(0),imgXCenter(0),imgYCenter(0),imgZCenter(0),
//	dx(0.9765625),dz(0.625),XN(512),YN(512),ZN(64)
//    {
//		xds.resize(DNU);
//		yds.resize(DNU);
//		zds.resize(DNV);
//
//		hangs.resize(PN);
//		hzpos.resize(PN);
//
//		mask.resize(XN*YN);
//
//		hprj.resize(DNU * DNV * PN);
//		hvol.resize(XN * YN * ZN);
//
//		float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//		float curBeta = 0;
//		for (int ii = 0; ii != DNU; ++ii)
//		{
//			curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//			xds[ii] = sinf(curBeta) * sdd;
//			yds[ii] = sid - cosf(curBeta) * sdd;
//		}
//
//		for (int ii = 0; ii != DNV; ++ii)
//		{
//			zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//		}
//
//		for (int ii = 0; ii != PN; ++ii)
//		{
//			hangs[ii] = ii * TWOPI / static_cast<float>(984);
//			hzpos[ii] = 0; //Cone beam CT
//		}
//
//		for (int i = 0; i != YN; ++i)
//		{
//			for (int j = 0; j != XN; ++j)
//			{
//				if (sqrtf(powf((i - YN / 2.0 + 0.5) / (YN / 2.0), 2.0) + powf((j - XN / 2.0 + 0.5) / (XN / 2.0), 2.0)) < 0.98)
//				{
//					mask[i * XN + j] = 1;
//				}
//				else
//				{
//					mask[i * XN + j] = 0;
//				}
//			}
//		}
//
//    }
//
//
//
//};
//
//
//
//void CGreconBigPatient20151124()
//{
//	//Define the parameters
//	float sid = 538.5200193125;
//	float sdd = 946.745971679688;
//
//	float x0 = 0.0f;
//	float y0 = sid;
//	float z0 = 0.0f;
//
//	int DNU = 888;
//	int DNV = 64;
//	int PN = 8139;
//
//	float imgXCenter = 0.0;
//	float imgYCenter = 0.0;
//	float imgZCenter = 0.0;
//
//	int XN = 600;
//	int YN = 600;
//	int ZN = 64 * 10;
//
//	float dx(500.0f / 512.0f);
//	float dz = 0.625;
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//
//	//Generate the positions of the detectors
//	float col_size = 1.02390003204346;
//	float row_size = 1.0963;
//
//	float col_offset = -1.25;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(984) -39.0939140319824 / 180.0 * 3.14159265358 + 3.14159265258;
//		hzPos[ii] = (static_cast<float>(ii)-4070.0) / 984.0 * 0.625 * 63.0;
//	}
//
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != YN; ++i)
//	{
//		for (int j = 0; j != XN; ++j)
//		{
//			if (sqrtf(powf((i - YN / 2.0 + 0.5) / (YN / 2.0), 2.0) + powf((j - XN / 2.0 + 0.5) / (XN / 2.0), 2.0)) < 0.98)
//			{
//				mask[i * XN + j] = 1;
//			}
//			else
//			{
//				mask[i * XN + j] = 0;
//			}
//		}
//	}
//	float* hprj = new float[DNU * DNV * PN];
//
//	//Load the projection
//	std::ifstream fin("BigPatientProj.prj",std::ios::binary);
//	if(!fin.is_open())
//	{
//		std::cerr<<"Cannot open the file\n";
//		exit(-1);
//	}
//	fin.read((char*)hprj,sizeof(float)* DNU * DNV * PN);
//	fin.close();
//
//
//
//	//Calculate B
//	thrust::host_vector<float> X(XN * YN * ZN, 0);
//	thrust::host_vector<float> B(XN * YN * ZN, 0);
//	thrust::host_vector<float> R(XN * YN * ZN, 0);
//	thrust::host_vector<float> nextR(XN * YN * ZN, 0);
//
//	thrust::host_vector<float> Delta(XN * YN * ZN, 0);
//	thrust::host_vector<float> TDelta(DNU * DNV *PN,0);
//	thrust::host_vector<float> MDelta(XN *YN*ZN, 0);
//
//	DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//				imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//				XN, YN, ZN, &B[0], hprj, dx, dz, mask, 0, 0, 0);
//
//
//	R = B;
//	Delta = R;
//	double alpha,beta;
//	double alpha1,alpha2;
//	double beta1, beta2;
//	clock_t start = clock();
//	for(int i = 0; i != 2; ++i)
//	{
//		//1 Calculate the scalar alpha
//		// (1) calculate MDelta
//		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, &Delta[0], &TDelta[0], dx, dz, mask, 0, 0);
//		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//					imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//					XN, YN, ZN, &MDelta[0], &TDelta[0], dx, dz, mask, 0, 0,0);
//		// (2)
//		alpha1 = thrust::inner_product(R.begin(),R.end(),R.begin(),0.0);
//		alpha2 = thrust::inner_product(MDelta.begin(),MDelta.end(),Delta.begin(),0.0);
//		std::cout<<"alpha1 = "<<alpha1<<std::endl;
//		std::cout<<"alpha2 = "<<alpha2<<std::endl;
//		alpha = alpha1 / alpha2;
//		thrust::transform(X.begin(),X.end(),Delta.begin(),X.begin(),CGop<float>(alpha));
//		thrust::transform(R.begin(),R.end(),MDelta.begin(),nextR.begin(),CGop<float>(-alpha));
//		beta1 = thrust::inner_product(nextR.begin(),nextR.end(),nextR.begin(),0.0);
//		beta2 = thrust::inner_product(R.begin(),R.end(),R.begin(),0.0);
//		beta = beta1 / beta2;
//		thrust::transform(nextR.begin(),nextR.end(),Delta.begin(),Delta.begin(),CGop<float>(beta));
//		R = nextR;
//		std::cout<<i<<std::endl;
//
//	}
//	clock_t end = clock();
//	std::cout<<"Total Time is "<<(static_cast<double>(end) - static_cast<double>(start)) / static_cast<double>(CLOCKS_PER_SEC)<<" seconds\n";
//	std::string name = "CGreconBigPatient.raw";
//	std::ofstream fou3(name.c_str(), std::ios::binary);
//	fou3.write((char*) &X[0], sizeof(float) * XN * YN * ZN);
//	fou3.close();
//}
//
//
//
//void reconBigPatient20151121()
//{
//	//Define the parameters
//	float sid = 538.5200193125;
//	float sdd = 946.745971679688;
//
//	float x0 = 0.0f;
//	float y0 = sid;
//	float z0 = 0.0f;
//
//	int DNU = 888;
//	int DNV = 64;
//	int PN = 8139;
//
//	float imgXCenter = 0.0;
//	float imgYCenter = 0.0;
//	float imgZCenter = 0.0;
//
//	int XN = 700;
//	int YN = 700;
//	int ZN = 64 * 12;
//
//	float dx(500.0f / 512.0f);
//	float dz = 0.625;
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//
//	//Generate the positions of the detectors
//	float col_size = 1.02390003204346;
//	float row_size = 1.0963;
//
//	float col_offset = -1.25;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(984) -39.0939140319824 / 180.0 * 3.14159265358 + 3.14159265258;
//		hzPos[ii] = (static_cast<float>(ii)-4070.0) / 984.0 * 0.625 * 63.0;
//	}
//
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != YN; ++i)
//	{
//		for (int j = 0; j != XN; ++j)
//		{
//			if (sqrtf(powf((i - YN / 2.0 + 0.5) / (YN / 2.0), 2.0) + powf((j - XN / 2.0 + 0.5) / (XN / 2.0), 2.0)) < 0.98)
//			{
//				mask[i * XN + j] = 1;
//			}
//			else
//			{
//				mask[i * XN + j] = 0;
//			}
//		}
//	}
//	float* hprj = new float[DNU * DNV * PN];
//	float* hvol = new float[XN * YN * ZN];
//
//
//	//Load the projection
//	std::ifstream fin("BigPatientProj.prj",std::ios::binary);
//	if(!fin.is_open())
//	{
//		std::cerr<<"Cannot open the file\n";
//		exit(-1);
//	}
//	fin.read((char*)hprj,sizeof(float)* DNU * DNV * PN);
//	fin.close();
//
//
//
//	float* oneVol = new float[XN * YN * ZN];
//	float* onePrj = new float[DNU * DNV * PN];
//	float* rowSum = new float[DNU * DNV * PN];
//	float* colSum = new float[XN * YN * ZN];
//	float* lasImg = new float[XN * YN * ZN];
//
//	for (int ii = 0; ii != XN * YN * ZN; ++ii)
//	{
//		hvol[ii] = 0;
//		oneVol[ii] = 1.0;
//		colSum[ii] = 0;
//	}
//	for (int ii = 0; ii != DNU * DNV * PN; ++ii)
//	{
//		onePrj[ii] = 1.0;
//		rowSum[ii] = 0.0;
//	}
//
//	DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, oneVol, rowSum, dx, dz, mask, 0, 0);
//	std::ofstream fou1("rowSum.raw",std::ios::binary);
//	fou1.write((char*)rowSum,sizeof(float) * DNU * DNV * PN);
//	fou1.close();
//
//	DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, colSum, onePrj, dx, dz, mask, 0, 0, 0);
//	std::ofstream fou2("colSum.arw",std::ios::binary);
//	fou2.write((char*)colSum,sizeof(float)* XN*YN*ZN);
//	fou2.close();
//	float t1 = 1.0;
//	float t2 = 1.0;
//	int totN = XN * YN * ZN;
//	clock_t start = clock();
//	for (int ii = 0; ii != 800; ++ii)
//	{
//		memcpy(lasImg, hvol, sizeof(float) * totN);
//
//		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, hvol, onePrj, dx, dz, mask, 0, 0);
//
//		prjWeight(onePrj, hprj, rowSum, DNU*DNV*PN);
//
//		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, oneVol, onePrj, dx, dz, mask, 0, 0, 0);
//
//		bakWeight(oneVol, hvol, colSum, totN);
//
//		t2 = (1.0 + sqrtf(1.0 + 4.0 * t1 * t1)) / 2.0;
//		FISTA(lasImg, hvol, t1, t2, totN);
//		t1 = t2;
//
//
//		if ((ii + 1) % 10 == 0)
//		{
//			std::stringstream ss;
//			ss << ii;
//			std::string name = "BigPatientRecon" + ss.str() + ".raw";
//			std::ofstream fou3(name.c_str(), std::ios::binary);
//			fou3.write((char*) hvol, sizeof(float) * XN * YN * ZN);
//			fou3.close();
//		}
//		std::cout << "Big patient recon iteration # = "<<ii<<std::endl;
//	}
//	clock_t end = clock();
//	std::cout<<"The total time is "<<static_cast<double>(end - start)/static_cast<double>(CLOCKS_PER_SEC);
//}
//
//
//
//void reconBigPatient20151126_alreadyinGPU()
//{
//	//Define the parameters
//	float sid = 538.5200193125;
//	float sdd = 946.745971679688;
//
//	float x0 = 0.0f;
//	float y0 = sid;
//	float z0 = 0.0f;
//
//	int DNU = 888;
//	int DNV = 64;
//	int PN = 8139;
//
//	float imgXCenter = 0.0;
//	float imgYCenter = 0.0;
//	float imgZCenter = 0.0;
//
//	int XN = 700;
//	int YN = 700;
//	int ZN = 64 * 12;
//
//	float dx(500.0f / 512.0f);
//	float dz = 0.625;
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//
//	//Generate the positions of the detectors
//	float col_size = 1.02390003204346;
//	float row_size = 1.0963;
//
//	float col_offset = -1.25;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(984) -39.0939140319824 / 180.0 * 3.14159265358 + 3.14159265258;
//		hzPos[ii] = (static_cast<float>(ii)-4070.0) / 984.0 * 0.625 * 63.0;
//	}
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != YN; ++i)
//	{
//		for (int j = 0; j != XN; ++j)
//		{
//			if (sqrtf(powf((i - YN / 2.0 + 0.5) / (YN / 2.0), 2.0) + powf((j - XN / 2.0 + 0.5) / (XN / 2.0), 2.0)) < 0.98)
//			{
//				mask[i * XN + j] = 1;
//			}
//			else
//			{
//				mask[i * XN + j] = 0;
//			}
//		}
//	}
//
//
//	float* hprj = new float[DNU * DNV * PN];
//	float* hvol = new float[XN * YN * ZN];
//
//
//	//Load the projection
//	std::ifstream fin("BigPatientProj.prj",std::ios::binary);
//	if(!fin.is_open())
//	{
//		std::cerr<<"Cannot open the file\n";
//		exit(-1);
//	}
//	fin.read((char*)hprj,sizeof(float)* DNU * DNV * PN);
//	fin.close();
//
//
//	float* oneVol = new float[XN * YN * ZN];
//	float* onePrj = new float[DNU * DNV * PN];
//	//float* rowSum = new float[DNU * DNV * PN];
//	//float* colSum = new float[XN * YN * ZN];
//	//float* lasImg = new float[XN * YN * ZN];
//
//	for (int ii = 0; ii != XN * YN * ZN; ++ii)
//	{
//		//hvol[ii] = 0;
//		oneVol[ii] = 1.0;
//		//colSum[ii] = 0;
//	}
//	for (int ii = 0; ii != DNU * DNV * PN; ++ii)
//	{
//		onePrj[ii] = 1.0;
//		//rowSum[ii] = 0.0;
//	}
//
//	thrust::device_vector<float> d_xds(xds,xds+DNU);
//	thrust::device_vector<float> d_yds(yds,yds+DNU);
//	thrust::device_vector<float> d_zds(zds,zds+DNV);
//	thrust::device_vector<float> d_angs(hangs,hangs+PN);
//	thrust::device_vector<float> d_zPos(hzPos,hzPos+PN);
//	thrust::device_vector<byte> d_mask(mask,mask+XN*YN);
//
//	thrust::device_vector<float> d_oneVol(oneVol,oneVol + XN * YN * ZN);
//	thrust::device_vector<float> d_rowSum(DNU*DNV*PN,0);
//	DD3Proj_gpu_alreadyinGPU(x0, y0, z0, DNU, DNV, d_xds, d_yds, d_zds,
//		imgXCenter, imgYCenter, imgZCenter, d_angs, d_zPos, PN,
//		XN, YN, ZN, d_oneVol, d_rowSum, dx, dz, d_mask, 0, 0);
//	thrust::host_vector<float> h_oneVol = d_oneVol;
//	thrust::host_vector<float> h_rowSum = d_rowSum;
//	d_oneVol.clear();
//	d_rowSum.clear();
//
//	thrust::device_vector<float> d_onePrj(onePrj,onePrj+DNU*DNV*PN);
//	thrust::device_vector<float> d_colSum(XN *YN*ZN,0);
//	DD3Back_gpu_alreadyinGPU(x0, y0, z0, DNU, DNV, d_xds, d_yds, d_zds,
//		imgXCenter, imgYCenter, imgZCenter, d_angs, d_zPos, PN,
//		XN, YN, ZN, d_colSum, d_onePrj, dx, dz, d_mask, 0, 0, 0);
//	thrust::host_vector<float> h_colSum = d_colSum;
//	thrust::host_vector<float> h_onePrj = d_onePrj;
//	d_colSum.clear();
//	d_onePrj.clear();
//
//
//	thrust::device_vector<float> d_prj(hprj,hprj+DNU*DNV*PN);
//	thrust::device_vector<float> d_vol(XN*YN*ZN,0);
//	thrust::device_vector<float> d_lasImg(XN*YN*ZN,0);
//
////
//
////	float t1 = 1.0;
////	float t2 = 1.0;
////	int totN = XN * YN * ZN;
////	clock_t start = clock();
////	for (int ii = 0; ii != 800; ++ii)
////	{
////		d_lasImg = d_vol;
////
////		//memcpy(lasImg, hvol, sizeof(float) * totN);
////
////		DD3Proj_gpu_alreadyinGPU(x0, y0, z0, DNU, DNV, d_xds, d_yds, d_zds,
////			imgXCenter, imgYCenter, imgZCenter, d_angs, d_zPos, PN,
////			XN, YN, ZN, d_vol, d_onePrj, dx, dz, d_mask, 0, 0);
////
////		prjWeight_v2<float>(d_onePrj, d_prj, d_rowSum);
////
////		DD3Back_gpu_alreadyinGPU(x0, y0, z0, DNU, DNV, d_xds, d_yds, d_zds,
////			imgXCenter, imgYCenter, imgZCenter, d_angs, d_zPos, PN,
////			XN, YN, ZN, d_oneVol, d_onePrj, dx, dz, d_mask, 0, 0, 0);
////
////		bakWeight_v2<float>(d_oneVol, d_vol, d_colSum);
////
////		t2 = (1.0 + sqrtf(1.0 + 4.0 * t1 * t1)) / 2.0;
////		FISTA<float>(d_lasImg, d_vol, t1, t2);
////		t1 = t2;
////
////
////		if ((ii + 1) % 10 == 0)
////		{
////			std::stringstream ss;
////			ss << ii;
////			thrust::host_vector<float> hvol = d_vol;
////
////			std::string name = "BigPatientReconAlreadyinGPU" + ss.str() + ".raw";
////			std::ofstream fou3(name.c_str(), std::ios::binary);
////			fou3.write((char*) &hvol[0], sizeof(float) * XN * YN * ZN);
////			fou3.close();
////			hvol.clear();
////
////		}
////		std::cout << "Big patient recon iteration # = "<<ii<<std::endl;
////	}
////	clock_t end = clock();
////	std::cout<<"The total time is "<<static_cast<double>(end - start)/static_cast<double>(CLOCKS_PER_SEC);
//}
//
//
//
//void realPatientrecon20151126_alreadyinGPU()
//{
//	//Define the parameters
//	float sid = 538.5200193125;
//	float sdd = 946.745971679688;
//
//	float x0 = 0.0f;
//	float y0 = sid;
//	float z0 = 0.0f;
//
//	int DNU = 888;
//	int DNV = 64;
//	int PN = 2201;
//
//	float imgXCenter = 0.0;
//	float imgYCenter = 0.0;
//	float imgZCenter = 0.0;
//
//	int XN = 512;
//	int YN = 512;
//	int ZN = 64;
//
//	float dx(500.0f / 512.0f);
//	float dz = 0.625;
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//
//	//Generate the positions of the detectors
//	float col_size = 1.02390003204346;
//	float row_size = 1.0963;
//
//	float col_offset = -1.25;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(PN);// -39.0939140319824 / 180.0 * 3.14159265358 + 3.14159265258;
//		hzPos[ii] = 0;//(static_cast<float>(ii)) / 984.0 * 0.625 * 63.0;
//	}
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != YN; ++i)
//	{
//		for (int j = 0; j != XN; ++j)
//		{
//			if (sqrtf(powf((i - YN / 2.0 + 0.5) / (YN / 2.0), 2.0) + powf((j - XN / 2.0 + 0.5) / (XN / 2.0), 2.0)) < 0.98)
//			{
//				mask[i * XN + j] = 1;
//			}
//			else
//			{
//				mask[i * XN + j] = 0;
//			}
//		}
//	}
//
//
//	float* hprj = new float[DNU * DNV * PN];
//	float* hvol = new float[XN * YN * ZN];
//
//
//	//Load the projection
//	std::ifstream fin("patientReal.raw",std::ios::binary);
//	if(!fin.is_open())
//	{
//		std::cerr<<"Cannot open the file\n";
//		exit(-1);
//	}
//	fin.read((char*)hvol,sizeof(float)* XN*YN*ZN);
//	fin.close();
//
//
//	DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, hvol, hprj, dx, dz, mask, 0, 0);
//
//
//	float* oneVol = new float[XN * YN * ZN];
//	float* onePrj = new float[DNU * DNV * PN];
//
//
//	for (int ii = 0; ii != XN * YN * ZN; ++ii)
//	{
//		//hvol[ii] = 0;
//		oneVol[ii] = 1.0;
//		//colSum[ii] = 0;
//	}
//	for (int ii = 0; ii != DNU * DNV * PN; ++ii)
//	{
//		onePrj[ii] = 1.0;
//		//rowSum[ii] = 0.0;
//	}
//
//	thrust::device_vector<float> d_xds(xds,xds+DNU);
//	thrust::device_vector<float> d_yds(yds,yds+DNU);
//	thrust::device_vector<float> d_zds(zds,zds+DNV);
//	thrust::device_vector<float> d_angs(hangs,hangs+PN);
//	thrust::device_vector<float> d_zPos(hzPos,hzPos+PN);
//	thrust::device_vector<byte> d_mask(mask,mask+XN*YN);
//
//	thrust::device_vector<float> d_oneVol(oneVol,oneVol + XN * YN * ZN);
//	thrust::device_vector<float> d_rowSum(DNU*DNV*PN,0);
//
//	thrust::device_vector<float> d_onePrj(onePrj,onePrj+DNU*DNV*PN);
//	thrust::device_vector<float> d_colSum(XN *YN*ZN,0);
//
//	thrust::device_vector<float> d_prj(hprj,hprj+DNU*DNV*PN);
//	thrust::device_vector<float> d_vol(XN*YN*ZN,0);
//	thrust::device_vector<float> d_lasImg(XN*YN*ZN,0);
//
//
//	DD3Proj_gpu_alreadyinGPU(x0, y0, z0, DNU, DNV, d_xds, d_yds, d_zds,
//		imgXCenter, imgYCenter, imgZCenter, d_angs, d_zPos, PN,
//		XN, YN, ZN, d_oneVol, d_rowSum, dx, dz, d_mask, 0, 0);
//	DD3Back_gpu_alreadyinGPU(x0, y0, z0, DNU, DNV, d_xds, d_yds, d_zds,
//		imgXCenter, imgYCenter, imgZCenter, d_angs, d_zPos, PN,
//		XN, YN, ZN, d_colSum, d_onePrj, dx, dz, d_mask, 0, 0, 0);
//
//
//
//
//
//	float t1 = 1.0;
//	float t2 = 1.0;
//	//int totN = XN * YN * ZN;
//	clock_t start = clock();
//	for (int ii = 0; ii != 800; ++ii)
//	{
//		d_lasImg = d_vol;
//
//		//memcpy(lasImg, hvol, sizeof(float) * totN);
//
//		DD3Proj_gpu_alreadyinGPU(x0, y0, z0, DNU, DNV, d_xds, d_yds, d_zds,
//			imgXCenter, imgYCenter, imgZCenter, d_angs, d_zPos, PN,
//			XN, YN, ZN, d_vol, d_onePrj, dx, dz, d_mask, 0, 0);
//
//		prjWeight_v2<float>(d_onePrj, d_prj, d_rowSum);
//
//		DD3Back_gpu_alreadyinGPU(x0, y0, z0, DNU, DNV, d_xds, d_yds, d_zds,
//			imgXCenter, imgYCenter, imgZCenter, d_angs, d_zPos, PN,
//			XN, YN, ZN, d_oneVol, d_onePrj, dx, dz, d_mask, 0, 0, 0);
//
//		bakWeight_v2<float>(d_oneVol, d_vol, d_colSum);
//
//		t2 = (1.0 + sqrtf(1.0 + 4.0 * t1 * t1)) / 2.0;
//		FISTA<float>(d_lasImg, d_vol, t1, t2, -2000, 4000);
//		t1 = t2;
//
//
//		if ((ii + 1) % 10 == 0)
//		{
//			std::stringstream ss;
//			ss << ii;
//			thrust::host_vector<float> hvol = d_vol;
//
//			std::string name = "realPatientReconAlreadyinGPU" + ss.str() + ".raw";
//			std::ofstream fou3(name.c_str(), std::ios::binary);
//			fou3.write((char*) &hvol[0], sizeof(float) * XN * YN * ZN);
//			fou3.close();
//			hvol.clear();
//
//
//		}
//		std::cout << "Iteration # = "<<ii<<std::endl;
//	}
//	clock_t end = clock();
//	std::cout<<"The total time is "<<static_cast<double>(end - start)/static_cast<double>(CLOCKS_PER_SEC);
//}
//
//
//// With OS-technique
//void reconBigPatient20151123()
//{
//	//Define the parameters
//	float sid = 538.5200193125;
//	float sdd = 946.745971679688;
//
//	float x0 = 0.0f;
//	float y0 = sid;
//	float z0 = 0.0f;
//
//	int DNU = 888;
//	int DNV = 64;
//	int PN = 8139;
//
//	float imgXCenter = 0.0;
//	float imgYCenter = 0.0;
//	float imgZCenter = 0.0;
//
//	int XN = 700;
//	int YN = 700;
//	int ZN = 64 * 12;
//
//	float dx(500.0f / 512.0f);
//	float dz = 0.625;
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//
//	//Generate the positions of the detectors
//	float col_size = 1.02390003204346;
//	float row_size = 1.0963;
//
//	float col_offset = -1.25;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(984) -39.0939140319824 / 180.0 * 3.14159265358 + 3.14159265258;
//		hzPos[ii] = (static_cast<float>(ii)-4070.0) / 984.0 * 0.625 * 63.0;
//	}
//
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != YN; ++i)
//	{
//		for (int j = 0; j != XN; ++j)
//		{
//			if (sqrtf(powf((i - YN / 2.0 + 0.5) / (YN / 2.0), 2.0) + powf((j - XN / 2.0 + 0.5) / (XN / 2.0), 2.0)) < 0.98)
//			{
//				mask[i * XN + j] = 1;
//			}
//			else
//			{
//				mask[i * XN + j] = 0;
//			}
//		}
//	}
//	float* hprj = new float[DNU * DNV * PN];
//	float* hvol = new float[XN * YN * ZN];
//
//
//	//Load the projection
//	std::ifstream fin("BigPatientProj.prj",std::ios::binary);
//	if(!fin.is_open())
//	{
//		std::cerr<<"Cannot open the file\n";
//		exit(-1);
//	}
//	fin.read((char*)hprj,sizeof(float)* DNU * DNV * PN);
//	fin.close();
//
//	float* oneVol = new float[XN * YN * ZN];
//	float* onePrj = new float[DNU * DNV * PN];
//	float* rowSum = new float[DNU * DNV * PN];
//	float* colSum = new float[XN * YN * ZN];
//	float* lasImg = new float[XN * YN * ZN];
//
//	for (int ii = 0; ii != XN * YN * ZN; ++ii)
//	{
//		hvol[ii] = 0;
//		oneVol[ii] = 1.0;
//		colSum[ii] = 0;
//	}
//	for (int ii = 0; ii != DNU * DNV * PN; ++ii)
//	{
//		onePrj[ii] = 1.0;
//		rowSum[ii] = 0.0;
//	}
//
//	int osNum = 20;
//	int prjPerOS = PN / osNum;
//	int* SPN = new int[osNum];
//	float** shangs = new float*[osNum];
//	float** shzPos = new float*[osNum];
//	float** shprj = new float*[osNum];
//	float** shrowSum = new float*[osNum];
//	for(int i = 0; i != osNum-1; ++i)
//	{
//		SPN[i] = prjPerOS;
//		shangs[i] = new float[SPN[i]];
//		shzPos[i] = new float[SPN[i]];
//		shprj[i] = new float[DNU * DNV * SPN[i]];
//		shrowSum[i] = new float[DNU * DNV * SPN[i]];
//
//		memcpy(shangs[i],hangs + i * prjPerOS, sizeof(float) * prjPerOS);
//		memcpy(shzPos[i],hzPos + i * prjPerOS, sizeof(float) * prjPerOS);
//		memcpy(shprj[i], hprj + i * prjPerOS * DNU * DNV, sizeof(float) * prjPerOS * DNU * DNV);
//
//	}
//	SPN[osNum-1] = PN - prjPerOS * (osNum-1);
//	shangs[osNum-1] = new float[SPN[osNum-1]];
//	shzPos[osNum-1] = new float[SPN[osNum-1]];
//	shprj[osNum-1] = new float[DNU * DNV * SPN[osNum-1]];
//	shrowSum[osNum-1] = new float[DNU * DNV * SPN[osNum-1]];
//
//	memcpy(shangs[osNum-1], hangs + (osNum-1) * prjPerOS, sizeof(float) * (PN - prjPerOS * osNum));
//	memcpy(shzPos[osNum-1], hzPos + (osNum-1) * prjPerOS, sizeof(float) * (PN - prjPerOS * osNum));
//	memcpy(shprj[osNum-1], hprj + (osNum-1) * prjPerOS * DNU * DNV, sizeof(float) * (PN - prjPerOS * osNum) * DNU * DNV);
//
//	for(int i = 0; i != osNum; ++i)
//	{
//		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, shangs[i], shzPos[i], SPN[i],
//			XN, YN, ZN, oneVol, shrowSum[i], dx, dz, mask, 0, 0);
//		std::cout<<"os "<<i<<std::endl;
//	}
////
////
////
//
////
////
////	DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
////		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
////		XN, YN, ZN, colSum, onePrj, dx, dz, mask, 0, 0, 0);
////
////	float t1 = 1.0;
////	float t2 = 1.0;
////	int totN = XN * YN * ZN;
////	clock_t start = clock();
////	for (int ii = 0; ii != 800; ++ii)
////	{
////		memcpy(lasImg, hvol, sizeof(float) * totN);
////
////		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
////			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
////			XN, YN, ZN, hvol, onePrj, dx, dz, mask, 0, 0);
////
////		prjWeight(onePrj, hprj, rowSum, DNU*DNV*PN);
////
////		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
////			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
////			XN, YN, ZN, oneVol, onePrj, dx, dz, mask, 0, 0, 0);
////
////		bakWeight(oneVol, hvol, colSum, totN);
////
////		t2 = (1.0 + sqrtf(1.0 + 4.0 * t1 * t1)) / 2.0;
////		FISTA(lasImg, hvol, t1, t2, totN);
////		t1 = t2;
////
////
////		if ((ii + 1) % 10 == 0)
////		{
////			std::stringstream ss;
////			ss << ii;
////			std::string name = "BigPatientRecon" + ss.str() + ".raw";
////			std::ofstream fou3(name.c_str(), std::ios::binary);
////			fou3.write((char*) hvol, sizeof(float) * XN * YN * ZN);
////			fou3.close();
////		}
////		std::cout << "Big patient recon iteration # = "<<ii<<std::endl;
////	}
////	clock_t end = clock();
////	std::cout<<"The total time is "<<static_cast<double>(end - start)/static_cast<double>(CLOCKS_PER_SEC);
//}
//
//
//
//
//
//// With OS-technique
//void reconBigPatient20151123_v2()
//{
//	//Define the parameters
//	float sid = 538.5200193125;
//	float sdd = 946.745971679688;
//
//	float x0 = 0.0f;
//	float y0 = sid;
//	float z0 = 0.0f;
//
//	int DNU = 888;
//	int DNV = 64;
//	int PN = 8139;
//
//	float imgXCenter = 0.0;
//	float imgYCenter = 0.0;
//	float imgZCenter = 0.0;
//
//	int XN = 512;
//	int YN = 512;
//	int ZN = 512;
//
//	float dx(600.0f / 512.0f);
//	float dz = 600.0/512.0f;
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//
//	//Generate the positions of the detectors
//	float col_size = 1.02390003204346;
//	float row_size = 1.0963;
//
//	float col_offset = -1.25;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(984) - 3.14159265358979 * 8139.0 / 984.0;// -39.0939140319824 / 180.0 * 3.14159265358 + 3.14159265258;
//		hzPos[ii] = (static_cast<float>(ii)-4070.0) / 984.0 * 0.625 * 63.0;
//	}
//
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != YN; ++i)
//	{
//		for (int j = 0; j != XN; ++j)
//		{
//			if (sqrtf(powf((i - YN / 2.0 + 0.5) / (YN / 2.0), 2.0) + powf((j - XN / 2.0 + 0.5) / (XN / 2.0), 2.0)) < 0.98)
//			{
//				mask[i * XN + j] = 1;
//			}
//			else
//			{
//				mask[i * XN + j] = 0;
//			}
//		}
//	}
//
//
//	float* hprj = new float[DNU * DNV * PN];
//	float* hvol = new float[XN * YN * ZN];
//
//
//	//Load the projection
//	std::ifstream fin("BigPatientProj.prj",std::ios::binary);
//	if(!fin.is_open())
//	{
//		std::cerr<<"Cannot open the file\n";
//		exit(-1);
//	}
//	fin.read((char*)hprj,sizeof(float)* DNU * DNV * PN);
//	fin.close();
//
//	int osNum = 20; // Divide the projection views into 20 subsets
//	int viewPersubSet = PN / osNum; // How many views are in one subset
//	int res = PN - viewPersubSet * osNum; //How many views are more in the last subset
//	int viewInLastSubSet = res + viewPersubSet;
//
//	int* SPN = new int[osNum];
//	float** shangs = new float*[osNum];
//	float** shzPos = new float*[osNum];
//	float** shprj = new float*[osNum];
//	float** srowSum = new float*[osNum];
//	float* oneVol = new float[XN * YN * ZN];
//	//float* colSum = new float[XN * YN * ZN];
//	float** onePrj = new float*[osNum];
//	float** colSum = new float*[osNum];
//
//	for(int i = 0; i != XN * YN * ZN; ++i)
//	{
//		oneVol[i] = 1.0;
//		//colSum[i] = 0;
//	}
//	for(int i = 0 ; i != osNum-1; ++i)
//	{
//		SPN[i] = viewPersubSet;
//	}
//	SPN[osNum - 1] = viewInLastSubSet;
//
//	for(int i = 0; i != osNum; ++i)
//	{
//		shangs[i] = new float[SPN[i]];
//		shzPos[i] = new float[SPN[i]];
//		shprj[i] = new float[SPN[i] * DNU * DNV];
//		srowSum[i] = new float[SPN[i] * DNU * DNV];
//		onePrj[i] = new float[SPN[i] * DNU * DNV];
//		colSum[i] = new	float[XN * YN * ZN];
//		for(int curIdx = 0; curIdx != viewPersubSet; ++curIdx)
//		{
//			shangs[i][curIdx] = hangs[curIdx * osNum + i];
//			shzPos[i][curIdx] = hzPos[curIdx * osNum + i];
//			memcpy(shprj[i] + curIdx * DNU * DNV, hprj + (curIdx * osNum + i) * DNU * DNV, sizeof(float) * DNU * DNV);
//
//		}
//
//	}
//	memcpy(shangs[osNum-1] + viewPersubSet, hangs + viewPersubSet * osNum, sizeof(float) * res);
//	memcpy(shzPos[osNum-1] + viewPersubSet, hzPos + viewPersubSet * osNum, sizeof(float) * res);
//	memcpy(shprj[osNum-1] + viewPersubSet * DNU * DNV, hprj + viewPersubSet * osNum * DNU * DNV, sizeof(float) * DNU * DNV * res);
//
//	float* allOneproj = new float[DNU * DNV * PN];
//	for(int i =0; i != DNU * DNV * PN; ++i)
//	{
//		allOneproj[i] = 1.0;
//	}
//	for(int i = 0; i != osNum; ++i)
//	{
//		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, shangs[i], shzPos[i], SPN[i],
//			XN, YN, ZN, oneVol, srowSum[i], dx, dz, mask, 0, 0);
//		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, shangs[i], shzPos[i], SPN[i],
//			XN, YN, ZN, colSum[i], allOneproj, dx, dz, mask, 0, 0, 0);
//
//		std::cout<<"os "<<i<<std::endl;
//	}
//
//	float* lasImg = new float[XN * YN * ZN];
//	float* tempVol = new float[XN * YN * ZN];
//
//	//float* colSum = new float[XN * YN * ZN];
//
//
//	float t1 = 1.0;
//	float t2 = 1.0;
//	int totN = XN * YN * ZN;
//	clock_t start = clock();
//	for (int ii = 0; ii != 30; ++ii)
//	{
//		memcpy(lasImg, hvol, sizeof(float) * totN);
//		for(int subIdx = 0; subIdx != osNum;++subIdx)
//		{
//			DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//				imgXCenter, imgYCenter, imgZCenter, shangs[subIdx], shzPos[subIdx], SPN[subIdx],
//				XN, YN, ZN, hvol, onePrj[subIdx], dx, dz, mask, 0, 0);
//			prjWeight(onePrj[subIdx], shprj[subIdx], srowSum[subIdx], DNU*DNV*SPN[subIdx]);
//
//
//			DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//				imgXCenter, imgYCenter, imgZCenter, shangs[subIdx], shzPos[subIdx], SPN[subIdx],
//				XN, YN, ZN, oneVol, onePrj[subIdx], dx, dz, mask, 0, 0, 0);
//			bakWeight(oneVol, hvol, colSum[subIdx], totN);
//			std::cout<<"OS is "<<subIdx<<std::endl;
//
//		}
//
//		t2 = (1.0 + sqrtf(1.0 + 4.0 * t1 * t1)) / 2.0;
//		FISTA(lasImg, hvol, t1, t2, totN);
//		t1 = t2;
//
//		std::stringstream ss;
//		ss<<ii;
//		std::string name = "BigPatientReconWithOSSART20151218_" + ss.str() + ".raw";
//		std::ofstream fu(name.c_str(), std::ios::binary);
//		fu.write((char*) hvol, sizeof(float) * XN * YN * ZN);
//		fu.close();
//
//		std::cout << "Big patient recon iteration # = "<<ii<<std::endl;
//	}
//	clock_t end = clock();
//	std::cout<<"The total time is "<<static_cast<double>(end - start)/static_cast<double>(CLOCKS_PER_SEC);
//
//
//}
//
//struct Recon
//{
//	float sid;
//	float sdd;
//	float x0;
//	float y0;
//	float z0;
//	int DNU;
//	int DNV;
//	int PN;
//	float imgXCenter;
//	float imgYCenter;
//	float imgZCenter;
//	int XN;
//	int YN;
//	int ZN;
//	float dx;
//	float dz;
//	thrust::host_vector<float> x_ds;
//	thrust::host_vector<float> y_ds;
//	thrust::host_vector<float> z_ds;
//	float col_size;
//	float row_size;
//	float col_offset;
//	float row_offset;
//	thrust::host_vector<float> angs;
//	thrust::host_vector<float> zPos;
//	thrust::host_vector<byte> msk;
//	void SART_recon(std::string ProjDataFileName, std::string ReconFileName, int osNum, int totalIterNum);
//	void CG_recon(std::string ProjDataFileName, std::string ReconFileName, int totalIterNum);
//
//};
//
//void Recon::SART_recon(std::string ProjDataFileName, std::string ReconFileName, int osNum, int totalIterNum)
//{
//	float* xds = &x_ds[0];
//	float* yds = &y_ds[0];
//	float* zds = &z_ds[0];
//
//	float* hangs = &angs[0];
//	float* hzPos = &zPos[0];
//	byte* mask = &msk[0];
//
//	float* hprj = new float[DNU * DNV * PN];
//	float* hvol = new float[XN * YN * ZN];
//
//	//Load the projection
//	std::ifstream fin(ProjDataFileName.c_str(),std::ios::binary);
//	if(!fin.is_open())
//	{
//		std::cerr<<"Cannot open the file\n";
//		exit(-1);
//	}
//	fin.read((char*)hprj,sizeof(float)* DNU * DNV * PN);
//	fin.close();
//
//	int viewPersubSet = PN / osNum; // How many views are in one subset
//	int res = PN - viewPersubSet * osNum; //How many views are more in the last subset
//	int viewInLastSubSet = res + viewPersubSet;
//
//	int* SPN = new int[osNum];
//	float** shangs = new float*[osNum];
//	float** shzPos = new float*[osNum];
//	float** shprj = new float*[osNum];
//	float** srowSum = new float*[osNum];
//	float* oneVol = new float[XN * YN * ZN];
//	//float* colSum = new float[XN * YN * ZN];
//	float** onePrj = new float*[osNum];
//	float** colSum = new float*[osNum];
//
//	for(int i = 0; i != XN * YN * ZN; ++i)
//	{
//		oneVol[i] = 1.0;
//		//colSum[i] = 0;
//	}
//	for(int i = 0 ; i != osNum-1; ++i)
//	{
//		SPN[i] = viewPersubSet;
//	}
//	SPN[osNum - 1] = viewInLastSubSet;
//
//	for(int i = 0; i != osNum; ++i)
//	{
//		shangs[i] = new float[SPN[i]];
//		shzPos[i] = new float[SPN[i]];
//		shprj[i] = new float[SPN[i] * DNU * DNV];
//		srowSum[i] = new float[SPN[i] * DNU * DNV];
//		onePrj[i] = new float[SPN[i] * DNU * DNV];
//		colSum[i] = new	float[XN * YN * ZN];
//		for(int curIdx = 0; curIdx != viewPersubSet; ++curIdx)
//		{
//			shangs[i][curIdx] = hangs[curIdx * osNum + i];
//			shzPos[i][curIdx] = hzPos[curIdx * osNum + i];
//			memcpy(shprj[i] + curIdx * DNU * DNV, hprj + (curIdx * osNum + i) * DNU * DNV, sizeof(float) * DNU * DNV);
//
//		}
//
//	}
//
//	memcpy(shangs[osNum-1] + viewPersubSet, hangs + viewPersubSet * osNum, sizeof(float) * res);
//	memcpy(shzPos[osNum-1] + viewPersubSet, hzPos + viewPersubSet * osNum, sizeof(float) * res);
//	memcpy(shprj[osNum-1] + viewPersubSet * DNU * DNV, hprj + viewPersubSet * osNum * DNU * DNV, sizeof(float) * DNU * DNV * res);
//
//	float* allOneproj = new float[DNU * DNV * PN];
//	for(int i =0; i != DNU * DNV * PN; ++i)
//	{
//		allOneproj[i] = 1.0;
//	}
//	for(int i = 0; i != osNum; ++i)
//	{
//		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, shangs[i], shzPos[i], SPN[i],
//			XN, YN, ZN, oneVol, srowSum[i], dx, dz, mask, 0, 0);
//		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, shangs[i], shzPos[i], SPN[i],
//			XN, YN, ZN, colSum[i], allOneproj, dx, dz, mask, 0, 0, 0);
//
//		std::cout<<"os "<<i<<std::endl;
//	}
//
//	float* lasImg = new float[XN * YN * ZN];
//	float* tempVol = new float[XN * YN * ZN];
//
//	float t1 = 1.0;
//	float t2 = 1.0;
//	int totN = XN * YN * ZN;
//	clock_t start = clock();
//	for (int ii = 0; ii != totalIterNum; ++ii)
//	{
//		memcpy(lasImg, hvol, sizeof(float) * totN);
//		for(int subIdx = 0; subIdx != osNum;++subIdx)
//		{
//			DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//				imgXCenter, imgYCenter, imgZCenter, shangs[subIdx], shzPos[subIdx], SPN[subIdx],
//				XN, YN, ZN, hvol, onePrj[subIdx], dx, dz, mask, 0, 0);
//			prjWeight(onePrj[subIdx], shprj[subIdx], srowSum[subIdx], DNU*DNV*SPN[subIdx]);
//
//
//			DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//				imgXCenter, imgYCenter, imgZCenter, shangs[subIdx], shzPos[subIdx], SPN[subIdx],
//				XN, YN, ZN, oneVol, onePrj[subIdx], dx, dz, mask, 0, 0, 0);
//
//			bakWeight(oneVol, hvol, colSum[subIdx], totN);
//			std::cout<<"OS is "<<subIdx<<std::endl;
//
//		}
//
//		t2 = (1.0 + sqrtf(1.0 + 4.0 * t1 * t1)) / 2.0;
//		FISTA(lasImg, hvol, t1, t2, totN);
//		t1 = t2;
//
//		std::cout << "iteration # = "<<ii<<std::endl;
//	}
//	clock_t end = clock();
//
//	std::ofstream fou3(ReconFileName.c_str(), std::ios::binary);
//	fou3.write((char*) hvol, sizeof(float) * XN * YN * ZN);
//	fou3.close();
//
//	std::cout<<"The total time is "<<static_cast<double>(end - start)/static_cast<double>(CLOCKS_PER_SEC);
//	delete[] SPN;
//	for(int i = 0; i != osNum;++i)
//	{
//		delete[] shangs[i];
//		delete[] shzPos[i];
//		delete[] shprj[i];
//		delete[] srowSum[i];
//		delete[] onePrj[i];
//		delete[] colSum[i];
//	}
//	delete[] shangs;
//	delete[] shzPos;
//	delete[] shprj;
//	delete[] srowSum;
//	delete[] onePrj;
//	delete[] colSum;
//	delete[] allOneproj;
//	delete[] lasImg;
//	delete[] tempVol;
//	delete[] hprj;
//	delete[] hvol;
//
//}
//
//void SART_recon(float sid, float sdd, float x0, float y0, float z0, int DNU, int DNV, int PN,
//		float imgXCenter, float imgYCenter, float imgZCenter, int XN, int YN, int ZN, float dx, float dz,
//		float* xds, float* yds, float* zds, float col_size, float row_size,
//		float col_offset, float row_offset, float* hangs, float* hzPos, byte* mask,
//		const std::string ProjDataFileName, const std::string ReconFileName,
//		int osNum, int totalIterNum)
//{
//	float* hprj = new float[DNU * DNV * PN];
//	float* hvol = new float[XN * YN * ZN];
//
//	//Load the projection
//	std::ifstream fin(ProjDataFileName.c_str(),std::ios::binary);
//	if(!fin.is_open())
//	{
//		std::cerr<<"Cannot open the file\n";
//		exit(-1);
//	}
//	fin.read((char*)hprj,sizeof(float)* DNU * DNV * PN);
//	fin.close();
//
//	int viewPersubSet = PN / osNum; // How many views are in one subset
//	int res = PN - viewPersubSet * osNum; //How many views are more in the last subset
//	int viewInLastSubSet = res + viewPersubSet;
//
//	int* SPN = new int[osNum];
//	float** shangs = new float*[osNum];
//	float** shzPos = new float*[osNum];
//	float** shprj = new float*[osNum];
//	float** srowSum = new float*[osNum];
//	float* oneVol = new float[XN * YN * ZN];
//	//float* colSum = new float[XN * YN * ZN];
//	float** onePrj = new float*[osNum];
//	float** colSum = new float*[osNum];
//
//	for(int i = 0; i != XN * YN * ZN; ++i)
//	{
//		oneVol[i] = 1.0;
//		//colSum[i] = 0;
//	}
//	for(int i = 0 ; i != osNum-1; ++i)
//	{
//		SPN[i] = viewPersubSet;
//	}
//	SPN[osNum - 1] = viewInLastSubSet;
//
//	for(int i = 0; i != osNum; ++i)
//	{
//		shangs[i] = new float[SPN[i]];
//		shzPos[i] = new float[SPN[i]];
//		shprj[i] = new float[SPN[i] * DNU * DNV];
//		srowSum[i] = new float[SPN[i] * DNU * DNV];
//		onePrj[i] = new float[SPN[i] * DNU * DNV];
//		colSum[i] = new	float[XN * YN * ZN];
//		for(int curIdx = 0; curIdx != viewPersubSet; ++curIdx)
//		{
//			shangs[i][curIdx] = hangs[curIdx * osNum + i];
//			shzPos[i][curIdx] = hzPos[curIdx * osNum + i];
//			memcpy(shprj[i] + curIdx * DNU * DNV, hprj + (curIdx * osNum + i) * DNU * DNV, sizeof(float) * DNU * DNV);
//
//		}
//
//	}
//
//	memcpy(shangs[osNum-1] + viewPersubSet, hangs + viewPersubSet * osNum, sizeof(float) * res);
//	memcpy(shzPos[osNum-1] + viewPersubSet, hzPos + viewPersubSet * osNum, sizeof(float) * res);
//	memcpy(shprj[osNum-1] + viewPersubSet * DNU * DNV, hprj + viewPersubSet * osNum * DNU * DNV, sizeof(float) * DNU * DNV * res);
//
//	float* allOneproj = new float[DNU * DNV * PN];
//	for(int i =0; i != DNU * DNV * PN; ++i)
//	{
//		allOneproj[i] = 1.0;
//	}
//	for(int i = 0; i != osNum; ++i)
//	{
//		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, shangs[i], shzPos[i], SPN[i],
//			XN, YN, ZN, oneVol, srowSum[i], dx, dz, mask, 0, 0);
//		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, shangs[i], shzPos[i], SPN[i],
//			XN, YN, ZN, colSum[i], allOneproj, dx, dz, mask, 0, 0, 0);
//
//		std::cout<<"os "<<i<<std::endl;
//	}
//
//	float* lasImg = new float[XN * YN * ZN];
//	float* tempVol = new float[XN * YN * ZN];
//
//	float t1 = 1.0;
//	float t2 = 1.0;
//	int totN = XN * YN * ZN;
//	clock_t start = clock();
//	for (int ii = 0; ii != totalIterNum; ++ii)
//	{
//		memcpy(lasImg, hvol, sizeof(float) * totN);
//		for(int subIdx = 0; subIdx != osNum;++subIdx)
//		{
//			DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//				imgXCenter, imgYCenter, imgZCenter, shangs[subIdx], shzPos[subIdx], SPN[subIdx],
//				XN, YN, ZN, hvol, onePrj[subIdx], dx, dz, mask, 0, 0);
//			prjWeight(onePrj[subIdx], shprj[subIdx], srowSum[subIdx], DNU*DNV*SPN[subIdx]);
//
//
//			DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//				imgXCenter, imgYCenter, imgZCenter, shangs[subIdx], shzPos[subIdx], SPN[subIdx],
//				XN, YN, ZN, oneVol, onePrj[subIdx], dx, dz, mask, 0, 0, 0);
//
//			bakWeight(oneVol, hvol, colSum[subIdx], totN);
//			std::cout<<"OS is "<<subIdx<<std::endl;
//
//		}
//
//		t2 = (1.0 + sqrtf(1.0 + 4.0 * t1 * t1)) / 2.0;
//		FISTA(lasImg, hvol, t1, t2, totN);
//		t1 = t2;
//
//		std::cout << "iteration # = "<<ii<<std::endl;
//	}
//	clock_t end = clock();
//
//	std::ofstream fou3(ReconFileName.c_str(), std::ios::binary);
//	fou3.write((char*) hvol, sizeof(float) * XN * YN * ZN);
//	fou3.close();
//
//	std::cout<<"The total time is "<<static_cast<double>(end - start)/static_cast<double>(CLOCKS_PER_SEC);
//	delete[] SPN;
//	for(int i = 0; i != osNum;++i)
//	{
//		delete[] shangs[i];
//		delete[] shzPos[i];
//		delete[] shprj[i];
//		delete[] srowSum[i];
//		delete[] onePrj[i];
//		delete[] colSum[i];
//	}
//	delete[] shangs;
//	delete[] shzPos;
//	delete[] shprj;
//	delete[] srowSum;
//	delete[] onePrj;
//	delete[] colSum;
//	delete[] allOneproj;
//	delete[] lasImg;
//	delete[] tempVol;
//	delete[] hprj;
//	delete[] hvol;
//
//}
//
//
//void reconBigPatientFunction(const std::string& ProjectionName,
//		const std::string& VolumeName,
//		const int osNum = 20,
//		const int iterNum = 30,
//		const bool outputMedRes = true,
//		const float sid=538.5200193125, const float sdd=946.745971679699,
//		const int DNU = 888, const int DNV = 64, const int PN = 8139,
//		const float imgXCenter = 0.0f, const float imgYCenter = 0.0f, const float imgZCenter = 0.0f,
//		const int XN = 512, const int YN = 512, const int ZN = 512,
//		const float dx = 1.171875, const float dz = 1.171875,
//		const float col_size = 1.0239f, const float row_size = 1.0963f,
//		const float col_offset = -1.25f, const float row_offset = 0.0f,
//		const float start_view = 0, const int view_per_rot = 984, const float PITCH = 63.0)
//{
//	//float PITCH = 63.0;
//	const float x0 = 0.0f;
//	const float y0 = sid;
//	const float z0 = 0.0f;
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//
//	//Generate the positions of the detectors
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = start_view + ii * TWOPI / static_cast<float>(view_per_rot);
//		hzPos[ii] = (static_cast<float>(ii) - static_cast<float>(PN) * 0.5) / static_cast<float>(view_per_rot) * 0.625 * PITCH;
//	}
//
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != YN; ++i)
//	{
//		for (int j = 0; j != XN; ++j)
//		{
//			if (sqrtf(powf((i - YN / 2.0 + 0.5) / (YN / 2.0), 2.0) + powf((j - XN / 2.0 + 0.5) / (XN / 2.0), 2.0)) < 0.98)
//			{
//				mask[i * XN + j] = 1;
//			}
//			else
//			{
//				mask[i * XN + j] = 0;
//			}
//		}
//	}
//
//
//	float* hprj = new float[DNU * DNV * PN];
//	float* hvol = new float[XN * YN * ZN];
//
//
//	//Load the projection
//	std::ifstream fin(ProjectionName.c_str(),std::ios::binary);
//	if(!fin.is_open())
//	{
//		std::cerr<<"Cannot open the file\n";
//		exit(-1);
//	}
//	fin.read((char*)hprj,sizeof(float)* DNU * DNV * PN);
//	fin.close();
//
//	//int osNum = 20; // Divide the projection views into 20 subsets
//	int viewPersubSet = PN / osNum; // How many views are in one subset
//	int res = PN - viewPersubSet * osNum; //How many views are more in the last subset
//	int viewInLastSubSet = res + viewPersubSet;
//
//	int* SPN = new int[osNum];
//	float** shangs = new float*[osNum];
//	float** shzPos = new float*[osNum];
//	float** shprj = new float*[osNum];
//	float** srowSum = new float*[osNum];
//	float* oneVol = new float[XN * YN * ZN];
//	//float* colSum = new float[XN * YN * ZN];
//	float** onePrj = new float*[osNum];
//	float** colSum = new float*[osNum];
//
//	for(int i = 0; i != XN * YN * ZN; ++i)
//	{
//		oneVol[i] = 1.0;
//		//colSum[i] = 0;
//	}
//	for(int i = 0 ; i != osNum-1; ++i)
//	{
//		SPN[i] = viewPersubSet;
//	}
//	SPN[osNum - 1] = viewInLastSubSet;
//
//	for(int i = 0; i != osNum; ++i)
//	{
//		shangs[i] = new float[SPN[i]];
//		shzPos[i] = new float[SPN[i]];
//		shprj[i] = new float[SPN[i] * DNU * DNV];
//		srowSum[i] = new float[SPN[i] * DNU * DNV];
//		onePrj[i] = new float[SPN[i] * DNU * DNV];
//		colSum[i] = new	float[XN * YN * ZN];
//		for(int curIdx = 0; curIdx != viewPersubSet; ++curIdx)
//		{
//			shangs[i][curIdx] = hangs[curIdx * osNum + i];
//			shzPos[i][curIdx] = hzPos[curIdx * osNum + i];
//			memcpy(shprj[i] + curIdx * DNU * DNV, hprj + (curIdx * osNum + i) * DNU * DNV, sizeof(float) * DNU * DNV);
//
//		}
//
//	}
//	memcpy(shangs[osNum-1] + viewPersubSet, hangs + viewPersubSet * osNum, sizeof(float) * res);
//	memcpy(shzPos[osNum-1] + viewPersubSet, hzPos + viewPersubSet * osNum, sizeof(float) * res);
//	memcpy(shprj[osNum-1] + viewPersubSet * DNU * DNV, hprj + viewPersubSet * osNum * DNU * DNV, sizeof(float) * DNU * DNV * res);
//
//	float* allOneproj = new float[DNU * DNV * PN];
//	for(int i =0; i != DNU * DNV * PN; ++i)
//	{
//		allOneproj[i] = 1.0;
//	}
//	for(int i = 0; i != osNum; ++i)
//	{
//		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, shangs[i], shzPos[i], SPN[i],
//			XN, YN, ZN, oneVol, srowSum[i], dx, dz, mask, 0, 0);
//		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, shangs[i], shzPos[i], SPN[i],
//			XN, YN, ZN, colSum[i], allOneproj, dx, dz, mask, 0, 0, 0);
//
//		std::cout<<"os "<<i<<std::endl;
//	}
//
//	float* lasImg = new float[XN * YN * ZN];
//	float* tempVol = new float[XN * YN * ZN];
//
//	//float* colSum = new float[XN * YN * ZN];
//
//
//	float t1 = 1.0;
//	float t2 = 1.0;
//	int totN = XN * YN * ZN;
//	clock_t start = clock();
//	for (int ii = 0; ii != iterNum; ++ii)
//	{
//		memcpy(lasImg, hvol, sizeof(float) * totN);
//		for(int subIdx = 0; subIdx != osNum;++subIdx)
//		{
//			DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//				imgXCenter, imgYCenter, imgZCenter, shangs[subIdx], shzPos[subIdx], SPN[subIdx],
//				XN, YN, ZN, hvol, onePrj[subIdx], dx, dz, mask, 0, 0);
//			prjWeight(onePrj[subIdx], shprj[subIdx], srowSum[subIdx], DNU*DNV*SPN[subIdx]);
//
//
//			DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//				imgXCenter, imgYCenter, imgZCenter, shangs[subIdx], shzPos[subIdx], SPN[subIdx],
//				XN, YN, ZN, oneVol, onePrj[subIdx], dx, dz, mask, 0, 0, 0);
//			bakWeight(oneVol, hvol, colSum[subIdx], totN);
//			std::cout<<"OS is "<<subIdx<<std::endl;
//
//		}
//
//		t2 = (1.0 + sqrtf(1.0 + 4.0 * t1 * t1)) / 2.0;
//		FISTA(lasImg, hvol, t1, t2, totN);
//		t1 = t2;
//
//		if(outputMedRes)
//		{
//			std::stringstream ss;
//			ss<<ii;
//			std::string name = VolumeName + "_" + ss.str() + ".raw";
//			std::ofstream fu(name.c_str(), std::ios::binary);
//			fu.write((char*) hvol, sizeof(float) * XN * YN * ZN);
//			fu.close();
//		}
//		std::cout << "Iteration # = "<<ii<<std::endl;
//	}
//	clock_t end = clock();
//
//	std::string name = VolumeName + ".raw";
//	std::ofstream fu(name.c_str(), std::ios::binary);
//	fu.write((char*) hvol, sizeof(float) * XN * YN * ZN);
//	fu.close();
//
//	std::cout<<"The total time is "<<static_cast<double>(end - start)/static_cast<double>(CLOCKS_PER_SEC);
//}
//
//
//#define EPSILON 1.0e-9
//
//template<typename T>
//__global__ void _GradientOfTV_ker(T* f, T* d, const T coef, const int  L, const int  W, const int H)
//{
//	const int idx = threadIdx.x + blockDim.x * blockIdx.x;
//	const int idy = threadIdx.y + blockDim.y * blockIdx.y;
//	const int idz = threadIdx.z + blockDim.z * blockIdx.z;
//	if (idx < L && idy < W && idz < H)
//	{
//		const int  curId = (idz * W + idy) * L + idx;
//		const T f_ijk = f[(idz * W + idy) * L + idx];
//		const T f_i1jk = f[(idz * W + idy) * L + ((idx + 1) % L)];
//		const T f_ij1k = f[(idz * W + ((idy + 1) % W)) * L + idx];
//		const T f_ijk1 = f[(((idz + 1) % H) * W + idy) * L + idx];
//		const T f_i_1jk = f[(idz * W + idy) * L + ((idx + L - 1) % L)];
//		const T f_i_1j1k = f[(idz * W + ((idy + 1) % W)) * L + ((idx + L - 1) % L)];
//		const T f_i_1jk1 = f[(((idz + 1) % H) * W + idy) * L + ((idx + L - 1) % L)];
//		const T f_ij_1k = f[(idz * W + ((idy + W - 1) % W)) * L + idx];
//		const T f_i1j_1k = f[(idz * W + ((idy + W - 1) % W)) * L + ((idx + 1) % L)];
//		const T f_ij_1k1 = f[(((idz + 1) % H) * W + ((idy + W - 1) % W)) * L + idx];
//		const T f_ijk_1 = f[(((idz + H - 1) % H) * W + idy) * L + idx];
//		const T f_i1jk_1 = f[(((idz + H - 1) % H) * W + idy) * L + ((idx + 1) % L)];
//		const T f_ij1k_1 = f[(((idz + H - 1) % H) * W + ((idy + 1) % W)) * L + idx];
//
//		const T dom1 = 1.0 / sqrt((f_i1jk - f_ijk)*(f_i1jk - f_ijk) + (f_ij1k - f_ijk)*(f_ij1k - f_ijk) + (f_ijk1 - f_ijk)*(f_ijk1 - f_ijk) + EPSILON);
//		const T dom2 = 1.0 / sqrt((f_ijk - f_i_1jk)*(f_ijk - f_i_1jk) + (f_i_1j1k - f_i_1jk)*(f_i_1j1k - f_i_1jk) + (f_i_1jk1 - f_i_1jk)*(f_i_1jk1 - f_i_1jk) + EPSILON);
//		const T dom3 = 1.0 / sqrt((f_i1j_1k - f_ij_1k)*(f_i1j_1k - f_ij_1k) + (f_ijk - f_ij_1k)*(f_ijk - f_ij_1k) + (f_ij_1k1 - f_ij_1k)*(f_ij_1k1 - f_ij_1k) + EPSILON);
//		const T dom4 = 1.0 / sqrt((f_i1jk_1 - f_ijk_1)*(f_i1jk_1 - f_ijk_1) + (f_ij1k_1 - f_ijk_1)*(f_ij1k_1 - f_ijk_1) + (f_ijk - f_ijk_1)*(f_ijk - f_ijk_1) + EPSILON);
//
//		d[curId] = ((3.0 * f_ijk - f_i1jk - f_ij1k - f_ijk1) * dom1 + (f_ijk - f_i_1jk) * dom2 + (f_ijk - f_ij_1k) * dom3 + (f_ijk - f_ijk_1) * dom4) * coef;
//		return;
//	}
//}
//
//
//
//
//
//void GradientOfTV_GPU(float* f, float* d, const float coef, const int & L, const int & W, const int & H, const dim3& blk, const dim3& gid)
//{
//	_GradientOfTV_ker<float> << <gid, blk >> >(f, d, coef, L, W, H);
//}
////
////template<typename T>
////struct UpdateGongHao
////{
////	T coef;
////	UpdateGongHao(const T& c):coef(c){}
////	__host__ __device__ T operator()(const T& f, const T& d)
////	{
////		return f - d * coef;
////	}
////};
//
//__global__ void updateGongHao(float* f, float* d, float coef, int LEN)
//{
//	int i = threadIdx.x + blockIdx.x * blockDim.x;
//	if(i < LEN)
//	{
//		f[i] = f[i] - d[i] * coef;
//	}
//}
//
//void gonghaoTV(float* cpuF, const float coef, const int iterNum, const int L, const int W, const int H)
//{
//	int ll = L * W * H;
//	int S = sizeof(float) * L * W * H;
//	float* gpuF;
//	float* gpuD;
//	cudaMalloc((void**)&gpuF, S);
//	cudaMalloc((void**)&gpuD, S);
//	cudaMemcpy(gpuF, cpuF, S,cudaMemcpyHostToDevice);
//	cudaMemset(gpuD, 0, S);
//	dim3 blk(8,8,8);
//	dim3 gid(
//			(L+blk.x-1)/blk.x,
//			(W+blk.y-1)/blk.y,
//			(H+blk.z-1)/blk.z);
//	dim3 blk2(1024);
//	dim3 gid2((ll + blk2.x - 1) / blk2.x);
//
//	for (int i = 0; i != iterNum; ++i)
//	{
//		GradientOfTV_GPU(gpuF, gpuD, 1.0, L, W,  H, blk, gid);
//		updateGongHao<<<gid2,blk2>>>(gpuF, gpuD, coef, ll);
//		//thrust::transform(gpuF,gpuF + ll, gpuD, gpuF, UpdateGongHao<float>(coef));
//		cudaMemset(gpuD, 0, S);
//	}
//	cudaMemcpy(cpuF,gpuF, S, cudaMemcpyDeviceToHost);
//	cudaFree(gpuF);
//	cudaFree(gpuD);
//
//}
//
//// The detector is a panel detector
//void reconGongHaoDataFunction(const std::string& ProjectionName,
//		const std::string& VolumeName,
//		const int osNum = 20,
//		const int iterNum = 30,
//		const bool outputMedRes = true,
//		const float sid = 522.0373, const float sdd = 674.6248,
//		const int DNU = 725, const int DNV = 253, const int PN = 720,
//		const float imgXCenter = 0.0f, const float imgYCenter = 0.0f, const float imgZCenter = 0.0f,
//		const int XN = 512, const int YN = 512, const int ZN = 512,
//		const float dx = 0.3095, const float dz = 0.3095,
//		const float col_size = 0.2, const float row_size = 0.2,
//		const float col_offset = 1.3477, const float row_offset = 0.0f,
//		const float start_view = 0, const int view_per_rot = 720, const float PITCH = 0)
//{
//
//	const float x0 = 0.0f;
//	const float y0 = sid;
//	const float z0 = 0.0f;
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//
//	//Generate the positions of the detectors
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		xds[ii] = (ii - (DNU - 1.0) * 0.5 + col_offset) * col_size;
//		yds[ii] = sid - sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = start_view + ii * TWOPI / static_cast<float>(view_per_rot);
//		hzPos[ii] = (static_cast<float>(ii) - static_cast<float>(PN) * 0.5) / static_cast<float>(view_per_rot) * 0.625 * PITCH;
//	}
//
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != YN; ++i)
//	{
//		for (int j = 0; j != XN; ++j)
//		{
//			if (sqrtf(powf((i - YN / 2.0 + 0.5) / (YN / 2.0), 2.0) + powf((j - XN / 2.0 + 0.5) / (XN / 2.0), 2.0)) < 0.98)
//			{
//				mask[i * XN + j] = 1;
//			}
//			else
//			{
//				mask[i * XN + j] = 0;
//			}
//		}
//	}
//
//
//	float* hprj = new float[DNU * DNV * PN];
//	float* hvol = new float[XN * YN * ZN];
//
//
//	//Load the projection
//	std::ifstream fin(ProjectionName.c_str(),std::ios::binary);
//	if(!fin.is_open())
//	{
//		std::cerr<<"Cannot open the file\n";
//		exit(-1);
//	}
//	fin.read((char*)hprj,sizeof(float)* DNU * DNV * PN);
//	fin.close();
//
//	int viewPersubSet = PN / osNum; // How many views are in one subset
//	int res = PN - viewPersubSet * osNum; //How many views are more in the last subset
//	int viewInLastSubSet = res + viewPersubSet;
//
//	int* SPN = new int[osNum];
//	float** shangs = new float*[osNum];
//	float** shzPos = new float*[osNum];
//	float** shprj = new float*[osNum];
//	float** srowSum = new float*[osNum];
//	float* oneVol = new float[XN * YN * ZN];
//	//float* colSum = new float[XN * YN * ZN];
//	float** onePrj = new float*[osNum];
//	float** colSum = new float*[osNum];
//
//	for(int i = 0; i != XN * YN * ZN; ++i)
//	{
//		oneVol[i] = 1.0;
//		//colSum[i] = 0;
//	}
//	for(int i = 0 ; i != osNum-1; ++i)
//	{
//		SPN[i] = viewPersubSet;
//	}
//	SPN[osNum - 1] = viewInLastSubSet;
//
//	for(int i = 0; i != osNum; ++i)
//	{
//		shangs[i] = new float[SPN[i]];
//		shzPos[i] = new float[SPN[i]];
//		shprj[i] = new float[SPN[i] * DNU * DNV];
//		srowSum[i] = new float[SPN[i] * DNU * DNV];
//		onePrj[i] = new float[SPN[i] * DNU * DNV];
//		colSum[i] = new	float[XN * YN * ZN];
//		for(int curIdx = 0; curIdx != viewPersubSet; ++curIdx)
//		{
//			shangs[i][curIdx] = hangs[curIdx * osNum + i];
//			shzPos[i][curIdx] = hzPos[curIdx * osNum + i];
//			memcpy(shprj[i] + curIdx * DNU * DNV, hprj + (curIdx * osNum + i) * DNU * DNV, sizeof(float) * DNU * DNV);
//
//		}
//
//	}
//	memcpy(shangs[osNum-1] + viewPersubSet, hangs + viewPersubSet * osNum, sizeof(float) * res);
//	memcpy(shzPos[osNum-1] + viewPersubSet, hzPos + viewPersubSet * osNum, sizeof(float) * res);
//	memcpy(shprj[osNum-1] + viewPersubSet * DNU * DNV, hprj + viewPersubSet * osNum * DNU * DNV, sizeof(float) * DNU * DNV * res);
//
//	float* allOneproj = new float[DNU * DNV * PN];
//	for(int i =0; i != DNU * DNV * PN; ++i)
//	{
//		allOneproj[i] = 1.0;
//	}
//	for(int i = 0; i != osNum; ++i)
//	{
//		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, shangs[i], shzPos[i], SPN[i],
//			XN, YN, ZN, oneVol, srowSum[i], dx, dz, mask, 0, 0);
//
//		DD3_panel_gpu_back_branchless_sat2d(x0, y0, z0,
//		    DNU, DNV,
//		    xds, yds, zds,
//		    imgXCenter, imgYCenter, imgZCenter,
//		    shangs[i], shzPos[i], SPN[i],
//		    XN,  YN,  ZN,
//		    colSum[i], allOneproj, dx, dz, mask,0,0);
//
//		std::cout<<"os "<<i<<std::endl;
//	}
//
//	float* lasImg = new float[XN * YN * ZN];
//	float* tempVol = new float[XN * YN * ZN];
//
//	float t1 = 1.0;
//	float t2 = 1.0;
//	int totN = XN * YN * ZN;
//	clock_t start = clock();
//	for (int ii = 0; ii != iterNum; ++ii)
//	{
//		memcpy(lasImg, hvol, sizeof(float) * totN);
//		for(int subIdx = 0; subIdx != osNum;++subIdx)
//		{
//			DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//				imgXCenter, imgYCenter, imgZCenter, shangs[subIdx], shzPos[subIdx], SPN[subIdx],
//				XN, YN, ZN, hvol, onePrj[subIdx], dx, dz, mask, 0, 0);
//			prjWeight(onePrj[subIdx], shprj[subIdx], srowSum[subIdx], DNU*DNV*SPN[subIdx]);
//
//
//			DD3_panel_gpu_back_branchless_sat2d(x0, y0, z0, DNU, DNV, xds, yds, zds,
//				imgXCenter, imgYCenter, imgZCenter, shangs[subIdx], shzPos[subIdx], SPN[subIdx],
//				XN, YN, ZN, oneVol, onePrj[subIdx], dx, dz, mask, 0, 0);
//
//			bakWeight(oneVol, hvol, colSum[subIdx], totN);
//			std::cout<<"OS is "<<subIdx<<std::endl;
//
//		}
//
//		t2 = (1.0 + sqrtf(1.0 + 4.0 * t1 * t1)) / 2.0;
//		FISTA(lasImg, hvol, t1, t2, totN);
//		t1 = t2;
//
//		// TV minimizations
//		if(ii > 60)
//		{
//			gonghaoTV(hvol, 0.0004 * pow(0.999, ii), 5, ZN,XN,YN);
//		}
//
//
//		if(outputMedRes)
//		{
//			std::stringstream ss;
//			ss<<ii;
//			std::string name = VolumeName + "_" + ss.str() + ".raw";
//			std::ofstream fu(name.c_str(), std::ios::binary);
//			fu.write((char*) hvol, sizeof(float) * XN * YN * ZN);
//			fu.close();
//		}
//		std::cout << "Iteration # = "<<ii<<std::endl;
//	}
//	clock_t end = clock();
//
//	std::string name = VolumeName + ".raw";
//	std::ofstream fu(name.c_str(), std::ios::binary);
//	fu.write((char*) hvol, sizeof(float) * XN * YN * ZN);
//	fu.close();
//
//	std::cout<<"The total time is "<<static_cast<double>(end - start)/static_cast<double>(CLOCKS_PER_SEC);
//}
//
//
//void reconBigPatientFunctionCG(const std::string& ProjectionName,
//		const std::string& VolumeName,
//		const int iterNum = 30,
//		const bool outputMedRes = false,
//		const float sid=538.5200193125, const float sdd=946.745971679699,
//		const int DNU = 888, const int DNV = 64, const int PN = 8139,
//		const float imgXCenter = 0.0f, const float imgYCenter = 0.0f, const float imgZCenter = 0.0f,
//		const int XN = 512, const int YN = 512, const int ZN = 512,
//		const float dx = 1.171875, const float dz = 1.171875,
//		const float col_size = 1.0239f, const float row_size = 1.0963f,
//		const float col_offset = -1.25f, const float row_offset = 0.0f,
//		const float start_view = 0, const int view_per_rot = 984,
//		const float PITCH = 63.0)
//{
//	//Define the parameters
//	float x0 = 0.0f;
//	float y0 = sid;
//	float z0 = 0.0f;
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//
//	//Generate the positions of the detectors
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = start_view + ii * TWOPI / static_cast<float>(view_per_rot) * 6;
//		hzPos[ii] = (static_cast<float>(ii)-static_cast<float>(PN) * 0.5)
//				/ static_cast<float>(view_per_rot) * 0.625 * PITCH * 6;
//	}
//
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != YN; ++i)
//	{
//		for (int j = 0; j != XN; ++j)
//		{
//			if (sqrtf(powf((i - YN / 2.0 + 0.5) / (YN / 2.0), 2.0) + powf((j - XN / 2.0 + 0.5) / (XN / 2.0), 2.0)) < 0.98)
//			{
//				mask[i * XN + j] = 1.0;
//			}
//			else
//			{
//				mask[i * XN + j] = 1.0;
//			}
//		}
//	}
//	float* hprj = new float[DNU * DNV * PN];
//	//float* hvol = new float[XN * YN * ZN];
//
//
//	//Load the projection
//	std::ifstream fin(ProjectionName.c_str(),std::ios::binary);
//	if(!fin.is_open())
//	{
//		std::cerr<<"Cannot open the file\n";
//		exit(-1);
//	}
//	fin.read((char*)hprj,sizeof(float)* DNU * DNV * PN);
//	fin.close();
//
//
//
//	//Calculate B
//	thrust::host_vector<float> X(XN * YN * ZN, 0);
//	thrust::host_vector<float> B(XN * YN * ZN, 0);
//	thrust::host_vector<float> R(XN * YN * ZN, 0);
//	thrust::host_vector<float> nextR(XN * YN * ZN, 0);
//
//	thrust::host_vector<float> Delta(XN * YN * ZN, 0);
//	thrust::host_vector<float> TDelta(DNU * DNV *PN,0);
//	thrust::host_vector<float> MDelta(XN *YN*ZN, 0);
//
//	DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//				imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//				XN, YN, ZN, &B[0], hprj, dx, dz, mask, 0, 0, 0);
//
//
//	R = B;
//	Delta = R;
//	double alpha,beta;
//	double alpha1,alpha2;
//	double beta1, beta2;
//	clock_t start = clock();
//	for(int i = 0; i != iterNum; ++i)
//	{
//		//1 Calculate the scalar alpha
//		// (1) calculate MDelta
//		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, &Delta[0], &TDelta[0], dx, dz, mask, 0, 0);
//		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//					imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//					XN, YN, ZN, &MDelta[0], &TDelta[0], dx, dz, mask, 0, 0, 0);
//		// (2)
//		alpha1 = thrust::inner_product(R.begin(),R.end(),R.begin(),0.0);
//		alpha2 = thrust::inner_product(MDelta.begin(),MDelta.end(),Delta.begin(),0.0);
//		std::cout<<"alpha1 = "<<alpha1<<std::endl;
//		std::cout<<"alpha2 = "<<alpha2<<std::endl;
//		alpha = alpha1 / alpha2;
//		thrust::transform(X.begin(),X.end(),Delta.begin(),X.begin(),CGop<float>(alpha));
//		thrust::transform(R.begin(),R.end(),MDelta.begin(),nextR.begin(),CGop<float>(-alpha));
//		beta1 = thrust::inner_product(nextR.begin(),nextR.end(),nextR.begin(),0.0);
//		beta2 = thrust::inner_product(R.begin(),R.end(),R.begin(),0.0);
//		beta = beta1 / beta2;
//		thrust::transform(nextR.begin(),nextR.end(),Delta.begin(),Delta.begin(),CGop<float>(beta));
//		R = nextR;
//		std::cout<<i<<std::endl;
//
//		if(outputMedRes)
//		{
//			std::stringstream ss;
//			ss<<i;
//			std::string name = VolumeName + "_" + ss.str() + ".raw";
//			std::ofstream fu(name.c_str(), std::ios::binary);
//			fu.write((char*) &X[0], sizeof(float) * XN * YN * ZN);
//			fu.close();
//		}
//
//	}
//	clock_t end = clock();
//	std::cout<<"Total Time is "<<(static_cast<double>(end) - static_cast<double>(start)) / static_cast<double>(CLOCKS_PER_SEC)<<" seconds\n";
//	std::ofstream fou3(VolumeName.c_str(), std::ios::binary);
//	fou3.write((char*) &X[0], sizeof(float) * XN * YN * ZN);
//	fou3.close();
//}
//
//
//struct helicalCT_Geo
//{
//public:
//	thrust::host_vector<float> initImg;
//	thrust::host_vector<float> hangs;
//	thrust::host_vector<float> hzPos;
//	thrust::host_vector<byte> mask;
//	float sid;
//	float sdd;
//	int DNU;
//	int DNV;
//	int PN;
//	float imgXCenter;
//	float imgYCenter;
//	float imgZCenter;
//	int XN;
//	int YN;
//	int ZN;
//	float dx;
//	float dz;
//	float col_size;
//	float row_size;
//	float col_offset;
//	float row_offset;
//	float start_view;
//	float view_per_rot;
//	float PITCH;
//};
//
//
//template<typename T>
//struct minus :public thrust::binary_function<T,T,T>
//{
//	__host__ __device__ T operator()(const T& a, const T& b)
//	{
//		return a-b;
//	}
//};
//
//template<typename T>
//thrust::host_vector<T> operator-(const thrust::host_vector<T>& A, const thrust::host_vector<T>& B)
//{
//	thrust::host_vector<T> C(A.size(),0);
//	thrust::transform(A.begin(),A.end(),B.begin(),C.begin(),minus<T>());
//	return C;
//}
//
//void reconBigPatientFunctionCG(
//		const std::string& VolumeName,
//		thrust::host_vector<float>& X, // initial guess and the result
//		thrust::host_vector<float>& B, // backprojected of the subprojection
//		thrust::host_vector<float>& hangs,
//		thrust::host_vector<float>& hzPos,
//		thrust::host_vector<byte>& mask,
//		thrust::host_vector<float>& xds,
//		thrust::host_vector<float>& yds,
//		thrust::host_vector<float>& zds,
//		const int iterNum = 30,
//		const bool outputMedRes = false,
//		const float x0 = 0, const float y0 = 600, const float z0 = 0,
//		const int DNU = 888, const int DNV = 64, const int PN = 8139,
//		const float imgXCenter = 0.0f, const float imgYCenter = 0.0f, const float imgZCenter = 0.0f,
//		const int XN = 512, const int YN = 512, const int ZN = 512,
//		const float dx = 1.171875, const float dz = 1.171875)
//{
//	//Calculate B
//	thrust::host_vector<float> R(XN * YN * ZN, 0);
//	thrust::host_vector<float> nextR(XN * YN * ZN, 0);
//
//	thrust::host_vector<float> Delta(XN * YN * ZN, 0);
//	thrust::host_vector<float> TDelta(DNU * DNV *PN,0);
//	thrust::host_vector<float> MDelta(XN *YN*ZN, 0);
//
//
//	DD3Proj_gpu(x0, y0, z0, DNU, DNV, &xds[0], &yds[0], &zds[0],
//		imgXCenter, imgYCenter, imgZCenter, &hangs[0], &hzPos[0], PN,
//		XN, YN, ZN, &X[0], &TDelta[0], dx, dz, &mask[0], 0, 0);
//	DD3Back_gpu(x0, y0, z0, DNU, DNV,&xds[0], &yds[0], &zds[0],
//				imgXCenter, imgYCenter, imgZCenter, &hangs[0], &hzPos[0], PN,
//				XN, YN, ZN, &MDelta[0], &TDelta[0], dx, dz, &mask[0], 0, 0, 0);
//
//	R = B - MDelta;
//	Delta = R;
//	double alpha,beta;
//	double alpha1,alpha2;
//	double beta1, beta2;
//	clock_t start = clock();
//	for(int i = 0; i != iterNum; ++i)
//	{
//		//1 Calculate the scalar alpha
//		// (1) calculate MDelta
//		DD3Proj_gpu(x0, y0, z0, DNU, DNV, &xds[0], &yds[0], &zds[0],
//			imgXCenter, imgYCenter, imgZCenter, &hangs[0], &hzPos[0], PN,
//			XN, YN, ZN, &Delta[0], &TDelta[0], dx, dz, &mask[0], 0, 0);
//		DD3Back_gpu(x0, y0, z0, DNU, DNV, &xds[0], &yds[0], &zds[0],
//					imgXCenter, imgYCenter, imgZCenter, &hangs[0], &hzPos[0], PN,
//					XN, YN, ZN, &MDelta[0], &TDelta[0], dx, dz, &mask[0], 0, 0, 0);
//		// (2)
//		alpha1 = thrust::inner_product(R.begin(),R.end(),R.begin(),0.0);
//		alpha2 = thrust::inner_product(MDelta.begin(),MDelta.end(),Delta.begin(),0.0);
//		std::cout<<"alpha1 = "<<alpha1<<std::endl;
//		std::cout<<"alpha2 = "<<alpha2<<std::endl;
//		alpha = alpha1 / alpha2;
//		thrust::transform(X.begin(),X.end(),Delta.begin(),X.begin(),CGop<float>(alpha));
//		thrust::transform(R.begin(),R.end(),MDelta.begin(),nextR.begin(),CGop<float>(-alpha));
//		beta1 = thrust::inner_product(nextR.begin(),nextR.end(),nextR.begin(),0.0);
//		beta2 = thrust::inner_product(R.begin(),R.end(),R.begin(),0.0);
//		beta = beta1 / beta2;
//		thrust::transform(nextR.begin(),nextR.end(),Delta.begin(),Delta.begin(),CGop<float>(beta));
//		R = nextR;
//		std::cout<<i<<std::endl;
//
//		if(outputMedRes)
//		{
//			std::stringstream ss;
//			ss<<i;
//			std::string name = VolumeName + "_" + ss.str() + ".raw";
//			std::ofstream fu(name.c_str(), std::ios::binary);
//			fu.write((char*) &X[0], sizeof(float) * XN * YN * ZN);
//			fu.close();
//		}
//
//	}
//	clock_t end = clock();
//	std::cout<<"Total Time is "<<(static_cast<double>(end) - static_cast<double>(start)) / static_cast<double>(CLOCKS_PER_SEC)<<" seconds\n";
//	std::ofstream fou3(VolumeName.c_str(), std::ios::binary);
//	fou3.write((char*) &X[0], sizeof(float) * XN * YN * ZN);
//	fou3.close();
//}
//
//
//enum CG_Method{
//	CG_FR,
//	CG_PRP,
//	CG_CW,
//	CG_DI,
//	CG_DY};
////还没开始写
//void reconBigPatientFunctionCG(const std::string& ProjectionName,
//		const std::string& VolumeName,
//		CG_Method cg_method = CG_FR,
//		const int iterNum = 30,
//		const bool outputMedRes = false,
//		const float sid=538.5200193125, const float sdd=946.745971679699,
//		const int DNU = 888, const int DNV = 64, const int PN = 8139,
//		const float imgXCenter = 0.0f, const float imgYCenter = 0.0f, const float imgZCenter = 0.0f,
//		const int XN = 512, const int YN = 512, const int ZN = 512,
//		const float dx = 1.171875, const float dz = 1.171875,
//		const float col_size = 1.0239f, const float row_size = 1.0963f,
//		const float col_offset = -1.25f, const float row_offset = 0.0f,
//		const float start_view = 0, const int view_per_rot = 984,
//		const float PITCH = 63.0)
//{
//	//Define the parameters
//	float x0 = 0.0f;
//	float y0 = sid;
//	float z0 = 0.0f;
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//
//	//Generate the positions of the detectors
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = start_view + ii * TWOPI / static_cast<float>(view_per_rot);
//		hzPos[ii] = (static_cast<float>(ii)-static_cast<float>(PN) * 0.5)
//				/ static_cast<float>(view_per_rot) * 0.625 * PITCH;
//	}
//
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != YN; ++i)
//	{
//		for (int j = 0; j != XN; ++j)
//		{
//			if (sqrtf(powf((i - YN / 2.0 + 0.5) / (YN / 2.0), 2.0) + powf((j - XN / 2.0 + 0.5) / (XN / 2.0), 2.0)) < 0.98)
//			{
//				mask[i * XN + j] = 1;
//			}
//			else
//			{
//				mask[i * XN + j] = 0;
//			}
//		}
//	}
//	float* hprj = new float[DNU * DNV * PN];
//	//float* hvol = new float[XN * YN * ZN];
//
//
//	//Load the projection
//	std::ifstream fin(ProjectionName.c_str(),std::ios::binary);
//	if(!fin.is_open())
//	{
//		std::cerr<<"Cannot open the file\n";
//		exit(-1);
//	}
//	fin.read((char*)hprj,sizeof(float)* DNU * DNV * PN);
//	fin.close();
//
//
//
//	//Calculate B
//	thrust::host_vector<float> X(XN * YN * ZN, 0);
//	thrust::host_vector<float> B(XN * YN * ZN, 0);
//	thrust::host_vector<float> R(XN * YN * ZN, 0);
//	thrust::host_vector<float> nextR(XN * YN * ZN, 0);
//
//	thrust::host_vector<float> Delta(XN * YN * ZN, 0);
//	thrust::host_vector<float> TDelta(DNU * DNV *PN,0);
//	thrust::host_vector<float> MDelta(XN *YN*ZN, 0);
//
//	DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//				imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//				XN, YN, ZN, &B[0], hprj, dx, dz, mask, 0, 0, 0);
//
//
//	R = B;
//	Delta = R;
//	double alpha,beta;
//	double alpha1,alpha2;
//	double beta1, beta2;
//	clock_t start = clock();
//	for(int i = 0; i != iterNum; ++i)
//	{
//		//1 Calculate the scalar alpha
//		// (1) calculate MDelta
//		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, &Delta[0], &TDelta[0], dx, dz, mask, 0, 0);
//		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//					imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//					XN, YN, ZN, &MDelta[0], &TDelta[0], dx, dz, mask, 0, 0,0);
//		// (2)
//		alpha1 = thrust::inner_product(R.begin(),R.end(),R.begin(),0.0);
//		alpha2 = thrust::inner_product(MDelta.begin(),MDelta.end(),Delta.begin(),0.0);
//		std::cout<<"alpha1 = "<<alpha1<<std::endl;
//		std::cout<<"alpha2 = "<<alpha2<<std::endl;
//		alpha = alpha1 / alpha2;
//		thrust::transform(X.begin(),X.end(),Delta.begin(),X.begin(),CGop<float>(alpha));
//		thrust::transform(R.begin(),R.end(),MDelta.begin(),nextR.begin(),CGop<float>(-alpha));
//		beta1 = thrust::inner_product(nextR.begin(),nextR.end(),nextR.begin(),0.0);
//		beta2 = thrust::inner_product(R.begin(),R.end(),R.begin(),0.0);
//		beta = beta1 / beta2;
//		thrust::transform(nextR.begin(),nextR.end(),Delta.begin(),Delta.begin(),CGop<float>(beta));
//		R = nextR;
//		std::cout<<i<<std::endl;
//
//		if(outputMedRes)
//		{
//			std::stringstream ss;
//			ss<<i;
//			std::string name = VolumeName + "_" + ss.str() + ".raw";
//			std::ofstream fu(name.c_str(), std::ios::binary);
//			fu.write((char*) &X[0], sizeof(float) * XN * YN * ZN);
//			fu.close();
//		}
//
//	}
//	clock_t end = clock();
//	std::cout<<"Total Time is "<<(static_cast<double>(end) - static_cast<double>(start)) / static_cast<double>(CLOCKS_PER_SEC)<<" seconds\n";
//	std::ofstream fou3(VolumeName.c_str(), std::ios::binary);
//	fou3.write((char*) &X[0], sizeof(float) * XN * YN * ZN);
//	fou3.close();
//}
//
//
//void reconSheppLogan20151108()
//{
//	float sid = 538.52f;
//	float sdd = 946.75f;
//
//	float x0(0.0f);
//	float y0(sid);
//	float z0(0.0f);
//
//	int DNU(888);
//	int DNV(64);
//	int PN(2201);
//	float imgXCenter(0);
//	float imgYCenter(0);
//	float imgZCenter(0);
//
//	int XN(512);
//	int YN(512);
//	int ZN(64); //mouse volume is larger
//
//
//	float dx(500.0f / 512.0f);
//	float dz(0.625);
//
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//	//Generate the positions of the detectors
//	float col_size = 1.0239;
//	float row_size = 1.0963;
//
//	float col_offset = 0;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//	imgXCenter = 0;
//	imgYCenter = 0;
//	imgZCenter = 0;
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(PN);
//		hzPos[ii] = 0;// (ii - PN / 2) * 0.0015;
//	}
//
//	float* hvol = new float[XN * YN * ZN];
//	float* hprj = new float[DNU * DNV * PN]; // Each time we only perform on projection
//
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != YN; ++i)
//	{
//		for (int j = 0; j != XN; ++j)
//		{
//			if (sqrtf(powf((i - YN / 2.0 + 0.5) / (YN / 2.0), 2.0) + powf((j - XN / 2.0 + 0.5) / (XN / 2.0), 2.0)) < 0.98)
//			{
//				mask[i * XN + j] = 1;
//			}
//			else
//			{
//				mask[i * XN + j] = 0;
//			}
//		}
//	}
//
//	//Generate projection
//	std::ifstream fin("sheppLoganReal.raw",std::ios::binary);
//	fin.read((char*)hvol,sizeof(float)* XN * YN * ZN);
//	fin.close();
//	DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, hvol, hprj, dx, dz, mask, 0, 0);
//
//	float* oneVol = new float[XN * YN * ZN];
//	float* onePrj = new float[DNU * DNV * PN];
//	float* rowSum = new float[DNU * DNV * PN];
//	float* colSum = new float[XN * YN * ZN];
//	float* lasImg = new float[XN * YN * ZN];
//
//	for (int ii = 0; ii != XN * YN * ZN; ++ii)
//	{
//		hvol[ii] = 0;
//		oneVol[ii] = 1.0;
//		colSum[ii] = 0;
//	}
//	for (int ii = 0; ii != DNU * DNV * PN; ++ii)
//	{
//		onePrj[ii] = 1.0;
//		rowSum[ii] = 0.0;
//	}
//
//	DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, oneVol, rowSum, dx, dz, mask, 0, 0);
//
//	DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, colSum, onePrj, dx, dz, mask, 0, 0, 0);
//	float t1 = 1.0;
//	float t2 = 1.0;
//	int totN = XN * YN * ZN;
//	for (int ii = 0; ii != 800; ++ii)
//	{
//		memcpy(lasImg, hvol, sizeof(float) * totN);
//
//		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, hvol, onePrj, dx, dz, mask, 0, 0);
//
//		prjWeight(onePrj, hprj, rowSum, DNU*DNV*PN);
//
//		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, oneVol, onePrj, dx, dz, mask, 0, 0, 0);
//
//		bakWeight(oneVol, hvol, colSum, totN);
//
//		t2 = (1.0 + sqrtf(1.0 + 4.0 * t1 * t1)) / 2.0;
//		FISTA(lasImg, hvol, t1, t2, totN);
//		t1 = t2;
//
//
//		if ((ii + 1) % 20 == 0)
//		{
//			std::stringstream ss;
//			ss << ii;
//			std::string name = "SheppRecon" + ss.str() + ".raw";
//			std::ofstream fou3(name.c_str(), std::ios::binary);
//			fou3.write((char*) hvol, sizeof(float) * XN * YN * ZN);
//			fou3.close();
//		}
//		std::cout << "iteration # = "<<ii<<std::endl;
//	}
//
//}
//
//
//void reconRealPatient20151111()
//{
//	float sid = 538.52;
//	float sdd = 946.75;
//
//	float x0(0.0f);
//	float y0(sid);
//	float z0(0.0f);
//
//	int DNU(888);
//	int DNV(64);
//	int PN(2200);
//	float imgXCenter(0);
//	float imgYCenter(0);
//	float imgZCenter(0);
//
//	int XN(512);
//	int YN(512);
//	int ZN(64); //mouse volume is larger
//
//
//	float dx(500.0f / 512.0f);
//	float dz(0.625);
//
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//	//Generate the positions of the detectors
//	float col_size = 1.0239;
//	float row_size = 1.0963;
//
//
//	float col_offset = -1.25;
//	float row_offset = 0;
//
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//	imgXCenter = 0;
//	imgYCenter = 0;
//	imgZCenter = 0;
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = ii * TWOPI / static_cast<float>(PN);
//		hzPos[ii] = 0;// (ii - PN / 2) * 0.0015;
//	}
//
//	float* hvol = new float[XN * YN * ZN];
//	float* hprj = new float[DNU * DNV * PN]; // Each time we only perform on projection
//
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != YN; ++i)
//	{
//		for (int j = 0; j != XN; ++j)
//		{
//			if (sqrtf(powf((i - YN / 2.0 + 0.5) / (YN / 2.0), 2.0) + powf((j - XN / 2.0 + 0.5) / (XN / 2.0), 2.0)) < 0.98)
//			{
//				mask[i * XN + j] = 1;
//			}
//			else
//			{
//				mask[i * XN + j] = 0;
//			}
//		}
//	}
//
//	//Generate projection
//	std::ifstream fin("RealPatientProj1.raw",std::ios::binary);
//	fin.read((char*)hprj,sizeof(float)* DNU * DNV * PN);
//	fin.close();
//
//
//
//	float* oneVol = new float[XN * YN * ZN];
//	float* onePrj = new float[DNU * DNV * PN];
//	float* rowSum = new float[DNU * DNV * PN];
//	float* colSum = new float[XN * YN * ZN];
//	float* lasImg = new float[XN * YN * ZN];
//
//	for (int ii = 0; ii != XN * YN * ZN; ++ii)
//	{
//		hvol[ii] = 0;
//		oneVol[ii] = 1.0;
//		colSum[ii] = 0;
//	}
//	for (int ii = 0; ii != DNU * DNV * PN; ++ii)
//	{
//		onePrj[ii] = 1.0;
//		rowSum[ii] = 0.0;
//	}
//
//	DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, oneVol, rowSum, dx, dz, mask, 0, 0);
//
//	DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//		imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//		XN, YN, ZN, colSum, onePrj, dx, dz, mask, 0, 0, 0);
//	float t1 = 1.0;
//	float t2 = 1.0;
//	int totN = XN * YN * ZN;
//	for (int ii = 0; ii != 800; ++ii)
//	{
//		memcpy(lasImg, hvol, sizeof(float) * totN);
//
//		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, hvol, onePrj, dx, dz, mask, 0, 0);
//
//		prjWeight(onePrj, hprj, rowSum, DNU*DNV*PN);
//
//		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
//			imgXCenter, imgYCenter, imgZCenter, hangs, hzPos, PN,
//			XN, YN, ZN, oneVol, onePrj, dx, dz, mask, 0, 0, 0);
//
//		bakWeight(oneVol, hvol, colSum, totN);
//
//		t2 = (1.0 + sqrtf(1.0 + 4.0 * t1 * t1)) / 2.0;
//		FISTA(lasImg, hvol, t1, t2, totN);
//		t1 = t2;
//
//
//		if ((ii + 1) % 20 == 0)
//		{
//			std::stringstream ss;
//			ss << ii;
//			std::string name = "newRealPatientRecon" + ss.str() + ".raw";
//			std::ofstream fou3(name.c_str(), std::ios::binary);
//			fou3.write((char*) hvol, sizeof(float) * XN * YN * ZN);
//			fou3.close();
//		}
//		std::cout << "iteration # = "<<ii<<std::endl;
//	}
//
//}
//
//
/////////////////////////////////////////////////////////////////
////void testSIEMENS_trainingData()
////{
////	const int DNU = 736;
////	const int DNV = 64;
////	const int PN = 6074;
////	const int XN = 512;
////	const int YN = 512;
////	const int ZN = 560;
////	const float sid = 592.1447;
////	const float sdd = 1085.6;
////	const float dx = 500.0 / 512.0;
////	const float dz = 1.0;
////	const float col_size = 1.2858;
////	const float row_size = 1.0947;
////	const float col_offset = 2.125; // May be changed;
////	const float row_offset = 1.0; // May be changed
////
////	const float imgXCenter = 0;
////	const float imgYCenter = 0;
////	const float imgZCenter = 0;
////
////	const float x0 = 0;
////	const float y0 = sid;
////	const float z0 = 0;
////
////	const int iterNum = 6;// Iteration #
////	const int osNum = 1; // eight subsets are applied for CG
////	const int cgIter = 2; // each subset, 2 iteration are involved
////
////	//Open 8 files for CG
////	std::string volName[osNum];
////
////
////	thrust::host_vector<float> subProj[osNum];
////	thrust::host_vector<float> B[osNum];
////	thrust::host_vector<float> hangs[osNum];
////	thrust::host_vector<float> hzPos[osNum];
////
////	thrust::host_vector<byte> mask(XN * YN, 1);
////
////	thrust::host_vector<float> xds(DNU, 0);
////	thrust::host_vector<float> yds(DNU, 0);
////	thrust::host_vector<float> zds(DNV, 0);
////
////	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
////	float curBeta = 0;
////	for (int ii = 0; ii != DNU; ++ii)
////	{
////		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
////		xds[ii] = sinf(curBeta) * sdd;
////		yds[ii] = sid - cosf(curBeta) * sdd;
////	}
////
////	for (int ii = 0; ii != DNV; ++ii)
////	{
////		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
////	}
////
////
////	std::string subBName[osNum];
////	for(int i = 0; i != osNum; ++i)
////	{
////		std::stringstream ss;
////		ss<<i;
////
////		subBName[i] = "subVolume" + ss.str() + ".raw";
////	}
////
////
////	for(int i = 0; i != osNum; ++i)
////	{
////
////		volName[i] = "vol_os_";
////		// Read sub projections
////		subProj[i].resize(DNU * DNV * PN, 0);
////		std::stringstream ss;
////		ss<<i;
////		std::string name = "testSIEMENS_subProj_" + ss.str() + ".prj";
////		std::ifstream fid(name.c_str(),std::ios::binary);
////		if(!fid.is_open())
////		{
////			std::cerr<<"Cannot open the projection\n";
////			exit(-1);
////		}
////		fid.read((char*)(&subProj[i][0]), sizeof(float) * DNU * DNV * PN);
////		fid.close();
////
////		hangs[i].resize(PN, 0);
////		hzPos[i].resize(PN, 0);
////
////		std::string haName = "subAngle_" + ss.str() + ".prj";
////		std::ifstream fanIn(haName.c_str(),std::ios::binary);
////		if(!fanIn.is_open())
////		{
////			std::cerr<<"Cannot open the angle files \n";
////			exit(-1);
////		}
////		fanIn.read((char*)(&hangs[i][0]), sizeof(float) * PN);
////		fanIn.close();
////
////		std::string hzName = "sub_ZPos" + ss.str() + ".prj";
////		std::ifstream fzpIn(hzName.c_str(),std::ios::binary);
////		if(!fzpIn.is_open())
////		{
////			std::cerr<<"Cannot open the zpos files\n";
////			exit(-1);
////		}
////		fzpIn.read((char*)(&hzPos[i][0]), sizeof(float) * PN);
////		fzpIn.close();
////
////
////		B[i].resize(XN * YN * ZN, 0);
////		std::ifstream fin(subBName[i].c_str(),std::ios::binary);
////		if(!fin.is_open())
////		{
////			//Backprojection
////			DD3Back_gpu(x0, y0, z0, DNU, DNV, &xds[0], &yds[0], &zds[0],
////						imgXCenter, imgYCenter, imgZCenter,
////						&(hangs[i][0]), &(hzPos[i][0]),
////						PN,	XN, YN, ZN, &(B[i][0]),
////						&(subProj[i][0]),
////						dx, dz, &mask[0], 0, 0, 0);
////			std::ofstream fou(subBName[i].c_str(),std::ios::binary);
////			fou.write((char*)(&B[i][0]), sizeof(float) * XN * YN * ZN);
////			fou.close();
////		}
////		else
////		{
////			fin.read((char*)(&B[i][0]), sizeof(float) * XN * YN * ZN);
////			fin.close();
////		}
////
////		std::cout<<i<<std::endl;
////	}
////
////
////
////
////	thrust::host_vector<float> X(XN * YN * ZN, 0);
////
////
////	for(int iterIdx = 0; iterIdx != iterNum; ++iterIdx)
////	{
////		for(int osIdx = 0; osIdx != osNum; ++osIdx)
////		{
////
////			reconBigPatientFunctionCG(volName[osIdx],X,
////					B[osIdx], hangs[osIdx], hzPos[osIdx],
////					mask, xds, yds, zds, cgIter, true,
////					x0, y0, z0, DNU,DNV,PN,imgXCenter, imgYCenter,imgZCenter,
////					XN,YN,ZN,dx,dz);
////			std::cout<<"osIDX = " << osIdx << std::endl;
////		}
////		std::cout<<"iterIDX = " << iterIdx << std::endl;
////	}
////	std::ofstream foutName("SIEMENS_L067.raw",std::ios::binary);
////	foutName.write((char*)(&X[0]), sizeof(float) * XN * YN * ZN);
////	foutName.close();
////}
////
//
//
////This is Right. Do NOT to change the materials here.
//void testSIEMENS_trainingData_smallData()
//{
//	const int DNU = 736;
//	const int DNV = 64;
//	const int PN = 8099;
//	const int XN = 512;
//	const int YN = 512;
//	const int ZN = 256;
//	const float sid = 592.1447;
//	const float sdd = 1085.6;
//	const float dx = 500.0 / 512.0;
//	const float dz = 1.0;
//	const float col_size = 1.2858;
//	const float row_size = 1.0947;
//	const float col_offset = 2.125; // May be changed;
//	const float row_offset = 1.0; // May be changed
//
//	const float imgXCenter = 0;
//	const float imgYCenter = 0;
//	const float imgZCenter = 0;
//
////	reconBigPatientFunctionCG("testSIEMENS.prj","testPartSIEMENS_20160318.raw",
////			40, true,sid, sdd, DNU, DNV, PN,
////			imgXCenter, imgYCenter, imgZCenter,
////			XN, YN, ZN, dx, dz, col_size, row_size,
////			col_offset, row_offset,
////			-4.0 * 3.1415926, 2304, 64 * 0.6);
//
//	reconBigPatientFunctionCG("subProj1in6.prj","testPartSIEMENS_20160328.raw",
//			40, true,sid, sdd, DNU, DNV, PN,
//			imgXCenter, imgYCenter, imgZCenter,
//			XN, YN, ZN, dx, dz, col_size, row_size,
//			col_offset, row_offset,
//			-4.0 * 3.1415926, 2304, 64 * 0.6);
//
//}
//
//#include <algorithm>
//
//void test3GPU_Proj(const std::string& ProjectionName,
//		const std::string& VolumeName,
//		const int iterNum = 30,
//		const bool outputMedRes = false,
//		const float sid=538.5200193125, const float sdd=946.745971679699,
//		const int DNU = 888, const int DNV = 64, const int PN = 8139,
//		const float imgXCenter = 0.0f, const float imgYCenter = 0.0f, const float imgZCenter = 0.0f,
//		const int XN = 512, const int YN = 512, const int ZN = 512,
//		const float dx = 1.171875, const float dz = 1.171875,
//		const float col_size = 1.0239f, const float row_size = 1.0963f,
//		const float col_offset = -1.25f, const float row_offset = 0.0f,
//		const float start_view = 0, const int view_per_rot = 984,
//		const float PITCH = 63.0)
//{
//	//Define the parameters
//	float x0 = 0.0f;
//	float y0 = sid;
//	float z0 = 0.0f;
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//
//	//Generate the positions of the detectors
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii)
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//
//	float* hangs = new float[PN];
//	float* hzPos = new float[PN];
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = start_view + ii * TWOPI / static_cast<float>(view_per_rot);
//		hzPos[ii] = (static_cast<float>(ii)-static_cast<float>(PN) * 0.5)
//				/ static_cast<float>(view_per_rot) * 64;
//	}
//
//
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != YN; ++i)
//	{
//		for (int j = 0; j != XN; ++j)
//		{
//			if (sqrtf(powf((i - YN / 2.0 + 0.5) / (YN / 2.0), 2.0) + powf((j - XN / 2.0 + 0.5) / (XN / 2.0), 2.0)) < 0.98)
//			{
//				mask[i * XN + j] = 1;
//			}
//			else
//			{
//				mask[i * XN + j] = 1;
//			}
//		}
//	}
//	float* hprj = new float[DNU * DNV * PN];
//	float* hvol = new float[XN * YN * ZN];
//	std::ifstream fin(VolumeName.c_str(),std::ios::binary);
//	for(int i = 0; i != XN * YN * ZN; ++i)
//	{
//		hvol[i] = 1.0;
//	}
//	fin.read((char*)hvol,sizeof(float) * XN * YN * ZN);
//	fin.close();
//
//	int PNN = 4700;
//	int RES = (PN - PNN) / 2;
//	int NXPNN = PNN + RES;
//	int startPN[3] = {0, PNN, NXPNN};
//	int startVOL[3] = {0, 258, 323};
//
//
//	DD3ProjHelical_3GPU(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
//		hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz, mask, 3, startPN);
//
//	const int TPN = DNU * DNV * PN;
//	const int TVN = XN * YN * ZN;
//
//	float* rowSum = new float[TPN];
//	float* colSum = new float[TVN];
//
//	float* oneVol = new float[TVN];
//	float* onePrj = new float[TPN];
//
//	omp_set_num_threads(32);
//#pragma omp parallel for
//	for(int i = 0; i < TPN; ++i)
//	{
//		rowSum[i] = 0;
//		onePrj[i] = 1.0f;
//	}
//
//#pragma omp parallel for
//	for(int i = 0; i < TVN; ++i)
//	{
//		colSum[i] = 0;
//		oneVol[i] = 1.0f;
//		hvol[i] = 0.0;
//	}
//
//	DD3ProjHelical_3GPU(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
//		hangs, hzPos, PN, XN, YN, ZN, oneVol, rowSum, dx, dz, mask, 3, startPN);
//
//	DD3BackHelical_3GPU(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
//			hangs, hzPos, PN, XN, YN, ZN, colSum, onePrj, dx, dz,mask, 3, startVOL);
//	//delete[] oneVol;
//	//delete[] onePrj;
//	float* lastV = new float[TVN];
//	float t0 = 1;
//	float t1 = 1;
//
//	for(size_t iterIdx = 0; iterIdx != 400; ++iterIdx)
//	{
//		std::copy(hvol,hvol+TVN,lastV);
//
//		//Projection
//		DD3ProjHelical_3GPU(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
//			hangs, hzPos, PN, XN, YN, ZN, hvol, onePrj, dx, dz, mask, 3, startPN);
//		//Weighting
//#pragma omp parallel for
//		for(int pI = 0; pI < TPN; ++pI)
//		{
//			onePrj[pI] = hprj[pI] - onePrj[pI];
//
//			if(abs(rowSum[pI]) < 1.0E-9)
//			{
//				onePrj[pI] = 0;
//			}
//			else
//			{
//				onePrj[pI] = onePrj[pI] / rowSum[pI];
//			}
//		}
//		//Backprojection
//		DD3BackHelical_3GPU(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
//			hangs, hzPos, PN, XN, YN, ZN, oneVol, onePrj, dx, dz,mask, 3, startVOL);
//
//		//Add and weighting
//#pragma omp parallel for
//		for(int vI = 0; vI < TVN; ++vI)
//		{
//			if(abs(colSum[vI]) < 1.0E-9)
//			{
//				hvol[vI] += 0;
//			}
//			else
//			{
//				hvol[vI] += (oneVol[vI] / colSum[vI]);
//			}
//
//		}
//
//		t1 = (1.0 + sqrt(1.0 + 4.0 * t0 * t0)) / 2.0;
//#pragma omp parallel for
//		for(int vI = 0; vI < TVN; ++vI)
//		{
//			hvol[vI] = hvol[vI] + (t0 - 1.0) / t1 * (hvol[vI] - lastV[vI]);
//		}
//
//		t0 = t1;
//		std::cout<<iterIdx<<"\n";
//
//		std::stringstream ss;
//		ss<<iterIdx;
//		std::string Name = "temp" + ss.str() + ".raw";
//		if(iterIdx%8==0)
//		{
//			std::ofstream fff(Name.c_str(), std::ios::binary);
//			fff.write((char*) &hvol[0], sizeof(float) * XN * YN * ZN);
//			fff.close();
//		}
//
//
//	}
//
//
//
//
//	std::ofstream fou3("TestVOL.raw", std::ios::binary);
//	fou3.write((char*) &hvol[0], sizeof(float) * XN * YN * ZN);
//	fou3.close();
//	delete[] lastV;
//	delete[] oneVol;
//	delete[] onePrj;
//	delete[] rowSum;
//	delete[] colSum;
//
//	delete[] hprj;
//	delete[] hvol;
//
//}
//
//
//
//void testFromZ0PositionScan(const std::string& ProjectionName,
//		const std::string& VolumeName,
//		const int iterNum = 30,
//		const bool outputMedRes = false,
//		const float sid=538.5200193125, const float sdd=946.745971679699,
//		const int DNU = 888, const int DNV = 64, const int PN = 8139,
//		const float imgXCenter = 0.0f, const float imgYCenter = 0.0f, const float imgZCenter = 0.0f,
//		const int XN = 512, const int YN = 512, const int ZN = 512,
//		const float dx = 1.171875, const float dz = 1.171875,
//		const float col_size = 1.0239f, const float row_size = 1.0963f,
//		const float col_offset = -1.25f, const float row_offset = 0.0f,
//		const float start_view = 0, const int view_per_rot = 984,
//		const float PITCH = 63.0)
//{
//	//Define the parameters
//	float x0 = 0.0f;
//	float y0 = sid;
//	float z0 = 0.0f;
//
//	float* xds = new float[DNU];
//	float* yds = new float[DNU];
//	float* zds = new float[DNV];
//
//	//Generate the positions of the detectors
//
//	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
//	float curBeta = 0;
//	for (int ii = 0; ii != DNU; ++ii) //This will be calculated in MATLAB
//	{
//		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
//		xds[ii] = sinf(curBeta) * sdd;
//		yds[ii] = sid - cosf(curBeta) * sdd;
//	}
//
//	for (int ii = 0; ii != DNV; ++ii)
//	{
//		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
//	}
//
//
//	float* hangs = new float[PN]; // This will be calculated in MATLAB
//	float* hzPos = new float[PN]; // This will be calculated in MATLAB
//
//	for (int ii = 0; ii != PN; ++ii)
//	{
//		hangs[ii] = start_view + ii * TWOPI / static_cast<float>(view_per_rot);
//		hzPos[ii] = static_cast<float>(ii) / static_cast<float>(view_per_rot) * 64; //pitch == 1
//	}
//
//	// This will be calculated in MATLAB
//	byte* mask = new byte[XN * YN];
//	for (int i = 0; i != YN; ++i)
//	{
//		for (int j = 0; j != XN; ++j)
//		{
//			if (sqrtf(powf((i - YN / 2.0 + 0.5) / (YN / 2.0), 2.0) + powf((j - XN / 2.0 + 0.5) / (XN / 2.0), 2.0)) < 0.98)
//			{
//				mask[i * XN + j] = 1;
//			}
//			else
//			{
//				mask[i * XN + j] = 1;
//			}
//		}
//	}
//
//
//	float* hprj = new float[DNU * DNV * PN];
//	float* hvol = new float[XN * YN * ZN];
//	std::ifstream fin(VolumeName.c_str(),std::ios::binary);
//	for(int i = 0; i != XN * YN * ZN; ++i)
//	{
//		hvol[i] = 1.0;
//	}
//	fin.read((char*)hvol,sizeof(float) * XN * YN * ZN);
//	fin.close();
//
//	int PNN = 4700;
//	int RES = (PN - PNN) / 2;
//	int NXPNN = PNN + RES;
//	int startPN[3] = {0, PNN, NXPNN};
//	int startVOL[3] = {0, 258, 323};
//
//
//	DD3ProjHelical_3GPU(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
//		hangs, hzPos, PN, XN, YN, ZN, hvol, hprj, dx, dz, mask, 3, startPN);
//
//	const int TPN = DNU * DNV * PN;
//	const int TVN = XN * YN * ZN;
//
//	float* rowSum = new float[TPN];
//	float* colSum = new float[TVN];
//
//	float* oneVol = new float[TVN];
//	float* onePrj = new float[TPN];
//
//	omp_set_num_threads(32);
//#pragma omp parallel for
//	for(int i = 0; i < TPN; ++i)
//	{
//		rowSum[i] = 0;
//		onePrj[i] = 1.0f;
//	}
//
//#pragma omp parallel for
//	for(int i = 0; i < TVN; ++i)
//	{
//		colSum[i] = 0;
//		oneVol[i] = 1.0f;
//		hvol[i] = 0.0;
//	}
//
//	DD3ProjHelical_3GPU(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
//		hangs, hzPos, PN, XN, YN, ZN, oneVol, rowSum, dx, dz, mask, 3, startPN);
//
//	DD3BackHelical_3GPU(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
//			hangs, hzPos, PN, XN, YN, ZN, colSum, onePrj, dx, dz,mask, 3, startVOL);
//	//delete[] oneVol;
//	//delete[] onePrj;
//	float* lastV = new float[TVN];
//	float t0 = 1;
//	float t1 = 1;
//
//	for(size_t iterIdx = 0; iterIdx != 400; ++iterIdx)
//	{
//		std::copy(hvol,hvol+TVN,lastV);
//
//		//Projection
//		DD3ProjHelical_3GPU(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
//			hangs, hzPos, PN, XN, YN, ZN, hvol, onePrj, dx, dz, mask, 3, startPN);
//		//Weighting
//#pragma omp parallel for
//		for(int pI = 0; pI < TPN; ++pI)
//		{
//			onePrj[pI] = hprj[pI] - onePrj[pI];
//
//			if(abs(rowSum[pI]) < 1.0E-9)
//			{
//				onePrj[pI] = 0;
//			}
//			else
//			{
//				onePrj[pI] = onePrj[pI] / rowSum[pI];
//			}
//		}
//		//Backprojection
//		DD3BackHelical_3GPU(x0, y0, z0, DNU, DNV, xds, yds, zds, imgXCenter, imgYCenter, imgZCenter,
//			hangs, hzPos, PN, XN, YN, ZN, oneVol, onePrj, dx, dz,mask, 3, startVOL);
//
//		//Add and weighting
//#pragma omp parallel for
//		for(int vI = 0; vI < TVN; ++vI)
//		{
//			if(abs(colSum[vI]) < 1.0E-9)
//			{
//				hvol[vI] += 0;
//			}
//			else
//			{
//				hvol[vI] += (oneVol[vI] / colSum[vI]);
//			}
//
//		}
//
//		t1 = (1.0 + sqrt(1.0 + 4.0 * t0 * t0)) / 2.0;
//#pragma omp parallel for
//		for(int vI = 0; vI < TVN; ++vI)
//		{
//			hvol[vI] = hvol[vI] + (t0 - 1.0) / t1 * (hvol[vI] - lastV[vI]);
//		}
//
//		t0 = t1;
//		std::cout<<iterIdx<<"\n";
//
//		std::stringstream ss;
//		ss<<iterIdx;
//		std::string Name = "temp" + ss.str() + ".raw";
//		if(iterIdx%8==0)
//		{
//			std::ofstream fff(Name.c_str(), std::ios::binary);
//			fff.write((char*) &hvol[0], sizeof(float) * XN * YN * ZN);
//			fff.close();
//		}
//
//
//	}
//
//	std::ofstream fou3("TestVOL.raw", std::ios::binary);
//	fou3.write((char*) &hvol[0], sizeof(float) * XN * YN * ZN);
//	fou3.close();
//	delete[] lastV;
//	delete[] oneVol;
//	delete[] onePrj;
//	delete[] rowSum;
//	delete[] colSum;
//
//	delete[] hprj;
//	delete[] hvol;
//
//}
//
//// Test the projection splitting
//void testSIEMENS_ProjSplit()
//{
//	const int DNU = 736;
//	const int DNV = 64;
//	//const int PN = 2034 * 4 * 5; //This is the largest for projection for 3 GPUs
//	const int PN = 2034 * 4;
//	const int XN = 512;
//	const int YN = 512;
//	const int ZN = 512;
//	const float sid = 592.1447;
//	const float sdd = 1085.6;
//	const float dx = 500.0 / 512.0;
//	const float dz = 1.0;
//	const float col_size = 1.2858;
//	const float row_size = 1.0947;
//	const float col_offset = 2.125; // May be changed;
//	const float row_offset = 1.0; // May be changed
//
//	const float imgXCenter = 0;
//	const float imgYCenter = 0;
//	const float imgZCenter = 0 + 1.0 * 256;
//
//	testFromZ0PositionScan("test3GPUProj.prj","test3GPUProj_vol.raw",
//			40, true,sid, sdd, DNU, DNV, PN,
//			imgXCenter, imgYCenter, imgZCenter,
//			XN, YN, ZN, dx, dz, col_size, row_size,
//			col_offset, row_offset,
//			-4.0 * 3.1415926, 2304, 64 * 0.6);
//
////	test3GPU_Proj("test3GPUProj.prj","test3GPUProj_vol.raw",
////			40, true,sid, sdd, DNU, DNV, PN,
////			imgXCenter, imgYCenter, imgZCenter,
////			XN, YN, ZN, dx, dz, col_size, row_size,
////			col_offset, row_offset,
////			-4.0 * 3.1415926, 2304, 64 * 0.6);
//
//
//}
//
//
////int main()
////{
//	//testSIEMENS_trainingData();
//
//	//testSIEMENS_trainingData_smallData();
//
////	testSIEMENS_ProjSplit();
//
//
////	std::string haoName = "HaoProj.prj";
////	std::string haoRecon = "HaoVol.raw";
////	reconGongHaoDataFunction(haoName,haoRecon,1,80,false,
////			52.20373,67.46248,
////			725, 253, 720,
////			0.0f, 0.0f, 0.0f,
////			512, 512, 200,
////			0.03095, 0.03095,
////			0.02, 0.02,
////			1.3477, 0.0f,
////			0, 720,  0);
////
//
////
////	std::string ProjectionName1 = "fake_BigPatProj1.prj"; // Only used for test the speed performance
////	std::string VolumeName11 = "VolumeOSSART1_";
////	std::string VolumeName12 = "VolumeCG1_";
////	std::string ProjectionName2 = "BigPatProj2.prj";
////	std::string VolumeName21 = "VolumeOSSART2_";
////	std::string VolumeName22 = "VolumeCG2_";
////
////	std::string ProjectionName3 = "BigPatProj3.prj";
////	std::string VolumeName31 = "VolumeOSSART3_";
////	std::string VolumeName32 = "VolumeCG3_";
///**
// * CG with Zeng Larry's method
// */
////	reconBigPatientFunctionCG(ProjectionName1, VolumeName12,
////			2, true,538.52, 946.75, 888, 64, 8139,
////			0.0f, 0.0f, 0.0f, 512, 512, 512,
////			1.171875,1.171875,
////			1.0239f,1.0963f,
////			-1.25,0,
////			-8139.0/984.0*3.1415926, 984, 63.0);
//
////
////	reconBigPatientFunctionCG(ProjectionName2, VolumeName22,
////			40, true,538.52, 946.75, 888, 64, 12950,
////			0.0f, 0.0f, 0.0f, 512, 512, 512,
////			1.171875, 1.171875,
////			1.0239f,1.0963f,
////			-1.25,0,
////			-12950.0/984.0*3.1415926, 984, 63.0);
////
////	reconBigPatientFunctionCG(ProjectionName3, VolumeName32,
////			40, true,538.52, 946.75, 888, 64, 6420,
////			0.0f, 0.0f, 0.0f, 512, 512, 512,
////			1.171875, 1.171875,
////			1.0239f,1.0963f,
////			-1.25,0,
////			-6420.0/984.0*3.1415926, 984, 88.0);
//
//
///**
// * OS-SART with 10 OS, iter = 30
// */
////
////	clock_t start = clock();
////
////	reconBigPatientFunction(ProjectionName1, VolumeName11, 1, 3, false,
////			538.52, 946.75,
////			888, 64, 8139,
////			0.0f, 0.0f,0.0f,
////			512, 512, 512,
////			1.171875, 1.171875,
////			1.0239f, 1.0963f,
////			-1.25f,0.0f,
////			-8139.0/984.0*3.1415926, 984, 63.0);
////	clock_t end = clock();
////	std::cout<<(static_cast<double>(end)-static_cast<double>(start)) / CLOCKS_PER_SEC;
//
//////	reconBigPatientFunction(ProjectionName2, VolumeName21, 20, 30, false,
//////				538.52, 946.75,
//////				888, 64, 12950,
//////				0.0f, 0.0f,0.0f,
//////				512, 512, 512,
//////				1.171875, 1.171875,
//////				1.0239f, 1.0963f,
//////				-1.25f,0.0f,
//////				-12950.0/984.0*3.1415926, 984, 63.0);
//////
//////	reconBigPatientFunction(ProjectionName3, VolumeName31, 20, 30, false,
//////				538.52, 946.75,
//////				888, 64, 6420,
//////				0.0f, 0.0f,0.0f,
//////				512, 512, 512,
//////				1.171875, 1.171875,
//////				1.0239f, 1.0963f,
//////				-1.25f,0.0f,
//////				-6420.0/984.0*3.1415926, 984, 88.0);
////
////
//////	//TestProjectionTimeWithNumOfGPUs();
//////	//projFORBILD();
//	//TestProjectionTime3();
//////	reconFORBILD20151105();
//////	reconSheppLogan20151108();
//////	//reconRealPatient20151111();
//////	reconRealPatient20151111();
//////	return 0;
////	//realPatientrecon20151126_alreadyinGPU(); //It has problem
////	//reconBigPatient20151126_alreadyinGPU();
////
////	//reconBigPatient20151123_v2(); //OS-SART reconstruction of the patient.
////
////
////	//CGreconBigPatient20151124(); //CG reconstruction of the patient
////	//testForbildHeadReconWithOffset();
////
////	//CG_reconRealPatient20151124();
////
//////	int SIZE1 = 2620/4;
////
//////	TestProjectionTimeWithNumOfGPUs(SIZE1);
//////
//////	do{
//////		std::cout<<"Number of Titan X is ";
//////		std::cin>>SIZE1;
//////		if(SIZE1==99999)
//////			break;
//////
//////	}while(SIZE1!=99999);
//////	std::cout<<"End";
//////
////
////
////
//////
//////	std::string sheppInput = "sheppLoganReal.raw";
//////	std::string sheppOutput = "newSheppRecon";
//////	reconDEMO20151115(sheppInput,sheppOutput);
//////	std::string forbildInput = "forbildReal.raw";
//////	std::string forbildOutput = "ForbildRecon";
//////	reconDEMO20151115(forbildInput,forbildOutput);
//////	std::string patientInput = "patientReal.raw";
//////	std::string patientOutput = "patientRecon";
//////	reconDEMO20151115(patientInput,patientOutput);
////
//////  TestProjectionTime20151117();
//
////}
//
