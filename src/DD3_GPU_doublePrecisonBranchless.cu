//#include "DD3_GPU_doublePrecisonBranchless.h"
//
//INLINE __device__ double bilerp(int2 v0, int2 v1, int2 v2, int2 v3, double t1, double t2)
//{
//	double v0_ = __hiloint2double(v0.y, v0.x);
//	double v1_ = __hiloint2double(v1.y, v1.x);
//	double v2_ = __hiloint2double(v2.y, v2.x);
//	double v3_ = __hiloint2double(v3.y, v3.x);
//
//	double vv0 = v0_ * (1.0 - t1) + v1_ * t1;
//	double vv1 = v2_ * (1.0 - t1) + v3_ * t1;
//	return vv0 * (1 - t2) + vv1 * t2;
//}
//
//// \brief Kernel function of Double precision based branchless DD with 2D SAT
//__global__ void DD3_gpu_proj_doubleprecisionbranchless_ker(
//	cudaTextureObject_t volTex1, // volume SAT in ZXY order
//	cudaTextureObject_t volTex2, // volume SAT in ZYX order
//	double* proj, // projection data
//	double3 s, // initial source position
//	const double3* __restrict cossinZT, // bind (cosine, sine, zshift)
//	const double* __restrict xds,  //
//	const double* __restrict yds,  //
//	const double* __restrict zds,  // detector cells center positions
//	const double* __restrict bxds,
//	const double* __restrict byds,
//	const double* __restrict bzds, // detector boundary positions
//	double3 objCntIdx, // object center index
//	double dx, double dz, // pixel size in xy plane and Z direction
//	int XN, int YN, // pixel # in XY plane
//	int DNU, int DNV, // detector cell # in xy plane and Z direction
//	int PN) // view #
//{
//	int detIdV = threadIdx.x + blockIdx.x * blockDim.x;
//	int detIdU = threadIdx.y + blockIdx.y * blockDim.y;
//	int angIdx = threadIdx.z + blockIdx.z * blockDim.z;
//
//	__shared__ double _xds[BLKY];
//	__shared__ double _yds[BLKY];
//	_xds[threadIdx.y] = xds[detIdU];
//	_yds[threadIdx.y] = yds[detIdU];
//	__syncthreads();
//
//	if (detIdU < DNU && detIdV < DNV && angIdx < PN)
//	{
//		double3 dir = cossinZT[angIdx];
//		double3 cursour = make_double3(
//			s.x * dir.x - s.y * dir.y,
//			s.x * dir.y + s.y * dir.x,
//			s.z + dir.z);
//		s = cossinZT[angIdx];
//		double summ = _xds[threadIdx.y] * s.x - _yds[threadIdx.y] * s.y;
//		double obj = _xds[threadIdx.y] * s.y + _yds[threadIdx.y] * s.x;
//		double realL = bxds[detIdU];
//		double realR = byds[detIdU];
//		double realU = bxds[detIdU + 1];
//		double realD = byds[detIdU + 1]; // intersection coordinates (mm); float2 is equv to (obj1,obj2) above
//		double2 curDetL = make_double2(
//			realL * s.x - realR * s.y,
//			realL * s.y + realR * s.x);
//
//		double2 curDetR = make_double2(
//			realU * s.x - realD * s.y,
//			realU * s.y + realD * s.x);
//		double4 curDet = make_double4(summ, obj, bzds[detIdV] + s.z, bzds[detIdV + 1] + s.z); //(center x, center y, lower z, upper z)
//
//		dir = normalize(make_double3(
//			summ,
//			obj,
//			zds[detIdV] + s.z) - cursour);
//
//		summ = 0; // to accumulate projection value
//		obj = 0; // slice location (mm) along the ray tracing direction TODO: is this variable needed?
//
//		double intersectLength, intersectHeight;
//		double invdz = 1.0 / dz;
//		double invdx = 1.0 / dx;
//
//
//		double factL(1.0f); // dy/dx for (0,pi/4)
//		double factR(1.0f);
//		double factU(1.0f);
//		double factD(1.0f);
//		double constVal = 0;
//
//		int crealD, crealR, crealU, crealL;
//		int frealD, frealR, frealU, frealL;
//
//		if (abs(s.x) <= abs(s.y))
//		{
//			summ = 0;
//			// a few book keeping variables
//			assert(curDetL.x != cursour.x);
//			assert(curDetR.x != cursour.x);
//			assert(curDet.x != cursour.x);
//
//			factL = (curDetL.y - cursour.y) / (curDetL.x - cursour.x);
//			factR = (curDetR.y - cursour.y) / (curDetR.x - cursour.x);
//			factU = (curDet.w - cursour.z) / (curDet.x - cursour.x);
//			factD = (curDet.z - cursour.z) / (curDet.x - cursour.x);
//
//			assert(dir.x != 0);
//			constVal = dx * dx * dz / (abs(dir.x));
//#pragma unroll
//			for (int ii = 0; ii < XN; ii++)
//			{
//				obj = (ii - objCntIdx.x) * dx;
//
//				realL = (obj - curDetL.x) * factL + curDetL.y;
//				realR = (obj - curDetR.x) * factR + curDetR.y;
//				realU = (obj - curDet.x) * factU + curDet.w;
//				realD = (obj - curDet.x) * factD + curDet.z;
//
//				intersectLength = realR - realL;
//				intersectHeight = realU - realD;
//				assert(intersectLength != 0 && intersectHeight != 0);
//
//				// 1D LUT to address inaccuracies in texture coordinates
//				realD = realD * invdz + objCntIdx.z + 1;
//				realR = realR * invdx + objCntIdx.y + 1;
//				realU = realU * invdz + objCntIdx.z + 1;
//				realL = realL * invdx + objCntIdx.y + 1;
//
//				crealD = ceil(realD);
//				crealR = ceil(realR);
//				crealU = ceil(realU);
//				crealL = ceil(realL);
//
//				frealD = floor(realD);
//				frealR = floor(realR);
//				frealU = floor(realU);
//				frealL = floor(realL);
//
//
//				summ +=
//					(bilerp(
//						tex3D<int2>(volTex2, frealD, frealL, ii + 0.5),
//						tex3D<int2>(volTex2, frealD, crealL, ii + 0.5),
//						tex3D<int2>(volTex2, crealD, frealL, ii + 0.5),
//						tex3D<int2>(volTex2, crealD, crealL, ii + 0.5),
//						realL - frealL, realD - frealD) +
//						bilerp(
//							tex3D<int2>(volTex2, frealU, frealR, ii + 0.5),
//							tex3D<int2>(volTex2, frealU, crealR, ii + 0.5),
//							tex3D<int2>(volTex2, crealU, frealR, ii + 0.5),
//							tex3D<int2>(volTex2, crealU, crealR, ii + 0.5),
//							realR - frealR, realU - frealU) -
//						bilerp(
//							tex3D<int2>(volTex2, frealD, frealR, ii + 0.5),
//							tex3D<int2>(volTex2, frealD, crealR, ii + 0.5),
//							tex3D<int2>(volTex2, crealD, frealR, ii + 0.5),
//							tex3D<int2>(volTex2, crealD, crealR, ii + 0.5),
//							realR - frealR, realD - frealD) -
//						bilerp(
//							tex3D<int2>(volTex2, frealU, frealL, ii + 0.5),
//							tex3D<int2>(volTex2, frealU, crealL, ii + 0.5),
//							tex3D<int2>(volTex2, crealU, frealL, ii + 0.5),
//							tex3D<int2>(volTex2, crealU, crealL, ii + 0.5),
//							realL - frealL, realU - frealU)) / (intersectLength * intersectHeight);
//			}
//			__syncthreads();
//			proj[(angIdx * DNU + detIdU) * DNV + detIdV] = summ * constVal;
//		}
//		else
//		{
//			summ = 0;
//			assert(curDetL.y - cursour.y);
//			assert(curDetR.y - cursour.y);
//			assert(curDet.y - cursour.y);
//
//			factL = (curDetL.x - cursour.x) / (curDetL.y - cursour.y);
//			factR = (curDetR.x - cursour.x) / (curDetR.y - cursour.y);
//			factU = (curDet.w - cursour.z) / (curDet.y - cursour.y);
//			factD = (curDet.z - cursour.z) / (curDet.y - cursour.y);
//
//			assert(dir.y);
//			constVal = dx * dx * dz / (abs(dir.y));
//#pragma unroll
//			for (int jj = 0; jj < YN; jj++)
//			{
//				obj = (jj - objCntIdx.y) * dx;
//				realL = (obj - curDetL.y) * factL + curDetL.x;
//				realR = (obj - curDetR.y) * factR + curDetR.x;
//				realU = (obj - curDet.y) * factU + curDet.w;
//				realD = (obj - curDet.y) * factD + curDet.z;
//
//				intersectLength = realR - realL;
//				intersectHeight = realU - realD;
//				assert(intersectLength != 0 && intersectHeight != 0);
//
//				realD = realD * invdz + objCntIdx.z + 1;
//				realR = realR * invdx + objCntIdx.x + 1;
//				realU = realU * invdz + objCntIdx.z + 1;
//				realL = realL * invdx + objCntIdx.x + 1;
//
//
//				crealD = ceil(realD);
//				crealR = ceil(realR);
//				crealU = ceil(realU);
//				crealL = ceil(realL);
//
//				frealD = floor(realD);
//				frealR = floor(realR);
//				frealU = floor(realU);
//				frealL = floor(realL);
//
//				summ += (bilerp(
//					tex3D<int2>(volTex1, frealD, frealL, jj + 0.5),
//					tex3D<int2>(volTex1, frealD, crealL, jj + 0.5),
//					tex3D<int2>(volTex1, crealD, frealL, jj + 0.5),
//					tex3D<int2>(volTex1, crealD, crealL, jj + 0.5),
//					realL - frealL, realD - frealD) +
//					bilerp(
//						tex3D<int2>(volTex1, frealU, frealR, jj + 0.5),
//						tex3D<int2>(volTex1, frealU, crealR, jj + 0.5),
//						tex3D<int2>(volTex1, crealU, frealR, jj + 0.5),
//						tex3D<int2>(volTex1, crealU, crealR, jj + 0.5),
//						realR - frealR, realU - frealU) -
//					bilerp(
//						tex3D<int2>(volTex1, frealD, frealR, jj + 0.5),
//						tex3D<int2>(volTex1, frealD, crealR, jj + 0.5),
//						tex3D<int2>(volTex1, crealD, frealR, jj + 0.5),
//						tex3D<int2>(volTex1, crealD, crealR, jj + 0.5),
//						realR - frealR, realD - frealD) -
//					bilerp(
//						tex3D<int2>(volTex1, frealU, frealL, jj + 0.5),
//						tex3D<int2>(volTex1, frealU, crealL, jj + 0.5),
//						tex3D<int2>(volTex1, crealU, frealL, jj + 0.5),
//						tex3D<int2>(volTex1, crealU, crealL, jj + 0.5),
//						realL - frealL, realU - frealU)) / (intersectLength * intersectHeight);
//			}
//			__syncthreads();
//			proj[(angIdx * DNU + detIdU) * DNV + detIdV] = summ * constVal;
//		}
//
//	}
//}
//
//// \brief C interface of Double Precision Branchless DD with 2D SAT
//void DD3_gpu_proj_doubleprecisionbranchless(
//	float x0, float y0, float z0,
//	int DNU, int DNV,
//	float* xds, float* yds, float* zds,
//	float imgXCenter, float imgYCenter, float imgZCenter,
//	float* hangs, float* hzPos, int PN,
//	int XN, int YN, int ZN,
//	float* vol, float* hprj,
//	float dx, float dz,
//	byte* mask, int gpunum) {
//	//Pre compute mask.*vol;
//	maskingVolume(vol, mask, XN, YN, ZN);
//
//	float* bxds = new float[DNU + 1];
//	float* byds = new float[DNU + 1];
//	float* bzds = new float[DNV + 1];
//	DD3Boundaries(DNU + 1, xds, bxds);
//	DD3Boundaries(DNU + 1, yds, byds);
//	DD3Boundaries(DNV + 1, zds, bzds);
//
//	CUDA_CHECK_RETURN(cudaSetDevice(gpunum));
//	CUDA_CHECK_RETURN(cudaDeviceReset());
//
//	cudaStream_t streams[4];
//	CUDA_CHECK_RETURN(cudaStreamCreate(&streams[0]));
//	CUDA_CHECK_RETURN(cudaStreamCreate(&streams[1]));
//	CUDA_CHECK_RETURN(cudaStreamCreate(&streams[2]));
//	CUDA_CHECK_RETURN(cudaStreamCreate(&streams[3]));
//
//	int TOTVN = XN * YN * ZN;
//	double objCntIdxX = (XN - 1.0) * 0.5 - imgXCenter / dx;
//	double objCntIdxY = (YN - 1.0) * 0.5 - imgYCenter / dx;
//	double objCntIdxZ = (ZN - 1.0) * 0.5 - imgZCenter / dz;
//
//	thrust::device_vector<float> in(vol, vol + TOTVN); // original img volume
//	thrust::device_vector<double> in_ZXY((ZN + 1) * (XN + 1) * YN, 0); //
//	thrust::device_vector<double> in_ZYX((ZN + 1) * (YN + 1) * XN, 0); // transposed img volume
//
//	dim3 blk(64, 16, 1);
//	dim3 gid(
//		(ZN + blk.x - 1) / blk.x,
//		(XN + blk.y - 1) / blk.y,
//		(YN + blk.z - 1) / blk.z);
//	// copy to original and transposed image volume with left- and top-side boarder padding to be consistent with SAT dimensions
//	naive_copyToTwoVolumes<float, double> << <gid, blk >> > (
//		thrust::raw_pointer_cast(&in[0]),
//		thrust::raw_pointer_cast(&in_ZXY[0]),
//		thrust::raw_pointer_cast(&in_ZYX[0]), XN, YN, ZN);
//	in.clear();
//
//	thrust::device_vector<double> in_ZXY_summ1((ZN + 1) * (XN + 1) * YN, 0);
//	thrust::device_vector<int2> in_ZXY_summ((ZN + 1) * (XN + 1) * YN);
//
//
//	blk.x = 64;							blk.y = 1;		blk.z = 1;
//	gid.x = (ZN + blk.x) / blk.x;		gid.y = 1;		gid.z = 1;
//
//	dim3 blk2(64);
//	dim3 gid2((YN + blk2.x) / blk2.x);
//	dim3 blk3(64);
//	dim3 gid3((XN + blk3.x) / blk3.x);
//
//	// compute SAT for the original img volume
//	for (int jj = 0; jj != YN; ++jj)
//	{
//		// for each Y slice
//		naive_herizontalIntegral << <gid, blk, 0, streams[0] >> > (
//			thrust::raw_pointer_cast(&in_ZXY[0]) + jj * (ZN + 1) * (XN + 1),
//			thrust::raw_pointer_cast(&in_ZXY_summ1[0]) + jj * (ZN + 1) * (XN + 1), XN + 1, ZN + 1);
//		naive_verticalIntegral << <gid2, blk2, 0, streams[0] >> > (
//			thrust::raw_pointer_cast(&in_ZXY_summ1[0]) + jj * (ZN + 1) * (XN + 1),
//			thrust::raw_pointer_cast(&in_ZXY_summ[0]) + jj * (ZN + 1) * (XN + 1), XN + 1, ZN + 1);
//	}
//	in_ZXY.clear();
//	in_ZXY_summ1.clear();
//
//
//	cudaArray* d_volumeArray1 = nullptr;
//	cudaTextureObject_t texObj1;
//
//	createTextureObject<int2>(texObj1, d_volumeArray1, ZN + 1, XN + 1, YN,
//		thrust::raw_pointer_cast(&in_ZXY_summ[0]),
//		cudaMemcpyDeviceToDevice, cudaAddressModeClamp, cudaFilterModePoint,
//		cudaReadModeElementType, false);
//	in_ZXY_summ.clear();
//
//	thrust::device_vector<double> in_ZYX_summ1((ZN + 1) * (YN + 1) * XN, 0); // SAT for the transposed img volume
//	thrust::device_vector<int2> in_ZYX_summ((ZN + 1) * (YN + 1) * XN);
//	// compute SAT for the transposed img volume
//	for (int ii = 0; ii != XN; ++ii)
//	{
//		// for each X slice
//		naive_herizontalIntegral << <gid, blk, 0, streams[1] >> > (
//			thrust::raw_pointer_cast(&in_ZYX[0]) + ii * (ZN + 1) * (YN + 1),
//			thrust::raw_pointer_cast(&in_ZYX_summ1[0]) + ii * (ZN + 1) * (YN + 1), YN + 1, ZN + 1);
//		naive_verticalIntegral << <gid3, blk3, 0, streams[1] >> > (
//			thrust::raw_pointer_cast(&in_ZYX_summ1[0]) + ii * (ZN + 1) * (YN + 1),
//			thrust::raw_pointer_cast(&in_ZYX_summ[0]) + ii * (ZN + 1) * (YN + 1), YN + 1, ZN + 1);
//	}
//	in_ZYX.clear();
//	in_ZYX_summ1.clear();
//
//	cudaArray* d_volumeArray2 = nullptr;
//	cudaTextureObject_t texObj2;
//
//	createTextureObject<int2>(texObj2, d_volumeArray2, ZN + 1, YN + 1, XN,
//		thrust::raw_pointer_cast(&in_ZYX_summ[0]),
//		cudaMemcpyDeviceToDevice, cudaAddressModeClamp, cudaFilterModePoint,
//		cudaReadModeElementType, false);
//	in_ZYX_summ.clear();
//
//	thrust::device_vector<double> prj(DNU * DNV * PN, 0);
//	thrust::device_vector<double> angs(hangs, hangs + PN);
//	thrust::device_vector<double> zPos(hzPos, hzPos + PN);
//	thrust::device_vector<double> d_xds(xds, xds + DNU);
//	thrust::device_vector<double> d_yds(yds, yds + DNU);
//	thrust::device_vector<double> d_zds(zds, zds + DNV);
//	thrust::device_vector<double> d_bxds(bxds, bxds + DNU + 1);
//	thrust::device_vector<double> d_byds(byds, byds + DNU + 1);
//	thrust::device_vector<double> d_bzds(bzds, bzds + DNV + 1);
//
//
//	// constant values for DD calculation
//	thrust::device_vector<double3> cossinZT(PN);
//	thrust::transform(
//		thrust::make_zip_iterator(thrust::make_tuple(angs.begin(), zPos.begin())),
//		thrust::make_zip_iterator(thrust::make_tuple(angs.end(), zPos.end())),
//		cossinZT.begin(), CTMBIR::ConstantForBackProjection<double>(x0, y0, z0));
//
//	//precalculate all constant values in CUDA
//	dim3 blkc(64, 16, 1);
//	dim3 gidc(
//		(DNV + blkc.x) / blkc.x,
//		(DNU + blkc.y) / blkc.y,
//		(PN + blkc.z - 1) / blkc.z);
//
//
//	//Configure BLOCKs for projection
//	blk.x = BLKX; // det row index
//	blk.y = BLKY; // det col index
//	blk.z = BLKZ; // view index
//	gid.x = (DNV + blk.x - 1) / blk.x;
//	gid.y = (DNU + blk.y - 1) / blk.y;
//	gid.z = (PN + blk.z - 1) / blk.z;
//
//	//Projection kernel
//	DD3_gpu_proj_doubleprecisionbranchless_ker << <gid, blk >> > (texObj1, texObj2,
//		thrust::raw_pointer_cast(&prj[0]), make_double3(x0, y0, z0),
//		thrust::raw_pointer_cast(&cossinZT[0]),
//		thrust::raw_pointer_cast(&d_xds[0]),
//		thrust::raw_pointer_cast(&d_yds[0]),
//		thrust::raw_pointer_cast(&d_zds[0]),
//		thrust::raw_pointer_cast(&d_bxds[0]),
//		thrust::raw_pointer_cast(&d_byds[0]),
//		thrust::raw_pointer_cast(&d_bzds[0]),
//		make_double3(objCntIdxX, objCntIdxY, objCntIdxZ),
//		dx, dz, XN, YN, DNU, DNV, PN);
//	thrust::copy(prj.begin(), prj.end(), hprj);
//
//	CUDA_CHECK_RETURN(cudaDestroyTextureObject(texObj1));
//	CUDA_CHECK_RETURN(cudaDestroyTextureObject(texObj2));
//
//	destroyTextureObject(texObj1, d_volumeArray1);
//	destroyTextureObject(texObj2, d_volumeArray2);
//	prj.clear();
//	angs.clear();
//	zPos.clear();
//
//	d_xds.clear();
//	d_yds.clear();
//	d_zds.clear();
//	d_bxds.clear();
//	d_byds.clear();
//	d_bzds.clear();
//
//
//	delete[] bxds;
//	delete[] byds;
//	delete[] bzds;
//}