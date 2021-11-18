#include "GenerateSummedAreaTable.h"
template<typename Ta, typename Tb>
static __global__ void naive_copyToTwoVolumes(Ta* in_ZXY,
	Tb* out_ZXY, Tb* out_ZYX,
	int XN, int YN, int ZN) {
	const int idz = threadIdx.x + blockIdx.x * blockDim.x; // Z dimension major
	const int idx = threadIdx.y + blockIdx.y * blockDim.y;
	const int idy = threadIdx.z + blockIdx.z * blockDim.z;
	if (idx < XN && idy < YN && idz < ZN) {
		int i = (idy * XN + idx) * ZN + idz;
		int ni = (idy * (XN + 1) + (idx + 1)) * (ZN + 1) + idz + 1;
		int nj = (idx * (YN + 1) + (idy + 1)) * (ZN + 1) + idz + 1;

		out_ZXY[ni] = in_ZXY[i];
		out_ZYX[nj] = in_ZXY[i];
	}
}

template<typename Ta, typename Tb>
static __global__ void naive_herizontalIntegral(Ta* in, Tb* out, int N, int ZN) {
	const int zi = threadIdx.x + blockIdx.x * blockDim.x;
	if (zi < ZN) {
		out[zi] = in[zi];
		for (int i = 1; i < N; ++i) {
			out[i * ZN + zi] = out[(i - 1) * ZN + zi]
				+ in[i * ZN + zi];
		}
	}
}

template<typename Ta, typename Tb>
static __global__ void naive_verticalIntegral(Ta* in, Tb* out, int N, int ZN) {
	int xyi = threadIdx.x + blockIdx.x * blockDim.x;
	if (xyi < N) {
		out[xyi * ZN] = in[xyi * ZN];
		for (int ii = 1; ii < ZN; ++ii)	{
			out[xyi * ZN + ii] = out[xyi * ZN + ii - 1]
				+ in[xyi * ZN + ii];
		}
	}
}

// template specialization for double precision into int2 datatype
template<>
static __global__ void naive_verticalIntegral(double*  in, int2* out, int N, int ZN) {
	int xyi = threadIdx.x + blockIdx.x * blockDim.x;
	if (xyi < N) {
		double temp = in[xyi * ZN];
		out[xyi * ZN] = make_int2(__double2loint(temp), __double2hiint(temp));
		double temp2 = 0;
		for (int ii = 1; ii < ZN; ++ii) {
			temp2 = temp + in[xyi * ZN + ii];
			out[xyi * ZN + ii] = make_int2(__double2loint(temp2), __double2hiint(temp2));
			temp = temp2;
		}
	}
}

template<typename T>
static __global__ void verticalIntegral(T* prj, int ZN, int N) {
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < N) {
		int currentHead = idx * ZN;
		for (int ii = 1; ii < ZN; ++ii) {
			prj[currentHead + ii] += prj[currentHead + ii - 1];
		}
	}
}

template<typename T>
static __global__ void horizontalIntegral(T* prj, int DNU, int DNV, int PN) {
	int idv = threadIdx.x + blockIdx.x * blockDim.x;
	int pIdx = threadIdx.y + blockIdx.y * blockDim.y;
	if (idv < DNV && pIdx < PN)	{
		int headPtr = pIdx * DNU * DNV + idv;
		for (int ii = 1; ii < DNU; ++ii) {
			prj[headPtr + ii * DNV] += prj[headPtr + (ii - 1) * DNV];
		}
	}
}


void genSAT_fof_Volume(float* hvol,
	thrust::device_vector<float>&ZXY,
	thrust::device_vector<float>&ZYX,
	int XN, int YN, int ZN) {
	const int siz = XN * YN * ZN;
	const int nsiz_ZXY = (ZN + 1) * (XN + 1) * YN;
	const int nsiz_ZYX = (ZN + 1) * (YN + 1) * XN;
	ZXY.resize(nsiz_ZXY);
	ZYX.resize(nsiz_ZYX);

	thrust::device_vector<float> vol(hvol, hvol + siz);

	dim3 blk(64, 16, 1);
	dim3 gid(
		(ZN + blk.x - 1) / blk.x,
		(XN + blk.y - 1) / blk.y,
		(YN + blk.z - 1) / blk.z);

	naive_copyToTwoVolumes << <gid, blk >> >(
		thrust::raw_pointer_cast(&vol[0]),
		thrust::raw_pointer_cast(&ZXY[0]),
		thrust::raw_pointer_cast(&ZYX[0]),
		XN, YN, ZN);

	vol.clear();
	const int nZN = ZN + 1;
	const int nXN = XN + 1;
	const int nYN = YN + 1;

	blk.x = 32;
	blk.y = 1;
	blk.z = 1;
	gid.x = (nXN * YN + blk.x - 1) / blk.x;
	gid.y = 1;
	gid.z = 1;
	verticalIntegral << <gid, blk >> >(
		thrust::raw_pointer_cast(&ZXY[0]),
		nZN, nXN * YN);

	blk.x = 64;
	blk.y = 16;
	blk.z = 1;
	gid.x = (nZN + blk.x - 1) / blk.x;
	gid.y = (YN + blk.y - 1) / blk.y;
	gid.z = 1;

	horizontalIntegral << <gid, blk >> >(
		thrust::raw_pointer_cast(&ZXY[0]),
		nXN, nZN, YN);

	blk.x = 32;
	blk.y = 1;
	blk.z = 1;
	gid.x = (nYN * XN + blk.x - 1) / blk.x;
	gid.y = 1;
	gid.z = 1;
	verticalIntegral << <gid, blk >> >(
		thrust::raw_pointer_cast(&ZYX[0]),
		nZN, nXN * YN);

	blk.x = 64;
	blk.y = 16;
	blk.z = 1;
	gid.x = (nZN + blk.x - 1) / blk.x;
	gid.y = (XN + blk.y - 1) / blk.y;
	gid.z = 1;

	horizontalIntegral << <gid, blk >> >(
		thrust::raw_pointer_cast(&ZYX[0]),
		nYN, nZN, XN);
}

thrust::pair<thrust::device_vector<float>, int> generateSumAreaTableImpl(thrust::device_vector<float>& vol,
	int XN, int YN, int ZN) {
	const int nZN = ZN + 1;
	const int nXN = XN + 1;
	const int nYN = YN + 1;
	const int nsiz_ZXY = nZN * nXN * YN;
	const int nsiz_ZYX = nZN * nYN * XN;

	thrust::device_vector<float> res(nsiz_ZXY + nsiz_ZYX, 0);
	float* pvol = thrust::raw_pointer_cast(&vol[0]);
	float* pZXY = thrust::raw_pointer_cast(&res[0]);
	float* pZYX = thrust::raw_pointer_cast(&res[nsiz_ZXY]);

	dim3 blk(64, 16, 1);
	dim3 gid(
		(ZN + blk.x - 1) / blk.x,
		(XN + blk.y - 1) / blk.y,
		(YN + blk.z - 1) / blk.z);

	naive_copyToTwoVolumes << <gid, blk >> > (pvol, pZXY, pZYX, XN, YN, ZN);

	blk.x = 32;
	blk.y = 1;
	blk.z = 1;
	gid.x = (nXN * YN + blk.x - 1) / blk.x;
	gid.y = 1;
	gid.z = 1;
	verticalIntegral << <gid, blk >> > (pZXY, nZN, nXN * YN);

	blk.x = 64;
	blk.y = 16;
	blk.z = 1;
	gid.x = (nZN + blk.x - 1) / blk.x;
	gid.y = (YN + blk.y - 1) / blk.y;
	gid.z = 1;

	horizontalIntegral << <gid, blk >> > (pZXY, nXN, nZN, YN);

	blk.x = 32;  blk.y = 1; blk.z = 1;
	gid.x = (nYN * XN + blk.x - 1) / blk.x;
	gid.y = 1;
	gid.z = 1;
	verticalIntegral << <gid, blk >> > (pZYX, nZN, nXN * YN);

	blk.x = 64;  blk.y = 16;  blk.z = 1;
	gid.x = (nZN + blk.x - 1) / blk.x;
	gid.y = (XN + blk.y - 1) / blk.y;
	gid.z = 1;

	horizontalIntegral << <gid, blk >> > (pZYX,	nYN, nZN, XN);

	return thrust::make_pair(res, nsiz_ZXY);
}

void genSAT_fof_Volume_alreadyinGPU(
	const thrust::device_vector<float>& vol,
	thrust::device_vector<float>&ZXY,
	thrust::device_vector<float>&ZYX,
	int XN, int YN, int ZN) {
	const int nsiz_ZXY = (ZN + 1) * (XN + 1) * YN;
	const int nsiz_ZYX = (ZN + 1) * (YN + 1) * XN;
	ZXY.resize(nsiz_ZXY);
	ZYX.resize(nsiz_ZYX);

	dim3 blk(64, 16, 1);
	dim3 gid(
		(ZN + blk.x - 1) / blk.x,
		(XN + blk.y - 1) / blk.y,
		(YN + blk.z - 1) / blk.z);

	naive_copyToTwoVolumes << <gid, blk >> >(
		thrust::raw_pointer_cast(&vol[0]),
		thrust::raw_pointer_cast(&ZXY[0]),
		thrust::raw_pointer_cast(&ZYX[0]),
		XN, YN, ZN);

	const int nZN = ZN + 1;
	const int nXN = XN + 1;
	const int nYN = YN + 1;

	blk.x = 32;
	blk.y = 1;
	blk.z = 1;
	gid.x = (nXN * YN + blk.x - 1) / blk.x;
	gid.y = 1;
	gid.z = 1;
	verticalIntegral << <gid, blk >> >(
		thrust::raw_pointer_cast(&ZXY[0]),
		nZN, nXN * YN);

	blk.x = 64;
	blk.y = 16;
	blk.z = 1;
	gid.x = (nZN + blk.x - 1) / blk.x;
	gid.y = (YN + blk.y - 1) / blk.y;
	gid.z = 1;

	horizontalIntegral << <gid, blk >> >(
		thrust::raw_pointer_cast(&ZXY[0]),
		nXN, nZN, YN);

	blk.x = 32;
	blk.y = 1;
	blk.z = 1;
	gid.x = (nYN * XN + blk.x - 1) / blk.x;
	gid.y = 1;
	gid.z = 1;
	verticalIntegral << <gid, blk >> >(
		thrust::raw_pointer_cast(&ZYX[0]),
		nZN, nXN * YN);

	blk.x = 64;
	blk.y = 16;
	blk.z = 1;
	gid.x = (nZN + blk.x - 1) / blk.x;
	gid.y = (XN + blk.y - 1) / blk.y;
	gid.z = 1;

	horizontalIntegral << <gid, blk >> >(
		thrust::raw_pointer_cast(&ZYX[0]),
		nYN, nZN, XN);
}