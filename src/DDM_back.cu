

template<typename T>
__global__ void DDM_ED_bproj_ker(const T* proj, T* img,
	const T S2O, const T O2D, const T objSizeX, const T objSizeY,
	const T detSize,
	const T detCntIdx, const int XN, const int YN, const int DN, const int PN,
	const T dd, const T dx, const T dy, const T hfXN, const T hfYN,
	const T S2D,
	const T* angs)
{
	const int ii = threadIdx.x + blockIdx.x * blockDim.x;
	const int jj = threadIdx.y + blockIdx.y * blockDim.y;
	if (ii < XN && jj < YN)
	{
		T summ = 0;
		for (int angIdx = 0; angIdx < PN; angIdx++)
		{
			T curang = angs[angIdx];
			T cosT = cos(curang);
			T sinT = sin(curang);
			T lefPx(0);
			T lefPy(0);
			T rghPx(0);
			T rghPy(0);
			T temp(0);
			if ((curang > PI * 0.25 && curang <= PI * 0.75) || (curang >= PI * 1.25 && curang < PI * 1.75))
			{
				lefPx = (ii - hfXN) * dx;
				lefPy = (jj - hfYN + 0.5) * dy;
				rghPx = (ii - hfXN + 1.0) * dx;
				rghPy = (jj - hfYN + 0.5) * dy;
			}
			else
			{
				lefPx = (ii - hfXN + 0.5) * dx;
				lefPy = (jj - hfYN + 1.0) * dy;
				rghPx = (ii - hfXN + 0.5) * dx;
				rghPy = (jj - hfYN) * dy;
			}

			T initObjX1 = lefPx * cosT + lefPy * sinT;
			T initObjY1 = -lefPx * sinT + lefPy * cosT;
			T initObjX2 = rghPx * cosT + rghPy * sinT;
			T initObjY2 = -rghPx * sinT + rghPy * cosT;

			T objYdetPosMin = initObjY1 * S2D / (S2O - initObjX1);
			T objYdetPosMax = initObjY2 * S2D / (S2O - initObjX2);

			if (objYdetPosMax < objYdetPosMin)
			{
				temp = objYdetPosMax;
				objYdetPosMax = objYdetPosMin;
				objYdetPosMin = temp;

			}
			int minDetIdx = floor(objYdetPosMax / (-dd) + detCntIdx);
			int maxDetIdx = ceil(objYdetPosMin / (-dd) + detCntIdx);

			if (minDetIdx > DN)
			{
				continue;
			}
			if (maxDetIdx < 0)
			{
				continue;
			}

			T objYaxisPosMin = initObjX1 * initObjY1 / (S2O - initObjX1) + initObjY1; //pixel端点在Y轴上的投影;
			T objYaxisPosMax = initObjX2 * initObjY2 / (S2O - initObjX2) + initObjY2;
			if (objYaxisPosMax < objYaxisPosMin)
			{
				temp = objYaxisPosMax;
				objYaxisPosMax = objYaxisPosMin;
				objYaxisPosMin = temp;

			}
			T objYaxisLength = abs(objYaxisPosMax - objYaxisPosMin);

			if (minDetIdx < 0)
			{
				minDetIdx = 0;
			}
			if (maxDetIdx >= DN)
			{
				maxDetIdx = DN;
			}

			for (int detIdx = minDetIdx; detIdx < maxDetIdx; ++detIdx)
			{
				T maxDetPos = (-(detIdx - detCntIdx) * dd) * S2O / S2D;
				T minDetPos = (-(detIdx + 1.0 - detCntIdx) * dd) * S2O / S2D;
				if (maxDetPos < minDetPos)
				{
					temp = minDetPos;
					minDetPos = maxDetPos;
					maxDetPos = temp;

				}
				T s = (-(detIdx + 0.5 - detCntIdx) * dd) * S2O / S2D;

				T ll = sqrt(S2O * S2O + s * s);

				T cosAng = abs(S2O / ll);
				summ += proj[angIdx * DN + detIdx] * intersectLength_device<T>(objYaxisPosMin, objYaxisPosMax, minDetPos, maxDetPos) / (objYaxisLength * cosAng);

			}

			if ((curang > PI * 0.25 && curang <= PI * 0.75) || (curang >= PI * 1.25 && curang < PI * 1.75))
			{
				summ *= dy;
			}
			else
			{
				summ *= dx;
			}

		}
		img[jj * XN + ii] = summ;
	}
}


template<typename T>
void DDM_ED_bproj_GPU_template(const T* proj, T* img,
	const T S2O, const T O2D, const T objSizeX, const T objSizeY,
	const T detSize,
	const T detCntIdx, const int XN, const int YN, const int DN, const int PN,
	const T dd, const T dx, const T dy, const T hfXN, const T hfYN,
	const T S2D,
	const T* angs, const dim3 blk, const dim3 gid)
{
	DDM_ED_bproj_ker<T> << <gid, blk >> >(proj, img, S2O, O2D, objSizeX, objSizeY,
		detSize, detCntIdx, XN, YN, DN, PN, dd, dx, dy, hfXN, hfYN, S2D, angs);
}


void DDM_ED_bproj_GPU(const double* proj, double* img,
	const double S2O, const double O2D, const double objSizeX, const double objSizeY,
	const double detSize, const double detCntIdx, const int XN, const int YN, const int DN, const int PN,
	const double dd, const double dx, const double dy, const double hfXN, const double hfYN,
	const double S2D,
	const double* angs, const dim3 blk, const dim3 gid)
{
	DDM_ED_bproj_GPU_template<double>(proj, img, S2O, O2D, objSizeX, objSizeY,
		detSize, detCntIdx, XN, YN, DN, PN,
		dd, dx, dy, hfXN, hfYN, S2D, angs, blk, gid);
}


void DDM_ED_bproj_GPU(const float* proj, float* img,
	const float S2O, const float O2D, const float objSizeX, const float objSizeY,
	const float detSize, const float detCntIdx, const int XN, const int YN, const int DN, const int PN,
	const float dd, const float dx, const float dy, const float hfXN, const float hfYN,
	const float S2D,
	const float* angs, const dim3 blk, const dim3 gid)
{
	DDM_ED_bproj_GPU_template<float>(proj, img, S2O, O2D, objSizeX, objSizeY,
		detSize, detCntIdx, XN, YN, DN, PN,
		dd, dx, dy, hfXN, hfYN, S2D, angs, blk, gid);
}


template<typename T>
void DDM_ED_bproj_GPU_template(const thrust::device_vector<T>& proj, thrust::device_vector<T>& img,
	const T S2O, const T O2D, const T objSizeX, const T objSizeY,
	const T detSize, const T detCntIdx, const int XN, const int YN, const int DN, const int PN,
	const T dd, const T dx, const T dy, const T hfXN, const T hfYN,
	const T S2D,
	const thrust::device_vector<T>& angs, const dim3 blk, const dim3 gid)
{
	const T* pProj = thrust::raw_pointer_cast(&proj[0]);
	T* pImg = thrust::raw_pointer_cast(&img[0]);
	const T* pAngs = thrust::raw_pointer_cast(&angs[0]);

	DDM_ED_bproj_GPU_template<T>(pProj, pImg, S2O, O2D, objSizeX, objSizeY,
		detSize, detCntIdx, XN, YN, DN, PN,
		dd, dx, dy, hfXN, hfYN, S2D, pAngs, blk, gid);
}

void DDM_ED_bproj_GPU(const thrust::device_vector<float>& proj, thrust::device_vector<float>& img,
	const float S2O, const float O2D, const float objSizeX, const float objSizeY,
	const float detSize, const float detCntIdx, const int XN, const int YN, const int DN, const int PN,
	const float dd, const float dx, const float dy, const float hfXN, const float hfYN,
	const float S2D,
	const thrust::device_vector<float>& angs, const dim3 blk, const dim3 gid)
{
	DDM_ED_bproj_GPU_template<float>(proj, img,
		S2O, O2D, objSizeX, objSizeY, detSize, detCntIdx,
		XN, YN, DN, PN, dd, dx, dy, hfXN, hfYN, S2D,
		angs, blk, gid);
}

void DDM_ED_bproj_GPU(const thrust::device_vector<double>& proj, thrust::device_vector<double>& img,
	const double S2O, const double O2D, const double objSizeX, const float objSizeY,
	const double detSize, const double detCntIdx, const int XN, const int YN, const int DN, const int PN,
	const double dd, const double dx, const double dy, const double hfXN, const double hfYN,
	const double S2D,
	const thrust::device_vector<double>& angs, const dim3 blk, const dim3 gid)
{
	DDM_ED_bproj_GPU_template<double>(proj, img,
		S2O, O2D, objSizeX, objSizeY, detSize, detCntIdx,
		XN, YN, DN, PN, dd, dx, dy, hfXN, hfYN, S2D,
		angs, blk, gid);
}



template<typename T>
__global__ void DDM3D_ED_bproj_GPU_template(const T* proj, T* vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detSizeU, const T detSizeV,
	const T detCntIdU, const T detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const T ddu, const T ddv, const T dx, const T dy, const T dz,
	const T hfXN, const T hfYN, const T hfZN, const T S2D,
	const T* angs)
{
	const int ii = threadIdx.x + blockIdx.x * blockDim.x;
	const int jj = threadIdx.y + blockIdx.y * blockDim.y;
	const int kk = threadIdx.z + blockIdx.z * blockDim.z;
	if (ii < XN && jj < YN  && kk < ZN)
	{
		T summ = 0;
		T temp = 0;
		for (int angIdx = 0; angIdx != PN; ++angIdx)
		{
			T curang = angs[angIdx];
			T cosT = cos(curang);
			T sinT = sin(curang);
			T Px = (ii - hfXN + 0.5) * dx;
			T Py = (jj - hfYN + 0.5) * dy;
			T Pz = (kk - hfZN + 0.5) * dz;

			T lefPx(0);
			T lefPy(0);
			//T lefPz(0);

			T rghPx(0);
			T rghPy(0);
			//T rghPz(0);

			T uppPx(0);
			T uppPy(0);
			T uppPz(0);

			T dowPx(0);
			T dowPy(0);
			T dowPz(0);

			if ((curang > PI * 0.25 && curang <= PI * 0.75) || (curang >= PI * 1.25 && curang < PI * 1.75))
			{
				lefPx = Px - 0.5 * dx;
				lefPy = Py;
				//lefPz = Pz;

				rghPx = Px + 0.5 * dx;
				rghPy = Py;
				//rghPz = Pz;

				uppPx = Px;
				uppPy = Py;
				uppPz = Pz + 0.5 * dz;

				dowPx = Px;
				dowPy = Py;
				dowPz = Pz - 0.5 * dz;

			}
			else
			{
				lefPx = Px;
				lefPy = Py - 0.5 * dy;
				//lefPz = Pz;

				rghPx = Px;
				rghPy = Py + 0.5 * dy;
				//rghPz = Pz;

				uppPx = Px;
				uppPy = Py;
				uppPz = Pz + 0.5 * dz;

				dowPx = Px;
				dowPy = Py;
				dowPz = Pz - 0.5 * dz;
			}

			T initObjlefx = lefPx * cosT + lefPy * sinT;
			T initObjlefy = -lefPx * sinT + lefPy * cosT;
			//T initObjlefz = lefPz;

			T initObjrghx = rghPx * cosT + rghPy * sinT;
			T initObjrghy = -rghPx * sinT + rghPy * cosT;
			//T initObjrghz = rghPz;

			T initObjuppx = uppPx * cosT + uppPy * sinT;
			//T initObjuppy = -uppPx * sinT + uppPy * cosT;
			T initObjuppz = uppPz;

			T initObjdowx = dowPx* cosT + dowPy * sinT;
			//T initObjdowy = -dowPx * sinT + dowPy * cosT;
			T initObjdowz = dowPz;

			T objYdetPosUMin = initObjlefy * S2D / (S2O - initObjlefx);
			T objYdetPosUMax = initObjrghy * S2D / (S2O - initObjrghx);
			T objYdetPosVMin = initObjdowz * S2D / (S2O - initObjdowx);
			T objYdetPosVMax = initObjuppz * S2D / (S2O - initObjuppx);

			if (objYdetPosUMin > objYdetPosUMax)
			{
				temp = objYdetPosUMin;
				objYdetPosUMin = objYdetPosUMax;
				objYdetPosUMax = temp;
				//std::swap(objYdetPosUMin, objYdetPosUMax);
			}
			if (objYdetPosVMin > objYdetPosVMax)
			{
				temp = objYdetPosVMin;
				objYdetPosVMin = objYdetPosVMax;
				objYdetPosVMax = temp;

			}
			int minDetUIdx = floor(objYdetPosUMin / ddu + detCntIdU) - 1;
			int maxDetUIdx = ceil(objYdetPosUMax / ddu + detCntIdU) + 1;
			int minDetVIdx = floor(objYdetPosVMin / ddv + detCntIdV) - 1;
			int maxDetVIdx = ceil(objYdetPosVMax / ddv + detCntIdV) + 1;

			if (minDetUIdx > DNU)
			{
				continue;
			}
			if (maxDetUIdx < 0)
			{
				continue;
			}
			if (minDetVIdx > DNV)
			{
				continue;
			}
			if (maxDetVIdx < 0)
			{
				continue;
			}

			T objYOZLength = objYdetPosUMax - objYdetPosUMin;
			T objYOZHeight = objYdetPosVMax - objYdetPosVMin;


			if (minDetUIdx < 0)
			{
				minDetUIdx = 0;
			}
			if (maxDetUIdx > DNU)
			{
				maxDetUIdx = DNU;
			}
			if (minDetVIdx < 0)
			{
				minDetVIdx = 0;
			}
			if (maxDetVIdx > DNV)
			{
				maxDetVIdx = DNV;
			}


			for (int detIdU = minDetUIdx; detIdU < maxDetUIdx; ++detIdU)
			{
				T minDetUPos = (detIdU - detCntIdU - 0.5) * ddu;// *S2O / S2D;
				T maxDetUPos = (detIdU - detCntIdU + 0.5) * ddu;// *S2O / S2D;

				T ll = intersectLength_device<T>(objYdetPosUMin, objYdetPosUMax, minDetUPos, maxDetUPos);
				if (ll > 0)
				{
					for (int detIdV = minDetVIdx; detIdV < maxDetVIdx; ++detIdV)
					{
						T minDetVPos = (detIdV - detCntIdV - 0.5) * ddv;// *S2O / S2D;
						T maxDetVPos = (detIdV - detCntIdV + 0.5) * ddv;// *S2O / S2D;

						T DU = (detIdU - detCntIdU) * ddu;
						T DV = (detIdV - detCntIdV) * ddv;
						T cosAlphacosGamma = S2D / sqrt(DU*DU + DV*DV + S2D*S2D);
						T mm = intersectLength_device<T>(objYdetPosVMin, objYdetPosVMax, minDetVPos, maxDetVPos);
						if (mm > 0)
						{
							summ += (proj[(angIdx * DNV + detIdV) * DNU + detIdU] * ll * mm / (objYOZLength * objYOZHeight * cosAlphacosGamma) * dx);
						}
						else
						{
							summ += 0;
						}

					}
				}

			}


		}
		vol[(kk * YN + jj) * XN + ii] = summ;
	}
}

void DDM3D_ED_bproj_GPU(const double* proj, double* vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeU, const double detSizeV,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double ddu, const double ddv, const double dx, const double dy, const double dz,
	const double hfXN, const double hfYN, const double hfZN, const double S2D,
	const double* angs, const dim3 blk, const dim3 gid)
{
	DDM3D_ED_bproj_GPU_template<double> << <gid, blk >> >(proj, vol,
		S2O, O2D, objSizeX, objSizeY, objSizeZ, detSizeU, detSizeV,
		detCntIdU, detCntIdV, XN, YN, ZN, DNU, DNV, PN,
		ddu, ddv, dx, dy, dz, hfXN, hfYN, hfZN, S2D, angs);
}

void DDM3D_ED_bproj_GPU(const float* proj, float* vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeU, const float detSizeV,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float ddu, const float ddv, const float dx, const float dy, const float dz,
	const float hfXN, const float hfYN, const float hfZN, const float S2D,
	const float* angs, const dim3 blk, const dim3 gid)
{
	DDM3D_ED_bproj_GPU_template<float> << <gid, blk >> >(proj, vol,
		S2O, O2D, objSizeX, objSizeY, objSizeZ, detSizeU, detSizeV,
		detCntIdU, detCntIdV, XN, YN, ZN, DNU, DNV, PN,
		ddu, ddv, dx, dy, dz, hfXN, hfYN, hfZN, S2D, angs);
}


template<typename T>
void DDM3D_ED_bproj_GPU_template(const thrust::device_vector<T>& proj, thrust::device_vector<T>& vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detSizeU, const T detSizeV,
	const T detCntIdU, const T detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const T ddu, const T ddv, const T dx, const T dy, const T dz,
	const T hfXN, const T hfYN, const T hfZN, const T S2D,
	const thrust::device_vector<T>& angs, const dim3 blk, const dim3 gid)
{
	const T* pProj = thrust::raw_pointer_cast(&proj[0]);
	T* pVol = thrust::raw_pointer_cast(&vol[0]);
	const T* pAngs = thrust::raw_pointer_cast(&angs[0]);

	DDM3D_ED_bproj_GPU_template<T> << <gid, blk >> >(pProj, pVol,
		S2O, O2D, objSizeX, objSizeY, objSizeZ, detSizeU, detSizeV,
		detCntIdU, detCntIdV, XN, YN, ZN, DNU, DNV, PN,
		ddu, ddv, dx, dy, dz, hfXN, hfYN, hfZN, S2D, pAngs);
}

void DDM3D_ED_bproj_GPU(const thrust::device_vector<double>& proj, thrust::device_vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeU, const double detSizeV,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double ddu, const double ddv, const double dx, const double dy, const double dz,
	const double hfXN, const double hfYN, const double hfZN, const double S2D,
	const thrust::device_vector<double>& angs, const dim3 blk, const dim3 gid)
{
	DDM3D_ED_bproj_GPU_template<double>(proj, vol, S2O, O2D,
		objSizeX, objSizeY, objSizeZ, detSizeU, detSizeV, detCntIdU, detCntIdV,
		XN, YN, ZN, DNU, DNV, PN, ddu, ddv, dx, dy, dz, hfXN, hfYN, hfZN, S2D,
		angs, blk, gid);
}


void DDM3D_ED_bproj_GPU(const thrust::device_vector<float>& proj, thrust::device_vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeU, const float detSizeV,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float ddu, const float ddv, const float dx, const float dy, const float dz,
	const float hfXN, const float hfYN, const float hfZN, const float S2D,
	const thrust::device_vector<float>& angs, const dim3 blk, const dim3 gid)
{
	DDM3D_ED_bproj_GPU_template<float>(proj, vol, S2O, O2D,
		objSizeX, objSizeY, objSizeZ, detSizeU, detSizeV, detCntIdU, detCntIdV,
		XN, YN, ZN, DNU, DNV, PN, ddu, ddv, dx, dy, dz, hfXN, hfYN, hfZN, S2D,
		angs, blk, gid);
}


template<typename T>
__global__ void DDM3D_EA_helical_bproj_GPU_template_ker(const T* proj, T* vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detArc, const T detSizeV,
	const T detCntIdU, const T detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const T dbeta, const T ddv, const T dx, const T dy, const T dz,
	const T hfXN, const T hfYN, const T hfZN, const T S2D,
	const T initZPos, const T pitch, const T* angs)
{
	const int ii = threadIdx.x + blockIdx.x * blockDim.x;
	const int jj = threadIdx.y + blockIdx.y * blockDim.y;
	const int kk = threadIdx.z + blockIdx.z * blockDim.z;

	if (ii < XN && jj < YN && kk < ZN)
	{
		T summ = 0;
		for (int angIdx = 0; angIdx < PN; ++angIdx)
		{
			T curang = angs[angIdx];
			T cosT = cos(curang);
			T sinT = sin(curang);

			T Px = (ii - hfXN + 0.5) * dx;
			T Py = (jj - hfYN + 0.5) * dy;
			T Pz = (kk - hfZN + 0.5) * dz;

			T lefPx(0);
			T lefPy(0);

			T rghPx(0);
			T rghPy(0);

			T uppPx(0);
			T uppPy(0);
			T uppPz(0);

			T dowPx(0);
			T dowPy(0);
			T dowPz(0);

			if ((curang > PI * 0.25 && curang <= PI * 0.75) || (curang >= PI * 1.25 && curang < PI * 1.75))
			{
				lefPx = Px - 0.5 * dx;
				lefPy = Py;

				rghPx = Px + 0.5 * dx;
				rghPy = Py;

				uppPx = Px;
				uppPy = Py;
				uppPz = Pz + 0.5 * dz;

				dowPx = Px;
				dowPy = Py;
				dowPz = Pz - 0.5 * dz;

			}
			else
			{
				lefPx = Px;
				lefPy = Py - 0.5 * dy;

				rghPx = Px;
				rghPy = Py + 0.5 * dy;

				uppPx = Px;
				uppPy = Py;
				uppPz = Pz + 0.5 * dz;

				dowPx = Px;
				dowPy = Py;
				dowPz = Pz - 0.5 * dz;

			}

			T initObjLefx = lefPx * cosT + lefPy * sinT;
			T initObjlefy = -lefPx * sinT + lefPy * cosT;

			T initObjrghx = rghPx * cosT + rghPy * sinT;
			T initObjrghy = -rghPx * sinT + rghPy * cosT;

			T initObjuppx = uppPx * cosT + uppPy * sinT;

			T initObjuppz = uppPz - (initZPos + pitch * angIdx);

			T initObjdowx = dowPx * cosT + dowPy * sinT;

			T initObjdowz = dowPz - (initZPos + pitch * angIdx);

			T objYdetPosUMinbeta = atan(initObjLefx / (S2O - initObjlefy));
			T objYdetPosUMaxbeta = atan(initObjrghx / (S2O - initObjrghy));
			if (objYdetPosUMinbeta > objYdetPosUMaxbeta)
			{
				T tt = objYdetPosUMaxbeta;
				objYdetPosUMaxbeta = objYdetPosUMinbeta;
				objYdetPosUMinbeta = tt;
			}

			T objYdetPosVMin = initObjdowz * S2D / (S2O - initObjdowx);
			T objYdetPosVMax = initObjuppz * S2D / (S2O - initObjuppx);
			if (objYdetPosVMin > objYdetPosVMax)
			{
				T tt = objYdetPosVMax;
				objYdetPosVMax = objYdetPosVMin;
				objYdetPosVMin = tt;
			}

			int minDetUIdx = floor(objYdetPosUMinbeta / dbeta + detCntIdU) - 1;
			int maxDetUIdx = ceil(objYdetPosUMaxbeta / dbeta + detCntIdU) + 1;

			int minDetVIdx = floor(objYdetPosVMin / ddv + detCntIdV) - 1;
			int maxDetVIdx = ceil(objYdetPosVMax / ddv + detCntIdV) + 1;

			if (minDetUIdx > DNU)
			{
				continue;
			}

			if (maxDetUIdx < 0)
			{
				continue;
			}

			if (minDetVIdx > DNV)
			{
				continue;
			}
			if (maxDetVIdx < 0)
			{
				continue;
			}

			T objYOZLengthbeta = abs(objYdetPosUMaxbeta - objYdetPosUMinbeta);//žÄ;
			T objYOZHeight = abs(objYdetPosVMax - objYdetPosVMin);

			if (minDetUIdx < 0)
			{
				minDetUIdx = 0;
			}
			if (maxDetUIdx > DNU)
			{
				maxDetUIdx = DNU;
			}
			if (minDetVIdx < 0)
			{
				minDetVIdx = 0;
			}
			if (maxDetVIdx > DNV)
			{
				maxDetVIdx = DNV;
			}

			for (int detIdU = minDetUIdx; detIdU < maxDetUIdx; ++detIdU)
			{
				T minDetUPosbeta = (detIdU - detCntIdU - 0.5) * dbeta;
				T maxDetUPosbeta = (detIdU - detCntIdU + 0.5) * dbeta;

				T ll = intersectLength_device<T>(objYdetPosUMinbeta, objYdetPosUMaxbeta, minDetUPosbeta, maxDetUPosbeta);
				if (ll > 0)
				{
					for (int detIdV = minDetVIdx; detIdV < maxDetVIdx; ++detIdV)
					{
						T minDetVPos = (detIdV - detCntIdV - 0.5) * ddv;
						T maxDetVPos = (detIdV - detCntIdV + 0.5) * ddv;

						T DU = abs(sin((detIdU - detCntIdU) * dbeta) * S2D);
						T DV = (detIdV - detCntIdV) * ddv;
						T cosAlphacosGamma = S2D / sqrt(DU * DU + DV * DV + S2D * S2D);
						T mm = intersectLength_device<T>(objYdetPosVMin, objYdetPosVMax, minDetVPos, maxDetVPos);
						if (mm > 0)
						{
							summ += (proj[(angIdx * DNV + detIdV) * DNU + detIdU] * ll * mm / (objYOZLengthbeta * objYOZHeight * cosAlphacosGamma) * dx);

						}
						else
						{
							summ += 0;
						}
					}
				}
			}

		}
		vol[(kk * YN + jj) * XN + ii] = summ;
	}
}


void DDM3D_EA_helical_bproj(const float* proj, float* vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float dbeta, const float ddv, const float dx, const float dy, const float dz,
	const float hfXN, const float hfYN, const float hfZN, const float S2D,
	const float initZPos, const float pitch,
	const float*  angs, const dim3 blk, const dim3 gid)
{
	DDM3D_EA_helical_bproj_GPU_template_ker<float> << <gid, blk >> >(proj, vol, S2O, O2D,
		objSizeX, objSizeY, objSizeZ, detArc, detSizeH, detCntIdU, detCntIdV,
		XN, YN, ZN, DNU, DNV, PN, dbeta, ddv, dx, dy, dz, hfXN, hfYN, hfZN, S2D, initZPos, pitch, angs);
}

void DDM3D_EA_helical_bproj(const double* proj, double* vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double dbeta, const double ddv, const double dx, const double dy, const double dz,
	const double hfXN, const double hfYN, const double hfZN, const double S2D,
	const double initZPos, const double pitch,
	const double*  angs, const dim3 blk, const dim3 gid)
{
	DDM3D_EA_helical_bproj_GPU_template_ker<double> << <gid, blk >> >(proj, vol, S2O, O2D,
		objSizeX, objSizeY, objSizeZ, detArc, detSizeH, detCntIdU, detCntIdV,
		XN, YN, ZN, DNU, DNV, PN, dbeta, ddv, dx, dy, dz, hfXN, hfYN, hfZN, S2D, initZPos, pitch, angs);
}




template<typename T>
__global__ void DDM3D_EA_helical_bproj_GPU_template(const T* proj, T* vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detArc, const T detSizeV,
	const T detCntIdU, const T detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const T dbeta, const T ddv, const T dx, const T dy, const T dz,
	const T hfXN, const T hfYN, const T hfZN, const T S2D,
	const T* zShifts, const T* angs)
{
	const int ii = threadIdx.x + blockIdx.x * blockDim.x;
	const int jj = threadIdx.y + blockIdx.y * blockDim.y;
	const int kk = threadIdx.z + blockIdx.z * blockDim.z;

	if (ii < XN && jj < YN && kk < ZN)
	{
		T summ = 0;
		for (int angIdx = 0; angIdx < PN; ++angIdx)
		{
			T curang = angs[angIdx];
			T cosT = cos(curang);
			T sinT = sin(curang);

			T Px = (ii - hfXN + 0.5) * dx;
			T Py = (jj - hfYN + 0.5) * dy;
			T Pz = (kk - hfZN + 0.5) * dz;

			T lefPx(0);
			T lefPy(0);

			T rghPx(0);
			T rghPy(0);

			T uppPx(0);
			T uppPy(0);
			T uppPz(0);

			T dowPx(0);
			T dowPy(0);
			T dowPz(0);

			if ((curang > PI * 0.25 && curang <= PI * 0.75) || (curang >= PI * 1.25 && curang < PI * 1.75))
			{
				lefPx = Px - 0.5 * dx;
				lefPy = Py;

				rghPx = Px + 0.5 * dx;
				rghPy = Py;

				uppPx = Px;
				uppPy = Py;
				uppPz = Pz + 0.5 * dz;

				dowPx = Px;
				dowPy = Py;
				dowPz = Pz - 0.5 * dz;

			}
			else
			{
				lefPx = Px;
				lefPy = Py - 0.5 * dy;

				rghPx = Px;
				rghPy = Py + 0.5 * dy;

				uppPx = Px;
				uppPy = Py;
				uppPz = Pz + 0.5 * dz;

				dowPx = Px;
				dowPy = Py;
				dowPz = Pz - 0.5 * dz;

			}

			T initObjLefx = lefPx * cosT + lefPy * sinT;
			T initObjlefy = -lefPx * sinT + lefPy * cosT;

			T initObjrghx = rghPx * cosT + rghPy * sinT;
			T initObjrghy = -rghPx * sinT + rghPy * cosT;

			T initObjuppx = uppPx * cosT + uppPy * sinT;

			T initObjuppz = uppPz - zShifts[angIdx];

			T initObjdowx = dowPx * cosT + dowPy * sinT;

			T initObjdowz = dowPz - zShifts[angIdx];

			T objYdetPosUMinbeta = atan(initObjLefx / (S2O - initObjlefy));
			T objYdetPosUMaxbeta = atan(initObjrghx / (S2O - initObjrghy));
			if (objYdetPosUMinbeta > objYdetPosUMaxbeta)
			{
				T tt = objYdetPosUMaxbeta;
				objYdetPosUMaxbeta = objYdetPosUMinbeta;
				objYdetPosUMinbeta = tt;
			}

			T objYdetPosVMin = initObjdowz * S2D / (S2O - initObjdowx);
			T objYdetPosVMax = initObjuppz * S2D / (S2O - initObjuppx);
			if (objYdetPosVMin > objYdetPosVMax)
			{
				T tt = objYdetPosVMax;
				objYdetPosVMax = objYdetPosVMin;
				objYdetPosVMin = tt;
			}

			int minDetUIdx = floor(objYdetPosUMinbeta / dbeta + detCntIdU) - 1;
			int maxDetUIdx = ceil(objYdetPosUMaxbeta / dbeta + detCntIdU) + 1;

			int minDetVIdx = floor(objYdetPosVMin / ddv + detCntIdV) - 1;
			int maxDetVIdx = ceil(objYdetPosVMax / ddv + detCntIdV) + 1;

			if (minDetUIdx > DNU)
			{
				continue;
			}

			if (maxDetUIdx < 0)
			{
				continue;
			}

			if (minDetVIdx > DNV)
			{
				continue;
			}
			if (maxDetVIdx < 0)
			{
				continue;
			}

			T objYOZLengthbeta = abs(objYdetPosUMaxbeta - objYdetPosUMinbeta);//žÄ;
			T objYOZHeight = abs(objYdetPosVMax - objYdetPosVMin);

			if (minDetUIdx < 0)
			{
				minDetUIdx = 0;
			}
			if (maxDetUIdx > DNU)
			{
				maxDetUIdx = DNU;
			}
			if (minDetVIdx < 0)
			{
				minDetVIdx = 0;
			}
			if (maxDetVIdx > DNV)
			{
				maxDetVIdx = DNV;
			}

			for (int detIdU = minDetUIdx; detIdU < maxDetUIdx; ++detIdU)
			{
				T minDetUPosbeta = (detIdU - detCntIdU - 0.5) * dbeta;
				T maxDetUPosbeta = (detIdU - detCntIdU + 0.5) * dbeta;

				T ll = intersectLength_device<T>(objYdetPosUMinbeta, objYdetPosUMaxbeta, minDetUPosbeta, maxDetUPosbeta);
				if (ll > 0)
				{
					for (int detIdV = minDetVIdx; detIdV < maxDetVIdx; ++detIdV)
					{
						T minDetVPos = (detIdV - detCntIdV - 0.5) * ddv;
						T maxDetVPos = (detIdV - detCntIdV + 0.5) * ddv;

						T DU = abs(sin((detIdU - detCntIdU) * dbeta) * S2D);
						T DV = (detIdV - detCntIdV) * ddv;
						T cosAlphacosGamma = S2D / sqrt(DU * DU + DV * DV + S2D * S2D);
						T mm = intersectLength_device<T>(objYdetPosVMin, objYdetPosVMax, minDetVPos, maxDetVPos);
						if (mm > 0)
						{
							summ += (proj[(angIdx * DNV + detIdV) * DNU + detIdU] * ll * mm / (objYOZLengthbeta * objYOZHeight * cosAlphacosGamma) * dx);

						}
						else
						{
							summ += 0;
						}
					}
				}
			}

		}
		vol[(kk * YN + jj) * XN + ii] = summ;
	}
}

template<typename T>
void DDM3D_EA_helical_bproj_template(const T* proj, T* vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detArc, const T detSizeH,
	const T detCntIdU, const T detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const T dbeta, const T ddv, const T dx, const T dy, const T dz,
	const T hfXN, const T hfYN, const T hfZN, const T S2D,
	const T* zShifts,
	const T*  angs, const dim3 blk, const dim3 gid)
{
	DDM3D_EA_helical_bproj_GPU_template<T> << <gid, blk >> >(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ,
		detArc, detSizeH, detCntIdU, detCntIdV, XN, YN, ZN, DNU, DNV, PN, dbeta, ddv, dx, dy, dz,
		hfXN, hfYN, hfZN, S2D, zShifts, angs);
}

void DDM3D_EA_helical_bproj(const float* proj, float* vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float dbeta, const float ddv, const float dx, const float dy, const float dz,
	const float hfXN, const float hfYN, const float hfZN, const float S2D,
	const float* zShifts,
	const float*  angs, const dim3 blk, const dim3 gid)
{
	DDM3D_EA_helical_bproj_template<float>(proj, vol, S2O, O2D,
		objSizeX, objSizeY, objSizeZ, detArc, detSizeH, detCntIdU,
		detCntIdV, XN, YN, ZN, DNU, DNV, PN, dbeta, ddv, dx, dy,
		dz, hfXN, hfYN, hfZN, S2D, zShifts, angs, blk, gid);
}
void DDM3D_EA_helical_bproj(const double* proj, double* vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double dbeta, const double ddv, const double dx, const double dy, const double dz,
	const double hfXN, const double hfYN, const double hfZN, const double S2D,
	const double* zShifts,
	const double*  angs, const dim3 blk, const dim3 gid)
{
	DDM3D_EA_helical_bproj_template<double>(proj, vol, S2O, O2D,
		objSizeX, objSizeY, objSizeZ, detArc, detSizeH, detCntIdU,
		detCntIdV, XN, YN, ZN, DNU, DNV, PN, dbeta, ddv, dx, dy,
		dz, hfXN, hfYN, hfZN, S2D, zShifts, angs, blk, gid);
}


template<typename T>
void DDM3D_EA_helical_bproj_template(const thrust::device_vector<T>& proj, thrust::device_vector<T>& vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detArc, const T detSizeH,
	const T detCntIdU, const T detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const T dbeta, const T ddv, const T dx, const T dy, const T dz,
	const T hfXN, const T hfYN, const T hfZN, const T S2D,
	const T initZPos, const T pitch,
	const thrust::device_vector<T>& angs, const dim3 blk, const dim3 gid)
{
	const T* pProj = thrust::raw_pointer_cast(&proj[0]);
	T* pVol = thrust::raw_pointer_cast(&vol[0]);
	const T* pAngs = thrust::raw_pointer_cast(&angs[0]);

	DDM3D_EA_helical_bproj_GPU_template_ker<T> << <gid, blk >> >(pProj, pVol, S2O, O2D,
		objSizeX, objSizeY, objSizeZ, detArc, detSizeH,
		detCntIdU, detCntIdV,
		XN, YN, ZN, DNU, DNV, PN,
		dbeta, ddv, dx, dy, dz,
		hfXN, hfYN, hfZN, S2D,
		initZPos, pitch, pAngs);
}


void DDM3D_EA_helical_bproj(const thrust::device_vector<float>& proj, thrust::device_vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float dbeta, const float ddv, const float dx, const float dy, const float dz,
	const float hfXN, const float hfYN, const float hfZN, const float S2D,
	const float initZPos, const float pitch,
	const thrust::device_vector<float>& angs, const dim3 blk, const dim3 gid)
{
	DDM3D_EA_helical_bproj_template<float>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ,
		detArc, detSizeH, detCntIdU, detCntIdV, XN, YN, ZN, DNU, DNV, PN, dbeta, ddv, dx, dy, dz,
		hfXN, hfYN, hfZN, S2D, initZPos, pitch, angs, blk, gid);
}
void DDM3D_EA_helical_bproj(const thrust::device_vector<double>& proj, thrust::device_vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double dbeta, const double ddv, const double dx, const double dy, const double dz,
	const double hfXN, const double hfYN, const double hfZN, const double S2D,
	const double initZPos, const double pitch,
	const thrust::device_vector<double>& angs, const dim3 blk, const dim3 gid)
{
	DDM3D_EA_helical_bproj_template<double>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ,
		detArc, detSizeH, detCntIdU, detCntIdV, XN, YN, ZN, DNU, DNV, PN, dbeta, ddv, dx, dy, dz,
		hfXN, hfYN, hfZN, S2D, initZPos, pitch, angs, blk, gid);
}

template<typename T>
void DDM3D_EA_helical_bproj_template(const thrust::device_vector<T>& proj, thrust::device_vector<T>& vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detArc, const T detSizeH,
	const T detCntIdU, const T detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const T dbeta, const T ddv, const T dx, const T dy, const T dz,
	const T hfXN, const T hfYN, const T hfZN, const T S2D,
	const thrust::device_vector<T>& zShifts,
	const thrust::device_vector<T>& angs, const dim3 blk, const dim3 gid)
{
	const T* pProj = thrust::raw_pointer_cast(&proj[0]);
	T* pVol = thrust::raw_pointer_cast(&vol[0]);
	const T* pZ = thrust::raw_pointer_cast(&zShifts[0]);
	const T* pAngs = thrust::raw_pointer_cast(&angs[0]);
	DDM3D_EA_helical_bproj_template<T>(pProj, pVol, S2O, O2D, objSizeX, objSizeY, objSizeZ,
		detArc, detSizeH, detCntIdU, detCntIdV, XN, YN, ZN, DNU, DNV, PN,
		dbeta, ddv, dx, dy, dz, hfXN, hfYN, hfZN, S2D, pZ, pAngs, blk, gid);
}

void DDM3D_EA_helical_bproj(const thrust::device_vector<float>& proj, thrust::device_vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detArc, const float detSizeH,
	const float detCntIdU, const float detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const float dbeta, const float ddv, const float dx, const float dy, const float dz,
	const float hfXN, const float hfYN, const float hfZN, const float S2D,
	const thrust::device_vector<float>& zShifts,
	const thrust::device_vector<float>& angs, const dim3 blk, const dim3 gid)
{
	DDM3D_EA_helical_bproj_template<float>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ,
		detArc, detSizeH, detCntIdU, detCntIdV, XN, YN, ZN, DNU, DNV, PN, dbeta, ddv, dx, dy, dz,
		hfXN, hfYN, hfZN, S2D, zShifts, angs, blk, gid);
}
void DDM3D_EA_helical_bproj(const thrust::device_vector<double>& proj, thrust::device_vector<double>&  vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detArc, const double detSizeH,
	const double detCntIdU, const double detCntIdV,
	const int XN, const int YN, const int ZN, const int DNU, const int DNV, const int PN,
	const double dbeta, const double ddv, const double dx, const double dy, const double dz,
	const double hfXN, const double hfYN, const double hfZN, const double S2D,
	const thrust::device_vector<double>& zShifts,
	const thrust::device_vector<double>& angs, const dim3 blk, const dim3 gid)
{
	DDM3D_EA_helical_bproj_template<double>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ,
		detArc, detSizeH, detCntIdU, detCntIdV, XN, YN, ZN, DNU, DNV, PN, dbeta, ddv, dx, dy, dz,
		hfXN, hfYN, hfZN, S2D, zShifts, angs, blk, gid);
}


