/*
 * COPYRIGHT NOTICE
 * COPYRIGHT (c) 2015, Wake Forest and UMass Lowell
 * All rights reserved
 *
 * @file SF_proj.cpp
 * @brief The CPU based SF projection/backprojection in conventional method
 *
 * @version 1.0
 * @author Rui Liu
 * @date May. 1, 2015
 *
 */



#include "SF_proj.hpp"

#include <vector>
#include <algorithm>
#include <cmath>
template<typename T>
inline T cal_Tau(const T x, const T y, const T cosBeta, const T sinBeta, const T S2O, const T S2D)
{
	T tp = x * cosBeta + y * sinBeta;
	T tv = -x * sinBeta + y * cosBeta;
	T ds = S2O - tv;
	return S2D * tp / ds;
}

template<typename T>
void sort_Tau(T* tau)
{
	int i = 0;
	int j = 0;
	T temp = 0;
	for (size_t i = 0; i < 4; i++)
	{
		for (size_t j = i + 1; j < 4; j++)
		{
			if (tau[i] > tau[j])
			{
				temp = tau[i];
				tau[i] = tau[j];
				tau[j] = temp;
			}
		}
	}
}


template<typename T>
T F1(const T sk, const T rs, const T* tau)
{
	const T s1 = sk - rs / 2.0;
	const T s2 = sk + rs / 2.0;
	//Calculate gamma 1


	T b1 = std::max(s1, tau[0]);
	T b2 = std::min(s2, tau[1]);
	T gamma1 = 1.0 / (2.0 * (tau[1] - tau[0])) * (std::pow(b2 - tau[0], 2.0) - std::pow(b1 - tau[0], 2.0));
	if (b2 > b1)
	{
		//gamma1 = gamma1;
	}
	else
	{
		gamma1 = 0;
	}

	b1 = std::max(s1, tau[1]);
	b2 = std::min(s2, tau[2]);
	T gamma2 = b2 - b1;
	if (b2 > b1)
	{
		//gamma2 = gamma2;
	}
	else
	{
		gamma2 = 0;
	}

	b1 = std::max(s1, tau[2]);
	b2 = std::min(s2, tau[3]);
	T gamma3 = 1.0 / (2.0 * (tau[3] - tau[2])) * (std::pow(b1 - tau[3], 2.0) - std::pow(b2 - tau[3], 2.0));
	if (b2 > b1)
	{
		//gamma3 = gamma3;
	}
	else
	{
		gamma3 = 0;
	}

	return gamma1 + gamma2 + gamma3;
}


template<typename T>
T F2(const T tl, const T rt, const T upp_t, const T dow_t)
{
	T minr = std::min(tl + rt / 2.0f, upp_t);
	T maxr = std::max(tl - rt / 2.0f, dow_t);
	if (minr > maxr)
	{
		return (minr - maxr) / rt;
	}
	else
	{
		return 0.0;
	}
}


template<typename T>
T F2_TT(const T tl, const T rt, const T* xi)
{
	const T s1 = tl - rt / 2.0;
	const T s2 = tl + rt / 2.0;


	T b1 = std::max(s1, xi[0]);
	T b2 = std::min(s2, xi[1]);
	T gamma1 = 1.0 / (2.0 * (xi[1] - xi[0])) * (pow(b2 - xi[0], 2.0) - pow(b1 - xi[0], 2.0));
	if (b2 > b1)
	{
		//gamma1 = gamma1;
	}
	else
	{
		gamma1 = 0;
	}

	b1 = std::max(s1, xi[1]);
	b2 = std::min(s2, xi[2]);
	T gamma2 = b2 - b1;
	if (b2 > b1)
	{
		//gamma2 = gamma2;
	}
	else
	{
		gamma2 = 0;
	}

	b1 = std::max(s1, xi[2]);
	b2 = std::min(s2, xi[3]);
	T gamma3 = 1.0 / (2.0 * (xi[3] - xi[2])) * (pow(b1 - xi[3], 2.0) - pow(b2 - xi[3], 2.0));
	if (b2 > b1)
	{
		//gamma3 = gamma3;
	}
	else
	{
		gamma3 = 0;
	}

	return gamma1 + gamma2 + gamma3;


}

//This is the implementation in A1
template<typename T>
void SFTRA1_3D_ED_proj_template(std::vector<T>& proj, const std::vector<T>& vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detSizeS, const T detSizeT,
	const T detCntIdS, const T detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<T>& angs, const int angIdx)
{
	const T beta = angs[angIdx];
	const T cosBeta = std::cos(beta);
	const T sinBeta = std::sin(beta);
	const T S2D = S2O + O2D;

	const T dx = objSizeX / XN;
	const T dy = objSizeY / YN;
	const T dz = objSizeZ / ZN;
	const T ds = detSizeS / DNS;
	const T dt = detSizeT / DNT;


	const T objCntIdxX = (XN - 1.0) / 2.0;
	const T objCntIdxY = (YN - 1.0) / 2.0;
	const T objCntIdxZ = (ZN - 1.0) / 2.0;

	T tau[4];
	T curobjX, curobjY, curobjZ;
	T sk = 0;
	T Store1[20];
	T upp_t, dow_t;
	T tl = 0;
	T f2val = 0;
	int eleN = 0;
	int minSIdx, maxSIdx, minTIdx, maxTIdx;
	for (size_t ii = 0; ii < XN; ii++)
	{
		for (size_t jj = 0; jj < YN; jj++)
		{
			curobjX = (ii - objCntIdxX) * dx;
			curobjY = (jj - objCntIdxY) * dy;
			tau[0] = cal_Tau<T>(curobjX - dx / 2.0, curobjY - dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[1] = cal_Tau<T>(curobjX - dx / 2.0, curobjY + dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[2] = cal_Tau<T>(curobjX + dx / 2.0, curobjY - dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[3] = cal_Tau<T>(curobjX + dx / 2.0, curobjY + dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			sort_Tau(tau);
			minSIdx = std::floor(tau[0] / ds + detCntIdS) - 1;
			maxSIdx = std::ceil(tau[3] / ds + detCntIdS) + 1;

			if (maxSIdx < 0 || minSIdx > DNS - 1)
			{
				continue;
			}

			if (minSIdx < 0)
			{
				minSIdx = 0;
			}
			if (maxSIdx > DNS - 1)
			{
				maxSIdx = DNS - 1;
			}
			eleN = 0;

			for (int sIdx = minSIdx; sIdx <= maxSIdx; sIdx++)
			{
				sk = (sIdx - detCntIdS) * ds;
				Store1[eleN++] = F1(sk, ds, tau);
			}


			//If use SF-TR-A2, calculate the l_phi_0
			for (size_t kk = 0; kk < ZN; kk++)
			{
				//计算上下两个点对应的坐标 t;
				dow_t = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-curobjX * sinBeta + curobjY * cosBeta));
				upp_t = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-curobjX * sinBeta + curobjY * cosBeta));
				//计算 T Idx
				minTIdx = std::floor(dow_t / dt + detCntIdT) - 1;
				maxTIdx = std::ceil(upp_t / dt + detCntIdT) + 1;
				if (maxTIdx < 0 || minTIdx > DNT - 1)
				{
					continue;
				}
				if (minTIdx < 0)
				{
					minTIdx = 0;
				}
				if (maxTIdx > DNT - 1)
				{
					maxTIdx = DNT - 1;
				}
				//std::cout << (maxTIdx - minTIdx) << " ";

				for (int tIdx = minTIdx; tIdx <= maxTIdx; tIdx++)
				{
					tl = (tIdx - detCntIdT) * dt;
					f2val = F2(tl, dt, upp_t, dow_t);
					for (int sIdx = minSIdx; sIdx <= maxSIdx; sIdx++)
					{
						proj[(angIdx * DNT + tIdx) * DNS + sIdx] +=
							vol[(kk * YN + jj) * XN + ii] * Store1[sIdx - minSIdx] * f2val;
					}
				}

			}
		}
	}

	T cur_s, cur_t;
	T curGamma;
	T curPhi;
	T curTheta;
	T l_theta, l_phi;
	// Scale all the projection using (36)
	for (int sIdx = 0; sIdx < DNS; sIdx++)
	{
		cur_s = (sIdx - detCntIdS) * ds;
		curGamma = atan(cur_s / S2D);
		curPhi = curGamma + beta;
		l_phi = dx / std::max(abs(cos(curPhi)), abs(sin(curPhi)));

		for (int tIdx = 0; tIdx < DNT; tIdx++)
		{
			cur_t = (tIdx - detCntIdT) * dt;
			curTheta = -atan(cur_t / hypot(cur_s, S2D));
			l_theta = 1.0 / abs(cos(curTheta));
			proj[(angIdx * DNT + tIdx) * DNS + sIdx] *= (l_phi * l_theta);

		}
	}
}

void SFTRA1_3D_ED_proj(std::vector<float>& proj, const std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx)
{
	SFTRA1_3D_ED_proj_template<float>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ, detSizeS, detSizeT,
		detCntIdS, detCntIdT, XN, YN, ZN, DNS, DNT, PN, angs, angIdx);
}


void SFTRA1_3D_ED_proj(std::vector<double>& proj, const std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx)
{
	SFTRA1_3D_ED_proj_template<double>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ, detSizeS, detSizeT,
		detCntIdS, detCntIdT, XN, YN, ZN, DNS, DNT, PN, angs, angIdx);
}





template<typename T>
void SFTRA1_3D_ED_bproj_template(const std::vector<T>& proj, std::vector<T>& vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detSizeS, const T detSizeT,
	const T detCntIdS, const T detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<T>& angs, const int angIdx)
{
	const T beta = angs[angIdx];
	const T cosBeta = cos(beta);
	const T sinBeta = sin(beta);
	const T S2D = S2O + O2D;

	const T dx = objSizeX / XN;
	const T dy = objSizeY / YN;
	const T dz = objSizeZ / ZN;
	const T ds = detSizeS / DNS;
	const T dt = detSizeT / DNT;

	const T objCntIdxX = (XN - 1.0) / 2.0;
	const T objCntIdxY = (YN - 1.0) / 2.0;
	const T objCntIdxZ = (ZN - 1.0) / 2.0;

	T tau[4];
	T curobjX, curobjY, curobjZ;
	T sk = 0;
	T Store1[20];
	T upp_t, dow_t;
	T tl = 0;
	T f2val = 0;
	int eleN = 0;
	int minSIdx, maxSIdx, minTIdx, maxTIdx;

	T cur_s, cur_t;
	T curGamma;
	T curPhi;
	T curTheta;
	T l_theta, l_phi;
	T summ = 0;
	for (size_t ii = 0; ii < XN; ii++)
	{
		for (size_t jj = 0; jj < YN; jj++)
		{
			curobjX = (ii - objCntIdxX) * dx;
			curobjY = (jj - objCntIdxY) * dy;
			tau[0] = cal_Tau<T>(curobjX - dx / 2.0, curobjY - dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[1] = cal_Tau<T>(curobjX - dx / 2.0, curobjY + dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[2] = cal_Tau<T>(curobjX + dx / 2.0, curobjY - dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[3] = cal_Tau<T>(curobjX + dx / 2.0, curobjY + dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			sort_Tau(tau);
			minSIdx = floor(tau[0] / ds + detCntIdS) - 1;
			maxSIdx = ceil(tau[3] / ds + detCntIdS) + 1;

			if (maxSIdx < 0 || minSIdx > DNS - 1)
			{
				continue;
			}

			if (minSIdx < 0)
			{
				minSIdx = 0;
			}
			if (maxSIdx > DNS - 1)
			{
				maxSIdx = DNS - 1;
			}
			eleN = 0;

			for (int sIdx = minSIdx; sIdx <= maxSIdx; sIdx++)
			{
				sk = (sIdx - detCntIdS) * ds;
				Store1[eleN++] = F1(sk, ds, tau);
			}


			//If use SF-TR-A2, calculate the l_phi_0
			for (size_t kk = 0; kk < ZN; kk++)
			{
				//计算上下两个点对应的坐标 t;
				dow_t = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-curobjX * sinBeta + curobjY * cosBeta));
				upp_t = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-curobjX * sinBeta + curobjY * cosBeta));
				//计算 T Idx
				minTIdx = floor(dow_t / dt + detCntIdT) - 1;
				maxTIdx = ceil(upp_t / dt + detCntIdT) + 1;
				if (maxTIdx < 0 || minTIdx > DNT - 1)
				{
					continue;
				}
				if (minTIdx < 0)
				{
					minTIdx = 0;
				}
				if (maxTIdx > DNT - 1)
				{
					maxTIdx = DNT - 1;
				}
				//std::cout << (maxTIdx - minTIdx) << " ";
				summ = 0;
				for (int tIdx = minTIdx; tIdx <= maxTIdx; tIdx++)
				{
					tl = (tIdx - detCntIdT) * dt;
					f2val = F2(tl, dt, upp_t, dow_t);

					cur_t = (tIdx - detCntIdT) * dt;
					curTheta = -atan(cur_t / hypot(cur_s, S2D));
					l_theta = 1.0 / abs(cos(curTheta));
					for (int sIdx = minSIdx; sIdx <= maxSIdx; sIdx++)
					{
						cur_s = (sIdx - detCntIdS) * ds;
						curGamma = atan(cur_s / S2D);
						curPhi = curGamma + beta;
						l_phi = dx / std::max(abs(cos(curPhi)), abs(sin(curPhi)));

						summ += proj[(angIdx * DNT + tIdx) * DNS + sIdx] * (l_phi * l_theta)  * Store1[sIdx - minSIdx] * f2val;

					}
				}
				vol[(kk * YN + jj) * XN + ii] = summ;
			}
		}
	}

}

void SFTRA1_3D_ED_bproj(const std::vector<float>& proj, std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx)
{
	SFTRA1_3D_ED_bproj_template<float>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ, detSizeS, detSizeT,
		detCntIdS, detCntIdT, XN, YN, ZN, DNS, DNT, PN, angs, angIdx);
}


void SFTRA1_3D_ED_bproj(const std::vector<double>& proj, std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx)
{
	SFTRA1_3D_ED_bproj_template<double>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ, detSizeS, detSizeT,
		detCntIdS, detCntIdT, XN, YN, ZN, DNS, DNT, PN, angs, angIdx);
}



template<typename T>
void SFTRA2_3D_ED_proj_template(std::vector<T>& proj, const std::vector<T>& vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detSizeS, const T detSizeT,
	const T detCntIdS, const T detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<T>& angs, const int angIdx)
{
	const T beta = angs[angIdx];
	const T cosBeta = std::cos(beta);
	const T sinBeta = std::sin(beta);
	const T S2D = S2O + O2D;

	const T dx = objSizeX / XN;
	const T dy = objSizeY / YN;
	const T dz = objSizeZ / ZN;
	const T ds = detSizeS / DNS;
	const T dt = detSizeT / DNT;


	const T objCntIdxX = (XN - 1.0) / 2.0;
	const T objCntIdxY = (YN - 1.0) / 2.0;
	const T objCntIdxZ = (ZN - 1.0) / 2.0;

	T tau[4];
	T curobjX, curobjY, curobjZ;
	T sk = 0;
	T Store1[20];
	T upp_t, dow_t;
	T tl = 0;
	T f2val = 0;
	int eleN = 0;
	T l_phi0 = 0;
	T curphi = 0;
	int minSIdx, maxSIdx, minTIdx, maxTIdx;
	for (size_t ii = 0; ii < XN; ii++)
	{
		for (size_t jj = 0; jj < YN; jj++)
		{
			curobjX = (ii - objCntIdxX) * dx;
			curobjY = (jj - objCntIdxY) * dy;
			tau[0] = cal_Tau<T>(curobjX - dx / 2.0, curobjY - dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[1] = cal_Tau<T>(curobjX - dx / 2.0, curobjY + dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[2] = cal_Tau<T>(curobjX + dx / 2.0, curobjY - dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[3] = cal_Tau<T>(curobjX + dx / 2.0, curobjY + dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			sort_Tau(tau);
			minSIdx = std::floor(tau[0] / ds + detCntIdS) - 1;
			maxSIdx = std::ceil(tau[3] / ds + detCntIdS) + 1;


			if (maxSIdx < 0 || minSIdx > DNS - 1)
			{
				continue;
			}

			if (minSIdx < 0)
			{
				minSIdx = 0;
			}
			if (maxSIdx > DNS - 1)
			{
				maxSIdx = DNS - 1;
			}
			eleN = 0;

			for (int sIdx = minSIdx; sIdx <= maxSIdx; sIdx++)
			{
				sk = (sIdx - detCntIdS) * ds;
				Store1[eleN++] = F1(sk, ds, tau);
			}

			curphi = beta + atan((curobjX * cosBeta + curobjY * sinBeta) / (S2O - (-curobjX * sinBeta + curobjY * cosBeta)));
			l_phi0 = dx / std::max(abs(cos(curphi)), abs(sin(curphi)));
			//If use SF-TR-A2, calculate the l_phi_0
			for (size_t kk = 0; kk < ZN; kk++)
			{
				//计算上下两个点对应的坐标 t;
				dow_t = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-curobjX * sinBeta + curobjY * cosBeta));
				upp_t = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-curobjX * sinBeta + curobjY * cosBeta));
				//计算 T Idx
				minTIdx = floor(dow_t / dt + detCntIdT) - 1;
				maxTIdx = ceil(upp_t / dt + detCntIdT) + 1;
				if (maxTIdx < 0 || minTIdx > DNT - 1)
				{
					continue;
				}
				if (minTIdx < 0)
				{
					minTIdx = 0;
				}
				if (maxTIdx > DNT - 1)
				{
					maxTIdx = DNT - 1;
				}
				//std::cout << (maxTIdx - minTIdx) << " ";

				for (int tIdx = minTIdx; tIdx <= maxTIdx; tIdx++)
				{
					tl = (tIdx - detCntIdT) * dt;
					f2val = F2(tl, dt, upp_t, dow_t);
					for (int sIdx = minSIdx; sIdx <= maxSIdx; sIdx++)
					{
						proj[(angIdx * DNT + tIdx) * DNS + sIdx] +=
							vol[(kk * YN + jj) * XN + ii] * Store1[sIdx - minSIdx] * f2val * l_phi0;
					}
				}

			}
		}
	}

	T cur_s, cur_t;

	T curTheta;
	T l_theta;
	// Scale all the projection using (36)
	for (int sIdx = 0; sIdx < DNS; sIdx++)
	{
		cur_s = (sIdx - detCntIdS) * ds;
		for (int tIdx = 0; tIdx < DNT; tIdx++)
		{
			cur_t = (tIdx - detCntIdT) * dt;
			curTheta = -std::atan(cur_t / hypot(cur_s, S2D));
			l_theta = 1.0 / abs(std::cos(curTheta));
			proj[(angIdx * DNT + tIdx) * DNS + sIdx] *= (l_theta);

		}
	}
}


void SFTRA2_3D_ED_proj(std::vector<float>& proj, const std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx)
{
	SFTRA2_3D_ED_proj_template<float>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ, detSizeS, detSizeT,
		detCntIdS, detCntIdT, XN, YN, ZN, DNS, DNT, PN, angs, angIdx);
}

void SFTRA2_3D_ED_proj(std::vector<double>& proj, const std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx)
{
	SFTRA2_3D_ED_proj_template<double>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ, detSizeS, detSizeT,
		detCntIdS, detCntIdT, XN, YN, ZN, DNS, DNT, PN, angs, angIdx);
}


template<typename T>
void SFTRA2_3D_ED_bproj_template(const std::vector<T>& proj, std::vector<T>& vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detSizeS, const T detSizeT,
	const T detCntIdS, const T detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<T>& angs, const int angIdx)
{
	const T beta = angs[angIdx];
	const T cosBeta = std::cos(beta);
	const T sinBeta = std::sin(beta);
	const T S2D = S2O + O2D;

	const T dx = objSizeX / XN;
	const T dy = objSizeY / YN;
	const T dz = objSizeZ / ZN;
	const T ds = detSizeS / DNS;
	const T dt = detSizeT / DNT;


	const T objCntIdxX = (XN - 1.0) / 2.0;
	const T objCntIdxY = (YN - 1.0) / 2.0;
	const T objCntIdxZ = (ZN - 1.0) / 2.0;

	T tau[4];
	T curobjX, curobjY, curobjZ;
	T sk = 0;
	T Store1[20];
	T upp_t, dow_t;
	T tl = 0;
	T f2val = 0;
	int eleN = 0;
	T l_phi0 = 0;
	T curphi = 0;
	int minSIdx, maxSIdx, minTIdx, maxTIdx;

	T cur_s, cur_t;

	T curTheta;
	T l_theta;
	T summ = 0;
	for (size_t ii = 0; ii < XN; ii++)
	{
		for (size_t jj = 0; jj < YN; jj++)
		{
			curobjX = (ii - objCntIdxX) * dx;
			curobjY = (jj - objCntIdxY) * dy;
			tau[0] = cal_Tau<T>(curobjX - dx / 2.0, curobjY - dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[1] = cal_Tau<T>(curobjX - dx / 2.0, curobjY + dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[2] = cal_Tau<T>(curobjX + dx / 2.0, curobjY - dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[3] = cal_Tau<T>(curobjX + dx / 2.0, curobjY + dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			sort_Tau(tau);
			minSIdx = floor(tau[0] / ds + detCntIdS) - 1;
			maxSIdx = ceil(tau[3] / ds + detCntIdS) + 1;

			if (maxSIdx < 0 || minSIdx > DNS - 1)
			{
				continue;
			}

			if (minSIdx < 0)
			{
				minSIdx = 0;
			}
			if (maxSIdx > DNS - 1)
			{
				maxSIdx = DNS - 1;
			}
			eleN = 0;

			for (int sIdx = minSIdx; sIdx <= maxSIdx; sIdx++)
			{
				sk = (sIdx - detCntIdS) * ds;
				Store1[eleN++] = F1(sk, ds, tau);
			}

			curphi = beta + std::atan((curobjX * cosBeta + curobjY * sinBeta) / (S2O - (-curobjX * sinBeta + curobjY * cosBeta)));
			l_phi0 = dx / std::max(abs(std::cos(curphi)), abs(std::sin(curphi)));
			//If use SF-TR-A2, calculate the l_phi_0
			for (size_t kk = 0; kk < ZN; kk++)
			{
				//计算上下两个点对应的坐标 t;
				dow_t = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-curobjX * sinBeta + curobjY * cosBeta));
				upp_t = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-curobjX * sinBeta + curobjY * cosBeta));
				//计算 T Idx
				minTIdx = floor(dow_t / dt + detCntIdT) - 1;
				maxTIdx = ceil(upp_t / dt + detCntIdT) + 1;
				if (maxTIdx < 0 || minTIdx > DNT - 1)
				{
					continue;
				}
				if (minTIdx < 0)
				{
					minTIdx = 0;
				}
				if (maxTIdx > DNT - 1)
				{
					maxTIdx = DNT - 1;
				}
				//std::cout << (maxTIdx - minTIdx) << " ";
				summ = 0;
				for (int tIdx = minTIdx; tIdx <= maxTIdx; tIdx++)
				{
					tl = (tIdx - detCntIdT) * dt;
					f2val = F2(tl, dt, upp_t, dow_t);
					for (int sIdx = minSIdx; sIdx <= maxSIdx; sIdx++)
					{
						cur_s = (sIdx - detCntIdS) * ds;
						cur_t = (tIdx - detCntIdT) * dt;

						curTheta = -std::atan(cur_t / hypot(cur_s, S2D));
						l_theta = 1.0 / abs(cos(curTheta));
						summ +=
							proj[(angIdx * DNT + tIdx) * DNS + sIdx] * l_theta
							* Store1[sIdx - minSIdx] * f2val * l_phi0;
					}
				}
				vol[(kk * YN + jj) * XN + ii] = summ;

			}
		}
	}

}



void SFTRA2_3D_ED_bproj(const std::vector<float>& proj, std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx)
{
	SFTRA2_3D_ED_bproj_template<float>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ, detSizeS, detSizeT,
		detCntIdS, detCntIdT, XN, YN, ZN, DNS, DNT, PN, angs, angIdx);
}


void SFTRA2_3D_ED_bproj(const std::vector<double>& proj, std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx)
{
	SFTRA2_3D_ED_bproj_template<double>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ, detSizeS, detSizeT,
		detCntIdS, detCntIdT, XN, YN, ZN, DNS, DNT, PN, angs, angIdx);
}


template<typename T>
void SFTTA1_3D_ED_proj_template(std::vector<T>& proj, const std::vector<T>& vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detSizeS, const T detSizeT,
	const T detCntIdS, const T detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<T>& angs, const int angIdx)
{
	const T beta = angs[angIdx];
	const T cosBeta = cos(beta);
	const T sinBeta = sin(beta);
	const T S2D = S2O + O2D;

	const T dx = objSizeX / XN;
	const T dy = objSizeY / YN;
	const T dz = objSizeZ / ZN;
	const T ds = detSizeS / DNS;
	const T dt = detSizeT / DNT;


	const T objCntIdxX = (XN - 1.0) / 2.0;
	const T objCntIdxY = (YN - 1.0) / 2.0;
	const T objCntIdxZ = (ZN - 1.0) / 2.0;

	T tau[4];
	T xi_dow[4], xi_upp[4], xi[4];
	T curobjX, curobjY, curobjZ;
	T sk = 0;
	T Store1[20];

	T upp_t, dow_t;
	T tl = 0;
	T f2val = 0;
	int eleN = 0;
	int minSIdx, maxSIdx, minTIdx, maxTIdx;
	for (size_t ii = 0; ii < XN; ii++)
	{
		for (size_t jj = 0; jj < YN; jj++)
		{
			curobjX = (ii - objCntIdxX) * dx;
			curobjY = (jj - objCntIdxY) * dy;
			tau[0] = cal_Tau<T>(curobjX - dx / 2.0, curobjY - dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[1] = cal_Tau<T>(curobjX - dx / 2.0, curobjY + dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[2] = cal_Tau<T>(curobjX + dx / 2.0, curobjY - dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[3] = cal_Tau<T>(curobjX + dx / 2.0, curobjY + dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			sort_Tau(tau);
			minSIdx = floor(tau[0] / ds + detCntIdS) - 1;
			maxSIdx = ceil(tau[3] / ds + detCntIdS) + 1;

			if (maxSIdx < 0 || minSIdx > DNS - 1)
			{
				continue;
			}

			if (minSIdx < 0)
			{
				minSIdx = 0;
			}
			if (maxSIdx > DNS - 1)
			{
				maxSIdx = DNS - 1;
			}
			eleN = 0;

			for (int sIdx = minSIdx; sIdx <= maxSIdx; sIdx++)
			{
				sk = (sIdx - detCntIdS) * ds;
				Store1[eleN++] = F1(sk, ds, tau);
			}


			//If use SF-TR-A2, calculate the l_phi_0
			for (size_t kk = 0; kk < ZN; kk++)
			{
				//计算上下两个点对应的坐标 t;
				dow_t = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-curobjX * sinBeta + curobjY * cosBeta));
				upp_t = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-curobjX * sinBeta + curobjY * cosBeta));
				//计算 T Idx
				minTIdx = floor(dow_t / dt + detCntIdT) - 1;
				maxTIdx = ceil(upp_t / dt + detCntIdT) + 1;
				if (maxTIdx < 0 || minTIdx > DNT - 1)
				{
					continue;
				}
				if (minTIdx < 0)
				{
					minTIdx = 0;
				}
				if (maxTIdx > DNT - 1)
				{
					maxTIdx = DNT - 1;
				}

				xi_dow[0] = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-(curobjX - 0.5 * dx) * sinBeta + (curobjY - 0.5 * dy) * cosBeta));
				xi_dow[1] = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-(curobjX - 0.5 * dx) * sinBeta + (curobjY + 0.5 * dy) * cosBeta));
				xi_dow[2] = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-(curobjX + 0.5 * dx) * sinBeta + (curobjY - 0.5 * dy) * cosBeta));
				xi_dow[3] = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-(curobjX + 0.5 * dx) * sinBeta + (curobjY + 0.5 * dy) * cosBeta));
				sort_Tau(xi_dow);
				xi_upp[0] = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-(curobjX - 0.5 * dx) * sinBeta + (curobjY - 0.5 * dy) * cosBeta));
				xi_upp[1] = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-(curobjX - 0.5 * dx) * sinBeta + (curobjY + 0.5 * dy) * cosBeta));
				xi_upp[2] = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-(curobjX + 0.5 * dx) * sinBeta + (curobjY - 0.5 * dy) * cosBeta));
				xi_upp[3] = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-(curobjX + 0.5 * dx) * sinBeta + (curobjY + 0.5 * dy) * cosBeta));
				sort_Tau(xi_upp);
				xi[0] = xi_dow[0];
				xi[1] = xi_dow[3];
				xi[2] = xi_upp[0];
				xi[3] = xi_upp[3];

				for (int tIdx = minTIdx; tIdx <= maxTIdx; tIdx++)
				{
					tl = (tIdx - detCntIdT) * dt;
					f2val = F2_TT(tl, dt, xi);
					for (int sIdx = minSIdx; sIdx <= maxSIdx; sIdx++)
					{
						proj[(angIdx * DNT + tIdx) * DNS + sIdx] +=
							vol[(kk * YN + jj) * XN + ii] * Store1[sIdx - minSIdx] * f2val;
					}
				}

			}
		}
	}

	T cur_s, cur_t;
	T curGamma;
	T curPhi;
	T curTheta;
	T l_theta, l_phi;
	// Scale all the projection using (36)
	for (int sIdx = 0; sIdx < DNS; sIdx++)
	{
		cur_s = (sIdx - detCntIdS) * ds;
		curGamma = atan(cur_s / S2D);
		curPhi = curGamma + beta;
		l_phi = dx / std::max(abs(cos(curPhi)), abs(sin(curPhi)));

		for (int tIdx = 0; tIdx < DNT; tIdx++)
		{
			cur_t = (tIdx - detCntIdT) * dt;
			curTheta = -atan(cur_t / hypot(cur_s, S2D));
			l_theta = 1.0 / abs(cos(curTheta));
			proj[(angIdx * DNT + tIdx) * DNS + sIdx] *= (l_phi * l_theta);

		}
	}
}


void SFTTA1_3D_ED_proj(std::vector<float>& proj, const std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx)
{
	SFTTA1_3D_ED_proj_template<float>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ, detSizeS, detSizeT,
		detCntIdS, detCntIdT, XN, YN, ZN, DNS, DNT, PN, angs, angIdx);
}

void SFTTA1_3D_ED_proj(std::vector<double>& proj, const std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx)
{
	SFTTA1_3D_ED_proj_template<double>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ, detSizeS, detSizeT,
		detCntIdS, detCntIdT, XN, YN, ZN, DNS, DNT, PN, angs, angIdx);
}



template<typename T>
void SFTTA1_3D_ED_bproj_template(const std::vector<T>& proj, std::vector<T>& vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detSizeS, const T detSizeT,
	const T detCntIdS, const T detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<T>& angs, const int angIdx)
{
	const T beta = angs[angIdx];
	const T cosBeta = cos(beta);
	const T sinBeta = sin(beta);
	const T S2D = S2O + O2D;

	const T dx = objSizeX / XN;
	const T dy = objSizeY / YN;
	const T dz = objSizeZ / ZN;
	const T ds = detSizeS / DNS;
	const T dt = detSizeT / DNT;


	const T objCntIdxX = (XN - 1.0) / 2.0;
	const T objCntIdxY = (YN - 1.0) / 2.0;
	const T objCntIdxZ = (ZN - 1.0) / 2.0;

	T tau[4];
	T xi_dow[4], xi_upp[4], xi[4];
	T curobjX, curobjY, curobjZ;
	T sk = 0;
	T Store1[20];

	T upp_t, dow_t;
	T tl = 0;
	T f2val = 0;
	int eleN = 0;
	int minSIdx, maxSIdx, minTIdx, maxTIdx;

	T cur_s, cur_t;
	T curGamma;
	T curPhi;
	T curTheta;
	T l_theta, l_phi;
	T summ = 0;

	for (size_t ii = 0; ii < XN; ii++)
	{
		for (size_t jj = 0; jj < YN; jj++)
		{
			curobjX = (ii - objCntIdxX) * dx;
			curobjY = (jj - objCntIdxY) * dy;
			tau[0] = cal_Tau<T>(curobjX - dx / 2.0, curobjY - dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[1] = cal_Tau<T>(curobjX - dx / 2.0, curobjY + dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[2] = cal_Tau<T>(curobjX + dx / 2.0, curobjY - dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[3] = cal_Tau<T>(curobjX + dx / 2.0, curobjY + dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			sort_Tau(tau);
			minSIdx = floor(tau[0] / ds + detCntIdS) - 1;
			maxSIdx = ceil(tau[3] / ds + detCntIdS) + 1;

			if (maxSIdx < 0 || minSIdx > DNS - 1)
			{
				continue;
			}

			if (minSIdx < 0)
			{
				minSIdx = 0;
			}
			if (maxSIdx > DNS - 1)
			{
				maxSIdx = DNS - 1;
			}
			eleN = 0;

			for (int sIdx = minSIdx; sIdx <= maxSIdx; sIdx++)
			{
				sk = (sIdx - detCntIdS) * ds;
				Store1[eleN++] = F1(sk, ds, tau);
			}


			//If use SF-TR-A2, calculate the l_phi_0
			for (size_t kk = 0; kk < ZN; kk++)
			{
				//计算上下两个点对应的坐标 t;
				dow_t = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-curobjX * sinBeta + curobjY * cosBeta));
				upp_t = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-curobjX * sinBeta + curobjY * cosBeta));
				//计算 T Idx
				minTIdx = floor(dow_t / dt + detCntIdT) - 1;
				maxTIdx = ceil(upp_t / dt + detCntIdT) + 1;
				if (maxTIdx < 0 || minTIdx > DNT - 1)
				{
					continue;
				}
				if (minTIdx < 0)
				{
					minTIdx = 0;
				}
				if (maxTIdx > DNT - 1)
				{
					maxTIdx = DNT - 1;
				}

				xi_dow[0] = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-(curobjX - 0.5 * dx) * sinBeta + (curobjY - 0.5 * dy) * cosBeta));
				xi_dow[1] = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-(curobjX - 0.5 * dx) * sinBeta + (curobjY + 0.5 * dy) * cosBeta));
				xi_dow[2] = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-(curobjX + 0.5 * dx) * sinBeta + (curobjY - 0.5 * dy) * cosBeta));
				xi_dow[3] = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-(curobjX + 0.5 * dx) * sinBeta + (curobjY + 0.5 * dy) * cosBeta));
				sort_Tau(xi_dow);
				xi_upp[0] = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-(curobjX - 0.5 * dx) * sinBeta + (curobjY - 0.5 * dy) * cosBeta));
				xi_upp[1] = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-(curobjX - 0.5 * dx) * sinBeta + (curobjY + 0.5 * dy) * cosBeta));
				xi_upp[2] = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-(curobjX + 0.5 * dx) * sinBeta + (curobjY - 0.5 * dy) * cosBeta));
				xi_upp[3] = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-(curobjX + 0.5 * dx) * sinBeta + (curobjY + 0.5 * dy) * cosBeta));
				sort_Tau(xi_upp);
				xi[0] = xi_dow[0];
				xi[1] = xi_dow[3];
				xi[2] = xi_upp[0];
				xi[3] = xi_upp[3];

				summ = 0;
				for (int tIdx = minTIdx; tIdx <= maxTIdx; tIdx++)
				{
					tl = (tIdx - detCntIdT) * dt;
					f2val = F2_TT(tl, dt, xi);
					cur_t = (tIdx - detCntIdT) * dt;
					curTheta = -atan(cur_t / hypot(cur_s, S2D));
					l_theta = 1.0 / abs(cos(curTheta));

					for (int sIdx = minSIdx; sIdx <= maxSIdx; sIdx++)
					{
						cur_s = (sIdx - detCntIdS) * ds;
						curGamma = atan(cur_s / S2D);
						curPhi = curGamma + beta;
						l_phi = dx / std::max(abs(cos(curPhi)), abs(sin(curPhi)));
						summ += proj[(angIdx * DNT + tIdx) * DNS + sIdx] * (l_phi * l_theta) * Store1[sIdx - minSIdx] * f2val;
					}
				}
				vol[(kk * YN + jj) * XN + ii] = summ;
			}
		}
	}
}



void SFTTA1_3D_ED_bproj(const std::vector<float>& proj, std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx)
{
	SFTTA1_3D_ED_bproj_template<float>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ, detSizeS, detSizeT,
		detCntIdS, detCntIdT, XN, YN, ZN, DNS, DNT, PN, angs, angIdx);
}


void SFTTA1_3D_ED_bproj(const std::vector<double>& proj, std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx)
{
	SFTTA1_3D_ED_bproj_template<double>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ, detSizeS, detSizeT,
		detCntIdS, detCntIdT, XN, YN, ZN, DNS, DNT, PN, angs, angIdx);
}


template<typename T>
void SFTTA2_3D_ED_proj_template(std::vector<T>& proj, const std::vector<T>& vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detSizeS, const T detSizeT,
	const T detCntIdS, const T detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<T>& angs, const int angIdx)
{
	const T beta = angs[angIdx];
	const T cosBeta = cos(beta);
	const T sinBeta = sin(beta);
	const T S2D = S2O + O2D;

	const T dx = objSizeX / XN;
	const T dy = objSizeY / YN;
	const T dz = objSizeZ / ZN;
	const T ds = detSizeS / DNS;
	const T dt = detSizeT / DNT;


	const T objCntIdxX = (XN - 1.0) / 2.0;
	const T objCntIdxY = (YN - 1.0) / 2.0;
	const T objCntIdxZ = (ZN - 1.0) / 2.0;

	T tau[4];
	T curobjX, curobjY, curobjZ;
	T sk = 0;
	T Store1[20];
	T upp_t, dow_t;
	T tl = 0;
	T f2val = 0;
	int eleN = 0;
	T l_phi0 = 0;
	T curphi = 0;
	T xi_dow[4], xi_upp[4], xi[4];
	int minSIdx, maxSIdx, minTIdx, maxTIdx;
	for (size_t ii = 0; ii < XN; ii++)
	{
		for (size_t jj = 0; jj < YN; jj++)
		{
			curobjX = (ii - objCntIdxX) * dx;
			curobjY = (jj - objCntIdxY) * dy;
			tau[0] = cal_Tau<T>(curobjX - dx / 2.0, curobjY - dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[1] = cal_Tau<T>(curobjX - dx / 2.0, curobjY + dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[2] = cal_Tau<T>(curobjX + dx / 2.0, curobjY - dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[3] = cal_Tau<T>(curobjX + dx / 2.0, curobjY + dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			sort_Tau(tau);
			minSIdx = floor(tau[0] / ds + detCntIdS) - 1;
			maxSIdx = ceil(tau[3] / ds + detCntIdS) + 1;

			if (maxSIdx < 0 || minSIdx > DNS - 1)
			{
				continue;
			}

			if (minSIdx < 0)
			{
				minSIdx = 0;
			}
			if (maxSIdx > DNS - 1)
			{
				maxSIdx = DNS - 1;
			}
			eleN = 0;

			for (int sIdx = minSIdx; sIdx <= maxSIdx; sIdx++)
			{
				sk = (sIdx - detCntIdS) * ds;
				Store1[eleN++] = F1(sk, ds, tau);
			}

			curphi = beta + atan((curobjX * cosBeta + curobjY * sinBeta) / (S2O - (-curobjX * sinBeta + curobjY * cosBeta)));
			l_phi0 = dx / std::max(abs(cos(curphi)), abs(sin(curphi)));
			//If use SF-TR-A2, calculate the l_phi_0
			for (size_t kk = 0; kk < ZN; kk++)
			{
				//计算上下两个点对应的坐标 t;
				dow_t = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-curobjX * sinBeta + curobjY * cosBeta));
				upp_t = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-curobjX * sinBeta + curobjY * cosBeta));
				//计算 T Idx
				minTIdx = floor(dow_t / dt + detCntIdT) - 1;
				maxTIdx = ceil(upp_t / dt + detCntIdT) + 1;
				if (maxTIdx < 0 || minTIdx > DNT - 1)
				{
					continue;
				}
				if (minTIdx < 0)
				{
					minTIdx = 0;
				}
				if (maxTIdx > DNT - 1)
				{
					maxTIdx = DNT - 1;
				}
				//std::cout << (maxTIdx - minTIdx) << " ";

				xi_dow[0] = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-(curobjX - 0.5 * dx) * sinBeta + (curobjY - 0.5 * dy) * cosBeta));
				xi_dow[1] = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-(curobjX - 0.5 * dx) * sinBeta + (curobjY + 0.5 * dy) * cosBeta));
				xi_dow[2] = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-(curobjX + 0.5 * dx) * sinBeta + (curobjY - 0.5 * dy) * cosBeta));
				xi_dow[3] = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-(curobjX + 0.5 * dx) * sinBeta + (curobjY + 0.5 * dy) * cosBeta));
				sort_Tau(xi_dow);
				xi_upp[0] = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-(curobjX - 0.5 * dx) * sinBeta + (curobjY - 0.5 * dy) * cosBeta));
				xi_upp[1] = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-(curobjX - 0.5 * dx) * sinBeta + (curobjY + 0.5 * dy) * cosBeta));
				xi_upp[2] = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-(curobjX + 0.5 * dx) * sinBeta + (curobjY - 0.5 * dy) * cosBeta));
				xi_upp[3] = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-(curobjX + 0.5 * dx) * sinBeta + (curobjY + 0.5 * dy) * cosBeta));
				sort_Tau(xi_upp);
				xi[0] = xi_dow[0];
				xi[1] = xi_dow[3];
				xi[2] = xi_upp[0];
				xi[3] = xi_upp[3];


				for (int tIdx = minTIdx; tIdx <= maxTIdx; tIdx++)
				{
					tl = (tIdx - detCntIdT) * dt;
					f2val = F2_TT(tl, dt, xi);
					for (int sIdx = minSIdx; sIdx <= maxSIdx; sIdx++)
					{
						proj[(angIdx * DNT + tIdx) * DNS + sIdx] +=
							vol[(kk * YN + jj) * XN + ii] * Store1[sIdx - minSIdx] * f2val * l_phi0;
					}
				}

			}
		}
	}

	T cur_s, cur_t;

	T curTheta;
	T l_theta;
	// Scale all the projection using (36)
	for (int sIdx = 0; sIdx < DNS; sIdx++)
	{
		cur_s = (sIdx - detCntIdS) * ds;
		for (int tIdx = 0; tIdx < DNT; tIdx++)
		{
			cur_t = (tIdx - detCntIdT) * dt;
			curTheta = -atan(cur_t / hypot(cur_s, S2D));
			l_theta = 1.0 / abs(cos(curTheta));
			proj[(angIdx * DNT + tIdx) * DNS + sIdx] *= (l_theta);

		}
	}
}



void SFTTA2_3D_ED_proj(std::vector<float>& proj, const std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx)
{
	SFTTA2_3D_ED_proj_template<float>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ, detSizeS, detSizeT,
		detCntIdS, detCntIdT, XN, YN, ZN, DNS, DNT, PN, angs, angIdx);
}

void SFTTA2_3D_ED_proj(std::vector<double>& proj, const std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx)
{
	SFTTA2_3D_ED_proj_template<double>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ, detSizeS, detSizeT,
		detCntIdS, detCntIdT, XN, YN, ZN, DNS, DNT, PN, angs, angIdx);
}



template<typename T>
void SFTTA2_3D_ED_bproj_template(const std::vector<T>& proj, std::vector<T>& vol,
	const T S2O, const T O2D,
	const T objSizeX, const T objSizeY, const T objSizeZ,
	const T detSizeS, const T detSizeT,
	const T detCntIdS, const T detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<T>& angs, const int angIdx)
{
	const T beta = angs[angIdx];
	const T cosBeta = cos(beta);
	const T sinBeta = sin(beta);
	const T S2D = S2O + O2D;

	const T dx = objSizeX / XN;
	const T dy = objSizeY / YN;
	const T dz = objSizeZ / ZN;
	const T ds = detSizeS / DNS;
	const T dt = detSizeT / DNT;


	const T objCntIdxX = (XN - 1.0) / 2.0;
	const T objCntIdxY = (YN - 1.0) / 2.0;
	const T objCntIdxZ = (ZN - 1.0) / 2.0;

	T tau[4];
	T curobjX, curobjY, curobjZ;
	T sk = 0;
	T Store1[20];
	T upp_t, dow_t;
	T tl = 0;
	T f2val = 0;
	int eleN = 0;
	T l_phi0 = 0;
	T curphi = 0;
	T xi_dow[4], xi_upp[4], xi[4];

	T cur_s, cur_t;

	T curTheta;
	T l_theta;

	T summ = 0;
	int minSIdx, maxSIdx, minTIdx, maxTIdx;
	for (size_t ii = 0; ii < XN; ii++)
	{
		for (size_t jj = 0; jj < YN; jj++)
		{
			curobjX = (ii - objCntIdxX) * dx;
			curobjY = (jj - objCntIdxY) * dy;
			tau[0] = cal_Tau<T>(curobjX - dx / 2.0, curobjY - dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[1] = cal_Tau<T>(curobjX - dx / 2.0, curobjY + dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[2] = cal_Tau<T>(curobjX + dx / 2.0, curobjY - dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			tau[3] = cal_Tau<T>(curobjX + dx / 2.0, curobjY + dy / 2.0, cosBeta, sinBeta, S2O, S2D);
			sort_Tau(tau);
			minSIdx = std::floor(tau[0] / ds + detCntIdS) - 1;
			maxSIdx = std::ceil(tau[3] / ds + detCntIdS) + 1;

			if (maxSIdx < 0 || minSIdx > DNS - 1)
			{
				continue;
			}

			if (minSIdx < 0)
			{
				minSIdx = 0;
			}
			if (maxSIdx > DNS - 1)
			{
				maxSIdx = DNS - 1;
			}
			eleN = 0;

			for (int sIdx = minSIdx; sIdx <= maxSIdx; sIdx++)
			{
				sk = (sIdx - detCntIdS) * ds;
				Store1[eleN++] = F1(sk, ds, tau);
			}

			curphi = beta + atan((curobjX * cosBeta + curobjY * sinBeta) / (S2O - (-curobjX * sinBeta + curobjY * cosBeta)));
			l_phi0 = dx / std::max(abs(cos(curphi)), abs(sin(curphi)));
			//If use SF-TR-A2, calculate the l_phi_0
			for (size_t kk = 0; kk < ZN; kk++)
			{
				//计算上下两个点对应的坐标 t;
				dow_t = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-curobjX * sinBeta + curobjY * cosBeta));
				upp_t = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-curobjX * sinBeta + curobjY * cosBeta));
				//计算 T Idx
				minTIdx = std::floor(dow_t / dt + detCntIdT) - 1;
				maxTIdx = std::ceil(upp_t / dt + detCntIdT) + 1;
				if (maxTIdx < 0 || minTIdx > DNT - 1)
				{
					continue;
				}
				if (minTIdx < 0)
				{
					minTIdx = 0;
				}
				if (maxTIdx > DNT - 1)
				{
					maxTIdx = DNT - 1;
				}

				xi_dow[0] = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-(curobjX - 0.5 * dx) * sinBeta + (curobjY - 0.5 * dy) * cosBeta));
				xi_dow[1] = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-(curobjX - 0.5 * dx) * sinBeta + (curobjY + 0.5 * dy) * cosBeta));
				xi_dow[2] = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-(curobjX + 0.5 * dx) * sinBeta + (curobjY - 0.5 * dy) * cosBeta));
				xi_dow[3] = (kk - objCntIdxZ - 0.5) * dz * S2D / (S2O - (-(curobjX + 0.5 * dx) * sinBeta + (curobjY + 0.5 * dy) * cosBeta));
				sort_Tau(xi_dow);
				xi_upp[0] = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-(curobjX - 0.5 * dx) * sinBeta + (curobjY - 0.5 * dy) * cosBeta));
				xi_upp[1] = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-(curobjX - 0.5 * dx) * sinBeta + (curobjY + 0.5 * dy) * cosBeta));
				xi_upp[2] = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-(curobjX + 0.5 * dx) * sinBeta + (curobjY - 0.5 * dy) * cosBeta));
				xi_upp[3] = (kk - objCntIdxZ + 0.5) * dz * S2D / (S2O - (-(curobjX + 0.5 * dx) * sinBeta + (curobjY + 0.5 * dy) * cosBeta));
				sort_Tau(xi_upp);
				xi[0] = xi_dow[0];
				xi[1] = xi_dow[3];
				xi[2] = xi_upp[0];
				xi[3] = xi_upp[3];

				summ = 0;
				for (int tIdx = minTIdx; tIdx <= maxTIdx; tIdx++)
				{
					tl = (tIdx - detCntIdT) * dt;
					f2val = F2_TT(tl, dt, xi);
					for (int sIdx = minSIdx; sIdx <= maxSIdx; sIdx++)
					{
						cur_s = (sIdx - detCntIdS) * ds;
						cur_t = (tIdx - detCntIdT) * dt;
						curTheta = -std::atan(cur_t / hypot(cur_s, S2D));
						l_theta = 1.0 / abs(std::cos(curTheta));
						summ += proj[(angIdx * DNT + tIdx) * DNS + sIdx] * (l_theta) * Store1[sIdx - minSIdx] * f2val * l_phi0;
					}
				}
				vol[(kk * YN + jj) * XN + ii] = summ;
			}
		}
	}
}



void SFTTA2_3D_ED_bproj(const std::vector<float>& proj, std::vector<float>& vol,
	const float S2O, const float O2D,
	const float objSizeX, const float objSizeY, const float objSizeZ,
	const float detSizeS, const float detSizeT,
	const float detCntIdS, const float detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<float>& angs, const int angIdx)
{
	SFTTA2_3D_ED_bproj_template<float>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ, detSizeS, detSizeT,
		detCntIdS, detCntIdT, XN, YN, ZN, DNS, DNT, PN, angs, angIdx);
}


void SFTTA2_3D_ED_bproj(const std::vector<double>& proj, std::vector<double>& vol,
	const double S2O, const double O2D,
	const double objSizeX, const double objSizeY, const double objSizeZ,
	const double detSizeS, const double detSizeT,
	const double detCntIdS, const double detCntIdT,
	const int XN, const int YN, const int ZN, const int DNS, const int DNT, const int PN,
	const std::vector<double>& angs, const int angIdx)
{
	SFTTA2_3D_ED_bproj_template<double>(proj, vol, S2O, O2D, objSizeX, objSizeY, objSizeZ, detSizeS, detSizeT,
		detCntIdS, detCntIdT, XN, YN, ZN, DNS, DNT, PN, angs, angIdx);
}
