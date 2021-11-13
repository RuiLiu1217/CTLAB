#include "differentialFuncs.hpp"
#include "utilities.hpp"


template<typename T>
void _Dx_CPU(T* f, T* d, const T coef, cuint L, cuint W)
{
	unsigned int idx(0), idy(0), curid(0), lasid(0);
	for (idy = 0; idy != W; ++idy)
	{
		for (idx = 0; idx != L; ++idx)
		{
			curid = idy * L + idx;
			lasid = idy * L + ((idx + L - 1) % L);
			d[curid] = (f[curid] - f[lasid]) * coef;
		}
	}
}


template<typename T>
void _D2x_CPU(T* f, T* d, const T coef, cuint L, cuint W)
{
	unsigned int idx(0), idy(0), curid(0), lasid(0), las2id(0);
	for (idy = 0; idy != W; ++idy)
	{
		for (idx = 0; idx != L; ++idx)
		{
			curid = idy * L + idx;
			lasid = idy * L + ((idx + L - 1) % L);
			las2id = idy * L + ((idx + L - 2) % L);
			d[curid] = (f[curid] + f[las2id] - f[lasid] * 2.0f) * coef;
		}
	}
}


template<typename T>
void _Dx_CPU(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	unsigned int idx(0), idy(0), idz(0), curid(0), lasid(0);
	for (idz = 0; idz != H; ++idz)
	{
		for (idy = 0; idy != W; ++idy)
		{
			for (idx = 0; idx != L; ++idx)
			{
				curid = (idz * W + idy) * L + idx;
				lasid = (idz * W + idy) * L + ((idx + L - 1) % L);
				d[curid] = (f[curid] - f[lasid]) * coef;
			}
		}
	}
}


template<typename T>
void _D2x_CPU(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	unsigned int idx(0), idy(0), idz(0), curid(0), lasid(0), las2id(0);
	for (idz = 0; idz != H; ++idz)
	{
		for (idy = 0; idy != W; ++idy)
		{
			for (idx = 0; idx != L; ++idx)
			{
				curid = (idz * W + idy) * L + idx;
				lasid = (idz * W + idy) * L + ((idx + L - 1) % L);
				las2id = (idz * W + idy) * L + ((idx + L - 2) % L);
				d[curid] = (f[curid] + f[las2id] - f[lasid] * 2.0f) * coef;
			}
		}
	}
}




template<typename T>
void _Dxt_CPU(T* f, T* d, const T coef, cuint L, cuint W)
{
	unsigned int idx(0), idy(0), curid(0), nexid(0);
	for (idy = 0; idy != W; ++idy)
	{
		for (idx = 0; idx != L; ++idx)
		{
			curid = idy * L + idx;
			nexid = idy * L + ((idx + 1) % L);
			d[curid] = (f[curid] - f[nexid]) * coef;
		}
	}
}


template<typename T>
void _D2xt_CPU(T* f, T* d, const T coef, cuint L, cuint W)
{
	unsigned int idx(0), idy(0), curid(0), nexid(0), nex2id(0);
	for (idy = 0; idy != W; ++idy)
	{
		for (idx = 0; idx != L; ++idx)
		{
			curid = idy * L + idx;
			nexid = idy * L + ((idx + 1) % L);
			nex2id = idy * L + ((idx + 2) % L);
			d[curid] = (f[curid] + f[nex2id] - f[nexid] * 2.0f) * coef;
		}
	}
}



template<typename T>
void _Dxt_CPU(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	unsigned int idx(0), idy(0), idz(0), curid(0), nexid(0);
	for (idz = 0; idz != H; ++idz)
	{
		for (idy = 0; idy != W; ++idy)
		{
			for (idx = 0; idx != L; ++idx)
			{
				curid = (idz * W + idy) * L + idx;
				nexid = (idz * W + idy) * L + ((idx + 1) % L);
				d[curid] = (f[curid] - f[nexid]) * coef;
			}
		}
	}
}

template<typename T>
void _D2xt_CPU(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	unsigned int idx(0), idy(0), idz(0), curid(0), nexid(0), nex2id(0);
	for (idz = 0; idz != H; ++idz)
	{
		for (idy = 0; idy != W; ++idy)
		{
			for (idx = 0; idx != L; ++idx)
			{
				curid = (idz * W + idy) * L + idx;
				nexid = (idz * W + idy) * L + ((idx + 1) % L);
				nex2id = (idz * W + idy) * L + ((idx + 2) % L);
				d[curid] = (f[curid] + f[nex2id] - f[nexid] * 2.0f) * coef;
			}
		}
	}
}



template<typename T>
void _Dy_CPU(T* f, T* d, const T coef, cuint L, cuint W)
{
	unsigned int idx(0), idy(0), curid(0), lasid(0);
	for (idy = 0; idy != W; ++idy)
	{
		for (idx = 0; idx != L; ++idx)
		{
			curid = idy * L + idx;
			lasid = ((idy + W - 1) % W) * L + idx;
			d[curid] = (f[curid] - f[lasid]) * coef;
		}
	}
}


template<typename T>
void _D2y_CPU(T* f, T* d, const T coef, cuint L, cuint W)
{
	unsigned int idx(0), idy(0), curid(0), lasid(0), las2id(0);
	for (idy = 0; idy != W; ++idy)
	{
		for (idx = 0; idx != L; ++idx)
		{
			curid = idy * L + idx;
			lasid = ((idy + W - 1) % W) * L + idx;
			las2id = ((idy + W - 2) % W) * L + idx;
			d[curid] = (f[curid] + f[las2id] - f[lasid] * 2.0f) * coef;

		}
	}
}

template<typename T>
void _Dy_CPU(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	unsigned int idx(0), idy(0), idz(0), curid(0), lasid(0);
	for (idz = 0; idz != H; ++idz)
	{
		for (idy = 0; idy != W; ++idy)
		{
			for (idx = 0; idx != L; ++idx)
			{
				curid = (idz * W + idy) * L + idx;
				lasid = (idz * W + ((idy + W - 1) % W)) * L + idx;
				d[curid] = (f[curid] - f[lasid]) * coef;
			}
		}
	}
}


template<typename T>
void _D2y_CPU(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	unsigned int idx(0), idy(0), idz(0), curid(0), lasid(0), las2id(0);
	for (idz = 0; idz != H; ++idz)
	{
		for (idy = 0; idy != W; ++idy)
		{
			for (idx = 0; idx != L; ++idx)
			{
				curid = (idz * W + idy) * L + idx;
				lasid = (idz * W + ((idy + W - 1) % W)) * L + idx;
				las2id = (idz * W + ((idy + W - 2) % W)) * L + idx;
				d[curid] = (f[curid] + f[las2id] - f[lasid] * 2.0f) * coef;
			}
		}
	}
}



template<typename T>
void _Dyt_CPU(T* f, T* d, const T coef, cuint L, cuint W)
{
	unsigned int idx(0), idy(0), curid(0), nexid(0);
	for (idy = 0; idy != W; ++idy)
	{
		for (idx = 0; idx != L; ++idx)
		{
			curid = idy * L + idx;
			nexid = ((idy + 1) % W) * L + idx;
			d[curid] = (f[curid] - f[nexid]) * coef;
		}
	}
}


template<typename T>
void _D2yt_CPU(T* f, T* d, const T coef, cuint L, cuint W)
{
	unsigned int idx(0), idy(0), curid(0), nexid(0), nex2id(0);
	for (idy = 0; idy != W; ++idy)
	{
		for (idx = 0; idx != L; ++idx)
		{
			curid = idy * L + idx;
			nexid = ((idy + 1) % W) * L + idx;
			nex2id = ((idy + 2) % W) * L + idx;
			d[curid] = (f[curid] + f[nex2id] - f[nexid] * 2.0f) * coef;
		}
	}
}


template<typename T>
void _Dyt_CPU(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	unsigned int idx(0), idy(0), idz(0), curid(0), nexid(0);
	for (idz = 0; idz != H; ++idz)
	{
		for (idy = 0; idy != W; ++idy)
		{
			for (idx = 0; idx != L; ++idx)
			{
				curid = (idz * W + idy) * L + idx;
				nexid = (idz * W + ((idy + 1) % W)) * L + idx;
				d[curid] = (f[curid] - f[nexid]) * coef;
			}
		}
	}
}

template<typename T>
void _D2yt_CPU(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	unsigned int idx(0), idy(0), idz(0), curid(0), nexid(0), nex2id(0);
	for (idz = 0; idz != H; ++idz)
	{
		for (idy = 0; idy != W; ++idy)
		{
			for (idx = 0; idx != L; ++idx)
			{
				curid = (idz * W + idy) * L + idx;
				nexid = (idz * W + ((idy + 1) % W)) * L + idx;
				nex2id = (idz * W + ((idy + 2) % W)) * L + idx;
				d[curid] = (f[curid] + f[nex2id] - f[nexid] * 2.0f) * coef;
			}
		}
	}
}


template<typename T>
void _Dz_CPU(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	unsigned int idx(0), idy(0), idz(0), curid(0), lasid(0);
	for (idz = 0; idz != H; ++idz)
	{
		for (idy = 0; idy != W; ++idy)
		{
			for (idx = 0; idx != L; ++idx)
			{
				curid = (idz * W + idy) * L + idx;
				lasid = (((idz + H - 1) % H) * W + idy) * L + idx;
				d[curid] = (f[curid] - f[lasid]) * coef;
			}
		}
	}

}



template<typename T>
void _D2z_CPU(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	unsigned int idx(0), idy(0), idz(0), curid(0), lasid(0), las2id(0);
	for (idz = 0; idz != H; ++idz)
	{
		for (idy = 0; idy != W; ++idy)
		{
			for (idx = 0; idx != L; ++idx)
			{
				curid = (idz * W + idy) * L + idx;
				lasid = (((idz + H - 1) % H) * W + idy) * L + idx;
				las2id = (((idz + H - 2) % H) * W + idy) * L + idx;
				d[curid] = (f[curid] + f[las2id] - f[lasid] * 2.0f) * coef;
			}
		}
	}
}


template<typename T>
void _Dzt_CPU(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	unsigned int idx(0), idy(0), idz(0), curid(0), nexid(0);
	for (idz = 0; idz != H; ++idz)
	{
		for (idy = 0; idy != W; ++idy)
		{
			for (idx = 0; idx != L; ++idx)
			{
				curid = (idz * W + idy) * L + idx;
				nexid = (((idz + 1) % H) * W + idy) * L + idx;
				d[curid] = (f[curid] - f[nexid]) * coef;
			}
		}
	}
}

template<typename T>
void _D2zt_CPU(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	unsigned int idx(0), idy(0), idz(0), curid(0), nexid(0), nex2id(0);
	for (idz = 0; idz != H; ++idz)
	{
		for (idy = 0; idy != W; ++idy)
		{
			for (idx = 0; idx != L; ++idx)
			{
				curid = (idz * W + idy) * L + idx;
				nexid = (((idz + 1) % H) * W + idy) * L + idx;
				nex2id = (((idz + 2) % H) * W + idy) * L + idx;
				d[curid] = (f[curid] + f[nex2id] - f[nexid] * 2.0f) * coef;
			}
		}
	}
}


template<typename T>
void _Laplacian_CPU(T* f, T* l, const T coef, cuint L, cuint W)
{
	unsigned int idx(0), idy(0), curid(0), lasxd(0), nexxd(0), lasyd(0), nexyd(0);
	for (idy = 0; idy != W; ++idy)
	{
		for (idx = 0; idx != L; ++idx)
		{
			curid = idy * L + idx;
			lasxd = idy * L + ((idx + L - 1) % L);
			nexxd = idy * L + ((idx + 1) % L);
			lasyd = ((idy + W - 1) % W) * L + idx;
			nexyd = ((idy + 1) % W) * L + idx;
			l[curid] = (4.0f * f[curid] - f[lasxd] - f[lasyd] - f[nexxd] - f[nexyd]) * coef;
		}
	}
}

template<typename T>
void _Laplacian_CPU(T* f, T* l, const T coef, cuint L, cuint W, cuint H)
{
	unsigned int idx(0), idy(0), idz(0),
		curid(0), lasxd(0), lasyd(0), laszd(0), nexxd(0), nexyd(0), nexzd(0);
	for (idz = 0; idz != H; ++idz)
	{
		for (idy = 0; idy != W; ++idy)
		{
			for (idx = 0; idy != L; ++idx)
			{
				curid = (idz * W + idy) * L + idx;
				lasxd = (idz * W + idy) * L + ((idx + L - 1) % L);
				nexxd = (idz * W + idy) * L + ((idx + 1) % L);
				lasyd = (idz * W + ((idy + W - 1) % W)) * L + idx;
				nexyd = (idz * W + ((idy + 1) % W)) * L + idx;
				laszd = (((idz + H - 1) % H) * W + idy) * L + idx;
				nexzd = (((idz + 1) % H) * W + idy) * L + idx;
				l[curid] = (6.0f * f[curid] - f[lasxd] - f[lasyd] - f[nexxd] - f[nexyd] - f[laszd] - f[nexzd]) * coef;
			}
		}
	}
}


template<typename T>
void _DiscreteGradientTrans_CPU(T* f, T* d, const T coef, cuint L, cuint W)
{
	unsigned int idx(0), idy(0), curid(0), nexxd(0), nexyd(0);
	T difx(0.0), dify(0.0);
	for (idy = 0; idy != W; ++idy)
	{
		for (idx = 0; idx != L; ++idx)
		{
			curid = idy * L + idx;
			nexxd = idy * L + ((idx + 1) % L);
			nexyd = ((idy + 1) % W) * L + idx;
			difx = f[nexxd] - f[curid];
			dify = f[nexyd] - f[curid];
			d[curid] = std::sqrt(difx * difx + dify * dify) * coef;
		}
	}
}



template<typename T>
void _DiscreteGradientTrans_CPU(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	unsigned int idx(0), idy(0), idz(0), curid(0), nexxd(0), nexyd(0), nexzd(0);
	T difx(0), dify(0), difz(0);
	for (idz = 0; idz != H; ++idz)
	{
		for (idy = 0; idy != W; ++idy)
		{
			for (idx = 0; idx != L; ++idx)
			{
				curid = (idz * W + idy) * L + idx;
				nexxd = (idz * W + idy) * L + ((idx + 1) % L);
				nexyd = (idz * W + ((idy + 1) % W)) * L + idx;
				nexzd = (((idz + 1) % H) * W + idy) * L + idx;
				difx = f[nexxd] - f[curid];
				dify = f[nexyd] - f[curid];
				difz = f[nexzd] - f[curid];
				d[curid] = std::sqrt(difx * difx + dify * dify + difz * difz) * coef;
			}
		}
	}
}


template<typename T>
void _GradientOfTV_CPU(T* f, T* d, const T coef, cuint L, cuint W)
{

	unsigned int idx(0), idy(0), curid(0);
	T fij(0), fi1j(0), fi_1j(0), fij1(0),
		fij_1(0), fi_1j1(0), fi1j_1(0), dom1(0), dom2(0), dom3(0);
	for (idy = 0; idy != W; ++idy)
	{
		for (idx = 0; idx != L; ++idx)
		{
			curid = idy * L + idx;
			fij = f[curid];
			fi1j = f[idy * L + (idx + 1) % L];
			fi_1j = f[idy * L + (idx + L - 1) % L];
			fij1 = f[((idy + 1) % W) * L + idx];
			fij_1 = f[((idy + W - 1) % W) * L + idx];
			fi_1j1 = f[((idy + 1) % W) * L + (idx + L - 1) % L];
			fi1j_1 = f[((idy + W - 1) % W) * L + (idx + 1) % L];
			dom1 = 1.0f / (std::sqrt((fi1j - fij)*(fi1j - fij) + (fij1 - fij) * (fij1 - fij) + _EPSILON));
			dom2 = 1.0f / (std::sqrt((fij - fi_1j)*(fij - fi_1j) + (fi_1j1 - fi_1j)*(fi_1j1 - fi_1j) + _EPSILON));
			dom3 = 1.0f / (std::sqrt((fi1j_1 - fij_1) * (fi1j_1 - fij_1) + (fij - fij_1)*(fij - fij_1) + _EPSILON));
			d[curid] = ((2.0 * fij - fi1j - fij1) * dom1 + (fij - fi_1j) * dom2 + (fij - fij_1) * dom3) * coef;
		}
	}
}



template<typename T>
void _GradientOfTV_CPU(T* f, T* d, const T coef, cuint L, cuint W, cuint H)
{
	unsigned int idx(0), idy(0), idz(0), curId(0);
	T f_ijk(0), f_i1jk(0), f_ij1k(0), f_ijk1(0), f_i_1jk(0), f_i_1j1k(0),
		f_i_1jk1(0), f_ij_1k(0), f_i1j_1k(0), f_ij_1k1(0),
		f_ijk_1(0), f_i1jk_1(0), f_ij1k_1(0),
		dom1(0), dom2(0), dom3(0), dom4(0);
	for (idz = 0; idz != H; ++idz)
	{
		for (idy = 0; idy != W; ++idy)
		{
			for (idx = 0; idx != L; ++idx)
			{
				curId = (idz * W + idy) * L + idx;
				f_ijk = f[curId];
				f_i1jk = f[(idz * W + idy) * L + ((idx + 1) % L)];
				f_ij1k = f[(idz * W + ((idy + 1) % W)) * L + idx];
				f_ijk1 = f[(((idz + 1) % H) * W + idy) * L + idx];
				f_i_1jk = f[(idz * W + idy) * L + ((idx + L - 1) % L)];
				f_i_1j1k = f[(idz * W + ((idy + 1) % W)) * L + ((idx + L - 1) % L)];
				f_i_1jk1 = f[(((idz + 1) % H) * W + idy) * L + ((idx + L - 1) % L)];
				f_ij_1k = f[(idz * W + ((idy + W - 1) % W)) * L + idx];
				f_i1j_1k = f[(idz * W + ((idy + W - 1) % W)) * L + ((idx + 1) % L)];
				f_ij_1k1 = f[(((idz + 1) % H) * W + ((idy + W - 1) % W)) * L + idx];
				f_ijk_1 = f[(((idz + H - 1) % H) * W + idy) * L + idx];
				f_i1jk_1 = f[(((idz + H - 1) % H) * W + idy) * L + ((idx + 1) % L)];
				f_ij1k_1 = f[(((idz + H - 1) % H) * W + ((idy + 1) % W)) * L + idx];

				dom1 = 1.0f / std::sqrt((f_i1jk - f_ijk)*(f_i1jk - f_ijk) + (f_ij1k - f_ijk)*(f_ij1k - f_ijk) + (f_ijk1 - f_ijk)*(f_ijk1 - f_ijk) + _EPSILON);
				dom2 = 1.0f / std::sqrt((f_ijk - f_i_1jk)*(f_ijk - f_i_1jk) + (f_i_1j1k - f_i_1jk)*(f_i_1j1k - f_i_1jk) + (f_i_1jk1 - f_i_1jk)*(f_i_1jk1 - f_i_1jk) + _EPSILON);
				dom3 = 1.0f / std::sqrt((f_i1j_1k - f_ij_1k)*(f_i1j_1k - f_ij_1k) + (f_ijk - f_ij_1k)*(f_ijk - f_ij_1k) + (f_ij_1k1 - f_ij_1k)*(f_ij_1k1 - f_ij_1k) + _EPSILON);
				dom4 = 1.0f / std::sqrt((f_i1jk_1 - f_ijk_1)*(f_i1jk_1 - f_ijk_1) + (f_ij1k_1 - f_ijk_1)*(f_ij1k_1 - f_ijk_1) + (f_ijk - f_ijk_1)*(f_ijk - f_ijk_1) + _EPSILON);

				d[curId] = ((3.0f * f_ijk - f_i1jk - f_ij1k - f_ijk1) * dom1 + (f_ijk - f_i_1jk) * dom2 + (f_ijk - f_ij_1k) * dom3 + (f_ijk - f_ijk_1) * dom4) * coef;
				return;
			}
		}
	}
}







void Dx_CPU(float* f, float* d, const float coef, cuint L, cuint W)
{
	_Dx_CPU<float>(f, d, coef, L, W);
}

void Dx_CPU(double* f, double* d, const double coef, cuint L, cuint W)
{
	_Dx_CPU<double>(f, d, coef, L, W);
}

void Dx_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H)
{
	_Dx_CPU<float>(f, d, coef, L, W, H);
}

void Dx_CPU(double* f, double* d, const double coef, cuint L, cuint W, cuint H)
{
	_Dx_CPU<double>(f, d, coef, L, W, H);
}

void D2x_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H)
{
	_D2x_CPU<float>(f, d, coef, L, W, H);
}

void D2x_CPU(double* f, double* d, const double coef, cuint L, cuint W, cuint H)
{
	_D2x_CPU<double>(f, d, coef, L, W, H);
}

void D2x_CPU(float* f, float* d, const float coef, cuint L, cuint W)
{
	_D2x_CPU<float>(f, d, coef, L, W);
}

void D2x_CPU(double* f, double* d, const double coef, cuint L, cuint W)
{
	_D2x_CPU<double>(f, d, coef, L, W);
}


void Dxt_CPU(float* f, float* d, const float coef, cuint L, cuint W)
{
	_Dxt_CPU<float>(f, d, coef, L, W);
}

void Dxt_CPU(double* f, double* d, const double coef, cuint L, cuint W)
{
	_Dxt_CPU<double>(f, d, coef, L, W);
}

void Dxt_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H)
{
	_Dxt_CPU<float>(f, d, coef, L, W, H);
}

void Dxt_CPU(double* f, double* d, const double coef, cuint L, cuint W, cuint H)
{
	_Dxt_CPU<double>(f, d, coef, L, W, H);
}

void D2xt_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H)
{
	_D2xt_CPU<float>(f, d, coef, L, W, H);
}

void D2xt_CPU(double* f, double* d, const double coef, cuint L, cuint W, cuint H)
{
	_D2xt_CPU<double>(f, d, coef, L, W, H);
}

void D2xt_CPU(float* f, float* d, const float coef, cuint L, cuint W)
{
	_D2xt_CPU<float>(f, d, coef, L, W);
}

void D2xt_CPU(double* f, double* d, const double coef, cuint L, cuint W)
{
	_D2xt_CPU<double>(f, d, coef, L, W);
}




void Dy_CPU(float* f, float* d, const float coef, cuint L, cuint W)
{
	_Dy_CPU<float>(f, d, coef, L, W);
}

void Dy_CPU(double* f, double* d, const double coef, cuint L, cuint W)
{
	_Dy_CPU<double>(f, d, coef, L, W);
}

void Dy_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H)
{
	_Dy_CPU<float>(f, d, coef, L, W, H);
}

void Dy_CPU(double* f, double* d, const double coef, cuint L, cuint W, cuint H)
{
	_Dy_CPU<double>(f, d, coef, L, W, H);
}

void D2y_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H)
{
	_D2y_CPU<float>(f, d, coef, L, W, H);
}

void D2y_CPU(double* f, double* d, const double coef, cuint L, cuint W, cuint H)
{
	_D2y_CPU<double>(f, d, coef, L, W, H);
}

void D2y_CPU(float* f, float* d, const float coef, cuint L, cuint W)
{
	_D2y_CPU<float>(f, d, coef, L, W);
}

void D2y_CPU(double* f, double* d, const double coef, cuint L, cuint W)
{
	_D2y_CPU<double>(f, d, coef, L, W);
}


void Dyt_CPU(float* f, float* d, const float coef, cuint L, cuint W)
{
	_Dyt_CPU<float>(f, d, coef, L, W);
}

void Dyt_CPU(double* f, double* d, const double coef, cuint L, cuint W)
{
	_Dyt_CPU<double>(f, d, coef, L, W);
}

void Dyt_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H)
{
	_Dyt_CPU<float>(f, d, coef, L, W, H);
}

void Dyt_CPU(double* f, double* d, const double coef, cuint L, cuint W, cuint H)
{
	_Dyt_CPU<double>(f, d, coef, L, W, H);
}

void D2yt_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H)
{
	_D2yt_CPU<float>(f, d, coef, L, W, H);
}

void D2yt_CPU(double* f, double* d, const double coef, cuint L, cuint W, cuint H)
{
	_D2yt_CPU<double>(f, d, coef, L, W, H);
}

void D2yt_CPU(float* f, float* d, const float coef, cuint L, cuint W)
{
	_D2yt_CPU<float>(f, d, coef, L, W);
}

void D2yt_CPU(double* f, double* d, const double coef, cuint L, cuint W)
{
	_D2yt_CPU<double>(f, d, coef, L, W);
}



void Dz_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H)
{
	_Dz_CPU<float>(f, d, coef, L, W, H);
}

void Dz_CPU(double* f, double* d, const double coef, cuint L, cuint W, cuint H)
{
	_Dz_CPU<double>(f, d, coef, L, W, H);
}

void D2z_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H)
{
	_D2z_CPU<float>(f, d, coef, L, W, H);
}

void D2z_CPU(double* f, double* d, const double coef, cuint L, cuint W, cuint H)
{
	_D2z_CPU<double>(f, d, coef, L, W, H);
}

void Dzt_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H)
{
	_Dzt_CPU<float>(f, d, coef, L, W, H);
}

void Dzt_CPU(double* f, double* d, const double coef, cuint L, cuint W, cuint H)
{
	_Dzt_CPU<double>(f, d, coef, L, W, H);
}

void D2zt_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H)
{
	_D2zt_CPU<float>(f, d, coef, L, W, H);
}

void D2zt_CPU(double* f, double* d, const double coef, cuint L, cuint W, cuint H)
{
	_D2zt_CPU<double>(f, d, coef, L, W, H);
}



void Laplacian_CPU(float*f, float* l, const float coef, cuint L, cuint W)
{
	_Laplacian_CPU<float>(f, l, coef, L, W);
}

void Laplacian_CPU(double*f, double* l, const double coef, cuint L, cuint W)
{
	_Laplacian_CPU<double>(f, l, coef, L, W);
}

void Laplacian_CPU(float* f, float* l, const float coef, cuint L, cuint W, cuint H)
{
	_Laplacian_CPU<float>(f, l, coef, L, W, H);
}

void Laplacian_CPU(double* f, double* l, const double coef, cuint L, cuint W, cuint H)
{
	_Laplacian_CPU<double>(f, l, coef, L, W, H);
}


void DiscreteGradientTrans_CPU(float* f, float* d, const float coef, cuint L, cuint W)
{
	_DiscreteGradientTrans_CPU<float>(f, d, coef, L, W);
}
void DiscreteGradientTrans_CPU(double* f, double* d, const double coef, cuint L, cuint W)
{
	_DiscreteGradientTrans_CPU<double>(f, d, coef, L, W);
}

void DiscreteGradientTrans_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H)
{
	_DiscreteGradientTrans_CPU<float>(f, d, coef, L, W, H);
}
void DiscreteGradientTrans_CPU(double* f, double* d, const double coef, cuint L, cuint W, cuint H)
{
	_DiscreteGradientTrans_CPU<double>(f, d, coef, L, W, H);
}

void GradientOfTV_CPU(float* f, float* d, const float coef, cuint L, cuint W)
{
	_GradientOfTV_CPU<float>(f, d, coef, L, W);
}
void GradientOfTV_CPU(double* f, double* d, const double coef, cuint L, cuint W)
{
	_GradientOfTV_CPU<double>(f, d, coef, L, W);
}


void GradientOfTV_CPU(float* f, float* d, const float coef, cuint L, cuint W, cuint H)
{
	_GradientOfTV_CPU<float>(f, d, coef, L, W, H);
}

void GradientOfTV_CPU(double* f, double* d, const double coef, cuint L, cuint W, cuint H)
{
	_GradientOfTV_CPU<double>(f, d, coef, L, W, H);
}



template<typename T>
void _invDiscreteGradientTransform_CPU(T* f, T* d, T* r, const T omega, cuint L, cuint W)
{
	unsigned int idx(0), idy(0), curid(0), lasxd(0), nexxd(0), lasyd(0), nexyd(0);
	T fa(0), fb(0), fc(0);
	for (idy = 0; idy != W; ++idy)
	{
		for (idx = 0; idx != L; ++idx)
		{
			curid = idy * L + idx;
			lasxd = idy * L + (idx + L - 1) % L;
			nexxd = idy * L + (idx + 1) % L;
			lasyd = ((idy + W - 1) % W) * L + idx;
			nexyd = ((idy + 1) % W) * L + idx;

			if (d[curid] < omega)
			{
				fa = (2.0f * f[curid] + f[nexxd] + f[nexyd]) * 0.25f;
			}
			else
			{
				fa = f[curid] - omega * (2.0f * f[curid] - f[nexxd] - f[nexyd]) / (4.0f * d[curid]);
			}
			if (d[lasxd] < omega)
			{
				fb = (f[curid] + f[lasxd]) * 0.5f;
			}
			else
			{
				fb = f[curid] - omega * (f[curid] - f[lasxd]) * 0.5f / d[lasxd];
			}
			if (d[lasyd] < omega)
			{
				fc = (f[curid] + f[lasyd]) * 0.5f;
			}
			else
			{
				fc = f[curid] - omega * (f[curid] - f[lasyd]) * 0.5f / d[lasyd];
			}
			r[curid] = 0.25f * (2.0f * fa + fb + fc);
		}
	}
}

void invDiscreteGradientTransform_CPU(float* f, float* d, float* r, const float omega, cuint L, cuint W)
{
	_invDiscreteGradientTransform_CPU<float>(f, d, r, omega, L, W);
}
void invDiscreteGradientTransform_CPU(double* f, double* d, double* r, const double omega, cuint L, cuint W)
{
	_invDiscreteGradientTransform_CPU<double>(f, d, r, omega, L, W);
}


template<typename T>
void _invDiscreteGradientTransform_CPU(T* f, T* d, T* r, const T omega, cuint L, cuint W, cuint H)
{
	unsigned int idx(0), idy(0), idz(0), curId(0),
		nexxd(0), lasxd(0), nexyd(0), lasyd(0), nexzd(0), laszd(0);
	T fa(0), fb(0), fc(0), fd(0);
	for (idz = 0; idz != H; ++idz)
	{
		for (idy = 0; idy != W; ++idy)
		{
			for (idx = 0; idx != L; ++idx)
			{
				curId = (idz * W + idy) * L + idx;
				nexxd = (idz * W + idy) * L + ((idx + 1) % L);
				lasxd = (idz * W + idy) * L + ((idx + L - 1) % L);
				nexyd = (idz * W + ((idy + 1) % W)) * L + idx;
				lasyd = (idz * W + ((idy + W - 1) % W)) * L + idx;
				nexzd = (((idz + 1) % H) * W + idy) * L + idx;
				laszd = (((idz + H - 1) % H) * W + idy) * L + idx;
				if (d[curId] < omega)
				{
					fa = (3.0f * f[curId] + f[nexxd] + f[nexyd] + f[nexzd]) * 0.1666666666666667f;
				}
				else
				{
					fa = f[curId] - omega *(3.0f * f[curId] - f[nexxd] - f[nexyd] - f[nexzd]) * 0.1666666666666667f / d[curId];
				}
				if (d[lasxd] < omega)
				{
					fb = (f[curId] + f[lasxd]) * 0.5f;
				}
				else
				{
					fb = f[curId] - omega * (f[curId] - f[lasxd]) * 0.5f / d[lasxd];
				}
				if (d[lasyd] < omega)
				{
					fc = (f[curId] + f[lasyd]) * 0.5f;
				}
				else
				{
					fc = f[curId] - omega * (f[curId] - f[lasyd]) * 0.5f / d[lasyd];
				}
				if (d[laszd] < omega)
				{
					fd = (f[curId] + f[laszd]) * 0.5f;
				}
				else
				{
					fd = f[curId] - omega * (f[curId] - f[laszd]) * 0.5f / d[laszd];
				}
				r[curId] = 0.1666666666666667f * (3.0f * fa + fb + fc + fd);
			}
		}
	}
}
void invDiscreteGradientTransform_CPU(float* f, float* d, float* r, const float omega, cuint L, cuint W, cuint H)
{
	_invDiscreteGradientTransform_CPU<float>(f, d, r, omega, L, W, H);
}
void invDiscreteGradientTransform_CPU(double* f, double* d, double* r, const double omega, cuint L, cuint W, cuint H)
{
	_invDiscreteGradientTransform_CPU<double>(f, d, r, omega, L, W, H);
}




template<typename T>
void _rawToUnchar(
	std::vector<T>& inData,
	std::vector<unsigned char>& ouData)
{

	auto result = std::minmax_element(inData.begin(), inData.end());
	unsigned int minPos(result.first - inData.begin());
	unsigned int maxPos(result.second - inData.begin());
	T minData = inData[minPos];
	T maxData = inData[maxPos];
	ouData.clear();
	for (unsigned int i = 0; i != inData.size(); ++i)
	{
		ouData.push_back(static_cast<unsigned char>(static_cast<double>(inData[i] - minData) / (maxData - minData) * 255.0));
	}
}


void rawToUnchar(std::vector<float>& inData, std::vector<unsigned char>& ouData)
{
	_rawToUnchar<float>(inData, ouData);
}
void rawToUnchar(std::vector<double>& inData, std::vector<unsigned char>& ouData)
{
	_rawToUnchar<double>(inData, ouData);
}
void rawToUnchar(float* inData, const unsigned int len, std::vector<unsigned char>& ouData)
{
	float minData = inData[0];
	float maxData = inData[0];
	for (unsigned int idx = 0; idx != len; ++idx)
	{
		if (inData[idx] < minData)
		{
			minData = inData[idx];
		}

		if (inData[idx] > maxData)
		{
			maxData = inData[idx];
		}
	}
	ouData.clear();
	if (maxData != minData)
	{
		for (unsigned int idx = 0; idx != len; ++idx)
		{
			ouData.push_back(static_cast<unsigned char>(static_cast<double>(inData[idx] - minData) / (maxData - minData) * 255.0));
		}
	}
	else
	{
		for (unsigned int idx = 0; idx != len; ++idx)
		{
			ouData.push_back(0);
		}
	}
}


void rawToUnchar(double* inData, const unsigned int len, std::vector<unsigned char>& ouData)
{
	double minData = inData[0];
	double maxData = inData[0];
	for (unsigned int idx = 0; idx != len; ++idx)
	{
		if (inData[idx] < minData)
		{
			minData = inData[idx];
		}

		if (inData[idx] > maxData)
		{
			maxData = inData[idx];
		}
	}
	ouData.clear();
	if (maxData != minData)
	{
		for (unsigned int idx = 0; idx != len; ++idx)
		{
			ouData.push_back(static_cast<unsigned char>(static_cast<double>(inData[idx] - minData) / (maxData - minData) * 255.0));
		}
	}
	else
	{
		for (unsigned int idx = 0; idx != len; ++idx)
		{
			ouData.push_back(0);
		}
	}
}


