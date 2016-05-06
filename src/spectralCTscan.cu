#include "utilities.hpp"

#include "spectralCTscan.hpp"
enum{ZERO=0,ONE=1,NUM_OF_CYLINDER=17,NUM_OF_ENERGYCHANNEL=801};

static const double PI = 3.14159265358979323846264;
namespace SPECTRALCT{
	template<typename T>
	struct T3
	{
		T x;
		T y;
		T z;
		__host__ __device__ T3(){ x = 0; y = 0; z = 0; }
		__host__ __device__ T3(const T3<T>& v)
		{
			x = v.x;
			y = v.y;
			z = v.z;
		}
	};



	template<typename T>
	struct Cylinder
	{
		T A[NUM_OF_CYLINDER]; // Materials
		T R[NUM_OF_CYLINDER];
		T posx[NUM_OF_CYLINDER];
		T posy[NUM_OF_CYLINDER];
		T height;
		Cylinder()
		{
			height = 1;
			for (unsigned int i = 0; i < NUM_OF_CYLINDER; i++)
			{
				A[i] = i + 1;
			}
			R[0] = 0.45 * 2; posx[0] = 0; posy[0] = 0;
			R[1] = 0.075 * 2; posx[1] = 0.55; posy[1] = 0;
			R[2] = 0.075 * 2; posx[2] = 0.275; posy[2] = -0.4763;
			R[3] = 0.075 * 2; posx[3] = -0.275; posy[3] = -0.4763;
			R[4] = 0.075 * 2; posx[4] = -0.55; posy[4] = 0;
			R[5] = 0.075 * 2; posx[5] = -0.275; posy[5] = 0.4763;
			R[6] = 0.075 * 2; posx[6] = 0.275; posy[6] = 0.4763;
			R[7] = 0.04 * 2; posx[7] = 0.275; posy[7] = 0;
			R[8] = 0.035 * 2; posx[8] = 0.1375; posy[8] = -0.23815;
			R[9] = 0.03 * 2; posx[9] = -0.1375; posy[9] = -0.23815;
			R[10] = 0.025 * 2; posx[10] = -0.275; posy[10] = 0;
			R[11] = 0.02 * 2; posx[11] = -0.1375; posy[11] = 0.23815;
			R[12] = 0.015 * 2; posx[12] = 0.1375; posy[12] = 0.23815;
			R[13] = 0.01 * 2; posx[13] = 0; posy[13] = -0.1;
			R[14] = 0.0075 * 2; posx[14] = 0.1; posy[14] = 0;
			R[15] = 0.005 * 2; posx[15] = 0; posy[15] = 0.1;
			R[16] = 0.0025 * 2; posx[16] = -0.1; posy[16] = 0;

		}
	};

	template<typename T>
	inline __host__ __device__ bool intersectLengthWithCyl(const T sourx, const T soury,
		const T sourz, const T dirx, const T diry, const T dirz,
		const T ra, const T rb, // center position of the cylinder
		const T R, //radius of the cylinder
		const T botz, // minumpos of the cylinder,
		const T topz, // maximupos of the cylinder
		T& intersectlength)
	{
		const T A = (dirx * dirx + diry * diry);
		const T B = 2.0 * dirx * (sourx - ra) + 2.0 * diry * (soury - rb);
		const T C = (sourx - ra) * (sourx - ra) + (soury - rb) * (soury - rb) - R*R;
		const T Delta = B * B - 4.0 * A * C;
		if (Delta > 0)
		{
			T t1 = (-B - sqrt(Delta)) / (2.0 * A);
			T t2 = (-B + sqrt(Delta)) / (2.0 * A);
			if (dirz == 0)
			{
				intersectlength = (t2 - t1);
				return (t2 > t1);
			}
			else
			{
				T t3 = (botz - sourz) / dirz;
				T t4 = (topz - sourz) / dirz;

				T tzmin = (t3 < t4) ? t3 : t4;
				T tzmax = (t3 > t4) ? t3 : t4;
				T tinter1 = (tzmin > t1) ? tzmin : t1;
				T tinter2 = (tzmax < t2) ? tzmax : t2;
				if (tinter2 > tinter1)
				{
					intersectlength = tinter2 - tinter1;
					return tinter2 > tinter1;
				}
				else
				{
					intersectlength = 0;
					return false;
				}
			}
		}
		else
		{
			intersectlength = 0;
			return false;
		}
	}

	template<typename T>
	__global__ void projCylinder_ker(T* proj, const T sourx, const T soury, const T sourz,
		const T detcntx, const T detcnty, const T detcntz, const T ux, const T uy, const T uz,
		const T vx, const T vy, const T vz, const T detustp, const T detvstp, const int DNi, const int DNj,
		const T* A, const T* R, const T height, const T* posx, const T* posy)
	{
		const int i = threadIdx.x + blockIdx.x * blockDim.x;
		const int j = threadIdx.y + blockIdx.y * blockDim.y;

		//__shared__ T sA[NUM_OF_CYLINDER];
		__shared__ T sR[NUM_OF_CYLINDER];
		__shared__ T sposx[NUM_OF_CYLINDER];
		__shared__ T sposy[NUM_OF_CYLINDER];
		for (size_t i = 0; i < NUM_OF_CYLINDER; i++)
		{
			//sA[i] = A[i];
			sR[i] = R[i];
			sposx[i] = posx[i];
			sposy[i] = posy[i];
		}
		__syncthreads();
		if (i < DNi && j < DNj)
		{
			const T uu = (detustp * (i - DNi * 0.5 + 0.5));
			const T vv = (detvstp * (j - DNj * 0.5 + 0.5));
			const T detpx = detcntx + ux * uu + vx * vv;
			const T detpy = detcnty + uy * uu + vy * vv;
			const T detpz = detcntz + uz * uu + vz * vv;
			T dirx = detpx - sourx;
			T diry = detpy - soury;
			T dirz = detpz - sourz;
			const T ll = sqrt(dirx * dirx + diry * diry + dirz * dirz);
			dirx /= ll;
			diry /= ll;
			dirz /= ll;
			T summ = 0;
			T intersectlength = 0;
			const T botz = -height * 0.5;
			const T topz = height * 0.5;
			bool intersect = intersectLengthWithCyl<T>(sourx, soury, sourz, dirx, diry, dirz,
				sposx[0], sposy[0], sR[0], botz, topz, intersectlength);

			if (!intersect)
			{
				proj[j * DNi + i] = 0;
				return;
			}
			else
			{
				summ = intersectlength * A[0];

				for (size_t ii = 1; ii < NUM_OF_CYLINDER; ii++)
				{
					if (intersectLengthWithCyl<T>(sourx, soury, sourz, dirx, diry, dirz,
						sposx[ii], sposy[ii], sR[ii], botz, topz, intersectlength))
					{
						summ += intersectlength * (A[ii] - A[0]);
					}
				}
				proj[j * DNi + i] = summ;
			}
		}
	}

	template<typename T>
	void projCylinder(T* proj, const T sourx, const T soury, const T sourz,
		const T detcntx, const T detcnty, const T detcntz, const T ux, const T uy, const T uz,
		const T vx, const T vy, const T vz, const T detustp, const T detvstp, const int DNi, const int DNj,
		const T* A, const T* R, const T height, const T* posx, const T* posy, const dim3 blk, const dim3 gid)
	{
		projCylinder_ker<T> << <gid, blk >> >(proj, sourx, soury, sourz, detcntx, detcnty, detcntz, ux, uy, uz, vx, vy, vz, detustp, detvstp, DNi, DNj, A, R, height, posx, posy);
	}




	extern "C"
		void conebeamScanCylinder()
	{
		Cylinder<double> cyl;
		double initx = 0.0;
		double inity = 5.0;
		double initz = 0.0;
		double sourx, soury, sourz;
		double initdetcntx = 0;
		double initdetcnty = -5.0;
		double initdetcntz = 0.0;
		double detcntx, detcnty, detcntz;

		double initux = 1;	double inituy = 0;	double inituz = 0;
		double ux, uy, uz;
		double initvx = 0;	double initvy = 0;	double initvz = 1;
		double vx, vy, vz;
		const double detustp = 0.01;
		const double detvstp = 0.01;
		const int DNi = 512;
		const int DNj = 512;
		const int angN = 360;
		const double angStp = 3.14159265358979323846264 * 2.0 / angN;

		double* proj;
		checkCudaErrors(cudaMalloc((void**) &proj, sizeof(double) * DNi * DNj * angN));
		checkCudaErrors(cudaMemset(proj, 0, sizeof(double) * DNi * DNj * angN));

		double* A;
		checkCudaErrors(cudaMalloc((void**) &A, sizeof(double) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(A, cyl.A, sizeof(double) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		double* R;
		checkCudaErrors(cudaMalloc((void**) &R, sizeof(double) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(R, cyl.R, sizeof(double) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		double* posx;
		checkCudaErrors(cudaMalloc((void**) &posx, sizeof(double) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(posx, cyl.posx, sizeof(double) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		double* posy;
		checkCudaErrors(cudaMalloc((void**) &posy, sizeof(double) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(posy, cyl.posy, sizeof(double) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		double height = cyl.height;

		dim3 blk(16, 16);
		dim3 gid((DNi + blk.x - 1) / blk.x,
			(DNj + blk.y - 1) / blk.y);
		double cosT, sinT;
		for (unsigned int i = 0; i != angN; ++i)
		{
			cosT = cos(i * angStp);
			sinT = sin(i * angStp);

			sourx = initx * cosT - inity * sinT;
			soury = initx * sinT + inity * cosT;
			sourz = initz;
			detcntx = initdetcntx * cosT - initdetcnty * sinT;
			detcnty = initdetcntx * sinT + initdetcnty * cosT;
			detcntz = initdetcntz;

			ux = initux * cosT - inituy * sinT;
			uy = initux * sinT + inituy * cosT;
			uz = inituz;

			vx = initvx * cosT - initvy * sinT;
			vy = initvx * sinT + initvy * cosT;
			vz = initvz;


			projCylinder<double>(proj + i * DNi * DNj, sourx, soury, sourz,
				detcntx, detcnty, detcntz, ux, uy, uz, vx, vy, vz, detustp, detvstp,
				DNi, DNj, A, R, height, posx, posy, blk, gid);
			/*projPhantom<double>(proj + i * DNi * DNj, sourx, soury, 0, detcntx, detcnty, 0,
			ux, uy, 0, 0, 0, 1, detustp, detvstp, DNi, DNj, attenCoefs, haxis,
			core, rots, blk, gid);*/
			std::cout << i << "\n";
		}

		double* hproj = new double[DNi * DNj * angN];
		checkCudaErrors(cudaMemcpy(hproj, proj, sizeof(double) * DNi * DNj * angN, cudaMemcpyDeviceToHost));

		std::ofstream fou("res.raw", std::ios::binary);
		fou.write((char*) hproj, sizeof(double) * DNi * DNj * angN);
		fou.close();
	}



	template<typename T>
	void spiralScanCylinder_template(
		const T initx,
		const T inity,
		const T initz,
		const T initdetcntx,
		const T initdetcnty,
		const T initdetcntz,
		const T pitch,
		const T detustp,
		const T detvstp,
		const int DNi,
		const int DNj,
		const int angN,
		const std::vector<T>& angSeries,
		const std::string& ProjFileName)
	{
		Cylinder<T> cyl;
		T sourx, soury, sourz;
		T detcntx, detcnty, detcntz;
		T initux = 1;	T inituy = 0;	T inituz = 0;
		T ux, uy, uz;
		T initvx = 0;	T initvy = 0;	T initvz = 1;
		T vx, vy, vz;

		T* proj;
		checkCudaErrors(cudaMalloc((void**) &proj, sizeof(T) * DNi * DNj * angN));
		checkCudaErrors(cudaMemset(proj, 0, sizeof(T) * DNi * DNj * angN));

		T* A;
		checkCudaErrors(cudaMalloc((void**) &A, sizeof(T) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(A, cyl.A, sizeof(T) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		T* R;
		checkCudaErrors(cudaMalloc((void**) &R, sizeof(T) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(R, cyl.R, sizeof(T) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		T* posx;
		checkCudaErrors(cudaMalloc((void**) &posx, sizeof(T) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(posx, cyl.posx, sizeof(T) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		T* posy;
		checkCudaErrors(cudaMalloc((void**) &posy, sizeof(T) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(posy, cyl.posy, sizeof(T) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		T height = cyl.height;

		dim3 blk(16, 16);
		dim3 gid((DNi + blk.x - 1) / blk.x,
			(DNj + blk.y - 1) / blk.y);
		T cosT, sinT;
		for (unsigned int i = 0; i != angN; ++i)
		{
			cosT = cos(angSeries[i]);
			sinT = sin(angSeries[i]);

			sourx = initx * cosT - inity * sinT;
			soury = initx * sinT + inity * cosT;
			sourz = initz + (angSeries[i] * pitch);
			detcntx = initdetcntx * cosT - initdetcnty * sinT;
			detcnty = initdetcntx * sinT + initdetcnty * cosT;
			detcntz = initdetcntz + (angSeries[i] * pitch);

			ux = initux * cosT - inituy * sinT;
			uy = initux * sinT + inituy * cosT;
			uz = inituz;

			vx = initvx * cosT - initvy * sinT;
			vy = initvx * sinT + initvy * cosT;
			vz = initvz;


			projCylinder<T>(proj + i * DNi * DNj, sourx, soury, sourz,
				detcntx, detcnty, detcntz, ux, uy, uz, vx, vy, vz, detustp, detvstp,
				DNi, DNj, A, R, height, posx, posy, blk, gid);
			std::cout << ".";
		}

		T* hproj = new T[DNi * DNj * angN];
		checkCudaErrors(cudaMemcpy(hproj, proj, sizeof(T) * DNi * DNj * angN, cudaMemcpyDeviceToHost));

		std::ofstream fou(ProjFileName.c_str(), std::ios::binary);
		fou.write((char*) hproj, sizeof(T) * DNi * DNj * angN);
		fou.close();
	}

	void spiralScanCylinder(
		const double initx,
		const double inity,
		const double initz,
		const double initdetcntx,
		const double initdetcnty,
		const double initdetcntz,
		const double pitch,
		const double detustp,
		const double detvstp,
		const int DNi,
		const int DNj,
		const int angN,
		const std::vector<double>& angSeries,
		const std::string& ProjFileName)
	{
		spiralScanCylinder_template<double>(initx, inity, initz, initdetcntx, initdetcnty, initdetcntz, pitch, detustp, detvstp, DNi, DNj, angN, angSeries, ProjFileName);
	}

	void spiralScanCylinder(
		const float initx,
		const float inity,
		const float initz,
		const float initdetcntx,
		const float initdetcnty,
		const float initdetcntz,
		const float pitch,
		const float detustp,
		const float detvstp,
		const int DNi,
		const int DNj,
		const int angN,
		const std::vector<float>& angSeries,
		const std::string& ProjFileName)
	{
		spiralScanCylinder_template<float>(initx, inity, initz, initdetcntx, initdetcnty, initdetcntz, pitch, detustp, detvstp, DNi, DNj, angN, angSeries, ProjFileName);
	}




	extern "C"
		void spiralbeamScanCylinder()
	{
		Cylinder<double> cyl;
		double initx = 0.0;
		double inity = 5.0;
		double initz = -1.0;
		double sourx, soury, sourz;
		double initdetcntx = 0;
		double initdetcnty = -5.0;
		double initdetcntz = -1.0;
		double detcntx, detcnty, detcntz;
		double pitch = 0.2;
		double initux = 1;	double inituy = 0;	double inituz = 0;
		double ux, uy, uz;
		double initvx = 0;	double initvy = 0;	double initvz = 1;
		double vx, vy, vz;
		const double detustp = 0.01;
		const double detvstp = 0.01;
		const int DNi = 512;
		const int DNj = 512;
		const int angN = 360;
		const double angStp = 3.14159265358979323846264 * 2.0 / angN;

		double* proj;
		checkCudaErrors(cudaMalloc((void**) &proj, sizeof(double) * DNi * DNj * angN));
		checkCudaErrors(cudaMemset(proj, 0, sizeof(double) * DNi * DNj * angN));

		double* A;
		checkCudaErrors(cudaMalloc((void**) &A, sizeof(double) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(A, cyl.A, sizeof(double) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		double* R;
		checkCudaErrors(cudaMalloc((void**) &R, sizeof(double) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(R, cyl.R, sizeof(double) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		double* posx;
		checkCudaErrors(cudaMalloc((void**) &posx, sizeof(double) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(posx, cyl.posx, sizeof(double) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		double* posy;
		checkCudaErrors(cudaMalloc((void**) &posy, sizeof(double) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(posy, cyl.posy, sizeof(double) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		double height = cyl.height;

		dim3 blk(16, 16);
		dim3 gid((DNi + blk.x - 1) / blk.x,
			(DNj + blk.y - 1) / blk.y);
		double cosT, sinT;
		for (unsigned int i = 0; i != angN; ++i)
		{
			cosT = cos(i * angStp);
			sinT = sin(i * angStp);

			sourx = initx * cosT - inity * sinT;
			soury = initx * sinT + inity * cosT;
			sourz = initz + (i * angStp * pitch);
			detcntx = initdetcntx * cosT - initdetcnty * sinT;
			detcnty = initdetcntx * sinT + initdetcnty * cosT;
			detcntz = initdetcntz + (i * angStp * pitch);

			ux = initux * cosT - inituy * sinT;
			uy = initux * sinT + inituy * cosT;
			uz = inituz;

			vx = initvx * cosT - initvy * sinT;
			vy = initvx * sinT + initvy * cosT;
			vz = initvz;


			projCylinder<double>(proj + i * DNi * DNj, sourx, soury, sourz,
				detcntx, detcnty, detcntz, ux, uy, uz, vx, vy, vz, detustp, detvstp,
				DNi, DNj, A, R, height, posx, posy, blk, gid);
			/*projPhantom<double>(proj + i * DNi * DNj, sourx, soury, 0, detcntx, detcnty, 0,
			ux, uy, 0, 0, 0, 1, detustp, detvstp, DNi, DNj, attenCoefs, haxis,
			core, rots, blk, gid);*/
			std::cout << i << "\n";
		}

		double* hproj = new double[DNi * DNj * angN];
		checkCudaErrors(cudaMemcpy(hproj, proj, sizeof(double) * DNi * DNj * angN, cudaMemcpyDeviceToHost));

		std::ofstream fou("res.raw", std::ios::binary);
		fou.write((char*) hproj, sizeof(double) * DNi * DNj * angN);
		fou.close();
	}









	template<typename T>
	struct CylinderV2
	{
		int A[NUM_OF_CYLINDER]; // Materials
		T R[NUM_OF_CYLINDER];
		T posx[NUM_OF_CYLINDER];
		T posy[NUM_OF_CYLINDER];
		T height;
		T* materials;

		CylinderV2()
		{
			materials = new T[7 * NUM_OF_ENERGYCHANNEL];
			height = 1;
			A[0] = 0;		A[1] = 1;		A[2] = 2;
			A[3] = 3;		A[4] = 4;		A[5] = 5;
			A[6] = 6;		A[7] = 1;		A[8] = 2;
			A[9] = 3;		A[10] = 4;		A[11] = 5;
			A[12] = 6;		A[13] = 5;		A[14] = 5;
			A[15] = 5;		A[16] = 5;

			R[0] = 0.45 * 2; posx[0] = 0; posy[0] = 0;
			R[1] = 0.075 * 2; posx[1] = 0.55; posy[1] = 0;
			R[2] = 0.075 * 2; posx[2] = 0.275; posy[2] = -0.4763;
			R[3] = 0.075 * 2; posx[3] = -0.275; posy[3] = -0.4763;
			R[4] = 0.075 * 2; posx[4] = -0.55; posy[4] = 0;
			R[5] = 0.075 * 2; posx[5] = -0.275; posy[5] = 0.4763;
			R[6] = 0.075 * 2; posx[6] = 0.275; posy[6] = 0.4763;
			R[7] = 0.04 * 2; posx[7] = 0.275; posy[7] = 0;
			R[8] = 0.035 * 2; posx[8] = 0.1375; posy[8] = -0.23815;
			R[9] = 0.03 * 2; posx[9] = -0.1375; posy[9] = -0.23815;
			R[10] = 0.025 * 2; posx[10] = -0.275; posy[10] = 0;
			R[11] = 0.02 * 2; posx[11] = -0.1375; posy[11] = 0.23815;
			R[12] = 0.015 * 2; posx[12] = 0.1375; posy[12] = 0.23815;
			R[13] = 0.01 * 2; posx[13] = 0; posy[13] = -0.1;
			R[14] = 0.0075 * 2; posx[14] = 0.1; posy[14] = 0;
			R[15] = 0.005 * 2; posx[15] = 0; posy[15] = 0.1;
			R[16] = 0.0025 * 2; posx[16] = -0.1; posy[16] = 0;

			//Read attenuations 
			std::ifstream fin0("NewMaterial_NUM_OF_ENERGYCHANNELenergies_0.raw", std::ios::binary);
			fin0.read((char*) (materials + 0 * NUM_OF_ENERGYCHANNEL), sizeof(T) * NUM_OF_ENERGYCHANNEL);
			fin0.close();
			std::ifstream fin1("NewMaterial_NUM_OF_ENERGYCHANNELenergies_1.raw", std::ios::binary);
			fin1.read((char*) (materials + 1 * NUM_OF_ENERGYCHANNEL), sizeof(T) * NUM_OF_ENERGYCHANNEL);
			fin1.close();
			std::ifstream fin2("NewMaterial_NUM_OF_ENERGYCHANNELenergies_2.raw", std::ios::binary);
			fin2.read((char*) (materials + 2 * NUM_OF_ENERGYCHANNEL), sizeof(T) * NUM_OF_ENERGYCHANNEL);
			fin2.close();
			std::ifstream fin3("NewMaterial_NUM_OF_ENERGYCHANNELenergies_3.raw", std::ios::binary);
			fin3.read((char*) (materials + 3 * NUM_OF_ENERGYCHANNEL), sizeof(T) * NUM_OF_ENERGYCHANNEL);
			fin3.close();
			std::ifstream fin4("NewMaterial_NUM_OF_ENERGYCHANNELenergies_4.raw", std::ios::binary);
			fin4.read((char*) (materials + 4 * NUM_OF_ENERGYCHANNEL), sizeof(T) * NUM_OF_ENERGYCHANNEL);
			fin4.close();
			std::ifstream fin5("NewMaterial_NUM_OF_ENERGYCHANNELenergies_5.raw", std::ios::binary);
			fin5.read((char*) (materials + 5 * NUM_OF_ENERGYCHANNEL), sizeof(T) * NUM_OF_ENERGYCHANNEL);
			fin5.close();
			std::ifstream fin6("NewMaterial_NUM_OF_ENERGYCHANNELenergies_6.raw", std::ios::binary);
			fin6.read((char*) (materials + 6 * NUM_OF_ENERGYCHANNEL), sizeof(T) * NUM_OF_ENERGYCHANNEL);
			fin6.close();
		}
		~CylinderV2()
		{
			delete [] materials;
		}
	};





	template<typename T>
	__global__ void projCylinderMultiEnergies_ker(T* proj, const T sourx, const T soury, const T sourz,
		const T detcntx, const T detcnty, const T detcntz, const T ux, const T uy, const T uz,
		const T vx, const T vy, const T vz, const T detustp, const T detvstp, const int DNi, const int DNj,
		const int* A, const T* materials, const int channelsNum, const T* R, const T height, const T* posx, const T* posy)
	{
		const int i = threadIdx.x + blockIdx.x * blockDim.x;
		const int j = threadIdx.y + blockIdx.y * blockDim.y;
		const int energIdx = threadIdx.z + blockIdx.z * blockDim.z;

		//__shared__ T sA[NUM_OF_CYLINDER];
		__shared__ T sR[NUM_OF_CYLINDER];
		__shared__ T sposx[NUM_OF_CYLINDER];
		__shared__ T sposy[NUM_OF_CYLINDER];
		for (size_t i = 0; i < NUM_OF_CYLINDER; i++)
		{
			//sA[i] = A[i];
			sR[i] = R[i];
			sposx[i] = posx[i];
			sposy[i] = posy[i];
		}
		__syncthreads();
		if (i < DNi && j < DNj && energIdx < channelsNum)
		{
			const T uu = (detustp * (i - DNi * 0.5 + 0.5));
			const T vv = (detvstp * (j - DNj * 0.5 + 0.5));
			const T detpx = detcntx + ux * uu + vx * vv;
			const T detpy = detcnty + uy * uu + vy * vv;
			const T detpz = detcntz + uz * uu + vz * vv;
			T dirx = detpx - sourx;
			T diry = detpy - soury;
			T dirz = detpz - sourz;
			const T ll = sqrt(dirx * dirx + diry * diry + dirz * dirz);
			dirx /= ll;
			diry /= ll;
			dirz /= ll;
			T summ; //Changes here

			T intersectlength = 0;
			const T botz = -height * 0.5;
			const T topz = height * 0.5;
			bool intersect = intersectLengthWithCyl<T>(sourx, soury, sourz, dirx, diry, dirz,
				sposx[0], sposy[0], sR[0], botz, topz, intersectlength);

			if (!intersect)
			{
				proj[(energIdx * DNj + j) * DNi + i] = 0;
				return;
			}
			else
			{
				summ = intersectlength * materials[channelsNum * A[0] + energIdx];
				for (size_t ii = 1; ii < NUM_OF_CYLINDER; ii++)
				{
					if (intersectLengthWithCyl<T>(sourx, soury, sourz, dirx, diry, dirz,
						sposx[ii], sposy[ii], sR[ii], botz, topz, intersectlength))
					{
						summ += intersectlength * (materials[channelsNum * A[ii] + energIdx] - materials[channelsNum * A[0] + energIdx]);
					}
				}
				proj[(energIdx * DNj + j) * DNi + i] = summ;
			}
		}
	}

	template<typename T>
	void projCylinderMultiEnergies(T* proj, const T sourx, const T soury, const T sourz,
		const T detcntx, const T detcnty, const T detcntz, const T ux, const T uy, const T uz,
		const T vx, const T vy, const T vz, const T detustp, const T detvstp, const int DNi, const int DNj,
		const int* A, const T* materials, const int channelsNum, const T* R, const T height, const T* posx, const T* posy, const dim3 blk, const dim3 gid)
	{
		projCylinderMultiEnergies_ker<T> << <gid, blk >> >(proj, sourx, soury, sourz, detcntx, detcnty, detcntz, ux, uy, uz, vx, vy, vz, detustp, detvstp, DNi, DNj,
			A, materials, channelsNum, R, height, posx, posy);
	}


	template<typename T>
	__global__ void projCylinderMultiEnergies_ker(T* proj, const T sourx, const T soury, const T sourz,
		const T detcntx, const T detcnty, const T detcntz, const T ux, const T uy, const T uz,
		const T vx, const T vy, const T vz, const T detustp, const T detvstp, const int DNi, const int DNj,
		const T* photonNum,
		const int* A, const T* materials, const int channelsNum, const T* R, const T height, const T* posx, const T* posy)
	{
		const int i = threadIdx.x + blockIdx.x * blockDim.x;
		const int j = threadIdx.y + blockIdx.y * blockDim.y;
		const int energIdx = threadIdx.z + blockIdx.z * blockDim.z;

		//__shared__ T sA[NUM_OF_CYLINDER];
		__shared__ T sR[NUM_OF_CYLINDER];
		__shared__ T sposx[NUM_OF_CYLINDER];
		__shared__ T sposy[NUM_OF_CYLINDER];
		for (size_t i = 0; i < NUM_OF_CYLINDER; i++)
		{
//			sA[i] = A[i];
			sR[i] = R[i];
			sposx[i] = posx[i];
			sposy[i] = posy[i];
		}
		__syncthreads();
		if (i < DNi && j < DNj && energIdx < channelsNum)
		{
			const T uu = (detustp * (i - DNi * 0.5 + 0.5));
			const T vv = (detvstp * (j - DNj * 0.5 + 0.5));
			const T detpx = detcntx + ux * uu + vx * vv;
			const T detpy = detcnty + uy * uu + vy * vv;
			const T detpz = detcntz + uz * uu + vz * vv;
			T dirx = detpx - sourx;
			T diry = detpy - soury;
			T dirz = detpz - sourz;
			const T ll = sqrt(dirx * dirx + diry * diry + dirz * dirz);
			dirx /= ll;
			diry /= ll;
			dirz /= ll;
			T summ; //Changes here

			T intersectlength = 0;
			const T botz = -height * 0.5;
			const T topz = height * 0.5;
			bool intersect = intersectLengthWithCyl<T>(sourx, soury, sourz, dirx, diry, dirz,
				sposx[0], sposy[0], sR[0], botz, topz, intersectlength);

			if (!intersect)
			{
				proj[(energIdx * DNj + j) * DNi + i] = photonNum[energIdx];
				return;
			}
			else
			{
				summ = intersectlength * materials[channelsNum * A[0] + energIdx];
				for (size_t ii = 1; ii < NUM_OF_CYLINDER; ii++)
				{
					if (intersectLengthWithCyl<T>(sourx, soury, sourz, dirx, diry, dirz,
						sposx[ii], sposy[ii], sR[ii], botz, topz, intersectlength))
					{
						summ += intersectlength * (materials[channelsNum * A[ii] + energIdx] - materials[channelsNum * A[0] + energIdx]);
					}
				}
				proj[(energIdx * DNj + j) * DNi + i] = photonNum[energIdx] * exp(-summ);
			}
		}
	}

	template<typename T>
	void projCylinderMultiEnergies(T* proj, const T sourx, const T soury, const T sourz,
		const T detcntx, const T detcnty, const T detcntz, const T ux, const T uy, const T uz,
		const T vx, const T vy, const T vz, const T detustp, const T detvstp, const int DNi, const int DNj, const T* photonNum,
		const int* A, const T* materials, const int channelsNum, const T* R, const T height, const T* posx, const T* posy, const dim3 blk, const dim3 gid)
	{
		projCylinderMultiEnergies_ker<T> << <gid, blk >> >(proj, sourx, soury, sourz, detcntx, detcnty, detcntz, ux, uy, uz, vx, vy, vz, detustp, detvstp, DNi, DNj, photonNum,
			A, materials, channelsNum, R, height, posx, posy);
	}

	template<typename T>
	struct divideConstValue_functor
	{
		T value;
		divideConstValue_functor(const T& v) :value(v){}
		__host__ __device__ T operator()(const T& inpu)
		{
			if (value != 0)
			{
				return inpu / value;
			}
			else
			{
				return 0;
			}
		}
	};







	// \B4\D3NUM_OF_ENERGYCHANNEL\B8\F6ͨ\B5\C0\BAϳ\C9һ\B8\F6ͨ\B5\C0\A3\AC\B2\A2\CA\E4\B3\F68\B8\F6ͨ\B5\C0\B5\C4ֵ;
	// already done
	void conebeamScanCylinderMultiEnergiesV4(
		cublasHandle_t handle, const double initx, const double inity, const double initz,
		const double initdetcntx, const double initdetcnty, const double initdetcntz,
		const double detustp, const double detvstp,
		const int DNi, const int DNj, const int angN, const double angStp)
	{
		const double ONE = 1.0;
		const double ZERO = 0.0;
		cublasCreate(&handle);
		CylinderV2<double> cyl;
		double sourx, soury, sourz;
		double detcntx, detcnty, detcntz;

		double initux = 1;	double inituy = 0;	double inituz = 0;
		double ux, uy, uz;
		double initvx = 0;	double initvy = 0;	double initvz = 1;
		double vx, vy, vz;
		const int EnergiesChannels = NUM_OF_ENERGYCHANNEL; // We totally define NUM_OF_ENERGYCHANNEL channels;
		const int matrixRowNum = DNi * DNj;
		const int matrixColNum = EnergiesChannels;


		//Read the spectrum file.
		double* spRatio = new double[EnergiesChannels];
		std::ifstream fin("spectrum.raw", std::ios::binary);
		if (!fin.is_open())
		{
			std::cout << "Cannot open the spectrum data\n";
			std::cout << "We will use uniform distribution\n";
			for (size_t i = 0; i < EnergiesChannels; i++)
			{
				spRatio[i] = 1.0 / static_cast<double>(EnergiesChannels);
			}
		}
		else
		{
			fin.read((char*) spRatio, sizeof(double) * EnergiesChannels);
		}
		double* channelRatio[8];
		for (size_t i = 0; i < 8; i++)
		{
			channelRatio[i] = new double[EnergiesChannels];
			for (size_t j = 0; j < EnergiesChannels; j++)
			{
				channelRatio[i][j] = 0;
			}
		}
		std::copy(spRatio, spRatio + 110, channelRatio[0]);
		std::copy(spRatio + 110, spRatio + 210, channelRatio[1] + 110);
		std::copy(spRatio + 210, spRatio + 310, channelRatio[2] + 210);
		std::copy(spRatio + 310, spRatio + 410, channelRatio[3] + 310);
		std::copy(spRatio + 410, spRatio + 510, channelRatio[4] + 410);
		std::copy(spRatio + 510, spRatio + 610, channelRatio[5] + 510);
		std::copy(spRatio + 610, spRatio + 710, channelRatio[6] + 610);
		std::copy(spRatio + 710, spRatio + NUM_OF_ENERGYCHANNEL, channelRatio[7] + 710);


		const double sumSpectrum = thrust::reduce(spRatio, spRatio + EnergiesChannels);
		double sumChannelSpectrum[8] = { 0 };
		sumChannelSpectrum[0] = thrust::reduce(channelRatio[0], channelRatio[0] + EnergiesChannels);
		sumChannelSpectrum[1] = thrust::reduce(channelRatio[1], channelRatio[1] + EnergiesChannels);
		sumChannelSpectrum[2] = thrust::reduce(channelRatio[2], channelRatio[2] + EnergiesChannels);
		sumChannelSpectrum[3] = thrust::reduce(channelRatio[3], channelRatio[3] + EnergiesChannels);
		sumChannelSpectrum[4] = thrust::reduce(channelRatio[4], channelRatio[4] + EnergiesChannels);
		sumChannelSpectrum[5] = thrust::reduce(channelRatio[5], channelRatio[5] + EnergiesChannels);
		sumChannelSpectrum[6] = thrust::reduce(channelRatio[6], channelRatio[6] + EnergiesChannels);
		sumChannelSpectrum[7] = thrust::reduce(channelRatio[7], channelRatio[7] + EnergiesChannels);

		double* dspRatio;
		checkCudaErrors(cudaMalloc((void**) &dspRatio, sizeof(double) * EnergiesChannels));
		checkCudaErrors(cudaMemcpy(dspRatio, spRatio, sizeof(double) *EnergiesChannels, cudaMemcpyHostToDevice));
		double* multiSpectrumProj = new double[matrixRowNum];
		double* dmultiSpectrumProj;
		checkCudaErrors(cudaMalloc((void**) &dmultiSpectrumProj, sizeof(double) * matrixRowNum));
		double* channelProj[8];
		double* dchannelProj[8];
		double* dchannelRatio[8];
		for (size_t i = 0; i < 8; i++)
		{
			channelProj[i] = new double[matrixRowNum];
			checkCudaErrors(cudaMalloc((void**) &dchannelProj[i], sizeof(double) * matrixRowNum));

			checkCudaErrors(cudaMalloc((void**) &dchannelRatio[i], sizeof(double) * EnergiesChannels));
			checkCudaErrors(cudaMemcpy(dchannelRatio[i], channelRatio[i], sizeof(double) * EnergiesChannels, cudaMemcpyHostToDevice));
		}


		double* proj; //All energies channels
		checkCudaErrors(cudaMalloc((void**) &proj, sizeof(double) * DNi * DNj * EnergiesChannels));
		checkCudaErrors(cudaMemset(proj, 0, sizeof(double) * DNi * DNj * EnergiesChannels));

		int* A; //Different materials
		checkCudaErrors(cudaMalloc((void**) &A, sizeof(int) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(A, cyl.A, sizeof(int) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		double* R;
		checkCudaErrors(cudaMalloc((void**) &R, sizeof(double) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(R, cyl.R, sizeof(double) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		double* posx;
		checkCudaErrors(cudaMalloc((void**) &posx, sizeof(double) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(posx, cyl.posx, sizeof(double) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		double* posy;
		checkCudaErrors(cudaMalloc((void**) &posy, sizeof(double) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(posy, cyl.posy, sizeof(double) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		double* materials;
		checkCudaErrors(cudaMalloc((void**) &materials, sizeof(double) * 7 * NUM_OF_ENERGYCHANNEL));
		checkCudaErrors(cudaMemcpy(materials, cyl.materials, sizeof(double) * 7 * NUM_OF_ENERGYCHANNEL, cudaMemcpyHostToDevice));

		double height = cyl.height;

		dim3 blk(8, 8, 4);
		dim3 gid((DNi + blk.x - 1) / blk.x,
			(DNj + blk.y - 1) / blk.y,
			(EnergiesChannels + blk.z - 1) / blk.z);
		double cosT, sinT;
		double* hostproj = new double[DNi * DNj * EnergiesChannels];
		for (unsigned int i = 0; i != angN; ++i)
		{
			cosT = cos(i * angStp);
			sinT = sin(i * angStp);

			sourx = initx * cosT - inity * sinT;
			soury = initx * sinT + inity * cosT;
			sourz = initz;
			detcntx = initdetcntx * cosT - initdetcnty * sinT;
			detcnty = initdetcntx * sinT + initdetcnty * cosT;
			detcntz = initdetcntz;

			ux = initux * cosT - inituy * sinT;
			uy = initux * sinT + inituy * cosT;
			uz = inituz;

			vx = initvx * cosT - initvy * sinT;
			vy = initvx * sinT + initvy * cosT;
			vz = initvz;

			projCylinderMultiEnergies<double>(proj, sourx, soury, sourz,
				detcntx, detcnty, detcntz, ux, uy, uz, vx, vy, vz,
				detustp, detvstp, DNi, DNj, A, materials, EnergiesChannels, R, height, posx, posy, blk, gid);
			//cudaMemcpy(hostproj, proj, sizeof(double) * DNi * DNj * EnergiesChannels, cudaMemcpyDeviceToHost);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dspRatio, 1, &ZERO, dmultiSpectrumProj, 1);

			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[0], 1, &ZERO, dchannelProj[0], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[1], 1, &ZERO, dchannelProj[1], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[2], 1, &ZERO, dchannelProj[2], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[3], 1, &ZERO, dchannelProj[3], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[4], 1, &ZERO, dchannelProj[4], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[5], 1, &ZERO, dchannelProj[5], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[6], 1, &ZERO, dchannelProj[6], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[7], 1, &ZERO, dchannelProj[7], 1);
			checkCudaErrors(cudaMemcpy(multiSpectrumProj, dmultiSpectrumProj, sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[0], dchannelProj[0], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[1], dchannelProj[1], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[2], dchannelProj[2], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[3], dchannelProj[3], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[4], dchannelProj[4], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[5], dchannelProj[5], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[6], dchannelProj[6], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[7], dchannelProj[7], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));


			thrust::transform(multiSpectrumProj, multiSpectrumProj + matrixRowNum, multiSpectrumProj, divideConstValue_functor<double>(sumSpectrum));
			thrust::transform(channelProj[0], channelProj[0] + matrixRowNum, channelProj[0], divideConstValue_functor<double>(sumChannelSpectrum[0]));
			thrust::transform(channelProj[1], channelProj[1] + matrixRowNum, channelProj[1], divideConstValue_functor<double>(sumChannelSpectrum[1]));
			thrust::transform(channelProj[2], channelProj[2] + matrixRowNum, channelProj[2], divideConstValue_functor<double>(sumChannelSpectrum[2]));
			thrust::transform(channelProj[3], channelProj[3] + matrixRowNum, channelProj[3], divideConstValue_functor<double>(sumChannelSpectrum[3]));
			thrust::transform(channelProj[4], channelProj[4] + matrixRowNum, channelProj[4], divideConstValue_functor<double>(sumChannelSpectrum[4]));
			thrust::transform(channelProj[5], channelProj[5] + matrixRowNum, channelProj[5], divideConstValue_functor<double>(sumChannelSpectrum[5]));
			thrust::transform(channelProj[6], channelProj[6] + matrixRowNum, channelProj[6], divideConstValue_functor<double>(sumChannelSpectrum[6]));
			thrust::transform(channelProj[7], channelProj[7] + matrixRowNum, channelProj[7], divideConstValue_functor<double>(sumChannelSpectrum[7]));

			std::stringstream ss;
			ss << i;
			/*std::string ProjName = "G:\\MultiEnergies\\ProjAngle" + ss.str() + ".raw";
			std::ofstream fou(ProjName.c_str(), std::ios::binary);
			fou.write((char*) hostproj, sizeof(double) * DNi * DNj * EnergiesChannels);
			fou.close();*/
			std::string chProjName0 = "Channel0Proj" + ss.str() + ".raw";
			std::string chProjName1 = "Channel1Proj" + ss.str() + ".raw";
			std::string chProjName2 = "Channel2Proj" + ss.str() + ".raw";
			std::string chProjName3 = "Channel3Proj" + ss.str() + ".raw";
			std::string chProjName4 = "Channel4Proj" + ss.str() + ".raw";
			std::string chProjName5 = "Channel5Proj" + ss.str() + ".raw";
			std::string chProjName6 = "Channel6Proj" + ss.str() + ".raw";
			std::string chProjName7 = "Channel7Proj" + ss.str() + ".raw";

			std::string chProjName8 = "allChannelsProj" + ss.str() + ".raw";

			std::ofstream fou0(chProjName0.c_str(), std::ios::binary);
			fou0.write((char*) channelProj[0], sizeof(double) * matrixRowNum);
			fou0.close();

			std::ofstream fou1(chProjName1.c_str(), std::ios::binary);
			fou1.write((char*) channelProj[1], sizeof(double) * matrixRowNum);
			fou1.close();

			std::ofstream fou2(chProjName2.c_str(), std::ios::binary);
			fou2.write((char*) channelProj[2], sizeof(double) * matrixRowNum);
			fou2.close();

			std::ofstream fou3(chProjName3.c_str(), std::ios::binary);
			fou3.write((char*) channelProj[3], sizeof(double) * matrixRowNum);
			fou3.close();

			std::ofstream fou4(chProjName4.c_str(), std::ios::binary);
			fou4.write((char*) channelProj[4], sizeof(double) * matrixRowNum);
			fou4.close();

			std::ofstream fou5(chProjName5.c_str(), std::ios::binary);
			fou5.write((char*) channelProj[5], sizeof(double) * matrixRowNum);
			fou5.close();

			std::ofstream fou6(chProjName6.c_str(), std::ios::binary);
			fou6.write((char*) channelProj[6], sizeof(double) * matrixRowNum);
			fou6.close();

			std::ofstream fou7(chProjName7.c_str(), std::ios::binary);
			fou7.write((char*) channelProj[7], sizeof(double) * matrixRowNum);
			fou7.close();

			std::ofstream fou8(chProjName8.c_str(), std::ios::binary);
			fou8.write((char*) multiSpectrumProj, sizeof(double) * matrixRowNum);
			fou8.close();

			std::cout << i << "\n";
		}
	}




	template<typename T>
	struct HeadPhantom3D
	{
		int M[10]; //index of Different Materials for these ellipsoids
		T3<T> axis[10]; // axis length of these ellipsoids
		T3<T> core[10]; // core position of these ellipsoids
		T rots[10]; //rotation angle of these ellipsoids

		T* materials;
		HeadPhantom3D()
		{
			M[0] = 0;
			M[1] = 1;
			M[2] = 2;
			M[3] = 2;
			M[4] = 3;
			M[5] = 4;
			M[6] = 5;
			M[7] = 6;
			M[8] = 6;
			M[9] = 6;
			materials = new T[7 * NUM_OF_ENERGYCHANNEL];

			//Read attenuations 
			std::ifstream fin0("biologicalmaterial_NUM_OF_ENERGYCHANNELenergies_0.raw", std::ios::binary);
			std::ifstream fin1("biologicalmaterial_NUM_OF_ENERGYCHANNELenergies_1.raw", std::ios::binary);
			std::ifstream fin2("biologicalmaterial_NUM_OF_ENERGYCHANNELenergies_2.raw", std::ios::binary);
			std::ifstream fin3("biologicalmaterial_NUM_OF_ENERGYCHANNELenergies_3.raw", std::ios::binary);
			std::ifstream fin4("biologicalmaterial_NUM_OF_ENERGYCHANNELenergies_4.raw", std::ios::binary);
			std::ifstream fin5("biologicalmaterial_NUM_OF_ENERGYCHANNELenergies_5.raw", std::ios::binary);
			std::ifstream fin6("biologicalmaterial_NUM_OF_ENERGYCHANNELenergies_6.raw", std::ios::binary);

			if (!(fin0.is_open() && fin1.is_open() && fin2.is_open() && fin3.is_open() &&
				fin4.is_open() && fin5.is_open() && fin6.is_open()))
			{
				std::cout << "Cannot read material information\n";
				exit(-1);
			}

			fin0.read((char*) (materials + 0 * NUM_OF_ENERGYCHANNEL), sizeof(T) * NUM_OF_ENERGYCHANNEL);
			fin0.close();

			fin1.read((char*) (materials + 1 * NUM_OF_ENERGYCHANNEL), sizeof(T) * NUM_OF_ENERGYCHANNEL);
			fin1.close();

			fin2.read((char*) (materials + 2 * NUM_OF_ENERGYCHANNEL), sizeof(T) * NUM_OF_ENERGYCHANNEL);
			fin2.close();

			fin3.read((char*) (materials + 3 * NUM_OF_ENERGYCHANNEL), sizeof(T) * NUM_OF_ENERGYCHANNEL);
			fin3.close();

			fin4.read((char*) (materials + 4 * NUM_OF_ENERGYCHANNEL), sizeof(T) * NUM_OF_ENERGYCHANNEL);
			fin4.close();

			fin5.read((char*) (materials + 5 * NUM_OF_ENERGYCHANNEL), sizeof(T) * NUM_OF_ENERGYCHANNEL);
			fin5.close();

			fin6.read((char*) (materials + 6 * NUM_OF_ENERGYCHANNEL), sizeof(T) * NUM_OF_ENERGYCHANNEL);
			fin6.close();

			axis[0].x = 0.69; axis[0].y = 0.92; axis[0].z = 0.81;
			axis[1].x = 0.6624; axis[1].y = 0.874; axis[1].z = 0.780;
			axis[2].x = 0.11; axis[2].y = 0.31; axis[2].z = 0.22;
			axis[3].x = 0.16; axis[3].y = 0.41; axis[3].z = 0.28;
			axis[4].x = 0.21; axis[4].y = 0.25; axis[4].z = 0.41;
			axis[5].x = 0.046; axis[5].y = 0.046; axis[5].z = 0.05;
			axis[6].x = 0.046; axis[6].y = 0.046; axis[6].z = 0.046;
			axis[7].x = 0.046; axis[7].y = 0.023; axis[7].z = 0.05;
			axis[8].x = 0.023; axis[8].y = 0.023; axis[8].z = 0.02;
			axis[9].x = 0.023; axis[9].y = 0.046; axis[9].z = 0.02;

			core[0].x = 0; core[0].y = 0; core[0].z = 0;
			core[1].x = 0; core[1].y = -0.0184; core[1].z = 0;
			core[2].x = 0.22; core[2].y = 0; core[2].z = 0;
			core[3].x = -0.22; core[3].y = 0; core[3].z = 0;
			core[4].x = 0; core[4].y = 0.35; core[4].z = -0.15;
			core[5].x = 0; core[5].y = 0.1; core[5].z = 0.25;
			core[6].x = 0; core[6].y = -0.1; core[6].z = 0.25;
			core[7].x = -0.08; core[7].y = -0.605; core[7].z = 0;
			core[8].x = 0; core[8].y = -0.606; core[8].z = 0;
			core[9].x = 0.06; core[9].y = -0.605; core[9].z = 0;

			rots[0] = 0; rots[1] = 0; rots[2] = -18.0 / 180.0 * PI;
			rots[3] = 18.0 / 180.0 * PI; rots[4] = 0; rots[5] = 0;
			rots[6] = 0; rots[7] = 0; rots[8] = 0;
			rots[9] = 0;
		}
		~HeadPhantom3D()
		{
			delete [] materials;
		}

	};



	template<typename T>
	__global__ void projHeadPhantomMultiEnergies_ker(T* proj, const T sourx, const T soury, const T sourz,
		const T detcntx, const T detcnty, const T detcntz, const T ux, const T uy, const T uz,
		const T vx, const T vy, const T vz, const T detustp, const T detvstp, const int DNi, const int DNj,
		const int* M, const T* materials, const int channelsNum, const T* axisX, const T* axisY, const T* axisZ,
		const T* coreX, const T* coreY, const T* coreZ, const T* rots)
	{
		const int i = threadIdx.x + blockIdx.x * blockDim.x;
		const int j = threadIdx.y + blockIdx.y * blockDim.y;
		const int energIdx = threadIdx.z + blockIdx.z * blockDim.z;
		__shared__ T sax[10];
		__shared__ T say[10];
		__shared__ T saz[10];
		__shared__ T cox[10];
		__shared__ T coy[10];
		__shared__ T coz[10];
		//__shared__ T mat[10][100]; //It can be changed with configurations
		for (size_t i = 0; i < 10; i++)
		{
			sax[i] = axisX[i];
			say[i] = axisY[i];
			saz[i] = axisZ[i];
			cox[i] = coreX[i];
			coy[i] = coreY[i];
			coz[i] = coreZ[i];
			//mat[i][energIdx] = materials[channelsNum * M[i] + energIdx];
		}
		__syncthreads();

		if (i < DNi && j < DNj && energIdx < channelsNum)
		{
			const T uu = (detustp * (i - DNi * 0.5 + 0.5));
			const T vv = (detvstp * (j - DNj * 0.5 + 0.5));
			const T detpx = detcntx + ux * uu + vx * vv;
			const T detpy = detcnty + uy * uu + vy * vv;
			const T detpz = detcntz + uz * uu + vz * vv;
			T summ = 0;

			T sxsX = 0;		T spsX = 0;		T sxsY = 0;
			T spsY = 0;		T sxsZ = 0;		T spsZ = 0;
			T pcosT = 0;		T psinT = 0;
			T axisSquareX = 0;		T axisSquareY = 0;		T axisSquareZ = 0;
			T tempvar = 0;
			T dirx = 0;		T diry = 0;		T dirz = 0;
			T AA = 0;		T BB = 0;		T CC = 0;
			T leng;
			T mater0 = materials[channelsNum * M[0] + energIdx];
			T mater1 = materials[channelsNum * M[1] + energIdx];
			T mater2 = materials[channelsNum * M[2] + energIdx];
			T mater3 = materials[channelsNum * M[3] + energIdx];
			T mater4 = materials[channelsNum * M[4] + energIdx];
			T mater5 = materials[channelsNum * M[5] + energIdx];
			T mater6 = materials[channelsNum * M[6] + energIdx];
			T mater7 = materials[channelsNum * M[7] + energIdx];
			T mater8 = materials[channelsNum * M[8] + energIdx];
			T mater9 = materials[channelsNum * M[9] + energIdx];


			sxsX = sourx - cox[0];
			spsX = detpx - cox[0];
			sxsY = soury - coy[0];
			spsY = detpy - coy[0];
			sxsZ = sourz - coz[0];
			spsZ = detpz - coz[0];
			pcosT = cos(rots[0]);
			psinT = sin(rots[0]);
			axisSquareX = sax[0] * sax[0];
			axisSquareY = say[0] * say[0];
			axisSquareZ = saz[0] * saz[0];
			tempvar = sxsX;
			sxsX = tempvar * pcosT + sxsY * psinT;
			sxsY = -tempvar * psinT + sxsY * pcosT;
			tempvar = spsX;
			spsX = tempvar * pcosT + spsY * psinT;
			spsY = -tempvar * psinT + spsY * pcosT;
			dirx = sxsX - spsX;
			diry = sxsY - spsY;
			dirz = sxsZ - spsZ;
			tempvar = sqrt(dirx * dirx + diry * diry + dirz * dirz);
			dirx /= tempvar; diry /= tempvar; dirz /= tempvar;

			AA = pow(dirx, 2) / axisSquareX + pow(diry, 2) / axisSquareY + pow(dirz, 2) / axisSquareZ;
			BB = dirx * spsX / axisSquareX + diry * spsY / axisSquareY + dirz * spsZ / axisSquareZ;
			CC = pow(spsX, 2) / axisSquareX + pow(spsY, 2) / axisSquareY + pow(spsZ, 2) / axisSquareZ - 1;
			tempvar = BB * BB - AA * CC;
			if (tempvar <= 0)
			{
				proj[(energIdx * DNj + j) * DNi + i] = 0;
				return;
			}
			leng = sqrt(tempvar) / AA;
			summ += leng * mater0;



			sxsX = sourx - cox[1];
			spsX = detpx - cox[1];
			sxsY = soury - coy[1];
			spsY = detpy - coy[1];
			sxsZ = sourz - coz[1];
			spsZ = detpz - coz[1];
			pcosT = cos(rots[1]);
			psinT = sin(rots[1]);
			axisSquareX = sax[1] * sax[1];
			axisSquareY = say[1] * say[1];
			axisSquareZ = saz[1] * saz[1];
			tempvar = sxsX;
			sxsX = tempvar * pcosT + sxsY * psinT;
			sxsY = -tempvar * psinT + sxsY * pcosT;
			tempvar = spsX;
			spsX = tempvar * pcosT + spsY * psinT;
			spsY = -tempvar * psinT + spsY * pcosT;
			dirx = sxsX - spsX;
			diry = sxsY - spsY;
			dirz = sxsZ - spsZ;
			tempvar = sqrt(dirx * dirx + diry * diry + dirz * dirz);
			dirx /= tempvar; diry /= tempvar; dirz /= tempvar;

			AA = pow(dirx, 2) / axisSquareX + pow(diry, 2) / axisSquareY + pow(dirz, 2) / axisSquareZ;
			BB = dirx * spsX / axisSquareX + diry * spsY / axisSquareY + dirz * spsZ / axisSquareZ;
			CC = pow(spsX, 2) / axisSquareX + pow(spsY, 2) / axisSquareY + pow(spsZ, 2) / axisSquareZ - 1;
			tempvar = BB * BB - AA * CC;
			if (tempvar <= 0)
			{
				proj[(energIdx * DNj + j) * DNi + i] = summ;
				return;
			}
			leng = sqrt(tempvar) / AA;
			summ += leng * (mater1 - mater0);

			sxsX = sourx - cox[2];
			spsX = detpx - cox[2];
			sxsY = soury - coy[2];
			spsY = detpy - coy[2];
			sxsZ = sourz - coz[2];
			spsZ = detpz - coz[2];
			pcosT = cos(rots[2]);
			psinT = sin(rots[2]);
			axisSquareX = sax[2] * sax[2];
			axisSquareY = say[2] * say[2];
			axisSquareZ = saz[2] * saz[2];
			tempvar = sxsX;
			sxsX = tempvar * pcosT + sxsY * psinT;
			sxsY = -tempvar * psinT + sxsY * pcosT;
			tempvar = spsX;
			spsX = tempvar * pcosT + spsY * psinT;
			spsY = -tempvar * psinT + spsY * pcosT;
			dirx = sxsX - spsX;
			diry = sxsY - spsY;
			dirz = sxsZ - spsZ;
			tempvar = sqrt(dirx * dirx + diry * diry + dirz * dirz);
			dirx /= tempvar; diry /= tempvar; dirz /= tempvar;

			AA = pow(dirx, 2) / axisSquareX + pow(diry, 2) / axisSquareY + pow(dirz, 2) / axisSquareZ;
			BB = dirx * spsX / axisSquareX + diry * spsY / axisSquareY + dirz * spsZ / axisSquareZ;
			CC = pow(spsX, 2) / axisSquareX + pow(spsY, 2) / axisSquareY + pow(spsZ, 2) / axisSquareZ - 1;
			tempvar = BB * BB - AA * CC;
			if (tempvar <= 0)
			{
				summ += 0;
			}
			else
			{
				leng = sqrt(tempvar) / AA;
				summ += leng * (mater2 - mater1);
			}

			sxsX = sourx - cox[3];
			spsX = detpx - cox[3];
			sxsY = soury - coy[3];
			spsY = detpy - coy[3];
			sxsZ = sourz - coz[3];
			spsZ = detpz - coz[3];
			pcosT = cos(rots[3]);
			psinT = sin(rots[3]);
			axisSquareX = sax[3] * sax[3];
			axisSquareY = say[3] * say[3];
			axisSquareZ = saz[3] * saz[3];
			tempvar = sxsX;
			sxsX = tempvar * pcosT + sxsY * psinT;
			sxsY = -tempvar * psinT + sxsY * pcosT;
			tempvar = spsX;
			spsX = tempvar * pcosT + spsY * psinT;
			spsY = -tempvar * psinT + spsY * pcosT;
			dirx = sxsX - spsX;
			diry = sxsY - spsY;
			dirz = sxsZ - spsZ;
			tempvar = sqrt(dirx * dirx + diry * diry + dirz * dirz);
			dirx /= tempvar; diry /= tempvar; dirz /= tempvar;

			AA = pow(dirx, 2) / axisSquareX + pow(diry, 2) / axisSquareY + pow(dirz, 2) / axisSquareZ;
			BB = dirx * spsX / axisSquareX + diry * spsY / axisSquareY + dirz * spsZ / axisSquareZ;
			CC = pow(spsX, 2) / axisSquareX + pow(spsY, 2) / axisSquareY + pow(spsZ, 2) / axisSquareZ - 1;
			tempvar = BB * BB - AA * CC;
			if (tempvar <= 0)
			{
				summ += 0;
			}
			else
			{
				leng = sqrt(tempvar) / AA;
				summ += leng * (mater3 - mater1);
			}


			sxsX = sourx - cox[4];
			spsX = detpx - cox[4];
			sxsY = soury - coy[4];
			spsY = detpy - coy[4];
			sxsZ = sourz - coz[4];
			spsZ = detpz - coz[4];
			pcosT = cos(rots[4]);
			psinT = sin(rots[4]);
			axisSquareX = sax[4] * sax[4];
			axisSquareY = say[4] * say[4];
			axisSquareZ = saz[4] * saz[4];
			tempvar = sxsX;
			sxsX = tempvar * pcosT + sxsY * psinT;
			sxsY = -tempvar * psinT + sxsY * pcosT;
			tempvar = spsX;
			spsX = tempvar * pcosT + spsY * psinT;
			spsY = -tempvar * psinT + spsY * pcosT;
			dirx = sxsX - spsX;
			diry = sxsY - spsY;
			dirz = sxsZ - spsZ;
			tempvar = sqrt(dirx * dirx + diry * diry + dirz * dirz);
			dirx /= tempvar; diry /= tempvar; dirz /= tempvar;

			AA = pow(dirx, 2) / axisSquareX + pow(diry, 2) / axisSquareY + pow(dirz, 2) / axisSquareZ;
			BB = dirx * spsX / axisSquareX + diry * spsY / axisSquareY + dirz * spsZ / axisSquareZ;
			CC = pow(spsX, 2) / axisSquareX + pow(spsY, 2) / axisSquareY + pow(spsZ, 2) / axisSquareZ - 1;
			tempvar = BB * BB - AA * CC;
			if (tempvar <= 0)
			{
				summ += 0;
			}
			else
			{
				leng = sqrt(tempvar) / AA;
				summ += leng * (mater4 - mater1);
			}


			sxsX = sourx - cox[5];
			spsX = detpx - cox[5];
			sxsY = soury - coy[5];
			spsY = detpy - coy[5];
			sxsZ = sourz - coz[5];
			spsZ = detpz - coz[5];
			pcosT = cos(rots[5]);
			psinT = sin(rots[5]);
			axisSquareX = sax[5] * sax[5];
			axisSquareY = say[5] * say[5];
			axisSquareZ = saz[5] * saz[5];
			tempvar = sxsX;
			sxsX = tempvar * pcosT + sxsY * psinT;
			sxsY = -tempvar * psinT + sxsY * pcosT;
			tempvar = spsX;
			spsX = tempvar * pcosT + spsY * psinT;
			spsY = -tempvar * psinT + spsY * pcosT;
			dirx = sxsX - spsX;
			diry = sxsY - spsY;
			dirz = sxsZ - spsZ;
			tempvar = sqrt(dirx * dirx + diry * diry + dirz * dirz);
			dirx /= tempvar; diry /= tempvar; dirz /= tempvar;

			AA = pow(dirx, 2) / axisSquareX + pow(diry, 2) / axisSquareY + pow(dirz, 2) / axisSquareZ;
			BB = dirx * spsX / axisSquareX + diry * spsY / axisSquareY + dirz * spsZ / axisSquareZ;
			CC = pow(spsX, 2) / axisSquareX + pow(spsY, 2) / axisSquareY + pow(spsZ, 2) / axisSquareZ - 1;
			tempvar = BB * BB - AA * CC;
			if (tempvar <= 0)
			{
				summ += 0;
			}
			else
			{
				leng = sqrt(tempvar) / AA;
				summ += leng * (mater5 - mater1);
			}


			sxsX = sourx - cox[6];
			spsX = detpx - cox[6];
			sxsY = soury - coy[6];
			spsY = detpy - coy[6];
			sxsZ = sourz - coz[6];
			spsZ = detpz - coz[6];
			pcosT = cos(rots[6]);
			psinT = sin(rots[6]);
			axisSquareX = sax[6] * sax[6];
			axisSquareY = say[6] * say[6];
			axisSquareZ = saz[6] * saz[6];
			tempvar = sxsX;
			sxsX = tempvar * pcosT + sxsY * psinT;
			sxsY = -tempvar * psinT + sxsY * pcosT;
			tempvar = spsX;
			spsX = tempvar * pcosT + spsY * psinT;
			spsY = -tempvar * psinT + spsY * pcosT;
			dirx = sxsX - spsX;
			diry = sxsY - spsY;
			dirz = sxsZ - spsZ;
			tempvar = sqrt(dirx * dirx + diry * diry + dirz * dirz);
			dirx /= tempvar; diry /= tempvar; dirz /= tempvar;

			AA = pow(dirx, 2) / axisSquareX + pow(diry, 2) / axisSquareY + pow(dirz, 2) / axisSquareZ;
			BB = dirx * spsX / axisSquareX + diry * spsY / axisSquareY + dirz * spsZ / axisSquareZ;
			CC = pow(spsX, 2) / axisSquareX + pow(spsY, 2) / axisSquareY + pow(spsZ, 2) / axisSquareZ - 1;
			tempvar = BB * BB - AA * CC;
			if (tempvar <= 0)
			{
				summ += 0;
			}
			else
			{
				leng = sqrt(tempvar) / AA;
				summ += leng * (mater6 - mater1);
			}

			sxsX = sourx - cox[7];
			spsX = detpx - cox[7];
			sxsY = soury - coy[7];
			spsY = detpy - coy[7];
			sxsZ = sourz - coz[7];
			spsZ = detpz - coz[7];
			pcosT = cos(rots[7]);
			psinT = sin(rots[7]);
			axisSquareX = sax[7] * sax[7];
			axisSquareY = say[7] * say[7];
			axisSquareZ = saz[7] * saz[7];
			tempvar = sxsX;
			sxsX = tempvar * pcosT + sxsY * psinT;
			sxsY = -tempvar * psinT + sxsY * pcosT;
			tempvar = spsX;
			spsX = tempvar * pcosT + spsY * psinT;
			spsY = -tempvar * psinT + spsY * pcosT;
			dirx = sxsX - spsX;
			diry = sxsY - spsY;
			dirz = sxsZ - spsZ;
			tempvar = sqrt(dirx * dirx + diry * diry + dirz * dirz);
			dirx /= tempvar; diry /= tempvar; dirz /= tempvar;

			AA = pow(dirx, 2) / axisSquareX + pow(diry, 2) / axisSquareY + pow(dirz, 2) / axisSquareZ;
			BB = dirx * spsX / axisSquareX + diry * spsY / axisSquareY + dirz * spsZ / axisSquareZ;
			CC = pow(spsX, 2) / axisSquareX + pow(spsY, 2) / axisSquareY + pow(spsZ, 2) / axisSquareZ - 1;
			tempvar = BB * BB - AA * CC;
			if (tempvar <= 0)
			{
				summ += 0;
			}
			else
			{
				leng = sqrt(tempvar) / AA;
				summ += leng * (mater7 - mater1);
			}




			sxsX = sourx - cox[8];
			spsX = detpx - cox[8];
			sxsY = soury - coy[8];
			spsY = detpy - coy[8];
			sxsZ = sourz - coz[8];
			spsZ = detpz - coz[8];
			pcosT = cos(rots[8]);
			psinT = sin(rots[8]);
			axisSquareX = sax[8] * sax[8];
			axisSquareY = say[8] * say[8];
			axisSquareZ = saz[8] * saz[8];
			tempvar = sxsX;
			sxsX = tempvar * pcosT + sxsY * psinT;
			sxsY = -tempvar * psinT + sxsY * pcosT;
			tempvar = spsX;
			spsX = tempvar * pcosT + spsY * psinT;
			spsY = -tempvar * psinT + spsY * pcosT;
			dirx = sxsX - spsX;
			diry = sxsY - spsY;
			dirz = sxsZ - spsZ;
			tempvar = sqrt(dirx * dirx + diry * diry + dirz * dirz);
			dirx /= tempvar; diry /= tempvar; dirz /= tempvar;

			AA = pow(dirx, 2) / axisSquareX + pow(diry, 2) / axisSquareY + pow(dirz, 2) / axisSquareZ;
			BB = dirx * spsX / axisSquareX + diry * spsY / axisSquareY + dirz * spsZ / axisSquareZ;
			CC = pow(spsX, 2) / axisSquareX + pow(spsY, 2) / axisSquareY + pow(spsZ, 2) / axisSquareZ - 1;
			tempvar = BB * BB - AA * CC;
			if (tempvar <= 0)
			{
				summ += 0;
			}
			else
			{
				leng = sqrt(tempvar) / AA;
				summ += leng * (mater8 - mater1);
			}


			sxsX = sourx - cox[9];
			spsX = detpx - cox[9];
			sxsY = soury - coy[9];
			spsY = detpy - coy[9];
			sxsZ = sourz - coz[9];
			spsZ = detpz - coz[9];
			pcosT = cos(rots[9]);
			psinT = sin(rots[9]);
			axisSquareX = sax[9] * sax[9];
			axisSquareY = say[9] * say[9];
			axisSquareZ = saz[9] * saz[9];
			tempvar = sxsX;
			sxsX = tempvar * pcosT + sxsY * psinT;
			sxsY = -tempvar * psinT + sxsY * pcosT;
			tempvar = spsX;
			spsX = tempvar * pcosT + spsY * psinT;
			spsY = -tempvar * psinT + spsY * pcosT;
			dirx = sxsX - spsX;
			diry = sxsY - spsY;
			dirz = sxsZ - spsZ;
			tempvar = sqrt(dirx * dirx + diry * diry + dirz * dirz);
			dirx /= tempvar; diry /= tempvar; dirz /= tempvar;

			AA = pow(dirx, 2) / axisSquareX + pow(diry, 2) / axisSquareY + pow(dirz, 2) / axisSquareZ;
			BB = dirx * spsX / axisSquareX + diry * spsY / axisSquareY + dirz * spsZ / axisSquareZ;
			CC = pow(spsX, 2) / axisSquareX + pow(spsY, 2) / axisSquareY + pow(spsZ, 2) / axisSquareZ - 1;
			tempvar = BB * BB - AA * CC;
			if (tempvar <= 0)
			{
				summ += 0;
			}
			else
			{
				leng = sqrt(tempvar) / AA;
				summ += leng * (mater9 - mater1);
			}

			proj[(energIdx * DNj + j) * DNi + i] = summ;
		}
	}

	template<typename T>
	void projHeadPhantomMultiEnergies(T* proj, const T sourx, const T soury, const T sourz,
		const T detcntx, const T detcnty, const T detcntz, const T ux, const T uy, const T uz,
		const T vx, const T vy, const T vz, const T detustp, const T detvstp, const int DNi, const int DNj,
		const int* M, const T* materials, const int channelsNum, const T* axisX, const T* axisY, const T* axisZ,
		const T* coreX, const T* coreY, const T* coreZ, const T* rots, const dim3 blk, const dim3 gid)
	{
		projHeadPhantomMultiEnergies_ker<T> << <gid, blk >> >(proj, sourx, soury, sourz, detcntx, detcnty, detcntz,
			ux, uy, uz, vx, vy, vz, detustp, detvstp, DNi, DNj, M, materials, channelsNum, axisX, axisY, axisZ,
			coreX, coreY, coreZ, rots);
	}


	// already done
	void conebeamScanHeadPhantomMultiEnergies(
		cublasHandle_t handle, const double initx, const double inity, const double initz,
		const double initdetcntx, const double initdetcnty, const double initdetcntz,
		const double detustp, const double detvstp,
		const int DNi, const int DNj, const int angN, const double angStp)
	{
		const double ONE = 1.0;
		const double ZERO = 0.0;
		cublasCreate(&handle);
		HeadPhantom3D<double> cyl;
		//CylinderV2<double> cyl;
		double sourx, soury, sourz;
		double detcntx, detcnty, detcntz;

		double initux = 1;	double inituy = 0;	double inituz = 0;
		double ux, uy, uz;
		double initvx = 0;	double initvy = 0;	double initvz = 1;
		double vx, vy, vz;
		const int EnergiesChannels = NUM_OF_ENERGYCHANNEL; // We totally define NUM_OF_ENERGYCHANNEL channels;
		const int matrixRowNum = DNi * DNj;
		const int matrixColNum = EnergiesChannels;


		//Read the spectrum file.
		double* spRatio = new double[EnergiesChannels];
		std::ifstream fin("spectrum.raw", std::ios::binary);
		if (!fin.is_open())
		{
			std::cout << "Cannot open the spectrum data\n";
			std::cout << "We will use uniform distribution\n";
			for (size_t i = 0; i < EnergiesChannels; i++)
			{
				spRatio[i] = 1.0 / static_cast<double>(EnergiesChannels);
			}
		}
		else
		{
			fin.read((char*) spRatio, sizeof(double) * EnergiesChannels);
		}
		double* channelRatio[8];
		for (size_t i = 0; i < 8; i++)
		{
			channelRatio[i] = new double[EnergiesChannels];
			for (size_t j = 0; j < EnergiesChannels; j++)
			{
				channelRatio[i][j] = 0;
			}
		}
		std::copy(spRatio, spRatio + 110, channelRatio[0]);
		std::copy(spRatio + 110, spRatio + 210, channelRatio[1] + 110);
		std::copy(spRatio + 210, spRatio + 310, channelRatio[2] + 210);
		std::copy(spRatio + 310, spRatio + 410, channelRatio[3] + 310);
		std::copy(spRatio + 410, spRatio + 510, channelRatio[4] + 410);
		std::copy(spRatio + 510, spRatio + 610, channelRatio[5] + 510);
		std::copy(spRatio + 610, spRatio + 710, channelRatio[6] + 610);
		std::copy(spRatio + 710, spRatio + NUM_OF_ENERGYCHANNEL, channelRatio[7] + 710);


		const double sumSpectrum = thrust::reduce(spRatio, spRatio + EnergiesChannels);
		double sumChannelSpectrum[8] = { 0 };
		sumChannelSpectrum[0] = thrust::reduce(channelRatio[0], channelRatio[0] + EnergiesChannels);
		sumChannelSpectrum[1] = thrust::reduce(channelRatio[1], channelRatio[1] + EnergiesChannels);
		sumChannelSpectrum[2] = thrust::reduce(channelRatio[2], channelRatio[2] + EnergiesChannels);
		sumChannelSpectrum[3] = thrust::reduce(channelRatio[3], channelRatio[3] + EnergiesChannels);
		sumChannelSpectrum[4] = thrust::reduce(channelRatio[4], channelRatio[4] + EnergiesChannels);
		sumChannelSpectrum[5] = thrust::reduce(channelRatio[5], channelRatio[5] + EnergiesChannels);
		sumChannelSpectrum[6] = thrust::reduce(channelRatio[6], channelRatio[6] + EnergiesChannels);
		sumChannelSpectrum[7] = thrust::reduce(channelRatio[7], channelRatio[7] + EnergiesChannels);

		double* dspRatio;
		checkCudaErrors(cudaMalloc((void**) &dspRatio, sizeof(double) * EnergiesChannels));
		checkCudaErrors(cudaMemcpy(dspRatio, spRatio, sizeof(double) *EnergiesChannels, cudaMemcpyHostToDevice));
		double* multiSpectrumProj = new double[matrixRowNum];
		double* dmultiSpectrumProj;
		checkCudaErrors(cudaMalloc((void**) &dmultiSpectrumProj, sizeof(double) * matrixRowNum));
		double* channelProj[8];
		double* dchannelProj[8];
		double* dchannelRatio[8];
		for (size_t i = 0; i < 8; i++)
		{
			channelProj[i] = new double[matrixRowNum];
			checkCudaErrors(cudaMalloc((void**) &dchannelProj[i], sizeof(double) * matrixRowNum));

			checkCudaErrors(cudaMalloc((void**) &dchannelRatio[i], sizeof(double) * EnergiesChannels));
			checkCudaErrors(cudaMemcpy(dchannelRatio[i], channelRatio[i], sizeof(double) * EnergiesChannels, cudaMemcpyHostToDevice));
		}


		double* proj; //All energies channels
		checkCudaErrors(cudaMalloc((void**) &proj, sizeof(double) * DNi * DNj * EnergiesChannels));
		checkCudaErrors(cudaMemset(proj, 0, sizeof(double) * DNi * DNj * EnergiesChannels));

		int* M; //Different materials
		checkCudaErrors(cudaMalloc((void**) &M, sizeof(int) * 10));
		checkCudaErrors(cudaMemcpy(M, cyl.M, sizeof(int) * 10, cudaMemcpyHostToDevice));

		double* haxisX = new double[10];
		double* haxisY = new double[10];
		double* haxisZ = new double[10];
		double* hcoreX = new double[10];
		double* hcoreY = new double[10];
		double* hcoreZ = new double[10];
		double* hrots = new double[10];
		for (size_t i = 0; i < 10; i++)
		{
			haxisX[i] = cyl.axis[i].x;
			haxisY[i] = cyl.axis[i].y;
			haxisZ[i] = cyl.axis[i].z;
			hcoreX[i] = cyl.core[i].x;
			hcoreY[i] = cyl.core[i].y;
			hcoreZ[i] = cyl.core[i].z;
			hrots[i] = cyl.rots[i];
		}


		double* axisX;
		double* axisY;
		double* axisZ;
		double* coreX;
		double* coreY;
		double* coreZ;
		double* rots;
		checkCudaErrors(cudaMalloc((void**) &axisX, sizeof(double) * 10));
		checkCudaErrors(cudaMalloc((void**) &axisY, sizeof(double) * 10));
		checkCudaErrors(cudaMalloc((void**) &axisZ, sizeof(double) * 10));
		checkCudaErrors(cudaMalloc((void**) &coreX, sizeof(double) * 10));
		checkCudaErrors(cudaMalloc((void**) &coreY, sizeof(double) * 10));
		checkCudaErrors(cudaMalloc((void**) &coreZ, sizeof(double) * 10));
		checkCudaErrors(cudaMalloc((void**) &rots, sizeof(double) * 10));
		checkCudaErrors(cudaMemcpy(axisX, haxisX, sizeof(double) * 10, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(axisY, haxisY, sizeof(double) * 10, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(axisZ, haxisZ, sizeof(double) * 10, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(coreX, hcoreX, sizeof(double) * 10, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(coreY, hcoreY, sizeof(double) * 10, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(coreZ, hcoreZ, sizeof(double) * 10, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(rots, hrots, sizeof(double) * 10, cudaMemcpyHostToDevice));

		double* materials;
		checkCudaErrors(cudaMalloc((void**) &materials, sizeof(double) * 7 * NUM_OF_ENERGYCHANNEL));
		checkCudaErrors(cudaMemcpy(materials, cyl.materials, sizeof(double) * 7 * NUM_OF_ENERGYCHANNEL, cudaMemcpyHostToDevice));


		dim3 blk(4, 4, 2);
		dim3 gid((DNi + blk.x - 1) / blk.x,
			(DNj + blk.y - 1) / blk.y,
			(EnergiesChannels + blk.z - 1) / blk.z);
		double cosT, sinT;
		double* hostproj = new double[DNi * DNj * EnergiesChannels];
		for (unsigned int i = 0; i != angN; ++i)
		{
			cosT = cos(i * angStp);
			sinT = sin(i * angStp);

			sourx = initx * cosT - inity * sinT;
			soury = initx * sinT + inity * cosT;
			sourz = initz;
			detcntx = initdetcntx * cosT - initdetcnty * sinT;
			detcnty = initdetcntx * sinT + initdetcnty * cosT;
			detcntz = initdetcntz;

			ux = initux * cosT - inituy * sinT;
			uy = initux * sinT + inituy * cosT;
			uz = inituz;

			vx = initvx * cosT - initvy * sinT;
			vy = initvx * sinT + initvy * cosT;
			vz = initvz;


			projHeadPhantomMultiEnergies(proj, sourx, soury, sourz, detcntx, detcnty, detcntz,
				ux, uy, uz, vx, vy, vz, detustp, detvstp, DNi, DNj, M, materials, EnergiesChannels,
				axisX, axisY, axisZ, coreX, coreY, coreZ, rots, blk, gid);

			checkCudaErrors(cudaMemcpy(hostproj, proj, sizeof(double) * DNi * DNj * EnergiesChannels, cudaMemcpyDeviceToHost)); //if you want to save all the channels
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dspRatio, 1, &ZERO, dmultiSpectrumProj, 1);

			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[0], 1, &ZERO, dchannelProj[0], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[1], 1, &ZERO, dchannelProj[1], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[2], 1, &ZERO, dchannelProj[2], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[3], 1, &ZERO, dchannelProj[3], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[4], 1, &ZERO, dchannelProj[4], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[5], 1, &ZERO, dchannelProj[5], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[6], 1, &ZERO, dchannelProj[6], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[7], 1, &ZERO, dchannelProj[7], 1);
			checkCudaErrors(cudaMemcpy(multiSpectrumProj, dmultiSpectrumProj, sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[0], dchannelProj[0], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[1], dchannelProj[1], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[2], dchannelProj[2], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[3], dchannelProj[3], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[4], dchannelProj[4], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[5], dchannelProj[5], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[6], dchannelProj[6], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[7], dchannelProj[7], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));


			thrust::transform(multiSpectrumProj, multiSpectrumProj + matrixRowNum, multiSpectrumProj, divideConstValue_functor<double>(sumSpectrum));
			thrust::transform(channelProj[0], channelProj[0] + matrixRowNum, channelProj[0], divideConstValue_functor<double>(sumChannelSpectrum[0]));
			thrust::transform(channelProj[1], channelProj[1] + matrixRowNum, channelProj[1], divideConstValue_functor<double>(sumChannelSpectrum[1]));
			thrust::transform(channelProj[2], channelProj[2] + matrixRowNum, channelProj[2], divideConstValue_functor<double>(sumChannelSpectrum[2]));
			thrust::transform(channelProj[3], channelProj[3] + matrixRowNum, channelProj[3], divideConstValue_functor<double>(sumChannelSpectrum[3]));
			thrust::transform(channelProj[4], channelProj[4] + matrixRowNum, channelProj[4], divideConstValue_functor<double>(sumChannelSpectrum[4]));
			thrust::transform(channelProj[5], channelProj[5] + matrixRowNum, channelProj[5], divideConstValue_functor<double>(sumChannelSpectrum[5]));
			thrust::transform(channelProj[6], channelProj[6] + matrixRowNum, channelProj[6], divideConstValue_functor<double>(sumChannelSpectrum[6]));
			thrust::transform(channelProj[7], channelProj[7] + matrixRowNum, channelProj[7], divideConstValue_functor<double>(sumChannelSpectrum[7]));

			std::stringstream ss;
			ss << i;
			/*std::string ProjName = "G:\\MultiEnergies\\ProjAngle" + ss.str() + ".raw";
			std::ofstream fou(ProjName.c_str(), std::ios::binary);
			fou.write((char*) hostproj, sizeof(double) * DNi * DNj * EnergiesChannels);
			fou.close();*/ //if you want to save all the channels
			std::string chProjName0 = "Channel0Proj" + ss.str() + ".raw";
			std::string chProjName1 = "Channel1Proj" + ss.str() + ".raw";
			std::string chProjName2 = "Channel2Proj" + ss.str() + ".raw";
			std::string chProjName3 = "Channel3Proj" + ss.str() + ".raw";
			std::string chProjName4 = "Channel4Proj" + ss.str() + ".raw";
			std::string chProjName5 = "Channel5Proj" + ss.str() + ".raw";
			std::string chProjName6 = "Channel6Proj" + ss.str() + ".raw";
			std::string chProjName7 = "Channel7Proj" + ss.str() + ".raw";

			std::string chProjName8 = "allChannelsProj" + ss.str() + ".raw";

			std::ofstream fou0(chProjName0.c_str(), std::ios::binary);
			fou0.write((char*) channelProj[0], sizeof(double) * matrixRowNum);
			fou0.close();

			std::ofstream fou1(chProjName1.c_str(), std::ios::binary);
			fou1.write((char*) channelProj[1], sizeof(double) * matrixRowNum);
			fou1.close();

			std::ofstream fou2(chProjName2.c_str(), std::ios::binary);
			fou2.write((char*) channelProj[2], sizeof(double) * matrixRowNum);
			fou2.close();

			std::ofstream fou3(chProjName3.c_str(), std::ios::binary);
			fou3.write((char*) channelProj[3], sizeof(double) * matrixRowNum);
			fou3.close();

			std::ofstream fou4(chProjName4.c_str(), std::ios::binary);
			fou4.write((char*) channelProj[4], sizeof(double) * matrixRowNum);
			fou4.close();

			std::ofstream fou5(chProjName5.c_str(), std::ios::binary);
			fou5.write((char*) channelProj[5], sizeof(double) * matrixRowNum);
			fou5.close();

			std::ofstream fou6(chProjName6.c_str(), std::ios::binary);
			fou6.write((char*) channelProj[6], sizeof(double) * matrixRowNum);
			fou6.close();

			std::ofstream fou7(chProjName7.c_str(), std::ios::binary);
			fou7.write((char*) channelProj[7], sizeof(double) * matrixRowNum);
			fou7.close();

			std::ofstream fou8(chProjName8.c_str(), std::ios::binary);
			fou8.write((char*) multiSpectrumProj, sizeof(double) * matrixRowNum);
			fou8.close();

			std::cout << i << "\n";
		}
	}







	void helicalScanHeadPhantomMultiEnergies(
		cublasHandle_t handle, const double initx, const double inity, const double initz,
		const double initdetcntx, const double initdetcnty, const double initdetcntz, const double pitch,
		const double detustp, const double detvstp,
		const int DNi, const int DNj, const int angN, const double angStp)
	{
		const double ONE = 1.0;
		const double ZERO = 0.0;
		cublasCreate(&handle);
		HeadPhantom3D<double> cyl;
		//CylinderV2<double> cyl;
		double sourx, soury, sourz;
		double detcntx, detcnty, detcntz;

		double initux = 1;	double inituy = 0;	double inituz = 0;
		double ux, uy, uz;
		double initvx = 0;	double initvy = 0;	double initvz = 1;
		double vx, vy, vz;
		const int EnergiesChannels = NUM_OF_ENERGYCHANNEL; // We totally define NUM_OF_ENERGYCHANNEL channels;
		const int matrixRowNum = DNi * DNj;
		const int matrixColNum = EnergiesChannels;


		//Read the spectrum file.
		double* spRatio = new double[EnergiesChannels];
		std::ifstream fin("spectrum.raw", std::ios::binary);
		if (!fin.is_open())
		{
			std::cout << "Cannot open the spectrum data\n";
			std::cout << "We will use uniform distribution\n";
			for (size_t i = 0; i < EnergiesChannels; i++)
			{
				spRatio[i] = 1.0 / static_cast<double>(EnergiesChannels);
			}
		}
		else
		{
			fin.read((char*) spRatio, sizeof(double) * EnergiesChannels);
		}
		double* channelRatio[8];
		for (size_t i = 0; i < 8; i++)
		{
			channelRatio[i] = new double[EnergiesChannels];
			for (size_t j = 0; j < EnergiesChannels; j++)
			{
				channelRatio[i][j] = 0;
			}
		}
		std::copy(spRatio, spRatio + 110, channelRatio[0]);
		std::copy(spRatio + 110, spRatio + 210, channelRatio[1] + 110);
		std::copy(spRatio + 210, spRatio + 310, channelRatio[2] + 210);
		std::copy(spRatio + 310, spRatio + 410, channelRatio[3] + 310);
		std::copy(spRatio + 410, spRatio + 510, channelRatio[4] + 410);
		std::copy(spRatio + 510, spRatio + 610, channelRatio[5] + 510);
		std::copy(spRatio + 610, spRatio + 710, channelRatio[6] + 610);
		std::copy(spRatio + 710, spRatio + NUM_OF_ENERGYCHANNEL, channelRatio[7] + 710);


		const double sumSpectrum = thrust::reduce(spRatio, spRatio + EnergiesChannels);
		double sumChannelSpectrum[8] = { 0 };
		sumChannelSpectrum[0] = thrust::reduce(channelRatio[0], channelRatio[0] + EnergiesChannels);
		sumChannelSpectrum[1] = thrust::reduce(channelRatio[1], channelRatio[1] + EnergiesChannels);
		sumChannelSpectrum[2] = thrust::reduce(channelRatio[2], channelRatio[2] + EnergiesChannels);
		sumChannelSpectrum[3] = thrust::reduce(channelRatio[3], channelRatio[3] + EnergiesChannels);
		sumChannelSpectrum[4] = thrust::reduce(channelRatio[4], channelRatio[4] + EnergiesChannels);
		sumChannelSpectrum[5] = thrust::reduce(channelRatio[5], channelRatio[5] + EnergiesChannels);
		sumChannelSpectrum[6] = thrust::reduce(channelRatio[6], channelRatio[6] + EnergiesChannels);
		sumChannelSpectrum[7] = thrust::reduce(channelRatio[7], channelRatio[7] + EnergiesChannels);

		double* dspRatio;
		checkCudaErrors(cudaMalloc((void**) &dspRatio, sizeof(double) * EnergiesChannels));
		checkCudaErrors(cudaMemcpy(dspRatio, spRatio, sizeof(double) *EnergiesChannels, cudaMemcpyHostToDevice));
		double* multiSpectrumProj = new double[matrixRowNum];
		double* dmultiSpectrumProj;
		checkCudaErrors(cudaMalloc((void**) &dmultiSpectrumProj, sizeof(double) * matrixRowNum));
		double* channelProj[8];
		double* dchannelProj[8];
		double* dchannelRatio[8];
		for (size_t i = 0; i < 8; i++)
		{
			channelProj[i] = new double[matrixRowNum];
			checkCudaErrors(cudaMalloc((void**) &dchannelProj[i], sizeof(double) * matrixRowNum));

			checkCudaErrors(cudaMalloc((void**) &dchannelRatio[i], sizeof(double) * EnergiesChannels));
			checkCudaErrors(cudaMemcpy(dchannelRatio[i], channelRatio[i], sizeof(double) * EnergiesChannels, cudaMemcpyHostToDevice));
		}


		double* proj; //All energies channels
		checkCudaErrors(cudaMalloc((void**) &proj, sizeof(double) * DNi * DNj * EnergiesChannels));
		checkCudaErrors(cudaMemset(proj, 0, sizeof(double) * DNi * DNj * EnergiesChannels));

		int* M; //Different materials
		checkCudaErrors(cudaMalloc((void**) &M, sizeof(int) * 10));
		checkCudaErrors(cudaMemcpy(M, cyl.M, sizeof(int) * 10, cudaMemcpyHostToDevice));

		double* haxisX = new double[10];
		double* haxisY = new double[10];
		double* haxisZ = new double[10];
		double* hcoreX = new double[10];
		double* hcoreY = new double[10];
		double* hcoreZ = new double[10];
		double* hrots = new double[10];
		for (size_t i = 0; i < 10; i++)
		{
			haxisX[i] = cyl.axis[i].x;
			haxisY[i] = cyl.axis[i].y;
			haxisZ[i] = cyl.axis[i].z;
			hcoreX[i] = cyl.core[i].x;
			hcoreY[i] = cyl.core[i].y;
			hcoreZ[i] = cyl.core[i].z;
			hrots[i] = cyl.rots[i];
		}


		double* axisX;
		double* axisY;
		double* axisZ;
		double* coreX;
		double* coreY;
		double* coreZ;
		double* rots;
		checkCudaErrors(cudaMalloc((void**) &axisX, sizeof(double) * 10));
		checkCudaErrors(cudaMalloc((void**) &axisY, sizeof(double) * 10));
		checkCudaErrors(cudaMalloc((void**) &axisZ, sizeof(double) * 10));
		checkCudaErrors(cudaMalloc((void**) &coreX, sizeof(double) * 10));
		checkCudaErrors(cudaMalloc((void**) &coreY, sizeof(double) * 10));
		checkCudaErrors(cudaMalloc((void**) &coreZ, sizeof(double) * 10));
		checkCudaErrors(cudaMalloc((void**) &rots, sizeof(double) * 10));
		checkCudaErrors(cudaMemcpy(axisX, haxisX, sizeof(double) * 10, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(axisY, haxisY, sizeof(double) * 10, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(axisZ, haxisZ, sizeof(double) * 10, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(coreX, hcoreX, sizeof(double) * 10, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(coreY, hcoreY, sizeof(double) * 10, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(coreZ, hcoreZ, sizeof(double) * 10, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(rots, hrots, sizeof(double) * 10, cudaMemcpyHostToDevice));

		double* materials;
		checkCudaErrors(cudaMalloc((void**) &materials, sizeof(double) * 7 * NUM_OF_ENERGYCHANNEL));
		checkCudaErrors(cudaMemcpy(materials, cyl.materials, sizeof(double) * 7 * NUM_OF_ENERGYCHANNEL, cudaMemcpyHostToDevice));


		dim3 blk(4, 4, 2);
		dim3 gid((DNi + blk.x - 1) / blk.x,
			(DNj + blk.y - 1) / blk.y,
			(EnergiesChannels + blk.z - 1) / blk.z);
		double cosT, sinT;
		double* hostproj = new double[DNi * DNj * EnergiesChannels];
		for (unsigned int i = 0; i != angN; ++i)
		{
			cosT = cos(i * angStp);
			sinT = sin(i * angStp);

			sourx = initx * cosT - inity * sinT;
			soury = initx * sinT + inity * cosT;
			sourz = initz + i * angStp * pitch;

			detcntx = initdetcntx * cosT - initdetcnty * sinT;
			detcnty = initdetcntx * sinT + initdetcnty * cosT;
			detcntz = initdetcntz + i * angStp * pitch;

			ux = initux * cosT - inituy * sinT;
			uy = initux * sinT + inituy * cosT;
			uz = inituz;

			vx = initvx * cosT - initvy * sinT;
			vy = initvx * sinT + initvy * cosT;
			vz = initvz;


			projHeadPhantomMultiEnergies(proj, sourx, soury, sourz, detcntx, detcnty, detcntz,
				ux, uy, uz, vx, vy, vz, detustp, detvstp, DNi, DNj, M, materials, EnergiesChannels,
				axisX, axisY, axisZ, coreX, coreY, coreZ, rots, blk, gid);

			checkCudaErrors(cudaMemcpy(hostproj, proj, sizeof(double) * DNi * DNj * EnergiesChannels, cudaMemcpyDeviceToHost)); //if you want to save all the channels
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dspRatio, 1, &ZERO, dmultiSpectrumProj, 1);

			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[0], 1, &ZERO, dchannelProj[0], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[1], 1, &ZERO, dchannelProj[1], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[2], 1, &ZERO, dchannelProj[2], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[3], 1, &ZERO, dchannelProj[3], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[4], 1, &ZERO, dchannelProj[4], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[5], 1, &ZERO, dchannelProj[5], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[6], 1, &ZERO, dchannelProj[6], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[7], 1, &ZERO, dchannelProj[7], 1);

			checkCudaErrors(cudaMemcpy(multiSpectrumProj, dmultiSpectrumProj, sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[0], dchannelProj[0], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[1], dchannelProj[1], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[2], dchannelProj[2], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[3], dchannelProj[3], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[4], dchannelProj[4], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[5], dchannelProj[5], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[6], dchannelProj[6], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[7], dchannelProj[7], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));


			thrust::transform(multiSpectrumProj, multiSpectrumProj + matrixRowNum, multiSpectrumProj, divideConstValue_functor<double>(sumSpectrum));
			thrust::transform(channelProj[0], channelProj[0] + matrixRowNum, channelProj[0], divideConstValue_functor<double>(sumChannelSpectrum[0]));
			thrust::transform(channelProj[1], channelProj[1] + matrixRowNum, channelProj[1], divideConstValue_functor<double>(sumChannelSpectrum[1]));
			thrust::transform(channelProj[2], channelProj[2] + matrixRowNum, channelProj[2], divideConstValue_functor<double>(sumChannelSpectrum[2]));
			thrust::transform(channelProj[3], channelProj[3] + matrixRowNum, channelProj[3], divideConstValue_functor<double>(sumChannelSpectrum[3]));
			thrust::transform(channelProj[4], channelProj[4] + matrixRowNum, channelProj[4], divideConstValue_functor<double>(sumChannelSpectrum[4]));
			thrust::transform(channelProj[5], channelProj[5] + matrixRowNum, channelProj[5], divideConstValue_functor<double>(sumChannelSpectrum[5]));
			thrust::transform(channelProj[6], channelProj[6] + matrixRowNum, channelProj[6], divideConstValue_functor<double>(sumChannelSpectrum[6]));
			thrust::transform(channelProj[7], channelProj[7] + matrixRowNum, channelProj[7], divideConstValue_functor<double>(sumChannelSpectrum[7]));

			std::stringstream ss;
			ss << i;
			/*std::string ProjName = "G:\\MultiEnergies\\ProjAngle" + ss.str() + ".raw";
			std::ofstream fou(ProjName.c_str(), std::ios::binary);
			fou.write((char*) hostproj, sizeof(double) * DNi * DNj * EnergiesChannels);
			fou.close();*/ //if you want to save all the channels
			std::string chProjName0 = "Channel0Proj" + ss.str() + ".raw";
			std::string chProjName1 = "Channel1Proj" + ss.str() + ".raw";
			std::string chProjName2 = "Channel2Proj" + ss.str() + ".raw";
			std::string chProjName3 = "Channel3Proj" + ss.str() + ".raw";
			std::string chProjName4 = "Channel4Proj" + ss.str() + ".raw";
			std::string chProjName5 = "Channel5Proj" + ss.str() + ".raw";
			std::string chProjName6 = "Channel6Proj" + ss.str() + ".raw";
			std::string chProjName7 = "Channel7Proj" + ss.str() + ".raw";

			std::string chProjName8 = "allChannelsProj" + ss.str() + ".raw";

			std::ofstream fou0(chProjName0.c_str(), std::ios::binary);
			fou0.write((char*) channelProj[0], sizeof(double) * matrixRowNum);
			fou0.close();

			std::ofstream fou1(chProjName1.c_str(), std::ios::binary);
			fou1.write((char*) channelProj[1], sizeof(double) * matrixRowNum);
			fou1.close();

			std::ofstream fou2(chProjName2.c_str(), std::ios::binary);
			fou2.write((char*) channelProj[2], sizeof(double) * matrixRowNum);
			fou2.close();

			std::ofstream fou3(chProjName3.c_str(), std::ios::binary);
			fou3.write((char*) channelProj[3], sizeof(double) * matrixRowNum);
			fou3.close();

			std::ofstream fou4(chProjName4.c_str(), std::ios::binary);
			fou4.write((char*) channelProj[4], sizeof(double) * matrixRowNum);
			fou4.close();

			std::ofstream fou5(chProjName5.c_str(), std::ios::binary);
			fou5.write((char*) channelProj[5], sizeof(double) * matrixRowNum);
			fou5.close();

			std::ofstream fou6(chProjName6.c_str(), std::ios::binary);
			fou6.write((char*) channelProj[6], sizeof(double) * matrixRowNum);
			fou6.close();

			std::ofstream fou7(chProjName7.c_str(), std::ios::binary);
			fou7.write((char*) channelProj[7], sizeof(double) * matrixRowNum);
			fou7.close();

			std::ofstream fou8(chProjName8.c_str(), std::ios::binary);
			fou8.write((char*) multiSpectrumProj, sizeof(double) * matrixRowNum);
			fou8.close();

			std::cout << i << "\n";
		}
	}











	void maanscanCylinder()
	{
		Cylinder<double> cyl;
		double initx = 0.0;
		double inity = 5.0;
		double initz = 0.0;
		double sourx, soury, sourz;
		double initdetcntx = 0;
		double initdetcnty = -5.0;
		double initdetcntz = 0.0;


		double detcntx, detcnty, detcntz;

		double initux = 1;	double inituy = 0;	double inituz = 0;
		double ux, uy, uz;
		double initvx = 0;	double initvy = 0;	double initvz = 1;
		double vx, vy, vz;
		const double detustp = 0.01;
		const double detvstp = 0.01;
		const int DNi = 512;
		const int DNj = 512;
		const int angN = 360;
		const double angStp = 3.14159265358979323846264 * 2.0 / angN;

		double* proj;
		checkCudaErrors(cudaMalloc((void**) &proj, sizeof(double) * DNi * DNj * angN));
		checkCudaErrors(cudaMemset(proj, 0, sizeof(double) * DNi * DNj * angN));

		double* A;
		checkCudaErrors(cudaMalloc((void**) &A, sizeof(double) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(A, cyl.A, sizeof(double) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		double* R;
		checkCudaErrors(cudaMalloc((void**) &R, sizeof(double) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(R, cyl.R, sizeof(double) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		double* posx;
		checkCudaErrors(cudaMalloc((void**) &posx, sizeof(double) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(posx, cyl.posx, sizeof(double) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		double* posy;
		checkCudaErrors(cudaMalloc((void**) &posy, sizeof(double) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(posy, cyl.posy, sizeof(double) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		double height = cyl.height;

		dim3 blk(16, 16);
		dim3 gid((DNi + blk.x - 1) / blk.x,
			(DNj + blk.y - 1) / blk.y);
		double cosT, sinT, cos2T, sin2T;
		for (unsigned int i = 0; i != angN; ++i)
		{

			initx = 0.0;
			inity = 5.0;
			initz = 0.0;

			initdetcntx = 0;
			initdetcnty = -5.0;
			initdetcntz = 0.0;
			initux = 1;	inituy = 0;	inituz = 0;
			initvx = 0;	initvy = 0;	initvz = 1;


			cosT = cos(i * angStp);
			sinT = sin(i * angStp);
			cos2T = cos(i * angStp * 2);
			sin2T = sin(i * angStp * 2);

			sourx = initx;
			soury = inity * cos2T - initz * sin2T;
			sourz = inity * sin2T + initz * cos2T;

			initx = sourx * cosT - soury * sinT;
			inity = sourx * sinT + soury * cosT;
			initz = sourz;

			sourx = initx;
			soury = inity;
			sourz = initz;



			detcntx = initdetcntx;
			detcnty = initdetcnty * cos2T - initdetcntz * sin2T;
			detcntz = initdetcnty * sin2T + initdetcntz * cos2T;

			initdetcntx = detcntx * cosT - detcnty * sinT;
			initdetcnty = detcntx * sinT + detcnty * cosT;
			initdetcntz = detcntz;

			detcntx = initdetcntx;
			detcnty = initdetcnty;
			detcntz = initdetcntz;



			ux = initux;
			uy = inituy * cos2T - inituz * sin2T;
			uz = inituy * sin2T + inituz * cos2T;

			initux = ux * cosT - uy * sinT;
			inituy = ux * sinT + uy * cosT;
			inituz = uz;

			ux = initux;
			uy = inituy;
			uz = inituz;


			vx = initvx;
			vy = initvy * cos2T - initvz * sin2T;
			vz = initvy * sin2T + initvz * cos2T;

			initvx = vx * cosT - vy * sinT;
			initvy = vx * sinT + vy * cosT;
			initvz = vz;

			vx = initvx;
			vy = initvy;
			vz = initvz;


			projCylinder<double>(proj + i * DNi * DNj, sourx, soury, sourz,
				detcntx, detcnty, detcntz, ux, uy, uz, vx, vy, vz, detustp, detvstp,
				DNi, DNj, A, R, height, posx, posy, blk, gid);

			std::cout << i << "\n";
		}

		double* hproj = new double[DNi * DNj * angN];
		checkCudaErrors(cudaMemcpy(hproj, proj, sizeof(double) * DNi * DNj * angN, cudaMemcpyDeviceToHost));

		std::ofstream fou("maanScanningRes.raw", std::ios::binary);
		fou.write((char*) hproj, sizeof(double) * DNi * DNj * angN);
		fou.close();
	}




	// 从NUM_OF_ENERGYCHANNEL个通道合成一个通道，并输出8个通道的值;
	void conebeamScanCylinderMultiEnergiesNew(
		cublasHandle_t handle, const double initx, const double inity, const double initz,
		const double initdetcntx, const double initdetcnty, const double initdetcntz,
		const double detustp, const double detvstp,
		const int DNi, const int DNj, const int angN, const double angStp)
	{
		const double ONE = 1.0;
		const double ZERO = 0.0;
		cublasCreate(&handle);
		CylinderV2<double> cyl; //This will read another file
		double sourx, soury, sourz;
		double detcntx, detcnty, detcntz;

		double initux = 1;	double inituy = 0;	double inituz = 0;
		double ux, uy, uz;
		double initvx = 0;	double initvy = 0;	double initvz = 1;
		double vx, vy, vz;
		const int EnergiesChannels = NUM_OF_ENERGYCHANNEL; // We totally define NUM_OF_ENERGYCHANNEL channels;
		const int matrixRowNum = DNi * DNj;
		const int matrixColNum = EnergiesChannels;

		//Read the photon number file
		double* photonNum = new double[EnergiesChannels];
		std::ifstream fin("photonNum.raw", std::ios::binary);
		if (!fin.is_open()) { std::cout << "Cannot open the photon number file\n"; exit(-1); }
		fin.read((char*) photonNum, sizeof(double) * EnergiesChannels);
		fin.close();



		double* channelRatio[8];
		for (size_t i = 0; i < 8; i++)
		{
			channelRatio[i] = new double[EnergiesChannels];
			for (size_t j = 0; j < EnergiesChannels; j++)
			{
				channelRatio[i][j] = 0;
			}
		}
		//Use this configuration to guarantee same photons in each channel
		std::fill(channelRatio[0], channelRatio[0] + 110, 1.0);
		std::fill(channelRatio[1] + 110, channelRatio[1] + 145, 1.0);
		std::fill(channelRatio[2] + 145, channelRatio[2] + 180, 1.0);
		std::fill(channelRatio[3] + 180, channelRatio[3] + 221, 1.0);
		std::fill(channelRatio[4] + 221, channelRatio[4] + 275, 1.0);
		std::fill(channelRatio[5] + 275, channelRatio[5] + 319, 1.0);
		std::fill(channelRatio[6] + 319, channelRatio[6] + 419, 1.0);
		std::fill(channelRatio[7] + 419, channelRatio[7] + NUM_OF_ENERGYCHANNEL, 1.0);

		double* dphotonNum; //各个通道的光子数;
		checkCudaErrors(cudaMalloc((void**) &dphotonNum, sizeof(double) * EnergiesChannels));
		checkCudaErrors(cudaMemcpy(dphotonNum, photonNum, sizeof(double) *EnergiesChannels, cudaMemcpyHostToDevice));

		double* spRatio = new double[EnergiesChannels];
		for (unsigned int i = 0; i != EnergiesChannels; ++i) { spRatio[i] = 1.0; }
		double* dspRatio;
		checkCudaErrors(cudaMalloc((void**) &dspRatio, sizeof(double) * EnergiesChannels));
		checkCudaErrors(cudaMemcpy(dspRatio, spRatio, sizeof(double) * EnergiesChannels, cudaMemcpyHostToDevice));



		double* multiSpectrumProj = new double[matrixRowNum];
		double* dmultiSpectrumProj;
		checkCudaErrors(cudaMalloc((void**) &dmultiSpectrumProj, sizeof(double) * matrixRowNum));
		double* channelProj[8];
		double* dchannelProj[8];
		double* dchannelRatio[8];
		for (size_t i = 0; i < 8; i++)
		{
			channelProj[i] = new double[matrixRowNum];
			checkCudaErrors(cudaMalloc((void**) &dchannelProj[i], sizeof(double) * matrixRowNum));

			checkCudaErrors(cudaMalloc((void**) &dchannelRatio[i], sizeof(double) * EnergiesChannels));
			checkCudaErrors(cudaMemcpy(dchannelRatio[i], channelRatio[i], sizeof(double) * EnergiesChannels, cudaMemcpyHostToDevice));
		}


		double* proj; //All energies channels
		checkCudaErrors(cudaMalloc((void**) &proj, sizeof(double) * DNi * DNj * EnergiesChannels));
		checkCudaErrors(cudaMemset(proj, 0, sizeof(double) * DNi * DNj * EnergiesChannels));

		int* A; //Different materials
		checkCudaErrors(cudaMalloc((void**) &A, sizeof(int) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(A, cyl.A, sizeof(int) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		double* R;
		checkCudaErrors(cudaMalloc((void**) &R, sizeof(double) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(R, cyl.R, sizeof(double) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		double* posx;
		checkCudaErrors(cudaMalloc((void**) &posx, sizeof(double) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(posx, cyl.posx, sizeof(double) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		double* posy;
		checkCudaErrors(cudaMalloc((void**) &posy, sizeof(double) * NUM_OF_CYLINDER));
		checkCudaErrors(cudaMemcpy(posy, cyl.posy, sizeof(double) * NUM_OF_CYLINDER, cudaMemcpyHostToDevice));

		double* materials;
		checkCudaErrors(cudaMalloc((void**) &materials, sizeof(double) * 7 * NUM_OF_ENERGYCHANNEL));
		checkCudaErrors(cudaMemcpy(materials, cyl.materials, sizeof(double) * 7 * NUM_OF_ENERGYCHANNEL, cudaMemcpyHostToDevice));

		double height = cyl.height;

		dim3 blk(8, 8, 4);
		dim3 gid((DNi + blk.x - 1) / blk.x,
			(DNj + blk.y - 1) / blk.y,
			(EnergiesChannels + blk.z - 1) / blk.z);
		double cosT, sinT;
		double* hostproj = new double[DNi * DNj * EnergiesChannels];
		for (unsigned int i = 0; i != angN; ++i)
		{
			cosT = cos(i * angStp);
			sinT = sin(i * angStp);

			sourx = initx * cosT - inity * sinT;
			soury = initx * sinT + inity * cosT;
			sourz = initz;
			detcntx = initdetcntx * cosT - initdetcnty * sinT;
			detcnty = initdetcntx * sinT + initdetcnty * cosT;
			detcntz = initdetcntz;

			ux = initux * cosT - inituy * sinT;
			uy = initux * sinT + inituy * cosT;
			uz = inituz;

			vx = initvx * cosT - initvy * sinT;
			vy = initvx * sinT + initvy * cosT;
			vz = initvz;

			projCylinderMultiEnergies<double>(proj, sourx, soury, sourz,
				detcntx, detcnty, detcntz, ux, uy, uz, vx, vy, vz,
				detustp, detvstp, DNi, DNj, dphotonNum, A, materials, EnergiesChannels, R, height, posx, posy, blk, gid);


			// Get I = I0 * exp(-\mu)


			checkCudaErrors(cudaMemcpy(hostproj, proj, sizeof(double) * DNi * DNj * EnergiesChannels, cudaMemcpyDeviceToHost));
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dspRatio, 1, &ZERO, dmultiSpectrumProj, 1);

			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[0], 1, &ZERO, dchannelProj[0], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[1], 1, &ZERO, dchannelProj[1], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[2], 1, &ZERO, dchannelProj[2], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[3], 1, &ZERO, dchannelProj[3], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[4], 1, &ZERO, dchannelProj[4], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[5], 1, &ZERO, dchannelProj[5], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[6], 1, &ZERO, dchannelProj[6], 1);
			cublasDgemv(handle, CUBLAS_OP_N, matrixRowNum, matrixColNum, &ONE, proj, matrixRowNum, dchannelRatio[7], 1, &ZERO, dchannelProj[7], 1);
			checkCudaErrors(cudaMemcpy(multiSpectrumProj, dmultiSpectrumProj, sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[0], dchannelProj[0], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[1], dchannelProj[1], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[2], dchannelProj[2], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[3], dchannelProj[3], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[4], dchannelProj[4], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[5], dchannelProj[5], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[6], dchannelProj[6], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(channelProj[7], dchannelProj[7], sizeof(double) * matrixRowNum, cudaMemcpyDeviceToHost));

			std::stringstream ss;
			ss << i;
			/*std::string ProjName = "G:\\MultiEnergies\\ProjAngle" + ss.str() + ".raw";
			std::ofstream fou(ProjName.c_str(), std::ios::binary);
			fou.write((char*) hostproj, sizeof(double) * DNi * DNj * EnergiesChannels);
			fou.close();*/
			std::string chProjName0 = "photonNumProj" + ss.str() + ".raw";
			std::string chProjName1 = "photonNumProj" + ss.str() + ".raw";
			std::string chProjName2 = "photonNumProj" + ss.str() + ".raw";
			std::string chProjName3 = "photonNumProj" + ss.str() + ".raw";
			std::string chProjName4 = "photonNumProj" + ss.str() + ".raw";
			std::string chProjName5 = "photonNumProj" + ss.str() + ".raw";
			std::string chProjName6 = "photonNumProj" + ss.str() + ".raw";
			std::string chProjName7 = "photonNumProj" + ss.str() + ".raw";

			std::string chProjName8 = "G:\\NEW_MultiEnergies\\allChannels\\photonNumProj" + ss.str() + ".raw";

			std::ofstream fou0(chProjName0.c_str(), std::ios::binary);
			fou0.write((char*) channelProj[0], sizeof(double) * matrixRowNum);
			fou0.close();

			std::ofstream fou1(chProjName1.c_str(), std::ios::binary);
			fou1.write((char*) channelProj[1], sizeof(double) * matrixRowNum);
			fou1.close();

			std::ofstream fou2(chProjName2.c_str(), std::ios::binary);
			fou2.write((char*) channelProj[2], sizeof(double) * matrixRowNum);
			fou2.close();

			std::ofstream fou3(chProjName3.c_str(), std::ios::binary);
			fou3.write((char*) channelProj[3], sizeof(double) * matrixRowNum);
			fou3.close();

			std::ofstream fou4(chProjName4.c_str(), std::ios::binary);
			fou4.write((char*) channelProj[4], sizeof(double) * matrixRowNum);
			fou4.close();

			std::ofstream fou5(chProjName5.c_str(), std::ios::binary);
			fou5.write((char*) channelProj[5], sizeof(double) * matrixRowNum);
			fou5.close();

			std::ofstream fou6(chProjName6.c_str(), std::ios::binary);
			fou6.write((char*) channelProj[6], sizeof(double) * matrixRowNum);
			fou6.close();

			std::ofstream fou7(chProjName7.c_str(), std::ios::binary);
			fou7.write((char*) channelProj[7], sizeof(double) * matrixRowNum);
			fou7.close();

			std::ofstream fou8(chProjName8.c_str(), std::ios::binary);
			fou8.write((char*) multiSpectrumProj, sizeof(double) * matrixRowNum);
			fou8.close();

			std::cout << i << "\n";
		}
	}



	void testNewGeneratePhotonNumProj()
	{
		cublasHandle_t handle;
		cublasCreate(&handle);
		const double initx = 0;
		const double inity = 5.0;
		const double initz = 0;
		const double initdetcntx = 0;
		const double initdetcnty = -5.0;
		const double initdetcntz = 0;
		const double detustp = 0.008;
		const double detvstp = 0.008;
		const int DNi = 512;
		const int DNj = 512;
		const int angN = 360;
		const double angStp = 0.01745329251994329576923688888889 * 2;
		//const double pitch = 2;
		conebeamScanCylinderMultiEnergiesNew(handle, initx, inity, initz, initdetcntx, initdetcnty, initdetcntz,
			detustp, detvstp, DNi, DNj, angN, angStp);
	}

}

