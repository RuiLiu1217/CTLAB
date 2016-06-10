#include "DD3_GPU_recon.h"
#include "DD3_GPU_proj.h"
#include "DD3_GPU_back.h"
#include "FastIterativeShrinkageAlgorithm.h"
#include "SARTWeighting.h"


#include <sstream>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/transform.h>
#include <thrust/tuple.h>
#include <thrust/inner_product.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>


// Evaluate the projection/backprojection with all one input
void evalProjBack(
		bool forwardProjection,
		ForwardDDMethod forwardMethod = PROJ_BRANCHLESS,
		BackwardDDMethod backwardMethod = BACK_BRANCHLESS)
{
	float sid = 541.0f;
	float sdd = 949.0f;

	float x0(0.0f);
	float y0(sid);
	float z0(0.0f);


	int DNU(888);
	int DNV(64);
	int PN(1200);
	float imgXCenter(0);
	float imgYCenter(0);
	float imgZCenter(0);

	int XN(512);
	int YN(512);
	int ZN(64);


	float dx(500.0f / 512.0f);
	float dz(0.625);


	float* xds = new float[DNU];
	float* yds = new float[DNU];
	float* zds = new float[DNV];
	//Generate the positions of the detectors
	float col_size = 1.0239;
	float row_size = 1.0963;

	float col_offset = 0;
	float row_offset = 0;


	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
	float curBeta = 0;
	for (int ii = 0; ii != DNU; ++ii)
	{
		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
		xds[ii] = sinf(curBeta) * sdd;
		yds[ii] = sid - cosf(curBeta) * sdd;
	}

	for (int ii = 0; ii != DNV; ++ii)
	{
		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
	}

	imgXCenter = 0;
	imgYCenter = 0;
	imgZCenter = 0;

	float* hangs = new float[PN];
	float* hzPos = new float[PN];

	for (int ii = 0; ii != PN; ++ii)
	{
		hangs[ii] = ii * TWOPI / static_cast<float>(PN);
		hzPos[ii] = 0;// (ii - PN / 2) * 0.0015;
	}


	byte* mask = new byte[XN * YN];
	for (int i = 0; i != XN * YN; ++i)
	{
		mask[i] = 1;
	}

	float* hvol = new float[XN * YN * ZN];
	float* hprj = new float[DNU * DNV * PN];


	if(forwardProjection)
	{
		std::generate(hvol,hvol+XN*YN*ZN,[](){return 1.0f;});
		std::generate(hprj, hprj + DNU * DNV * PN,[](){return 0.0f;});
		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
			imgXCenter, imgYCenter, imgZCenter,
			hangs, hzPos, PN, XN, YN, ZN, hvol, hprj,
			dx, dz, mask, 0, forwardMethod);
		std::ofstream fou("evalProjResult.prj", std::ios::binary);
		fou.write((char*)hprj, sizeof(float) * DNU * DNV * PN);
		fou.close();
	}
	else
	{
		std::generate(hvol,hvol+XN*YN*ZN,[](){return 0.0f;});
		std::generate(hprj, hprj + DNU * DNV * PN,[](){return 1.0f;});
		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
			imgXCenter, imgYCenter, imgZCenter,
			hangs, hzPos, PN, XN, YN, ZN, hvol, hprj,
			dx, dz, mask, 0, 0, backwardMethod);
		std::ofstream fou("evalBackResult.raw", std::ios::binary);
		fou.write((char*)hvol, sizeof(float) * XN * YN * ZN);
		fou.close();
	}
}



// The OS-SART algorithm with one GPU (Titan X)
void OS_SART(thrust::host_vector<float>& reconImg, // The image volume to be reconstructed
		thrust::host_vector<float>& initImg, // The initial image volume
		thrust::host_vector<float>& hprj, // Projection data
		thrust::host_vector<float>& hangs, // projection views
		thrust::host_vector<float>& hzPos, // source positions in different views
		thrust::host_vector<byte>& mask, // image mask along in-plane direction
		const std::string& VolumeName, // The reconstructed volume name to be stored (Z, X, Y) order in float datatype
		const int osNum = 20, // # of OS
		const int iterNum = 30, // # of total iterations
		const bool outputMedRes = true, // output the intermediate results
		const float sid=538.5200193125, // source to iso-center distance
		const float sdd=946.745971679699, // source to detector distance
		const int DNU = 888, // # of detector elements along in-plane direction
		const int DNV = 64, // # of detector elements along bench moving direction
		const int PN = 8139, // # total view #
		const float imgXCenter = 0.0f, // volume center coordinate x
		const float imgYCenter = 0.0f, // volume center coordinate y
		const float imgZCenter = 0.0f, // volume center coordinate z
		const int XN = 512, // # of volume pixel along x
		const int YN = 512, // # of volume pixel along y
		const int ZN = 512, // # of volume pixel along z
		const float dx = 1.171875, // image pixel size along in-plane direction
		const float dz = 1.171875, // slice thickness
		const float col_size = 1.0239f, // detector elements size along in-plane direction
		const float row_size = 1.0963f, // detector elements size along bench moving direction
		const float col_offset = -1.25f, // #(rational) of detector elements that offset from the center of the iso-center along in-plane direction
		const float row_offset = 0.0f, // #(rational) of detector elements that offset from the center of the iso-center along bench moving direction
		const ForwardDDMethod forwardMethod = PROJ_BRANCHLESS, // Forward projection model
		const BackwardDDMethod backwardMethod = BACK_BRANCHLESS)
{
	// Source position
	const float x0 = 0.0f;
	const float y0 = sid;
	const float z0 = 0.0f;

	thrust::host_vector<float> xds(DNU, 0);
	thrust::host_vector<float> yds(DNU, 0);
	thrust::host_vector<float> zds(DNV, 0);


	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
	int ii = 0;
	thrust::generate(
			thrust::make_zip_iterator(thrust::make_tuple(xds.begin(),yds.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(xds.end(),yds.end())),
			[&ii, DNU, col_offset, stepTheta, sdd, sid](){
		auto res = thrust::make_tuple(sinf((ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta) * sdd,
				sid - cosf((ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta) * sdd);
		ii = ii + 1;
		return res;
	});

	ii = 0;
	thrust::generate(zds.begin(), zds.end(), [&ii, DNV, row_offset, row_size](){
		return (ii++ - (DNV - 1.0) * 0.5 + row_offset) * row_size;
	});

	thrust::host_vector<float> hvol = initImg;

	const int viewPersubSet = PN / osNum; // How many views are in one subset
	const int res = PN - viewPersubSet * osNum; //How many views are more in the last subset
	const int viewInLastSubSet = res + viewPersubSet;

	const thrust::host_vector<float>::size_type totVN = reconImg.size();
	const thrust::host_vector<float>::size_type totJN = hprj.size();
	const int DDNPV = DNU * DNV;

	int* SPN = new int[osNum];
	float** shangs = new float*[osNum];
	float** shzPos = new float*[osNum];
	float** shprj = new float*[osNum];
	float** srowSum = new float*[osNum];
	//float* oneVol = new float[totVN];

	float** onePrj = new float*[osNum];
	float** colSum = new float*[osNum];

	thrust::host_vector<float> oneVol(totVN, 1.0f);
//	for(thrust::host_vector<float>::size_type i = 0; i != totVN; ++i)
//	{
//		oneVol[i] = 1.0;
//	}
	for(int i = 0 ; i != osNum-1; ++i)
	{
		SPN[i] = viewPersubSet;
	}
	SPN[osNum - 1] = viewInLastSubSet;

	for(int i = 0; i != osNum; ++i)
	{
		shangs[i] = new float[SPN[i]];
		shzPos[i] = new float[SPN[i]];
		shprj[i] = new float[SPN[i] * DDNPV];
		srowSum[i] = new float[SPN[i] * DDNPV];
		onePrj[i] = new float[SPN[i] * DDNPV];
		colSum[i] = new	float[totVN];
		for(int curIdx = 0; curIdx != viewPersubSet; ++curIdx)
		{
			shangs[i][curIdx] = hangs[curIdx * osNum + i];
			shzPos[i][curIdx] = hzPos[curIdx * osNum + i];
			memcpy(shprj[i] + curIdx * DDNPV, &hprj[0] + (curIdx * osNum + i) * DDNPV, sizeof(float) * DDNPV);

		}

	}
	memcpy(shangs[osNum-1] + viewPersubSet, &hangs[0] + viewPersubSet * osNum, sizeof(float) * res);
	memcpy(shzPos[osNum-1] + viewPersubSet, &hzPos[0] + viewPersubSet * osNum, sizeof(float) * res);
	memcpy(shprj[osNum-1] + viewPersubSet * DDNPV, &hprj[0] + viewPersubSet * osNum * DDNPV, sizeof(float) * DDNPV * res);

	// Generate the row sum and col sum weights
	thrust::host_vector<float> allOneproj(totJN, 1.0);
	for(int i = 0; i != osNum; ++i)
	{
		DD3Proj_gpu(x0, y0, z0, DNU, DNV,
				&xds[0], &yds[0], &zds[0],
			imgXCenter, imgYCenter, imgZCenter, shangs[i], shzPos[i], SPN[i],
			XN, YN, ZN, &oneVol[0], srowSum[i], dx, dz, &mask[0], 0, forwardMethod);
		DD3Back_gpu(x0, y0, z0, DNU, DNV,
				&xds[0],&yds[0],&zds[0],
			imgXCenter, imgYCenter, imgZCenter, shangs[i], shzPos[i], SPN[i],
			XN, YN, ZN, colSum[i], &allOneproj[0], dx, dz, &mask[0], 0, 0, backwardMethod);
		std::cout<<"os "<<i<<std::endl;
	}
	allOneproj.clear();


	thrust::host_vector<float> lasImg(totVN, 0);

	float t1 = 1.0;
	float t2 = 1.0;
	int totN = static_cast<int>(totVN);
	clock_t start = clock();
	for (int ii = 0; ii != iterNum; ++ii)
	{
		thrust::copy(hvol.begin(),hvol.end(),lasImg.begin());
		for(int subIdx = 0; subIdx != osNum;++subIdx)
		{
			DD3Proj_gpu(x0, y0, z0, DNU, DNV,
					&xds[0],&yds[0],&zds[0],
				imgXCenter, imgYCenter, imgZCenter, shangs[subIdx], shzPos[subIdx], SPN[subIdx],
				XN, YN, ZN, &hvol[0], onePrj[subIdx], dx, dz, &mask[0], 0, forwardMethod);
			SART_Weighting<float>::Proj_ptr(onePrj[subIdx], shprj[subIdx], srowSum[subIdx], DNU*DNV*SPN[subIdx]);

			DD3Back_gpu(x0, y0, z0, DNU, DNV,
					&xds[0],&yds[0],&zds[0],
				imgXCenter, imgYCenter, imgZCenter, shangs[subIdx], shzPos[subIdx], SPN[subIdx],
				XN, YN, ZN, &oneVol[0], onePrj[subIdx], dx, dz, &mask[0], 0, 0, backwardMethod);
			SART_Weighting<float>::Back_ptr(&oneVol[0], &hvol[0], colSum[subIdx], totN);

			std::cout<<"OS is "<<subIdx<<std::endl;

		}

		t2 = (1.0f + sqrtf(1.0f + 4.0f * t1 * t1)) / 2.0f;
		FISTA<float>(lasImg, hvol, t1, t2, -2000.0f, 6000.0f);
		t1 = t2;

		if(outputMedRes && (ii%10 == 0))
		{
			std::stringstream ss;
			ss<<ii;
			std::string name = VolumeName + "_" + ss.str();
			std::ofstream fu(name.c_str(), std::ios::binary);
			fu.write((char*) &hvol[0], sizeof(float) * totVN);
			fu.close();
		}
		std::cout << "Iteration # = "<<ii<<std::endl;
	}
	clock_t end = clock();

	std::string name = VolumeName + ".raw";
	std::ofstream fu(name.c_str(), std::ios::binary);
	fu.write((char*) &hvol[0], sizeof(float) * totVN);
	fu.close();

	std::cout<<"The total time is "<<static_cast<double>(end - start)/static_cast<double>(CLOCKS_PER_SEC);
	lasImg.clear();
	xds.clear();
	yds.clear();
	zds.clear();

	for(int i = 0; i != osNum; ++i)
	{
		delete[] shangs[i];
		delete[] shzPos[i];
		delete[] shprj[i];
		delete[] srowSum[i];
		delete[] onePrj[i];
		delete[] colSum[i];
	}
	delete[] shangs;
	delete[] shzPos;
	delete[] shprj;
	delete[] srowSum;
	delete[] onePrj;
	delete[] colSum;

}


template<typename T>
struct CGop
{
	T alpha;
	CGop(T a):alpha(a){}
	__host__ __device__ T operator()(T x, T y)
	{
		return x + alpha * y;
	}
};


template<typename T>
double calculateAlpha(const thrust::host_vector<T>& R, const thrust::host_vector<T>& MDelta, const thrust::host_vector<T>& Delta)
{
	// Calculate alpha
	// TODO: Use only one transform to generate the result
	double alpha1 = thrust::inner_product(R.begin(),R.end(),R.begin(),0.0);
	double alpha2 = thrust::inner_product(MDelta.begin(),MDelta.end(),Delta.begin(),0.0);
	return alpha1 / alpha2;

//	//Float datatype overflow
//	typedef thrust::tuple<double, double> OutputType;
//	OutputType res(thrust::make_tuple<double,double>(0,0));
//
//	struct Product
//	{
//		typedef thrust::tuple<T, T, T> InputType;
//		typedef thrust::tuple<double, double> OutputType;
//		__host__ __device__ OutputType operator()(InputType input)
//		{
//			double R = thrust::get<0>(input);
//			double MDelta = thrust::get<1>(input);
//			double Delta = thrust::get<2>(input);
//			return thrust::make_tuple(R*R,MDelta * Delta);
//		}
//	};
//
//
//	struct Summ
//	{
//		typedef thrust::tuple<double, double> InputType;
//		typedef thrust::tuple<double, double> OutputType;
//		__host__ __device__ OutputType operator()(InputType left, InputType right)
//		{
//			double A1 = thrust::get<0>(left);
//			double B1 = thrust::get<1>(left);
//			double A2 = thrust::get<0>(right);
//			double B2 = thrust::get<1>(right);
//			return thrust::make_tuple(A1 + A2, B1 + B2);
//		}
//	};
//
//
//	//return aTa / bTC;
//	thrust::transform_reduce(
//			thrust::make_zip_iterator(thrust::make_tuple(R.begin(),MDelta.begin(),Delta.begin())),
//			thrust::make_zip_iterator(thrust::make_tuple(R.end(),MDelta.end(),Delta.end())),
//			Product(),
//			res,
//			Summ());
//	return thrust::get<0>(res) / thrust::get<1>(res);

}

template<typename T>
double calculateBeta(
		const thrust::host_vector<T>& X,
		const thrust::host_vector<T>& R,
		const thrust::host_vector<T>& nextR,
		const thrust::host_vector<T>& Delta,
		const CGMethod cgm = CG_FR)
{
	double beta1;
	double beta2;
	double beta;

	thrust::host_vector<T> betaTempVar(nextR.size(),0);

	//Calculate beta;
	switch(cgm)
	{
	case CG_FR:
		// TODO: USE ONE TRANSORM TO UPDATE IT (NOTE: This will not reduce the complexity of the calculation)
		beta1 = thrust::inner_product(nextR.begin(),nextR.end(),nextR.begin(),0.0);
		beta2 = thrust::inner_product(R.begin(),R.end(),R.begin(),0.0);
		beta = beta1 / beta2;
		break;
	case CG_PRP:
		// TODO: USE ONLY ONE TRANSFORM TO UPDATE IT (NOTE: This will not reduce the complexity of the calculation)
		thrust::transform(nextR.begin(),nextR.end(),R.begin(),betaTempVar.begin(),thrust::minus<float>());
		beta1 = thrust::inner_product(nextR.begin(),nextR.end(),betaTempVar.begin(),0.0);
		beta2 = thrust::inner_product(R.begin(),R.end(),R.begin(),0.0);
		beta = beta1 / beta2;
		break;
	case CG_CW:
		// TODO: USE ONLY ONE TRANSFORM TO UPDATE IT (NOTE: This will not reduce the complexity of the calculation)
		thrust::transform(nextR.begin(),nextR.end(),R.begin(),betaTempVar.begin(),thrust::minus<float>());\
		beta1 = thrust::inner_product(nextR.begin(), nextR.end(),betaTempVar.begin(),0.0);
		beta2 = thrust::inner_product(Delta.begin(), Delta.end(),betaTempVar.begin(),0.0);
		beta = beta1 / beta2;
		break;
	case CG_DI:
		// TODO: USE ONLY ONE TRANSFORM TO UPDATE IT (NOTE: This will not reduce the complexity of the calculation)
		beta1 = thrust::inner_product(nextR.begin(),nextR.end(),nextR.begin(),0.0);
		beta2 = thrust::inner_product(Delta.begin(),Delta.end(),R.begin(),0.0);
		beta = -beta1 / beta2;
		break;
	case CG_DY:
		// TODO: USE ONLY ONE TRANSFORM TO UPDATE IT (NOTE: This will not reduce the complexity of the calculation)
		beta1 = thrust::inner_product(nextR.begin(),nextR.end(),nextR.begin(),0.0);
		thrust::transform(nextR.begin(),nextR.end(),R.begin(),betaTempVar.begin(),thrust::minus<float>());
		beta2 = thrust::inner_product(Delta.begin(), Delta.end(),betaTempVar.begin(),0.0);
		beta = beta1 / beta2;
		break;
	default:
		// TODO: USE ONLY ONE TRANSFORM TO UPDATE IT (NOTE: This will not reduce the complexity of the calculation)
		beta1 = thrust::inner_product(nextR.begin(),nextR.end(),nextR.begin(),0.0);
		beta2 = thrust::inner_product(R.begin(),R.end(),R.begin(),0.0);
		beta = beta1 / beta2;
		break;
	}

	betaTempVar.clear();
	return beta;
}



void updateXandR(
		thrust::host_vector<float>& X,
		thrust::host_vector<float>& nextR,
		thrust::host_vector<float>& R,
		thrust::host_vector<float>& MDelta,
		thrust::host_vector<float>& Delta,
		const float apa)
{
	struct TransForm
	{
		float alpha;

		typedef thrust::tuple<float,float,float,float> InputDataType;
		typedef thrust::tuple<float,float> OutputDataType;

		TransForm(const float alp):alpha(alp){}

		__host__ __device__ OutputDataType operator()(InputDataType data)
		{
			float x = thrust::get<0>(data);
			float r = thrust::get<1>(data);
			float md = thrust::get<2>(data);
			float d = thrust::get<3>(data);

			float n1 = x + alpha * d;
			float n2 = r - alpha * md;
			return thrust::make_tuple<float,float>(n1,n2);
		}
	};


	thrust::transform(
			thrust::make_zip_iterator(thrust::make_tuple(X.begin(),R.begin(),MDelta.begin(),Delta.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(X.end(),R.end(),MDelta.end(),Delta.end())),
			thrust::make_zip_iterator(thrust::make_tuple(X.begin(),nextR.begin())),
			TransForm(apa));
}



// The CG algorithm with one GPU (Titan X)
void CG(thrust::host_vector<float>& reconImg, // The image volume to be reconstructed
		thrust::host_vector<float>& initImg, // The initial image volume
		thrust::host_vector<float>& hprj, // Projection data
		thrust::host_vector<float>& hangs, // projection views
		thrust::host_vector<float>& hzPos, // source positions in different views
		thrust::host_vector<byte>& mask, // image mask along in-plane direction
		const std::string& VolumeName, // The reconstructed volume name to be stored (Z, X, Y) order in float datatype
		const int osNum, // # of OS
		const int iterNum, // # of total iterations
		const bool outputMedRes, // output the intermediate results
		const float sid, // source to iso-center distance
		const float sdd, // source to detector distance
		const int DNU, // # of detector elements along in-plane direction
		const int DNV, // # of detector elements along bench moving direction
		const int PN, // # total view #
		const float imgXCenter, // volume center coordinate x
		const float imgYCenter, // volume center coordinate y
		const float imgZCenter, // volume center coordinate z
		const int XN, // # of volume pixel along x
		const int YN, // # of volume pixel along y
		const int ZN, // # of volume pixel along z
		const float dx, // image pixel size along in-plane direction
		const float dz, // slice thickness
		const float col_size, // detector elements size along in-plane direction
		const float row_size, // detector elements size along bench moving direction
		const float col_offset, // #(rational) of detector elements that offset from the center of the iso-center along in-plane direction
		const float row_offset, // #(rational) of detector elements that offset from the center of the iso-center along bench moving direction
		const ForwardDDMethod forwardMethod, // Forward projection model
		const BackwardDDMethod backwardMethod,
		const CGMethod cgm) // Back projection model
{
	// Define the parameters
	float x0 = 0.0f;
	float y0 = sid;
	float z0 = 0.0f;

	float* xds = new float[DNU];
	float* yds = new float[DNU];
	float* zds = new float[DNV];


	//Generate the positions of the detectors
	const float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
	float curBeta = 0;
	for (int ii = 0; ii != DNU; ++ii)
	{
		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
		xds[ii] = sinf(curBeta) * sdd;
		yds[ii] = sid - cosf(curBeta) * sdd;
	}

	for (int ii = 0; ii != DNV; ++ii)
	{
		zds[ii] = (ii - (DNV - 1.0) * 0.5 + row_offset) * row_size;
	}

	const int TOTVOLN(XN * YN * ZN);
	const int TOTPRJN(DNU * DNV * PN);

	// Calculate B
	typedef thrust::host_vector<float> hostvec;



	hostvec X(TOTVOLN, 0);
	hostvec B(TOTVOLN, 0);
	hostvec R(TOTVOLN, 0);
	hostvec nextR(TOTVOLN, 0);
	hostvec Delta(TOTVOLN, 0);
	hostvec TDelta(TOTPRJN, 0);
	hostvec MDelta(TOTVOLN, 0);


	DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
				imgXCenter, imgYCenter, imgZCenter, &hangs[0], &hzPos[0], PN,
				XN, YN, ZN, &B[0], &hprj[0], dx, dz, &mask[0], 0, 0, backwardMethod);

	DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
				imgXCenter, imgYCenter, imgZCenter, &hangs[0], &hzPos[0], PN,
				XN, YN, ZN, &initImg[0], &TDelta[0], dx, dz, &mask[0], 0, forwardMethod);

	DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
					imgXCenter, imgYCenter, imgZCenter, &hangs[0], &hzPos[0], PN,
					XN, YN, ZN, &R[0], &TDelta[0], dx, dz, &mask[0], 0, 0, backwardMethod);

	thrust::transform(B.begin(),B.end(),R.begin(),R.begin(),thrust::minus<float>()); //R = B - R; // Residital

	Delta = R;

	double alpha(1.0);
	double beta(1.0);

//	double beta1(1.0);
//	double beta2(1.0);
	//thrust::host_vector<float> betaTempVar(nextR.size(), 0);

	for(int i = 0; i != iterNum; ++i)
	{
		// Calculate the scalar alpha
		// (1) calculate MDelta
		DD3Proj_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
			imgXCenter, imgYCenter, imgZCenter, &hangs[0], &hzPos[0], PN,
			XN, YN, ZN, &Delta[0], &TDelta[0], dx, dz, &mask[0], 0, forwardMethod);
		DD3Back_gpu(x0, y0, z0, DNU, DNV, xds, yds, zds,
					imgXCenter, imgYCenter, imgZCenter, &hangs[0], &hzPos[0], PN,
					XN, YN, ZN, &MDelta[0], &TDelta[0], dx, dz, &mask[0], 0, 0, backwardMethod);

		// (2) Calculate alpha
		// TODO: Optimize this algorithm with only one loop
		alpha = calculateAlpha<float>(R, MDelta, Delta);

		updateXandR(X,nextR,R,MDelta,Delta,alpha);

		// (3) Calculate beta
		beta = calculateBeta<float>(X, R, nextR, Delta, cgm);
//
//		// (3) Calculate beta
//		switch(cgm)
//		{
//		case CG_FR:
//			// TODO: USE ONE TRANSORM TO UPDATE IT (NOTE: This will not reduce the complexity of the calculation)
//			beta1 = thrust::inner_product(nextR.begin(),nextR.end(),nextR.begin(),0.0);
//			beta2 = thrust::inner_product(R.begin(),R.end(),R.begin(),0.0);
//			beta = beta1 / beta2;
//			break;
//		case CG_PRP:
//			// TODO: USE ONLY ONE TRANSFORM TO UPDATE IT (NOTE: This will not reduce the complexity of the calculation)
//			thrust::transform(nextR.begin(),nextR.end(),R.begin(),betaTempVar.begin(),thrust::minus<float>());
//			beta1 = thrust::inner_product(nextR.begin(),nextR.end(),betaTempVar.begin(),0.0);
//			beta2 = thrust::inner_product(R.begin(),R.end(),R.begin(),0.0);
//			beta = beta1 / beta2;
//			break;
//		case CG_CW:
//			// TODO: USE ONLY ONE TRANSFORM TO UPDATE IT (NOTE: This will not reduce the complexity of the calculation)
//			thrust::transform(nextR.begin(),nextR.end(),R.begin(),betaTempVar.begin(),thrust::minus<float>());\
//			beta1 = thrust::inner_product(nextR.begin(), nextR.end(),betaTempVar.begin(),0.0);
//			beta2 = thrust::inner_product(Delta.begin(), Delta.end(),betaTempVar.begin(),0.0);
//			beta = beta1 / beta2;
//			break;
//		case CG_DI:
//			// TODO: USE ONLY ONE TRANSFORM TO UPDATE IT (NOTE: This will not reduce the complexity of the calculation)
//			beta1 = thrust::inner_product(nextR.begin(),nextR.end(),nextR.begin(),0.0);
//			beta2 = thrust::inner_product(Delta.begin(),Delta.end(),R.begin(),0.0);
//			beta = -beta1 / beta2;
//			break;
//		case CG_DY:
//			// TODO: USE ONLY ONE TRANSFORM TO UPDATE IT (NOTE: This will not reduce the complexity of the calculation)
//			beta1 = thrust::inner_product(nextR.begin(),nextR.end(),nextR.begin(),0.0);
//			thrust::transform(nextR.begin(),nextR.end(),R.begin(),betaTempVar.begin(),thrust::minus<float>());
//			beta2 = thrust::inner_product(Delta.begin(), Delta.end(),betaTempVar.begin(),0.0);
//			beta = beta1 / beta2;
//			break;
//		default:
//			// TODO: USE ONLY ONE TRANSFORM TO UPDATE IT (NOTE: This will not reduce the complexity of the calculation)
//			beta1 = thrust::inner_product(nextR.begin(),nextR.end(),nextR.begin(),0.0);
//			beta2 = thrust::inner_product(R.begin(),R.end(),R.begin(),0.0);
//			beta = beta1 / beta2;
//			break;
//		}

		// Update current Delta with nextR;
		thrust::transform(nextR.begin(),nextR.end(),Delta.begin(),Delta.begin(),CGop<float>(beta));
		R = nextR;

		if(outputMedRes)
		{
			std::stringstream ss;
			ss<<i;
			std::string name = VolumeName + "_" + ss.str() + ".raw";
			std::ofstream fu(name.c_str(), std::ios::binary);
			fu.write((char*) &X[0], sizeof(float) * XN * YN * ZN);
			fu.close();
		}

	}

	std::ofstream fou3(VolumeName.c_str(), std::ios::binary);
	fou3.write((char*) &X[0], sizeof(float) * XN * YN * ZN);
	fou3.close();

	delete[] xds;
	delete[] yds;
	delete[] zds;

	B.clear();
	R.clear();
	nextR.clear();
	Delta.clear();
	TDelta.clear();
	MDelta.clear();
	X.clear();
}


