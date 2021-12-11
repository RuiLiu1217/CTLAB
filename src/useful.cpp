//#include "useful.hpp"
//#include "utilities.hpp"
//#include "cudaCheckReturner.h"
///// \brief the rotation of the 2D vector according to cos(T) and sin(T)
///// \param p original vector
///// \param cosT cos(T)
///// \param sinT sin(T)
//inline __host__	__device__ float2 rotation(const float2& p, const float& cosT, const float& sinT)
//{
//	float2 curP;
//	curP.x = p.x * cosT - p.y * sinT;
//	curP.y = p.x * sinT + p.y * cosT;
//	return curP;
//}
///// \brief the rotation of the 3D vector according to cos(T) and sin(T)
///// \param p original vector
///// \param cosT cos(T)
///// \param sinT sin(T)
//inline __host__ __device__ float3 rotation(const float3& p, const float& cosT, const float& sinT)
//{
//	float3 curP;
//	curP.x = p.x * cosT - p.y * sinT;
//	curP.y = p.x * sinT + p.y * cosT;
//	curP.z = p.z;
//	return curP;
//}
//
//void generatingMatrix(
//	std::vector<int>& rowIdx,
//	std::vector<int>& colIdx,
//	std::vector<float>& weight,
//	std::vector<int>& times,
//	const bool angDir,
//	const FanEAGeo& FanGeo,
//	const Image& Img,
//	const unsigned int& sampNum)
//{
//	float2 MINO = make_float2(
//		-Img.m_Size.x / 2.0f + Img.m_Bias.x,
//		-Img.m_Size.y / 2.0f + Img.m_Bias.y);
//
//	float curAng = 0;
//	float cosT = 0;
//	float sinT = 0;
//	Ray2D ray;
//
//	//unsigned int detId;
//	float ang(0); //bias angle from the center of the Fan Beams
//
//	float2 curDetPos; //the current detector element position;
//	//float totWeight(0);
//
//	float smallDetStep = FanGeo.m_DetStp / sampNum; //ÏÂ²ÉÑùºóµÄdetector²œ³€;
//	float cntSmallDet = sampNum * 0.5f;
//	float realAng = 0;
//	unsigned int angIdx = 0;
//	unsigned int detIdx = 0;
//	unsigned int subDetIdx = 0;
//
//	int rowidx = 0;
//	for (angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
//	{
//		//Current rotation angle;
//		curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		cosT = cosf(curAng);
//		sinT = sinf(curAng);
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//		for (detIdx = 0; detIdx != FanGeo.m_DetN; ++detIdx)
//		{
//
//			rowidx = angIdx * FanGeo.m_DetN + detIdx;
//
//			//Calculate current Angle
//			if (angDir)
//			{
//				ang = ((float) detIdx - FanGeo.m_DetCntIdx + 0.5f) * FanGeo.m_DetStp;
//			}
//			else
//			{
//				ang = ((FanGeo.m_DetN - 1 - (float) detIdx) - FanGeo.m_DetCntIdx + 0.5f) * FanGeo.m_DetStp;
//			}
//
//			for (subDetIdx = 0; subDetIdx != sampNum; ++subDetIdx)
//			{
//				//correct ang
//				realAng = ang + (static_cast<float>(subDetIdx) -cntSmallDet + 0.5f) *smallDetStep;
//				// current detector element position;
//				curDetPos = rotation(make_float2(sinf(realAng) * FanGeo.m_S2D, -cosf(realAng) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//
//				// X-ray direction
//				ray.d = normalize(curDetPos - ray.o);
//
//				pushMatrix(ray.o.x, ray.o.y,
//					curDetPos.x, curDetPos.y,
//					MINO.x, MINO.y,
//					Img.m_Step.x, Img.m_Step.y,
//					Img.m_Reso.x, Img.m_Reso.y,
//					rowidx, rowIdx, colIdx, weight);
//			}
//		}
//	}
//}
//
//
//
//void generatingMatrix(
//	std::vector<int>& rowIdx,
//	std::vector<int>& colIdx,
//	std::vector<float>& weight,
//	std::vector<int>& times,
//	const bool angDir,
//	const FanEDGeo& FanGeo,
//	const Image& Img,
//	const unsigned int& sampNum)
//{
//	float2 MINO = make_float2(
//		-Img.m_Size.x / 2.0f + Img.m_Bias.x,
//		-Img.m_Size.y / 2.0f + Img.m_Bias.y);
//
//	float curAng = 0;
//	float cosT = 0;
//	float sinT = 0;
//	Ray2D ray;
//
//	//unsigned int detId;
//	float ang(0); //bias angle from the center of the Fan Beams
//
//	float2 curDetPos; //the current detector element position;
//	//float totWeight(0);
//
//	float smallDetStep = FanGeo.m_DetStp / sampNum; //ÏÂ²ÉÑùºóµÄdetector²œ³€;
//	float cntSmallDet = sampNum * 0.5f;
//	float realAng = 0;
//	unsigned int angIdx = 0;
//	unsigned int detIdx = 0;
//	unsigned int subDetIdx = 0;
//
//	int rowidx = 0;
//	for (angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
//	{
//		//Current rotation angle;
//		curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		cosT = cosf(curAng);
//		sinT = sinf(curAng);
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//		for (detIdx = 0; detIdx != FanGeo.m_DetN; ++detIdx)
//		{
//
//			rowidx = angIdx * FanGeo.m_DetN + detIdx;
//
//			//Calculate current Angle
//			if (angDir)
//			{
//				ang = ((float) detIdx - FanGeo.m_DetCntIdx + 0.5f) * FanGeo.m_DetStp;
//			}
//			else
//			{
//				ang = ((FanGeo.m_DetN - 1 - (float) detIdx) - FanGeo.m_DetCntIdx + 0.5f) * FanGeo.m_DetStp;
//			}
//
//			for (subDetIdx = 0; subDetIdx != sampNum; ++subDetIdx)
//			{
//				//correct ang
//				realAng = ang + (static_cast<float>(subDetIdx) -cntSmallDet + 0.5f) *smallDetStep;
//				// current detector element position;
//				//curDetPos = rotation(make_float2(sinf(realAng) * FanGeo.m_S2D, -cosf(realAng) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//				curDetPos = rotation(make_float2(realAng, -FanGeo.m_O2D), cosT, sinT);
//				// X-ray direction
//				ray.d = normalize(curDetPos - ray.o);
//
//				pushMatrix(ray.o.x, ray.o.y,
//					curDetPos.x, curDetPos.y,
//					MINO.x, MINO.y,
//					Img.m_Step.x, Img.m_Step.y,
//					Img.m_Reso.x, Img.m_Reso.y,
//					rowidx, rowIdx, colIdx, weight);
//			}
//		}
//	}
//}
//
//// Generate the projection matrix in COO format with the given
//// geometry configuration
//void GenerateProjectionMatrix(const FanEAGeo& FanGeo, const Image& Img)
//{
//	std::vector<int> rowIdx;
//	std::vector<int> colIdx;
//	std::vector<float> weight;
//	std::vector<int> times;
//
//	generatingMatrix(rowIdx, colIdx, weight, times, true, FanGeo, Img, 1);
//
//	size_t sizeOfVec = weight.size();
//
//	std::stringstream ss;
//	ss << sizeOfVec;
//	std::string rowIdxFile = "rowIdx" + ss.str() + ".max";
//	std::string colIdxFile = "colIdx" + ss.str() + ".max";
//	std::string wegFile = "weight" + ss.str() + ".max";
//	std::string timFile = "times" + ss.str() + ".max";
//	std::string confFile = "config_for_" + ss.str() + ".txt";
//
//	std::ofstream fid1(rowIdxFile.c_str(), std::ios::binary);
//	std::ofstream fid2(colIdxFile.c_str(), std::ios::binary);
//	std::ofstream fid3(wegFile.c_str(), std::ios::binary);
//	std::ofstream fid4(timFile.c_str(), std::ios::binary);
//	std::ofstream fid5(confFile.c_str());
//	if (!(fid5.is_open() && fid1.is_open() && fid2.is_open() && fid3.is_open() && fid4.is_open()))
//	{
//		std::cout << "Cannot open the file\n";
//	}
//	fid1.write((const char*) (&rowIdx[0]), sizeof(int) * rowIdx.size());
//	fid2.write((const char*) (&colIdx[0]), sizeof(int) * colIdx.size());
//	fid3.write((const char*) (&weight[0]), sizeof(float) * weight.size());
//	fid4.write((const char*) (&times[0]), sizeof(int) * times.size());
//	fid5 << "Image size " << Img.m_Reso.x << "x" << Img.m_Reso.y << "\n" <<
//		"Detector reso " << FanGeo.m_DetN << ";\n ViewNum " << FanGeo.m_ViwN << "\n";
//	fid1.close();
//	fid2.close();
//	fid3.close();
//	fid4.close();
//	fid5.close();
//}
//
//
//
//
//
//
//
//
//
//
//
//
////
////
////void GenerateProjectionMatrix(const FanEDGeo& FanGeo, const Image& Img)
////{
////	std::vector<int> rowIdx;
////	std::vector<int> colIdx;
////	std::vector<float> weight;
////	std::vector<int> times;
////
////	generatingMatrix(rowIdx, colIdx, weight, times, true, FanGeo, Img, 1);
////
////	size_t sizeOfVec = weight.size();
////
////	std::stringstream ss;
////	ss << sizeOfVec;
////	std::string rowIdxFile = "rowIdx" + ss.str() + ".max";
////	std::string colIdxFile = "colIdx" + ss.str() + ".max";
////	std::string wegFile = "weight" + ss.str() + ".max";
////	std::string timFile = "times" + ss.str() + ".max";
////	std::string confFile = "config_for_" + ss.str() + ".txt";
////
////	std::ofstream fid1(rowIdxFile.c_str(), std::ios::binary);
////	std::ofstream fid2(colIdxFile.c_str(), std::ios::binary);
////	std::ofstream fid3(wegFile.c_str(), std::ios::binary);
////	std::ofstream fid4(timFile.c_str(), std::ios::binary);
////	std::ofstream fid5(confFile.c_str());
////	if (!(fid5.is_open() && fid1.is_open() && fid2.is_open() && fid3.is_open() && fid4.is_open()))
////	{
////		std::cout << "Cannot open the file\n";
////	}
////	fid1.write((const char*) (&rowIdx[0]), sizeof(int) * rowIdx.size());
////	fid2.write((const char*) (&colIdx[0]), sizeof(int) * colIdx.size());
////	fid3.write((const char*) (&weight[0]), sizeof(float) * weight.size());
////	fid4.write((const char*) (&times[0]), sizeof(int) * times.size());
////	fid5 << "Image size " << Img.m_Reso.x << "x" << Img.m_Reso.y << "\n" <<
////		"Detector reso " << FanGeo.m_DetN << ";\n ViewNum " << FanGeo.m_ViwN << "\n";
////	fid1.close();
////	fid2.close();
////	fid3.close();
////	fid4.close();
////	fid5.close();
////}
//
//
//void DEMO6()
//{
//	//FanEAGeo FanGeo;
//	//Image Img;
//
//	//FanGeo.m_DetArc = 0.95928517242269f;
//	//FanGeo.m_DetN = 1024;
//	//FanGeo.m_DetCntIdx = FanGeo.m_DetN * 0.5f;
//	//FanGeo.m_DetStp = FanGeo.m_DetArc / FanGeo.m_DetN;
//	//FanGeo.m_O2D = 4.082259521484375e+02;
//
//	//FanGeo.m_S2O = 5.385200195312500e+02;
//	//FanGeo.m_S2D = FanGeo.m_S2O + FanGeo.m_O2D;
//	//FanGeo.m_ViwBeg = 0.0f;// (-2.668082275390625e+02 / 180.0* 3.14159265358979323846264);
//	//FanGeo.m_ViwN = 360;
//	//FanGeo.m_ViwStp = 3.14159265358979323846264f * 2.0f / FanGeo.m_ViwN;
//
//	//Img.m_Bias.x = 0.0f;
//	//Img.m_Bias.y = 0.0f;
//	//Img.m_Reso.x = Img.m_Reso.y = 512;
//	//Img.m_Size.x = Img.m_Size.y = static_cast<float>(4.484740011196460e+02);
//	//Img.m_Step.x = Img.m_Size.x / Img.m_Reso.x;
//	//Img.m_Step.y = Img.m_Size.y / Img.m_Reso.y;
//
//	//GenerateProjectionMatrix(FanGeo, Img);
//	//std::cout << "Demo 6 finished\n";
//}
//
//
//
//
////
////
////template<typename T>
////void GenerateProjectionMatrix_AIM_temp(const FanEAGeo& FanGeo, const Image& Img, const std::string& FileName)
////{
////	std::vector<T>coeff;
////	std::vector<int> rowIdx;
////	std::vector<int> colIdx;
////
////	int detIdx(0);
////	int angIdx(0);
////	int imgLIdx(0);
////	int imgWIdx(0);
////	int currowIdx(0);
////	int curcolIdx(0);
////
////	const T cntImgX = (static_cast<T>(Img.m_Reso.x) - 1) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
////	const T cntImgY = (static_cast<T>(Img.m_Reso.y) - 1) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
////	const T area = Img.m_Step.x * Img.m_Step.y;
////	T curAng = 0;
////	T cosT(0), sinT(0);
////	T sour[2] = { 0, 0 };
////	T SVA[3], SVB[3];
////	T grid[4][3];
////	T coef = 0;
////	T pdist = 0;
////	T dir[2], initDir[2];
////	T len = 0;
////
////	for (angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
////	{
////		curAng = FanGeo.m_ViwBeg + angIdx* FanGeo.m_ViwStp;
////		while (curAng < 0){ curAng += TWOPI; }
////		while (curAng > TWOPI){ curAng -= (TWOPI); }
////		cosT = cos(curAng);
////		sinT = sin(curAng);
////		sour[0] = -FanGeo.m_S2O * sinT;
////		sour[1] = FanGeo.m_S2O * cosT;
////
////		for (detIdx = 0; detIdx != FanGeo.m_DetN; ++detIdx)
////		{
////			currowIdx = angIdx * FanGeo.m_DetN + detIdx;
////			calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, detIdx);
////			for (imgWIdx = 0; imgWIdx != Img.m_Reso.y; ++imgWIdx)
////			{
////				for (imgLIdx = 0; imgLIdx != Img.m_Reso.x; ++imgLIdx)
////				{
////					curcolIdx = imgWIdx *Img.m_Reso.x + imgLIdx;
////					grid[0][0] = (imgLIdx - cntImgX - 0.5) * Img.m_Step.x;
////					grid[0][1] = (imgWIdx - cntImgY - 0.5) * Img.m_Step.y;
////
////					dir[0] = grid[0][0] - sour[0];
////					dir[1] = grid[0][1] - sour[1];
////					len = hypot(dir[0], dir[1]);
////					dir[0] /= len;
////					dir[1] /= len;
////					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
////					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
////					initDir[0] = dir[0] * cosT + dir[1] * sinT;
////					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
////					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
////
////
////					grid[1][0] = grid[0][0] + Img.m_Step.x;
////					grid[1][1] = grid[0][1];
////					dir[0] = grid[1][0] - sour[0];
////					dir[1] = grid[1][1] - sour[1];
////					len = hypot(dir[0], dir[1]);
////					dir[0] /= len;
////					dir[1] /= len;
////					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
////					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
////					initDir[0] = dir[0] * cosT + dir[1] * sinT;
////					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
////					grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
////
////					grid[2][0] = grid[1][0];
////					grid[2][1] = grid[1][1] + Img.m_Step.y;
////					dir[0] = grid[2][0] - sour[0];
////					dir[1] = grid[2][1] - sour[1];
////					len = hypot(dir[0], dir[1]);
////					dir[0] /= len;
////					dir[1] /= len;
////					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
////					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
////					initDir[0] = dir[0] * cosT + dir[1] * sinT;
////					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
////					grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
////
////					grid[3][0] = grid[0][0];
////					grid[3][1] = grid[2][1];
////					dir[0] = grid[3][0] - sour[0];
////					dir[1] = grid[3][1] - sour[1];
////					len = hypot(dir[0], dir[1]);
////					dir[0] /= len;
////					dir[1] /= len;
////					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
////					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
////					initDir[0] = dir[0] * cosT + dir[1] * sinT;
////					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
////					grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
////
////					//Sort grid
////					SortProjection<T>(grid);
////
////					pdist = hypot((imgLIdx - cntImgX)*Img.m_Step.x - sour[0], (imgWIdx - cntImgY)*Img.m_Step.y - sour[1]);
////					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
////					coef = coef / (pdist * fabs(FanGeo.m_DetStp));
////					if (!IS_ZERO<T>(coef))
////					{
////						coeff.push_back(coef);
////						rowIdx.push_back(currowIdx);
////						colIdx.push_back(curcolIdx);
////					}
////				}
////			}
////		}
////	}
////
////	std::stringstream ss;
////	ss << coeff.size();
////	std::string num(ss.str());
////
////	if (typeid(T) == typeid(float))
////	{
////		std::string coefFile = FileName + num + "f.cof";
////		std::string rowFile = FileName + num + "f.row";
////		std::string colFile = FileName + num + "f.col";
////		std::ofstream fouCoef(coefFile.c_str(), std::ios::binary);
////		fouCoef.write((char*) &(coeff[0]), sizeof(T) *coeff.size());
////		std::ofstream fouCoef2(rowFile.c_str(), std::ios::binary);
////		fouCoef2.write((char*) &(rowIdx[0]), sizeof(int) *rowIdx.size());
////		std::ofstream fouCoef3(colFile.c_str(), std::ios::binary);
////		fouCoef3.write((char*) &(colIdx[0]), sizeof(int) *colIdx.size());
////	}
////	if (typeid(T) == typeid(double))
////	{
////		std::string coefFile = FileName + num + "d.cof";
////		std::string rowFile = FileName + num + "d.row";
////		std::string colFile = FileName + num + "d.col";
////		std::ofstream fouCoef(coefFile.c_str(), std::ios::binary);
////		fouCoef.write((char*) &(coeff[0]), sizeof(T) *coeff.size());
////		std::ofstream fouCoef2(rowFile.c_str(), std::ios::binary);
////		fouCoef2.write((char*) &(rowIdx[0]), sizeof(int) *rowIdx.size());
////		std::ofstream fouCoef3(colFile.c_str(), std::ios::binary);
////		fouCoef3.write((char*) &(colIdx[0]), sizeof(int) *colIdx.size());
////	}
////
////	coeff.clear();
////	rowIdx.clear();
////	colIdx.clear();
////}
////void GenerateProjectionMatrixf_AIM(const FanEAGeo& FanGeo, const Image& Img, const std::string& FileName)
////{
////	GenerateProjectionMatrix_AIM_temp<float>(FanGeo, Img, FileName);
////}
////void GenerateProjectionMatrixd_AIM(const FanEAGeo& FanGeo, const Image& Img, const std::string& FileName)
////{
////	GenerateProjectionMatrix_AIM_temp<double>(FanGeo, Img, FileName);
////}
////
//
//
//template<typename T>
//void GenerateMatrix_OpenMP_AIM_temp(FanEAGeo FanGeo, Image Img, const std::string& FileName)
//{
//	std::vector<T> coeff[32];
//	std::vector<int> colIdx[32];
//	std::vector<int> rowIdx[32];
//
//
//	int threadID;
//
//	//omp_set_num_threads(30);
//#pragma omp parallel for 
//	for (threadID = 0; threadID < 32; threadID++)
//	{
//		int thID = omp_get_thread_num();
//		//std::cout << "current Thread : " << thID << std::endl;
//		int detIdx(0);
//		int angIdx(0);
//		int imgLIdx(0);
//		int imgWIdx(0);
//		int currowIdx(0);
//		int curcolIdx(0);
//		const T cntImgX = (static_cast<T>(Img.m_Reso.x) - 1) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
//		const T cntImgY = (static_cast<T>(Img.m_Reso.y) - 1) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
//		const T area = Img.m_Step.x * Img.m_Step.y;
//		T curAng = 0;
//		T cosT(0), sinT(0);
//		T sour[2] = { 0, 0 };
//		T SVA[3], SVB[3];
//		T grid[4][3];
//		//T summ = 0;
//		T coef = 0;
//		T pdist = 0;
//		T dir[2], initDir[2];
//		T len = 0;
//		int angN = FanGeo.m_ViwN / 30;
//		std::cout << "Thread " << thID << " has " << angN << " angles\n";
//		for (angIdx = 0; angIdx < angN; ++angIdx)
//		{
//			//std::cout << "I am Here " << thID << "\n";
//			curAng = FanGeo.m_ViwBeg + (angIdx + thID * angN)* FanGeo.m_ViwStp;
//			while (curAng < 0){ curAng += TWOPI; }
//			while (curAng > TWOPI){ curAng -= TWOPI; }
//			cosT = cos(curAng);
//			sinT = sin(curAng);
//			sour[0] = -FanGeo.m_S2O * sinT;
//			sour[1] = FanGeo.m_S2O * cosT;
//
//			for (detIdx = 0; detIdx != FanGeo.m_DetN; ++detIdx)
//			{
//				currowIdx = (angIdx + thID * angN) * FanGeo.m_DetN + detIdx;
//				calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, detIdx);
//				for (imgWIdx = 0; imgWIdx != Img.m_Reso.y; ++imgWIdx)
//				{
//					for (imgLIdx = 0; imgLIdx != Img.m_Reso.x; ++imgLIdx)
//					{
//						curcolIdx = imgWIdx *Img.m_Reso.x + imgLIdx;
//						grid[0][0] = (imgLIdx - cntImgX - 0.5) * Img.m_Step.x;
//						grid[0][1] = (imgWIdx - cntImgY - 0.5) * Img.m_Step.y;
//
//						dir[0] = grid[0][0] - sour[0];
//						dir[1] = grid[0][1] - sour[1];
//						len = hypot(dir[0], dir[1]);
//						dir[0] /= len;
//						dir[1] /= len;
//						dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//						dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//						initDir[0] = dir[0] * cosT + dir[1] * sinT;
//						initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//						grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//						grid[1][0] = grid[0][0] + Img.m_Step.x;
//						grid[1][1] = grid[0][1];
//						dir[0] = grid[1][0] - sour[0];
//						dir[1] = grid[1][1] - sour[1];
//						len = hypot(dir[0], dir[1]);
//						dir[0] /= len;
//						dir[1] /= len;
//						dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//						dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//						initDir[0] = dir[0] * cosT + dir[1] * sinT;
//						initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//						grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//						grid[2][0] = grid[1][0];
//						grid[2][1] = grid[1][1] + Img.m_Step.y;
//						dir[0] = grid[2][0] - sour[0];
//						dir[1] = grid[2][1] - sour[1];
//						len = hypot(dir[0], dir[1]);
//						dir[0] /= len;
//						dir[1] /= len;
//						dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//						dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//						initDir[0] = dir[0] * cosT + dir[1] * sinT;
//						initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//						grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//						grid[3][0] = grid[0][0];
//						grid[3][1] = grid[2][1];
//						dir[0] = grid[3][0] - sour[0];
//						dir[1] = grid[3][1] - sour[1];
//						len = hypot(dir[0], dir[1]);
//						dir[0] /= len;
//						dir[1] /= len;
//						dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//						dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//						initDir[0] = dir[0] * cosT + dir[1] * sinT;
//						initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//						grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//						//Sort grid
//						SortProjection<T>(grid);
//
//						pdist = hypot((imgLIdx - cntImgX)*Img.m_Step.x - sour[0], (imgWIdx - cntImgY)*Img.m_Step.y - sour[1]);
//						coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//						coef = coef / (pdist * fabs(FanGeo.m_DetStp));
//						if (!IS_ZERO<T>(coef))
//						{
//							coeff[thID].push_back(coef);
//							rowIdx[thID].push_back(currowIdx);
//							colIdx[thID].push_back(curcolIdx);
//						}
//
//					}
//				}
//			}
//		}
//		std::cout << "current Thread : " << thID << std::endl;
//		//std::cout << "The size of coeff " << thID << "is " << coeff[thID].size() << std::endl;
//		//std::cout << "The size of coeff " << thID << "is " << rowIdx[thID].size() << std::endl;
//		//std::cout << "The size of coeff " << thID << "is " << colIdx[thID].size() << std::endl;
//	}
//
//	size_t totV = 0;
//	for (unsigned int i = 0; i != 32; ++i)
//	{
//		totV += coeff[i].size();
//	}
//	std::stringstream ss;
//	ss << totV;
//	std::string num(ss.str());
//
//	if (typeid(T) == typeid(float))
//	{
//		std::string coefFile = FileName + num + "f.cof";
//		std::string rowFile = FileName + num + "f.row";
//		std::string colFile = FileName + num + "f.col";
//		std::ofstream fouCoef(coefFile.c_str(), std::ios::binary | std::ios::app);
//		std::ofstream fouCoef2(rowFile.c_str(), std::ios::binary | std::ios::app);
//		std::ofstream fouCoef3(colFile.c_str(), std::ios::binary | std::ios::app);
//		for (unsigned int id = 0; id != 32; ++id)
//		{
//			fouCoef.write((char*) &((coeff[id])[0]), sizeof(T) *coeff[id].size());
//			fouCoef2.write((char*) &((rowIdx[id])[0]), sizeof(int) *rowIdx[id].size());
//			fouCoef3.write((char*) &((colIdx[id])[0]), sizeof(int) *colIdx[id].size());
//		}
//	}
//	if (typeid(T) == typeid(double))
//	{
//		std::string coefFile = FileName + num + "d.cof";
//		std::string rowFile = FileName + num + "d.row";
//		std::string colFile = FileName + num + "d.col";
//		std::ofstream fouCoef(coefFile.c_str(), std::ios::binary | std::ios::app);
//		std::ofstream fouCoef2(rowFile.c_str(), std::ios::binary | std::ios::app);
//		std::ofstream fouCoef3(colFile.c_str(), std::ios::binary | std::ios::app);
//		for (unsigned int id = 0; id != 32; ++id)
//		{
//			fouCoef.write((char*) &(coeff[id][0]), sizeof(T) *coeff[id].size());
//			fouCoef2.write((char*) &(rowIdx[id][0]), sizeof(int) *rowIdx[id].size());
//			fouCoef3.write((char*) &(colIdx[id][0]), sizeof(int) *colIdx[id].size());
//		}
//	}
//}
//
//
//void GenerateMatrix_OpenMP_AIM_float(const FanEAGeo& FanGeo, const Image& Img, const std::string& FileName)
//{
//	GenerateMatrix_OpenMP_AIM_temp<float>(FanGeo, Img, FileName);
//}
//
//void GenerateMatrix_OpenMP_AIM_double(const FanEAGeo& FanGeo, const Image& Img, const std::string& FileName)
//{
//	GenerateMatrix_OpenMP_AIM_temp<double>(FanGeo, Img, FileName);
//}
//
//
//void genProjectionMatrix_AIMf(const FanEDGeo FanGeo, const Image Img)
//{
//	genProjectionMatrix_AIM_template<float>(FanGeo, Img);
//}
//void genProjectionMatrix_AIMd(const FanEDGeo FanGeo, const Image Img)
//{
//	genProjectionMatrix_AIM_template<double>(FanGeo, Img);
//}
//
//
//void genProjectionMatrix_AIMf(const FanEAGeo FanGeo, const Image Img)
//{
//	genProjectionMatrix_AIM_template<float>(FanGeo, Img);
//}
//
//void genProjectionMatrix_AIMd(const FanEAGeo FanGeo, const Image Img)
//{
//	genProjectionMatrix_AIM_template<double>(FanGeo, Img);
//}
