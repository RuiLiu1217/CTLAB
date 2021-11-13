#include "useful.hpp"
#include "utilities.hpp"

FanEAGeo::FanEAGeo(void)
{
	m_DetArc = 0.95928517242269f;
	m_DetN = 888;
	m_DetCntIdx = 444.75;
	m_DetStp = m_DetArc / m_DetN;
	m_O2D = 4.082259521484375e+02;

	m_S2O = 5.385200195312500e+02;
	m_S2D = m_S2O + m_O2D;
	m_ViwBeg = 0.0f;// (-2.668082275390625e+02 / 180.0* 3.14159265358979323846264);
	m_ViwN = 360;
	m_ViwStp = 3.14159265358979323846264f * 2.0f / m_ViwN;
}

FanEAGeo::FanEAGeo(const FanEAGeo& rhs)
{
	m_DetArc = rhs.m_DetArc;
	m_DetN = rhs.m_DetN;
	m_DetCntIdx = rhs.m_DetCntIdx;
	m_DetStp = rhs.m_DetStp;
	m_O2D = rhs.m_O2D;
	m_S2O = rhs.m_S2O;
	m_S2D = rhs.m_S2D;
	m_ViwBeg = rhs.m_ViwBeg;
	m_ViwN = rhs.m_ViwN;
	m_ViwStp = rhs.m_ViwStp;
}


FanEAGeo::FanEAGeo(const float S2O, const float O2D, const unsigned int ViwN,
	const float ViwBeg, const float ViwEnd, const float DetArc,
	const unsigned int DetN) :m_S2O(S2O), m_O2D(O2D), m_S2D(S2O + O2D),
	m_ViwN(ViwN), m_ViwBeg(ViwBeg), m_ViwStp((ViwEnd - ViwBeg) / float(ViwN)),
	m_DetArc(DetArc), m_DetN(DetN), m_DetStp(DetArc / DetN),
	m_DetCntIdx(DetN * 0.5f - 0.5f){}


FanEDGeo::FanEDGeo(void)
{
	m_S2O = 85.0f;
	m_O2D = 15.0f;
	m_S2D = 100.0f;

	m_ViwN = 360; // view number
	m_ViwBeg = 0.0f;
	m_ViwStp = 3.14159265358979323846264f * 2.0f / m_ViwN;

	m_DetSize = 30.0f; // Detector size;
	m_DetN = 888; // How many detector cells
	m_DetStp = m_DetSize / m_DetN; // detector cells size;
	m_DetCntIdx = m_DetN * 0.5f - 0.5f; //detector center index;
}

FanEDGeo::FanEDGeo(const FanEDGeo& rhs)
{
	m_S2O = rhs.m_S2O;
	m_O2D = rhs.m_O2D;
	m_S2D = rhs.m_S2D;

	m_ViwN = rhs.m_ViwN;
	m_ViwBeg = rhs.m_ViwBeg;
	m_ViwStp = rhs.m_ViwStp;

	m_DetSize = rhs.m_DetSize;
	m_DetN = rhs.m_DetN;
	m_DetStp = rhs.m_DetStp;
	m_DetCntIdx = rhs.m_DetCntIdx;
}


FanEDGeo::FanEDGeo(const float S2O, const float O2D, const unsigned int ViwN,
	const float ViwBeg, const float ViwEnd, const float DetSize,
	const unsigned int DetN) :m_S2O(S2O), m_O2D(O2D), m_S2D(S2O + O2D),
	m_ViwN(ViwN), m_ViwBeg(ViwBeg), m_ViwStp((ViwEnd - ViwBeg) / ViwN),
	m_DetSize(DetSize), m_DetN(DetN), m_DetStp(DetSize / DetN),
	m_DetCntIdx(DetN * 0.5f - 0.5f){}







ConeEAGeo::ConeEAGeo(void){

	m_DetCntIdx.x = 444.75;

	m_ViwBeg = 0.0f;// (-2.668082275390625e+02 / 180.0* 3.14159265358979323846264);


	m_S2O = 5.385200195312500e+02;
	m_O2D = 4.082259521484375e+02;
	m_S2D = m_S2O + m_O2D;

	m_ViwN = 2200;
	m_ViwBeg = 0.0f;
	m_ViwStp = 3.14159265358979323846264f * 2.0f / m_ViwN;

	m_DetArc = 0.95928517242269f;
	m_DetN = 888;
	m_DetStp = m_DetArc / m_DetN;

	m_DetHei = 64.0f;
	m_DetHN = 64;
	m_DetHStp = m_DetHei / m_DetHN;

	m_DetCntIdx.x = m_DetN * 0.5f - 0.5f;
	m_DetCntIdx.y = m_DetN * 0.5f - 0.5f;
}

ConeEAGeo::ConeEAGeo(const ConeEAGeo& rhs)
{
	m_DetCntIdx = rhs.m_DetCntIdx;
	m_DetArc = rhs.m_DetArc;
	m_DetHei = rhs.m_DetHei;
	m_DetHN = rhs.m_DetHN;
	m_DetHStp = rhs.m_DetHStp;
	m_DetN = rhs.m_DetN;
	m_DetStp = rhs.m_DetStp;
	m_O2D = rhs.m_O2D;
	m_S2D = rhs.m_S2D;
	m_S2O = rhs.m_S2O;
	m_ViwBeg = rhs.m_ViwBeg;
	m_ViwN = rhs.m_ViwN;
	m_ViwStp = rhs.m_ViwStp;
}

ConeEAGeo::ConeEAGeo(const float S2O, const float O2D, const unsigned int viwN,
	const float ViwBeg, const float ViwEnd, const float DetArc, const unsigned int DetN,
	const float DetHei, const unsigned int DetHN) :m_S2O(S2O), m_O2D(O2D),
	m_S2D(S2O + O2D), m_ViwN(viwN), m_ViwBeg(ViwBeg), m_ViwStp((ViwEnd - ViwBeg) / viwN),
	m_DetArc(DetArc), m_DetN(DetN), m_DetStp(DetArc / DetN), m_DetHei(DetHei),
	m_DetHN(DetHN), m_DetHStp(DetHei / DetHN)
{
	m_DetCntIdx.x = m_DetN * 0.5f - 0.5f;
	m_DetCntIdx.y = m_DetHN * 0.5f - 0.5f;
}









ConeEDGeo::ConeEDGeo()
{
	m_S2O = 85.0f;
	m_O2D = 15.0f;
	m_S2D = m_S2O + m_O2D;
	m_ViwN = 360;
	m_ViwBeg = 0.0f;
	m_ViwStp = 3.14159265358979323846264f * 2.0f / m_ViwN;

	m_DetSize = make_float2(34.0f, 17.0f);
	m_DetN = make_int2(1024, 512);
	m_DetStp = make_float2(m_DetSize.x / m_DetN.x, m_DetSize.y / m_DetN.y);
	m_DetCntIdx = make_float2(m_DetN.x * 0.5f - 0.5f, m_DetN.y * 0.5f - 0.5f);
}

ConeEDGeo::ConeEDGeo(const ConeEDGeo& rhs)
{
	m_DetCntIdx = rhs.m_DetCntIdx;
	m_DetN = rhs.m_DetN;
	m_DetSize = rhs.m_DetSize;
	m_DetStp = rhs.m_DetStp;
	m_O2D = rhs.m_O2D;
	m_S2D = rhs.m_S2D;
	m_S2O = rhs.m_S2O;
	m_ViwBeg = rhs.m_ViwBeg;
	m_ViwN = rhs.m_ViwN;
	m_ViwStp = rhs.m_ViwStp;
}


ConeEDGeo::ConeEDGeo(const float S2O, const float O2D, const  int ViwN,
	const float ViwBeg, const float ViwEnd, const float2 DetSize,
	const int2 DetN) :m_S2O(S2O), m_O2D(O2D), m_S2D(S2O + O2D),
	m_ViwN(ViwN), m_ViwBeg(ViwBeg), m_ViwStp((ViwEnd - ViwBeg) / ViwN),
	m_DetSize(DetSize), m_DetN(DetN)
{
	m_DetStp.x = m_DetSize.x / m_DetN.x;
	m_DetStp.y = m_DetSize.y / m_DetN.y;
	m_DetCntIdx.x = m_DetN.x * 0.5f - 0.5f;
	m_DetCntIdx.y = m_DetN.y * 0.5f - 0.5f;
}

Image::Image(void)
{
	m_Bias.x = 0.0f;
	m_Bias.y = 0.0f;  //���ƫ�Ƶĵ�λ����ʵ���?λ;
	m_Reso.x = m_Reso.y = 512;
	m_Size.x = m_Size.y = 400.0;// static_cast<float>(4.484740011196460e+02);
	m_Step.x = m_Size.x / m_Reso.x;
	m_Step.y = m_Size.y / m_Reso.y;
}

Image::Image(const Image& rhs)
{
	m_Bias = rhs.m_Bias;
	m_Reso = rhs.m_Reso;
	m_Size = rhs.m_Size;
	m_Step = rhs.m_Step;
}



Image::Image(
	const unsigned int resoL,
	const unsigned int resoW,
	const float sizeL,
	const float sizeW,
	const float BiasL,
	const float BiasW) :m_Reso(make_uint2(resoL, resoW)),
	m_Size(make_float2(sizeL, sizeW)),
	m_Step(make_float2(sizeL / resoL, sizeW / resoW)),
	m_Bias(make_float2(BiasL, BiasW)){}



Volume::Volume(void)
{
	m_Reso.x = m_Reso.y = 512;
	m_Reso.z = 256;
	m_Size.x = m_Size.y = 20.0f;
	m_Size.z = 10.0f;
	m_Step.x = m_Size.x / m_Reso.x;
	m_Step.y = m_Size.y / m_Reso.y;
	m_Step.z = m_Size.z / m_Reso.z;
	m_Bias = make_float3(0.0f, 0.0f, 0.0f);
}

Volume::Volume(const Volume& rhs)
{
	m_Bias = rhs.m_Bias;
	m_Reso = rhs.m_Reso;
	m_Size = rhs.m_Size;
	m_Step = rhs.m_Step;
}

Volume::Volume(
	const unsigned int resoL,
	const unsigned int resoW,
	const unsigned int resoH,
	const float sizeL,
	const float sizeW,
	const float sizeH,
	const float biasL,
	const float biasW,
	const float biasH) :m_Reso(make_uint3(resoL, resoW, resoH)),
	m_Size(make_float3(sizeL, sizeW, sizeH)),
	m_Step(make_float3(sizeL / resoL, sizeW / resoW, sizeH / resoH)),
	m_Bias(make_float3(biasL, biasW, biasH)){}



//! Transform the number to string
std::string num2String(int& num)
{
	std::stringstream ss;
	ss << num;
	return ss.str();
}
std::string num2String(unsigned int& num)
{
	std::stringstream ss;
	ss << num;
	return ss.str();
}
std::string num2String(const int& num)
{
	std::stringstream ss;
	ss << num;
	return ss.str();
}
std::string num2String(const unsigned int& num)
{
	std::stringstream ss;
	ss << num;
	return ss.str();
}


double3 crossProduct(const double3& a, const double3& b)
{
	double3 res;
	res.x = (a.y*b.z - a.z*b.y);
	res.y = (a.z*b.x - a.x*b.z);
	res.z = (a.x*b.y - a.y*b.x);
	return res;
}
float3 crossProduct(const float3& a, const float3& b)
{
	float3 res;
	res.x = (a.y*b.z - a.z*b.y);
	res.y = (a.z*b.x - a.x*b.z);
	res.z = (a.x*b.y - a.y*b.x);
	return res;
}


void toPolarCoords(const float& x, const float& y, float& r, float& t)
{
	if (IS_ZERO(x)) {
		if (IS_ZERO(y)) {
			r = t = 0;
		}
		else if (y > 0) {
			r = y; t = _PI_2;
		}
		else {
			r = -y; t = -_PI_2;
		}
	}
	else if (IS_ZERO(y)) {
		if (x > 0) {
			r = x; t = 0;
		}
		else {
			r = -x; t = _PI;
		}
	}
	else {
		r = std::hypotf(y, x);
		t = std::atan2(y, x);
	}
}
void toSphericalCoords(const float& x, const float& y, const float& z, float& r, float& s, float& t)
{
	if (IS_ZERO(x) && IS_ZERO(y)) {
		if (IS_ZERO(z)) {
			r = s = t = 0;
		}
		else if (z > 0) {
			r = z;
			s = t = 0;
		}
		else {
			r = -z;
			s = _PI;
			t = 0;
		}
	}
	else {
		toPolarCoords(x, y, r, s);
		if (IS_ZERO(z)) {
			t = _PI_2;
		}
		else {
			r = std::hypotf(z, r);
			t = std::atan2(z, r);
		}
	}
}

void generatingMatrix(
	std::vector<int>& rowIdx,
	std::vector<int>& colIdx,
	std::vector<float>& weight,
	std::vector<int>& times,
	const bool angDir,
	const FanEAGeo& FanGeo,
	const Image& Img,
	const unsigned int& sampNum)
{
	float2 MINO = make_float2(
		-Img.m_Size.x / 2.0f + Img.m_Bias.x,
		-Img.m_Size.y / 2.0f + Img.m_Bias.y);

	float curAng = 0;
	float cosT = 0;
	float sinT = 0;
	Ray2D ray;

	//unsigned int detId;
	float ang(0); //bias angle from the center of the Fan Beams

	float2 curDetPos; //the current detector element position;
	//float totWeight(0);

	float smallDetStep = FanGeo.m_DetStp / sampNum; //ÏÂ²ÉÑùºóµÄdetector²œ³€;
	float cntSmallDet = sampNum * 0.5f;
	float realAng = 0;
	unsigned int angIdx = 0;
	unsigned int detIdx = 0;
	unsigned int subDetIdx = 0;

	int rowidx = 0;
	for (angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
	{
		//Current rotation angle;
		curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
		cosT = cosf(curAng);
		sinT = sinf(curAng);
		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);

		for (detIdx = 0; detIdx != FanGeo.m_DetN; ++detIdx)
		{

			rowidx = angIdx * FanGeo.m_DetN + detIdx;

			//Calculate current Angle
			if (angDir)
			{
				ang = ((float) detIdx - FanGeo.m_DetCntIdx + 0.5f) * FanGeo.m_DetStp;
			}
			else
			{
				ang = ((FanGeo.m_DetN - 1 - (float) detIdx) - FanGeo.m_DetCntIdx + 0.5f) * FanGeo.m_DetStp;
			}

			for (subDetIdx = 0; subDetIdx != sampNum; ++subDetIdx)
			{
				//correct ang
				realAng = ang + (static_cast<float>(subDetIdx) -cntSmallDet + 0.5f) *smallDetStep;
				// current detector element position;
				curDetPos = rotation(make_float2(sinf(realAng) * FanGeo.m_S2D, -cosf(realAng) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);

				// X-ray direction
				ray.d = normalize(curDetPos - ray.o);

				pushMatrix(ray.o.x, ray.o.y,
					curDetPos.x, curDetPos.y,
					MINO.x, MINO.y,
					Img.m_Step.x, Img.m_Step.y,
					Img.m_Reso.x, Img.m_Reso.y,
					rowidx, rowIdx, colIdx, weight);
			}
		}
	}
}



void generatingMatrix(
	std::vector<int>& rowIdx,
	std::vector<int>& colIdx,
	std::vector<float>& weight,
	std::vector<int>& times,
	const bool angDir,
	const FanEDGeo& FanGeo,
	const Image& Img,
	const unsigned int& sampNum)
{
	float2 MINO = make_float2(
		-Img.m_Size.x / 2.0f + Img.m_Bias.x,
		-Img.m_Size.y / 2.0f + Img.m_Bias.y);

	float curAng = 0;
	float cosT = 0;
	float sinT = 0;
	Ray2D ray;

	//unsigned int detId;
	float ang(0); //bias angle from the center of the Fan Beams

	float2 curDetPos; //the current detector element position;
	//float totWeight(0);

	float smallDetStep = FanGeo.m_DetStp / sampNum; //ÏÂ²ÉÑùºóµÄdetector²œ³€;
	float cntSmallDet = sampNum * 0.5f;
	float realAng = 0;
	unsigned int angIdx = 0;
	unsigned int detIdx = 0;
	unsigned int subDetIdx = 0;

	int rowidx = 0;
	for (angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
	{
		//Current rotation angle;
		curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
		cosT = cosf(curAng);
		sinT = sinf(curAng);
		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);

		for (detIdx = 0; detIdx != FanGeo.m_DetN; ++detIdx)
		{

			rowidx = angIdx * FanGeo.m_DetN + detIdx;

			//Calculate current Angle
			if (angDir)
			{
				ang = ((float) detIdx - FanGeo.m_DetCntIdx + 0.5f) * FanGeo.m_DetStp;
			}
			else
			{
				ang = ((FanGeo.m_DetN - 1 - (float) detIdx) - FanGeo.m_DetCntIdx + 0.5f) * FanGeo.m_DetStp;
			}

			for (subDetIdx = 0; subDetIdx != sampNum; ++subDetIdx)
			{
				//correct ang
				realAng = ang + (static_cast<float>(subDetIdx) -cntSmallDet + 0.5f) *smallDetStep;
				// current detector element position;
				//curDetPos = rotation(make_float2(sinf(realAng) * FanGeo.m_S2D, -cosf(realAng) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
				curDetPos = rotation(make_float2(realAng, -FanGeo.m_O2D), cosT, sinT);
				// X-ray direction
				ray.d = normalize(curDetPos - ray.o);

				pushMatrix(ray.o.x, ray.o.y,
					curDetPos.x, curDetPos.y,
					MINO.x, MINO.y,
					Img.m_Step.x, Img.m_Step.y,
					Img.m_Reso.x, Img.m_Reso.y,
					rowidx, rowIdx, colIdx, weight);
			}
		}
	}
}

// Generate the projection matrix in COO format with the given
// geometry configuration
void GenerateProjectionMatrix(const FanEAGeo& FanGeo, const Image& Img)
{
	std::vector<int> rowIdx;
	std::vector<int> colIdx;
	std::vector<float> weight;
	std::vector<int> times;

	generatingMatrix(rowIdx, colIdx, weight, times, true, FanGeo, Img, 1);

	size_t sizeOfVec = weight.size();

	std::stringstream ss;
	ss << sizeOfVec;
	std::string rowIdxFile = "rowIdx" + ss.str() + ".max";
	std::string colIdxFile = "colIdx" + ss.str() + ".max";
	std::string wegFile = "weight" + ss.str() + ".max";
	std::string timFile = "times" + ss.str() + ".max";
	std::string confFile = "config_for_" + ss.str() + ".txt";

	std::ofstream fid1(rowIdxFile.c_str(), std::ios::binary);
	std::ofstream fid2(colIdxFile.c_str(), std::ios::binary);
	std::ofstream fid3(wegFile.c_str(), std::ios::binary);
	std::ofstream fid4(timFile.c_str(), std::ios::binary);
	std::ofstream fid5(confFile.c_str());
	if (!(fid5.is_open() && fid1.is_open() && fid2.is_open() && fid3.is_open() && fid4.is_open()))
	{
		std::cout << "Cannot open the file\n";
	}
	fid1.write((const char*) (&rowIdx[0]), sizeof(int) * rowIdx.size());
	fid2.write((const char*) (&colIdx[0]), sizeof(int) * colIdx.size());
	fid3.write((const char*) (&weight[0]), sizeof(float) * weight.size());
	fid4.write((const char*) (&times[0]), sizeof(int) * times.size());
	fid5 << "Image size " << Img.m_Reso.x << "x" << Img.m_Reso.y << "\n" <<
		"Detector reso " << FanGeo.m_DetN << ";\n ViewNum " << FanGeo.m_ViwN << "\n";
	fid1.close();
	fid2.close();
	fid3.close();
	fid4.close();
	fid5.close();
}














void GenerateProjectionMatrix(const FanEDGeo& FanGeo, const Image& Img)
{
	std::vector<int> rowIdx;
	std::vector<int> colIdx;
	std::vector<float> weight;
	std::vector<int> times;

	generatingMatrix(rowIdx, colIdx, weight, times, true, FanGeo, Img, 1);

	size_t sizeOfVec = weight.size();

	std::stringstream ss;
	ss << sizeOfVec;
	std::string rowIdxFile = "rowIdx" + ss.str() + ".max";
	std::string colIdxFile = "colIdx" + ss.str() + ".max";
	std::string wegFile = "weight" + ss.str() + ".max";
	std::string timFile = "times" + ss.str() + ".max";
	std::string confFile = "config_for_" + ss.str() + ".txt";

	std::ofstream fid1(rowIdxFile.c_str(), std::ios::binary);
	std::ofstream fid2(colIdxFile.c_str(), std::ios::binary);
	std::ofstream fid3(wegFile.c_str(), std::ios::binary);
	std::ofstream fid4(timFile.c_str(), std::ios::binary);
	std::ofstream fid5(confFile.c_str());
	if (!(fid5.is_open() && fid1.is_open() && fid2.is_open() && fid3.is_open() && fid4.is_open()))
	{
		std::cout << "Cannot open the file\n";
	}
	fid1.write((const char*) (&rowIdx[0]), sizeof(int) * rowIdx.size());
	fid2.write((const char*) (&colIdx[0]), sizeof(int) * colIdx.size());
	fid3.write((const char*) (&weight[0]), sizeof(float) * weight.size());
	fid4.write((const char*) (&times[0]), sizeof(int) * times.size());
	fid5 << "Image size " << Img.m_Reso.x << "x" << Img.m_Reso.y << "\n" <<
		"Detector reso " << FanGeo.m_DetN << ";\n ViewNum " << FanGeo.m_ViwN << "\n";
	fid1.close();
	fid2.close();
	fid3.close();
	fid4.close();
	fid5.close();
}


void DEMO6()
{
	FanEAGeo FanGeo;
	Image Img;

	FanGeo.m_DetArc = 0.95928517242269f;
	FanGeo.m_DetN = 1024;
	FanGeo.m_DetCntIdx = FanGeo.m_DetN * 0.5f;
	FanGeo.m_DetStp = FanGeo.m_DetArc / FanGeo.m_DetN;
	FanGeo.m_O2D = 4.082259521484375e+02;

	FanGeo.m_S2O = 5.385200195312500e+02;
	FanGeo.m_S2D = FanGeo.m_S2O + FanGeo.m_O2D;
	FanGeo.m_ViwBeg = 0.0f;// (-2.668082275390625e+02 / 180.0* 3.14159265358979323846264);
	FanGeo.m_ViwN = 360;
	FanGeo.m_ViwStp = 3.14159265358979323846264f * 2.0f / FanGeo.m_ViwN;

	Img.m_Bias.x = 0.0f;
	Img.m_Bias.y = 0.0f;
	Img.m_Reso.x = Img.m_Reso.y = 512;
	Img.m_Size.x = Img.m_Size.y = static_cast<float>(4.484740011196460e+02);
	Img.m_Step.x = Img.m_Size.x / Img.m_Reso.x;
	Img.m_Step.y = Img.m_Size.y / Img.m_Reso.y;

	GenerateProjectionMatrix(FanGeo, Img);
	std::cout << "Demo 6 finished\n";
}






template<typename T>
void GenerateProjectionMatrix_AIM_temp(const FanEAGeo& FanGeo, const Image& Img, const std::string& FileName)
{
	std::vector<T>coeff;
	std::vector<int> rowIdx;
	std::vector<int> colIdx;

	int detIdx(0);
	int angIdx(0);
	int imgLIdx(0);
	int imgWIdx(0);
	int currowIdx(0);
	int curcolIdx(0);

	const T cntImgX = (static_cast<T>(Img.m_Reso.x) - 1) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
	const T cntImgY = (static_cast<T>(Img.m_Reso.y) - 1) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
	const T area = Img.m_Step.x * Img.m_Step.y;
	T curAng = 0;
	T cosT(0), sinT(0);
	T sour[2] = { 0, 0 };
	T SVA[3], SVB[3];
	T grid[4][3];
	T coef = 0;
	T pdist = 0;
	T dir[2], initDir[2];
	T len = 0;

	for (angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
	{
		curAng = FanGeo.m_ViwBeg + angIdx* FanGeo.m_ViwStp;
		while (curAng < 0){ curAng += _TWOPI; }
		while (curAng > _TWOPI){ curAng -= (_TWOPI); }
		cosT = cos(curAng);
		sinT = sin(curAng);
		sour[0] = -FanGeo.m_S2O * sinT;
		sour[1] = FanGeo.m_S2O * cosT;

		for (detIdx = 0; detIdx != FanGeo.m_DetN; ++detIdx)
		{
			currowIdx = angIdx * FanGeo.m_DetN + detIdx;
			calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, detIdx);
			for (imgWIdx = 0; imgWIdx != Img.m_Reso.y; ++imgWIdx)
			{
				for (imgLIdx = 0; imgLIdx != Img.m_Reso.x; ++imgLIdx)
				{
					curcolIdx = imgWIdx *Img.m_Reso.x + imgLIdx;
					grid[0][0] = (imgLIdx - cntImgX - 0.5) * Img.m_Step.x;
					grid[0][1] = (imgWIdx - cntImgY - 0.5) * Img.m_Step.y;

					dir[0] = grid[0][0] - sour[0];
					dir[1] = grid[0][1] - sour[1];
					len = hypot(dir[0], dir[1]);
					dir[0] /= len;
					dir[1] /= len;
					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
					initDir[0] = dir[0] * cosT + dir[1] * sinT;
					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;


					grid[1][0] = grid[0][0] + Img.m_Step.x;
					grid[1][1] = grid[0][1];
					dir[0] = grid[1][0] - sour[0];
					dir[1] = grid[1][1] - sour[1];
					len = hypot(dir[0], dir[1]);
					dir[0] /= len;
					dir[1] /= len;
					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
					initDir[0] = dir[0] * cosT + dir[1] * sinT;
					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
					grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

					grid[2][0] = grid[1][0];
					grid[2][1] = grid[1][1] + Img.m_Step.y;
					dir[0] = grid[2][0] - sour[0];
					dir[1] = grid[2][1] - sour[1];
					len = hypot(dir[0], dir[1]);
					dir[0] /= len;
					dir[1] /= len;
					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
					initDir[0] = dir[0] * cosT + dir[1] * sinT;
					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
					grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

					grid[3][0] = grid[0][0];
					grid[3][1] = grid[2][1];
					dir[0] = grid[3][0] - sour[0];
					dir[1] = grid[3][1] - sour[1];
					len = hypot(dir[0], dir[1]);
					dir[0] /= len;
					dir[1] /= len;
					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
					initDir[0] = dir[0] * cosT + dir[1] * sinT;
					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
					grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

					//Sort grid
					SortProj<T>(grid);

					pdist = hypot((imgLIdx - cntImgX)*Img.m_Step.x - sour[0], (imgWIdx - cntImgY)*Img.m_Step.y - sour[1]);
					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
					coef = coef / (pdist * fabs(FanGeo.m_DetStp));
					if (!IS_ZERO<T>(coef))
					{
						coeff.push_back(coef);
						rowIdx.push_back(currowIdx);
						colIdx.push_back(curcolIdx);
					}
				}
			}
		}
	}

	std::stringstream ss;
	ss << coeff.size();
	std::string num(ss.str());

	if (typeid(T) == typeid(float))
	{
		std::string coefFile = FileName + num + "f.cof";
		std::string rowFile = FileName + num + "f.row";
		std::string colFile = FileName + num + "f.col";
		std::ofstream fouCoef(coefFile.c_str(), std::ios::binary);
		fouCoef.write((char*) &(coeff[0]), sizeof(T) *coeff.size());
		std::ofstream fouCoef2(rowFile.c_str(), std::ios::binary);
		fouCoef2.write((char*) &(rowIdx[0]), sizeof(int) *rowIdx.size());
		std::ofstream fouCoef3(colFile.c_str(), std::ios::binary);
		fouCoef3.write((char*) &(colIdx[0]), sizeof(int) *colIdx.size());
	}
	if (typeid(T) == typeid(double))
	{
		std::string coefFile = FileName + num + "d.cof";
		std::string rowFile = FileName + num + "d.row";
		std::string colFile = FileName + num + "d.col";
		std::ofstream fouCoef(coefFile.c_str(), std::ios::binary);
		fouCoef.write((char*) &(coeff[0]), sizeof(T) *coeff.size());
		std::ofstream fouCoef2(rowFile.c_str(), std::ios::binary);
		fouCoef2.write((char*) &(rowIdx[0]), sizeof(int) *rowIdx.size());
		std::ofstream fouCoef3(colFile.c_str(), std::ios::binary);
		fouCoef3.write((char*) &(colIdx[0]), sizeof(int) *colIdx.size());
	}

	coeff.clear();
	rowIdx.clear();
	colIdx.clear();
}
void GenerateProjectionMatrixf_AIM(const FanEAGeo& FanGeo, const Image& Img, const std::string& FileName)
{
	GenerateProjectionMatrix_AIM_temp<float>(FanGeo, Img, FileName);
}
void GenerateProjectionMatrixd_AIM(const FanEAGeo& FanGeo, const Image& Img, const std::string& FileName)
{
	GenerateProjectionMatrix_AIM_temp<double>(FanGeo, Img, FileName);
}



template<typename T>
void GenerateMatrix_OpenMP_AIM_temp(FanEAGeo FanGeo, Image Img, const std::string& FileName)
{
	std::vector<T> coeff[32];
	std::vector<int> colIdx[32];
	std::vector<int> rowIdx[32];


	int threadID;

	//omp_set_num_threads(30);
#pragma omp parallel for 
	for (threadID = 0; threadID < 32; threadID++)
	{
		int thID = omp_get_thread_num();
		//std::cout << "current Thread : " << thID << std::endl;
		int detIdx(0);
		int angIdx(0);
		int imgLIdx(0);
		int imgWIdx(0);
		int currowIdx(0);
		int curcolIdx(0);
		const T cntImgX = (static_cast<T>(Img.m_Reso.x) - 1) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
		const T cntImgY = (static_cast<T>(Img.m_Reso.y) - 1) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
		const T area = Img.m_Step.x * Img.m_Step.y;
		T curAng = 0;
		T cosT(0), sinT(0);
		T sour[2] = { 0, 0 };
		T SVA[3], SVB[3];
		T grid[4][3];
		//T summ = 0;
		T coef = 0;
		T pdist = 0;
		T dir[2], initDir[2];
		T len = 0;
		int angN = FanGeo.m_ViwN / 30;
		std::cout << "Thread " << thID << " has " << angN << " angles\n";
		for (angIdx = 0; angIdx < angN; ++angIdx)
		{
			//std::cout << "I am Here " << thID << "\n";
			curAng = FanGeo.m_ViwBeg + (angIdx + thID * angN)* FanGeo.m_ViwStp;
			while (curAng < 0){ curAng += _TWOPI; }
			while (curAng > _TWOPI){ curAng -= _TWOPI; }
			cosT = cos(curAng);
			sinT = sin(curAng);
			sour[0] = -FanGeo.m_S2O * sinT;
			sour[1] = FanGeo.m_S2O * cosT;

			for (detIdx = 0; detIdx != FanGeo.m_DetN; ++detIdx)
			{
				currowIdx = (angIdx + thID * angN) * FanGeo.m_DetN + detIdx;
				calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, detIdx);
				for (imgWIdx = 0; imgWIdx != Img.m_Reso.y; ++imgWIdx)
				{
					for (imgLIdx = 0; imgLIdx != Img.m_Reso.x; ++imgLIdx)
					{
						curcolIdx = imgWIdx *Img.m_Reso.x + imgLIdx;
						grid[0][0] = (imgLIdx - cntImgX - 0.5) * Img.m_Step.x;
						grid[0][1] = (imgWIdx - cntImgY - 0.5) * Img.m_Step.y;

						dir[0] = grid[0][0] - sour[0];
						dir[1] = grid[0][1] - sour[1];
						len = hypot(dir[0], dir[1]);
						dir[0] /= len;
						dir[1] /= len;
						dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
						dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
						initDir[0] = dir[0] * cosT + dir[1] * sinT;
						initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
						grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;


						grid[1][0] = grid[0][0] + Img.m_Step.x;
						grid[1][1] = grid[0][1];
						dir[0] = grid[1][0] - sour[0];
						dir[1] = grid[1][1] - sour[1];
						len = hypot(dir[0], dir[1]);
						dir[0] /= len;
						dir[1] /= len;
						dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
						dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
						initDir[0] = dir[0] * cosT + dir[1] * sinT;
						initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
						grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						grid[2][0] = grid[1][0];
						grid[2][1] = grid[1][1] + Img.m_Step.y;
						dir[0] = grid[2][0] - sour[0];
						dir[1] = grid[2][1] - sour[1];
						len = hypot(dir[0], dir[1]);
						dir[0] /= len;
						dir[1] /= len;
						dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
						dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
						initDir[0] = dir[0] * cosT + dir[1] * sinT;
						initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
						grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						grid[3][0] = grid[0][0];
						grid[3][1] = grid[2][1];
						dir[0] = grid[3][0] - sour[0];
						dir[1] = grid[3][1] - sour[1];
						len = hypot(dir[0], dir[1]);
						dir[0] /= len;
						dir[1] /= len;
						dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
						dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
						initDir[0] = dir[0] * cosT + dir[1] * sinT;
						initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
						grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						//Sort grid
						SortProj<T>(grid);

						pdist = hypot((imgLIdx - cntImgX)*Img.m_Step.x - sour[0], (imgWIdx - cntImgY)*Img.m_Step.y - sour[1]);
						coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
						coef = coef / (pdist * fabs(FanGeo.m_DetStp));
						if (!IS_ZERO<T>(coef))
						{
							coeff[thID].push_back(coef);
							rowIdx[thID].push_back(currowIdx);
							colIdx[thID].push_back(curcolIdx);
						}

					}
				}
			}
		}
		std::cout << "current Thread : " << thID << std::endl;
		//std::cout << "The size of coeff " << thID << "is " << coeff[thID].size() << std::endl;
		//std::cout << "The size of coeff " << thID << "is " << rowIdx[thID].size() << std::endl;
		//std::cout << "The size of coeff " << thID << "is " << colIdx[thID].size() << std::endl;
	}

	size_t totV = 0;
	for (unsigned int i = 0; i != 32; ++i)
	{
		totV += coeff[i].size();
	}
	std::stringstream ss;
	ss << totV;
	std::string num(ss.str());

	if (typeid(T) == typeid(float))
	{
		std::string coefFile = FileName + num + "f.cof";
		std::string rowFile = FileName + num + "f.row";
		std::string colFile = FileName + num + "f.col";
		std::ofstream fouCoef(coefFile.c_str(), std::ios::binary | std::ios::app);
		std::ofstream fouCoef2(rowFile.c_str(), std::ios::binary | std::ios::app);
		std::ofstream fouCoef3(colFile.c_str(), std::ios::binary | std::ios::app);
		for (unsigned int id = 0; id != 32; ++id)
		{
			fouCoef.write((char*) &((coeff[id])[0]), sizeof(T) *coeff[id].size());
			fouCoef2.write((char*) &((rowIdx[id])[0]), sizeof(int) *rowIdx[id].size());
			fouCoef3.write((char*) &((colIdx[id])[0]), sizeof(int) *colIdx[id].size());
		}
	}
	if (typeid(T) == typeid(double))
	{
		std::string coefFile = FileName + num + "d.cof";
		std::string rowFile = FileName + num + "d.row";
		std::string colFile = FileName + num + "d.col";
		std::ofstream fouCoef(coefFile.c_str(), std::ios::binary | std::ios::app);
		std::ofstream fouCoef2(rowFile.c_str(), std::ios::binary | std::ios::app);
		std::ofstream fouCoef3(colFile.c_str(), std::ios::binary | std::ios::app);
		for (unsigned int id = 0; id != 32; ++id)
		{
			fouCoef.write((char*) &(coeff[id][0]), sizeof(T) *coeff[id].size());
			fouCoef2.write((char*) &(rowIdx[id][0]), sizeof(int) *rowIdx[id].size());
			fouCoef3.write((char*) &(colIdx[id][0]), sizeof(int) *colIdx[id].size());
		}
	}
}


void GenerateMatrix_OpenMP_AIM_float(const FanEAGeo& FanGeo, const Image& Img, const std::string& FileName)
{
	GenerateMatrix_OpenMP_AIM_temp<float>(FanGeo, Img, FileName);
}

void GenerateMatrix_OpenMP_AIM_double(const FanEAGeo& FanGeo, const Image& Img, const std::string& FileName)
{
	GenerateMatrix_OpenMP_AIM_temp<double>(FanGeo, Img, FileName);
}


void genProjectionMatrix_AIMf(const FanEDGeo FanGeo, const Image Img)
{
	genProjectionMatrix_AIM_template<float>(FanGeo, Img);
}
void genProjectionMatrix_AIMd(const FanEDGeo FanGeo, const Image Img)
{
	genProjectionMatrix_AIM_template<double>(FanGeo, Img);
}


void genProjectionMatrix_AIMf(const FanEAGeo FanGeo, const Image Img)
{
	genProjectionMatrix_AIM_template<float>(FanGeo, Img);
}

void genProjectionMatrix_AIMd(const FanEAGeo FanGeo, const Image Img)
{
	genProjectionMatrix_AIM_template<double>(FanGeo, Img);
}






