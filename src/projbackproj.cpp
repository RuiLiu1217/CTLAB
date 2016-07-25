#include "utilities.hpp"
#include "projbackproj.hpp"

void bakproj_PIXEL_CPU_OPENMP(float* donp, float* dimg, float* dmsk, const FanEAGeo& FanGeo, const Image& Img)
{
	int idY = 0;
#pragma omp parallel for
	for (idY = 0; idY < Img.m_Reso.y; idY++)
	{
		int idX = 0;
		for (idX = 0; idX < Img.m_Reso.x; idX++)
		{

			float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
			//cur rotation angle;

			float curAng, cosT, sinT;
			Ray2D ray;
			float2 minBox, maxBox, curImg, curDet, initDet, sincosAng;
			float biasAng, tnear, tfar, weg, detID;
			int flrID, ceilID;
			float summ(0);
			float msk = dmsk[idY * Img.m_Reso.x + idX];
			for (unsigned int angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
			{
				curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
				cosT = cosf(curAng);
				sinT = sinf(curAng);
				//current source position, assume the initial position is on the positive Y axis;
				ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);

				minBox = make_float2(
					MINO.x + idX * Img.m_Step.x,
					MINO.y + idY * Img.m_Step.y);
				maxBox = minBox + make_float2(Img.m_Step.x, Img.m_Step.y);
				curImg = minBox + 0.5f * make_float2(Img.m_Step.x, Img.m_Step.y);
				ray.d = normalize(curImg - ray.o);


				//current detector element Cartesian coordinates;
				curDet = ray.o + ray.d * FanGeo.m_S2D;
				//original detector position;
				initDet = rotation(curDet, cosT, -sinT);

				//Bias angle;

				sincosAng.x = initDet.x / FanGeo.m_S2D; //sin(ang)
				sincosAng.y = (initDet.y - FanGeo.m_S2O) / (-FanGeo.m_S2D); //cos(ang)
				biasAng = atan(sincosAng.x / sincosAng.y);

				detID = biasAng / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
				if (detID < 0)
				{
					dimg[idY * Img.m_Reso.x + idX] += 0;
					continue;//return;
				}
				if (detID > FanGeo.m_DetN)
				{
					dimg[idY * Img.m_Reso.x + idX] += 0;
					continue;//return;
				}
				tnear = 0;
				tfar = 0;
				flrID = detID;
				ceilID = ceilf(detID);
				intersectBox(ray, minBox, maxBox, &tnear, &tfar);

				weg = (tfar - tnear); //no weighting;
				if (!IS_ZERO(flrID - ceilID))
				{
					summ += ((donp[angIdx * FanGeo.m_DetN + flrID] * (ceilID - detID) + donp[angIdx * FanGeo.m_DetN + ceilID] * (detID - flrID)) * weg);
				}
				else
				{
					summ += (donp[angIdx * FanGeo.m_DetN + flrID] * weg);
				}
			}
			dimg[idY * Img.m_Reso.x + idX] = summ * msk;
		}
	}
}
void bakproj_PIXEL_CPU(float* donp, float* dimg, float* dmsk, const FanEAGeo& FanGeo, const Image& Img)
{
	unsigned int idX = 0, idY = 0;

	for (idY = 0; idY < Img.m_Reso.y; idY++)
	{
		for (idX = 0; idX < Img.m_Reso.x; idX++)
		{

			float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
			//cur rotation angle;

			float curAng, cosT, sinT;
			Ray2D ray;
			float2 minBox, maxBox, curImg, curDet, initDet, sincosAng;
			float biasAng, tnear, tfar, weg, detID;
			int flrID, ceilID;
			float summ(0);
			float msk = dmsk[idY * Img.m_Reso.x + idX];
			for (unsigned int angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
			{
				curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
				cosT = cosf(curAng);
				sinT = sinf(curAng);
				//current source position, assume the initial position is on the positive Y axis;
				ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);

				minBox = make_float2(
					MINO.x + idX * Img.m_Step.x,
					MINO.y + idY * Img.m_Step.y);
				maxBox = minBox + make_float2(Img.m_Step.x, Img.m_Step.y);
				curImg = minBox + 0.5f * make_float2(Img.m_Step.x, Img.m_Step.y);
				ray.d = normalize(curImg - ray.o);


				//current detector element Cartesian coordinates;
				curDet = ray.o + ray.d * FanGeo.m_S2D;
				//original detector position;
				initDet = rotation(curDet, cosT, -sinT);

				//Bias angle;

				sincosAng.x = initDet.x / FanGeo.m_S2D; //sin(ang)
				sincosAng.y = (initDet.y - FanGeo.m_S2O) / (-FanGeo.m_S2D); //cos(ang)
				biasAng = atan(sincosAng.x / sincosAng.y);

				detID = biasAng / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
				if (detID < 0)
				{
					dimg[idY * Img.m_Reso.x + idX] += 0;
					continue;//return;
				}
				if (detID > FanGeo.m_DetN)
				{
					dimg[idY * Img.m_Reso.x + idX] += 0;
					continue;//return;
				}
				tnear = 0;
				tfar = 0;
				flrID = detID;
				ceilID = ceilf(detID);
				intersectBox(ray, minBox, maxBox, &tnear, &tfar);

				weg = (tfar - tnear); //no weighting;
				if (!IS_ZERO(flrID - ceilID))
				{
					summ += ((donp[angIdx * FanGeo.m_DetN + flrID] * (ceilID - detID) + donp[angIdx * FanGeo.m_DetN + ceilID] * (detID - flrID)) * weg);
				}
				else
				{
					summ += (donp[angIdx * FanGeo.m_DetN + flrID] * weg);
				}
			}
			dimg[idY * Img.m_Reso.x + idX] = summ * msk;
		}
	}
}


void bakproj_BOXED_CPU(float* dprj, float* dimg, const FanEAGeo& FanGeo, const Image& Img)
{
	unsigned int i = 0;
	unsigned int j = 0;
	float2 minBox, maxBox, initD;
	unsigned int angIdx = 0;
	Ray2D ray;
	int minID = 10000;
	int maxID = -10000;
	int curID = 0;
	float ang = 0.0f;

	float summ(0.0f), tnear(0.0f), tfar(0.0f);

	float cosT(0.0f);
	float sinT(0.0f);
	float detID(0.0f);
	int hit = 0;
	//float curAng(0.0f);
	float2 curDet;
	for (j = 0; j != Img.m_Reso.y; ++j)
	{
		for (i = 0; i != Img.m_Reso.x; ++i)
		{
			minBox = make_float2(
				-Img.m_Size.x * 0.5f + Img.m_Bias.x + i * Img.m_Step.x,
				-Img.m_Size.y * 0.5f + Img.m_Bias.y + j * Img.m_Step.y);
			maxBox = make_float2(
				-Img.m_Size.x * 0.5f + Img.m_Bias.x + (i + 1) * Img.m_Step.x,
				-Img.m_Size.y * 0.5f + Img.m_Bias.y + (j + 1) * Img.m_Step.y);
			summ = 0;
			for (angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
			{
				minID = 10000;
				maxID = -10000;

				cosT = cos(FanGeo.m_ViwBeg + (double) angIdx * FanGeo.m_ViwStp);
				sinT = sin(FanGeo.m_ViwBeg + (double) angIdx * FanGeo.m_ViwStp);
				ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT); //current ray position;

				ray.d = normalize(minBox - ray.o);
				initD = rotation(ray.o + FanGeo.m_S2D * ray.d, cosT, -sinT);
				detID = atanf(initD.x / (FanGeo.m_S2O - initD.y)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
				minID = static_cast<int>((detID < minID) ? detID : minID);
				maxID = static_cast<int>((detID > maxID) ? ceilf(detID) : maxID);


				ray.d = normalize(minBox + make_float2(Img.m_Step.x, 0) - ray.o);
				//initD = rotation(make_float2(ray.d.x * FanGeo.m_S2D,ray.d.y * FanGeo.m_S2D), cosT,-sinT); 
				initD = rotation(ray.o + FanGeo.m_S2D * ray.d, cosT, -sinT);
				detID = atanf(initD.x / (FanGeo.m_S2O - initD.y)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
				minID = static_cast<int>((detID < minID) ? detID : minID);
				maxID = static_cast<int>((detID > maxID) ? ceilf(detID) : maxID);


				ray.d = normalize(minBox + make_float2(0, Img.m_Step.y) - ray.o);
				initD = rotation(ray.o + FanGeo.m_S2D * ray.d, cosT, -sinT);
				detID = atanf(initD.x / (FanGeo.m_S2O - initD.y)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
				minID = static_cast<int>((detID < minID) ? detID : minID);
				maxID = static_cast<int>((detID > maxID) ? ceilf(detID) : maxID);


				ray.d = normalize(minBox + make_float2(Img.m_Step.x, Img.m_Step.y) - ray.o);
				initD = rotation(ray.o + FanGeo.m_S2D * ray.d, cosT, -sinT);
				detID = atanf(initD.x / (FanGeo.m_S2O - initD.y)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
				minID = static_cast<int>((detID < minID) ? detID : minID);
				maxID = static_cast<int>((detID > maxID) ? ceilf(detID) : maxID);


				if (minID != 10000){ minID = (minID < 0) ? 0 : minID; }
				if (maxID != -10000){ maxID = (maxID > static_cast<int>(FanGeo.m_DetN - 1)) ? static_cast<int>(FanGeo.m_DetN - 1) : maxID; }



				for (curID = minID; curID <= maxID; ++curID)
				{
					//���index�㵱ǰdetector λ��;
					ang = (curID - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
					curDet = rotation(make_float2(sin(static_cast<double>(ang)) * FanGeo.m_S2D, -cos(static_cast<double>(ang))*FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
					ray.d = normalize(curDet - ray.o);
					tnear = 0.0f; tfar = 0.0f;
					hit = intersectBox(ray, minBox, maxBox, &tnear, &tfar);
					if (hit)
					{
						//intersectLength += (tfar - tnear);
						summ += (tfar - tnear) * dprj[angIdx * FanGeo.m_DetN + curID];
					}
					//summ += weightVal;
				}

			}//end for angIdx
			dimg[j * Img.m_Reso.x + i] = summ;
		} // end for i
	} //end for j
}

void bakproj_BOXED_CPU_OPENMP(float* dprj, float* dimg, const FanEAGeo& FanGeo, const Image& Img)
{
	int j = 0;
#pragma omp parallel for
	for (j = 0; j < Img.m_Reso.y; ++j)
	{
		unsigned int i = 0;

		float2 minBox, maxBox, initD;
		unsigned int angIdx = 0;
		Ray2D ray;
		int minID = 10000;
		int maxID = -10000;
		int curID = 0;
		float ang = 0.0f;

		float summ(0.0f), tnear(0.0f), tfar(0.0f);

		float cosT(0.0f);
		float sinT(0.0f);
		float detID(0.0f);
		int hit = 0;
		//float curAng(0.0f);
		float2 curDet;
		for (i = 0; i != Img.m_Reso.x; ++i)
		{
			minBox = make_float2(
				-Img.m_Size.x * 0.5f + Img.m_Bias.x + i * Img.m_Step.x,
				-Img.m_Size.y * 0.5f + Img.m_Bias.y + j * Img.m_Step.y);
			maxBox = make_float2(
				-Img.m_Size.x * 0.5f + Img.m_Bias.x + (i + 1) * Img.m_Step.x,
				-Img.m_Size.y * 0.5f + Img.m_Bias.y + (j + 1) * Img.m_Step.y);
			summ = 0;
			for (angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
			{
				minID = 10000;
				maxID = -10000;

				cosT = cos(FanGeo.m_ViwBeg + (double) angIdx * FanGeo.m_ViwStp);
				sinT = sin(FanGeo.m_ViwBeg + (double) angIdx * FanGeo.m_ViwStp);
				ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT); //current ray position;

				ray.d = normalize(minBox - ray.o);
				initD = rotation(ray.o + FanGeo.m_S2D * ray.d, cosT, -sinT);
				detID = atanf(initD.x / (FanGeo.m_S2O - initD.y)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
				minID = static_cast<int>((detID < minID) ? detID : minID);
				maxID = static_cast<int>((detID > maxID) ? ceilf(detID) : maxID);


				ray.d = normalize(minBox + make_float2(Img.m_Step.x, 0) - ray.o);
				//initD = rotation(make_float2(ray.d.x * FanGeo.m_S2D,ray.d.y * FanGeo.m_S2D), cosT,-sinT); 
				initD = rotation(ray.o + FanGeo.m_S2D * ray.d, cosT, -sinT);
				detID = atanf(initD.x / (FanGeo.m_S2O - initD.y)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
				minID = static_cast<int>((detID < minID) ? detID : minID);
				maxID = static_cast<int>((detID > maxID) ? ceilf(detID) : maxID);


				ray.d = normalize(minBox + make_float2(0, Img.m_Step.y) - ray.o);
				initD = rotation(ray.o + FanGeo.m_S2D * ray.d, cosT, -sinT);
				detID = atanf(initD.x / (FanGeo.m_S2O - initD.y)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
				minID = static_cast<int>((detID < minID) ? detID : minID);
				maxID = static_cast<int>((detID > maxID) ? ceilf(detID) : maxID);


				ray.d = normalize(minBox + make_float2(Img.m_Step.x, Img.m_Step.y) - ray.o);
				initD = rotation(ray.o + FanGeo.m_S2D * ray.d, cosT, -sinT);
				detID = atanf(initD.x / (FanGeo.m_S2O - initD.y)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
				minID = static_cast<int>((detID < minID) ? detID : minID);
				maxID = static_cast<int>((detID > maxID) ? ceilf(detID) : maxID);


				if (minID != 10000){ minID = (minID < 0) ? 0 : minID; }
				if (maxID != -10000){ maxID = (maxID > static_cast<int>(FanGeo.m_DetN - 1)) ? static_cast<int>(FanGeo.m_DetN - 1) : maxID; }



				for (curID = minID; curID <= maxID; ++curID)
				{
					//���index�㵱ǰdetector λ��;
					ang = (curID - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
					curDet = rotation(make_float2(sin(static_cast<double>(ang)) * FanGeo.m_S2D, -cos(static_cast<double>(ang))*FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
					ray.d = normalize(curDet - ray.o);
					tnear = 0.0f; tfar = 0.0f;
					hit = intersectBox(ray, minBox, maxBox, &tnear, &tfar);
					if (hit)
					{
						//intersectLength += (tfar - tnear);
						summ += (tfar - tnear) * dprj[angIdx * FanGeo.m_DetN + curID];
					}
					//summ += weightVal;
				}

			}//end for angIdx
			dimg[j * Img.m_Reso.x + i] = summ;
		} // end for i
	} //end for j
}


void proj_CPU(float* dvol, float* draw, const ConeEAGeo ConeGeo, const Volume Vol)
{
	unsigned int detId(0), rowId(0), angIdx(0);
	const float3 MINO = make_float3(
		-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
		-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
		-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
	float curAng;// = ConeGeo.m_ViwBeg + ConeGeo.m_ViwStp * angIdx;
	float cosT;
	float sinT;
	Ray ray;
	float ang(0);
	float totWeight(0);
	float3 curDetPos;// = rotation(make_float3(

	for (angIdx = 0; angIdx < ConeGeo.m_ViwN; ++angIdx)
	{
		curAng = ConeGeo.m_ViwBeg + ConeGeo.m_ViwStp * angIdx;
		cosT = cos(static_cast<double>(curAng));
		sinT = sin(static_cast<double>(curAng));

		for (rowId = 0; rowId != ConeGeo.m_DetHN; ++rowId)
		{

			for (detId = 0; detId < ConeGeo.m_DetN; ++detId)
			{
				// judging two indices addressing method, from large to small or small to large
				ang = ((float) detId - ConeGeo.m_DetCntIdx.x) * ConeGeo.m_DetStp;

				// current detector element position
				curDetPos = rotation(make_float3(
					sin(static_cast<double>(ang)) * ConeGeo.m_S2D,
					-cos(static_cast<double>(ang)) * ConeGeo.m_S2D + ConeGeo.m_S2O,
					((float) rowId - ConeGeo.m_DetCntIdx.y) * ConeGeo.m_DetHStp),
					cosT, sinT);

				// x-ray direction;
				ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
				ray.d = normalize(curDetPos - ray.o);

				totWeight = 0;
				draw[(angIdx * ConeGeo.m_DetHN + rowId) * ConeGeo.m_DetN + detId] =
					calSiddonOneRayKer(ray.o.x, ray.o.y, ray.o.z, curDetPos.x, curDetPos.y, curDetPos.z, MINO.x, MINO.y, MINO.z,
					Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z, Vol.m_Reso.x, Vol.m_Reso.y, Vol.m_Reso.z, dvol, &totWeight);
			}
		}
	}

}


void proj_CPU_OPENMP(float* dvol, float* draw, const ConeEAGeo ConeGeo, const Volume Vol)
{

	int angIdx(0);

	const float3 MINO = make_float3(
		-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
		-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
		-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
#pragma omp parallel for
	for (angIdx = 0; angIdx < ConeGeo.m_ViwN; ++angIdx)
	{
		unsigned int detId(0), rowId(0);
		float curAng = ConeGeo.m_ViwBeg + ConeGeo.m_ViwStp * angIdx;
		float cosT = cos(static_cast<double>(curAng));
		float sinT = sin(static_cast<double>(curAng));
		Ray ray;
		float ang(0);
		float totWeight(0);
		float3 curDetPos;
		for (rowId = 0; rowId != ConeGeo.m_DetHN; ++rowId)
		{
			for (detId = 0; detId < ConeGeo.m_DetN; ++detId)
			{
				// judging two indices addressing method, from large to small or small to large
				ang = ((float) detId - ConeGeo.m_DetCntIdx.x) * ConeGeo.m_DetStp;

				// current detector element position
				curDetPos = rotation(make_float3(
					sin(static_cast<double>(ang)) * ConeGeo.m_S2D,
					-cos(static_cast<double>(ang)) * ConeGeo.m_S2D + ConeGeo.m_S2O,
					((float) rowId - ConeGeo.m_DetCntIdx.y) * ConeGeo.m_DetHStp),
					cosT, sinT);

				// x-ray direction;
				ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
				ray.d = normalize(curDetPos - ray.o);

				totWeight = 0;
				draw[(angIdx * ConeGeo.m_DetHN + rowId) * ConeGeo.m_DetN + detId] =
					calSiddonOneRayKer(ray.o.x, ray.o.y, ray.o.z, curDetPos.x, curDetPos.y, curDetPos.z, MINO.x, MINO.y, MINO.z,
					Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z, Vol.m_Reso.x, Vol.m_Reso.y, Vol.m_Reso.z, dvol, &totWeight);
			}
		}
	}

}



template<typename T>
void proj_AIM_CPU_temp(T* dprj, T* dimg, const FanEAGeo& FanGeo, const Image& Img)
{
	int angIdx = 0;
	const T cntImgX = (static_cast<T>(Img.m_Reso.x) - 1) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
	const T cntImgY = (static_cast<T>(Img.m_Reso.y) - 1) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
	const T area(Img.m_Step.x * Img.m_Step.y);
	const int extraIdx = 1;
	T curAng;// = FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp;

	T cosT;
	T sinT;
	T sour[2]; // = { -FanGeo.m_S2O * sinT, FanGeo.m_S2O * cosT };
	T SVA[3];
	T SVB[3];

	unsigned int sliceIdx;
	T summ;
	T coord;
	T minC, maxC;
	int maxIdx;
	int curIdx;
	T grid[4][3];
	T dir[2];
	T initDir[2];
	T len;
	T pdist, coef;

	T curDirAng;// = (detIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp + curAng;

	for (angIdx = 0; angIdx < FanGeo.m_ViwN; ++angIdx)
	{
		curAng = FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp;
		while (curAng < 0){ curAng += _TWOPI; }
		while (curAng> _TWOPI){ curAng -= _TWOPI; }
		cosT = cos(curAng);
		sinT = sin(curAng);
		sour[0] = -FanGeo.m_S2O * sinT;
		sour[1] = FanGeo.m_S2O * cosT;
		for (int detIdx = 0; detIdx != FanGeo.m_DetN; ++detIdx)
		{
			calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, detIdx);
			curDirAng = (detIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp + curAng;
			while (curDirAng < 0)	{ curDirAng += _TWOPI; }
			while (curDirAng > _TWOPI){ curDirAng -= _TWOPI; }
			if (curDirAng <= _PI_4 || curDirAng > _7PI_4) //��һ�����;
			{
				//����y����ֲ� ??����������;
				summ = 0;
				for (sliceIdx = 0; sliceIdx != Img.m_Reso.y; ++sliceIdx)
				{
					coord = (static_cast<T>(sliceIdx) -cntImgY - 0.5) * Img.m_Step.y;
					minC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
					maxC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
					if (maxC < minC)
					{
						pdist = minC;
						minC = maxC;
						maxC = pdist;
					}
					curIdx = int(minC / Img.m_Step.x + cntImgX) - extraIdx;
					maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + extraIdx;
					if (curIdx > static_cast<int>(Img.m_Reso.x - 1) || maxIdx < 0)
					{

						continue;
					}
					if (curIdx < 0)
					{
						curIdx = 0;
					}
					if (maxIdx > Img.m_Reso.x - 1)
					{
						maxIdx = Img.m_Reso.x - 1;
					}

					for (; curIdx <= maxIdx; curIdx++)
					{
						grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
						grid[0][1] = coord;
						grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
						grid[1][1] = coord;
						grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
						grid[2][1] = coord + Img.m_Step.y;
						grid[3][0] = (curIdx - cntImgX - 0.5)* Img.m_Step.x;
						grid[3][1] = coord + Img.m_Step.y;


						//�����ĸ����Ӧ��ǰ��det index
						dir[0] = grid[0][0] - sour[0];
						dir[1] = grid[0][1] - sour[1];
						len = hypot(dir[0], dir[1]);
						dir[0] /= len;
						dir[1] /= len;
						dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
						dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
						initDir[0] = dir[0] * cosT + dir[1] * sinT;
						initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
						grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index

						//�����ĸ����Ӧ��ǰ��det index
						dir[0] = grid[0][0] - sour[0];
						dir[1] = grid[0][1] - sour[1];
						len = hypot(dir[0], dir[1]);
						dir[0] /= len;
						dir[1] /= len;
						dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
						dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
						initDir[0] = dir[0] * cosT + dir[1] * sinT;
						initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
						grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index

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

						//����;
						SortProj<T>(grid);

						pdist = hypot((static_cast<T>(curIdx) -cntImgX)*static_cast<T>(Img.m_Step.x) - sour[0], (static_cast<T>(sliceIdx) -cntImgY)*static_cast<T>(Img.m_Step.y) - sour[1]);
						coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area) / (pdist * FanGeo.m_DetStp);
						summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
					}
				}
				dprj[angIdx* FanGeo.m_DetN + detIdx] = summ;

			}
			else if (curDirAng > _PI_4 && curDirAng <= _3PI_4)
			{
				summ = 0;
				for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
				{
					coord = (sliceIdx - cntImgX - 0.5)* Img.m_Step.x;
					//�������λ�ù����ཻ�����;
					minC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
					maxC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
					if (maxC < minC)
					{
						pdist = minC;
						minC = maxC;
						maxC = pdist;
					}
					curIdx = int(minC / Img.m_Step.y + cntImgY) - extraIdx;
					maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + extraIdx;

					if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0)
					{
						continue;
					}

					if (curIdx < 0)
					{
						curIdx = 0;
					}

					if (maxIdx > Img.m_Reso.y - 1)
					{
						maxIdx = Img.m_Reso.y - 1;
					}

					for (; curIdx <= maxIdx; ++curIdx)
					{
						//��grid
						grid[0][0] = coord;
						grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
						grid[1][0] = coord + Img.m_Step.x;
						grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
						grid[2][0] = coord + Img.m_Step.x;
						grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
						grid[3][0] = coord;
						grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;

						//�����ĸ����Ӧ��ǰ��det index
						dir[0] = grid[0][0] - sour[0];
						dir[1] = grid[0][1] - sour[1];
						len = hypot(dir[0], dir[1]);
						dir[0] /= len;
						dir[1] /= len;
						dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
						dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
						initDir[0] = dir[0] * cosT + dir[1] * sinT;
						initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
						grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index

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

						//����;
						SortProj<T>(grid);
						pdist = hypot((static_cast<T>(sliceIdx) -cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(curIdx) -cntImgY)*Img.m_Step.y - sour[1]);
						coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area) / (pdist * FanGeo.m_DetStp);
						summ += coef * dimg[curIdx * Img.m_Reso.x + sliceIdx];

					} //End for one slice
				}// End all slices
				dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
			}
			else if (curDirAng > _3PI_4 && curDirAng <= _5PI_4)
			{
				summ = 0;
				for (sliceIdx = 0; sliceIdx < Img.m_Reso.y; ++sliceIdx)
				{
					coord = (sliceIdx - cntImgY + 0.5)* Img.m_Step.y;
					//�������λ�ù����ཻ�����;
					maxC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
					minC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
					if (maxC < minC)
					{
						pdist = minC;
						minC = maxC;
						maxC = pdist;
					}
					curIdx = int(minC / Img.m_Step.x + cntImgX) - extraIdx;
					maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + extraIdx;

					if (curIdx > static_cast<int>(Img.m_Reso.x) || maxIdx < 0)
					{
						continue;
					}

					if (curIdx < 0)
					{
						curIdx = 0;
					}

					if (maxIdx > Img.m_Reso.x - 1)
					{
						maxIdx = Img.m_Reso.x - 1;
					}

					for (; curIdx <= maxIdx; ++curIdx)
					{
						//��grid
						grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
						grid[0][1] = coord - Img.m_Step.y;
						grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
						grid[1][1] = coord - Img.m_Step.y;
						grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
						grid[2][1] = coord;
						grid[3][0] = (curIdx - cntImgX - 0.5)* Img.m_Step.x;
						grid[3][1] = coord;

						//�����ĸ����Ӧ��ǰ��det index
						dir[0] = grid[0][0] - sour[0];
						dir[1] = grid[0][1] - sour[1];
						len = hypot(dir[0], dir[1]);
						dir[0] /= len;
						dir[1] /= len;
						dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
						dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
						initDir[0] = dir[0] * cosT + dir[1] * sinT;
						initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
						grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index

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

						//����;
						SortProj<T>(grid);
						pdist = hypot((static_cast<T>(curIdx) -cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(sliceIdx) -cntImgY)*Img.m_Step.y - sour[1]);
						coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area) / (pdist * FanGeo.m_DetStp);
						summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];

					} //End for one slice
				}// End all slices
				dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
			}
			else
			{
				summ = 0;
				for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
				{
					coord = (sliceIdx - cntImgX + 0.5)* Img.m_Step.x;
					//�������λ�ù����ཻ�����;
					maxC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
					minC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
					if (maxC < minC)
					{
						pdist = minC;
						minC = maxC;
						maxC = pdist;
					}
					curIdx = int(minC / Img.m_Step.y + cntImgY) - extraIdx;
					maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + extraIdx;

					if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0)
					{
						continue;
					}

					if (curIdx < 0)
					{
						curIdx = 0;
					}

					if (maxIdx > Img.m_Reso.y - 1)
					{
						maxIdx = Img.m_Reso.y - 1;
					}

					for (; curIdx <= maxIdx; ++curIdx)
					{
						//��grid
						grid[0][0] = coord - Img.m_Step.x;
						grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
						grid[1][0] = coord;
						grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
						grid[2][0] = coord - Img.m_Step.x;
						grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
						grid[3][0] = coord;
						grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;

						//�����ĸ����Ӧ��ǰ��det index
						dir[0] = grid[0][0] - sour[0];
						dir[1] = grid[0][1] - sour[1];
						len = hypot(dir[0], dir[1]);
						dir[0] /= len;
						dir[1] /= len;
						dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
						dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
						initDir[0] = dir[0] * cosT + dir[1] * sinT;
						initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
						grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index

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

						//����;
						SortProj<T>(grid);
						pdist = hypot((static_cast<T>(sliceIdx) -cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(curIdx) -cntImgY)*Img.m_Step.y - sour[1]);
						coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
						coef = coef / (pdist * FanGeo.m_DetStp);
						summ += coef * dimg[curIdx* Img.m_Reso.x + sliceIdx];
						//summ += *coef;
					} //End for one slice
				}// End all slices
				dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
			}
		}
	}
}
void proj_AIM_CPU(float* prj, float* img, const FanEAGeo& FanGeo, const Image& Img)
{
	proj_AIM_CPU_temp<float>(prj, img, FanGeo, Img);
}
void proj_AIM_CPU(double* prj, double* img, const FanEAGeo& FanGeo, const Image& Img)
{
	proj_AIM_CPU_temp<double>(prj, img, FanGeo, Img);
}

template<typename T>
void proj_AIM_CPU_OPENMP_temp(T* dprj, T* dimg, const FanEAGeo& FanGeo, const Image& Img)
{
	int angIdx = 0;
	const T cntImgX = (static_cast<T>(Img.m_Reso.x) - 1) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
	const T cntImgY = (static_cast<T>(Img.m_Reso.y) - 1) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
	const T area(Img.m_Step.x * Img.m_Step.y);
	const int extraIdx = 1;
#pragma omp parallel for
	for (angIdx = 0; angIdx < FanGeo.m_ViwN; ++angIdx)
	{
		T curAng = FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp;
		while (curAng < 0){ curAng += _TWOPI; }
		while (curAng > _TWOPI){ curAng -= _TWOPI; }
		T cosT = cos(curAng);
		T sinT = sin(curAng);
		T sour[2] = { -FanGeo.m_S2O * sinT, FanGeo.m_S2O * cosT };
		T SVA[3];
		T SVB[3];

		unsigned int sliceIdx;
		T summ;
		T coord;
		T minC, maxC;
		int maxIdx;
		int curIdx;
		T grid[4][3];
		T dir[2];
		T initDir[2];
		T len;
		T pdist, coef;
		curAng = FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp;
		cosT = cos(curAng);
		sinT = sin(curAng);

		sour[0] = -FanGeo.m_S2O * sinT;
		sour[1] = FanGeo.m_S2O * cosT;
		while (curAng < 0){ curAng += _TWOPI; }
		while (curAng > _TWOPI){ curAng -= _TWOPI; }
		T curDirAng;// = (detIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp + curAng;
		for (int detIdx = 0; detIdx != FanGeo.m_DetN; ++detIdx)
		{
			calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, detIdx);
			curDirAng = (detIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp + curAng;
			while (curDirAng < 0)	{ curDirAng += _TWOPI; }
			while (curDirAng > _TWOPI){ curDirAng -= _TWOPI; }
			if (curDirAng <= _PI_4 || curDirAng > _7PI_4)
			{
				//����y����ֲ� ??����������;
				summ = 0;
				for (sliceIdx = 0; sliceIdx != Img.m_Reso.y; ++sliceIdx)
				{
					coord = (static_cast<T>(sliceIdx) -cntImgY - 0.5) * Img.m_Step.y;
					minC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
					maxC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
					if (maxC < minC)
					{
						pdist = minC;
						minC = maxC;
						maxC = pdist;
					}
					curIdx = int(minC / Img.m_Step.x + cntImgX) - extraIdx;
					maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + extraIdx;
					if (curIdx > static_cast<int>(Img.m_Reso.x - 1) || maxIdx < 0)
					{

						continue;
					}
					if (curIdx < 0)
					{
						curIdx = 0;
					}
					if (maxIdx > Img.m_Reso.x - 1)
					{
						maxIdx = Img.m_Reso.x - 1;
					}

					for (; curIdx <= maxIdx; curIdx++)
					{
						grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
						grid[0][1] = coord;
						grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
						grid[1][1] = coord;
						grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
						grid[2][1] = coord + Img.m_Step.y;
						grid[3][0] = (curIdx - cntImgX - 0.5)* Img.m_Step.x;
						grid[3][1] = coord + Img.m_Step.y;


						//�����ĸ����Ӧ��ǰ��det index
						dir[0] = grid[0][0] - sour[0];
						dir[1] = grid[0][1] - sour[1];
						len = hypot(dir[0], dir[1]);
						dir[0] /= len;
						dir[1] /= len;
						dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
						dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
						initDir[0] = dir[0] * cosT + dir[1] * sinT;
						initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
						grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index

						//�����ĸ����Ӧ��ǰ��det index
						dir[0] = grid[0][0] - sour[0];
						dir[1] = grid[0][1] - sour[1];
						len = hypot(dir[0], dir[1]);
						dir[0] /= len;
						dir[1] /= len;
						dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
						dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
						initDir[0] = dir[0] * cosT + dir[1] * sinT;
						initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
						grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index

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

						//����;
						SortProj<T>(grid);

						pdist = hypot((static_cast<T>(curIdx) -cntImgX)*static_cast<T>(Img.m_Step.x) - sour[0], (static_cast<T>(sliceIdx) -cntImgY)*static_cast<T>(Img.m_Step.y) - sour[1]);
						coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area) / (pdist * FanGeo.m_DetStp);
						summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
					}
				}
				dprj[angIdx* FanGeo.m_DetN + detIdx] = summ;

			}
			else if (curDirAng > _PI_4 && curDirAng <= _3PI_4)
			{
				summ = 0;
				for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
				{
					coord = (sliceIdx - cntImgX - 0.5)* Img.m_Step.x;
					//�������λ�ù����ཻ�����;
					minC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
					maxC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
					if (maxC < minC)
					{
						pdist = minC;
						minC = maxC;
						maxC = pdist;
					}
					curIdx = int(minC / Img.m_Step.y + cntImgY) - extraIdx;
					maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + extraIdx;

					if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0)
					{
						continue;
					}

					if (curIdx < 0)
					{
						curIdx = 0;
					}

					if (maxIdx > Img.m_Reso.y - 1)
					{
						maxIdx = Img.m_Reso.y - 1;
					}

					for (; curIdx <= maxIdx; ++curIdx)
					{
						//��grid
						grid[0][0] = coord;
						grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
						grid[1][0] = coord + Img.m_Step.x;
						grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
						grid[2][0] = coord + Img.m_Step.x;
						grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
						grid[3][0] = coord;
						grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;

						//�����ĸ����Ӧ��ǰ��det index
						dir[0] = grid[0][0] - sour[0];
						dir[1] = grid[0][1] - sour[1];
						len = hypot(dir[0], dir[1]);
						dir[0] /= len;
						dir[1] /= len;
						dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
						dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
						initDir[0] = dir[0] * cosT + dir[1] * sinT;
						initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
						grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index

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

						//����;
						SortProj<T>(grid);
						pdist = hypot((static_cast<T>(sliceIdx) -cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(curIdx) -cntImgY)*Img.m_Step.y - sour[1]);
						coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area) / (pdist * FanGeo.m_DetStp);
						summ += coef * dimg[curIdx * Img.m_Reso.x + sliceIdx];

					} //End for one slice
				}// End all slices
				dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
			}
			else if (curDirAng > _3PI_4 && curDirAng <= _5PI_4)
			{
				summ = 0;
				for (sliceIdx = 0; sliceIdx < Img.m_Reso.y; ++sliceIdx)
				{
					coord = (sliceIdx - cntImgY + 0.5)* Img.m_Step.y;
					//�������λ�ù����ཻ�����;
					maxC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
					minC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
					if (maxC < minC)
					{
						pdist = minC;
						minC = maxC;
						maxC = pdist;
					}
					curIdx = int(minC / Img.m_Step.x + cntImgX) - extraIdx;
					maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + extraIdx;

					if (curIdx > static_cast<int>(Img.m_Reso.x) || maxIdx < 0)
					{
						continue;
					}

					if (curIdx < 0)
					{
						curIdx = 0;
					}

					if (maxIdx > Img.m_Reso.x - 1)
					{
						maxIdx = Img.m_Reso.x - 1;
					}

					for (; curIdx <= maxIdx; ++curIdx)
					{
						//��grid
						grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
						grid[0][1] = coord - Img.m_Step.y;
						grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
						grid[1][1] = coord - Img.m_Step.y;
						grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
						grid[2][1] = coord;
						grid[3][0] = (curIdx - cntImgX - 0.5)* Img.m_Step.x;
						grid[3][1] = coord;

						//�����ĸ����Ӧ��ǰ��det index
						dir[0] = grid[0][0] - sour[0];
						dir[1] = grid[0][1] - sour[1];
						len = hypot(dir[0], dir[1]);
						dir[0] /= len;
						dir[1] /= len;
						dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
						dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
						initDir[0] = dir[0] * cosT + dir[1] * sinT;
						initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
						grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index

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

						//����;
						SortProj<T>(grid);
						pdist = hypot((static_cast<T>(curIdx) -cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(sliceIdx) -cntImgY)*Img.m_Step.y - sour[1]);
						coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area) / (pdist * FanGeo.m_DetStp);
						summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];

					} //End for one slice
				}// End all slices
				dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
			}
			else
			{
				summ = 0;
				for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
				{
					coord = (sliceIdx - cntImgX + 0.5)* Img.m_Step.x;
					//�������λ�ù����ཻ�����;
					maxC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
					minC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
					if (maxC < minC)
					{
						pdist = minC;
						minC = maxC;
						maxC = pdist;
					}
					curIdx = int(minC / Img.m_Step.y + cntImgY) - extraIdx;
					maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + extraIdx;

					if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0)
					{
						continue;
					}

					if (curIdx < 0)
					{
						curIdx = 0;
					}

					if (maxIdx > Img.m_Reso.y - 1)
					{
						maxIdx = Img.m_Reso.y - 1;
					}

					for (; curIdx <= maxIdx; ++curIdx)
					{
						//��grid
						grid[0][0] = coord - Img.m_Step.x;
						grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
						grid[1][0] = coord;
						grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
						grid[2][0] = coord - Img.m_Step.x;
						grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
						grid[3][0] = coord;
						grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;

						//�����ĸ����Ӧ��ǰ��det index
						dir[0] = grid[0][0] - sour[0];
						dir[1] = grid[0][1] - sour[1];
						len = hypot(dir[0], dir[1]);
						dir[0] /= len;
						dir[1] /= len;
						dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
						dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
						initDir[0] = dir[0] * cosT + dir[1] * sinT;
						initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
						grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index

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

						//����;
						SortProj<T>(grid);
						pdist = hypot((static_cast<T>(sliceIdx) -cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(curIdx) -cntImgY)*Img.m_Step.y - sour[1]);
						coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
						coef = coef / (pdist * FanGeo.m_DetStp);
						summ += coef * dimg[curIdx* Img.m_Reso.x + sliceIdx];
						//summ += *coef;
					} //End for one slice
				}// End all slices
				dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
			}
		}
	}
}
void proj_AIM_CPU_OPENMP(float* prj, float* img, const FanEAGeo& FanGeo, const Image& Img)
{
	proj_AIM_CPU_OPENMP_temp<float>(prj, img, FanGeo, Img);
}
void proj_AIM_CPU_OPENMP(double* prj, double* img, const FanEAGeo& FanGeo, const Image& Img)
{
	proj_AIM_CPU_OPENMP_temp<double>(prj, img, FanGeo, Img);
}



template<typename T>
void bakproj_AIM_CPU_temp(T* dprj, T* dimg, cuint angIdx, const FanEAGeo& FanGeo, const Image& Img)
{
	unsigned int xi(0), yi(0);
	const T dx = Img.m_Step.x;
	const T dy = Img.m_Step.y;
	T cntImgX = static_cast<T>(Img.m_Reso.x - 1.0) * 0.5f;
	T cntImgY = static_cast<T>(Img.m_Reso.y - 1.0) * 0.5f;
	T grid[4][3]; //Storing all four pixels positions and the corresponding index of the detector.
	T cosT = cos(FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp);
	T sinT = sin(FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp);

	//Source Position
	T sour[2];
	sour[0] = -FanGeo.m_S2O*sinT;
	sour[1] = FanGeo.m_S2O*cosT;
	T dir[4][2];//direction of four lights
	T det[4][2]; //current detector position coordinate
	T initDet[4][2]; //original detector position coordinate
	T beta[4];
	int minDetIdx;//��Сdet index
	int maxDetIdx;//���det index
	int curDetIdx;
	T lengths; //Storing the length of the source;

	//Ϊ�������߱ߵļ���;
	T pangle;
	//T temp;
	T SVA[3], SVB[3];
	const T area = dx * dy;
	T pdist; //��Դ�����Ԫ�صľ���;
	T summ;
	T coef;
	unsigned int imgIdx;
	unsigned int prjIdx;
	T initX, initY;
	T curX, curY;
	T direct[2];
	T legth;
	for (yi = 0; yi != Img.m_Reso.y; ++yi)
	{
		for (xi = 0; xi != Img.m_Reso.x; ++xi)
		{
			/*for (yi = 69; yi != 70; ++yi)
			{
			for (xi = 126; xi != 127; ++xi)
			{*/
			imgIdx = yi * Img.m_Reso.x + xi;

			////���㵱ǰ���ص��ĸ�λ��
			//////////////////////////////////////////////////
			////            ^y
			////        P3  |---------| P2
			////            |         |
			////            |    +    |
			////            |         |
			////          P0|---------|----------->x
			////                       P1
			////Fetch the four points of the pixel and their projection positions
			grid[0][0] = (xi - cntImgX - 0.5f)*dx;
			grid[0][1] = (yi - cntImgY - 0.5f)*dy;

			grid[1][0] = (xi - cntImgX + 0.5f)*dx;
			grid[1][1] = (yi - cntImgY - 0.5f)*dy;

			grid[2][0] = (xi - cntImgX + 0.5f)*dx;
			grid[2][1] = (yi - cntImgY + 0.5f)*dy;

			grid[3][0] = (xi - cntImgX - 0.5f)*dx;
			grid[3][1] = (yi - cntImgY + 0.5f)*dy;


			dir[0][0] = grid[0][0] - sour[0];
			dir[0][1] = grid[0][1] - sour[1];
			dir[1][0] = grid[1][0] - sour[0];
			dir[1][1] = grid[1][1] - sour[1];
			dir[2][0] = grid[2][0] - sour[0];
			dir[2][1] = grid[2][1] - sour[1];
			dir[3][0] = grid[3][0] - sour[0];
			dir[3][1] = grid[3][1] - sour[1];

			lengths = sqrt(dir[0][0] * dir[0][0] + dir[0][1] * dir[0][1]);
			dir[0][0] /= lengths;
			dir[0][1] /= lengths;
			lengths = sqrt(dir[1][0] * dir[1][0] + dir[1][1] * dir[1][1]);
			dir[1][0] /= lengths;
			dir[1][1] /= lengths;
			lengths = sqrt(dir[2][0] * dir[2][0] + dir[2][1] * dir[2][1]);
			dir[2][0] /= lengths;
			dir[2][1] /= lengths;
			lengths = sqrt(dir[3][0] * dir[3][0] + dir[3][1] * dir[3][1]);
			dir[3][0] /= lengths;
			dir[3][1] /= lengths;


			////���㵱ǰ���Ӧ���ĸ�̽����λ�����
			det[0][0] = sour[0] + dir[0][0] * FanGeo.m_S2D;
			det[0][1] = sour[1] + dir[0][1] * FanGeo.m_S2D;
			det[1][0] = sour[0] + dir[1][0] * FanGeo.m_S2D;
			det[1][1] = sour[1] + dir[1][1] * FanGeo.m_S2D;
			det[2][0] = sour[0] + dir[2][0] * FanGeo.m_S2D;
			det[2][1] = sour[1] + dir[2][1] * FanGeo.m_S2D;
			det[3][0] = sour[0] + dir[3][0] * FanGeo.m_S2D;
			det[3][1] = sour[1] + dir[3][1] * FanGeo.m_S2D;


			////��������㷶Χ���������ƫ�ƽǶ�;
			////1 ���Ƚ�detector��ת��ȥ
			initDet[0][0] = det[0][0] * cosT + det[0][1] * sinT;
			initDet[0][1] = -det[0][0] * sinT + det[0][1] * cosT;
			initDet[1][0] = det[1][0] * cosT + det[1][1] * sinT;
			initDet[1][1] = -det[1][0] * sinT + det[1][1] * cosT;
			initDet[2][0] = det[2][0] * cosT + det[2][1] * sinT;
			initDet[2][1] = -det[2][0] * sinT + det[2][1] * cosT;
			initDet[3][0] = det[3][0] * cosT + det[3][1] * sinT;
			initDet[3][1] = -det[3][0] * sinT + det[3][1] * cosT;

			////2 ��������Ƕ�; �����ң���С����;
			beta[0] = atan(-initDet[0][0] / (initDet[0][1] - FanGeo.m_S2O));  //��������ʽ����Ҫ��;
			beta[1] = atan(-initDet[1][0] / (initDet[1][1] - FanGeo.m_S2O));
			beta[2] = atan(-initDet[2][0] / (initDet[2][1] - FanGeo.m_S2O));
			beta[3] = atan(-initDet[3][0] / (initDet[3][1] - FanGeo.m_S2O));
			//dimg[imgIdx] = beta[0];
			////3 ��ݽǶ���index;
			grid[0][2] = beta[0] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
			grid[1][2] = beta[1] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
			grid[2][2] = beta[2] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
			grid[3][2] = beta[3] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

			/////����;
			SortProj<T>(grid);

			////��С���index
			minDetIdx = int(grid[0][2]);
			maxDetIdx = ceil(grid[3][2]);

			minDetIdx = (minDetIdx < 0) ? 0 : minDetIdx;
			maxDetIdx = (maxDetIdx > static_cast<int>(FanGeo.m_DetN)) ? static_cast<int>(FanGeo.m_DetN) : maxDetIdx;

			pdist = hypot((xi - cntImgX)*dx - sour[0], (yi - cntImgY)*dy - sour[1]);
			//dimg[imgIdx] = pdist;
			summ = 0;
			for (curDetIdx = minDetIdx; curDetIdx < maxDetIdx; ++curDetIdx)
			{
				//�����߱߼���;
				pangle = (curDetIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp; //С�߱�ƫ�ƽ�
				initY = -cos(pangle) * FanGeo.m_S2D + FanGeo.m_S2O;
				initX = sin(pangle) * FanGeo.m_S2D;
				curX = initX * cosT - initY * sinT;
				curY = initX * sinT + initY * cosT;
				direct[0] = curX - sour[0];
				direct[1] = curY - sour[1];
				legth = sqrt(direct[0] * direct[0] + direct[1] * direct[1]);
				SVA[0] = direct[0] / legth;
				SVA[1] = direct[1] / legth;
				SVA[2] = static_cast<T>(curDetIdx);


				pangle = (curDetIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp + FanGeo.m_DetStp; //��߱�ƫ�ƽ�
				initY = -cos(pangle) * FanGeo.m_S2D + FanGeo.m_S2O;
				initX = sin(pangle) * FanGeo.m_S2D;
				curX = initX * cosT - initY * sinT;
				curY = initX * sinT + initY * cosT;
				direct[0] = curX - sour[0];
				direct[1] = curY - sour[1];
				legth = sqrt(direct[0] * direct[0] + direct[1] * direct[1]);
				SVB[0] = direct[0] / legth;
				SVB[1] = direct[1] / legth;
				SVB[2] = static_cast<T>(curDetIdx) +1;

				//Compute the weighting coefficient for a special projection data
				coef = ComputeCoefficient(grid, SVA, SVB, sour, area);
				coef = coef / (pdist * FanGeo.m_DetStp);
				prjIdx = angIdx* FanGeo.m_DetN + curDetIdx;
				summ += dprj[prjIdx] * coef;
			}
			dimg[imgIdx] = summ;
		}
	}
}
void bakproj_AIM_CPU(float* dprj, float* dimg, cuint angIdx, const FanEAGeo& FanGeo, const Image& Img)
{
	bakproj_AIM_CPU_temp<float>(dprj, dimg, angIdx, FanGeo, Img);
}
void bakproj_AIM_CPU(double* dprj, double* dimg, cuint angIdx, const FanEAGeo& FanGeo, const Image& Img)
{
	bakproj_AIM_CPU_temp<double>(dprj, dimg, angIdx, FanGeo, Img);
}



template<typename T>
void bakproj_AIM_CPU_OPENMP_temp(T* dprj, T* dimg, cuint angIdx, const FanEAGeo& FanGeo, const Image& Img)
{
	int yi(0);
#pragma omp parallel for
	for (yi = 0; yi < Img.m_Reso.y; ++yi)
	{
		unsigned int xi(0);
		const T dx = Img.m_Step.x;
		const T dy = Img.m_Step.y;
		T cntImgX = static_cast<T>(Img.m_Reso.x - 1.0) * 0.5f;
		T cntImgY = static_cast<T>(Img.m_Reso.y - 1.0) * 0.5f;
		T grid[4][3]; //Storing all four pixels positions and the corresponding index of the detector.
		T cosT = cos(FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp);
		T sinT = sin(FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp);

		//Source Position
		T sour[2];
		sour[0] = -FanGeo.m_S2O*sinT;
		sour[1] = FanGeo.m_S2O*cosT;
		T dir[4][2];//direction of four lights
		T det[4][2]; //current detector position coordinate
		T initDet[4][2]; //original detector position coordinate
		T beta[4];
		int minDetIdx;//��Сdet index
		int maxDetIdx;//���det index
		int curDetIdx;
		T lengths; //Storing the length of the source;

		//Ϊ�������߱ߵļ���;
		T pangle;
		//T temp;
		T SVA[3], SVB[3];
		const T area = dx * dy;
		T pdist; //��Դ�����Ԫ�صľ���;
		T summ;
		T coef;
		unsigned int imgIdx;
		unsigned int prjIdx;
		T initX, initY;
		T curX, curY;
		T direct[2];
		T legth;
		for (xi = 0; xi != Img.m_Reso.x; ++xi)
		{
			imgIdx = yi * Img.m_Reso.x + xi;

			////���㵱ǰ���ص��ĸ�λ��
			//////////////////////////////////////////////////
			////            ^y
			////        P3  |---------| P2
			////            |         |
			////            |    +    |
			////            |         |
			////          P0|---------|----------->x
			////                       P1
			////Fetch the four points of the pixel and their projection positions
			grid[0][0] = (xi - cntImgX - 0.5f)*dx;
			grid[0][1] = (yi - cntImgY - 0.5f)*dy;

			grid[1][0] = (xi - cntImgX + 0.5f)*dx;
			grid[1][1] = (yi - cntImgY - 0.5f)*dy;

			grid[2][0] = (xi - cntImgX + 0.5f)*dx;
			grid[2][1] = (yi - cntImgY + 0.5f)*dy;

			grid[3][0] = (xi - cntImgX - 0.5f)*dx;
			grid[3][1] = (yi - cntImgY + 0.5f)*dy;


			dir[0][0] = grid[0][0] - sour[0];
			dir[0][1] = grid[0][1] - sour[1];
			dir[1][0] = grid[1][0] - sour[0];
			dir[1][1] = grid[1][1] - sour[1];
			dir[2][0] = grid[2][0] - sour[0];
			dir[2][1] = grid[2][1] - sour[1];
			dir[3][0] = grid[3][0] - sour[0];
			dir[3][1] = grid[3][1] - sour[1];

			lengths = sqrt(dir[0][0] * dir[0][0] + dir[0][1] * dir[0][1]);
			dir[0][0] /= lengths;
			dir[0][1] /= lengths;
			lengths = sqrt(dir[1][0] * dir[1][0] + dir[1][1] * dir[1][1]);
			dir[1][0] /= lengths;
			dir[1][1] /= lengths;
			lengths = sqrt(dir[2][0] * dir[2][0] + dir[2][1] * dir[2][1]);
			dir[2][0] /= lengths;
			dir[2][1] /= lengths;
			lengths = sqrt(dir[3][0] * dir[3][0] + dir[3][1] * dir[3][1]);
			dir[3][0] /= lengths;
			dir[3][1] /= lengths;


			////���㵱ǰ���Ӧ���ĸ�̽����λ�����
			det[0][0] = sour[0] + dir[0][0] * FanGeo.m_S2D;
			det[0][1] = sour[1] + dir[0][1] * FanGeo.m_S2D;
			det[1][0] = sour[0] + dir[1][0] * FanGeo.m_S2D;
			det[1][1] = sour[1] + dir[1][1] * FanGeo.m_S2D;
			det[2][0] = sour[0] + dir[2][0] * FanGeo.m_S2D;
			det[2][1] = sour[1] + dir[2][1] * FanGeo.m_S2D;
			det[3][0] = sour[0] + dir[3][0] * FanGeo.m_S2D;
			det[3][1] = sour[1] + dir[3][1] * FanGeo.m_S2D;


			////��������㷶Χ���������ƫ�ƽǶ�;
			////1 ���Ƚ�detector��ת��ȥ
			initDet[0][0] = det[0][0] * cosT + det[0][1] * sinT;
			initDet[0][1] = -det[0][0] * sinT + det[0][1] * cosT;
			initDet[1][0] = det[1][0] * cosT + det[1][1] * sinT;
			initDet[1][1] = -det[1][0] * sinT + det[1][1] * cosT;
			initDet[2][0] = det[2][0] * cosT + det[2][1] * sinT;
			initDet[2][1] = -det[2][0] * sinT + det[2][1] * cosT;
			initDet[3][0] = det[3][0] * cosT + det[3][1] * sinT;
			initDet[3][1] = -det[3][0] * sinT + det[3][1] * cosT;

			////2 ��������Ƕ�; �����ң���С����;
			beta[0] = atan(-initDet[0][0] / (initDet[0][1] - FanGeo.m_S2O));  //��������ʽ����Ҫ��;
			beta[1] = atan(-initDet[1][0] / (initDet[1][1] - FanGeo.m_S2O));
			beta[2] = atan(-initDet[2][0] / (initDet[2][1] - FanGeo.m_S2O));
			beta[3] = atan(-initDet[3][0] / (initDet[3][1] - FanGeo.m_S2O));
			//dimg[imgIdx] = beta[0];
			////3 ��ݽǶ���index;
			grid[0][2] = beta[0] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
			grid[1][2] = beta[1] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
			grid[2][2] = beta[2] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
			grid[3][2] = beta[3] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

			/////����;
			SortProj<T>(grid);

			////��С���index
			minDetIdx = int(grid[0][2]);
			maxDetIdx = ceil(grid[3][2]);

			minDetIdx = (minDetIdx < 0) ? 0 : minDetIdx;
			maxDetIdx = (maxDetIdx > static_cast<int>(FanGeo.m_DetN)) ? static_cast<int>(FanGeo.m_DetN) : maxDetIdx;

			pdist = hypot((xi - cntImgX)*dx - sour[0], (yi - cntImgY)*dy - sour[1]);
			//dimg[imgIdx] = pdist;
			summ = 0;
			for (curDetIdx = minDetIdx; curDetIdx < maxDetIdx; ++curDetIdx)
			{
				calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, curDetIdx);

				//Compute the weighting coefficient for a special projection data
				coef = ComputeCoefficient(grid, SVA, SVB, sour, area);
				coef = coef / (pdist * FanGeo.m_DetStp);
				prjIdx = angIdx* FanGeo.m_DetN + curDetIdx;
				summ += dprj[prjIdx] * coef;
			}
			dimg[imgIdx] = summ;
		}
	}
}
void bakproj_AIM_CPU_OPENMP(float* dprj, float* dimg, cuint angIdx, const FanEAGeo& FanGeo, const Image& Img)
{
	bakproj_AIM_CPU_OPENMP_temp<float>(dprj, dimg, angIdx, FanGeo, Img);
}
void bakproj_AIM_CPU_OPENMP(double* dprj, double* dimg, cuint angIdx, const FanEAGeo& FanGeo, const Image& Img)
{
	bakproj_AIM_CPU_OPENMP_temp<double>(dprj, dimg, angIdx, FanGeo, Img);
}

template<typename T>
void bakproj_AIM_CPU_OPENMP_temp(T* dprj, T* dimg, const FanEAGeo& FanGeo, const Image& Img)
{
	int yi(0);
	const T dx = Img.m_Step.x;
	const T dy = Img.m_Step.y;
	const T cntImgX = static_cast<T>(Img.m_Reso.x - 1.0) * 0.5f;
	const T cntImgY = static_cast<T>(Img.m_Reso.y - 1.0) * 0.5f;

#pragma omp parallel for
	for (yi = 0; yi < Img.m_Reso.y; ++yi)
	{
		unsigned int angIdx;
		unsigned int xi(0);

		T grid[4][3]; //Storing all four pixels positions and the corresponding index of the detector.
		T cosT;// = cos(FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp);
		T sinT; // = sin(FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp);

		//Source Position
		T sour[2];

		T dir[4][2];//direction of four lights
		T det[4][2]; //current detector position coordinate
		T initDet[4][2]; //original detector position coordinate
		T beta[4];
		int minDetIdx;//��Сdet index
		int maxDetIdx;//���det index
		int curDetIdx;
		T lengths; //Storing the length of the source;

		//Ϊ�������߱ߵļ���;
		T pangle;
		//T temp;
		T SVA[3], SVB[3];
		const T area = dx * dy;
		T pdist; //��Դ�����Ԫ�صľ���;
		T summ;
		T coef;
		unsigned int imgIdx;
		unsigned int prjIdx;
		T initX, initY;
		T curX, curY;
		T direct[2];
		T legth;
		for (xi = 0; xi != Img.m_Reso.x; ++xi)
		{


			imgIdx = yi * Img.m_Reso.x + xi;
			summ = 0;
			for (angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
			{
				cosT = cos(FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp);
				sinT = sin(FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp);
				sour[0] = -FanGeo.m_S2O*sinT;
				sour[1] = FanGeo.m_S2O*cosT;
				////���㵱ǰ���ص��ĸ�λ��
				//////////////////////////////////////////////////
				////            ^y
				////        P3  |---------| P2
				////            |         |
				////            |    +    |
				////            |         |
				////          P0|---------|----------->x
				////                       P1
				////Fetch the four points of the pixel and their projection positions
				grid[0][0] = (xi - cntImgX - 0.5f)*dx;
				grid[0][1] = (yi - cntImgY - 0.5f)*dy;

				grid[1][0] = (xi - cntImgX + 0.5f)*dx;
				grid[1][1] = (yi - cntImgY - 0.5f)*dy;

				grid[2][0] = (xi - cntImgX + 0.5f)*dx;
				grid[2][1] = (yi - cntImgY + 0.5f)*dy;

				grid[3][0] = (xi - cntImgX - 0.5f)*dx;
				grid[3][1] = (yi - cntImgY + 0.5f)*dy;


				dir[0][0] = grid[0][0] - sour[0];
				dir[0][1] = grid[0][1] - sour[1];
				dir[1][0] = grid[1][0] - sour[0];
				dir[1][1] = grid[1][1] - sour[1];
				dir[2][0] = grid[2][0] - sour[0];
				dir[2][1] = grid[2][1] - sour[1];
				dir[3][0] = grid[3][0] - sour[0];
				dir[3][1] = grid[3][1] - sour[1];

				lengths = sqrt(dir[0][0] * dir[0][0] + dir[0][1] * dir[0][1]);
				dir[0][0] /= lengths;
				dir[0][1] /= lengths;
				lengths = sqrt(dir[1][0] * dir[1][0] + dir[1][1] * dir[1][1]);
				dir[1][0] /= lengths;
				dir[1][1] /= lengths;
				lengths = sqrt(dir[2][0] * dir[2][0] + dir[2][1] * dir[2][1]);
				dir[2][0] /= lengths;
				dir[2][1] /= lengths;
				lengths = sqrt(dir[3][0] * dir[3][0] + dir[3][1] * dir[3][1]);
				dir[3][0] /= lengths;
				dir[3][1] /= lengths;


				////���㵱ǰ���Ӧ���ĸ�̽����λ�����
				det[0][0] = sour[0] + dir[0][0] * FanGeo.m_S2D;
				det[0][1] = sour[1] + dir[0][1] * FanGeo.m_S2D;
				det[1][0] = sour[0] + dir[1][0] * FanGeo.m_S2D;
				det[1][1] = sour[1] + dir[1][1] * FanGeo.m_S2D;
				det[2][0] = sour[0] + dir[2][0] * FanGeo.m_S2D;
				det[2][1] = sour[1] + dir[2][1] * FanGeo.m_S2D;
				det[3][0] = sour[0] + dir[3][0] * FanGeo.m_S2D;
				det[3][1] = sour[1] + dir[3][1] * FanGeo.m_S2D;


				////��������㷶Χ���������ƫ�ƽǶ�;
				////1 ���Ƚ�detector��ת��ȥ
				initDet[0][0] = det[0][0] * cosT + det[0][1] * sinT;
				initDet[0][1] = -det[0][0] * sinT + det[0][1] * cosT;
				initDet[1][0] = det[1][0] * cosT + det[1][1] * sinT;
				initDet[1][1] = -det[1][0] * sinT + det[1][1] * cosT;
				initDet[2][0] = det[2][0] * cosT + det[2][1] * sinT;
				initDet[2][1] = -det[2][0] * sinT + det[2][1] * cosT;
				initDet[3][0] = det[3][0] * cosT + det[3][1] * sinT;
				initDet[3][1] = -det[3][0] * sinT + det[3][1] * cosT;

				////2 ��������Ƕ�; �����ң���С����;
				beta[0] = atan(-initDet[0][0] / (initDet[0][1] - FanGeo.m_S2O));  //��������ʽ����Ҫ��;
				beta[1] = atan(-initDet[1][0] / (initDet[1][1] - FanGeo.m_S2O));
				beta[2] = atan(-initDet[2][0] / (initDet[2][1] - FanGeo.m_S2O));
				beta[3] = atan(-initDet[3][0] / (initDet[3][1] - FanGeo.m_S2O));
				//dimg[imgIdx] = beta[0];
				////3 ��ݽǶ���index;
				grid[0][2] = beta[0] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
				grid[1][2] = beta[1] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
				grid[2][2] = beta[2] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
				grid[3][2] = beta[3] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

				/////����;
				SortProj<T>(grid);

				////��С���index
				minDetIdx = int(grid[0][2]);
				maxDetIdx = ceil(grid[3][2]);

				minDetIdx = (minDetIdx < 0) ? 0 : minDetIdx;
				maxDetIdx = (maxDetIdx > static_cast<int>(FanGeo.m_DetN)) ? static_cast<int>(FanGeo.m_DetN) : maxDetIdx;

				pdist = hypot((xi - cntImgX)*dx - sour[0], (yi - cntImgY)*dy - sour[1]);

				for (curDetIdx = minDetIdx; curDetIdx < maxDetIdx; ++curDetIdx)
				{
					calSVASVB(SVA, SVB, sour, cosT, sinT, FanGeo, Img, curDetIdx);

					//Compute the weighting coefficient for a special projection data
					coef = ComputeCoefficient(grid, SVA, SVB, sour, area);
					coef = coef / (pdist * FanGeo.m_DetStp);
					prjIdx = angIdx* FanGeo.m_DetN + curDetIdx;
					summ += dprj[prjIdx] * coef;
				}
			}
			dimg[imgIdx] = summ;
		}
	}
}
void bakproj_AIM_CPU_OPENMP(float* dprj, float* dimg, const FanEAGeo& FanGeo, const Image& Img)
{
	bakproj_AIM_CPU_OPENMP_temp(dprj, dimg, FanGeo, Img);
}
void bakproj_AIM_CPU_OPENMP(double* dprj, double* dimg, const FanEAGeo& FanGeo, const Image& Img)
{
	bakproj_AIM_CPU_OPENMP_temp(dprj, dimg, FanGeo, Img);
}




template<typename T>
void proj_AIM_CPU_temp(T* dprj, T* dimg, const FanEDGeo& FanGeo, const Image& Img)
{
	int angIdx = 0;
	int detIdx = 0;

	const T cntImgX = (static_cast<T>(Img.m_Reso.x) - 1) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
	const T cntImgY = (static_cast<T>(Img.m_Reso.y) - 1) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
	const T area = Img.m_Step.x * Img.m_Step.y;
	T curAng = 0;
	T cosT = 0;
	T sinT = 0;
	T sour[2] = { 0, 0 };
	T SVA[3], SVB[3];
	unsigned intsliceIdx;
	T summ, coord, minC, maxC;
	int maxIdx, curIdx;
	T grid[4][3];
	T dir[2];
	T initDir[2];
	T len, pdist, coef, curDirAng;
	T initP[2]; // initial pixel position
	T xPos;
	int sliceIdx;
	T pangle;

	const int extraV = 8;//if there is significant Dividing line, increase this value properly.

	for (angIdx = 0; angIdx < FanGeo.m_ViwN; angIdx++)
	{
		curAng = FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp;
		while (curAng < 0){ curAng += _TWOPI; }
		while (curAng> _TWOPI){ curAng -= _TWOPI; }
		cosT = cos(curAng);
		sinT = sin(curAng);
		sour[0] = -FanGeo.m_S2O * sinT;
		sour[1] = FanGeo.m_S2O * cosT;

		for (detIdx = 0; detIdx < FanGeo.m_DetN; ++detIdx)
		{
			calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, detIdx);
			pangle = acos(abs(SVA[0] * SVB[0] + SVA[1] * SVB[1]));


			curDirAng = atan2(((detIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp), FanGeo.m_O2D) + curAng;
			while (curDirAng < 0){ curDirAng += _TWOPI; }
			while (curDirAng> _TWOPI){ curDirAng -= _TWOPI; }

			if (curDirAng <= _PI_4 || curDirAng > _7PI_4)
			{
				summ = 0;
				for (sliceIdx = 0; sliceIdx < Img.m_Reso.y; ++sliceIdx)
				{
					coord = (static_cast<T>(sliceIdx) -cntImgY - 0.5) * Img.m_Step.y;
					minC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
					maxC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
					if (maxC < minC)
					{
						pdist = minC;
						minC = maxC;
						maxC = pdist;
					}

					curIdx = int(minC / Img.m_Step.x + cntImgX) - extraV;
					maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + extraV;

					if (curIdx > static_cast<int>(Img.m_Reso.x - 1) || maxIdx < 0)
					{
						continue;
					}
					if (curIdx < 0)
					{
						curIdx = 0;
					}
					if (maxIdx > Img.m_Reso.x - 1)
					{
						maxIdx = Img.m_Reso.x - 1;
					}


					for (; curIdx <= maxIdx; curIdx++)
					{
						grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
						grid[0][1] = coord;
						grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
						grid[1][1] = coord;
						grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
						grid[2][1] = coord + Img.m_Step.y;
						grid[3][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
						grid[3][1] = coord + Img.m_Step.y;


						//�����ĸ����Ӧ��ǰ��det index
						initP[0] = grid[0][0] * cosT + grid[0][1] * sinT;
						initP[1] = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[0][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[1][0] * cosT + grid[1][1] * sinT;
						initP[1] = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[1][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[2][0] * cosT + grid[2][1] * sinT;
						initP[1] = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[2][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[3][0] * cosT + grid[3][1] * sinT;
						initP[1] = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[3][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						//����;
						SortProj<T>(grid);
						pdist = hypot((static_cast<T>(curIdx) -cntImgX)*static_cast<T>(Img.m_Step.x) - sour[0], (static_cast<T>(sliceIdx) -cntImgY)*static_cast<T>(Img.m_Step.y) - sour[1]);
						coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
						coef = coef / (pdist*pangle);
						summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
					}
				}
				dprj[angIdx* FanGeo.m_DetN + detIdx] = summ;
				//continue;
			} //End case 1
			else if (curDirAng > _PI_4 && curDirAng <= _3PI_4)
			{
				summ = 0;
				for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
				{
					coord = (static_cast<T>(sliceIdx) -cntImgX - 0.5)* Img.m_Step.x;
					//�������λ�ù����ཻ�����;
					minC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
					maxC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
					if (maxC < minC)
					{
						pdist = minC;
						minC = maxC;
						maxC = pdist;
					}
					curIdx = int(minC / Img.m_Step.y + cntImgY) - extraV;
					maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + extraV;

					if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0)
					{
						continue;
					}
					if (curIdx < 0)
					{
						curIdx = 0;
					}

					if (maxIdx > Img.m_Reso.y - 1)
					{
						maxIdx = Img.m_Reso.y - 1;
					}

					for (; curIdx <= maxIdx; ++curIdx)
					{
						//��grid
						grid[0][0] = coord;
						grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
						grid[1][0] = coord + Img.m_Step.x;
						grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
						grid[2][0] = coord + Img.m_Step.x;
						grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
						grid[3][0] = coord;
						grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;


						//�����ĸ����Ӧ��ǰ��det index
						initP[0] = grid[0][0] * cosT + grid[0][1] * sinT;
						initP[1] = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[0][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[1][0] * cosT + grid[1][1] * sinT;
						initP[1] = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[1][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[2][0] * cosT + grid[2][1] * sinT;
						initP[1] = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[2][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[3][0] * cosT + grid[3][1] * sinT;
						initP[1] = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[3][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;


						//����;
						SortProj<T>(grid);
						pdist = hypot((static_cast<T>(curIdx) -cntImgX)*static_cast<T>(Img.m_Step.x) - sour[0], (static_cast<T>(sliceIdx) -cntImgY)*static_cast<T>(Img.m_Step.y) - sour[1]);
						coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
						coef = coef / (pdist*pangle);
						summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];

					} //End for one slice
				}// End all slices
				dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
				//continue;
			}// End case 2
			else if (curDirAng > _3PI_4 && curDirAng <= _5PI_4)
			{

				summ = 0;
				for (sliceIdx = 0; sliceIdx < Img.m_Reso.y; ++sliceIdx)
				{
					coord = (static_cast<T>(sliceIdx) -cntImgY + 0.5)* Img.m_Step.y;
					//�������λ�ù����ཻ�����;
					maxC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
					minC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
					if (maxC < minC)
					{
						pdist = minC;
						minC = maxC;
						maxC = pdist;
					}
					curIdx = int(minC / Img.m_Step.x + cntImgX) - extraV;
					maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + extraV;

					if (curIdx > static_cast<int>(Img.m_Reso.x) || maxIdx < 0)
					{
						continue;
					}

					if (curIdx < 0)
					{
						curIdx = 0;
					}

					if (maxIdx > Img.m_Reso.x - 1)
					{
						maxIdx = Img.m_Reso.x - 1;
					}

					for (; curIdx <= maxIdx; ++curIdx)
					{
						//��grid
						grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
						grid[0][1] = coord - Img.m_Step.y;
						grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
						grid[1][1] = coord - Img.m_Step.y;
						grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
						grid[2][1] = coord;
						grid[3][0] = (curIdx - cntImgX - 0.5)* Img.m_Step.x;
						grid[3][1] = coord;

						//�����ĸ����Ӧ��ǰ��det index
						initP[0] = grid[0][0] * cosT + grid[0][1] * sinT;
						initP[1] = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[0][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[1][0] * cosT + grid[1][1] * sinT;
						initP[1] = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[1][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[2][0] * cosT + grid[2][1] * sinT;
						initP[1] = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[2][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[3][0] * cosT + grid[3][1] * sinT;
						initP[1] = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[3][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;


						//����;
						SortProj<T>(grid);
						pdist = hypot((static_cast<T>(curIdx) -cntImgX)*static_cast<T>(Img.m_Step.x) - sour[0], (static_cast<T>(sliceIdx) -cntImgY)*static_cast<T>(Img.m_Step.y) - sour[1]);
						coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
						coef = coef / (pdist*pangle);
						summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];

					} //End for one slice
				}// End all slices
				dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
				//continue;
			}
			else
			{

				summ = 0;
				for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
				{
					coord = (static_cast<T>(sliceIdx) -cntImgX + 0.5)* Img.m_Step.x;
					//�������λ�ù����ཻ�����;
					maxC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
					minC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
					if (maxC < minC)
					{
						pdist = minC;
						minC = maxC;
						maxC = pdist;
					}
					curIdx = int(minC / Img.m_Step.y + cntImgY) - extraV;
					maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + extraV;

					if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0)
					{
						continue;
					}

					if (curIdx < 0)
					{
						curIdx = 0;
					}

					if (maxIdx > Img.m_Reso.y - 1)
					{
						maxIdx = Img.m_Reso.y - 1;
					}

					for (; curIdx <= maxIdx; ++curIdx)
					{
						grid[0][0] = coord - Img.m_Step.x;
						grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
						grid[1][0] = coord;
						grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
						grid[2][0] = coord - Img.m_Step.x;
						grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
						grid[3][0] = coord;
						grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;


						//�����ĸ����Ӧ��ǰ��det index
						initP[0] = grid[0][0] * cosT + grid[0][1] * sinT;
						initP[1] = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[0][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[1][0] * cosT + grid[1][1] * sinT;
						initP[1] = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[1][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[2][0] * cosT + grid[2][1] * sinT;
						initP[1] = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[2][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[3][0] * cosT + grid[3][1] * sinT;
						initP[1] = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[3][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;


						//����;
						SortProj<T>(grid);
						pdist = hypot((static_cast<T>(curIdx) -cntImgX)*static_cast<T>(Img.m_Step.x) - sour[0], (static_cast<T>(sliceIdx) -cntImgY)*static_cast<T>(Img.m_Step.y) - sour[1]);
						coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
						coef = coef / (pdist*pangle);
						summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
						//summ += *coef;
					} //End for one slice
				}// End all slices
				dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
				//continue;
			}
		}
	}
}
void proj_AIM_CPU(float* dprj, float* dimg, const FanEDGeo& FanGeo, const Image& Img)
{
	proj_AIM_CPU_temp<float>(dprj, dimg, FanGeo, Img);
}
void proj_AIM_CPU(double* dprj, double* dimg, const FanEDGeo& FanGeo, const Image& Img)
{
	proj_AIM_CPU_temp<double>(dprj, dimg, FanGeo, Img);
}


template<typename T>
void proj_AIM_CPU_OPENMP_temp(T* dprj, T* dimg, const FanEDGeo& FanGeo, const Image& Img)
{
	int angIdx = 0;

#pragma omp parallel for 
	for (angIdx = 0; angIdx < FanGeo.m_ViwN; angIdx++)
	{
		const T cntImgX = (static_cast<T>(Img.m_Reso.x) - 1) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
		const T cntImgY = (static_cast<T>(Img.m_Reso.y) - 1) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
		const T area = Img.m_Step.x * Img.m_Step.y;
		T curAng = 0;
		T cosT = 0;
		T sinT = 0;
		T sour[2] = { 0, 0 };
		T SVA[3], SVB[3];
		unsigned intsliceIdx;
		T summ, coord, minC, maxC;
		int maxIdx, curIdx;
		T grid[4][3];
		T dir[2];
		T initDir[2];
		T len, pdist, coef, curDirAng;
		T initP[2]; // initial pixel position
		T xPos;
		int sliceIdx;
		T pangle;

		int extraV = 8;
		for (int detIdx = 0; detIdx < FanGeo.m_DetN; ++detIdx)
		{
			curAng = FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp;
			while (curAng < 0){ curAng += _TWOPI; }
			while (curAng> _TWOPI){ curAng -= _TWOPI; }
			cosT = cos(curAng);
			sinT = sin(curAng);
			sour[0] = -FanGeo.m_S2O * sinT;
			sour[1] = FanGeo.m_S2O * cosT;

			calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, detIdx);
			pangle = acos(abs(SVA[0] * SVB[0] + SVA[1] * SVB[1]));


			curDirAng = atan2(((detIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp), FanGeo.m_O2D) + curAng;
			while (curDirAng < 0){ curDirAng += _TWOPI; }
			while (curDirAng> _TWOPI){ curDirAng -= _TWOPI; }
			//dprj[angIdx * FanGeo.m_DetN + detIdx] = curDirAng;
			//continue;
			if (curDirAng <= _PI_4 || curDirAng > _7PI_4)
			{
				summ = 0;
				for (sliceIdx = 0; sliceIdx < Img.m_Reso.y; ++sliceIdx)
				{
					coord = (static_cast<T>(sliceIdx) -cntImgY - 0.5) * Img.m_Step.y;
					minC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
					maxC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
					if (maxC < minC)
					{
						pdist = minC;
						minC = maxC;
						maxC = pdist;
					}
					curIdx = int(minC / Img.m_Step.x + cntImgX) - extraV;
					maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + extraV;

					if (curIdx > static_cast<int>(Img.m_Reso.x - 1) || maxIdx < 0)
					{
						continue;
					}
					if (curIdx < 0)
					{
						curIdx = 0;
					}
					if (maxIdx > Img.m_Reso.x - 1)
					{
						maxIdx = Img.m_Reso.x - 1;
					}

					for (; curIdx <= maxIdx; curIdx++)
					{
						grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
						grid[0][1] = coord;
						grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
						grid[1][1] = coord;
						grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
						grid[2][1] = coord + Img.m_Step.y;
						grid[3][0] = (curIdx - cntImgX - 0.5)* Img.m_Step.x;
						grid[3][1] = coord + Img.m_Step.y;


						//�����ĸ����Ӧ��ǰ��det index
						initP[0] = grid[0][0] * cosT + grid[0][1] * sinT;
						initP[1] = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[0][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[1][0] * cosT + grid[1][1] * sinT;
						initP[1] = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[1][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[2][0] * cosT + grid[2][1] * sinT;
						initP[1] = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[2][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[3][0] * cosT + grid[3][1] * sinT;
						initP[1] = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[3][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						//����;
						SortProj<T>(grid);
						pdist = hypot((static_cast<T>(curIdx) -cntImgX)*static_cast<T>(Img.m_Step.x) - sour[0], (static_cast<T>(sliceIdx) -cntImgY)*static_cast<T>(Img.m_Step.y) - sour[1]);
						coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
						coef = coef / (pdist*pangle);
						summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
					}
				}
				dprj[angIdx* FanGeo.m_DetN + detIdx] = summ;
			} //End case 1
			else if (curDirAng > _PI_4 && curDirAng <= _3PI_4)
			{
				summ = 0;
				for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
				{
					coord = (static_cast<T>(sliceIdx) -cntImgX - 0.5)* Img.m_Step.x;
					//�������λ�ù����ཻ�����;
					minC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
					maxC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
					if (maxC < minC)
					{
						pdist = minC;
						minC = maxC;
						maxC = pdist;
					}
					curIdx = int(minC / Img.m_Step.y + cntImgY) - extraV;
					maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + extraV;

					if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0)
					{
						continue;
					}
					if (curIdx < 0)
					{
						curIdx = 0;
					}

					if (maxIdx > Img.m_Reso.y - 1)
					{
						maxIdx = Img.m_Reso.y - 1;
					}

					for (; curIdx <= maxIdx; ++curIdx)
					{
						//��grid
						grid[0][0] = coord;
						grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
						grid[1][0] = coord + Img.m_Step.x;
						grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
						grid[2][0] = coord + Img.m_Step.x;
						grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
						grid[3][0] = coord;
						grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;


						//�����ĸ����Ӧ��ǰ��det index
						initP[0] = grid[0][0] * cosT + grid[0][1] * sinT;
						initP[1] = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[0][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[1][0] * cosT + grid[1][1] * sinT;
						initP[1] = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[1][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[2][0] * cosT + grid[2][1] * sinT;
						initP[1] = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[2][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[3][0] * cosT + grid[3][1] * sinT;
						initP[1] = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[3][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;


						//����;
						SortProj<T>(grid);
						pdist = hypot((static_cast<T>(curIdx) -cntImgX)*static_cast<T>(Img.m_Step.x) - sour[0], (static_cast<T>(sliceIdx) -cntImgY)*static_cast<T>(Img.m_Step.y) - sour[1]);
						coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
						coef = coef / (pdist*pangle);
						summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];

					} //End for one slice
				}// End all slices
				dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
			}// End case 2
			else if (curDirAng > _3PI_4 && curDirAng <= _5PI_4)
			{

				summ = 0;
				for (sliceIdx = 0; sliceIdx < Img.m_Reso.y; ++sliceIdx)
				{
					coord = (static_cast<T>(sliceIdx) -cntImgY + 0.5)* Img.m_Step.y;
					//�������λ�ù����ཻ�����;
					maxC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
					minC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
					if (maxC < minC)
					{
						pdist = minC;
						minC = maxC;
						maxC = pdist;
					}
					curIdx = int(minC / Img.m_Step.x + cntImgX) - extraV;
					maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + extraV;

					if (curIdx > static_cast<int>(Img.m_Reso.x) || maxIdx < 0)
					{
						continue;
					}

					if (curIdx < 0)
					{
						curIdx = 0;
					}

					if (maxIdx > Img.m_Reso.x - 1)
					{
						maxIdx = Img.m_Reso.x - 1;
					}

					for (; curIdx <= maxIdx; ++curIdx)
					{
						//��grid
						grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
						grid[0][1] = coord - Img.m_Step.y;
						grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
						grid[1][1] = coord - Img.m_Step.y;
						grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
						grid[2][1] = coord;
						grid[3][0] = (curIdx - cntImgX - 0.5)* Img.m_Step.x;
						grid[3][1] = coord;

						//�����ĸ����Ӧ��ǰ��det index
						initP[0] = grid[0][0] * cosT + grid[0][1] * sinT;
						initP[1] = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[0][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[1][0] * cosT + grid[1][1] * sinT;
						initP[1] = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[1][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[2][0] * cosT + grid[2][1] * sinT;
						initP[1] = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[2][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[3][0] * cosT + grid[3][1] * sinT;
						initP[1] = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[3][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;


						//����;
						SortProj<T>(grid);
						pdist = hypot((static_cast<T>(curIdx) -cntImgX)*static_cast<T>(Img.m_Step.x) - sour[0], (static_cast<T>(sliceIdx) -cntImgY)*static_cast<T>(Img.m_Step.y) - sour[1]);
						coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
						coef = coef / (pdist*pangle);
						summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];

					} //End for one slice
				}// End all slices
				dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
			}
			else
			{

				summ = 0;
				for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
				{
					coord = (static_cast<T>(sliceIdx) -cntImgX + 0.5)* Img.m_Step.x;
					//�������λ�ù����ཻ�����;
					maxC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
					minC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
					if (maxC < minC)
					{
						pdist = minC;
						minC = maxC;
						maxC = pdist;
					}
					curIdx = int(minC / Img.m_Step.y + cntImgY) - extraV;
					maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + extraV;

					if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0)
					{
						continue;
					}

					if (curIdx < 0)
					{
						curIdx = 0;
					}

					if (maxIdx > Img.m_Reso.y - 1)
					{
						maxIdx = Img.m_Reso.y - 1;
					}

					for (; curIdx <= maxIdx; ++curIdx)
					{
						grid[0][0] = coord - Img.m_Step.x;
						grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
						grid[1][0] = coord;
						grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
						grid[2][0] = coord - Img.m_Step.x;
						grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
						grid[3][0] = coord;
						grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;


						//�����ĸ����Ӧ��ǰ��det index
						initP[0] = grid[0][0] * cosT + grid[0][1] * sinT;
						initP[1] = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[0][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[1][0] * cosT + grid[1][1] * sinT;
						initP[1] = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[1][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[2][0] * cosT + grid[2][1] * sinT;
						initP[1] = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[2][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

						initP[0] = grid[3][0] * cosT + grid[3][1] * sinT;
						initP[1] = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
						xPos = -initP[0] * FanGeo.m_S2D / initP[1];
						grid[3][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;


						//����;
						SortProj<T>(grid);
						pdist = hypot((static_cast<T>(curIdx) -cntImgX)*static_cast<T>(Img.m_Step.x) - sour[0], (static_cast<T>(sliceIdx) -cntImgY)*static_cast<T>(Img.m_Step.y) - sour[1]);
						coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
						coef = coef / (pdist*pangle);
						summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
						//summ += *coef;
					} //End for one slice
				}// End all slices
				dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
			}
		}
	}
}
void proj_AIM_CPU_OPENMP(float* dprj, float* dimg, const FanEDGeo& FanGeo, const Image& Img)
{
	proj_AIM_CPU_OPENMP_temp<float>(dprj, dimg, FanGeo, Img);
}
void proj_AIM_CPU_OPENMP(double* dprj, double* dimg, const FanEDGeo& FanGeo, const Image& Img)
{
	proj_AIM_CPU_OPENMP_temp<double>(dprj, dimg, FanGeo, Img);
}



template<typename T>
void bakproj_AIM_CPU_temp(T* dprj, T*dimg, cuint angidx, const FanEDGeo& FanGeo, const Image& Img)
{
	const T cosT = cos(FanGeo.m_ViwBeg + angidx * FanGeo.m_ViwStp);
	const T sinT = sin(FanGeo.m_ViwBeg + angidx * FanGeo.m_ViwStp);
	const T cntImgX = static_cast<T>(Img.m_Reso.x - 1.0) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
	const T cntImgY = static_cast<T>(Img.m_Reso.y - 1.0) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
	const T area = Img.m_Step.x * Img.m_Step.y;

	T sour[2];
	sour[0] = -FanGeo.m_S2O * sinT;
	sour[1] = FanGeo.m_S2O * cosT;
	int xi(0), yi(0);
	T grid[4][3];
	int minDetIdx, maxDetIdx;
	T pdist;
	T SVA[3], SVB[3];
	T vv;
	T coef;
	T summ;
	for (yi = 0; yi != Img.m_Reso.y; ++yi)
	{
		for (xi = 0; xi != Img.m_Reso.x; ++xi)
		{
			summ = 0;
			grid[0][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
			grid[0][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
			pdist = grid[0][0] * cosT + grid[0][1] * sinT;
			vv = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
			grid[0][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

			grid[1][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
			grid[1][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
			pdist = grid[1][0] * cosT + grid[1][1] * sinT;
			vv = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
			grid[1][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

			grid[2][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
			grid[2][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
			pdist = grid[2][0] * cosT + grid[2][1] * sinT;
			vv = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
			grid[2][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

			grid[3][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
			grid[3][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
			pdist = grid[3][0] * cosT + grid[3][1] * sinT;
			vv = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
			grid[3][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

			SortProj<T>(grid);

			minDetIdx = static_cast<int>(grid[0][2]);
			maxDetIdx = (int(ceil(grid[3][2])));
			minDetIdx = (minDetIdx < 0) ? 0 : minDetIdx;
			maxDetIdx = (maxDetIdx > static_cast<int>(FanGeo.m_DetN)) ? static_cast<int>(FanGeo.m_DetN) : maxDetIdx;
			pdist = (hypot((xi - cntImgX)*Img.m_Step.x - sour[0], (yi - cntImgY)*Img.m_Step.y - sour[1]));

			for (; minDetIdx < maxDetIdx; ++minDetIdx)
			{
				calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, minDetIdx);
				vv = acos(abs(SVA[0] * SVB[0] + SVA[1] * SVB[1]));
				coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
				coef = coef / (pdist * vv);
				summ += dprj[angidx * FanGeo.m_DetN + minDetIdx] * coef;

			}
			dimg[yi* Img.m_Reso.x + xi] = summ;
		}//End for xi
	}// End for yi

}
void bakproj_AIM_CPU(float* dprj, float* dimg, cuint angIdx, const FanEDGeo& FanGeo, const Image& Img)
{
	bakproj_AIM_CPU_temp<float>(dprj, dimg, angIdx, FanGeo, Img);
}
void bakproj_AIM_CPU(double* dprj, double* dimg, cuint angIdx, const FanEDGeo& FanGeo, const Image& Img)
{
	bakproj_AIM_CPU_temp<double>(dprj, dimg, angIdx, FanGeo, Img);
}


template<typename T>
void bakproj_AIM_CPU_temp(T* dprj, T*dimg, const FanEDGeo& FanGeo, const Image& Img)
{
	T cosT;// = cos(FanGeo.m_ViwBeg + angidx * FanGeo.m_ViwStp);
	T sinT;// = sin(FanGeo.m_ViwBeg + angidx * FanGeo.m_ViwStp);
	const T cntImgX = static_cast<T>(Img.m_Reso.x - 1.0) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
	const T cntImgY = static_cast<T>(Img.m_Reso.y - 1.0) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
	const T area = Img.m_Step.x * Img.m_Step.y;

	T sour[2];
	//sour[0] = -FanGeo.m_S2O * sinT;
	//sour[1] = FanGeo.m_S2O * cosT;

	int xi(0), yi(0);
	T grid[4][3];

	int minDetIdx, maxDetIdx;
	T pdist;
	T SVA[3], SVB[3];
	T vv;
	T coef;
	T summ;
	int angidx = 0;
	for (yi = 0; yi != Img.m_Reso.y; ++yi)
	{
		for (xi = 0; xi != Img.m_Reso.x; ++xi)
		{

			summ = 0;
			for (angidx = 0; angidx != FanGeo.m_ViwN; ++angidx)
			{
				cosT = cos(FanGeo.m_ViwBeg + angidx * FanGeo.m_ViwStp);
				sinT = sin(FanGeo.m_ViwBeg + angidx * FanGeo.m_ViwStp);

				sour[0] = -FanGeo.m_S2O * sinT;
				sour[1] = FanGeo.m_S2O * cosT;

				grid[0][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
				grid[0][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
				pdist = grid[0][0] * cosT + grid[0][1] * sinT;
				vv = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
				grid[0][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

				grid[1][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
				grid[1][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
				pdist = grid[1][0] * cosT + grid[1][1] * sinT;
				vv = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
				grid[1][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

				grid[2][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
				grid[2][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
				pdist = grid[2][0] * cosT + grid[2][1] * sinT;
				vv = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
				grid[2][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

				grid[3][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
				grid[3][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
				pdist = grid[3][0] * cosT + grid[3][1] * sinT;
				vv = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
				grid[3][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;



				//����
				SortProj<T>(grid);

				minDetIdx = static_cast<int>(grid[0][2]);
				maxDetIdx = (int(ceil(grid[3][2])));
				minDetIdx = (minDetIdx < 0) ? 0 : minDetIdx;
				maxDetIdx = (maxDetIdx > static_cast<int>(FanGeo.m_DetN)) ? static_cast<int>(FanGeo.m_DetN) : maxDetIdx;
				pdist = (hypot((xi - cntImgX)*Img.m_Step.x - sour[0], (yi - cntImgY)*Img.m_Step.y - sour[1]));

				for (; minDetIdx < maxDetIdx; ++minDetIdx)
				{
					calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, minDetIdx);
					vv = acos(abs(SVA[0] * SVB[0] + SVA[1] * SVB[1]));
					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
					coef = coef / (pdist * vv);
					summ += dprj[angidx * FanGeo.m_DetN + minDetIdx] * coef;

				}
			}


			dimg[yi* Img.m_Reso.x + xi] = summ;
		}//End for xi
	}// End for yi

}
void bakproj_AIM_CPU(float* dprj, float* dimg, const FanEDGeo& FanGeo, const Image& Img)
{
	bakproj_AIM_CPU_temp<float>(dprj, dimg, FanGeo, Img);
}
void bakproj_AIM_CPU(double* dprj, double* dimg, const FanEDGeo& FanGeo, const Image& Img)
{
	bakproj_AIM_CPU_temp<double>(dprj, dimg, FanGeo, Img);
}




template<typename T>
void bakproj_AIM_CPU_OPENMP_temp(T* dprj, T*dimg, const FanEDGeo& FanGeo, const Image& Img)
{
	const T cntImgX = static_cast<T>(Img.m_Reso.x - 1.0) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
	const T cntImgY = static_cast<T>(Img.m_Reso.y - 1.0) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
	const T area = Img.m_Step.x * Img.m_Step.y;

	int yi(0);
#pragma omp parallel for
	for (yi = 0; yi < Img.m_Reso.y; yi++)
	{
		T cosT;// = cos(FanGeo.m_ViwBeg + angidx * FanGeo.m_ViwStp);
		T sinT;// = sin(FanGeo.m_ViwBeg + angidx * FanGeo.m_ViwStp);

		T sour[2];

		int xi(0);
		T grid[4][3];

		int minDetIdx, maxDetIdx;
		T pdist;
		T SVA[3], SVB[3];
		T vv;
		T coef;
		T summ;
		int angidx = 0;
		for (xi = 0; xi != Img.m_Reso.x; ++xi)
		{

			summ = 0;
			for (angidx = 0; angidx != FanGeo.m_ViwN; ++angidx)
			{
				cosT = cos(FanGeo.m_ViwBeg + angidx * FanGeo.m_ViwStp);
				sinT = sin(FanGeo.m_ViwBeg + angidx * FanGeo.m_ViwStp);

				sour[0] = -FanGeo.m_S2O * sinT;
				sour[1] = FanGeo.m_S2O * cosT;

				grid[0][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
				grid[0][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
				pdist = grid[0][0] * cosT + grid[0][1] * sinT;
				vv = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
				grid[0][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

				grid[1][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
				grid[1][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
				pdist = grid[1][0] * cosT + grid[1][1] * sinT;
				vv = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
				grid[1][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

				grid[2][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
				grid[2][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
				pdist = grid[2][0] * cosT + grid[2][1] * sinT;
				vv = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
				grid[2][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

				grid[3][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
				grid[3][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
				pdist = grid[3][0] * cosT + grid[3][1] * sinT;
				vv = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
				grid[3][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;


				//����
				SortProj<T>(grid);

				minDetIdx = static_cast<int>(grid[0][2]);
				maxDetIdx = (int(ceil(grid[3][2])));
				minDetIdx = (minDetIdx < 0) ? 0 : minDetIdx;
				maxDetIdx = (maxDetIdx > static_cast<int>(FanGeo.m_DetN)) ? static_cast<int>(FanGeo.m_DetN) : maxDetIdx;
				pdist = (hypot((xi - cntImgX)*Img.m_Step.x - sour[0], (yi - cntImgY)*Img.m_Step.y - sour[1]));

				for (; minDetIdx < maxDetIdx; ++minDetIdx)
				{
					calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, minDetIdx);
					vv = acos(abs(SVA[0] * SVB[0] + SVA[1] * SVB[1]));
					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
					coef = coef / (pdist * vv);
					summ += dprj[angidx * FanGeo.m_DetN + minDetIdx] * coef;

				}
			}


			dimg[yi* Img.m_Reso.x + xi] = summ;
		}//End for xi
	}// End for yi

}
void bakproj_AIM_CPU_OPENMP(float* dprj, float* dimg, const FanEDGeo& FanGeo, const Image& Img)
{
	bakproj_AIM_CPU_OPENMP_temp<float>(dprj, dimg, FanGeo, Img);
}

void bakproj_AIM_CPU_OPENMP(double* dprj, double* dimg, const FanEDGeo& FanGeo, const Image& Img)
{
	bakproj_AIM_CPU_OPENMP_temp<double>(dprj, dimg, FanGeo, Img);
}


template<typename T>
void bakproj_AIM_CPU_OPENMP_temp(T* dprj, T*dimg, cuint angidx, const FanEDGeo& FanGeo, const Image& Img)
{
	const T cosT = cos(FanGeo.m_ViwBeg + angidx * FanGeo.m_ViwStp);
	const T sinT = sin(FanGeo.m_ViwBeg + angidx * FanGeo.m_ViwStp);
	const T cntImgX = static_cast<T>(Img.m_Reso.x - 1.0) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
	const T cntImgY = static_cast<T>(Img.m_Reso.y - 1.0) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
	const T area = Img.m_Step.x * Img.m_Step.y;
	int yi;

#pragma omp parallel for
	for (yi = 0; yi < Img.m_Reso.y; yi++)
	{
		T sour[2];
		sour[0] = -FanGeo.m_S2O * sinT;
		sour[1] = FanGeo.m_S2O * cosT;

		int xi(0);
		T grid[4][3];

		int minDetIdx, maxDetIdx;
		T pdist;
		T SVA[3], SVB[3];
		T vv;
		T coef;
		T summ;
		for (xi = 0; xi != Img.m_Reso.x; ++xi)
		{
			summ = 0;
			grid[0][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
			grid[0][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
			pdist = grid[0][0] * cosT + grid[0][1] * sinT;
			vv = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
			grid[0][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

			grid[1][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
			grid[1][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
			pdist = grid[1][0] * cosT + grid[1][1] * sinT;
			vv = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
			grid[1][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

			grid[2][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
			grid[2][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
			pdist = grid[2][0] * cosT + grid[2][1] * sinT;
			vv = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
			grid[2][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

			grid[3][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
			grid[3][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
			pdist = grid[3][0] * cosT + grid[3][1] * sinT;
			vv = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
			grid[3][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;

			//����
			SortProj<T>(grid);

			minDetIdx = static_cast<int>(grid[0][2]);
			maxDetIdx = (int(ceil(grid[3][2])));
			minDetIdx = (minDetIdx < 0) ? 0 : minDetIdx;
			maxDetIdx = (maxDetIdx > static_cast<int>(FanGeo.m_DetN)) ? static_cast<int>(FanGeo.m_DetN) : maxDetIdx;
			pdist = (hypot((xi - cntImgX)*Img.m_Step.x - sour[0], (yi - cntImgY)*Img.m_Step.y - sour[1]));

			for (; minDetIdx < maxDetIdx; ++minDetIdx)
			{
				calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, minDetIdx);
				vv = acos(abs(SVA[0] * SVB[0] + SVA[1] * SVB[1]));
				coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
				coef = coef / (pdist * vv);
				summ += dprj[angidx * FanGeo.m_DetN + minDetIdx] * coef;

			}
			dimg[yi* Img.m_Reso.x + xi] = summ;
		}//End for xi
	}// End for yi

}
void bakproj_AIM_CPU_OPENMP(float* dprj, float* dimg, cuint angIdx, const FanEDGeo& FanGeo, const Image& Img)
{
	bakproj_AIM_CPU_OPENMP_temp<float>(dprj, dimg, angIdx, FanGeo, Img);
}
void bakproj_AIM_CPU_OPENMP(double* dprj, double* dimg, cuint angIdx, const FanEDGeo& FanGeo, const Image& Img)
{
	bakproj_AIM_CPU_OPENMP_temp<double>(dprj, dimg, angIdx, FanGeo, Img);
}
