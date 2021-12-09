//#include "projbackproj.hpp"
//#include "utilities.hpp"
//#include "useful.hpp"
//#include "cudaCheckReturner.h"
//
//#ifndef DEMO10_SLICESNUM
//#define DEMO10_SLICESNUM 64
//#endif
//
//#ifndef DEMO11_SLICESNUM
//#define DEMO11_SLICESNUM 22
//#endif
//
//
//
//
//
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
//
///// \brief the rotation of the 2D vector according to cos(T) and sin(T)
///// \param p original vector
///// \param cosT cos(T)
///// \param sinT sin(T)
//inline __host__	__device__ double2 rotation(const double2& p, const double& cosT, const double& sinT)
//{
//	double2 curP;
//	curP.x = p.x * cosT - p.y * sinT;
//	curP.y = p.x * sinT + p.y * cosT;
//	return curP;
//}
//
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
//
//
///// \brief the rotation of the 3D vector according to cos(T) and sin(T)
///// \param p original vector
///// \param cosT cos(T)
///// \param sinT sin(T)
//inline __host__ __device__ double3 rotation(const double3& p, const double& cosT, const double& sinT)
//{
//	double3 curP;
//	curP.x = p.x * cosT - p.y * sinT;
//	curP.y = p.x * sinT + p.y * cosT;
//	curP.z = p.z;
//	return curP;
//}
//
//
//__global__ void _proj_Ker(float* dimg, float* donp, float* draw, cuint angIdx, const bool angDir, const FanEAGeo FanGeo, const Image Img)
//{
//	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
//	if (detId < FanGeo.m_DetN)
//	{
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//
//		//Current rotation angle;
//		float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//Current source position, assuming initial position is on the positive Y axis
//		Ray2D ray;
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//		float ang(0); //bias angle from the center of the Fan Beams
//
//		float2 curDetPos; //the current detector element position;
//		float totWeight(0);
//
//
//		// judging two indices addressing method, from large to small or small to large
//		if (angDir)
//		{
//			ang = ((float)detId - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//		}
//		else
//		{
//			ang = ((FanGeo.m_DetN - 1 - (float)detId) - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//		}
//		// current detector element position
//		curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//
//		// x-ray direction;
//		ray.d = normalize(curDetPos - ray.o);
//
//		totWeight = 0;
//		donp[detId] = draw[angIdx * FanGeo.m_DetN + detId] - calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
//			MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg, &totWeight);
//		if (!IS_ZERO(totWeight))
//		{
//			donp[detId] /= totWeight;
//		}
//		else
//		{
//			donp[detId] = 0;
//		}
//	}
//}
////////////////////////////////////////////////////////////////////////////
//// This function is the projection difference with weighting;
//// dimg: image to be projected
//// donp: one angle projection data difference with weighting
//// draw: raw projection data
//// angIdx: angle index
//// angDir: detector from positive to negative or verse
//// FanGeo: Fan Beam geometry
//// Img: Image configuration;
//// blk: it is better to be 1024, 1
//// gid: correspondingly configurated
////////////////////////////////////////////////////////////////////////////
//void proj(float* dimg, float* donp, float* draw, cuint angIdx, const bool angDir, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker << <gid, blk >> >(dimg, donp, draw, angIdx, angDir, FanGeo, Img);
//}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//__global__ void _proj_Ker(float* dimg, float* donp, float* draw, cuint angIdx, const bool angDir, const FanEDGeo FanGeo, const Image Img)
//{
//	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
//	if (detId < FanGeo.m_DetN)
//	{
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//
//		//Current rotation angle;
//		float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//Current source position, assuming initial position is on the positive Y axis
//		Ray2D ray;
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//		float ang(0); //bias angle from the center of the Fan Beams
//
//		float2 curDetPos; //the current detector element position;
//		float totWeight(0);
//
//
//		// judging two indices addressing method, from large to small or small to large
//		if (angDir)
//		{
//			ang = ((float)detId - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//		}
//		else
//		{
//			ang = ((FanGeo.m_DetN - 1 - (float)detId) - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//		}
//		// current detector element position
//		//curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//		curDetPos = rotation(make_float2(ang, -FanGeo.m_O2D), cosT, sinT);
//
//		// x-ray direction;
//		ray.d = normalize(curDetPos - ray.o);
//
//		////Only consider the case inside the FOV
//		//if (fabsf(ray.d.x * ray.o.y - ray.d.y * ray.o.x) >= (Img.m_Size.x * 0.5))
//		//{
//		//	donp[detId] = 0;
//		//	return;
//		//}
//		totWeight = 0;
//		donp[detId] = draw[angIdx * FanGeo.m_DetN + detId] - calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
//			MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg, &totWeight);
//		if (!IS_ZERO(totWeight))
//		{
//			donp[detId] /= totWeight;
//		}
//		else
//		{
//			donp[detId] = 0;
//		}
//	}
//}
////////////////////////////////////////////////////////////////////////////
//// This function is the projection difference with weighting;
//// ���������һ������ͶӰ�в�ļ�Ȩ;
//// dimg: image to be projected
//// donp: one angle projection data difference with weighting
//// draw: raw projection data
//// angIdx: angle index
//// angDir: detector from positive to negative or verse
//// FanGeo: Fan Beam geometry
//// Img: Image configuration;
//// blk: it is better to be 1024, 1
//// gid: correspondingly configurated
////////////////////////////////////////////////////////////////////////////
//void proj(float* dimg, float* donp, float* draw, cuint angIdx, const bool angDir, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker << <gid, blk >> >(dimg, donp, draw, angIdx, angDir, FanGeo, Img);
//}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//__global__ void _proj_Ker(float* dimg, float* donp, float* draw, cuint angIdx, const bool angDir, const ConeEAGeo ConeGeo, const Volume Vol)
//{
//	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint lineId = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detId < ConeGeo.m_DetN && lineId < ConeGeo.m_DetHN)
//	{
//		float3 MINO = make_float3(-Vol.m_Size.x * 0.5f + Vol.m_Bias.x, -Vol.m_Size.y * 0.5f + Vol.m_Bias.y, -Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//
//
//		//Current rotation angle;
//		float curAng = ConeGeo.m_ViwBeg + ConeGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//Current source position, assuming initial position is on the positive Y axis
//		Ray ray;
//		ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
//
//		float ang(0); //bias angle from the center of the Fan Beams
//
//		float3 curDetPos; //the current detector element position;
//		float totWeight(0);
//
//
//		// judging two indices addressing method, from large to small or small to large
//		if (angDir)
//		{
//			ang = ((float)detId - ConeGeo.m_DetCntIdx.x) * ConeGeo.m_DetStp;
//		}
//		else
//		{
//			ang = ((ConeGeo.m_DetN - 1 - (float)detId) - ConeGeo.m_DetCntIdx.x) * ConeGeo.m_DetStp;
//		}
//		// current detector element position
//		curDetPos = rotation(make_float3(
//			sinf(ang) * ConeGeo.m_S2D,
//			-cosf(ang) * ConeGeo.m_S2D + ConeGeo.m_S2O,
//			(lineId - ConeGeo.m_DetCntIdx.y) * ConeGeo.m_DetHStp), cosT, sinT);
//
//		// x-ray direction;
//		ray.d = normalize(curDetPos - ray.o);
//
//		////Only consider the case inside the FOV
//		//if (fabsf(ray.d.x * ray.o.y - ray.d.y * ray.o.x) >= (Img.m_Size.x * 0.5))
//		//{
//		//	donp[detId] = 0;
//		//	return;
//		//}
//		totWeight = 0;
//		donp[lineId * ConeGeo.m_DetN + detId] = draw[(angIdx * ConeGeo.m_DetHN + lineId) * ConeGeo.m_DetN + detId] -
//			calSiddonOneRayKer(ray.o.x, ray.o.y, ray.o.z, curDetPos.x, curDetPos.y, curDetPos.z, MINO.x, MINO.y, MINO.z,
//			Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z, Vol.m_Reso.x, Vol.m_Reso.y, Vol.m_Reso.z, dimg, &totWeight);
//		if (!IS_ZERO(totWeight))
//		{
//			donp[detId] /= totWeight;
//		}
//		else
//		{
//			donp[detId] = 0;
//		}
//	}
//}
//void proj(float* dvol, float* donp, float* draw, cuint angIdx, const bool angDir, const ConeEAGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker << <gid, blk >> >(dvol, donp, draw, angIdx, angDir, ConeGeo, Vol);
//}
//
//
//
//
//
//
//__global__ void _proj_Ker(float* dvol, float* donp, float* draw, cuint angIdx, const bool angDir, const ConeEDGeo ConeGeo, const Volume Vol)
//{
//	cuint detx = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint detz = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detx < ConeGeo.m_DetN.x && detz < ConeGeo.m_DetN.y)
//	{
//		float3 MINO = make_float3(
//			-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//			-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//			-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//		float curAng = ConeGeo.m_ViwBeg + angIdx * ConeGeo.m_ViwStp;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//		float3 curDet = rotation(
//			make_float3(
//			(detx - ConeGeo.m_DetCntIdx.x) * ConeGeo.m_DetStp.x,
//			-ConeGeo.m_O2D,
//			(detz - ConeGeo.m_DetCntIdx.y) * ConeGeo.m_DetStp.y
//			)
//			, cosT, sinT);
//		Ray ray;
//		ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
//		ray.d = normalize(curDet - ray.o);
//		float totWeight(0.0f);
//		donp[detz * ConeGeo.m_DetN.x + detx] = draw[(angIdx * ConeGeo.m_DetN.y + detz) * ConeGeo.m_DetN.x + detx]
//			- calSiddonOneRayKer(ray.o.x, ray.o.y, ray.o.z, curDet.x, curDet.y, curDet.z, MINO.x, MINO.y, MINO.z,
//			Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z, Vol.m_Reso.x, Vol.m_Reso.y, Vol.m_Reso.z,
//			dvol, &totWeight);
//		if (!IS_ZERO(totWeight))
//		{
//			donp[detz * ConeGeo.m_DetN.x + detx] /= totWeight;
//		}
//		else
//		{
//			donp[detz * ConeGeo.m_DetN.x + detx] = 0;
//		}
//	}
//}
//void proj(float* dvol, float* donp, float* draw, cuint angIdx, const bool angDir, const ConeEDGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker << <gid, blk >> >(dvol, donp, draw, angIdx, angDir, ConeGeo, Vol);
//}
//
//
//
//__global__ void _proj_Ker_V2(float* dvol, float* donp, float* draw, cuint angIdx, const bool angDir, const ConeEDGeo ConeGeo, const Volume Vol)
//{
//	cuint detx = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint detz = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detx < ConeGeo.m_DetN.x && detz < ConeGeo.m_DetN.y)
//	{
//		float3 MINO = make_float3(
//			-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//			-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//			-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//		float curAng = ConeGeo.m_ViwBeg + angIdx * ConeGeo.m_ViwStp;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//		float3 curDet = rotation(
//			make_float3(
//			(detx - ConeGeo.m_DetCntIdx.x) * ConeGeo.m_DetStp.x,
//			-ConeGeo.m_O2D,
//			(detz - ConeGeo.m_DetCntIdx.y) * ConeGeo.m_DetStp.y
//			)
//			, cosT, sinT);
//		Ray ray;
//		ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
//		ray.d = normalize(curDet - ray.o);
//		float totWeight(0.0f);
//		donp[detz * ConeGeo.m_DetN.x + detx] = draw[detz * ConeGeo.m_DetN.x + detx]
//			- calSiddonOneRayKer(ray.o.x, ray.o.y, ray.o.z, curDet.x, curDet.y, curDet.z, MINO.x, MINO.y, MINO.z,
//			Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z, Vol.m_Reso.x, Vol.m_Reso.y, Vol.m_Reso.z,
//			dvol, &totWeight);
//		if (!IS_ZERO(totWeight))
//		{
//			donp[detz * ConeGeo.m_DetN.x + detx] /= totWeight;
//		}
//		else
//		{
//			donp[detz * ConeGeo.m_DetN.x + detx] = 0;
//		}
//	}
//}
//void proj_V2(float* dvol, float* donp, float* draw, cuint angIdx, const bool angDir, const ConeEDGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker_V2 << <gid, blk >> >(dvol, donp, draw, angIdx, angDir, ConeGeo, Vol);
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
//__global__ void _proj_Ker(float* dimg, float* donp, float* draw, const bool angDir, const FanEAGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint prjIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detId < FanGeo.m_DetN && prjIdx < numPerSubSet)
//	{
//		cuint angIdx = prjIdx * subSetNum + curSubSetIdx;
//
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//
//		//Current rotation angle;
//		float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//Current source position, assuming initial position is on the positive Y axis
//		Ray2D ray;
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//		float ang(0); //bias angle from the center of the Fan Beams
//
//		float2 curDetPos; //the current detector element position;
//		float totWeight(0);
//
//
//		// judging two indices addressing method, from large to small or small to large
//		if (angDir)
//		{
//			ang = ((float)detId - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//		}
//		else
//		{
//			ang = ((FanGeo.m_DetN - 1 - (float)detId) - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//		}
//		// current detector element position
//		curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//
//		// x-ray direction;
//		ray.d = normalize(curDetPos - ray.o);
//
//		totWeight = 0;
//		donp[prjIdx * FanGeo.m_DetN + detId] = draw[angIdx * FanGeo.m_DetN + detId] - calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
//			MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg, &totWeight);
//		if (!IS_ZERO(totWeight))
//		{
//			donp[prjIdx * FanGeo.m_DetN + detId] /= totWeight;
//		}
//		else
//		{
//			donp[prjIdx * FanGeo.m_DetN + detId] = 0;
//		}
//	}
//}
////////////////////////////////////////////////////////////////////////////
//// This function is the projection difference with weighting;
//// �����������subset �е�ͶӰ�в�ļ�Ȩ;
//// dimg: image to be projected
//// donp: projection data difference with weighting inside the subset
//// draw: raw projection data
//// angDir: detector from positive to negative or verse
//// FanGeo: Fan Beam geometry
//// Img: Image configuration;
//// numPerSubSet: how many angles are included in the subset
//// subSetNum: How many subsets are divided from the original projection;
//// curSubSetIdx: current subset index
//// blk: it is better to be 1024, 1
//// gid: correspondingly configurated
////////////////////////////////////////////////////////////////////////////
//void proj(float* dimg, float* donp, float* draw, const bool angDir, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker << <gid, blk >> >(dimg, donp, draw, angDir, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx);
//}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//__global__ void _proj_Ker(float* dimg, float* donp, float* draw, const bool angDir, const FanEDGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint prjIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detId < FanGeo.m_DetN && prjIdx < numPerSubSet)
//	{
//		cuint angIdx = prjIdx * subSetNum + curSubSetIdx;
//
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//
//		//Current rotation angle;
//		float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//Current source position, assuming initial position is on the positive Y axis
//		Ray2D ray;
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//		float ang(0); //bias angle from the center of the Fan Beams
//
//		float2 curDetPos; //the current detector element position;
//		float totWeight(0);
//
//
//		// judging two indices addressing method, from large to small or small to large
//		if (angDir)
//		{
//			ang = ((float)detId - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//		}
//		else
//		{
//			ang = ((FanGeo.m_DetN - 1 - (float)detId) - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//		}
//		// current detector element position
//		//curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//		curDetPos = rotation(make_float2(ang, -FanGeo.m_O2D), cosT, sinT);
//		// x-ray direction;
//		ray.d = normalize(curDetPos - ray.o);
//
//		//Only consider the case inside the FOV
//		//if (fabsf(ray.d.x * ray.o.y - ray.d.y * ray.o.x) >= (Img.m_Size.x * 0.5))
//		//{
//		//	donp[detId] = 0;
//		//	return;
//		//}
//		totWeight = 0;
//		donp[prjIdx * FanGeo.m_DetN + detId] = draw[angIdx * FanGeo.m_DetN + detId] - calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
//			MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg, &totWeight);
//		if (!IS_ZERO(totWeight))
//		{
//			donp[prjIdx * FanGeo.m_DetN + detId] /= totWeight;
//		}
//		else
//		{
//			donp[prjIdx * FanGeo.m_DetN + detId] = 0;
//		}
//	}
//}
//void proj(float* dimg, float* donp, float* draw, const bool angDir, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker << <gid, blk >> >(dimg, donp, draw, angDir, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx);
//}
//
//
//
//
////
////
////__global__ void _proj_Ker(float* dvol, float* dprj,
////		const float3* __restrict__ cossinZT,
////		const float* __restrict__ xds,
////		const float* __restrict__ yds,
////		const float* __restrict__ zds,
////		const float* __restrict__ bxds,
////		const float* __restrict__ byds,
////		const float* __restrict__ bzds,
////		float3 objCntIdx,
////		float dx, float dz,
////		int XN, int YN, int ZN, int DNU, int DNV, int PN)
////{
////
////	int detIdV = threadIdx.y + blockIdx.y * blockDim.y; // rowIdx of the detector
////	int detIdU = threadIdx.x + blockIdx.x * blockDim.x; //channel index of the detector
////	int angIdx = threadIdx.z + blockIdx.z * blockDim.z;
////	__shared__ float _xds[BLKY];
////	__shared__ float _yds[BLKY];
////	__shared__ float _zds[BLKX];
////
////	_xds[threadIdx.y] = xds[detIdU];
////	_yds[threadIdx.y] = yds[detIdU];
////	_zds[threadIdx.x] = zds[detIdV];
////
////	__syncthreads();
////	if(detIdV < DNV && detIdU < DNU && angIdx < PN)
////	{
////		float3 MINO = make_float3(
////				(-objCntIdx.x - 0.5) * dx,
////				(-objCntIdx.y - 0.5) * dx,
////				(-objCntIdx.z - 0.5) * dz);
////		float cosT = cossinZT[angIdx].x;
////		float sinT = cossinZT[angIdx].y;
////		float zShift = cossinZT[angIdx].z;
////		float3 cursour = make_float3(
////			s.x * cosT - s.y * sinT,
////			s.x * sinT + s.y * cosT,
////			s.z + zShift);
////		float3 curDet = make_float3(
////				_xds[threadIdx.y] * cosT - _yds[threadIdx.y] * sinT,
////				_xds[threadIdx.y] * sinT + _yds[threadIdx.y] * cosT,
////				_zds[threadIdx.x] + zshift);
////
////		float3 dir = normalize(curDet - cursour);
////
////		float totWeight = 0;
////		dprj[(angIdx * DNU + detIdU) * DNV + detIdV] =
////				calSiddonOneRayKer_Zfirst(cursour.x, cursour.y, cursour.z,
////						curDet.x, curDet.y curDet.z,
////						MINO.x, MINO.x, MINO.z,
////						dx,dx,dz,XN,YN,ZN, dvol, &totWeight);
////
////
////	}
////}
////
//
//
//
//
//__global__ void _proj_Ker(float* dvol, float* donp, float* draw, const bool angDir, const ConeEAGeo ConeGeo, const Volume Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint rowId = threadIdx.y + blockIdx.y * blockDim.y;
//
//	cuint prjIdx = threadIdx.z + blockIdx.z * blockDim.z;
//	if (detId < ConeGeo.m_DetN && rowId < ConeGeo.m_DetHN && prjIdx < numPerSubSet)
//	{
//		cuint angIdx = prjIdx * subSetNum + curSubSetIdx;
//
//		float3 MINO = make_float3(-Vol.m_Size.x * 0.5f + Vol.m_Bias.x, -Vol.m_Size.y * 0.5f + Vol.m_Bias.y, -Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//
//		//Current rotation angle;
//		float curAng = ConeGeo.m_ViwBeg + ConeGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//Current source position, assuming initial position is on the positive Y axis
//		Ray ray;
//		ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
//
//		float ang(0); //bias angle from the center of the Fan Beams
//
//
//		float totWeight(0);
//
//		// judging two indices addressing method, from large to small or small to large
//		if (angDir)
//		{
//			ang = ((float)detId - ConeGeo.m_DetCntIdx.x) * ConeGeo.m_DetStp;
//		}
//		else
//		{
//			ang = ((ConeGeo.m_DetN - 1 - (float)detId) - ConeGeo.m_DetCntIdx.x) * ConeGeo.m_DetStp;
//		}
//		// current detector element position
//		float3 curDetPos = rotation(make_float3(
//			sinf(ang) * ConeGeo.m_S2D,
//			-cosf(ang) * ConeGeo.m_S2D + ConeGeo.m_S2O,
//			(rowId - ConeGeo.m_DetCntIdx.y) * ConeGeo.m_DetHStp), cosT, sinT);
//
//		// x-ray direction;
//		ray.d = normalize(curDetPos - ray.o);
//
//		//Only consider the case inside the FOV
//		//if (fabsf(ray.d.x * ray.o.y - ray.d.y * ray.o.x) >= (Img.m_Size.x * 0.5))
//		//{
//		//	donp[detId] = 0;
//		//	return;
//		//}
//		totWeight = 0;
//		donp[(prjIdx * ConeGeo.m_DetHN + rowId) * ConeGeo.m_DetN + detId] = draw[(angIdx * ConeGeo.m_DetHN + rowId) * ConeGeo.m_DetN + detId] -
//			calSiddonOneRayKer(ray.o.x, ray.o.y, ray.o.z, curDetPos.x, curDetPos.y, curDetPos.z, MINO.x, MINO.y, MINO.z,
//			Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z, Vol.m_Reso.x, Vol.m_Reso.y, Vol.m_Reso.z, dvol, &totWeight);
//
//		if (!IS_ZERO(totWeight))
//		{
//			donp[(prjIdx * ConeGeo.m_DetHN + rowId) * ConeGeo.m_DetN + detId] /= totWeight;
//		}
//		else
//		{
//			donp[(prjIdx * ConeGeo.m_DetHN + rowId) * ConeGeo.m_DetN + detId] = 0;
//		}
//	}
//}
//void proj(float* dvol, float* donp, float* draw, const bool angDir, const ConeEAGeo& ConeGeo, const Volume& Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker << <gid, blk >> >(dvol, donp, draw, angDir, ConeGeo, Vol, numPerSubSet, subSetNum, curSubSetIdx);
//}
//
//
//
//
//
////__global__ void _proj_Ker(float* dvol, float* donp, float* draw, const bool angDir, const ConeEDGeo& ConeGeo, const Volume& Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)//(float* dvol, float* draw, const bool angDir, const ConeEDGeo ConeGeo, const Volume Vol)
////{
////	cuint detidx = threadIdx.x + blockIdx.x * blockDim.x;
////	cuint detidz = threadIdx.y + blockIdx.y * blockDim.y;
////	cuint prjIdx = threadIdx.z + blockIdx.z * blockDim.z;
////	if (detidx < ConeGeo.m_DetN.x && detidz < ConeGeo.m_DetN.y && prjIdx < numPerSubSet)
////	{
////		float angIdx = curSubSetIdx * numPerSubSet + prjIdx;
////		float3 MINO = make_float3(
////			-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
////			-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
////			-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
////		float cosT = cosf(ConeGeo.m_ViwBeg + angIdx * ConeGeo.m_ViwStp);
////		float sinT = sinf(ConeGeo.m_ViwBeg + angIdx * ConeGeo.m_ViwStp);
////		float3 curDet = rotation(
////			make_float3(
////			(detidx - ConeGeo.m_DetCntIdx.x ) * ConeGeo.m_DetStp.x,
////			-ConeGeo.m_O2D,
////			(detidz - ConeGeo.m_DetCntIdx.y ) * ConeGeo.m_DetStp.y), cosT, sinT);
////		float totWeight(0.0f);
////		float3 curSor = rotation(
////			make_float3(0,ConeGeo.m_S2O,0),cosT,sinT);
////		//donp[(angIdx * ConeGeo.m_DetN.y + detidz) * ConeGeo.m_DetN.x + detidx] = draw[(angIdx * ConeGeo.m_DetN.y + detidz) * ConeGeo.m_DetN.x + detidx] = calSiddonOneRayKer(
////		//	curSor.x, curSor.y, curSor.z, curDet.x, curDet.y, curDet.z, MINO.x, MINO.y, MINO.z,
////		//Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z,Vol.m_Reso.x, Vol.m_Reso.y, Vol.m_Reso.z,
////		//	dvol,&totWeight);
////	}
////}
////
////void proj(float* dvol, float* donp, float* draw, const bool angDir, const ConeEDGeo& ConeGeo, const Volume& Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
////{
////	_proj_Ker<<<gid,blk>>>(dvol, draw, angDir, ConeGeo, Vol, numPerSubSet, subSetNum, curSubSetIdx);
////}
//
//
//
//
//__global__ void _proj_Ker(float* dimg, float* donp, float* draw, const FanEAGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint prjIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detId < FanGeo.m_DetN && prjIdx < numPerSubSet)
//	{
//		cuint angIdx = prjIdx * subSetNum + curSubSetIdx;
//
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//
//		//Current rotation angle;
//		float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//Current source position, assuming initial position is on the positive Y axis
//		Ray2D ray;
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//		float ang(0); //bias angle from the center of the Fan Beams
//
//		float2 curDetPos; //the current detector element position;
//		float totWeight(0);
//
//
//		// judging two indices addressing method, from large to small or small to large
//		ang = ((float)detId - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//		// current detector element position
//		curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//
//		// x-ray direction;
//		ray.d = normalize(curDetPos - ray.o);
//
//		totWeight = 0;
//		donp[prjIdx * FanGeo.m_DetN + detId] = draw[angIdx * FanGeo.m_DetN + detId] - calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
//			MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg, &totWeight);
//		if (!IS_ZERO(totWeight))
//		{
//			donp[prjIdx * FanGeo.m_DetN + detId] /= totWeight;
//		}
//		else
//		{
//			donp[prjIdx * FanGeo.m_DetN + detId] = 0;
//		}
//	}
//}
////////////////////////////////////////////////////////////////////////////
//// This function is the projection difference with weighting;
//// �����������subset �е�ͶӰ�в�ļ�Ȩ; ����������,ֻ��Ϊ������;
//// dimg: image to be projected
//// donp: projection data difference with weighting inside the subset
//// draw: raw projection data
//// FanGeo: Fan Beam geometry
//// Img: Image configuration;
//// numPerSubSet: how many angles are included in the subset
//// subSetNum: How many subsets are divided from the original projection;
//// curSubSetIdx: current subset index
//// blk: it is better to be 1024, 1
//// gid: correspondingly configurated
////////////////////////////////////////////////////////////////////////////
//void proj(float* dimg, float* donp, float* draw, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker << <gid, blk >> >(dimg, donp, draw, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx);
//}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//__global__ void _proj_Ker(float* dimg, float* donp, float* draw, const FanEDGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint prjIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detId < FanGeo.m_DetN && prjIdx < numPerSubSet)
//	{
//
//		cuint angIdx = prjIdx * subSetNum + curSubSetIdx;
//
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//
//		//Current rotation angle;
//		float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//Current source position, assuming initial position is on the positive Y axis
//		Ray2D ray;
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//		float ang(0); //bias angle from the center of the Fan Beams
//
//		float2 curDetPos; //the current detector element position;
//		float totWeight(0);
//
//
//		// judging two indices addressing method, from large to small or small to large
//
//		ang = ((float)detId - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//		// current detector element position
//		//curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//		curDetPos = rotation(make_float2(ang, -FanGeo.m_O2D), cosT, sinT);
//		// x-ray direction;
//		ray.d = normalize(curDetPos - ray.o);
//
//		totWeight = 0;
//		donp[prjIdx * FanGeo.m_DetN + detId] = draw[angIdx * FanGeo.m_DetN + detId] - calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
//			MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg, &totWeight);
//		if (!IS_ZERO(totWeight))
//		{
//			donp[prjIdx * FanGeo.m_DetN + detId] /= totWeight;
//		}
//		else
//		{
//			donp[prjIdx * FanGeo.m_DetN + detId] = 0;
//		}
//	}
//}
//void proj(float* dimg, float* donp, float* draw, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker << <gid, blk >> >(dimg, donp, draw, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx);
//}
//
//
//
//
//
//__global__ void _proj_Ker(float* dvol, float* donp, float* draw, const ConeEAGeo ConeGeo, const Volume Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint rowId = threadIdx.y + blockIdx.y * blockDim.y;
//	cuint prjIdx = threadIdx.z + blockIdx.z * blockDim.z;
//	if (detId < ConeGeo.m_DetN && rowId < ConeGeo.m_DetHN && prjIdx < numPerSubSet)
//	{
//		cuint angIdx = prjIdx * subSetNum + curSubSetIdx;
//
//		float3 MINO = make_float3(
//			-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//			-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//			-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//
//		//Current rotation angle;
//		float curAng = ConeGeo.m_ViwBeg + ConeGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//Current source position, assuming initial position is on the positive Y axis
//		Ray ray;
//		ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
//
//		float ang(0); //bias angle from the center of the Fan Beams
//
//		float3 curDetPos; //the current detector element position;
//		float totWeight(0);
//
//
//		// judging two indices addressing method, from large to small or small to large
//		ang = ((float)detId - ConeGeo.m_DetCntIdx.x) * ConeGeo.m_DetStp;
//
//		// current detector element position
//		curDetPos = rotation(make_float3(
//			sinf(ang) * ConeGeo.m_S2D,
//			-cosf(ang) * ConeGeo.m_S2D + ConeGeo.m_S2O,
//			((float)rowId - ConeGeo.m_DetCntIdx.y) * ConeGeo.m_DetHStp
//			), cosT, sinT);
//
//		// x-ray direction;
//		ray.d = normalize(curDetPos - ray.o);
//
//		totWeight = 0;
//		donp[(prjIdx * ConeGeo.m_DetHN + rowId) * ConeGeo.m_DetN + detId] = draw[(angIdx * ConeGeo.m_DetHN + rowId) * ConeGeo.m_DetN + detId] -
//			calSiddonOneRayKer(ray.o.x, ray.o.y, ray.o.z, curDetPos.x, curDetPos.y, curDetPos.z, MINO.x, MINO.y, MINO.z,
//			Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z, Vol.m_Reso.x, Vol.m_Reso.y, Vol.m_Reso.z, dvol, &totWeight);
//		//calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg, &totWeight);
//		if (!IS_ZERO(totWeight))
//		{
//			donp[(prjIdx * ConeGeo.m_DetHN + rowId) * ConeGeo.m_DetN + detId] /= totWeight;
//		}
//		else
//		{
//			donp[(prjIdx * ConeGeo.m_DetHN + rowId) * ConeGeo.m_DetN + detId] = 0;
//		}
//	}
//}
//void proj(float* dvol, float* donp, float* draw, const ConeEAGeo& ConeGeo, const Volume& Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker << <gid, blk >> >(dvol, donp, draw, ConeGeo, Vol, numPerSubSet, subSetNum, curSubSetIdx);
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
//__global__ void _proj_Ker(float* dimg, float* draw, const bool angDir, const FanEAGeo FanGeo, const Image Img)
//{
//	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint prjIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detId < FanGeo.m_DetN && prjIdx < FanGeo.m_ViwN)
//	{
//		cuint angIdx = prjIdx;
//
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//
//		//Current rotation angle;
//		float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//Current source position, assuming initial position is on the positive Y axis
//		Ray2D ray;
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//		float ang(0); //bias angle from the center of the Fan Beams
//
//		float2 curDetPos; //the current detector element position;
//		float totWeight(0);
//
//
//		// judging two indices addressing method, from large to small or small to large
//		if (angDir)
//		{
//			ang = ((float)detId - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//		}
//		else
//		{
//			ang = ((FanGeo.m_DetN - 1 - (float)detId) - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//		}
//		// current detector element position
//		curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//
//		// x-ray direction;
//		ray.d = normalize(curDetPos - ray.o);
//
//		totWeight = 0;
//		draw[prjIdx * FanGeo.m_DetN + detId] = calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
//			MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg, &totWeight);
//	}
//}
////////////////////////////////////////////////////////////////////////////
//// This function is the projection difference with weighting;
//// ��������������нǶȵ�ͶӰ;
//// dimg: image to be projected
//// draw: projection data
//// angDir: detector from positive to negative or verse
//// FanGeo: Fan Beam geometry
//// Img: Image configuration;
//// blk: it is better to be 1024, 1
//// gid: correspondingly configurated
////////////////////////////////////////////////////////////////////////////
//void proj(float* dimg, float* draw, const bool angDir, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker << <gid, blk >> >(dimg, draw, angDir, FanGeo, Img);
//}
//
//
//
//
//void proj_CPU(float*dimg, float* dprj, const FanEAGeo& FanGeo, const Image& Img)
//{
//	unsigned int detId(0), prjIdx(0);
//	for (prjIdx = 0; prjIdx < FanGeo.m_ViwN; prjIdx++)
//	{
//		for (detId = 0; detId < FanGeo.m_DetN; detId++)
//		{
//			cuint angIdx = prjIdx;
//
//			float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//
//			//Current rotation angle;
//			float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//			float cosT = cosf(curAng);
//			float sinT = sinf(curAng);
//
//			//Current source position, assuming initial position is on the positive Y axis
//			Ray2D ray;
//			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//			float ang(0); //bias angle from the center of the Fan Beams
//
//			float2 curDetPos; //the current detector element position;
//			float totWeight(0);
//
//
//			// judging two indices addressing method, from large to small or small to large
//			ang = ((float)detId - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//			// current detector element position
//			curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//
//			// x-ray direction;
//			ray.d = normalize(curDetPos - ray.o);
//
//			totWeight = 0;
//			dprj[prjIdx * FanGeo.m_DetN + detId] = calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
//				MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg, &totWeight);
//		}
//
//	}
//}
//
//
//
//
//void proj_CPU_OPENMP(float*dimg, float* dprj, const FanEAGeo& FanGeo, const Image& Img)
//{
//	int prjIdx(0);
//#pragma omp parallel for
//	for (prjIdx = 0; prjIdx < FanGeo.m_ViwN; prjIdx++)
//	{
//		unsigned int detId(0);
//		for (detId = 0; detId < FanGeo.m_DetN; detId++)
//		{
//			cuint angIdx = prjIdx;
//
//			float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//
//			//Current rotation angle;
//			float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//			float cosT = cosf(curAng);
//			float sinT = sinf(curAng);
//
//			//Current source position, assuming initial position is on the positive Y axis
//			Ray2D ray;
//			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//			float ang(0); //bias angle from the center of the Fan Beams
//
//			float2 curDetPos; //the current detector element position;
//			float totWeight(0);
//
//
//			// judging two indices addressing method, from large to small or small to large
//			ang = ((float)detId - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//			// current detector element position
//			curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//
//			// x-ray direction;
//			ray.d = normalize(curDetPos - ray.o);
//
//			totWeight = 0;
//			dprj[prjIdx * FanGeo.m_DetN + detId] = calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
//				MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg, &totWeight);
//		}
//
//	}
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
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//__global__ void _proj_Ker(float* dimg, float* draw, const bool angDir, const FanEDGeo FanGeo, const Image Img)
//{
//	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint prjIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detId < FanGeo.m_DetN && prjIdx < FanGeo.m_ViwN)
//	{
//		cuint angIdx = prjIdx;
//
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//
//		//Current rotation angle;
//		float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//Current source position, assuming initial position is on the positive Y axis
//		Ray2D ray;
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//		float ang(0); //bias angle from the center of the Fan Beams
//
//		float2 curDetPos; //the current detector element position;
//		float totWeight(0);
//
//
//		// judging two indices addressing method, from large to small or small to large
//		if (angDir)
//		{
//			ang = ((float)detId - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//		}
//		else
//		{
//			ang = ((FanGeo.m_DetN - 1 - (float)detId) - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//		}
//		// current detector element position
//		//curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//		curDetPos = rotation(make_float2(ang, -FanGeo.m_O2D), cosT, sinT);
//		// x-ray direction;
//		ray.d = normalize(curDetPos - ray.o);
//
//		totWeight = 0;
//		draw[prjIdx * FanGeo.m_DetN + detId] = calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
//			MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg, &totWeight);
//	}
//}
//void proj(float* dimg, float* draw, const bool angDir, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker << <gid, blk >> >(dimg, draw, angDir, FanGeo, Img);
//}
//
//
//__global__ void _proj_Ker(float* dimg, float* ddif, float* draw, const FanEDGeo FanGeo, const Image Img)
//{
//	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint prjIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detId < FanGeo.m_DetN && prjIdx < FanGeo.m_ViwN)
//	{
//		cuint angIdx = prjIdx;
//
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//
//		//Current rotation angle;
//		float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//Current source position, assuming initial position is on the positive Y axis
//		Ray2D ray;
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//		float ang(0); //bias angle from the center of the Fan Beams
//
//		float2 curDetPos; //the current detector element position;
//		float totWeight(0);
//
//
//		// judging two indices addressing method, from large to small or small to large
//		ang = ((float)detId - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//
//		// current detector element position
//		//curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//		curDetPos = rotation(make_float2(ang, -FanGeo.m_O2D), cosT, sinT);
//		// x-ray direction;
//		ray.d = normalize(curDetPos - ray.o);
//
//		totWeight = 0;
//		float val = draw[prjIdx * FanGeo.m_DetN + detId] - calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
//			MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg, &totWeight);
//		if (!IS_ZERO(totWeight))
//		{
//			ddif[prjIdx * FanGeo.m_DetN + detId] = val / totWeight;
//		}
//		else
//		{
//			ddif[prjIdx * FanGeo.m_DetN + detId] = 0;
//		}
//	}
//}
//void proj(float* dimg, float* ddif, float* draw, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker << <gid, blk >> >(dimg, ddif, draw, FanGeo, Img);
//}
////
////template<typename T>
////void proj_CPU_temp(T* dimg, T* ddif, T* draw, const FanEDGeo& FanGeo, const Image& Img)
////{
////	unsigned int detId = 0;
////	unsigned int prjIdx = 0;
////	for (prjIdx = 0; prjIdx < FanGeo.m_ViwN; prjIdx++)
////	{
////		for (detId = 0; detId < FanGeo.m_DetN; detId++)
////		{
////			cuint angIdx = prjIdx;
////
////			float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
////
////			//Current rotation angle;
////			float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
////			float cosT = cosf(curAng);
////			float sinT = sinf(curAng);
////
////			//Current source position, assuming initial position is on the positive Y axis
////			Ray2D ray;
////			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
////
////			float ang(0); //bias angle from the center of the Fan Beams
////
////			float2 curDetPos; //the current detector element position;
////			float totWeight(0);
////
////
////			// judging two indices addressing method, from large to small or small to large
////			ang = ((float) detId - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
////
////			// current detector element position
////			//curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
////			curDetPos = rotation(make_float2(ang, -FanGeo.m_O2D), cosT, sinT);
////			// x-ray direction;
////			ray.d = normalize(curDetPos - ray.o);
////
////			totWeight = 0;
////			float val = draw[prjIdx * FanGeo.m_DetN + detId] - calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
////				MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg, &totWeight);
////			if (!IS_ZERO(totWeight))
////			{
////				ddif[prjIdx * FanGeo.m_DetN + detId] = val / totWeight;
////			}
////			else
////			{
////				ddif[prjIdx * FanGeo.m_DetN + detId] = 0;
////			}
////		}
////	}
////}
//
//void proj_CPU(float* dimg, float* ddif, float* draw, const FanEDGeo& FanGeo, const Image& Img)
//{
//	unsigned int detId = 0;
//	unsigned int prjIdx = 0;
//	for (prjIdx = 0; prjIdx < FanGeo.m_ViwN; prjIdx++)
//	{
//		for (detId = 0; detId < FanGeo.m_DetN; detId++)
//		{
//			cuint angIdx = prjIdx;
//
//			float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//
//			//Current rotation angle;
//			float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//			float cosT = cosf(curAng);
//			float sinT = sinf(curAng);
//
//			//Current source position, assuming initial position is on the positive Y axis
//			Ray2D ray;
//			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//			float ang(0); //bias angle from the center of the Fan Beams
//
//			float2 curDetPos; //the current detector element position;
//			float totWeight(0);
//
//
//			// judging two indices addressing method, from large to small or small to large
//			ang = ((float)detId - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//
//			// current detector element position
//			//curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//			curDetPos = rotation(make_float2(ang, -FanGeo.m_O2D), cosT, sinT);
//			// x-ray direction;
//			ray.d = normalize(curDetPos - ray.o);
//
//			totWeight = 0;
//			float val = draw[prjIdx * FanGeo.m_DetN + detId] - calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
//				MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg, &totWeight);
//			if (!IS_ZERO(totWeight))
//			{
//				ddif[prjIdx * FanGeo.m_DetN + detId] = val / totWeight;
//			}
//			else
//			{
//				ddif[prjIdx * FanGeo.m_DetN + detId] = 0;
//			}
//		}
//	}
//}
//
//void proj_CPU_OPENMP(float* dimg, float* ddif, float* draw, const FanEDGeo& FanGeo, const Image& Img)
//{
//
//	int prjIdx = 0;
//	CUDA_CHECK_RETURN(cudaMalloc((void**)&draw[0], sizeof(float) * FanGeo.m_DetN *FanGeo.m_ViwN));
//#pragma omp parallel for
//	for (prjIdx = 0; prjIdx < FanGeo.m_ViwN; prjIdx++)
//	{
//		unsigned int detId = 0;
//		for (detId = 0; detId < FanGeo.m_DetN; detId++)
//		{
//			cuint angIdx = prjIdx;
//
//			float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//
//			//Current rotation angle;
//			float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//			float cosT = cosf(curAng);
//			float sinT = sinf(curAng);
//
//			//Current source position, assuming initial position is on the positive Y axis
//			Ray2D ray;
//			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//			float ang(0); //bias angle from the center of the Fan Beams
//
//			float2 curDetPos; //the current detector element position;
//			float totWeight(0);
//
//
//			// judging two indices addressing method, from large to small or small to large
//			ang = ((float)detId - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//
//			// current detector element position
//			//curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//			curDetPos = rotation(make_float2(ang, -FanGeo.m_O2D), cosT, sinT);
//			// x-ray direction;
//			ray.d = normalize(curDetPos - ray.o);
//
//			totWeight = 0;
//			float val = draw[prjIdx * FanGeo.m_DetN + detId] - calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
//				MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg, &totWeight);
//			if (!IS_ZERO(totWeight))
//			{
//				ddif[prjIdx * FanGeo.m_DetN + detId] = val / totWeight;
//			}
//			else
//			{
//				ddif[prjIdx * FanGeo.m_DetN + detId] = 0;
//			}
//		}
//	}
//}
//
////void proj_CPU(double* dimg, double* ddiff, double* ddraw, const FanEDGeo& FanGeo, const Image& Img)
////{
////	proj_CPU_temp<double>(dimg, ddiff, ddraw, FanGeo, Img);
////}
//
//
//
////
////__global__ void _proj_Ker(float* dimg,float* dcor, float* draw, const FanEDGeo FanGeo, const Image Img)
////{
////	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
////	cuint prjIdx = threadIdx.y + blockIdx.y * blockDim.y;
////	if (detId < FanGeo.m_DetN && prjIdx < FanGeo.m_ViwN)
////	{
////		cuint angIdx = prjIdx;
////
////		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
////
////		//Current rotation angle;
////		float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
////		float cosT = cosf(curAng);
////		float sinT = sinf(curAng);
////
////		//Current source position, assuming initial position is on the positive Y axis
////		Ray2D ray;
////		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
////
////		float ang(0); //bias angle from the center of the Fan Beams
////
////		float2 curDetPos; //the current detector element position;
////		float totWeight(0);
////
////
////		// judging two indices addressing method, from large to small or small to large
////		ang = ((float) detId - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
////		// current detector element position
////		//curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
////		curDetPos = rotation(make_float2(ang,-FanGeo.m_O2D),cosT,sinT);
////		// x-ray direction;
////		ray.d = normalize(curDetPos - ray.o);
////
////		totWeight = 0;
////		float val = draw[prjIdx * FanGeo.m_DetN + detId] - calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
////			MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg, &totWeight);
////		if (IS_ZERO(totWeight))
////		{
////			dcor[prjIdx * FanGeo.m_DetN + detId] = 0;
////		}
////		else
////		{
////			dcor[prjIdx * FanGeo.m_DetN + detId] = val / totWeight;
////		}
////	}
////}
////void proj(float* dimg, float* dcor, float* draw, const FanEDGeo& FanGeo, const Image& Img,const dim3& blk, const dim3& gid)
////{	_proj_Ker<<<gid,blk>>>(dimg,dcor, draw,FanGeo,Img);}
//
//
//
//
//
//
//__global__ void _proj_Ker(float* dvol, float* draw, const bool angDir, const ConeEAGeo ConeGeo, const Volume Vol)
//{
//	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint rowId = threadIdx.y + blockIdx.y * blockDim.y;
//
//	cuint prjIdx = threadIdx.z + blockIdx.z * blockDim.z;
//	if (detId < ConeGeo.m_DetN && rowId && ConeGeo.m_DetHN && prjIdx < ConeGeo.m_ViwN)
//	{
//		cuint angIdx = prjIdx;
//
//		float3 MINO = make_float3(
//			-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//			-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//			-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//
//		//Current rotation angle;
//		float curAng = ConeGeo.m_ViwBeg + ConeGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//Current source position, assuming initial position is on the positive Y axis
//		Ray ray;
//		ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
//
//		float ang(0); //bias angle from the center of the Fan Beams
//
//		float3 curDetPos; //the current detector element position;
//		float totWeight(0);
//
//
//		// judging two indices addressing method, from large to small or small to large
//		if (angDir)
//		{
//			ang = ((float)detId - ConeGeo.m_DetCntIdx.x) * ConeGeo.m_DetStp;
//		}
//		else
//		{
//			ang = ((ConeGeo.m_DetN - 1 - (float)detId) - ConeGeo.m_DetCntIdx.x) * ConeGeo.m_DetStp;
//		}
//		// current detector element position
//		curDetPos = rotation(make_float3(
//			sinf(ang) * ConeGeo.m_S2D,
//			-cosf(ang) * ConeGeo.m_S2D + ConeGeo.m_S2O,
//			((float)rowId - ConeGeo.m_DetCntIdx.y) * ConeGeo.m_DetHStp),
//			cosT, sinT);
//
//		// x-ray direction;
//		ray.d = normalize(curDetPos - ray.o);
//
//		totWeight = 0;
//		draw[(prjIdx * ConeGeo.m_DetHN + rowId) * ConeGeo.m_DetN + detId] =
//			calSiddonOneRayKer(ray.o.x, ray.o.y, ray.o.z, curDetPos.x, curDetPos.y, curDetPos.z, MINO.x, MINO.y, MINO.z,
//			Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z, Vol.m_Reso.x, Vol.m_Reso.y, Vol.m_Reso.z, dvol, &totWeight);
//	}
//}
//void proj(float* dvol, float* draw, const bool angDir, const ConeEAGeo& FanGeo, const Volume& Img, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker << <gid, blk >> >(dvol, draw, angDir, FanGeo, Img);
//}
//
//
//
//
//__global__ void _proj_Ker(float* dvol, float* draw, const bool angDir, const ConeEDGeo ConeGeo, const Volume Vol)
//{
//	cuint detidx = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint detidz = threadIdx.y + blockIdx.y * blockDim.y;
//	cuint angIdx = threadIdx.z + blockIdx.z * blockDim.z;
//	if (detidx < ConeGeo.m_DetN.x && detidz < ConeGeo.m_DetN.y && angIdx < ConeGeo.m_ViwN)
//	{
//		float3 MINO = make_float3(
//			-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//			-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//			-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//		float cosT = cosf(ConeGeo.m_ViwBeg + angIdx * ConeGeo.m_ViwStp);
//		float sinT = sinf(ConeGeo.m_ViwBeg + angIdx * ConeGeo.m_ViwStp);
//		float3 curDet = rotation(
//			make_float3(
//			(detidx - ConeGeo.m_DetCntIdx.x) * ConeGeo.m_DetStp.x,
//			-ConeGeo.m_O2D,
//			(detidz - ConeGeo.m_DetCntIdx.y) * ConeGeo.m_DetStp.y), cosT, sinT);
//		float totWeight(0.0f);
//		float3 curSor = rotation(
//			make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
//		draw[(angIdx * ConeGeo.m_DetN.y + detidz) * ConeGeo.m_DetN.x + detidx] = calSiddonOneRayKer(
//			curSor.x, curSor.y, curSor.z, curDet.x, curDet.y, curDet.z, MINO.x, MINO.y, MINO.z,
//			Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z, Vol.m_Reso.x, Vol.m_Reso.y, Vol.m_Reso.z,
//			dvol, &totWeight);
//	}
//}
//void proj(float* dvol, float* draw, const bool angDir, const ConeEDGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker << <gid, blk >> >(dvol, draw, angDir, ConeGeo, Vol);
//}
//
//
//
//
//
//
//
//__global__ void _proj_Ker(float* dimg, float* draw, const FanEAGeo FanGeo, const Image Img)
//{
//	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint prjIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detId < FanGeo.m_DetN && prjIdx < FanGeo.m_ViwN)
//	{
//		cuint angIdx = prjIdx;
//
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//
//		//Current rotation angle;
//		float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//Current source position, assuming initial position is on the positive Y axis
//		Ray2D ray;
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//		float ang(0); //bias angle from the center of the Fan Beams
//
//		float2 curDetPos; //the current detector element position;
//		float totWeight(0);
//
//
//		// judging two indices addressing method, from large to small or small to large
//
//		ang = ((float)detId - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//		// current detector element position
//		curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//
//		// x-ray direction;
//		ray.d = normalize(curDetPos - ray.o);
//
//		totWeight = 0;
//		draw[prjIdx * FanGeo.m_DetN + detId] = calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
//			MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg, &totWeight);
//	}
//}
////////////////////////////////////////////////////////////////////////////
//// This function is the projection difference with weighting;
//// ��������������нǶȵ�ͶӰ; ����������ֻ��������;
//// dimg: image to be projected
//// draw: projection data
//// FanGeo: Fan Beam geometry
//// Img: Image configuration;
//// blk: it is better to be 1024, 1
//// gid: correspondingly configurated
////////////////////////////////////////////////////////////////////////////
//void proj(float* dimg, float* draw, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker << <gid, blk >> >(dimg, draw, FanGeo, Img);
//}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//__global__ void _proj_Ker(float* dimg, float* draw, const FanEDGeo FanGeo, const Image Img)
//{
//	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint prjIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detId < FanGeo.m_DetN && prjIdx < FanGeo.m_ViwN)
//	{
//		cuint angIdx = prjIdx;
//
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//
//		//Current rotation angle;
//		float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//Current source position, assuming initial position is on the positive Y axis
//		Ray2D ray;
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//		float ang(0); //bias angle from the center of the Fan Beams
//
//		float2 curDetPos; //the current detector element position;
//		float totWeight(0);
//
//
//		// judging two indices addressing method, from large to small or small to large
//
//		ang = ((float)detId - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//		// current detector element position
//		//curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//		curDetPos = rotation(make_float2(ang, -FanGeo.m_O2D), cosT, sinT);
//		// x-ray direction;
//		ray.d = normalize(curDetPos - ray.o);
//
//		totWeight = 0;
//		draw[prjIdx * FanGeo.m_DetN + detId] = calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
//			MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg, &totWeight);
//	}
//}
//void proj(float* dimg, float* draw, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker << <gid, blk >> >(dimg, draw, FanGeo, Img);
//}
//
//
//
//
//
//
//
//__global__ void _proj_Ker(float* dvol, float* draw, const ConeEAGeo ConeGeo, const Volume Vol)
//{
//	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint rowId = threadIdx.y + blockIdx.y * blockDim.y;
//
//	cuint angIdx = threadIdx.z + blockIdx.z * blockDim.z;
//	if (detId < ConeGeo.m_DetN && rowId < ConeGeo.m_DetHN && angIdx < ConeGeo.m_ViwN)
//	{
//
//		const float3 MINO = make_float3(
//			-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//			-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//			-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//
//		float curAng = ConeGeo.m_ViwBeg + ConeGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//
//		float ang = ((float)detId - ConeGeo.m_DetCntIdx.x) * ConeGeo.m_DetStp;
//
//		// current detector element position
//		float3 curDetPos = rotation(make_float3(
//			sinf(ang) * ConeGeo.m_S2D,
//			-cosf(ang) * ConeGeo.m_S2D + ConeGeo.m_S2O,
//			((float)rowId - ConeGeo.m_DetCntIdx.y) * ConeGeo.m_DetHStp),
//			cosT, sinT);
//
//		// x-ray direction;
//		Ray ray;
//		ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
//		ray.d = normalize(curDetPos - ray.o);
//
//		float totWeight = 0;
//		draw[(angIdx * ConeGeo.m_DetHN + rowId) * ConeGeo.m_DetN + detId] =
//			calSiddonOneRayKer(ray.o.x, ray.o.y, ray.o.z, curDetPos.x, curDetPos.y, curDetPos.z, MINO.x, MINO.y, MINO.z,
//			Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z, Vol.m_Reso.x, Vol.m_Reso.y, Vol.m_Reso.z, dvol, &totWeight);
//
//
//
//		//float3 MINO = make_float3(
//		//	-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//		//	-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//		//	-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//
//		////Current rotation angle;
//		//float curAng = ConeGeo.m_ViwBeg + ConeGeo.m_ViwStp * angIdx;
//		//float cosT = cosf(curAng);
//		//float sinT = sinf(curAng);
//
//		//float ang = ((float) detId - ConeGeo.m_DetCntIdx.x) * ConeGeo.m_DetStp; //bias angle from the center of the Fan Beams
//
//		//// current detector element position
//		//float3 curDetPos = rotation(make_float3(
//		//	sinf(ang) * ConeGeo.m_S2D,
//		//	-cosf(ang) * ConeGeo.m_S2D + ConeGeo.m_S2O,
//		//	((float)rowId - ConeGeo.m_DetCntIdx.y) * ConeGeo.m_DetHStp),
//		//	cosT, sinT);
//
//		//// x-ray direction;
//		////Current source position, assuming initial position is on the positive Y axis
//		//Ray ray;
//		//ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
//		//ray.d = normalize(curDetPos - ray.o);
//
//		//float totWeight(0);
//		//draw[(angIdx * ConeGeo.m_DetHN + rowId) * ConeGeo.m_DetN + detId] = 
//		//	calSiddonOneRayKer(ray.o.x, ray.o.y, ray.o.z, curDetPos.x, curDetPos.y, curDetPos.z, MINO.x, MINO.y, MINO.z,
//		//	Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z, Vol.m_Reso.x, Vol.m_Reso.y, Vol.m_Reso.z, dvol, &totWeight);
//	}
//}
//void proj(float* dvol, float* draw, const ConeEAGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker << <gid, blk >> >(dvol, draw, ConeGeo, Vol);
//}
//
//
//
//
//
//
//
//__global__ void _proj_Ker(float* dvol, float* dprj, const ConeEDGeo ConeGeo, const Volume Vol)
//{
//	cuint detIdx = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint detIdz = threadIdx.y + blockIdx.y * blockDim.y;
//	cuint angIdx = threadIdx.z + blockIdx.z * blockDim.z;
//	if (detIdx < ConeGeo.m_DetN.x && detIdz < ConeGeo.m_DetN.y && angIdx < ConeGeo.m_ViwN)
//	{
//		float3 MINO = make_float3(
//			-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//			-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//			-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//
//		//current rotation
//		float cosT = cosf(ConeGeo.m_ViwBeg + (float)angIdx * ConeGeo.m_ViwStp);
//		float sinT = sinf(ConeGeo.m_ViwBeg + (float)angIdx * ConeGeo.m_ViwStp);
//
//		Ray ray;
//		ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
//		float3 curDet = rotation(make_float3(
//			((float)detIdx - ConeGeo.m_DetCntIdx.x) * ConeGeo.m_DetStp.x,
//			-ConeGeo.m_O2D,
//			((float)detIdz - ConeGeo.m_DetCntIdx.y) * ConeGeo.m_DetStp.y), cosT, sinT);
//		ray.d = normalize(curDet - ray.o);
//		float totWeight(0.0f);
//		dprj[(angIdx * ConeGeo.m_DetN.y + detIdz) * ConeGeo.m_DetN.x + detIdx] =
//			calSiddonOneRayKer(
//			ray.o.x, ray.o.y, ray.o.z,
//			curDet.x, curDet.y, curDet.z,
//			MINO.x, MINO.y, MINO.z,
//			Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z,
//			Vol.m_Reso.x, Vol.m_Reso.y, Vol.m_Reso.z, dvol, &totWeight);
//	}
//}
//void proj(float* dvol, float* dprj, const ConeEDGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker << <gid, blk >> >(dvol, dprj, ConeGeo, Vol);
//}
//
//
//
//
//
//__global__ void _proj_Ker(float* dimg, float* dcor, float* draw, const FanEAGeo FanGeo, const Image Img)
//{
//	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint prjIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detId < FanGeo.m_DetN && prjIdx < FanGeo.m_ViwN)
//	{
//		cuint angIdx = prjIdx;
//
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//
//		//Current rotation angle;
//		float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//Current source position, assuming initial position is on the positive Y axis
//		Ray2D ray;
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//		float ang(0); //bias angle from the center of the Fan Beams
//
//		float2 curDetPos; //the current detector element position;
//		float totWeight(0);
//
//
//		// judging two indices addressing method, from large to small or small to large
//
//		ang = ((float)detId - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//		// current detector element position
//		curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//
//		// x-ray direction;
//		ray.d = normalize(curDetPos - ray.o);
//
//		totWeight = 0;
//		dcor[prjIdx * FanGeo.m_DetN + detId] = draw[prjIdx * FanGeo.m_DetN + detId] - calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
//			MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg, &totWeight);
//		if (IS_ZERO(totWeight))
//		{
//			dcor[prjIdx * FanGeo.m_DetN + detId] = 0;
//		}
//		else
//		{
//			dcor[prjIdx * FanGeo.m_DetN + detId] /= totWeight;
//		}
//	}
//}
////////////////////////////////////////////////////////////////////////////
//// This function is the projection difference with weighting;
//// ��������������нǶȵ�ͶӰ; ����������ֻ��������;
//// dimg: image to be projected
//// draw: projection data
//// FanGeo: Fan Beam geometry
//// Img: Image configuration;
//// blk: it is better to be 1024, 1
//// gid: correspondingly configurated
////////////////////////////////////////////////////////////////////////////
//void proj(float* dimg, float* dcor, float* draw, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker << <gid, blk >> >(dimg, dcor, draw, FanGeo, Img);
//}
//
//
//
//
//
//
//__global__ void proj_Ker(float* dimg, float* dcor, float* draw, const FanEAGeo FanGeo, const Image Img, cuint sliceNum)
//{
//	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint angIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detId < FanGeo.m_DetN && angIdx < FanGeo.m_ViwN)
//	{
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//
//		//Current rotation angle;
//		float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//Current source position, assuming initial position is on the positive Y axis
//		Ray2D ray;
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//		//bias angle from the center of the Fan Beams
//		float ang = ((float)detId - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//		// current detector element position
//		float2 curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//
//		// x-ray direction;
//		ray.d = normalize(curDetPos - ray.o);
//
//		int sliceId = 0;
//		float totWeight(0);
//		for (sliceId = 0; sliceId != sliceNum; ++sliceId)
//		{
//			totWeight = 0;
//			curAng = draw[(sliceId * FanGeo.m_ViwN + angIdx) * FanGeo.m_DetN + detId] - calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
//				MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg + (sliceId * Img.m_Reso.x * Img.m_Reso.y) // Calculate the different slices in all loops
//				, &totWeight);
//			if (!IS_ZERO(totWeight))
//			{
//				dcor[(sliceId * FanGeo.m_ViwN + angIdx) * FanGeo.m_DetN + detId] = curAng / totWeight;
//			}
//			else
//			{
//				dcor[(sliceId * FanGeo.m_ViwN + angIdx) * FanGeo.m_DetN + detId] = 0;
//			}
//		}
//
//	}
//}
////��ϵ�к����ʾͬ�������configuration���ڲ�ͬ����ͬʱ��ͶӰ�ĺ���Ӧ����patientRec �Ǹ�64����ع�
//void proj(float* dimg, float* dcor, float* draw, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid, cuint sliceNum)
//{
//	proj_Ker << <gid, blk >> >(dimg, dcor, draw, FanGeo, Img, sliceNum);
//}
//
//
//
//
//
//
//
//__global__ void proj_Ker(float* dimg, float* donp, float* draw, const FanEAGeo FanGeo, const Image Img, cuint sliceNum, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint prjIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detId < FanGeo.m_DetN && prjIdx < numPerSubSet)
//	{
//		cuint angIdx = prjIdx * subSetNum + curSubSetIdx;
//
//		float2 MINO = make_float2(-Img.m_Size.x / 2.0f + Img.m_Bias.x, -Img.m_Size.y / 2.0f + Img.m_Bias.y);
//
//		//Current rotation angle;
//		float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//Current source position, assuming initial position is on the positive Y axis
//		Ray2D ray;
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//		float ang(0); //bias angle from the center of the Fan Beams
//
//		float2 curDetPos; //the current detector element position;
//		float totWeight(0);
//
//
//		// judging two indices addressing method, from large to small or small to large
//		ang = ((float)detId - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//		// current detector element position
//		curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//
//		// x-ray direction;
//		ray.d = normalize(curDetPos - ray.o);
//
//		//Only consider the case inside the FOV
//
//		int sliceId = 0;
//
//		for (sliceId = 0; sliceId != sliceNum; ++sliceId)
//		{
//			totWeight = 0;
//			curAng = draw[(sliceId * FanGeo.m_ViwN + angIdx) * FanGeo.m_DetN + detId] - calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
//				MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg + (sliceId * Img.m_Reso.x * Img.m_Reso.y), &totWeight);
//			if (!IS_ZERO(totWeight))
//			{
//				donp[(sliceId * numPerSubSet + prjIdx) * FanGeo.m_DetN + detId] = curAng / totWeight;
//			}
//			else
//			{
//				donp[(sliceId * numPerSubSet + prjIdx) * FanGeo.m_DetN + detId] = 0;
//			}
//		}
//
//	}
//}
//void proj(float* dimg, float* donp, float* draw, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid, cuint sliceNum, const cudaStream_t& streams)
//{
//	proj_Ker << <gid, blk, 0, streams >> >(dimg, donp, draw, FanGeo, Img, sliceNum, numPerSubSet, subSetNum, curSubSetIdx);
//}
//
//
//
//template<typename T>
//__global__ void _proj_AIM_Ker(T* dimg, T* dprj, T* draw, const FanEAGeo FanGeo, const Image Img, cuint sliceNum, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	int detIdx = threadIdx.x + blockIdx.x * blockDim.x;
//	int prjIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detIdx < FanGeo.m_DetN && prjIdx < numPerSubSet)
//	{
//		int angIdx = prjIdx * subSetNum + curSubSetIdx;
//		const T cntImgX = (static_cast<T>(Img.m_Reso.x) - 1) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
//		const T cntImgY = (static_cast<T>(Img.m_Reso.y) - 1) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
//		const T area = Img.m_Step.x * Img.m_Step.y;
//		T curAng = FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp;
//		while (curAng < 0){ curAng += (TWOPI); }
//		while (curAng > TWOPI){ curAng -= (TWOPI); }
//		T cosT = cos(curAng);
//		T sinT = sin(curAng);
//		T sour[2] = { -FanGeo.m_S2O * sinT, FanGeo.m_S2O * cosT };			
//		T SVA[3];
//		T SVB[3];
//		calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, detIdx);
//		unsigned int sliceIdx;
//		T summ[22]; //Usually used by the real Patient volume reconstruction
//		T coord;
//		T minC, maxC;
//		int maxIdx;
//		int curIdx;
//		T grid[4][3];
//		T dir[2];
//		T initDir[2];
//		T len;
//		T pdist, coef;
//		T curDirAng = (detIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp + curAng;
//
//		while (curDirAng < 0){ curDirAng += (TWOPI); }
//		while (curDirAng > TWOPI){ curDirAng -= (TWOPI); }
//		int levelIdx; // index different slices in Z direction
//		T weight = 0;
//		if (curDirAng <= PI / 4.0 || curDirAng > PI * 1.75)
//		{
//			for (levelIdx = 0; levelIdx != sliceNum; ++levelIdx)
//			{
//				summ[levelIdx] = 0;
//			}
//			weight = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.y; ++sliceIdx)
//			{
//				coord = (static_cast<T>(sliceIdx)-cntImgY - 0.5) * Img.m_Step.y;
//				minC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
//				maxC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.x + cntImgX) - 1;
//				maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + 1;
//				if (curIdx > static_cast<int>(Img.m_Reso.x - 1) || maxIdx < 0)
//				{
//
//					continue;
//				}
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//				if (maxIdx > Img.m_Reso.x - 1)
//				{
//					maxIdx = Img.m_Reso.x - 1;
//				}
//
//				for (; curIdx <= maxIdx; curIdx++)
//				{
//					grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
//					grid[0][1] = coord;
//					grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[1][1] = coord;
//					grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[2][1] = coord + Img.m_Step.y;
//					grid[3][0] = (curIdx - cntImgX - 0.5)* Img.m_Step.x;
//					grid[3][1] = coord + Img.m_Step.y;
//
//
//					//�����ĸ����Ӧ��ǰ��det index
//					dir[0] = grid[0][0] - sour[0];
//					dir[1] = grid[0][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//					//�����ĸ����Ӧ��ǰ��det index
//					dir[0] = grid[0][0] - sour[0];
//					dir[1] = grid[0][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//					dir[0] = grid[1][0] - sour[0];
//					dir[1] = grid[1][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					dir[0] = grid[2][0] - sour[0];
//					dir[1] = grid[2][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					dir[0] = grid[3][0] - sour[0];
//					dir[1] = grid[3][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					//����;
//					SortProjection<T>(grid);
//
//					pdist = hypot((static_cast<T>(curIdx)-cntImgX)*static_cast<T>(Img.m_Step.x) - sour[0], (static_cast<T>(sliceIdx)-cntImgY)*static_cast<T>(Img.m_Step.y) - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area) / (pdist * FanGeo.m_DetStp);
//
//					//�޸�;
//					weight += coef;
//					for (unsigned int levelIdx = 0; levelIdx != sliceNum; ++levelIdx)
//					{
//						summ[levelIdx] += coef *dimg[(levelIdx * Img.m_Reso.y + sliceIdx) * Img.m_Reso.x + curIdx];
//					}
//					//summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
//				}
//			}
//			if (!IS_ZERO<T>(weight))
//			{
//				for (levelIdx = 0; levelIdx < sliceNum; levelIdx++)
//				{
//					dprj[(levelIdx * numPerSubSet + prjIdx) * FanGeo.m_DetN + detIdx] = (draw[(levelIdx * FanGeo.m_ViwN + angIdx) * FanGeo.m_DetN + detIdx] - summ[levelIdx]) / weight;
//				}
//
//			}
//			else
//			{
//				for (levelIdx = 0; levelIdx < sliceNum; levelIdx++)
//				{
//					dprj[(levelIdx * numPerSubSet + prjIdx) * FanGeo.m_DetN + detIdx] = 0;
//				}
//			}
//			return;
//		} //End case 1
//		else if (curDirAng > PI * 0.25 && curDirAng <= PI * 0.75)
//		{
//			for (levelIdx = 0; levelIdx != sliceNum; ++levelIdx)
//			{
//				summ[levelIdx] = 0;
//			}
//			weight = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
//			{
//				coord = (sliceIdx - cntImgX - 0.5)* Img.m_Step.x;
//				
//				minC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
//				maxC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.y + cntImgY) - 1;
//				maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + 1;
//
//				if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0)
//				{
//					continue;
//				}
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//
//				if (maxIdx > Img.m_Reso.y - 1)
//				{
//					maxIdx = Img.m_Reso.y - 1;
//				}
//
//				for (; curIdx <= maxIdx; ++curIdx)
//				{
//					//��grid
//					grid[0][0] = coord;
//					grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[1][0] = coord + Img.m_Step.x;
//					grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[2][0] = coord + Img.m_Step.x;
//					grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//					grid[3][0] = coord;
//					grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//
//					//�����ĸ����Ӧ��ǰ��det index
//					dir[0] = grid[0][0] - sour[0];
//					dir[1] = grid[0][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//					dir[0] = grid[1][0] - sour[0];
//					dir[1] = grid[1][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					dir[0] = grid[2][0] - sour[0];
//					dir[1] = grid[2][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					dir[0] = grid[3][0] - sour[0];
//					dir[1] = grid[3][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					//����;
//					SortProjection<T>(grid);
//					pdist = hypot((static_cast<T>(sliceIdx)-cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(curIdx)-cntImgY)*Img.m_Step.y - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area) / (pdist * FanGeo.m_DetStp);
//					weight += coef;
//					for (levelIdx = 0; levelIdx != sliceNum; ++levelIdx)
//					{
//						summ[levelIdx] += coef * dimg[((levelIdx * Img.m_Reso.y + curIdx) * Img.m_Reso.x) + sliceIdx];
//					}
//
//
//				} //End for one slice
//			}// End all slices
//			if (!IS_ZERO<T>(weight))
//			{
//				for (levelIdx = 0; levelIdx < sliceNum; levelIdx++)
//				{
//					dprj[(levelIdx * numPerSubSet + prjIdx) * FanGeo.m_DetN + detIdx] = (draw[(levelIdx * FanGeo.m_ViwN + angIdx) * FanGeo.m_DetN + detIdx] - summ[levelIdx]) / weight;
//				}
//			}
//			else
//			{
//				for (levelIdx = 0; levelIdx < sliceNum; levelIdx++)
//				{
//					dprj[(levelIdx * numPerSubSet + prjIdx) * FanGeo.m_DetN + detIdx] = 0;
//				}
//			}
//			return;
//		}// End case 2
//		else if (curDirAng > PI * 0.75 && curDirAng <= PI * 1.25)
//		{
//			for (levelIdx = 0; levelIdx != sliceNum; ++levelIdx)
//			{
//				summ[levelIdx] = 0;
//			}
//			weight = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.y; ++sliceIdx)
//			{
//				coord = (sliceIdx - cntImgY + 0.5)* Img.m_Step.y;
//				//�������λ�ù����ཻ�����;
//				maxC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
//				minC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.x + cntImgX) - 1;
//				maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + 1;
//
//				if (curIdx > static_cast<int>(Img.m_Reso.x) || maxIdx < 0)
//				{
//					continue;
//				}
//
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//
//				if (maxIdx > Img.m_Reso.x - 1)
//				{
//					maxIdx = Img.m_Reso.x - 1;
//				}
//
//				for (; curIdx <= maxIdx; ++curIdx)
//				{
//					//��grid
//					grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
//					grid[0][1] = coord - Img.m_Step.y;
//					grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[1][1] = coord - Img.m_Step.y;
//					grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[2][1] = coord;
//					grid[3][0] = (curIdx - cntImgX - 0.5)* Img.m_Step.x;
//					grid[3][1] = coord;
//
//					//�����ĸ����Ӧ��ǰ��det index
//					dir[0] = grid[0][0] - sour[0];
//					dir[1] = grid[0][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//					dir[0] = grid[1][0] - sour[0];
//					dir[1] = grid[1][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					dir[0] = grid[2][0] - sour[0];
//					dir[1] = grid[2][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					dir[0] = grid[3][0] - sour[0];
//					dir[1] = grid[3][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					//����;
//					SortProjection<T>(grid);
//					pdist = hypot((static_cast<T>(curIdx)-cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(sliceIdx)-cntImgY)*Img.m_Step.y - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area) / (pdist * FanGeo.m_DetStp);
//					weight += coef;
//					for (levelIdx = 0; levelIdx != sliceNum; ++levelIdx)
//					{
//						summ[levelIdx] += coef * dimg[(levelIdx*Img.m_Reso.y + sliceIdx)* Img.m_Reso.x + curIdx];
//					}
//					//summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
//
//				} //End for one slice
//			}// End all slices
//			if (!IS_ZERO<T>(weight))
//			{
//				for (levelIdx = 0; levelIdx < sliceNum; levelIdx++)
//				{
//					dprj[(levelIdx * numPerSubSet + prjIdx) * FanGeo.m_DetN + detIdx] = (draw[(levelIdx * FanGeo.m_ViwN + angIdx) * FanGeo.m_DetN + detIdx] - summ[levelIdx]) / weight;
//				}
//			}
//			else
//			{
//				for (levelIdx = 0; levelIdx < sliceNum; levelIdx++)
//				{
//					dprj[(levelIdx * numPerSubSet + prjIdx) * FanGeo.m_DetN + detIdx] = 0;
//				}
//			}
//			return;
//		}
//		else //if (curDirAng > 1.25 * PI &&curDirAng <= 1.75*PI)
//		{
//			for (levelIdx = 0; levelIdx != sliceNum; ++levelIdx)
//			{
//				summ[levelIdx] = 0;
//			}
//			weight = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
//			{
//				coord = (sliceIdx - cntImgX + 0.5)* Img.m_Step.x;
//				//�������λ�ù����ཻ�����;
//				maxC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
//				minC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.y + cntImgY) - 1;
//				maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + 1;
//
//				if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0){ continue; }
//				if (curIdx < 0){ curIdx = 0; }
//				if (maxIdx > Img.m_Reso.y - 1)	{ maxIdx = Img.m_Reso.y - 1; }
//				for (; curIdx <= maxIdx; ++curIdx)
//				{
//					//��grid
//					grid[0][0] = coord - Img.m_Step.x;
//					grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[1][0] = coord;
//					grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[2][0] = coord - Img.m_Step.x;
//					grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//					grid[3][0] = coord;
//					grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//
//					//�����ĸ����Ӧ��ǰ��det index
//					dir[0] = grid[0][0] - sour[0];
//					dir[1] = grid[0][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//					dir[0] = grid[1][0] - sour[0];
//					dir[1] = grid[1][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					dir[0] = grid[2][0] - sour[0];
//					dir[1] = grid[2][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					dir[0] = grid[3][0] - sour[0];
//					dir[1] = grid[3][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					//����;
//					SortProjection<T>(grid);
//					pdist = hypot((static_cast<T>(sliceIdx)-cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(curIdx)-cntImgY)*Img.m_Step.y - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//					coef = coef / (pdist * FanGeo.m_DetStp);
//					weight += coef;
//					for (levelIdx = 0; levelIdx != sliceNum; ++levelIdx)
//					{
//						summ[levelIdx] += coef * dimg[(levelIdx * Img.m_Reso.y + curIdx) * Img.m_Reso.x + sliceIdx];
//					}
//				} //End for one slice
//			}// End all slices
//			if (!IS_ZERO<T>(weight))
//			{
//				for (levelIdx = 0; levelIdx < sliceNum; levelIdx++)
//				{
//					dprj[(levelIdx * numPerSubSet + prjIdx) * FanGeo.m_DetN + detIdx] = (draw[(levelIdx * FanGeo.m_ViwN + angIdx) * FanGeo.m_DetN + detIdx] - summ[levelIdx]) / weight;
//				}
//			}
//			else
//			{
//				for (levelIdx = 0; levelIdx < sliceNum; levelIdx++)
//				{
//					dprj[(levelIdx * numPerSubSet + prjIdx) * FanGeo.m_DetN + detIdx] = 0;
//				}
//			}
//			return;
//		}
//	}
//}
//
//
//template<typename T>
//void proj_AIM_temp(T* dimg, T* donp, T* draw, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid, cuint sliceNum, const cudaStream_t& streams)
//{
//	_proj_AIM_Ker<T> << <gid, blk, 0, streams >> >(dimg, donp, draw, FanGeo, Img, sliceNum, numPerSubSet, subSetNum, curSubSetIdx);
//}
//void proj_AIM(float* dimg, float* donp, float* draw, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid, cuint sliceNum, const cudaStream_t& streams)
//{
//	proj_AIM_temp<float>(dimg, donp, draw, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx, blk, gid, sliceNum, streams);
//}
//void proj_AIM(double* dimg, double* donp, double* draw, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid, cuint sliceNum, const cudaStream_t& streams)
//{
//	proj_AIM_temp<double>(dimg, donp, draw, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx, blk, gid, sliceNum, streams);
//}
//
//
//
//
//template<typename T>
//__global__ void _proj_AIMslow_ker(T* dprj, T* dimg, const FanEAGeo FanGeo, const Image Img)
//{
//	cuint detIdx = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint angIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detIdx < FanGeo.m_DetN && angIdx < FanGeo.m_ViwN)
//	{
//		const T cntImgX = (static_cast<T>(Img.m_Reso.x) - 1) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
//		const T cntImgY = (static_cast<T>(Img.m_Reso.y) - 1) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
//		const T area = Img.m_Step.x * Img.m_Step.y;
//		T curAng = FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp;
//		while (curAng < 0){ curAng += (TWOPI); }
//		while (curAng > TWOPI){ curAng -= (TWOPI); }
//		T cosT = cos(curAng);
//		T sinT = sin(curAng);
//		T sour[2] = { -FanGeo.m_S2O * sinT, FanGeo.m_S2O * cosT };
//		T SVA[3];
//		T SVB[3];
//		calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, detIdx);
//		T grid[4][3];
//		unsigned int imgLIdx = 0, imgWIdx = 0;
//		T summ = 0;
//		T coef = 0;
//		T pdist = 0;
//		T dir[2], initDir[2];
//		T len = 0;
//		for (imgWIdx = 0; imgWIdx < Img.m_Reso.y; imgWIdx++)
//		{
//			for (imgLIdx = 0; imgLIdx < Img.m_Reso.x; imgLIdx++)
//			{
//				grid[0][0] = (imgLIdx - cntImgX - 0.5) * Img.m_Step.x;
//				grid[0][1] = (imgWIdx - cntImgY - 0.5) * Img.m_Step.y;
//
//				dir[0] = grid[0][0] - sour[0];
//				dir[1] = grid[0][1] - sour[1];
//				len = hypot(dir[0], dir[1]);
//				dir[0] /= len;
//				dir[1] /= len;
//				dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//				dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//				initDir[0] = dir[0] * cosT + dir[1] * sinT;
//				initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//				grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//				grid[1][0] = grid[0][0] + Img.m_Step.x;
//				grid[1][1] = grid[0][1];
//				dir[0] = grid[1][0] - sour[0];
//				dir[1] = grid[1][1] - sour[1];
//				len = hypot(dir[0], dir[1]);
//				dir[0] /= len;
//				dir[1] /= len;
//				dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//				dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//				initDir[0] = dir[0] * cosT + dir[1] * sinT;
//				initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//				grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//				grid[2][0] = grid[1][0];
//				grid[2][1] = grid[1][1] + Img.m_Step.y;
//				dir[0] = grid[2][0] - sour[0];
//				dir[1] = grid[2][1] - sour[1];
//				len = hypot(dir[0], dir[1]);
//				dir[0] /= len;
//				dir[1] /= len;
//				dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//				dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//				initDir[0] = dir[0] * cosT + dir[1] * sinT;
//				initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//				grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//				grid[3][0] = grid[0][0];
//				grid[3][1] = grid[2][1];
//				dir[0] = grid[3][0] - sour[0];
//				dir[1] = grid[3][1] - sour[1];
//				len = hypot(dir[0], dir[1]);
//				dir[0] /= len;
//				dir[1] /= len;
//				dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//				dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//				initDir[0] = dir[0] * cosT + dir[1] * sinT;
//				initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//				grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//				//Sort grid
//				SortProjection<T>(grid);
//
//				pdist = hypot((imgLIdx - cntImgX)*Img.m_Step.x - sour[0], (imgWIdx - cntImgY)*Img.m_Step.y - sour[1]);
//				coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//				coef = coef / (pdist * FanGeo.m_DetStp);
//				summ += coef * dimg[imgWIdx* Img.m_Reso.x + imgLIdx];
//			}
//		}
//		dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
//	}
//}
//template<typename T>
//void proj_AIMslow_temp(T* dprj, T* dimg, const FanEAGeo& FanGeo, const Image&Img, const dim3& blk, const dim3& gid)
//{
//	_proj_AIMslow_ker<T> << <gid, blk >> >(dprj, dimg, FanGeo, Img);
//}
//void proj_AIMslow(float* dprj, float* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	proj_AIMslow_temp<float>(dprj, dimg, FanGeo, Img, blk, gid);
//}
//
//
//void proj_AIMslow(double* dprj, double* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	proj_AIMslow_temp<double>(dprj, dimg, FanGeo, Img, blk, gid);
//}
//
//
//
//
//template<typename T>
//__global__ void _proj_AIM_ker(T* dprj, T* dimg, const FanEAGeo FanGeo, const Image Img)
//{
//	int detIdx = threadIdx.x + blockIdx.x * blockDim.x;
//	int angIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detIdx < FanGeo.m_DetN && angIdx < FanGeo.m_ViwN)
//	{
//		const T cntImgX = (static_cast<T>(Img.m_Reso.x) - 1) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
//		const T cntImgY = (static_cast<T>(Img.m_Reso.y) - 1) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
//		const T area = Img.m_Step.x * Img.m_Step.y;
//		T curAng = FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp;
//		while (curAng < 0){ curAng += (TWOPI); }
//		while (curAng > TWOPI){ curAng -= (TWOPI); }
//		T cosT = cos(curAng);
//		T sinT = sin(curAng);
//		T sour[2] = { -FanGeo.m_S2O * sinT, FanGeo.m_S2O * cosT };
//		T SVA[3];
//		T SVB[3];
//		calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, detIdx);
//		unsigned int sliceIdx;
//		T summ;
//		T coord;
//		T minC, maxC;
//		int maxIdx;
//		int curIdx;
//		T grid[4][3];
//		T dir[2];
//		T initDir[2];
//		T len;
//		T pdist, coef;
//		T curDirAng = (detIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp + curAng;
//		while (curDirAng < 0){ curDirAng += (TWOPI); }
//		while (curDirAng > TWOPI){ curDirAng -= (TWOPI); }
//
//		if (curDirAng <= PI * 0.25 || curDirAng > PI * 1.75)
//		{
//			summ = 0;
//
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.y; ++sliceIdx)
//			{
//				coord = (static_cast<T>(sliceIdx)-cntImgY - 0.5) * Img.m_Step.y;
//				minC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
//				maxC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.x + cntImgX) - 1;
//				maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + 1;
//				if (curIdx > static_cast<int>(Img.m_Reso.x - 1) || maxIdx < 0)
//				{
//
//					continue;
//				}
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//				if (maxIdx > Img.m_Reso.x - 1)
//				{
//					maxIdx = Img.m_Reso.x - 1;
//				}
//
//				for (; curIdx <= maxIdx; curIdx++)
//				{
//					grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
//					grid[0][1] = coord;
//					grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[1][1] = coord;
//					grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[2][1] = coord + Img.m_Step.y;
//					grid[3][0] = (curIdx - cntImgX - 0.5)* Img.m_Step.x;
//					grid[3][1] = coord + Img.m_Step.y;
//
//
//					//�����ĸ����Ӧ��ǰ��det index
//					dir[0] = grid[0][0] - sour[0];
//					dir[1] = grid[0][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//					//�����ĸ����Ӧ��ǰ��det index
//					dir[0] = grid[0][0] - sour[0];
//					dir[1] = grid[0][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//					dir[0] = grid[1][0] - sour[0];
//					dir[1] = grid[1][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					dir[0] = grid[2][0] - sour[0];
//					dir[1] = grid[2][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					dir[0] = grid[3][0] - sour[0];
//					dir[1] = grid[3][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					SortProjection<T>(grid);
//
//					pdist = hypot((static_cast<T>(curIdx)-cntImgX)*static_cast<T>(Img.m_Step.x) - sour[0], (static_cast<T>(sliceIdx)-cntImgY)*static_cast<T>(Img.m_Step.y) - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area) / (pdist * FanGeo.m_DetStp);
//					summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
//				}
//			}
//			dprj[angIdx* FanGeo.m_DetN + detIdx] = summ;
//			return;
//		} //End case 1
//		else if (curDirAng > PI * 0.25 && curDirAng <= PI * 0.75)
//		{
//			summ = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
//			{
//				coord = (sliceIdx - cntImgX - 0.5)* Img.m_Step.x;
//				//�������λ�ù����ཻ�����;
//				minC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
//				maxC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.y + cntImgY) - 1;
//				maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + 1;
//
//				if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0)
//				{
//					continue;
//				}
//
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//
//				if (maxIdx > Img.m_Reso.y - 1)
//				{
//					maxIdx = Img.m_Reso.y - 1;
//				}
//
//				for (; curIdx <= maxIdx; ++curIdx)
//				{
//					//��grid
//					grid[0][0] = coord;
//					grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[1][0] = coord + Img.m_Step.x;
//					grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[2][0] = coord + Img.m_Step.x;
//					grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//					grid[3][0] = coord;
//					grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//
//					//�����ĸ����Ӧ��ǰ��det index
//					dir[0] = grid[0][0] - sour[0];
//					dir[1] = grid[0][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//					dir[0] = grid[1][0] - sour[0];
//					dir[1] = grid[1][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					dir[0] = grid[2][0] - sour[0];
//					dir[1] = grid[2][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					dir[0] = grid[3][0] - sour[0];
//					dir[1] = grid[3][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					SortProjection<T>(grid);
//					pdist = hypot((static_cast<T>(sliceIdx)-cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(curIdx)-cntImgY)*Img.m_Step.y - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area) / (pdist * FanGeo.m_DetStp);
//					summ += coef * dimg[curIdx * Img.m_Reso.x + sliceIdx];
//
//				} //End for one slice
//			}// End all slices
//			dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
//			return;
//		}// End case 2
//		else if (curDirAng > PI * 0.75 && curDirAng <= PI * 1.25)
//		{
//			summ = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.y; ++sliceIdx)
//			{
//				coord = (sliceIdx - cntImgY + 0.5)* Img.m_Step.y;
//				//�������λ�ù����ཻ�����;
//				maxC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
//				minC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.x + cntImgX) - 1;
//				maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + 1;
//
//				if (curIdx > static_cast<int>(Img.m_Reso.x) || maxIdx < 0)
//				{
//					continue;
//				}
//
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//
//				if (maxIdx > Img.m_Reso.x - 1)
//				{
//					maxIdx = Img.m_Reso.x - 1;
//				}
//
//				for (; curIdx <= maxIdx; ++curIdx)
//				{
//					//��grid
//					grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
//					grid[0][1] = coord - Img.m_Step.y;
//					grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[1][1] = coord - Img.m_Step.y;
//					grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[2][1] = coord;
//					grid[3][0] = (curIdx - cntImgX - 0.5)* Img.m_Step.x;
//					grid[3][1] = coord;
//
//					//�����ĸ����Ӧ��ǰ��det index
//					dir[0] = grid[0][0] - sour[0];
//					dir[1] = grid[0][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//					dir[0] = grid[1][0] - sour[0];
//					dir[1] = grid[1][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					dir[0] = grid[2][0] - sour[0];
//					dir[1] = grid[2][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					dir[0] = grid[3][0] - sour[0];
//					dir[1] = grid[3][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					SortProjection<T>(grid);
//					pdist = hypot((static_cast<T>(curIdx)-cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(sliceIdx)-cntImgY)*Img.m_Step.y - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area) / (pdist * FanGeo.m_DetStp);
//					summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
//
//				} //End for one slice
//			}// End all slices
//			dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
//			return;
//		}
//		else //if (curDirAng > 1.25 * PI &&curDirAng <= 1.75*PI)
//		{
//			summ = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
//			{
//				coord = (sliceIdx - cntImgX + 0.5)* Img.m_Step.x;
//				//�������λ�ù����ཻ�����;
//				maxC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
//				minC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.y + cntImgY) - 1;
//				maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + 1;
//
//				if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0)
//				{
//					continue;
//				}
//
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//
//				if (maxIdx > Img.m_Reso.y - 1)
//				{
//					maxIdx = Img.m_Reso.y - 1;
//				}
//
//				for (; curIdx <= maxIdx; ++curIdx)
//				{
//					//��grid
//					grid[0][0] = coord - Img.m_Step.x;
//					grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[1][0] = coord;
//					grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[2][0] = coord - Img.m_Step.x;
//					grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//					grid[3][0] = coord;
//					grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//
//					//�����ĸ����Ӧ��ǰ��det index
//					dir[0] = grid[0][0] - sour[0];
//					dir[1] = grid[0][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//					dir[0] = grid[1][0] - sour[0];
//					dir[1] = grid[1][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					dir[0] = grid[2][0] - sour[0];
//					dir[1] = grid[2][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					dir[0] = grid[3][0] - sour[0];
//					dir[1] = grid[3][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					//����;
//					SortProjection<T>(grid);
//					pdist = hypot((static_cast<T>(sliceIdx)-cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(curIdx)-cntImgY)*Img.m_Step.y - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//					coef = coef / (pdist * FanGeo.m_DetStp);
//					summ += coef * dimg[curIdx* Img.m_Reso.x + sliceIdx];
//					//summ += *coef;
//				} //End for one slice
//			}// End all slices
//			dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
//			return;
//		}
//	}
//}
//template<typename T>
//inline void proj_AIM_temp(T* dprj, T* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_proj_AIM_ker<T> << <gid, blk >> >(dprj, dimg, FanGeo, Img);
//}
//void proj_AIM(float* dprj, float* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3&blk, const dim3& gid)
//{
//	proj_AIM_temp<float>(dprj, dimg, FanGeo, Img, blk, gid);
//}
//void proj_AIM(double* dprj, double* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3&blk, const dim3& gid)
//{
//	proj_AIM_temp<double>(dprj, dimg, FanGeo, Img, blk, gid);
//}
//
//
//
//
//
//
//
//
////����ȷ
//template<typename T>
//__global__ void _proj_AIM_ker(T* dprj, T* draw, T* dimg, const FanEAGeo FanGeo, const Image Img)
//{
//	int detIdx = threadIdx.x + blockIdx.x * blockDim.x;
//	int angIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detIdx < FanGeo.m_DetN && angIdx < FanGeo.m_ViwN)
//	{
//		const T cntImgX = (static_cast<T>(Img.m_Reso.x) - 1) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
//		const T cntImgY = (static_cast<T>(Img.m_Reso.y) - 1) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
//		const T area = Img.m_Step.x * Img.m_Step.y;
//		T curAng = FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp;
//		while (curAng < 0){ curAng += (TWOPI); }
//		while (curAng > TWOPI){ curAng -= (TWOPI); }
//		T cosT = cos(curAng);
//		T sinT = sin(curAng);
//		T sour[2] = { -FanGeo.m_S2O * sinT, FanGeo.m_S2O * cosT };
//		T SVA[3];
//		T SVB[3];
//		calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, detIdx);
//		unsigned int sliceIdx;
//		T summ;
//		T coord;
//		T minC, maxC;
//		int maxIdx;
//		int curIdx;
//		T grid[4][3];
//		T dir[2];
//		T initDir[2];
//		T len;
//		T pdist, coef;
//		T curDirAng = (detIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp + curAng;
//
//		while (curDirAng < 0){ curDirAng += (TWOPI); }
//		while (curDirAng > TWOPI){ curDirAng -= (TWOPI); }
//
//		T weight = 0;
//		if (curDirAng <= PI * 0.25 || curDirAng > PI * 1.75)
//		{
//			summ = 0;
//			weight = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.y; ++sliceIdx)
//			{
//				coord = (static_cast<T>(sliceIdx)-cntImgY - 0.5) * Img.m_Step.y;
//				minC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
//				maxC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.x + cntImgX) - 1;
//				maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + 1;
//				if (curIdx > static_cast<int>(Img.m_Reso.x - 1) || maxIdx < 0)
//				{
//
//					continue;
//				}
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//				if (maxIdx > Img.m_Reso.x - 1)
//				{
//					maxIdx = Img.m_Reso.x - 1;
//				}
//
//				for (; curIdx <= maxIdx; curIdx++)
//				{
//					grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
//					grid[0][1] = coord;
//					grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[1][1] = coord;
//					grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[2][1] = coord + Img.m_Step.y;
//					grid[3][0] = (curIdx - cntImgX - 0.5)* Img.m_Step.x;
//					grid[3][1] = coord + Img.m_Step.y;
//
//
//					//�����ĸ����Ӧ��ǰ��det index
//					dir[0] = grid[0][0] - sour[0];
//					dir[1] = grid[0][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//					//�����ĸ����Ӧ��ǰ��det index
//					dir[0] = grid[0][0] - sour[0];
//					dir[1] = grid[0][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//					dir[0] = grid[1][0] - sour[0];
//					dir[1] = grid[1][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					dir[0] = grid[2][0] - sour[0];
//					dir[1] = grid[2][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					dir[0] = grid[3][0] - sour[0];
//					dir[1] = grid[3][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					SortProjection<T>(grid);
//
//					pdist = hypot((static_cast<T>(curIdx)-cntImgX)*static_cast<T>(Img.m_Step.x) - sour[0], (static_cast<T>(sliceIdx)-cntImgY)*static_cast<T>(Img.m_Step.y) - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area) / (pdist * FanGeo.m_DetStp);
//					weight += coef;
//					summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
//				}
//			}
//			if (!IS_ZERO<T>(weight))
//			{
//				dprj[angIdx* FanGeo.m_DetN + detIdx] = (draw[angIdx* FanGeo.m_DetN + detIdx] - summ) / weight;
//			}
//			else
//			{
//				dprj[angIdx* FanGeo.m_DetN + detIdx] = 0;
//			}
//
//			return;
//		} //End case 1
//		else if (curDirAng > PI * 0.25 && curDirAng <= PI * 0.75)
//		{
//			summ = 0;
//			weight = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
//			{
//				coord = (sliceIdx - cntImgX - 0.5)* Img.m_Step.x;
//				//�������λ�ù����ཻ�����;
//				minC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
//				maxC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.y + cntImgY) - 1;
//				maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + 1;
//
//				if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0)
//				{
//					continue;
//				}
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//
//				if (maxIdx > Img.m_Reso.y - 1)
//				{
//					maxIdx = Img.m_Reso.y - 1;
//				}
//
//				for (; curIdx <= maxIdx; ++curIdx)
//				{
//					//��grid
//					grid[0][0] = coord;
//					grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[1][0] = coord + Img.m_Step.x;
//					grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[2][0] = coord + Img.m_Step.x;
//					grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//					grid[3][0] = coord;
//					grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//
//					//�����ĸ����Ӧ��ǰ��det index
//					dir[0] = grid[0][0] - sour[0];
//					dir[1] = grid[0][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//					dir[0] = grid[1][0] - sour[0];
//					dir[1] = grid[1][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					dir[0] = grid[2][0] - sour[0];
//					dir[1] = grid[2][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					dir[0] = grid[3][0] - sour[0];
//					dir[1] = grid[3][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					SortProjection<T>(grid);
//					pdist = hypot((static_cast<T>(sliceIdx)-cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(curIdx)-cntImgY)*Img.m_Step.y - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area) / (pdist * FanGeo.m_DetStp);
//					weight += coef;
//					summ += coef * dimg[curIdx * Img.m_Reso.x + sliceIdx];
//
//				} //End for one slice
//			}// End all slices
//			if (!IS_ZERO<T>(weight))
//			{
//				dprj[angIdx* FanGeo.m_DetN + detIdx] = (draw[angIdx* FanGeo.m_DetN + detIdx] - summ) / weight;
//			}
//			else
//			{
//				dprj[angIdx* FanGeo.m_DetN + detIdx] = 0;
//			}
//			return;
//		}// End case 2
//		else if (curDirAng > PI * 0.75 && curDirAng <= PI * 1.25)
//		{
//			summ = 0;
//			weight = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.y; ++sliceIdx)
//			{
//				coord = (sliceIdx - cntImgY + 0.5)* Img.m_Step.y;
//				//�������λ�ù����ཻ�����;
//				maxC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
//				minC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.x + cntImgX) - 1;
//				maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + 1;
//
//				if (curIdx > static_cast<int>(Img.m_Reso.x) || maxIdx < 0)
//				{
//					continue;
//				}
//
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//
//				if (maxIdx > Img.m_Reso.x - 1)
//				{
//					maxIdx = Img.m_Reso.x - 1;
//				}
//
//				for (; curIdx <= maxIdx; ++curIdx)
//				{
//					//��grid
//					grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
//					grid[0][1] = coord - Img.m_Step.y;
//					grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[1][1] = coord - Img.m_Step.y;
//					grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[2][1] = coord;
//					grid[3][0] = (curIdx - cntImgX - 0.5)* Img.m_Step.x;
//					grid[3][1] = coord;
//
//					//�����ĸ����Ӧ��ǰ��det index
//					dir[0] = grid[0][0] - sour[0];
//					dir[1] = grid[0][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//					dir[0] = grid[1][0] - sour[0];
//					dir[1] = grid[1][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					dir[0] = grid[2][0] - sour[0];
//					dir[1] = grid[2][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					dir[0] = grid[3][0] - sour[0];
//					dir[1] = grid[3][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					SortProjection<T>(grid);
//					pdist = hypot((static_cast<T>(curIdx)-cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(sliceIdx)-cntImgY)*Img.m_Step.y - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area) / (pdist * FanGeo.m_DetStp);
//					weight += coef;
//					summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
//
//				} //End for one slice
//			}// End all slices
//			if (!IS_ZERO<T>(weight))
//			{
//				dprj[angIdx* FanGeo.m_DetN + detIdx] = (draw[angIdx* FanGeo.m_DetN + detIdx] - summ) / weight;
//			}
//			else
//			{
//				dprj[angIdx* FanGeo.m_DetN + detIdx] = 0;
//			}
//			return;
//		}
//		else //if (curDirAng > 1.25 * PI &&curDirAng <= 1.75*PI)
//		{
//			summ = 0;
//			weight = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
//			{
//				coord = (sliceIdx - cntImgX + 0.5)* Img.m_Step.x;
//				//�������λ�ù����ཻ�����;
//				maxC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
//				minC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.y + cntImgY) - 1;
//				maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + 1;
//
//				if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0)
//				{
//					continue;
//				}
//
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//
//				if (maxIdx > Img.m_Reso.y - 1)
//				{
//					maxIdx = Img.m_Reso.y - 1;
//				}
//
//				for (; curIdx <= maxIdx; ++curIdx)
//				{
//					//��grid
//					grid[0][0] = coord - Img.m_Step.x;
//					grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[1][0] = coord;
//					grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[2][0] = coord - Img.m_Step.x;
//					grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//					grid[3][0] = coord;
//					grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//
//					//�����ĸ����Ӧ��ǰ��det index
//					dir[0] = grid[0][0] - sour[0];
//					dir[1] = grid[0][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//					dir[0] = grid[1][0] - sour[0];
//					dir[1] = grid[1][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					dir[0] = grid[2][0] - sour[0];
//					dir[1] = grid[2][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					dir[0] = grid[3][0] - sour[0];
//					dir[1] = grid[3][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					SortProjection<T>(grid);
//					pdist = hypot((static_cast<T>(sliceIdx)-cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(curIdx)-cntImgY)*Img.m_Step.y - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//					coef = coef / (pdist * FanGeo.m_DetStp);
//					weight += coef;
//					summ += coef * dimg[curIdx* Img.m_Reso.x + sliceIdx];
//					//summ += *coef;
//				} //End for one slice
//			}// End all slices
//			if (!IS_ZERO<T>(weight))
//			{
//				dprj[angIdx* FanGeo.m_DetN + detIdx] = (draw[angIdx* FanGeo.m_DetN + detIdx] - summ) / weight;
//			}
//			else
//			{
//				dprj[angIdx* FanGeo.m_DetN + detIdx] = 0;
//			}
//			return;
//		}
//	}
//}
//template<typename T>
//void proj_AIM_temp(T* dprj, T* draw, T* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_proj_AIM_ker<T> << <gid, blk >> >(dprj, draw, dimg, FanGeo, Img);
//}
//void proj_AIM(float* dprj, float* draw, float* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3&blk, const dim3& gid)
//{
//	proj_AIM_temp<float>(dprj, draw, dimg, FanGeo, Img, blk, gid);
//}
//
//void proj_AIM(double* dprj, double* draw, double* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3&blk, const dim3& gid)
//{
//	proj_AIM_temp<double>(dprj, draw, dimg, FanGeo, Img, blk, gid);
//}
//
//
//
//template<typename T>
//__global__ void _proj_AIM_ker(T* dimg, T* dprj, T* draw, const FanEAGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	int detIdx = threadIdx.x + blockIdx.x * blockDim.x;
//	int prjIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detIdx < FanGeo.m_DetN && prjIdx < numPerSubSet)
//	{
//		int angIdx = prjIdx * subSetNum + curSubSetIdx;
//		const T cntImgX = (static_cast<T>(Img.m_Reso.x) - 1) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
//		const T cntImgY = (static_cast<T>(Img.m_Reso.y) - 1) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
//		const T area = Img.m_Step.x * Img.m_Step.y;
//		T curAng = FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp;
//		while (curAng < 0){ curAng += (TWOPI); }
//		while (curAng > TWOPI){ curAng -= (TWOPI); }
//		T cosT = cos(curAng);
//		T sinT = sin(curAng);
//		T sour[2] = { -FanGeo.m_S2O * sinT, FanGeo.m_S2O * cosT };
//		T SVA[3];
//		T SVB[3];
//		calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, detIdx);
//		unsigned int sliceIdx;
//		T summ;
//		T coord;
//		T minC, maxC;
//		int maxIdx;
//		int curIdx;
//		T grid[4][3];
//		T dir[2];
//		T initDir[2];
//		T len;
//		T pdist, coef;
//		T curDirAng = (detIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp + curAng;
//
//		while (curDirAng < 0){ curDirAng += (TWOPI); }
//		while (curDirAng > TWOPI){ curDirAng -= (TWOPI); }
//
//		T weight = 0;
//		if (curDirAng <= PI * 0.25 || curDirAng > PI * 1.75)
//		{
//			summ = 0;
//			weight = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.y; ++sliceIdx)
//			{
//				coord = (static_cast<T>(sliceIdx)-cntImgY - 0.5) * Img.m_Step.y;
//				minC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
//				maxC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.x + cntImgX) - 1;
//				maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + 1;
//				if (curIdx > static_cast<int>(Img.m_Reso.x - 1) || maxIdx < 0)
//				{
//
//					continue;
//				}
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//				if (maxIdx > Img.m_Reso.x - 1)
//				{
//					maxIdx = Img.m_Reso.x - 1;
//				}
//
//				for (; curIdx <= maxIdx; curIdx++)
//				{
//					grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
//					grid[0][1] = coord;
//					grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[1][1] = coord;
//					grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[2][1] = coord + Img.m_Step.y;
//					grid[3][0] = (curIdx - cntImgX - 0.5)* Img.m_Step.x;
//					grid[3][1] = coord + Img.m_Step.y;
//
//
//					//�����ĸ����Ӧ��ǰ��det index
//					dir[0] = grid[0][0] - sour[0];
//					dir[1] = grid[0][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//					//�����ĸ����Ӧ��ǰ��det index
//					dir[0] = grid[0][0] - sour[0];
//					dir[1] = grid[0][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//					dir[0] = grid[1][0] - sour[0];
//					dir[1] = grid[1][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					dir[0] = grid[2][0] - sour[0];
//					dir[1] = grid[2][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					dir[0] = grid[3][0] - sour[0];
//					dir[1] = grid[3][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					SortProjection<T>(grid);
//
//					pdist = hypot((static_cast<T>(curIdx)-cntImgX)*static_cast<T>(Img.m_Step.x) - sour[0], (static_cast<T>(sliceIdx)-cntImgY)*static_cast<T>(Img.m_Step.y) - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area) / (pdist * FanGeo.m_DetStp);
//					weight += coef;
//					summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
//				}
//			}
//			if (!IS_ZERO<T>(weight))
//			{
//				dprj[prjIdx* FanGeo.m_DetN + detIdx] = (draw[angIdx* FanGeo.m_DetN + detIdx] - summ) / weight;
//			}
//			else
//			{
//				dprj[prjIdx* FanGeo.m_DetN + detIdx] = 0;
//			}
//			return;
//		} //End case 1
//		else if (curDirAng > PI * 0.25 && curDirAng <= PI * 0.75)
//		{
//			summ = 0;
//			weight = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
//			{
//				coord = (sliceIdx - cntImgX - 0.5)* Img.m_Step.x;
//				//�������λ�ù����ཻ�����;
//				minC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
//				maxC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.y + cntImgY) - 1;
//				maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + 1;
//
//				if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0)
//				{
//					continue;
//				}
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//
//				if (maxIdx > Img.m_Reso.y - 1)
//				{
//					maxIdx = Img.m_Reso.y - 1;
//				}
//
//				for (; curIdx <= maxIdx; ++curIdx)
//				{
//					//��grid
//					grid[0][0] = coord;
//					grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[1][0] = coord + Img.m_Step.x;
//					grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[2][0] = coord + Img.m_Step.x;
//					grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//					grid[3][0] = coord;
//					grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//
//					//�����ĸ����Ӧ��ǰ��det index
//					dir[0] = grid[0][0] - sour[0];
//					dir[1] = grid[0][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//					dir[0] = grid[1][0] - sour[0];
//					dir[1] = grid[1][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					dir[0] = grid[2][0] - sour[0];
//					dir[1] = grid[2][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					dir[0] = grid[3][0] - sour[0];
//					dir[1] = grid[3][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					SortProjection<T>(grid);
//					pdist = hypot((static_cast<T>(sliceIdx)-cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(curIdx)-cntImgY)*Img.m_Step.y - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area) / (pdist * FanGeo.m_DetStp);
//					weight += coef;
//					summ += coef * dimg[curIdx * Img.m_Reso.x + sliceIdx];
//
//				} //End for one slice
//			}// End all slices
//			if (!IS_ZERO<T>(weight))
//			{
//				dprj[prjIdx* FanGeo.m_DetN + detIdx] = (draw[angIdx* FanGeo.m_DetN + detIdx] - summ) / weight;
//			}
//			else
//			{
//				dprj[prjIdx* FanGeo.m_DetN + detIdx] = 0;
//			}
//			return;
//		}// End case 2
//		else if (curDirAng > PI * 0.75 && curDirAng <= PI * 1.25)
//		{
//			summ = 0;
//			weight = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.y; ++sliceIdx)
//			{
//				coord = (sliceIdx - cntImgY + 0.5)* Img.m_Step.y;
//				//�������λ�ù����ཻ�����;
//				maxC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
//				minC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.x + cntImgX) - 1;
//				maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + 1;
//
//				if (curIdx > static_cast<int>(Img.m_Reso.x) || maxIdx < 0)
//				{
//					continue;
//				}
//
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//
//				if (maxIdx > Img.m_Reso.x - 1)
//				{
//					maxIdx = Img.m_Reso.x - 1;
//				}
//
//				for (; curIdx <= maxIdx; ++curIdx)
//				{
//					//��grid
//					grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
//					grid[0][1] = coord - Img.m_Step.y;
//					grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[1][1] = coord - Img.m_Step.y;
//					grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[2][1] = coord;
//					grid[3][0] = (curIdx - cntImgX - 0.5)* Img.m_Step.x;
//					grid[3][1] = coord;
//
//					//�����ĸ����Ӧ��ǰ��det index
//					dir[0] = grid[0][0] - sour[0];
//					dir[1] = grid[0][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//					dir[0] = grid[1][0] - sour[0];
//					dir[1] = grid[1][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					dir[0] = grid[2][0] - sour[0];
//					dir[1] = grid[2][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					dir[0] = grid[3][0] - sour[0];
//					dir[1] = grid[3][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					SortProjection<T>(grid);
//					pdist = hypot((static_cast<T>(curIdx)-cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(sliceIdx)-cntImgY)*Img.m_Step.y - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area) / (pdist * FanGeo.m_DetStp);
//					weight += coef;
//					summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
//
//				} //End for one slice
//			}// End all slices
//			if (!IS_ZERO<T>(weight))
//			{
//				dprj[prjIdx* FanGeo.m_DetN + detIdx] = (draw[angIdx* FanGeo.m_DetN + detIdx] - summ) / weight;
//			}
//			else
//			{
//				dprj[prjIdx* FanGeo.m_DetN + detIdx] = 0;
//			}
//			return;
//		}
//		else //if (curDirAng > 1.25 * PI &&curDirAng <= 1.75*PI)
//		{
//			summ = 0;
//			weight = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
//			{
//				coord = (sliceIdx - cntImgX + 0.5)* Img.m_Step.x;
//				//�������λ�ù����ཻ�����;
//				maxC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
//				minC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.y + cntImgY) - 1;
//				maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + 1;
//
//				if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0)
//				{
//					continue;
//				}
//
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//
//				if (maxIdx > Img.m_Reso.y - 1)
//				{
//					maxIdx = Img.m_Reso.y - 1;
//				}
//
//				for (; curIdx <= maxIdx; ++curIdx)
//				{
//					//��grid
//					grid[0][0] = coord - Img.m_Step.x;
//					grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[1][0] = coord;
//					grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[2][0] = coord - Img.m_Step.x;
//					grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//					grid[3][0] = coord;
//					grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//
//					//�����ĸ����Ӧ��ǰ��det index
//					dir[0] = grid[0][0] - sour[0];
//					dir[1] = grid[0][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//					dir[0] = grid[1][0] - sour[0];
//					dir[1] = grid[1][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					dir[0] = grid[2][0] - sour[0];
//					dir[1] = grid[2][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					dir[0] = grid[3][0] - sour[0];
//					dir[1] = grid[3][1] - sour[1];
//					len = hypot(dir[0], dir[1]);
//					dir[0] /= len;
//					dir[1] /= len;
//					dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//					dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//					initDir[0] = dir[0] * cosT + dir[1] * sinT;
//					initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//					grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					SortProjection<T>(grid);
//					pdist = hypot((static_cast<T>(sliceIdx)-cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(curIdx)-cntImgY)*Img.m_Step.y - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//					coef = coef / (pdist * FanGeo.m_DetStp);
//					weight += coef;
//					summ += coef * dimg[curIdx* Img.m_Reso.x + sliceIdx];
//					//summ += *coef;
//				} //End for one slice
//			}// End all slices
//			if (!IS_ZERO<T>(weight))
//			{
//				dprj[prjIdx* FanGeo.m_DetN + detIdx] = (draw[angIdx* FanGeo.m_DetN + detIdx] - summ) / weight;
//			}
//			else
//			{
//				dprj[prjIdx* FanGeo.m_DetN + detIdx] = 0;
//			}
//			return;
//		}
//	}
//
//}
//
//template<typename T>
//void proj_AIM_temp(T* dimg, T* donp, T* draw, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_proj_AIM_ker<T> << <gid, blk >> >(dimg, donp, draw, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx);
//}
//
//void proj_AIM(float* dimg, float* donp, float* draw, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	proj_AIM_temp<float>(dimg, donp, draw, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx, blk, gid);
//}
//void proj_AIM(double* dimg, double* donp, double* draw, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	proj_AIM_temp<double>(dimg, donp, draw, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx, blk, gid);
//}
//
//
//
//
//
//
//
////ÿ��block��һ���Ƕȵ�ͶӰ�����нǶȵ�ͶӰ�����ܵ�block����,ÿ���̶߳�Ӧһ��detector�߽�;
//template<typename T>
//__global__ void _proj_AIM2_ker(T* dprj, T* dimg, const FanEAGeo FanGeo, const Image Img)
//{
//	const int detIdx = threadIdx.x;
//	const int angIdx = blockIdx.x;
//
//
//	T SVA[3];
//	T SVB[3];
//	T curAng = FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp;
//	while (curAng < 0){ curAng += (TWOPI); }
//	while (curAng > (TWOPI)){ curAng -= (TWOPI); }
//
//	int sliceIdx;
//	T summ;
//	T coord;
//
//	int maxIdx;
//	int curIdx;
//	T grid[4][3];
//	T dir[2];
//	T initDir[2];
//	T len;
//	T coef;
//	const T area = Img.m_Step.x * Img.m_Step.y;
//	const T cosT = cos(curAng);
//	const T sinT = sin(curAng);
//	const T cntImgX = (static_cast<T>(Img.m_Reso.x) - 1) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
//	const T cntImgY = (static_cast<T>(Img.m_Reso.y) - 1) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
//	T sour[2];
//	sour[0] = -FanGeo.m_S2O * sinT;
//	sour[1] = FanGeo.m_S2O * cosT;
//	__shared__ T SV[1024][3];
//
//	SVA[0] = (detIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//	SVA[1] = -cos(SVA[0]) * FanGeo.m_S2D + FanGeo.m_S2O;
//	SVA[2] = sin(SVA[0]) * FanGeo.m_S2D;
//	SVB[0] = SVA[2] * cosT - SVA[1] * sinT;
//	SVB[1] = SVA[2] * sinT + SVA[1] * cosT;
//	SVA[2] = SVB[0] - sour[0];
//	SVA[1] = SVB[1] - sour[1];
//	SVB[2] = sqrt(SVA[2] * SVA[2] + SVA[1] * SVA[1]);
//	SV[detIdx][0] = SVA[2] / SVB[2];
//	SV[detIdx][1] = SVA[1] / SVB[2];
//	SV[detIdx][2] = static_cast<T>(detIdx);
//
//	__syncthreads();
//
//	if (detIdx >= FanGeo.m_DetN) return;
//
//	SVA[0] = SV[detIdx][0];
//	SVA[1] = SV[detIdx][1];
//	SVA[2] = SV[detIdx][2];
//	SVB[0] = SV[detIdx + 1][0];
//	SVB[1] = SV[detIdx + 1][1];
//	SVB[2] = SV[detIdx + 1][2];
//
//
//	if (curAng <= (PI * 0.25) || curAng > (PI * 1.75))
//	{
//		summ = 0;
//
//		for (sliceIdx = 0; sliceIdx < Img.m_Reso.y; sliceIdx++)
//		{
//			coord = (sliceIdx - cntImgY - 0.5) * Img.m_Step.y;
//			initDir[0] = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
//			initDir[1] = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
//
//			curIdx = int(initDir[0] / Img.m_Step.x + cntImgX - 1);
//			maxIdx = ceilf(initDir[1] / Img.m_Step.x + cntImgX + 1);
//			if (maxIdx >Img.m_Reso.x)
//			{
//				continue;
//			}
//			if (curIdx < 0)
//			{
//				continue;
//			}
//			if (curIdx < 0)
//			{
//				curIdx = 0;
//			}
//			if (maxIdx > Img.m_Reso.x)
//			{
//				maxIdx = Img.m_Reso.x;
//			}
//
//			for (; curIdx < maxIdx; curIdx++)
//			{
//				grid[0][0] = (curIdx - cntImgX) * Img.m_Step.x;
//				grid[0][1] = coord;
//				grid[1][0] = (curIdx - cntImgX + 1) * Img.m_Step.x;
//				grid[1][1] = coord;
//				grid[2][0] = (curIdx - cntImgX + 1) * Img.m_Step.x;
//				grid[2][1] = coord + Img.m_Step.y;
//				grid[3][0] = (curIdx - cntImgX)* Img.m_Step.x;
//				grid[3][1] = coord + Img.m_Step.y;
//
//
//				//�����ĸ����Ӧ��ǰ��det index
//				dir[0] = grid[0][0] - sour[0];
//				dir[1] = grid[0][1] - sour[1];
//				len = hypot(dir[0], dir[1]);
//				dir[0] /= len;
//				dir[1] /= len;
//				dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//				dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//				initDir[0] = dir[0] * cosT + dir[1] * sinT;
//				initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//				grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//				//�����ĸ����Ӧ��ǰ��det index
//				dir[0] = grid[0][0] - sour[0];
//				dir[1] = grid[0][1] - sour[1];
//				len = hypot(dir[0], dir[1]);
//				dir[0] /= len;
//				dir[1] /= len;
//				dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//				dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//				initDir[0] = dir[0] * cosT + dir[1] * sinT;
//				initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//				grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//				dir[0] = grid[1][0] - sour[0];
//				dir[1] = grid[1][1] - sour[1];
//				len = hypot(dir[0], dir[1]);
//				dir[0] /= len;
//				dir[1] /= len;
//				dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//				dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//				initDir[0] = dir[0] * cosT + dir[1] * sinT;
//				initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//				grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//				dir[0] = grid[2][0] - sour[0];
//				dir[1] = grid[2][1] - sour[1];
//				len = hypot(dir[0], dir[1]);
//				dir[0] /= len;
//				dir[1] /= len;
//				dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//				dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//				initDir[0] = dir[0] * cosT + dir[1] * sinT;
//				initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//				grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//				dir[0] = grid[3][0] - sour[0];
//				dir[1] = grid[3][1] - sour[1];
//				len = hypot(dir[0], dir[1]);
//				dir[0] /= len;
//				dir[1] /= len;
//				dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//				dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//				initDir[0] = dir[0] * cosT + dir[1] * sinT;
//				initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//				grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//				SortProjection<T>(grid);
//
//				len = hypot((curIdx - cntImgX)*Img.m_Step.x - sour[0], (sliceIdx - cntImgY)*Img.m_Step.y - sour[1]);
//				coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//				coef = coef / (len * FanGeo.m_DetStp);
//				summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
//			}
//		}
//		dprj[angIdx* FanGeo.m_DetN + detIdx] = summ;
//	} //End case 1
//	else if (curAng > (PI * 0.25) && curAng <= (PI * 0.75))
//	{
//		summ = 0;
//		for (sliceIdx = 0; sliceIdx != Img.m_Reso.x; ++sliceIdx)
//		{
//			coord = (sliceIdx - cntImgX - 0.5)* Img.m_Step.x;
//			initDir[0] = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
//			initDir[1] = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
//
//			curIdx = int(initDir[0] / Img.m_Step.y + cntImgY - 1);
//			maxIdx = ceil(initDir[1] / Img.m_Step.y + cntImgY + 1);
//			if (curIdx < 0)
//			{
//				curIdx = 0;
//			}
//			if (curIdx > Img.m_Reso.y)
//			{
//				continue;
//				//minyIdx = Img.m_Reso.y;
//			}
//			if (maxIdx < 0)
//			{
//				continue;
//				//maxyIdx = 0;
//			}
//			if (maxIdx > Img.m_Reso.y)
//			{
//				maxIdx = Img.m_Reso.y;
//			}
//
//			for (; curIdx < maxIdx; ++curIdx)
//			{
//				//��grid
//				grid[0][0] = coord;
//				grid[0][1] = (curIdx - cntImgY) * Img.m_Step.y;
//				grid[1][0] = coord + Img.m_Step.x;
//				grid[1][1] = (curIdx - cntImgY) * Img.m_Step.y;
//				grid[2][0] = coord + Img.m_Step.x;
//				grid[2][1] = (curIdx - cntImgY + 1) * Img.m_Step.y;
//				grid[3][0] = coord;
//				grid[3][1] = (curIdx - cntImgY + 1) * Img.m_Step.y;
//
//				//�����ĸ����Ӧ��ǰ��det index
//				dir[0] = grid[0][0] - sour[0];
//				dir[1] = grid[0][1] - sour[1];
//				len = hypot(dir[0], dir[1]);
//				dir[0] /= len;
//				dir[1] /= len;
//				dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//				dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//				initDir[0] = dir[0] * cosT + dir[1] * sinT;
//				initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//				grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//				dir[0] = grid[1][0] - sour[0];
//				dir[1] = grid[1][1] - sour[1];
//				len = hypot(dir[0], dir[1]);
//				dir[0] /= len;
//				dir[1] /= len;
//				dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//				dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//				initDir[0] = dir[0] * cosT + dir[1] * sinT;
//				initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//				grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//				dir[0] = grid[2][0] - sour[0];
//				dir[1] = grid[2][1] - sour[1];
//				len = hypot(dir[0], dir[1]);
//				dir[0] /= len;
//				dir[1] /= len;
//				dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//				dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//				initDir[0] = dir[0] * cosT + dir[1] * sinT;
//				initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//				grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//				dir[0] = grid[3][0] - sour[0];
//				dir[1] = grid[3][1] - sour[1];
//				len = hypot(dir[0], dir[1]);
//				dir[0] /= len;
//				dir[1] /= len;
//				dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//				dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//				initDir[0] = dir[0] * cosT + dir[1] * sinT;
//				initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//				grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//				SortProjection<T>(grid);
//				len = hypot((sliceIdx - cntImgX)*Img.m_Step.x - sour[0], (curIdx - cntImgY)*Img.m_Step.y - sour[1]);
//				coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//				coef = coef / (len * FanGeo.m_DetStp);
//				summ += coef * dimg[curIdx * Img.m_Reso.x + sliceIdx];
//
//			} //End for one slice
//		}// End all slices
//		dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
//	}// End case 2
//	else if (curAng > (PI * 0.75) && curAng <= (PI * 1.25))
//	{
//		summ = 0;
//		for (sliceIdx = 0; sliceIdx != Img.m_Reso.y; ++sliceIdx)
//		{
//			coord = (sliceIdx - cntImgY + 0.5)* Img.m_Step.y;
//			//�������λ�ù����ཻ�����;
//			initDir[1] = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
//			initDir[0] = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
//
//			curIdx = int(initDir[0] / Img.m_Step.x + cntImgX - 1);
//			maxIdx = ceilf(initDir[1] / Img.m_Step.x + cntImgX + 1);
//			if (curIdx < 0)
//			{
//				curIdx = 0;
//			}
//			if (curIdx > Img.m_Reso.x)
//			{
//				continue;
//
//			}
//			if (maxIdx < 0)
//			{
//				continue;
//			}
//			if (maxIdx > Img.m_Reso.x)
//			{
//				maxIdx = Img.m_Reso.x;
//			}
//
//			for (; curIdx < maxIdx; ++curIdx)
//			{
//				//��grid
//				grid[0][0] = (curIdx - cntImgX) * Img.m_Step.x;
//				grid[0][1] = coord;
//				grid[1][0] = (curIdx - cntImgX + 1) * Img.m_Step.x;
//				grid[1][1] = coord;
//				grid[2][0] = (curIdx - cntImgX + 1) * Img.m_Step.x;
//				grid[2][1] = coord + Img.m_Step.y;
//				grid[3][0] = (curIdx - cntImgX)* Img.m_Step.x;
//				grid[3][1] = coord + Img.m_Step.y;
//
//				//�����ĸ����Ӧ��ǰ��det index
//				dir[0] = grid[0][0] - sour[0];
//				dir[1] = grid[0][1] - sour[1];
//				len = hypot(dir[0], dir[1]);
//				dir[0] /= len;
//				dir[1] /= len;
//				dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//				dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//				initDir[0] = dir[0] * cosT + dir[1] * sinT;
//				initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//				grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//				dir[0] = grid[1][0] - sour[0];
//				dir[1] = grid[1][1] - sour[1];
//				len = hypot(dir[0], dir[1]);
//				dir[0] /= len;
//				dir[1] /= len;
//				dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//				dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//				initDir[0] = dir[0] * cosT + dir[1] * sinT;
//				initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//				grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//				dir[0] = grid[2][0] - sour[0];
//				dir[1] = grid[2][1] - sour[1];
//				len = hypot(dir[0], dir[1]);
//				dir[0] /= len;
//				dir[1] /= len;
//				dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//				dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//				initDir[0] = dir[0] * cosT + dir[1] * sinT;
//				initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//				grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//				dir[0] = grid[3][0] - sour[0];
//				dir[1] = grid[3][1] - sour[1];
//				len = hypot(dir[0], dir[1]);
//				dir[0] /= len;
//				dir[1] /= len;
//				dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//				dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//				initDir[0] = dir[0] * cosT + dir[1] * sinT;
//				initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//				grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//				SortProjection<T>(grid);
//				len = hypot((curIdx - cntImgX)*Img.m_Step.x - sour[0], (sliceIdx - cntImgY)*Img.m_Step.y - sour[1]);
//				coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//				coef = coef / (len * FanGeo.m_DetStp);
//				summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
//
//			} //End for one slice
//		}// End all slices
//		dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
//	}
//	else
//	{
//		summ = 0;
//		for (sliceIdx = 0; sliceIdx != Img.m_Reso.x; ++sliceIdx)
//		{
//			coord = (sliceIdx - cntImgX + 0.5)* Img.m_Step.x;
//			//�������λ�ù����ཻ�����;
//			initDir[1] = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
//			initDir[0] = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
//
//			curIdx = int(initDir[0] / Img.m_Step.y + cntImgY - 1);
//			maxIdx = ceilf(initDir[1] / Img.m_Step.y + cntImgY + 1);
//			if (curIdx < 0)
//			{
//				curIdx = 0;
//			}
//			if (curIdx > Img.m_Reso.y)
//			{
//				continue;
//			}
//			if (maxIdx < 0)
//			{
//				continue;
//			}
//			if (maxIdx > Img.m_Reso.y)
//			{
//				maxIdx = Img.m_Reso.y;
//			}
//
//			for (; curIdx < maxIdx; ++curIdx)
//			{
//				//��grid
//				grid[0][0] = coord;
//				grid[0][1] = (curIdx - cntImgY) * Img.m_Step.y;
//				grid[1][0] = coord + Img.m_Step.x;
//				grid[1][1] = (curIdx - cntImgY) * Img.m_Step.y;
//				grid[2][0] = coord + Img.m_Step.x;
//				grid[2][1] = (curIdx - cntImgY + 1) * Img.m_Step.y;
//				grid[3][0] = coord;
//				grid[3][1] = (curIdx - cntImgY + 1) * Img.m_Step.y;
//
//				//�����ĸ����Ӧ��ǰ��det index
//				dir[0] = grid[0][0] - sour[0];
//				dir[1] = grid[0][1] - sour[1];
//				len = hypot(dir[0], dir[1]);
//				dir[0] /= len;
//				dir[1] /= len;
//				dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//				dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//				initDir[0] = dir[0] * cosT + dir[1] * sinT;
//				initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//				grid[0][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx; //��index
//
//				dir[0] = grid[1][0] - sour[0];
//				dir[1] = grid[1][1] - sour[1];
//				len = hypot(dir[0], dir[1]);
//				dir[0] /= len;
//				dir[1] /= len;
//				dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//				dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//				initDir[0] = dir[0] * cosT + dir[1] * sinT;
//				initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//				grid[1][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//				dir[0] = grid[2][0] - sour[0];
//				dir[1] = grid[2][1] - sour[1];
//				len = hypot(dir[0], dir[1]);
//				dir[0] /= len;
//				dir[1] /= len;
//				dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//				dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//				initDir[0] = dir[0] * cosT + dir[1] * sinT;
//				initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//				grid[2][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//				dir[0] = grid[3][0] - sour[0];
//				dir[1] = grid[3][1] - sour[1];
//				len = hypot(dir[0], dir[1]);
//				dir[0] /= len;
//				dir[1] /= len;
//				dir[0] = sour[0] + FanGeo.m_S2D * dir[0];
//				dir[1] = sour[1] + FanGeo.m_S2D * dir[1];
//				initDir[0] = dir[0] * cosT + dir[1] * sinT;
//				initDir[1] = dir[0] * (-sinT) + dir[1] * cosT;
//				grid[3][2] = atan(-initDir[0] / (initDir[1] - FanGeo.m_S2O)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//				SortProjection<T>(grid);
//				len = hypot((sliceIdx - cntImgX)*Img.m_Step.x - sour[0], (curIdx - cntImgY)*Img.m_Step.y - sour[1]);
//				coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//				coef = coef / (len * FanGeo.m_DetStp);
//				summ += coef * dimg[curIdx* Img.m_Reso.x + sliceIdx];
//				//summ += *coef;
//			} //End for one slice
//		}// End all slices
//		dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
//	}
//}
//
//template<typename T>
//void proj_AIM2_temp(T* dprj, T* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_proj_AIM2_ker<T> << < gid, blk >> >(dprj, dimg, FanGeo, Img);
//}
//
//void proj_AIM2(float* dprj, float* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	proj_AIM2_temp<float>(dprj, dimg, FanGeo, Img, blk, gid);
//}
//
//void proj_AIM2(double* dprj, double* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	proj_AIM2_temp<double>(dprj, dimg, FanGeo, Img, blk, gid);
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
//template<typename T>
//__global__ void _proj_AIM_ker(T* dprj, T* dimg, const FanEDGeo FanGeo, const Image Img)
//{
//	int detIdx = threadIdx.x + blockIdx.x * blockDim.x;
//	int angIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detIdx < FanGeo.m_DetN &&angIdx < FanGeo.m_ViwN)
//	{
//		const int extraIdx = 8;
//		const T cntImgX = (static_cast<T>(Img.m_Reso.x) - 1) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
//		const T cntImgY = (static_cast<T>(Img.m_Reso.y) - 1) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
//		const T area = Img.m_Step.x * Img.m_Step.y;
//		T curAng = FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp;
//		while (curAng < 0){ curAng += (TWOPI); }
//		while (curAng > TWOPI){ curAng -= (TWOPI); }
//		T cosT = cos(curAng);
//		T sinT = sin(curAng);
//		T sour[2] = { -FanGeo.m_S2O * sinT, FanGeo.m_S2O * cosT };
//		T SVA[3], SVB[3];
//		calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, detIdx);
//		T vv = acos(abs(SVA[0] * SVB[0] + SVA[1] * SVB[1]));
//		//vv = asin(sqrt(1.0 - vv * vv));
//
//		unsigned int sliceIdx;
//		T summ;
//		T coord;
//		T minC, maxC;
//		int maxIdx;
//		int curIdx;
//		T grid[4][3];
//
//		T pdist, coef;
//		T curDirAng = atan(((detIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp) / FanGeo.m_O2D) + curAng;
//		while (curDirAng < 0){ curDirAng += (TWOPI); }
//		while (curDirAng > TWOPI){ curDirAng -= (TWOPI); }
//		T initP[2]; // initial pixel position
//		T xPos;
//		if (curDirAng <= PI * 0.25 || curDirAng > PI * 1.75)
//		{
//			summ = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.y; ++sliceIdx)
//			{
//				coord = (static_cast<T>(sliceIdx)-cntImgY - 0.5) * Img.m_Step.y;
//				minC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
//				maxC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.x + cntImgX) - extraIdx;
//				maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + extraIdx;
//				if (curIdx > static_cast<int>(Img.m_Reso.x - 1) || maxIdx < 0)
//				{
//					continue;
//				}
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//				if (maxIdx > Img.m_Reso.x - 1)
//				{
//					maxIdx = Img.m_Reso.x - 1;
//				}
//
//				for (; curIdx <= maxIdx; curIdx++)
//				{
//					grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
//					grid[0][1] = coord;
//					grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[1][1] = coord;
//					grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[2][1] = coord + Img.m_Step.y;
//					grid[3][0] = (curIdx - cntImgX - 0.5)* Img.m_Step.x;
//					grid[3][1] = coord + Img.m_Step.y;
//
//
//					initP[0] = grid[0][0] * cosT + grid[0][1] * sinT;
//					initP[1] = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[0][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[1][0] * cosT + grid[1][1] * sinT;
//					initP[1] = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[1][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[2][0] * cosT + grid[2][1] * sinT;
//					initP[1] = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[2][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[3][0] * cosT + grid[3][1] * sinT;
//					initP[1] = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[3][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					SortProjection<T>(grid);
//
//					pdist = hypot((static_cast<T>(curIdx)-cntImgX)*static_cast<T>(Img.m_Step.x) - sour[0], (static_cast<T>(sliceIdx)-cntImgY)*static_cast<T>(Img.m_Step.y) - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//					coef = coef / (pdist* vv);
//					summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
//				}
//			}
//			dprj[angIdx* FanGeo.m_DetN + detIdx] = summ;
//			return;
//		} //End case 1
//		else if (curDirAng > PI * 0.25 && curDirAng <= PI * 0.75)
//		{
//			summ = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
//			{
//				coord = (static_cast<T>(sliceIdx)-cntImgX - 0.5)* Img.m_Step.x;
//				//�������λ�ù����ཻ�����;
//				minC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
//				maxC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.y + cntImgY) - extraIdx;
//				maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + extraIdx;
//
//				if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0)
//				{
//					continue;
//				}
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//
//				if (maxIdx > Img.m_Reso.y - 1)
//				{
//					maxIdx = Img.m_Reso.y - 1;
//				}
//
//				for (; curIdx <= maxIdx; ++curIdx)
//				{
//					//��grid
//					grid[0][0] = coord;
//					grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[1][0] = coord + Img.m_Step.x;
//					grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[2][0] = coord + Img.m_Step.x;
//					grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//					grid[3][0] = coord;
//					grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//
//
//					//�����ĸ����Ӧ��ǰ��det index
//					initP[0] = grid[0][0] * cosT + grid[0][1] * sinT;
//					initP[1] = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[0][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[1][0] * cosT + grid[1][1] * sinT;
//					initP[1] = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[1][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[2][0] * cosT + grid[2][1] * sinT;
//					initP[1] = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[2][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[3][0] * cosT + grid[3][1] * sinT;
//					initP[1] = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[3][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//
//					SortProjection<T>(grid);
//					pdist = hypot((static_cast<T>(sliceIdx)-cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(curIdx)-cntImgY)*Img.m_Step.y - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//					coef = coef / (pdist* vv);
//					summ += coef * dimg[curIdx * Img.m_Reso.x + sliceIdx];
//
//				} //End for one slice
//			}// End all slices
//			dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
//			return;
//		}// End case 2
//		else if (curDirAng > PI * 0.75 && curDirAng <= PI * 1.25)
//		{
//			summ = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.y; ++sliceIdx)
//			{
//				coord = (static_cast<T>(sliceIdx)-cntImgY + 0.5)* Img.m_Step.y;
//				//�������λ�ù����ཻ�����;
//				maxC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
//				minC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.x + cntImgX) - extraIdx;
//				maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + extraIdx;
//
//				if (curIdx > static_cast<int>(Img.m_Reso.x) || maxIdx < 0)
//				{
//					continue;
//				}
//
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//
//				if (maxIdx > Img.m_Reso.x - 1)
//				{
//					maxIdx = Img.m_Reso.x - 1;
//				}
//
//				for (; curIdx <= maxIdx; ++curIdx)
//				{
//					grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
//					grid[0][1] = coord - Img.m_Step.y;
//					grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[1][1] = coord - Img.m_Step.y;
//					grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[2][1] = coord;
//					grid[3][0] = (curIdx - cntImgX - 0.5)* Img.m_Step.x;
//					grid[3][1] = coord;
//
//
//					//�����ĸ����Ӧ��ǰ��det index
//					initP[0] = grid[0][0] * cosT + grid[0][1] * sinT;
//					initP[1] = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[0][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[1][0] * cosT + grid[1][1] * sinT;
//					initP[1] = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[1][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[2][0] * cosT + grid[2][1] * sinT;
//					initP[1] = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[2][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[3][0] * cosT + grid[3][1] * sinT;
//					initP[1] = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[3][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//
//					SortProjection<T>(grid);
//					pdist = hypot((static_cast<T>(curIdx)-cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(sliceIdx)-cntImgY)*Img.m_Step.y - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//					coef = coef / (pdist* vv);
//					summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
//
//				} //End for one slice
//			}// End all slices
//			dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
//			return;
//		}
//		else //if (curDirAng > 1.25 * PI &&curDirAng <= 1.75*PI)
//		{
//			summ = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
//			{
//				coord = (static_cast<T>(sliceIdx)-cntImgX + 0.5)* Img.m_Step.x;
//				//�������λ�ù����ཻ�����;
//				maxC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
//				minC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.y + cntImgY) - extraIdx;
//				maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + extraIdx;
//
//				if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0)
//				{
//					continue;
//				}
//
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//
//				if (maxIdx > Img.m_Reso.y - 1)
//				{
//					maxIdx = Img.m_Reso.y - 1;
//				}
//
//				for (; curIdx <= maxIdx; ++curIdx)
//				{
//					grid[0][0] = coord - Img.m_Step.x;
//					grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[1][0] = coord;
//					grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[2][0] = coord - Img.m_Step.x;
//					grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//					grid[3][0] = coord;
//					grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//
//					initP[0] = grid[0][0] * cosT + grid[0][1] * sinT;
//					initP[1] = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[0][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[1][0] * cosT + grid[1][1] * sinT;
//					initP[1] = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[1][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[2][0] * cosT + grid[2][1] * sinT;
//					initP[1] = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[2][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[3][0] * cosT + grid[3][1] * sinT;
//					initP[1] = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[3][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					SortProjection<T>(grid);
//					pdist = hypot((static_cast<T>(sliceIdx)-cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(curIdx)-cntImgY)*Img.m_Step.y - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//					coef = coef / (pdist* vv);
//					summ += coef * dimg[curIdx* Img.m_Reso.x + sliceIdx];
//				} //End for one slice
//			}// End all slices
//			dprj[angIdx * FanGeo.m_DetN + detIdx] = summ;
//			return;
//		}
//	}
//}
//template<typename T>
//void proj_AIM_temp(T* dprj, T* dimg, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid){ _proj_AIM_ker<T> << <gid, blk >> >(dprj, dimg, FanGeo, Img); }
//void proj_AIM(float* dprj, float* dimg, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid){ proj_AIM_temp<float>(dprj, dimg, FanGeo, Img, blk, gid); }
//void proj_AIM(double* dprj, double* dimg, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid){ proj_AIM_temp<double>(dprj, dimg, FanGeo, Img, blk, gid); }
//
//template<typename T>
//void proj_AIM_3GPU_temp(T* dprj, T* dimg, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	cuint imgReso = Img.m_Reso.x * Img.m_Reso.y;
//	cuint prjReso = FanGeo.m_ViwN * FanGeo.m_DetN;
//	cuint thrPrjReso = prjReso / 3;
//	cudaStream_t streams[3];
//
//	FanEDGeo nFanGeo(FanGeo);
//	nFanGeo.m_ViwN = nFanGeo.m_ViwN / 3;
//
//	T* dImg[3];
//	T* dPrj[3];
//	dim3 nblk = blk;
//	dim3 ngid(
//		(FanGeo.m_DetN + nblk.x - 1) / nblk.x,
//		(FanGeo.m_ViwN / 3 + nblk.y - 1) / nblk.y);
//
//	for (size_t i = 0; i < 3; i++)
//	{
//		checkCudaErrors(cudaSetDevice(i));
//		checkCudaErrors(cudaStreamCreate(&streams[i]));
//		checkCudaErrors(cudaMalloc((void**)&dImg[i], sizeof(T) * imgReso));
//		checkCudaErrors(cudaMalloc((void**)&dPrj[i], sizeof(T) * thrPrjReso));
//		checkCudaErrors(cudaMemcpy(dImg[i], dimg, sizeof(T) * imgReso, cudaMemcpyHostToDevice));
//		checkCudaErrors(cudaMemcpy(dPrj[i], dprj + i* thrPrjReso, sizeof(T) * thrPrjReso, cudaMemcpyDeviceToDevice));
//		proj_AIM(dPrj[i], dImg[i], nFanGeo, Img, nblk, ngid);
//	}
//
//
//}
//
//template<typename T>
//__global__ void _proj_AIM_ker(T* dprj, T* draw, T* dimg, const FanEDGeo FanGeo, const Image Img)
//{
//	int detIdx = threadIdx.x + blockIdx.x * blockDim.x;
//	int angIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detIdx < FanGeo.m_DetN &&angIdx < FanGeo.m_ViwN)
//	{
//		const int extraIdx = 8;
//		const T cntImgX = (static_cast<T>(Img.m_Reso.x) - 1) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
//		const T cntImgY = (static_cast<T>(Img.m_Reso.y) - 1) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
//		const T area = Img.m_Step.x * Img.m_Step.y;
//		T curAng = FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp;
//		while (curAng < 0){ curAng += (TWOPI); }
//		while (curAng > TWOPI){ curAng -= (TWOPI); }
//		T cosT = cos(curAng);
//		T sinT = sin(curAng);
//		T sour[2] = { -FanGeo.m_S2O * sinT, FanGeo.m_S2O * cosT };
//		T SVA[3], SVB[3];
//		calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, detIdx);
//		T vv = acos(abs(SVA[0] * SVB[0] + SVA[1] * SVB[1]));
//		//vv = asin(sqrt(1.0 - vv * vv));
//
//		unsigned int sliceIdx;
//		T summ;
//		T coord;
//		T minC, maxC;
//		int maxIdx;
//		int curIdx;
//		T grid[4][3];
//
//		T pdist, coef;
//		T curDirAng = atan(((detIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp) / FanGeo.m_O2D) + curAng;
//		while (curDirAng < 0){ curDirAng += (TWOPI); }
//		while (curDirAng > TWOPI){ curDirAng -= (TWOPI); }
//		T initP[2]; // initial pixel position
//		T xPos;
//		T weight;
//		if (curDirAng <= PI * 0.25 || curDirAng > PI * 1.75)
//		{
//			summ = 0;
//			weight = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.y; ++sliceIdx)
//			{
//				coord = (static_cast<T>(sliceIdx)-cntImgY - 0.5) * Img.m_Step.y;
//				minC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
//				maxC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.x + cntImgX) - extraIdx;
//				maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + extraIdx;
//				if (curIdx > static_cast<int>(Img.m_Reso.x - 1) || maxIdx < 0)
//				{
//					continue;
//				}
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//				if (maxIdx > Img.m_Reso.x - 1)
//				{
//					maxIdx = Img.m_Reso.x - 1;
//				}
//
//				for (; curIdx <= maxIdx; curIdx++)
//				{
//					grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
//					grid[0][1] = coord;
//					grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[1][1] = coord;
//					grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[2][1] = coord + Img.m_Step.y;
//					grid[3][0] = (curIdx - cntImgX - 0.5)* Img.m_Step.x;
//					grid[3][1] = coord + Img.m_Step.y;
//
//
//					initP[0] = grid[0][0] * cosT + grid[0][1] * sinT;
//					initP[1] = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[0][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[1][0] * cosT + grid[1][1] * sinT;
//					initP[1] = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[1][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[2][0] * cosT + grid[2][1] * sinT;
//					initP[1] = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[2][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[3][0] * cosT + grid[3][1] * sinT;
//					initP[1] = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[3][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					//����;
//					SortProjection<T>(grid);
//
//					pdist = hypot((static_cast<T>(curIdx)-cntImgX)*static_cast<T>(Img.m_Step.x) - sour[0], (static_cast<T>(sliceIdx)-cntImgY)*static_cast<T>(Img.m_Step.y) - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//					coef = coef / (pdist* vv);
//					weight += coef;
//					summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
//				}
//			}
//			if (!IS_ZERO(weight))
//			{
//				dprj[angIdx* FanGeo.m_DetN + detIdx] = (draw[angIdx * FanGeo.m_DetN + detIdx] - summ) / weight;
//			}
//			else
//			{
//				dprj[angIdx* FanGeo.m_DetN + detIdx] = 0;
//			}
//			return;
//		} //End case 1
//		else if (curDirAng > PI * 0.25 && curDirAng <= PI * 0.75)
//		{
//			summ = 0;
//			weight = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
//			{
//				coord = (static_cast<T>(sliceIdx)-cntImgX - 0.5)* Img.m_Step.x;
//				//�������λ�ù����ཻ�����;
//				minC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
//				maxC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.y + cntImgY) - extraIdx;
//				maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + extraIdx;
//
//				if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0)
//				{
//					continue;
//				}
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//
//				if (maxIdx > Img.m_Reso.y - 1)
//				{
//					maxIdx = Img.m_Reso.y - 1;
//				}
//
//				for (; curIdx <= maxIdx; ++curIdx)
//				{
//					//��grid
//					grid[0][0] = coord;
//					grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[1][0] = coord + Img.m_Step.x;
//					grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[2][0] = coord + Img.m_Step.x;
//					grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//					grid[3][0] = coord;
//					grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//
//
//					//�����ĸ����Ӧ��ǰ��det index
//					initP[0] = grid[0][0] * cosT + grid[0][1] * sinT;
//					initP[1] = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[0][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[1][0] * cosT + grid[1][1] * sinT;
//					initP[1] = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[1][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[2][0] * cosT + grid[2][1] * sinT;
//					initP[1] = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[2][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[3][0] * cosT + grid[3][1] * sinT;
//					initP[1] = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[3][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					//����;
//					SortProjection<T>(grid);
//					pdist = hypot((static_cast<T>(sliceIdx)-cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(curIdx)-cntImgY)*Img.m_Step.y - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//					coef = coef / (pdist* vv);
//					weight += coef;
//					summ += coef * dimg[curIdx * Img.m_Reso.x + sliceIdx];
//
//				} //End for one slice
//			}// End all slices
//			if (!IS_ZERO(weight))
//			{
//				dprj[angIdx* FanGeo.m_DetN + detIdx] = (draw[angIdx * FanGeo.m_DetN + detIdx] - summ) / weight;
//			}
//			else
//			{
//				dprj[angIdx* FanGeo.m_DetN + detIdx] = 0;
//			}
//			return;
//		}// End case 2
//		else if (curDirAng > PI * 0.75 && curDirAng <= PI * 1.25)
//		{
//			summ = 0;
//			weight = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.y; ++sliceIdx)
//			{
//				coord = (static_cast<T>(sliceIdx)-cntImgY + 0.5)* Img.m_Step.y;
//				//�������λ�ù����ཻ�����;
//				maxC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
//				minC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.x + cntImgX) - extraIdx;
//				maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + extraIdx;
//
//				if (curIdx > static_cast<int>(Img.m_Reso.x) || maxIdx < 0)
//				{
//					continue;
//				}
//
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//
//				if (maxIdx > Img.m_Reso.x - 1)
//				{
//					maxIdx = Img.m_Reso.x - 1;
//				}
//
//				for (; curIdx <= maxIdx; ++curIdx)
//				{
//					grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
//					grid[0][1] = coord - Img.m_Step.y;
//					grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[1][1] = coord - Img.m_Step.y;
//					grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[2][1] = coord;
//					grid[3][0] = (curIdx - cntImgX - 0.5)* Img.m_Step.x;
//					grid[3][1] = coord;
//
//
//					//�����ĸ����Ӧ��ǰ��det index
//					initP[0] = grid[0][0] * cosT + grid[0][1] * sinT;
//					initP[1] = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[0][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[1][0] * cosT + grid[1][1] * sinT;
//					initP[1] = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[1][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[2][0] * cosT + grid[2][1] * sinT;
//					initP[1] = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[2][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[3][0] * cosT + grid[3][1] * sinT;
//					initP[1] = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[3][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					//����;
//					SortProjection<T>(grid);
//					pdist = hypot((static_cast<T>(curIdx)-cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(sliceIdx)-cntImgY)*Img.m_Step.y - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//					coef = coef / (pdist* vv);
//					weight += coef;
//					summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
//
//				} //End for one slice
//			}// End all slices
//			if (!IS_ZERO(weight))
//			{
//				dprj[angIdx* FanGeo.m_DetN + detIdx] = (draw[angIdx * FanGeo.m_DetN + detIdx] - summ) / weight;
//			}
//			else
//			{
//				dprj[angIdx* FanGeo.m_DetN + detIdx] = 0;
//			}
//			return;
//		}
//		else //if (curDirAng > 1.25 * PI &&curDirAng <= 1.75*PI)
//		{
//			summ = 0;
//			weight = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
//			{
//				coord = (static_cast<T>(sliceIdx)-cntImgX + 0.5)* Img.m_Step.x;
//				//�������λ�ù����ཻ�����;
//				maxC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
//				minC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.y + cntImgY) - extraIdx;
//				maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + extraIdx;
//
//				if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0)
//				{
//					continue;
//				}
//
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//
//				if (maxIdx > Img.m_Reso.y - 1)
//				{
//					maxIdx = Img.m_Reso.y - 1;
//				}
//
//				for (; curIdx <= maxIdx; ++curIdx)
//				{
//					grid[0][0] = coord - Img.m_Step.x;
//					grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[1][0] = coord;
//					grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[2][0] = coord - Img.m_Step.x;
//					grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//					grid[3][0] = coord;
//					grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//
//					initP[0] = grid[0][0] * cosT + grid[0][1] * sinT;
//					initP[1] = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[0][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[1][0] * cosT + grid[1][1] * sinT;
//					initP[1] = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[1][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[2][0] * cosT + grid[2][1] * sinT;
//					initP[1] = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[2][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[3][0] * cosT + grid[3][1] * sinT;
//					initP[1] = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[3][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					SortProjection<T>(grid);
//					pdist = hypot((static_cast<T>(sliceIdx)-cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(curIdx)-cntImgY)*Img.m_Step.y - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//					coef = coef / (pdist* vv);
//					weight += coef;
//					summ += coef * dimg[curIdx* Img.m_Reso.x + sliceIdx];
//				} //End for one slice
//			}// End all slices
//			if (!IS_ZERO(weight))
//			{
//				dprj[angIdx* FanGeo.m_DetN + detIdx] = (draw[angIdx * FanGeo.m_DetN + detIdx] - summ) / weight;
//			}
//			else
//			{
//				dprj[angIdx* FanGeo.m_DetN + detIdx] = 0;
//			}
//			return;
//		}
//	}
//}
//
//template<typename T>
//void proj_AIM_temp(T* dprj, T* draw, T* dimg, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid){ _proj_AIM_ker<T> << <gid, blk >> >(dprj, draw, dimg, FanGeo, Img); }
//
//void proj_AIM(float* dprj, float* draw, float* dimg, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	proj_AIM_temp<float>(dprj, draw, dimg, FanGeo, Img, blk, gid);
//}
//void proj_AIM(double* dprj, double* draw, double* dimg, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	proj_AIM_temp<double>(dprj, draw, dimg, FanGeo, Img, blk, gid);
//}
//
//
//
//
//
//
//
//
//template<typename T>
//__global__ void _proj_AIM_ker(T* dprj, T* draw, T* dimg, const FanEDGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	int detIdx = threadIdx.x + blockIdx.x * blockDim.x;
//	int prjIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detIdx < FanGeo.m_DetN && prjIdx < numPerSubSet)
//	{
//		int angIdx = prjIdx * subSetNum + curSubSetIdx;
//		const int extraIdx = 8;
//		const T cntImgX = (static_cast<T>(Img.m_Reso.x) - 1) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
//		const T cntImgY = (static_cast<T>(Img.m_Reso.y) - 1) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
//		const T area = Img.m_Step.x * Img.m_Step.y;
//		T curAng = FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp;
//		while (curAng < 0){ curAng += (TWOPI); }
//		while (curAng > TWOPI){ curAng -= (TWOPI); }
//		T cosT = cos(curAng);
//		T sinT = sin(curAng);
//		T sour[2] = { -FanGeo.m_S2O * sinT, FanGeo.m_S2O * cosT };
//		T SVA[3], SVB[3];
//		calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, detIdx);
//		T vv = acos(abs(SVA[0] * SVB[0] + SVA[1] * SVB[1]));
//		//vv = asin(sqrt(1.0 - vv * vv));
//
//		unsigned int sliceIdx;
//		T summ;
//		T coord;
//		T minC, maxC;
//		int maxIdx;
//		int curIdx;
//		T grid[4][3];
//
//		T pdist, coef;
//		T curDirAng = atan(((detIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp) / FanGeo.m_O2D) + curAng;
//		while (curDirAng < 0){ curDirAng += (TWOPI); }
//		while (curDirAng > TWOPI){ curDirAng -= (TWOPI); }
//		T initP[2]; // initial pixel position
//		T xPos;
//		T weight;
//		if (curDirAng <= PI * 0.25 || curDirAng > PI * 1.75)
//		{
//			summ = 0;
//			weight = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.y; ++sliceIdx)
//			{
//				coord = (static_cast<T>(sliceIdx)-cntImgY - 0.5) * Img.m_Step.y;
//				minC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
//				maxC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.x + cntImgX) - extraIdx;
//				maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + extraIdx;
//				if (curIdx > static_cast<int>(Img.m_Reso.x - 1) || maxIdx < 0)
//				{
//					continue;
//				}
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//				if (maxIdx > Img.m_Reso.x - 1)
//				{
//					maxIdx = Img.m_Reso.x - 1;
//				}
//
//				for (; curIdx <= maxIdx; curIdx++)
//				{
//					grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
//					grid[0][1] = coord;
//					grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[1][1] = coord;
//					grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[2][1] = coord + Img.m_Step.y;
//					grid[3][0] = (curIdx - cntImgX - 0.5)* Img.m_Step.x;
//					grid[3][1] = coord + Img.m_Step.y;
//
//
//					initP[0] = grid[0][0] * cosT + grid[0][1] * sinT;
//					initP[1] = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[0][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[1][0] * cosT + grid[1][1] * sinT;
//					initP[1] = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[1][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[2][0] * cosT + grid[2][1] * sinT;
//					initP[1] = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[2][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[3][0] * cosT + grid[3][1] * sinT;
//					initP[1] = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[3][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					//����;
//					SortProjection<T>(grid);
//
//					pdist = hypot((static_cast<T>(curIdx)-cntImgX)*static_cast<T>(Img.m_Step.x) - sour[0], (static_cast<T>(sliceIdx)-cntImgY)*static_cast<T>(Img.m_Step.y) - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//					coef = coef / (pdist* vv);
//					weight += coef;
//					summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
//				}
//			}
//			if (!IS_ZERO(weight))
//			{
//				dprj[prjIdx* FanGeo.m_DetN + detIdx] = (draw[angIdx * FanGeo.m_DetN + detIdx] - summ) / weight;
//			}
//			else
//			{
//				dprj[prjIdx* FanGeo.m_DetN + detIdx] = 0;
//			}
//			return;
//		} //End case 1
//		else if (curDirAng > PI * 0.25 && curDirAng <= PI * 0.75)
//		{
//			summ = 0;
//			weight = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
//			{
//				coord = (static_cast<T>(sliceIdx)-cntImgX - 0.5)* Img.m_Step.x;
//				//�������λ�ù����ཻ�����;
//				minC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
//				maxC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.y + cntImgY) - extraIdx;
//				maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + extraIdx;
//
//				if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0)
//				{
//					continue;
//				}
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//
//				if (maxIdx > Img.m_Reso.y - 1)
//				{
//					maxIdx = Img.m_Reso.y - 1;
//				}
//
//				for (; curIdx <= maxIdx; ++curIdx)
//				{
//					//��grid
//					grid[0][0] = coord;
//					grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[1][0] = coord + Img.m_Step.x;
//					grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[2][0] = coord + Img.m_Step.x;
//					grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//					grid[3][0] = coord;
//					grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//
//
//					//�����ĸ����Ӧ��ǰ��det index
//					initP[0] = grid[0][0] * cosT + grid[0][1] * sinT;
//					initP[1] = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[0][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[1][0] * cosT + grid[1][1] * sinT;
//					initP[1] = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[1][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[2][0] * cosT + grid[2][1] * sinT;
//					initP[1] = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[2][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[3][0] * cosT + grid[3][1] * sinT;
//					initP[1] = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[3][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					//����;
//					SortProjection<T>(grid);
//					pdist = hypot((static_cast<T>(sliceIdx)-cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(curIdx)-cntImgY)*Img.m_Step.y - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//					coef = coef / (pdist* vv);
//					weight += coef;
//					summ += coef * dimg[curIdx * Img.m_Reso.x + sliceIdx];
//
//				} //End for one slice
//			}// End all slices
//			if (!IS_ZERO(weight))
//			{
//				dprj[prjIdx* FanGeo.m_DetN + detIdx] = (draw[angIdx * FanGeo.m_DetN + detIdx] - summ) / weight;
//			}
//			else
//			{
//				dprj[prjIdx* FanGeo.m_DetN + detIdx] = 0;
//			}
//			return;
//		}// End case 2
//		else if (curDirAng > PI * 0.75 && curDirAng <= PI * 1.25)
//		{
//			summ = 0;
//			weight = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.y; ++sliceIdx)
//			{
//				coord = (static_cast<T>(sliceIdx)-cntImgY + 0.5)* Img.m_Step.y;
//				//�������λ�ù����ཻ�����;
//				maxC = sour[0] + SVA[0] * (coord - sour[1]) / SVA[1];
//				minC = sour[0] + SVB[0] * (coord - sour[1]) / SVB[1];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.x + cntImgX) - extraIdx;
//				maxIdx = ceil(maxC / Img.m_Step.x + cntImgX) + extraIdx;
//
//				if (curIdx > static_cast<int>(Img.m_Reso.x) || maxIdx < 0)
//				{
//					continue;
//				}
//
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//
//				if (maxIdx > Img.m_Reso.x - 1)
//				{
//					maxIdx = Img.m_Reso.x - 1;
//				}
//
//				for (; curIdx <= maxIdx; ++curIdx)
//				{
//					grid[0][0] = (curIdx - cntImgX - 0.5) * Img.m_Step.x;
//					grid[0][1] = coord - Img.m_Step.y;
//					grid[1][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[1][1] = coord - Img.m_Step.y;
//					grid[2][0] = (curIdx - cntImgX + 0.5) * Img.m_Step.x;
//					grid[2][1] = coord;
//					grid[3][0] = (curIdx - cntImgX - 0.5)* Img.m_Step.x;
//					grid[3][1] = coord;
//
//
//					//�����ĸ����Ӧ��ǰ��det index
//					initP[0] = grid[0][0] * cosT + grid[0][1] * sinT;
//					initP[1] = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[0][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[1][0] * cosT + grid[1][1] * sinT;
//					initP[1] = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[1][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[2][0] * cosT + grid[2][1] * sinT;
//					initP[1] = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[2][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[3][0] * cosT + grid[3][1] * sinT;
//					initP[1] = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[3][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//
//					//����;
//					SortProjection<T>(grid);
//					pdist = hypot((static_cast<T>(curIdx)-cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(sliceIdx)-cntImgY)*Img.m_Step.y - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//					coef = coef / (pdist* vv);
//					weight += coef;
//					summ += coef * dimg[sliceIdx* Img.m_Reso.x + curIdx];
//
//				} //End for one slice
//			}// End all slices
//			if (!IS_ZERO(weight))
//			{
//				dprj[prjIdx* FanGeo.m_DetN + detIdx] = (draw[angIdx * FanGeo.m_DetN + detIdx] - summ) / weight;
//			}
//			else
//			{
//				dprj[prjIdx* FanGeo.m_DetN + detIdx] = 0;
//			}
//			return;
//		}
//		else //if (curDirAng > 1.25 * PI &&curDirAng <= 1.75*PI)
//		{
//			summ = 0;
//			weight = 0;
//			for (sliceIdx = 0; sliceIdx < Img.m_Reso.x; ++sliceIdx)
//			{
//				coord = (static_cast<T>(sliceIdx)-cntImgX + 0.5)* Img.m_Step.x;
//				//�������λ�ù����ཻ�����;
//				maxC = sour[1] + SVA[1] * (coord - sour[0]) / SVA[0];
//				minC = sour[1] + SVB[1] * (coord - sour[0]) / SVB[0];
//				if (maxC < minC)
//				{
//					pdist = minC;
//					minC = maxC;
//					maxC = pdist;
//				}
//				curIdx = int(minC / Img.m_Step.y + cntImgY) - extraIdx;
//				maxIdx = ceil(maxC / Img.m_Step.y + cntImgY) + extraIdx;
//
//				if (curIdx > static_cast<int>(Img.m_Reso.y) || maxIdx < 0)
//				{
//					continue;
//				}
//
//				if (curIdx < 0)
//				{
//					curIdx = 0;
//				}
//
//				if (maxIdx > Img.m_Reso.y - 1)
//				{
//					maxIdx = Img.m_Reso.y - 1;
//				}
//
//				for (; curIdx <= maxIdx; ++curIdx)
//				{
//					grid[0][0] = coord - Img.m_Step.x;
//					grid[0][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[1][0] = coord;
//					grid[1][1] = (curIdx - cntImgY - 0.5) * Img.m_Step.y;
//					grid[2][0] = coord - Img.m_Step.x;
//					grid[2][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//					grid[3][0] = coord;
//					grid[3][1] = (curIdx - cntImgY + 0.5) * Img.m_Step.y;
//
//					initP[0] = grid[0][0] * cosT + grid[0][1] * sinT;
//					initP[1] = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[0][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[1][0] * cosT + grid[1][1] * sinT;
//					initP[1] = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[1][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[2][0] * cosT + grid[2][1] * sinT;
//					initP[1] = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[2][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					initP[0] = grid[3][0] * cosT + grid[3][1] * sinT;
//					initP[1] = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
//					xPos = -initP[0] * FanGeo.m_S2D / initP[1];
//					grid[3][2] = xPos / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//					SortProjection<T>(grid);
//					pdist = hypot((static_cast<T>(sliceIdx)-cntImgX)*Img.m_Step.x - sour[0], (static_cast<T>(curIdx)-cntImgY)*Img.m_Step.y - sour[1]);
//					coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//					coef = coef / (pdist* vv);
//					weight += coef;
//					summ += coef * dimg[curIdx* Img.m_Reso.x + sliceIdx];
//				} //End for one slice
//			}// End all slices
//			if (!IS_ZERO(weight))
//			{
//				dprj[prjIdx* FanGeo.m_DetN + detIdx] = (draw[angIdx * FanGeo.m_DetN + detIdx] - summ) / weight;
//			}
//			else
//			{
//				dprj[prjIdx* FanGeo.m_DetN + detIdx] = 0;
//			}
//			return;
//		}
//	}
//}
//
//
//template<typename T>
//void proj_AIM_temp(T* dprj, T* draw, T* dimg, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_proj_AIM_ker<T> << < gid, blk >> >(dprj, draw, dimg, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx, blk, gid);
//}
//void proj_AIM(float* dprj, float* draw, float* dimg, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	proj_AIM_temp<float>(dprj, draw, dimg, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx, blk, gid);
//}
//void proj_AIM(double* dprj, double* draw, double* dimg, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	proj_AIM_temp<double>(dprj, draw, dimg, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx, blk, gid);
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
////////////////////////////////////////////////////////////////////////////
//// Back projection processes
////////////////////////////////////////////////////////////////////////////
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dimg, cuint angIdx, const bool angDir, const bool FOV, const FanEAGeo FanGeo, const Image Img)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	if (idX < Img.m_Reso.x && idY < Img.m_Reso.y)
//	{
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//		//cur rotation angle;
//		float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//		//current source position, assume the initial position is on the positive Y axis;
//		Ray2D ray;
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//		float2 minBox, maxBox, curImg, curDet, initDet, sincosAng;
//		float biasAng, tnear, tfar, weg, detID;
//		int flrID, ceilID;
//
//
//		minBox = make_float2(
//			MINO.x + idX * Img.m_Step.x,
//			MINO.y + idY * Img.m_Step.y);
//		maxBox = minBox + make_float2(Img.m_Step.x, Img.m_Step.y);
//		curImg = minBox + 0.5f * make_float2(Img.m_Step.x, Img.m_Step.y);
//		ray.d = normalize(curImg - ray.o);
//
//
//		//current detector element Cartesian coordinates;
//		curDet = ray.o + ray.d * FanGeo.m_S2D;
//		//original detector position;
//		initDet = rotation(curDet, cosT, -sinT);
//
//		//Bias angle;
//
//		sincosAng.x = initDet.x / FanGeo.m_S2D; //sin(ang)
//		sincosAng.y = (initDet.y - FanGeo.m_S2O) / (-FanGeo.m_S2D); //cos(ang)
//		biasAng = atan(sincosAng.x / sincosAng.y);
//
//		//judging the postive or negative angle direction;
//		if (angDir)
//		{
//			detID = biasAng / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//		}
//		else
//		{
//			detID = FanGeo.m_DetN - 1 - (biasAng / FanGeo.m_DetStp + FanGeo.m_DetCntIdx);
//		}
//		if (detID < 0)
//		{
//			dimg[idY * Img.m_Reso.x + idX] += 0;
//			return;
//		}
//		if (detID > FanGeo.m_DetN)
//		{
//			dimg[idY * Img.m_Reso.x + idX] += 0;
//			return;
//		}
//		tnear = 0;
//		tfar = 0;
//		flrID = detID;
//		ceilID = ceilf(detID);
//		intersectBox(ray, minBox, maxBox, &tnear, &tfar);
//
//		weg = (tfar - tnear); //no weighting;
//		if (!IS_ZERO(flrID - ceilID))
//		{
//			dimg[idY * Img.m_Reso.x + idX] += ((donp[flrID] * (ceilID - detID) + donp[ceilID] * (detID - flrID)) * weg);
//		}
//		else
//		{
//			dimg[idY * Img.m_Reso.x + idX] += donp[flrID] * weg;
//		}
//	}
//}
////////////////////////////////////////////////////////////////////////////
//// Backprojection from one angle projection data
//// donp: projection data from one angle;
//// dimg: image;
//// angIdx: angle index;
//// angDir: positive or negative direction;
//// FOV: do we only calculate inside the FOV;
//// FanGeo: Fan Beam Geometry
//// Img: Image configuration
//// blk: 32,32
//// gid: configured
////////////////////////////////////////////////////////////////////////////
//void bakproj_PIXEL(float* donp, float* dimg, cuint angIdx, const bool angDir, const bool FOV, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dimg, angIdx, angDir, FOV, FanGeo, Img);
//}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dimg, cuint angIdx, const bool angDir, const bool FOV, const FanEDGeo FanGeo, const Image Img)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	if (idX < Img.m_Reso.x && idY < Img.m_Reso.y)
//	{
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//		//cur rotation angle;
//		float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//current pixel position
//		float2 curPix = rotation(make_float2(MINO.x + Img.m_Step.x * idX, MINO.y + Img.m_Step.y * idY), cosT, -sinT); //��ת���λ��;
//		//�����Ӧ��detidx;
//
//		float2 minBox, maxBox, curImg;
//		float tnear(0), tfar(0), weg(0), detID(0);
//		int flrID(0), ceilID(0);
//
//		//���ﲻ���� angDir; ���Ǵ�С�����;
//		detID = FanGeo.m_S2D * curPix.x / (curPix.y - FanGeo.m_S2O) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//		if (detID < 0 || detID > FanGeo.m_DetN - 1)
//		{
//			dimg[idY * Img.m_Reso.x + idX] += 0;
//			return;
//		}
//		flrID = detID;
//		ceilID = ceilf(detID);
//
//		Ray2D ray;
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//		minBox = make_float2(
//			MINO.x + idX * Img.m_Step.x,
//			MINO.y + idY * Img.m_Step.y);
//		maxBox = minBox + make_float2(Img.m_Step.x, Img.m_Step.y);
//		curImg = minBox + 0.5f * make_float2(Img.m_Step.x, Img.m_Step.y);
//		ray.d = normalize(curImg - ray.o);
//		intersectBox(ray, minBox, maxBox, &tnear, &tfar); //�ཻ����;
//
//		weg = (tfar - tnear); //no weighting;
//		if (!IS_ZERO(flrID - ceilID))
//		{
//			dimg[idY * Img.m_Reso.x + idX] += ((donp[flrID] * (ceilID - detID) + donp[ceilID] * (detID - flrID)) * weg);
//		}
//		else
//		{
//			dimg[idY * Img.m_Reso.x + idX] += donp[flrID] * weg;
//		}
//	}
//}
//void bakproj_PIXEL(float* donp, float* dimg, cuint angIdx, const bool angDir, const bool FOV, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dimg, angIdx, angDir, FOV, FanGeo, Img);
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
//
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dvol, cuint angIdx, const bool angDir, const bool FOV, const ConeEAGeo ConeGeo, const Volume Vol)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	cuint idZ = threadIdx.z + blockDim.z * blockIdx.z;
//
//	if (idX < Vol.m_Reso.x && idY < Vol.m_Reso.y && idZ < Vol.m_Reso.z)
//	{
//		float3 MINO = make_float3(
//			-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//			-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//			-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//		//cur rotation angle;
//		float curAng = ConeGeo.m_ViwBeg + ConeGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//		//current source position, assume the initial position is on the positive Y axis;
//		Ray ray;
//		ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
//		float3 minBox, maxBox, curImg, curDet, initDet, sincosAng;
//		float biasAng, tnear, tfar, weg;
//		float2 detID;
//
//		minBox = make_float3(
//			MINO.x + idX * Vol.m_Step.x,
//			MINO.y + idY * Vol.m_Step.y,
//			MINO.z + idZ * Vol.m_Step.z);
//
//		maxBox = minBox + make_float3(Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z);
//		curImg = minBox + 0.5f * make_float3(Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z);
//		ray.d = normalize(curImg - ray.o);
//
//		//current detector element Cartesian coordinates;
//		curDet = ray.o + ray.d * ConeGeo.m_S2D;
//		//original detector position;
//		initDet = rotation(curDet, cosT, -sinT);  //�����λ��, ���ڿ�ʼ��Ӧ��index;
//
//		//Bias angle;
//		sincosAng.x = initDet.x / ConeGeo.m_S2D; //sin(ang)
//		sincosAng.y = (initDet.y - ConeGeo.m_S2O) / (-ConeGeo.m_S2D); //cos(ang)
//		biasAng = atan(sincosAng.x / sincosAng.y);
//
//		//judging the postive or negative angle direction;
//		if (angDir)
//		{
//			detID.x = biasAng / ConeGeo.m_DetStp + ConeGeo.m_DetCntIdx.x;
//		}
//		else
//		{
//			detID.x = ConeGeo.m_DetN - 1 - (biasAng / ConeGeo.m_DetStp + ConeGeo.m_DetCntIdx.x);
//		}
//
//		detID.y = (initDet.z / ConeGeo.m_DetHStp + ConeGeo.m_DetCntIdx.y);
//
//		if (detID.x < 0 || detID.x >(ConeGeo.m_DetN - 1) || detID.y < 0 || detID.y >(ConeGeo.m_DetHN - 1))
//		{
//			dvol[(idZ * Vol.m_Reso.z + idY) * Vol.m_Reso.x + idX] += 0;
//			return;
//		}
//
//		int minX = detID.x;
//		int maxX = ceil(detID.x);
//		int minZ = detID.y;
//		int maxZ = ceil(detID.y);
//
//		float Qmm = donp[minZ * ConeGeo.m_DetN + minX];
//		float Qmb = donp[maxZ * ConeGeo.m_DetN + minX];
//		float Qbm = donp[minZ * ConeGeo.m_DetN + maxX];
//		float Qbb = donp[maxZ * ConeGeo.m_DetN + maxX];
//
//		float val = _bilinearInterpolation_ker<float>(detID.x, detID.y, minX, maxX, minZ, maxZ, Qmm, Qmb, Qbm, Qbb);
//
//		tnear = 0;
//		tfar = 0;
//		intersectBox(ray, minBox, maxBox, &tnear, &tfar);
//		weg = (tfar - tnear); //no weighting;
//
//		dvol[(idZ * Vol.m_Reso.y + idY) * Vol.m_Reso.x + idX] += (weg * val);
//	}
//}
//void bakproj_PIXEL(float* donp, float* dvol, cuint angIdx, const bool angDir, const bool FOV, const ConeEAGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dvol, angIdx, angDir, FOV, ConeGeo, Vol);
//}
//
//
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dvol, float* dmsk, cuint angIdx, const bool angDir, const ConeEAGeo ConeGeo, const Volume Vol)//(float* donp, float* dvol, cuint angIdx, const bool angDir, const bool FOV, const ConeEAGeo ConeGeo, const Volume Vol)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	cuint idZ = threadIdx.z + blockDim.z * blockIdx.z;
//
//	if (idX < Vol.m_Reso.x && idY < Vol.m_Reso.y && idZ < Vol.m_Reso.z)
//	{
//		float3 MINO = make_float3(
//			-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//			-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//			-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//		//cur rotation angle;
//		float curAng = ConeGeo.m_ViwBeg + ConeGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//		//current source position, assume the initial position is on the positive Y axis;
//		Ray ray;
//		ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
//		float3 minBox, maxBox, curImg, curDet, initDet, sincosAng;
//		float biasAng, tnear, tfar, weg;
//		float2 detID;
//
//		minBox = make_float3(
//			MINO.x + idX * Vol.m_Step.x,
//			MINO.y + idY * Vol.m_Step.y,
//			MINO.z + idZ * Vol.m_Step.z);
//
//		maxBox = minBox + make_float3(Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z);
//		curImg = minBox + 0.5f * make_float3(Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z);
//		ray.d = normalize(curImg - ray.o);
//
//		//current detector element Cartesian coordinates;
//		curDet = ray.o + ray.d * ConeGeo.m_S2D;
//		//original detector position;
//		initDet = rotation(curDet, cosT, -sinT);  //�����λ��, ���ڿ�ʼ��Ӧ��index;
//
//		//Bias angle;
//		sincosAng.x = initDet.x / ConeGeo.m_S2D; //sin(ang)
//		sincosAng.y = (initDet.y - ConeGeo.m_S2O) / (-ConeGeo.m_S2D); //cos(ang)
//		biasAng = atan(sincosAng.x / sincosAng.y);
//
//		//judging the postive or negative angle direction;
//		if (angDir)
//		{
//			detID.x = biasAng / ConeGeo.m_DetStp + ConeGeo.m_DetCntIdx.x;
//		}
//		else
//		{
//			detID.x = ConeGeo.m_DetN - 1 - (biasAng / ConeGeo.m_DetStp + ConeGeo.m_DetCntIdx.x);
//		}
//
//		detID.y = (initDet.z / ConeGeo.m_DetHStp + ConeGeo.m_DetCntIdx.y);
//
//		if (detID.x < 0 || detID.x >(ConeGeo.m_DetN - 1) || detID.y < 0 || detID.y >(ConeGeo.m_DetHN - 1))
//		{
//			dvol[(idZ * Vol.m_Reso.z + idY) * Vol.m_Reso.x + idX] += 0;
//			return;
//		}
//
//		int minX = detID.x;
//		int maxX = ceil(detID.x);
//		int minZ = detID.y;
//		int maxZ = ceil(detID.y);
//
//		float Qmm = donp[minZ * ConeGeo.m_DetN + minX];
//		float Qmb = donp[maxZ * ConeGeo.m_DetN + minX];
//		float Qbm = donp[minZ * ConeGeo.m_DetN + maxX];
//		float Qbb = donp[maxZ * ConeGeo.m_DetN + maxX];
//
//		float val = _bilinearInterpolation_ker<float>(detID.x, detID.y, minX, maxX, minZ, maxZ, Qmm, Qmb, Qbm, Qbb);
//
//		tnear = 0;
//		tfar = 0;
//		intersectBox(ray, minBox, maxBox, &tnear, &tfar);
//		weg = (tfar - tnear); //no weighting;
//
//		dvol[(idZ * Vol.m_Reso.y + idY) * Vol.m_Reso.x + idX] += (weg * val * dmsk[idY * Vol.m_Reso.x + idX]);
//	}
//}
//void bakproj_PIXEL(float* donp, float* dvol, float* dmsk, cuint angIdx, const bool angDir, const ConeEAGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dvol, dmsk, angIdx, angDir, ConeGeo, Vol);
//}
//
//
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dimg, const bool angDir, const bool FOV, const FanEAGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	if (idX < Img.m_Reso.x && idY < Img.m_Reso.y)
//	{
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//		//cur rotation angle;
//		unsigned int angIdx = 0;
//		float curAng, cosT, sinT;
//		Ray2D ray;
//		float2 minBox, maxBox, curImg, curDet, initDet, sincosAng;
//		float biasAng, tnear, tfar, weg, detID;
//		int flrID, ceilID;
//		float summ(0);
//		for (unsigned int prjIdx = 0; prjIdx != numPerSubSet; ++prjIdx)
//		{
//			angIdx = prjIdx * subSetNum + curSubSetIdx;
//			curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//			//current source position, assume the initial position is on the positive Y axis;
//			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//			minBox = make_float2(
//				MINO.x + idX * Img.m_Step.x,
//				MINO.y + idY * Img.m_Step.y);
//			maxBox = minBox + make_float2(Img.m_Step.x, Img.m_Step.y);
//			curImg = minBox + 0.5f * make_float2(Img.m_Step.x, Img.m_Step.y);
//			ray.d = normalize(curImg - ray.o);
//
//
//			//current detector element Cartesian coordinates;
//			curDet = ray.o + ray.d * FanGeo.m_S2D;
//			//original detector position;
//			initDet = rotation(curDet, cosT, -sinT);
//
//			//Bias angle;
//
//			sincosAng.x = initDet.x / FanGeo.m_S2D; //sin(ang)
//			sincosAng.y = (initDet.y - FanGeo.m_S2O) / (-FanGeo.m_S2D); //cos(ang)
//			biasAng = atan(sincosAng.x / sincosAng.y);
//
//			//judging the postive or negative angle direction;
//			if (angDir)
//			{
//				detID = biasAng / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			}
//			else
//			{
//				detID = FanGeo.m_DetN - 1 - (biasAng / FanGeo.m_DetStp + FanGeo.m_DetCntIdx);
//			}
//			if (detID < 0)
//			{
//				dimg[idY * Img.m_Reso.x + idX] += 0;
//				continue;
//			}
//			if (detID > FanGeo.m_DetN)
//			{
//				dimg[idY * Img.m_Reso.x + idX] += 0;
//				continue;
//			}
//			tnear = 0;
//			tfar = 0;
//			flrID = detID;
//			ceilID = ceilf(detID);
//			intersectBox(ray, minBox, maxBox, &tnear, &tfar);
//
//			weg = (tfar - tnear); //no weighting;
//			if (!IS_ZERO(flrID - ceilID))
//			{
//				summ += ((donp[prjIdx * FanGeo.m_DetN + flrID] * (ceilID - detID) + donp[prjIdx * FanGeo.m_DetN + ceilID] * (detID - flrID)) * weg);
//			}
//			else
//			{
//				summ += donp[prjIdx * FanGeo.m_DetN + flrID] * weg;
//			}
//		}
//		dimg[idY * Img.m_Reso.x + idX] = summ;
//	}
//}
////////////////////////////////////////////////////////////////////////////
//// Back projection inside the subset with given subset index
//// donp: the projection data inside the subset 
//// dimg: image 
//// angDir: we use positive or negative index direction
//// FOV: do we only calculate inside the FOV
//// FanGeo: Fan Beam Geometry
//// Img: Image configuration
//// numPerSubSet: how many angles are there inside the subset 
//// subSetNum: how many subsets 
//// curSubSetIdx: current subset index;
//// blk: 32,32
//// gid: configured
////////////////////////////////////////////////////////////////////////////
//void bakproj_PIXEL(float* donp, float* dimg, const bool angDir, const bool FOV, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dimg, angDir, FOV, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx);
//}
//
//
//
//
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dimg, const bool angDir, const bool FOV, const FanEDGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	if (idX < Img.m_Reso.x && idY < Img.m_Reso.y)
//	{
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//		//cur rotation angle;
//		unsigned int angIdx = 0;
//		float curAng(0), cosT(0), sinT(0), summ(0);
//		Ray2D ray;
//		float2 minBox, maxBox, curImg, curPix;
//		float tnear(0), tfar(0), weg(0), detID(0);
//		int flrID, ceilID;
//
//		for (unsigned int prjIdx = 0; prjIdx != numPerSubSet; ++prjIdx)
//		{
//			angIdx = prjIdx * subSetNum + curSubSetIdx;
//			curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//
//			//current pixel position
//			curPix = rotation(make_float2(MINO.x + Img.m_Step.x * idX, MINO.y + Img.m_Step.y * idY), cosT, -sinT); //��ת���λ��;
//
//			//�����Ӧ��detidx;			//���ﲻ���� angDir; ���Ǵ�С�����;
//			detID = FanGeo.m_S2D * curPix.x / (curPix.y - FanGeo.m_S2O) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			if (detID < 0 || detID > FanGeo.m_DetN - 1)
//			{
//				dimg[idY * Img.m_Reso.x + idX] += 0;
//				continue;
//			}
//			flrID = detID;
//			ceilID = ceilf(detID);
//
//			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//			minBox = make_float2(
//				MINO.x + idX * Img.m_Step.x,
//				MINO.y + idY * Img.m_Step.y);
//			maxBox = minBox + make_float2(Img.m_Step.x, Img.m_Step.y);
//			curImg = minBox + 0.5f * make_float2(Img.m_Step.x, Img.m_Step.y);
//			ray.d = normalize(curImg - ray.o);
//
//			tnear = 0;
//			tfar = 0;
//			intersectBox(ray, minBox, maxBox, &tnear, &tfar); //�ཻ����;
//
//			weg = (tfar - tnear); //no weighting;
//			if (!IS_ZERO(flrID - ceilID))
//			{
//				summ += ((donp[prjIdx * FanGeo.m_DetN + flrID] * (ceilID - detID) + donp[prjIdx * FanGeo.m_DetN + ceilID] * (detID - flrID)) * weg);
//			}
//			else
//			{
//				summ += donp[prjIdx * FanGeo.m_DetN + flrID] * weg;
//			}
//		}
//		dimg[idY * Img.m_Reso.x + idX] = summ;
//	}
//}
//void bakproj_PIXEL(float* donp, float* dimg, const bool angDir, const bool FOV, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dimg, angDir, FOV, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx);
//}
//
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dvol, const bool angDir, const bool FOV, const ConeEAGeo ConeGeo, const Volume Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)//(float* donp, float* dvol, cuint angIdx, const bool angDir, const bool FOV, const ConeEAGeo ConeGeo, const Volume Vol)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	cuint idZ = threadIdx.z + blockDim.z * blockIdx.z;
//
//	if (idX < Vol.m_Reso.x && idY < Vol.m_Reso.y && idZ < Vol.m_Reso.z)
//	{
//		float3 MINO = make_float3(
//			-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//			-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//			-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//
//		unsigned int prjIdx = 0;
//		unsigned int angIdx = 0;
//		float curAng(0), cosT(0), sinT(0), biasAng(0), tnear(0),
//			tfar(0), weg(0), Qmm(0), Qmb(0), Qbm(0), Qbb(0), val(0), summ(0);
//		float3 curDet, initDet, sincosAng;
//		float2 detID;
//
//		float3 minBox = make_float3(
//			MINO.x + idX * Vol.m_Step.x,
//			MINO.y + idY * Vol.m_Step.y,
//			MINO.z + idZ * Vol.m_Step.z);
//		float3 maxBox = minBox + make_float3(Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z);
//		float3 curImg = minBox + 0.5f * make_float3(Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z);
//		int minX(0), minZ(0), maxX(0), maxZ(0);
//		Ray ray;
//		for (prjIdx = 0; prjIdx != numPerSubSet; ++prjIdx)
//		{
//			angIdx = prjIdx * subSetNum + curSubSetIdx;
//			curAng = ConeGeo.m_ViwBeg + angIdx * ConeGeo.m_ViwStp;
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//
//			ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
//			ray.d = normalize(curImg - ray.o);
//
//			//current detector element Cartesian coordinates;
//			curDet = ray.o + ray.d * ConeGeo.m_S2D;
//			//original detector position;
//			initDet = rotation(curDet, cosT, -sinT);  //�����λ��, ���ڿ�ʼ��Ӧ��index;
//
//			//Bias angle;
//			sincosAng.x = initDet.x / ConeGeo.m_S2D; //sin(ang)
//			sincosAng.y = (initDet.y - ConeGeo.m_S2O) / (-ConeGeo.m_S2D); //cos(ang)
//			biasAng = atan(sincosAng.x / sincosAng.y);
//
//			//judging the postive or negative angle direction;
//			if (angDir)
//			{
//				detID.x = biasAng / ConeGeo.m_DetStp + ConeGeo.m_DetCntIdx.x;
//			}
//			else
//			{
//				detID.x = ConeGeo.m_DetN - 1 - (biasAng / ConeGeo.m_DetStp + ConeGeo.m_DetCntIdx.x);
//			}
//
//			detID.y = (initDet.z / ConeGeo.m_DetHStp + ConeGeo.m_DetCntIdx.y);
//
//			if (detID.x < 0 || detID.x >(ConeGeo.m_DetN - 1) || detID.y < 0 || detID.y >(ConeGeo.m_DetHN - 1))
//			{
//				dvol[(idZ * Vol.m_Reso.z + idY) * Vol.m_Reso.x + idX] += 0;
//				continue;
//			}
//
//			minX = detID.x;
//			maxX = ceil(detID.x);
//			minZ = detID.y;
//			maxZ = ceil(detID.y);
//
//			Qmm = donp[(prjIdx * ConeGeo.m_DetHN + minZ) * ConeGeo.m_DetN + minX];
//			Qmb = donp[(prjIdx * ConeGeo.m_DetHN + maxZ) * ConeGeo.m_DetN + minX];
//			Qbm = donp[(prjIdx * ConeGeo.m_DetHN + minZ) * ConeGeo.m_DetN + maxX];
//			Qbb = donp[(prjIdx * ConeGeo.m_DetHN + maxZ) * ConeGeo.m_DetN + maxX];
//
//			val = _bilinearInterpolation_ker<float>(detID.x, detID.y, minX, maxX, minZ, maxZ, Qmm, Qmb, Qbm, Qbb);
//
//			tnear = 0;			tfar = 0;
//			intersectBox(ray, minBox, maxBox, &tnear, &tfar);
//			weg = (tfar - tnear); //no weighting;
//			summ += (weg * val);
//		}
//		dvol[(idZ * Vol.m_Reso.y + idY) * Vol.m_Reso.x + idX] += summ;
//	}
//}
//void bakproj_PIXEL(float* donp, float* dvol, const bool angDir, const bool FOV, const ConeEAGeo& ConeGeo, const Volume& Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dvol, angDir, FOV, ConeGeo, Vol, numPerSubSet, subSetNum, curSubSetIdx);
//}
//
//
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dvol, float* dmsk, const ConeEAGeo ConeGeo, const Volume Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	cuint idZ = threadIdx.z + blockDim.z * blockIdx.z;
//
//	if (idX < Vol.m_Reso.x && idY < Vol.m_Reso.y && idZ < Vol.m_Reso.z)
//	{
//		float3 MINO = make_float3(
//			-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//			-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//			-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//
//		unsigned int prjIdx = 0;
//		unsigned int angIdx = 0;
//		float curAng(0), cosT(0), sinT(0), biasAng(0), tnear(0),
//			tfar(0), weg(0), Qmm(0), Qmb(0), Qbm(0), Qbb(0), val(0), summ(0);
//		float3 curDet, initDet, sincosAng;
//		float2 detID;
//
//		float3 minBox = make_float3(
//			MINO.x + idX * Vol.m_Step.x,
//			MINO.y + idY * Vol.m_Step.y,
//			MINO.z + idZ * Vol.m_Step.z);
//		float3 maxBox = minBox + make_float3(Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z);
//		float3 curImg = minBox + 0.5f * make_float3(Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z);
//		int minX(0), minZ(0), maxX(0), maxZ(0);
//		Ray ray;
//		float mask = dmsk[idY * Vol.m_Reso.x + idX];
//		for (prjIdx = 0; prjIdx != numPerSubSet; ++prjIdx)
//		{
//			angIdx = prjIdx * subSetNum + curSubSetIdx;
//			curAng = ConeGeo.m_ViwBeg + angIdx * ConeGeo.m_ViwStp;
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//
//			ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
//			ray.d = normalize(curImg - ray.o);
//
//			//current detector element Cartesian coordinates;
//			curDet = ray.o + ray.d * ConeGeo.m_S2D;
//			//original detector position;
//			initDet = rotation(curDet, cosT, -sinT);  //�����λ��, ���ڿ�ʼ��Ӧ��index;
//
//			//Bias angle;
//			sincosAng.x = initDet.x / ConeGeo.m_S2D; //sin(ang)
//			sincosAng.y = (initDet.y - ConeGeo.m_S2O) / (-ConeGeo.m_S2D); //cos(ang)
//			biasAng = atan(sincosAng.x / sincosAng.y);
//
//			//judging the postive or negative angle direction;
//			detID.x = biasAng / ConeGeo.m_DetStp + ConeGeo.m_DetCntIdx.x;
//
//			detID.y = (initDet.z / ConeGeo.m_DetHStp + ConeGeo.m_DetCntIdx.y);
//
//			if (detID.x < 0 || detID.x >(ConeGeo.m_DetN - 1) || detID.y < 0 || detID.y >(ConeGeo.m_DetHN - 1))
//			{
//				dvol[(idZ * Vol.m_Reso.z + idY) * Vol.m_Reso.x + idX] += 0;
//				continue;
//			}
//
//			minX = detID.x;
//			maxX = ceil(detID.x);
//			minZ = detID.y;
//			maxZ = ceil(detID.y);
//
//			Qmm = donp[(prjIdx * ConeGeo.m_DetHN + minZ) * ConeGeo.m_DetN + minX];
//			Qmb = donp[(prjIdx * ConeGeo.m_DetHN + maxZ) * ConeGeo.m_DetN + minX];
//			Qbm = donp[(prjIdx * ConeGeo.m_DetHN + minZ) * ConeGeo.m_DetN + maxX];
//			Qbb = donp[(prjIdx * ConeGeo.m_DetHN + maxZ) * ConeGeo.m_DetN + maxX];
//
//			val = _bilinearInterpolation_ker<float>(detID.x, detID.y, minX, maxX, minZ, maxZ, Qmm, Qmb, Qbm, Qbb);
//
//			tnear = 0;			tfar = 0;
//			intersectBox(ray, minBox, maxBox, &tnear, &tfar);
//			weg = (tfar - tnear); //no weighting;
//			summ += (weg * val);
//		}
//		dvol[(idZ * Vol.m_Reso.y + idY) * Vol.m_Reso.x + idX] += (summ * mask);
//	}
//}
//void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, const ConeEAGeo& ConeGeo, const Volume& Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dimg, dmsk, ConeGeo, Vol, numPerSubSet, subSetNum, curSubSetIdx);
//}
//
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dvol, float* dmsk, const bool angDir, const ConeEAGeo ConeGeo, const Volume Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	cuint idZ = threadIdx.z + blockDim.z * blockIdx.z;
//
//	if (idX < Vol.m_Reso.x && idY < Vol.m_Reso.y && idZ < Vol.m_Reso.z)
//	{
//		float3 MINO = make_float3(
//			-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//			-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//			-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//
//		unsigned int prjIdx = 0;
//		unsigned int angIdx = 0;
//		float curAng(0), cosT(0), sinT(0), biasAng(0), tnear(0),
//			tfar(0), weg(0), Qmm(0), Qmb(0), Qbm(0), Qbb(0), val(0), summ(0);
//		float3 curDet, initDet, sincosAng;
//		float2 detID;
//
//		float3 minBox = make_float3(
//			MINO.x + idX * Vol.m_Step.x,
//			MINO.y + idY * Vol.m_Step.y,
//			MINO.z + idZ * Vol.m_Step.z);
//		float3 maxBox = minBox + make_float3(Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z);
//		float3 curImg = minBox + 0.5f * make_float3(Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z);
//		int minX(0), minZ(0), maxX(0), maxZ(0);
//		Ray ray;
//		float mask = dmsk[idY * Vol.m_Reso.x + idX];
//		for (prjIdx = 0; prjIdx != numPerSubSet; ++prjIdx)
//		{
//			angIdx = prjIdx * subSetNum + curSubSetIdx;
//			curAng = ConeGeo.m_ViwBeg + angIdx * ConeGeo.m_ViwStp;
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//
//			ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
//			ray.d = normalize(curImg - ray.o);
//
//			//current detector element Cartesian coordinates;
//			curDet = ray.o + ray.d * ConeGeo.m_S2D;
//			//original detector position;
//			initDet = rotation(curDet, cosT, -sinT);  //�����λ��, ���ڿ�ʼ��Ӧ��index;
//
//			//Bias angle;
//			sincosAng.x = initDet.x / ConeGeo.m_S2D; //sin(ang)
//			sincosAng.y = (initDet.y - ConeGeo.m_S2O) / (-ConeGeo.m_S2D); //cos(ang)
//			biasAng = atan(sincosAng.x / sincosAng.y);
//
//			//judging the postive or negative angle direction;
//			if (angDir)
//			{
//				detID.x = biasAng / ConeGeo.m_DetStp + ConeGeo.m_DetCntIdx.x;
//			}
//			else
//			{
//				detID.x = ConeGeo.m_DetN - 1 - (biasAng / ConeGeo.m_DetStp + ConeGeo.m_DetCntIdx.x);
//			}
//
//			detID.y = (initDet.z / ConeGeo.m_DetHStp + ConeGeo.m_DetCntIdx.y);
//
//			if (detID.x < 0 || detID.x >(ConeGeo.m_DetN - 1) || detID.y < 0 || detID.y >(ConeGeo.m_DetHN - 1))
//			{
//				dvol[(idZ * Vol.m_Reso.z + idY) * Vol.m_Reso.x + idX] += 0;
//				continue;
//			}
//
//			minX = detID.x;
//			maxX = ceil(detID.x);
//			minZ = detID.y;
//			maxZ = ceil(detID.y);
//
//			Qmm = donp[(prjIdx * ConeGeo.m_DetHN + minZ) * ConeGeo.m_DetN + minX];
//			Qmb = donp[(prjIdx * ConeGeo.m_DetHN + maxZ) * ConeGeo.m_DetN + minX];
//			Qbm = donp[(prjIdx * ConeGeo.m_DetHN + minZ) * ConeGeo.m_DetN + maxX];
//			Qbb = donp[(prjIdx * ConeGeo.m_DetHN + maxZ) * ConeGeo.m_DetN + maxX];
//
//			val = _bilinearInterpolation_ker<float>(detID.x, detID.y, minX, maxX, minZ, maxZ, Qmm, Qmb, Qbm, Qbb);
//
//			tnear = 0;			tfar = 0;
//			intersectBox(ray, minBox, maxBox, &tnear, &tfar);
//			weg = (tfar - tnear); //no weighting;
//			summ += (weg * val);
//		}
//		dvol[(idZ * Vol.m_Reso.y + idY) * Vol.m_Reso.x + idX] += (summ * mask);
//	}
//}
//void bakproj_PIXEL(float* donp, float* dvol, float* dmsk, const bool angDir, const ConeEAGeo& ConeGeo, const Volume& Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dvol, dmsk, angDir, ConeGeo, Vol, numPerSubSet, subSetNum, curSubSetIdx);
//}
//
//
//
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dimg, const FanEAGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx) //;(float* donp, float* dimg, const FanEAGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx) 
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	if (idX < Img.m_Reso.x && idY < Img.m_Reso.y)
//	{
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//		//cur rotation angle;
//		unsigned int angIdx = 0;
//		float curAng, cosT, sinT;
//		Ray2D ray;
//		float2 minBox, maxBox, curImg, curDet, initDet, sincosAng;
//		float biasAng, tnear, tfar, weg, detID;
//		int flrID, ceilID;
//		float summ(0);
//		for (unsigned int prjIdx = 0; prjIdx != numPerSubSet; ++prjIdx)
//		{
//			angIdx = prjIdx * subSetNum + curSubSetIdx;
//			curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//			//current source position, assume the initial position is on the positive Y axis;
//			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//			minBox = make_float2(
//				MINO.x + idX * Img.m_Step.x,
//				MINO.y + idY * Img.m_Step.y);
//			maxBox = minBox + make_float2(Img.m_Step.x, Img.m_Step.y);
//			curImg = minBox + 0.5f * make_float2(Img.m_Step.x, Img.m_Step.y);
//			ray.d = normalize(curImg - ray.o);
//
//
//			//current detector element Cartesian coordinates;
//			curDet = ray.o + ray.d * FanGeo.m_S2D;
//			//original detector position;
//			initDet = rotation(curDet, cosT, -sinT);
//
//			//Bias angle;
//
//			sincosAng.x = initDet.x / FanGeo.m_S2D; //sin(ang)
//			sincosAng.y = (initDet.y - FanGeo.m_S2O) / (-FanGeo.m_S2D); //cos(ang)
//			biasAng = atan(sincosAng.x / sincosAng.y);
//
//			//judging the postive or negative angle direction;
//			detID = biasAng / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			if (detID < 0)
//			{
//				dimg[idY * Img.m_Reso.x + idX] += 0;
//				continue;
//			}
//			if (detID > FanGeo.m_DetN)
//			{
//				dimg[idY * Img.m_Reso.x + idX] += 0;
//				continue;
//			}
//			tnear = 0;
//			tfar = 0;
//			flrID = detID;
//			ceilID = ceilf(detID);
//			intersectBox(ray, minBox, maxBox, &tnear, &tfar);
//
//			weg = (tfar - tnear); //no weighting;
//			if (!IS_ZERO(flrID - ceilID))
//			{
//				summ += ((donp[prjIdx * FanGeo.m_DetN + flrID] * (ceilID - detID) + donp[prjIdx * FanGeo.m_DetN + ceilID] * (detID - flrID)) * weg);
//			}
//			else
//			{
//				summ += donp[prjIdx * FanGeo.m_DetN + flrID] * weg;
//			}
//		}
//		dimg[idY * Img.m_Reso.x + idX] = summ;
//	}
//}
////////////////////////////////////////////////////////////////////////////
//// Backprojection from the subset with given subset index
//// donp: the projection data from subset
//// dimg: image
//// FanGeo: Fan Beam Geometry
//// Img: Image configuration
//// numPerSubSet: how many angles are there inside the subset 
//// subSetNum: how many subsets 
//// curSubSetIdx: current subset index;
//// blk: 32,32
//// gid: configured
////////////////////////////////////////////////////////////////////////////
//void bakproj_PIXEL(float* donp, float* dimg, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dimg, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx);
//}
//
//
////
////__global__ void _bakproj_PIXEL_Ker(float* donp, float* dimg, float* dmask, const FanEDGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
////{
////	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
////	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
////	if (idX < Img.m_Reso.x && idY < Img.m_Reso.y)
////	{
////
////		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
////		//cur rotation angle;
////		//unsigned int angIdx = 0;
////		float curAng(0), cosT(0), sinT(0), summ(0),mask(dmask[idY * Img.m_Reso.x + idX]);
////		Ray2D ray;
////		float2 minBox, maxBox, curImg, curPix;
////		float tnear(0), tfar(0), weg(0), detID(0);
////		int flrID, ceilID;
////		unsigned int angIdx(0);
////
////		for (unsigned int prjIdx = 0; prjIdx != numPerSubSet; ++prjIdx)
////		{
////			angIdx = prjIdx * subSetNum + curSubSetIdx;
////			curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
////			cosT = cosf(curAng);
////			sinT = sinf(curAng);
////
////			//current pixel position
////			curPix = rotation(make_float2(
////				MINO.x + Img.m_Step.x * (idX + 0.5f),
////				MINO.y + Img.m_Step.y * (idY + 0.5f)),cosT,-sinT); //��ת���λ��;
////			//���λ����
////
////			//�����Ӧ��detidx;			//���ﲻ���� angDir; ���Ǵ�С�����;
////			detID = curPix.x / (FanGeo.m_S2O - curPix.y) * FanGeo.m_S2D / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;//FanGeo.m_S2D * curPix.x / (curPix.y - FanGeo.m_S2O) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
////			if (detID < 0 || detID > FanGeo.m_DetN - 1)
////			{
////				dimg[idY * Img.m_Reso.x + idX] += 0.0f;
////				return;
////			}
////			flrID = detID;
////			ceilID = ceilf(detID);
////
////			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
////			minBox = make_float2(
////				MINO.x + idX * Img.m_Step.x,
////				MINO.y + idY * Img.m_Step.y);
////			maxBox = minBox + make_float2(Img.m_Step.x, Img.m_Step.y);
////			curImg = minBox + 0.5f * make_float2(Img.m_Step.x, Img.m_Step.y);
////			ray.d = normalize(curImg - ray.o);
////
////			tnear = 0;
////			tfar = 0;
////			intersectBox(ray, minBox, maxBox, &tnear, &tfar); //�ཻ����;
////
////			weg = (tfar - tnear); //no weighting;
////			if (!IS_ZERO(flrID - ceilID))
////			{
////				//dimg[idY * Img.m_Reso.x + idX] += 1.0f;
////				summ += ((donp[prjIdx * FanGeo.m_DetN + flrID] * (ceilID - detID) + donp[prjIdx * FanGeo.m_DetN + ceilID] * (detID - flrID)) * weg);
////			}
////			else
////			{
////				//dimg[idY * Img.m_Reso.x + idX] += 100.0f;
////				summ += donp[prjIdx * FanGeo.m_DetN + flrID] * weg;
////			}
////		}
////		dimg[idY * Img.m_Reso.x + idX] = (summ * mask);
////	}
////}
////void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
////{	_bakproj_PIXEL_Ker<<<gid,blk>>>(donp, dimg, dmsk, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx); }
//
//
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dimg, float* dmask, const ConeEDGeo ConeGeo, const Volume Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	cuint idZ = threadIdx.z + blockDim.z * blockIdx.z;
//	if (idX < Vol.m_Reso.x && idY < Vol.m_Reso.y && idZ < Vol.m_Reso.z)
//	{
//
//		float3 MINO = make_float3(
//			-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//			-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//			-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//		//cur rotation angle;
//		//unsigned int angIdx = 0;
//		float cosT(0), sinT(0), summ(0), mask(dmask[idY * Vol.m_Reso.x + idX]);
//		Ray ray;
//		float3 minBox, maxBox, curImg, curPix;
//		float tnear(0), tfar(0), weg(0);
//		float2 detID;
//		int2 flrID, ceilID;
//		unsigned int angIdx(0);
//		float Qmm(0), Qmb(0), Qbm(0), Qbb(0);
//		for (unsigned int prjIdx = 0; prjIdx != numPerSubSet; ++prjIdx)
//		{
//			angIdx = prjIdx * subSetNum + curSubSetIdx;
//			cosT = cosf(ConeGeo.m_ViwBeg + ConeGeo.m_ViwStp * angIdx);
//			sinT = sinf(ConeGeo.m_ViwBeg + ConeGeo.m_ViwStp * angIdx);
//
//			//current pixel position
//			curPix = rotation(make_float3(
//				MINO.x + Vol.m_Step.x * (idX + 0.5f),
//				MINO.y + Vol.m_Step.y * (idY + 0.5f),
//				MINO.z + Vol.m_Step.z * (idZ + 0.5f)), cosT, -sinT); //��ת���λ��;
//
//			detID.x = curPix.x / (ConeGeo.m_S2O - curPix.y) * ConeGeo.m_S2D / ConeGeo.m_DetStp.x + ConeGeo.m_DetCntIdx.x;
//			detID.y = curPix.z / (ConeGeo.m_S2O - curPix.y) * ConeGeo.m_S2D / ConeGeo.m_DetStp.y + ConeGeo.m_DetCntIdx.y;
//
//			if (detID.x < 0 || detID.y < 0 || detID.x > ConeGeo.m_DetN.x - 1 || detID.y > ConeGeo.m_DetN.y - 1)
//			{
//				dimg[(idZ * Vol.m_Reso.y + idY) * Vol.m_Reso.x + idX] += 0.0f;
//				continue;
//			}
//			flrID.x = detID.x;
//			flrID.y = detID.y;
//			ceilID.x = ceilf(detID.x);
//			ceilID.y = ceilf(detID.y);
//
//			Qmm = donp[(prjIdx * ConeGeo.m_DetN.y + flrID.y) * ConeGeo.m_DetN.x + flrID.x];
//			Qmb = donp[(prjIdx * ConeGeo.m_DetN.y + ceilID.y) * ConeGeo.m_DetN.x + flrID.x];
//			Qbm = donp[(prjIdx * ConeGeo.m_DetN.y + flrID.y) * ConeGeo.m_DetN.x + ceilID.x];
//			Qbb = donp[(prjIdx * ConeGeo.m_DetN.y + ceilID.y) * ConeGeo.m_DetN.x + ceilID.x];
//
//
//			ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
//			minBox = make_float3(
//				MINO.x + idX * Vol.m_Step.x,
//				MINO.y + idY * Vol.m_Step.y,
//				MINO.z + idZ * Vol.m_Step.z);
//
//			maxBox = minBox + make_float3(Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z);
//			curImg = minBox + 0.5f * make_float3(Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z);
//			ray.d = normalize(curImg - ray.o);
//
//			tnear = 0;
//			tfar = 0;
//			intersectBox(ray, minBox, maxBox, &tnear, &tfar); //�ཻ����;
//
//			weg = (tfar - tnear); //no weighting;
//
//
//			summ += weg * _bilinearInterpolation_ker<float>((float)detID.x, (float)detID.y,
//				(float)flrID.x, (float)ceilID.x, (float)flrID.y, (float)ceilID.y, Qmm, Qmb, Qbm, Qbb);
//		}
//		dimg[(idZ * Vol.m_Reso.y + idY) * Vol.m_Reso.x + idX] = (summ * mask);
//	}
//}
//void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, const ConeEDGeo& ConeGeo, const Volume& Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dimg, dmsk, ConeGeo, Vol, numPerSubSet, subSetNum, curSubSetIdx);
//}
//
//
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dimg, const ConeEDGeo& ConeGeo, const Volume& Vol, const int angIdx)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	cuint idZ = threadIdx.z + blockDim.z * blockIdx.z;
//	if (idX < Vol.m_Reso.x && idY < Vol.m_Reso.y && idZ < Vol.m_Reso.z)
//	{
//
//		float3 MINO = make_float3(
//			-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//			-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//			-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//		//cur rotation angle;
//		//unsigned int angIdx = 0;
//		float cosT(0), sinT(0), summ(0);
//		Ray ray;
//		float3 minBox, maxBox, curImg, curPix;
//		float tnear(0), tfar(0), weg(0);
//		float2 detID;
//		int2 flrID, ceilID;
//		unsigned int angIdx(0);
//		float Qmm(0), Qmb(0), Qbm(0), Qbb(0);
//		cosT = cosf(ConeGeo.m_ViwBeg + ConeGeo.m_ViwStp * angIdx);
//		sinT = sinf(ConeGeo.m_ViwBeg + ConeGeo.m_ViwStp * angIdx);
//
//		//current pixel position
//		curPix = rotation(make_float3(
//			MINO.x + Vol.m_Step.x * (idX + 0.5f),
//			MINO.y + Vol.m_Step.y * (idY + 0.5f),
//			MINO.z + Vol.m_Step.z * (idZ + 0.5f)), cosT, -sinT); //��ת���λ��;
//
//		detID.x = curPix.x / (ConeGeo.m_S2O - curPix.y) * ConeGeo.m_S2D / ConeGeo.m_DetStp.x + ConeGeo.m_DetCntIdx.x;
//		detID.y = curPix.z / (ConeGeo.m_S2O - curPix.y) * ConeGeo.m_S2D / ConeGeo.m_DetStp.y + ConeGeo.m_DetCntIdx.y;
//
//
//		flrID.x = detID.x;
//		flrID.y = detID.y;
//		ceilID.x = ceilf(detID.x);
//		ceilID.y = ceilf(detID.y);
//		if (flrID.x < 0 || flrID.x > ConeGeo.m_DetN.x - 1 ||
//			flrID.y < 0 || flrID.y > ConeGeo.m_DetN.y - 1 ||
//			ceilID.x < 0 || ceilID.x > ConeGeo.m_DetN.x - 1 ||
//			ceilID.y < 0 || ceilID.y > ConeGeo.m_DetN.y - 1)
//		{
//			dimg[(idZ * Vol.m_Reso.y + idY) * Vol.m_Reso.x + idX] += 0.0f;
//			return;
//		}
//		Qmm = donp[flrID.y * ConeGeo.m_DetN.x + flrID.x];
//		Qmb = donp[ceilID.y * ConeGeo.m_DetN.x + flrID.x];
//		Qbm = donp[flrID.y * ConeGeo.m_DetN.x + ceilID.x];
//		Qbb = donp[ceilID.y * ConeGeo.m_DetN.x + ceilID.x];
//
//
//		ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
//		minBox = make_float3(
//			MINO.x + idX * Vol.m_Step.x,
//			MINO.y + idY * Vol.m_Step.y,
//			MINO.z + idZ * Vol.m_Step.z);
//
//		maxBox = minBox + make_float3(Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z);
//		curImg = minBox + 0.5f * make_float3(Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z);
//		ray.d = normalize(curImg - ray.o);
//
//		tnear = 0;
//		tfar = 0;
//		intersectBox(ray, minBox, maxBox, &tnear, &tfar); //�ཻ����;
//
//		weg = (tfar - tnear); //no weighting;
//
//		summ = weg * _bilinearInterpolation_ker<float>((float)detID.x, (float)detID.y,
//			(float)flrID.x, (float)ceilID.x, (float)flrID.y, (float)ceilID.y, Qmm, Qmb, Qbm, Qbb);
//		dimg[(idZ * Vol.m_Reso.y + idY) * Vol.m_Reso.x + idX] += summ;
//	}
//
//}
//void back_proj_DEMO18v3(float* donp, float* dimg, const ConeEDGeo& ConeGeo, const Volume& Vol, const int angIdx, const dim3& blockSize, const dim3& gridSize)
//{
//	_bakproj_PIXEL_Ker << <gridSize, blockSize >> >(donp, dimg, ConeGeo, Vol, angIdx);
//}
//
//
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dvol, const ConeEAGeo ConeGeo, const Volume Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	cuint idZ = threadIdx.z + blockDim.z * blockIdx.z;
//
//	if (idX < Vol.m_Reso.x && idY < Vol.m_Reso.y && idZ < Vol.m_Reso.z)
//	{
//		float3 MINO = make_float3(
//			-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//			-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//			-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//
//		unsigned int prjIdx = 0;
//		unsigned int angIdx = 0;
//		float curAng(0), cosT(0), sinT(0), biasAng(0), tnear(0),
//			tfar(0), weg(0), Qmm(0), Qmb(0), Qbm(0), Qbb(0), val(0), summ(0);
//		float3 curDet, initDet, sincosAng;
//		float2 detID;
//
//		float3 minBox = make_float3(
//			MINO.x + idX * Vol.m_Step.x,
//			MINO.y + idY * Vol.m_Step.y,
//			MINO.z + idZ * Vol.m_Step.z);
//		float3 maxBox = minBox + make_float3(Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z);
//		float3 curImg = minBox + 0.5f * make_float3(Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z);
//		int minX(0), minZ(0), maxX(0), maxZ(0);
//		Ray ray;
//		for (prjIdx = 0; prjIdx != numPerSubSet; ++prjIdx)
//		{
//			angIdx = prjIdx * subSetNum + curSubSetIdx;
//			curAng = ConeGeo.m_ViwBeg + angIdx * ConeGeo.m_ViwStp;
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//
//			ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
//			ray.d = normalize(curImg - ray.o);
//
//			//current detector element Cartesian coordinates;
//			curDet = ray.o + ray.d * ConeGeo.m_S2D;
//			//original detector position;
//			initDet = rotation(curDet, cosT, -sinT);  //�����λ��, ���ڿ�ʼ��Ӧ��index;
//
//			//Bias angle;
//			sincosAng.x = initDet.x / ConeGeo.m_S2D; //sin(ang)
//			sincosAng.y = (initDet.y - ConeGeo.m_S2O) / (-ConeGeo.m_S2D); //cos(ang)
//			biasAng = atan(sincosAng.x / sincosAng.y);
//
//			//judging the postive or negative angle direction;
//			detID.x = biasAng / ConeGeo.m_DetStp + ConeGeo.m_DetCntIdx.x;
//
//			detID.y = (initDet.z / ConeGeo.m_DetHStp + ConeGeo.m_DetCntIdx.y);
//
//			if (detID.x < 0 || detID.x >(ConeGeo.m_DetN - 1) || detID.y < 0 || detID.y >(ConeGeo.m_DetHN - 1))
//			{
//				dvol[(idZ * Vol.m_Reso.z + idY) * Vol.m_Reso.x + idX] += 0;
//				continue;
//			}
//
//			minX = detID.x;
//			maxX = ceil(detID.x);
//			minZ = detID.y;
//			maxZ = ceil(detID.y);
//
//			Qmm = donp[(prjIdx * ConeGeo.m_DetHN + minZ) * ConeGeo.m_DetN + minX];
//			Qmb = donp[(prjIdx * ConeGeo.m_DetHN + maxZ) * ConeGeo.m_DetN + minX];
//			Qbm = donp[(prjIdx * ConeGeo.m_DetHN + minZ) * ConeGeo.m_DetN + maxX];
//			Qbb = donp[(prjIdx * ConeGeo.m_DetHN + maxZ) * ConeGeo.m_DetN + maxX];
//
//			val = _bilinearInterpolation_ker<float>(detID.x, detID.y, minX, maxX, minZ, maxZ, Qmm, Qmb, Qbm, Qbb);
//
//			tnear = 0;			tfar = 0;
//			intersectBox(ray, minBox, maxBox, &tnear, &tfar);
//			weg = (tfar - tnear); //no weighting;
//
//			summ += (weg * val);
//
//		}
//		dvol[(idZ * Vol.m_Reso.y + idY) * Vol.m_Reso.x + idX] += summ;
//	}
//}
//void bakproj_PIXEL(float* donp, float* dvol, const ConeEAGeo& ConeGeo, const Volume& Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dvol, ConeGeo, Vol, numPerSubSet, subSetNum, curSubSetIdx);
//}
//
//
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dimg, const FanEAGeo FanGeo, const Image Img)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	if (idX < Img.m_Reso.x && idY < Img.m_Reso.y)
//	{
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//		//cur rotation angle;
//		float curAng, cosT, sinT;
//		Ray2D ray;
//		float2 minBox, maxBox, curImg, curDet, initDet, sincosAng;
//		float biasAng, tnear, tfar, weg, detID;
//		int flrID, ceilID;
//		float summ(0);
//		for (unsigned int angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
//		{
//			curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//			//current source position, assume the initial position is on the positive Y axis;
//			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//			minBox = make_float2(
//				MINO.x + idX * Img.m_Step.x,
//				MINO.y + idY * Img.m_Step.y);
//			maxBox = minBox + make_float2(Img.m_Step.x, Img.m_Step.y);
//			curImg = minBox + 0.5f * make_float2(Img.m_Step.x, Img.m_Step.y);
//			ray.d = normalize(curImg - ray.o);
//
//
//			//current detector element Cartesian coordinates;
//			curDet = ray.o + ray.d * FanGeo.m_S2D;
//			//original detector position;
//			initDet = rotation(curDet, cosT, -sinT);
//
//			//Bias angle;
//
//			sincosAng.x = initDet.x / FanGeo.m_S2D; //sin(ang)
//			sincosAng.y = (initDet.y - FanGeo.m_S2O) / (-FanGeo.m_S2D); //cos(ang)
//			biasAng = atan(sincosAng.x / sincosAng.y);
//
//			//judging the postive or negative angle direction;
//			detID = biasAng / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			if (detID < 0)
//			{
//				dimg[idY * Img.m_Reso.x + idX] += 0;//summ += 0;
//				continue;
//			}
//			if (detID > FanGeo.m_DetN - 1)
//			{
//				dimg[idY * Img.m_Reso.x + idX] += 0;
//				continue;
//			}
//			tnear = 0;
//			tfar = 0;
//			flrID = detID;
//			ceilID = ceilf(detID);
//			intersectBox(ray, minBox, maxBox, &tnear, &tfar);
//
//			weg = (tfar - tnear); //no weighting;
//			if (!IS_ZERO(flrID - ceilID))
//			{
//				summ += (
//					(
//					donp[angIdx * FanGeo.m_DetN + flrID] * (ceilID - detID) +
//					donp[angIdx * FanGeo.m_DetN + ceilID] * (detID - flrID)
//					) * weg);
//			}
//			else
//			{
//				summ += donp[angIdx * FanGeo.m_DetN + flrID] * weg;
//			}
//		}
//		dimg[idY * Img.m_Reso.x + idX] = summ;
//	}
//}
////////////////////////////////////////////////////////////////////////////
//// Backprojection from all the angles
//// donp: projection data
//// dimg: image;
//// FanGeo: Fan Beam Geometry
//// Img: Image configuration
//// blk: 32,32
//// gid: configured
////////////////////////////////////////////////////////////////////////////
//void bakproj_PIXEL(float* donp, float* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dimg, FanGeo, Img);
//}
//
//
//
//
//
//
//
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dimg, const FanEDGeo FanGeo, const Image Img)//;(float* donp, float* dimg, const FanEDGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	if (idX < Img.m_Reso.x && idY < Img.m_Reso.y)
//	{
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//		//cur rotation angle;
//		//unsigned int angIdx = 0;
//		float curAng(0), cosT(0), sinT(0), summ(0);
//		Ray2D ray;
//		float2 minBox, maxBox, curImg, curPix;
//		float tnear(0), tfar(0), weg(0), detID(0);
//		int flrID, ceilID;
//
//		for (unsigned int angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
//		{
//			//angIdx = prjIdx * subSetNum + curSubSetIdx;
//			curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//
//			//current pixel position
//			curPix = rotation(make_float2(MINO.x + Img.m_Step.x * idX, MINO.y + Img.m_Step.y * idY), cosT, -sinT); //��ת���λ��;
//
//			//�����Ӧ��detidx;			//���ﲻ���� angDir; ���Ǵ�С�����;
//			detID = FanGeo.m_S2D * curPix.x / (curPix.y - FanGeo.m_S2O) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			if (detID < 0 || detID > FanGeo.m_DetN - 1)
//			{
//				dimg[idY * Img.m_Reso.x + idX] += 0;
//				continue;
//			}
//			flrID = detID;
//			ceilID = ceilf(detID);
//
//			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//			minBox = make_float2(
//				MINO.x + idX * Img.m_Step.x,
//				MINO.y + idY * Img.m_Step.y);
//			maxBox = minBox + make_float2(Img.m_Step.x, Img.m_Step.y);
//			curImg = minBox + 0.5f * make_float2(Img.m_Step.x, Img.m_Step.y);
//			ray.d = normalize(curImg - ray.o);
//
//			tnear = 0;
//			tfar = 0;
//			intersectBox(ray, minBox, maxBox, &tnear, &tfar); //�ཻ����;
//
//			weg = (tfar - tnear); //no weighting;
//			if (!IS_ZERO(flrID - ceilID))
//			{
//				summ += ((donp[angIdx * FanGeo.m_DetN + flrID] * (ceilID - detID) + donp[angIdx * FanGeo.m_DetN + ceilID] * (detID - flrID)) * weg);
//			}
//			else
//			{
//				summ += donp[angIdx * FanGeo.m_DetN + flrID] * weg;
//			}
//		}
//		dimg[idY * Img.m_Reso.x + idX] = summ;
//	}
//}
//void bakproj_PIXEL(float* donp, float* dimg, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dimg, FanGeo, Img);
//}
//
//
//
//
//
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dvol, const ConeEAGeo ConeGeo, const Volume Vol)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	cuint idZ = threadIdx.z + blockDim.z * blockIdx.z;
//
//	if (idX < Vol.m_Reso.x && idY < Vol.m_Reso.y && idZ < Vol.m_Reso.z)
//	{
//		float3 MINO = make_float3(
//			-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//			-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//			-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//
//		unsigned int prjIdx = 0;
//		unsigned int angIdx = 0;
//		float curAng(0), cosT(0), sinT(0), biasAng(0), tnear(0),
//			tfar(0), weg(0), Qmm(0), Qmb(0), Qbm(0), Qbb(0), val(0), summ(0);
//		float3 curDet, initDet, sincosAng;
//		float2 detID;
//
//		float3 minBox = make_float3(
//			MINO.x + idX * Vol.m_Step.x,
//			MINO.y + idY * Vol.m_Step.y,
//			MINO.z + idZ * Vol.m_Step.z);
//		float3 maxBox = minBox + make_float3(Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z);
//		float3 curImg = minBox + 0.5f * make_float3(Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z);
//		int minX(0), minZ(0), maxX(0), maxZ(0);
//		Ray ray;
//		for (angIdx = 0; angIdx != ConeGeo.m_ViwN; ++prjIdx)
//		{
//			curAng = ConeGeo.m_ViwBeg + angIdx * ConeGeo.m_ViwStp;
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//
//			ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
//			ray.d = normalize(curImg - ray.o);
//
//			//current detector element Cartesian coordinates;
//			curDet = ray.o + ray.d * ConeGeo.m_S2D;
//			//original detector position;
//			initDet = rotation(curDet, cosT, -sinT);  //�����λ��, ���ڿ�ʼ��Ӧ��index;
//
//			//Bias angle;
//			sincosAng.x = initDet.x / ConeGeo.m_S2D; //sin(ang)
//			sincosAng.y = (initDet.y - ConeGeo.m_S2O) / (-ConeGeo.m_S2D); //cos(ang)
//			biasAng = atan(sincosAng.x / sincosAng.y);
//
//			//judging the postive or negative angle direction;
//			detID.x = biasAng / ConeGeo.m_DetStp + ConeGeo.m_DetCntIdx.x;
//
//			detID.y = (initDet.z / ConeGeo.m_DetHStp + ConeGeo.m_DetCntIdx.y);
//
//			if (detID.x < 0 || detID.x >(ConeGeo.m_DetN - 1) || detID.y < 0 || detID.y >(ConeGeo.m_DetHN - 1))
//			{
//				dvol[(idZ * Vol.m_Reso.z + idY) * Vol.m_Reso.x + idX] += 0;
//				continue;
//			}
//
//			minX = detID.x;
//			maxX = ceil(detID.x);
//			minZ = detID.y;
//			maxZ = ceil(detID.y);
//
//			Qmm = donp[(angIdx * ConeGeo.m_DetHN + minZ) * ConeGeo.m_DetN + minX];
//			Qmb = donp[(angIdx * ConeGeo.m_DetHN + maxZ) * ConeGeo.m_DetN + minX];
//			Qbm = donp[(angIdx * ConeGeo.m_DetHN + minZ) * ConeGeo.m_DetN + maxX];
//			Qbb = donp[(angIdx * ConeGeo.m_DetHN + maxZ) * ConeGeo.m_DetN + maxX];
//
//			val = _bilinearInterpolation_ker<float>(detID.x, detID.y, minX, maxX, minZ, maxZ, Qmm, Qmb, Qbm, Qbb);
//
//			tnear = 0;			tfar = 0;
//			intersectBox(ray, minBox, maxBox, &tnear, &tfar);
//			weg = (tfar - tnear); //no weighting;
//
//			summ += (weg * val);
//
//		}
//		dvol[(idZ * Vol.m_Reso.y + idY) * Vol.m_Reso.x + idX] += summ;
//	}
//}
//void bakproj_PIXEL(float* donp, float* dvol, const ConeEAGeo& CoeGeo, const Volume& Vol, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dvol, CoeGeo, Vol);
//}
//
//
//
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dvol, float* dmsk, const ConeEAGeo ConeGeo, const Volume Vol)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	cuint idZ = threadIdx.z + blockDim.z * blockIdx.z;
//
//	if (idX < Vol.m_Reso.x && idY < Vol.m_Reso.y && idZ < Vol.m_Reso.z)
//	{
//		float3 MINO = make_float3(
//			-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//			-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//			-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//
//		unsigned int prjIdx = 0;
//		unsigned int angIdx = 0;
//		float curAng(0), cosT(0), sinT(0), biasAng(0), tnear(0),
//			tfar(0), weg(0), Qmm(0), Qmb(0), Qbm(0), Qbb(0), val(0), summ(0);
//		float3 curDet, initDet, sincosAng;
//		float2 detID;
//
//		float3 minBox = make_float3(
//			MINO.x + idX * Vol.m_Step.x,
//			MINO.y + idY * Vol.m_Step.y,
//			MINO.z + idZ * Vol.m_Step.z);
//		float3 maxBox = minBox + make_float3(Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z);
//		float3 curImg = minBox + 0.5f * make_float3(Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z);
//		int minX(0), minZ(0), maxX(0), maxZ(0);
//		Ray ray;
//		float mask = dmsk[idY * Vol.m_Reso.x + idX];
//		for (angIdx = 0; angIdx != ConeGeo.m_ViwN; ++prjIdx)
//		{
//			curAng = ConeGeo.m_ViwBeg + angIdx * ConeGeo.m_ViwStp;
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//
//			ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
//			ray.d = normalize(curImg - ray.o);
//
//			//current detector element Cartesian coordinates;
//			curDet = ray.o + ray.d * ConeGeo.m_S2D;
//			//original detector position;
//			initDet = rotation(curDet, cosT, -sinT);  //�����λ��, ���ڿ�ʼ��Ӧ��index;
//
//			//Bias angle;
//			sincosAng.x = initDet.x / ConeGeo.m_S2D; //sin(ang)
//			sincosAng.y = (initDet.y - ConeGeo.m_S2O) / (-ConeGeo.m_S2D); //cos(ang)
//			biasAng = atan(sincosAng.x / sincosAng.y);
//
//			//judging the postive or negative angle direction;
//			detID.x = biasAng / ConeGeo.m_DetStp + ConeGeo.m_DetCntIdx.x;
//
//			detID.y = (initDet.z / ConeGeo.m_DetHStp + ConeGeo.m_DetCntIdx.y);
//
//			if (detID.x < 0 || detID.x >(ConeGeo.m_DetN - 1) || detID.y < 0 || detID.y >(ConeGeo.m_DetHN - 1))
//			{
//				dvol[(idZ * Vol.m_Reso.z + idY) * Vol.m_Reso.x + idX] += 0;
//				continue;
//			}
//
//			minX = detID.x;
//			maxX = ceil(detID.x);
//			minZ = detID.y;
//			maxZ = ceil(detID.y);
//
//			Qmm = donp[(angIdx * ConeGeo.m_DetHN + minZ) * ConeGeo.m_DetN + minX];
//			Qmb = donp[(angIdx * ConeGeo.m_DetHN + maxZ) * ConeGeo.m_DetN + minX];
//			Qbm = donp[(angIdx * ConeGeo.m_DetHN + minZ) * ConeGeo.m_DetN + maxX];
//			Qbb = donp[(angIdx * ConeGeo.m_DetHN + maxZ) * ConeGeo.m_DetN + maxX];
//
//			val = _bilinearInterpolation_ker<float>(detID.x, detID.y, minX, maxX, minZ, maxZ, Qmm, Qmb, Qbm, Qbb);
//
//			tnear = 0;			tfar = 0;
//			intersectBox(ray, minBox, maxBox, &tnear, &tfar);
//			weg = (tfar - tnear); //no weighting;
//
//			summ += (weg * val);
//
//		}
//		dvol[(idZ * Vol.m_Reso.y + idY) * Vol.m_Reso.x + idX] += (summ * mask);
//	}
//}
//void bakproj_PIXEL(float* donp, float* dvol, float* dmask, const ConeEAGeo& CoeGeo, const Volume& Vol, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dvol, dmask, CoeGeo, Vol);
//}
//
//
//
//
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dimg, float* dmsk, cuint angIdx, const bool angDir, const FanEAGeo FanGeo, const Image Img)//;(float* donp, float* dimg, cuint angIdx, const bool angDir, const bool FOV, const FanEAGeo FanGeo, const Image Img)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	if (idX < Img.m_Reso.x && idY < Img.m_Reso.y)
//	{
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y  * 0.5f + Img.m_Bias.y);
//		//cur rotation angle;
//		float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//		//current source position, assume the initial position is on the positive Y axis;
//		Ray2D ray;
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//		float2 minBox, maxBox, curImg, curDet, initDet, sincosAng;
//		float biasAng, tnear, tfar, weg, detID;
//		int flrID, ceilID;
//
//
//		minBox = make_float2(
//			MINO.x + idX * Img.m_Step.x,
//			MINO.y + idY * Img.m_Step.y);
//		maxBox = minBox + make_float2(Img.m_Step.x, Img.m_Step.y);
//		curImg = minBox + 0.5f * make_float2(Img.m_Step.x, Img.m_Step.y);
//		ray.d = normalize(curImg - ray.o);
//
//		//current detector element Cartesian coordinates;
//		curDet = ray.o + ray.d * FanGeo.m_S2D;
//		//original detector position;
//		initDet = rotation(curDet, cosT, -sinT);
//
//		//Bias angle;
//
//		sincosAng.x = initDet.x / FanGeo.m_S2D; //sin(ang)
//		sincosAng.y = (initDet.y - FanGeo.m_S2O) / (-FanGeo.m_S2D); //cos(ang)
//		biasAng = atan(sincosAng.x / sincosAng.y);
//
//		//judging the postive or negative angle direction;
//		if (angDir)
//		{
//			detID = biasAng / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//		}
//		else
//		{
//			detID = FanGeo.m_DetN - 1 - (biasAng / FanGeo.m_DetStp + FanGeo.m_DetCntIdx);
//		}
//		if (detID < 0)
//		{
//			dimg[idY * Img.m_Reso.x + idX] += 0;
//			return;
//		}
//		if (detID > FanGeo.m_DetN)
//		{
//			dimg[idY * Img.m_Reso.x + idX] += 0;
//			return;
//		}
//		tnear = 0;
//		tfar = 0;
//		flrID = detID;
//		ceilID = ceilf(detID);
//		intersectBox(ray, minBox, maxBox, &tnear, &tfar);
//
//		weg = (tfar - tnear); //no weighting;
//		if (!IS_ZERO(flrID - ceilID))
//		{
//			dimg[idY * Img.m_Reso.x + idX] += ((donp[flrID] * (ceilID - detID) + donp[ceilID] * (detID - flrID)) * weg * dmsk[idY * Img.m_Reso.x + idX]);
//		}
//		else
//		{
//			dimg[idY * Img.m_Reso.x + idX] += (donp[flrID] * weg * dmsk[idY * Img.m_Reso.x + idX]);
//		}
//	}
//}
////////////////////////////////////////////////////////////////////////////
//// Backprojection from one angle 
//// donp: projection data from one angle;
//// dimg: image;
//// dmsk: mask;
//// angIdx: current angle index;
//// angDir: index from positive or negative direction;
//// FanGeo: Fan Beam Geometry
//// Img: Image configuration
//// blk: 32,32
//// gid: configured
////////////////////////////////////////////////////////////////////////////
//void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, cuint angIdx, const bool angDir, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dimg, dmsk, angIdx, angDir, FanGeo, Img);
//}
//
//
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dimg, float* dmsk, cuint angIdx, const bool angDir, const FanEDGeo FanGeo, const Image Img)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	if (idX < Img.m_Reso.x && idY < Img.m_Reso.y)
//	{
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//		//cur rotation angle;
//		float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//current pixel position
//		float2 curPix = rotation(make_float2(MINO.x + Img.m_Step.x * idX, MINO.y + Img.m_Step.y * idY), cosT, -sinT); //��ת���λ��;
//		//�����Ӧ��detidx;
//
//		float2 minBox, maxBox, curImg;
//		float tnear(0), tfar(0), weg(0), detID(0);
//		int flrID(0), ceilID(0);
//
//		//���ﲻ���� angDir; ���Ǵ�С�����;
//		detID = FanGeo.m_S2D * curPix.x / (curPix.y - FanGeo.m_S2O) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//		if (detID < 0 || detID > FanGeo.m_DetN - 1)
//		{
//			dimg[idY * Img.m_Reso.x + idX] += 0;
//			return;
//		}
//		flrID = detID;
//		ceilID = ceilf(detID);
//
//		Ray2D ray;
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//		minBox = make_float2(
//			MINO.x + idX * Img.m_Step.x,
//			MINO.y + idY * Img.m_Step.y);
//		maxBox = minBox + make_float2(Img.m_Step.x, Img.m_Step.y);
//		curImg = minBox + 0.5f * make_float2(Img.m_Step.x, Img.m_Step.y);
//		ray.d = normalize(curImg - ray.o);
//		intersectBox(ray, minBox, maxBox, &tnear, &tfar); //�ཻ����;
//
//		weg = (tfar - tnear); //no weighting;
//		if (!IS_ZERO(flrID - ceilID))
//		{
//			dimg[idY * Img.m_Reso.x + idX] += ((donp[flrID] * (ceilID - detID) + donp[ceilID] * (detID - flrID)) * weg * dmsk[idY * Img.m_Reso.x + idX]);
//		}
//		else
//		{
//			dimg[idY * Img.m_Reso.x + idX] += (donp[flrID] * weg * dmsk[idY * Img.m_Reso.x + idX]);
//		}
//	}
//}
//void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, cuint angIdx, const bool angDir, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dimg, dmsk, angIdx, angDir, FanGeo, Img);
//}
//
//
//
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dimg, float* dmsk, const bool angDir, const FanEAGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)//(float* donp, float* dimg, const bool angDir, const bool FOV, const FanEAGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	if (idX < Img.m_Reso.x && idY < Img.m_Reso.y)
//	{
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//		//cur rotation angle;
//		unsigned int angIdx = 0;
//		float curAng, cosT, sinT;
//		Ray2D ray;
//		float2 minBox, maxBox, curImg, curDet, initDet, sincosAng;
//		float biasAng, tnear, tfar, weg, detID;
//		int flrID, ceilID;
//		float summ(0);
//		float masks(dmsk[idY * Img.m_Reso.x + idX]);
//		for (unsigned int prjIdx = 0; prjIdx != numPerSubSet; ++prjIdx)
//		{
//			angIdx = prjIdx * subSetNum + curSubSetIdx;
//			curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//			//current source position, assume the initial position is on the positive Y axis;
//			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//			minBox = make_float2(
//				MINO.x + idX * Img.m_Step.x,
//				MINO.y + idY * Img.m_Step.y);
//			maxBox = minBox + make_float2(Img.m_Step.x, Img.m_Step.y);
//			curImg = minBox + 0.5f * make_float2(Img.m_Step.x, Img.m_Step.y);
//			ray.d = normalize(curImg - ray.o);
//
//
//			//current detector element Cartesian coordinates;
//			curDet = ray.o + ray.d * FanGeo.m_S2D;
//			//original detector position;
//			initDet = rotation(curDet, cosT, -sinT);
//
//			//Bias angle;
//
//			sincosAng.x = initDet.x / FanGeo.m_S2D; //sin(ang)
//			sincosAng.y = (initDet.y - FanGeo.m_S2O) / (-FanGeo.m_S2D); //cos(ang)
//			biasAng = atan(sincosAng.x / sincosAng.y);
//
//			//judging the postive or negative angle direction;
//			if (angDir)
//			{
//				detID = biasAng / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			}
//			else
//			{
//				detID = FanGeo.m_DetN - 1 - (biasAng / FanGeo.m_DetStp + FanGeo.m_DetCntIdx);
//			}
//			if (detID < 0)
//			{
//				dimg[idY * Img.m_Reso.x + idX] += 0;
//				continue;
//			}
//			if (detID > FanGeo.m_DetN)
//			{
//				dimg[idY * Img.m_Reso.x + idX] += 0;
//				continue;
//			}
//			tnear = 0;
//			tfar = 0;
//			flrID = detID;
//			ceilID = ceilf(detID);
//			intersectBox(ray, minBox, maxBox, &tnear, &tfar);
//
//			weg = (tfar - tnear); //no weighting;
//			if (!IS_ZERO(flrID - ceilID))
//			{
//				summ += ((donp[prjIdx * FanGeo.m_DetN + flrID] * (ceilID - detID) + donp[prjIdx * FanGeo.m_DetN + ceilID] * (detID - flrID)) * weg);
//			}
//			else
//			{
//				summ += (donp[prjIdx * FanGeo.m_DetN + flrID] * weg);
//			}
//		}
//		dimg[idY * Img.m_Reso.x + idX] = (summ * masks);
//	}
//}
////////////////////////////////////////////////////////////////////////////
//// Backprojection from subset angles
//// donp: projection inside the subsets;
//// dimg: image;
//// dmsk: mask;
//// angDir: index direction from positive or negative;
//// FanGeo: Fan Beam Geometry
//// Img: Image configuration
//// numPerSubSet: How many projections are inside the subset;
//// subSetNum: How many subsets are distributed;
//// curSubSetIdx: current subset Index
//// blk: 32,32
//// gid: configured
////////////////////////////////////////////////////////////////////////////
//void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, const bool angDir, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dimg, dmsk, angDir, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx);
//}
//
//
//
//
//
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dimg, float* dmsk, const bool angDir, const FanEDGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)//(float* donp, float* dimg, const bool angDir, const bool FOV, const FanEDGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	if (idX < Img.m_Reso.x && idY < Img.m_Reso.y)
//	{
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//		//cur rotation angle;
//		unsigned int angIdx = 0;
//		float curAng(0), cosT(0), sinT(0), summ(0), mask(dmsk[idY * Img.m_Reso.x + idX]);
//		Ray2D ray;
//		float2 minBox, maxBox, curImg, curPix;
//		float tnear(0), tfar(0), weg(0), detID(0);
//		int flrID, ceilID;
//
//		for (unsigned int prjIdx = 0; prjIdx != numPerSubSet; ++prjIdx)
//		{
//			angIdx = prjIdx * subSetNum + curSubSetIdx;
//			curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//
//			//current pixel position
//			curPix = rotation(make_float2(MINO.x + Img.m_Step.x * idX, MINO.y + Img.m_Step.y * idY), cosT, -sinT); //��ת���λ��;
//
//			//�����Ӧ��detidx;			//���ﲻ���� angDir; ���Ǵ�С�����;
//			detID = FanGeo.m_S2D * curPix.x / (curPix.y - FanGeo.m_S2O) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			if (detID < 0 || detID > FanGeo.m_DetN - 1)
//			{
//				dimg[idY * Img.m_Reso.x + idX] += 0;
//				continue;
//			}
//			flrID = detID;
//			ceilID = ceilf(detID);
//
//			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//			minBox = make_float2(
//				MINO.x + idX * Img.m_Step.x,
//				MINO.y + idY * Img.m_Step.y);
//			maxBox = minBox + make_float2(Img.m_Step.x, Img.m_Step.y);
//			curImg = minBox + 0.5f * make_float2(Img.m_Step.x, Img.m_Step.y);
//			ray.d = normalize(curImg - ray.o);
//
//			tnear = 0;
//			tfar = 0;
//			intersectBox(ray, minBox, maxBox, &tnear, &tfar); //�ཻ����;
//
//			weg = (tfar - tnear); //no weighting;
//			if (!IS_ZERO(flrID - ceilID))
//			{
//				summ += ((donp[prjIdx * FanGeo.m_DetN + flrID] * (ceilID - detID) + donp[prjIdx * FanGeo.m_DetN + ceilID] * (detID - flrID)) * weg);
//			}
//			else
//			{
//				summ += (donp[prjIdx * FanGeo.m_DetN + flrID] * weg);
//			}
//		}
//		dimg[idY * Img.m_Reso.x + idX] = (summ * mask);
//	}
//}
//void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, const bool angDir, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dimg, dmsk, angDir, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx);
//}
//
//
//
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dimg, float* dmsk, const FanEAGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	if (idX < Img.m_Reso.x && idY < Img.m_Reso.y)
//	{
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//		//cur rotation angle;
//		unsigned int angIdx = 0;
//		float curAng, cosT, sinT;
//		Ray2D ray;
//		float2 minBox, maxBox, curImg, curDet, initDet, sincosAng;
//		float biasAng, tnear, tfar, weg, detID;
//		int flrID, ceilID;
//		float summ(0);
//		float mask(dmsk[idY * Img.m_Reso.x + idX]);
//		for (unsigned int prjIdx = 0; prjIdx != numPerSubSet; ++prjIdx)
//		{
//			angIdx = prjIdx * subSetNum + curSubSetIdx;
//			curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//			//current source position, assume the initial position is on the positive Y axis;
//			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//			minBox = make_float2(
//				MINO.x + idX * Img.m_Step.x,
//				MINO.y + idY * Img.m_Step.y);
//			maxBox = minBox + make_float2(Img.m_Step.x, Img.m_Step.y);
//			curImg = minBox + 0.5f * make_float2(Img.m_Step.x, Img.m_Step.y);
//			ray.d = normalize(curImg - ray.o);
//
//
//			//current detector element Cartesian coordinates;
//			curDet = ray.o + ray.d * FanGeo.m_S2D;
//			//original detector position;
//			initDet = rotation(curDet, cosT, -sinT);
//
//			//Bias angle;
//
//			sincosAng.x = initDet.x / FanGeo.m_S2D; //sin(ang)
//			sincosAng.y = (initDet.y - FanGeo.m_S2O) / (-FanGeo.m_S2D); //cos(ang)
//			biasAng = atan(sincosAng.x / sincosAng.y);
//
//			detID = biasAng / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//			if (detID < 0)
//			{
//				dimg[idY * Img.m_Reso.x + idX] += 0;
//				return;
//			}
//			if (detID > FanGeo.m_DetN)
//			{
//				dimg[idY * Img.m_Reso.x + idX] += 0;
//				return;
//			}
//			tnear = 0;
//			tfar = 0;
//			flrID = detID;
//			ceilID = ceilf(detID);
//			intersectBox(ray, minBox, maxBox, &tnear, &tfar);
//
//			weg = (tfar - tnear); //no weighting;
//			if (!IS_ZERO(flrID - ceilID))
//			{
//				summ += ((donp[prjIdx * FanGeo.m_DetN + flrID] * (ceilID - detID) + donp[prjIdx * FanGeo.m_DetN + ceilID] * (detID - flrID)) * weg);
//			}
//			else
//			{
//				summ += (donp[prjIdx * FanGeo.m_DetN + flrID] * weg);
//			}
//		}
//		dimg[idY * Img.m_Reso.x + idX] = (summ * mask);
//	}
//}
////////////////////////////////////////////////////////////////////////////
//// Backprojection from subset angles
//// donp: projection inside the subsets;
//// dimg: image;
//// dmsk: mask;
//// FanGeo: Fan Beam Geometry
//// Img: Image configuration
//// numPerSubSet: How many projections are inside the subset;
//// subSetNum: How many subsets are distributed;
//// curSubSetIdx: current subset Index
//// blk: 32,32
//// gid: configured
////////////////////////////////////////////////////////////////////////////
//void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dimg, dmsk, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx);
//}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dimg, float* dmask, const FanEDGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	if (idX < Img.m_Reso.x && idY < Img.m_Reso.y)
//	{
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//		//cur rotation angle;
//		unsigned int angIdx = 0;
//		float curAng(0), cosT(0), sinT(0), summ(0), mask(dmask[idY * Img.m_Reso.x + idX]);
//		Ray2D ray;
//		float2 minBox, maxBox, curImg, curPix;
//		float tnear(0), tfar(0), weg(0), detID(0);
//		int flrID, ceilID;
//
//		for (unsigned int prjIdx = 0; prjIdx != numPerSubSet; ++prjIdx)
//		{
//			angIdx = prjIdx * subSetNum + curSubSetIdx;
//			curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//
//			//current pixel position
//			curPix = rotation(make_float2(MINO.x + Img.m_Step.x * idX, MINO.y + Img.m_Step.y * idY), cosT, -sinT); //��ת���λ��;
//
//			//�����Ӧ��detidx;			//���ﲻ���� angDir; ���Ǵ�С�����;
//			detID = FanGeo.m_S2D * curPix.x / (curPix.y - FanGeo.m_S2O) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			if (detID < 0 || detID > FanGeo.m_DetN - 1)
//			{
//				dimg[idY * Img.m_Reso.x + idX] += 0;
//				return;
//			}
//			flrID = detID;
//			ceilID = ceilf(detID);
//
//			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//			minBox = make_float2(
//				MINO.x + idX * Img.m_Step.x,
//				MINO.y + idY * Img.m_Step.y);
//			maxBox = minBox + make_float2(Img.m_Step.x, Img.m_Step.y);
//			curImg = minBox + 0.5f * make_float2(Img.m_Step.x, Img.m_Step.y);
//			ray.d = normalize(curImg - ray.o);
//
//			tnear = 0;
//			tfar = 0;
//			intersectBox(ray, minBox, maxBox, &tnear, &tfar); //�ཻ����;
//
//			weg = (tfar - tnear); //no weighting;
//			if (!IS_ZERO(flrID - ceilID))
//			{
//				summ += ((donp[prjIdx * FanGeo.m_DetN + flrID] * (ceilID - detID) + donp[prjIdx * FanGeo.m_DetN + ceilID] * (detID - flrID)) * weg);
//			}
//			else
//			{
//				summ += donp[prjIdx * FanGeo.m_DetN + flrID] * weg;
//			}
//		}
//		dimg[idY * Img.m_Reso.x + idX] = (summ * mask);
//	}
//}
//void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dimg, dmsk, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx);
//}
//
//
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dimg, float* dmsk, const FanEAGeo FanGeo, const Image Img)//; (float* donp, float* dimg, float* dmsk, const FanEAGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	if (idX < Img.m_Reso.x && idY < Img.m_Reso.y)
//	{
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//		//cur rotation angle;
//
//		float curAng, cosT, sinT;
//		Ray2D ray;
//		float2 minBox, maxBox, curImg, curDet, initDet, sincosAng;
//		float biasAng, tnear, tfar, weg, detID;
//		int flrID, ceilID;
//		float summ(0);
//		float msk = dmsk[idY * Img.m_Reso.x + idX];
//		for (unsigned int angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
//		{
//			curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//			//current source position, assume the initial position is on the positive Y axis;
//			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//			minBox = make_float2(
//				MINO.x + idX * Img.m_Step.x,
//				MINO.y + idY * Img.m_Step.y);
//			maxBox = minBox + make_float2(Img.m_Step.x, Img.m_Step.y);
//			curImg = minBox + 0.5f * make_float2(Img.m_Step.x, Img.m_Step.y);
//			ray.d = normalize(curImg - ray.o);
//
//
//			//current detector element Cartesian coordinates;
//			curDet = ray.o + ray.d * FanGeo.m_S2D;
//			//original detector position;
//			initDet = rotation(curDet, cosT, -sinT);
//
//			//Bias angle;
//
//			sincosAng.x = initDet.x / FanGeo.m_S2D; //sin(ang)
//			sincosAng.y = (initDet.y - FanGeo.m_S2O) / (-FanGeo.m_S2D); //cos(ang)
//			biasAng = atan(sincosAng.x / sincosAng.y);
//
//			detID = biasAng / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			if (detID < 0)
//			{
//				dimg[idY * Img.m_Reso.x + idX] += 0;
//				return;
//			}
//			if (detID > FanGeo.m_DetN)
//			{
//				dimg[idY * Img.m_Reso.x + idX] += 0;
//				return;
//			}
//			tnear = 0;
//			tfar = 0;
//			flrID = detID;
//			ceilID = ceilf(detID);
//			intersectBox(ray, minBox, maxBox, &tnear, &tfar);
//
//			weg = (tfar - tnear); //no weighting;
//			if (!IS_ZERO(flrID - ceilID))
//			{
//				summ += ((donp[angIdx * FanGeo.m_DetN + flrID] * (ceilID - detID) + donp[angIdx * FanGeo.m_DetN + ceilID] * (detID - flrID)) * weg);
//			}
//			else
//			{
//				summ += (donp[angIdx * FanGeo.m_DetN + flrID] * weg);
//			}
//		}
//		dimg[idY * Img.m_Reso.x + idX] = summ * msk;
//	}
//}
////////////////////////////////////////////////////////////////////////////
//// Backprojection from all angles
//// donp: projection inside the subsets;
//// dimg: image;
//// dmsk: mask;
//// FanGeo: Fan Beam Geometry
//// Img: Image configuration
//// blk: 32,32
//// gid: configured
////////////////////////////////////////////////////////////////////////////
//void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dimg, dmsk, FanGeo, Img);
//}
//
//
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//__global__ void _bakproj_PIXEL_Ker(float* donp, float* dimg, float* dmsk, const FanEDGeo FanGeo, const Image Img)//;(float* donp, float* dimg, const FanEDGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	if (idX < Img.m_Reso.x && idY < Img.m_Reso.y)
//	{
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//		//cur rotation angle;
//		//unsigned int angIdx = 0;
//		float curAng(0), cosT(0), sinT(0), summ(0), mask(dmsk[idY * Img.m_Reso.x + idX]);
//		Ray2D ray;
//		float2 minBox, maxBox, curImg, curPix;
//		float tnear(0), tfar(0), weg(0), detID(0);
//		int flrID, ceilID;
//
//		for (unsigned int angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
//		{
//			//angIdx = prjIdx * subSetNum + curSubSetIdx;
//			curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//
//			//current pixel position
//			curPix = rotation(make_float2(MINO.x + Img.m_Step.x * idX, MINO.y + Img.m_Step.y * idY), cosT, -sinT); //��ת���λ��;
//
//			//�����Ӧ��detidx;			//���ﲻ���� angDir; ���Ǵ�С�����;
//			detID = FanGeo.m_S2D * curPix.x / (curPix.y - FanGeo.m_S2O) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			if (detID < 0 || detID > FanGeo.m_DetN - 1)
//			{
//				dimg[idY * Img.m_Reso.x + idX] += 0;
//				return;
//			}
//			flrID = detID;
//			ceilID = ceilf(detID);
//
//			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//			minBox = make_float2(
//				MINO.x + idX * Img.m_Step.x,
//				MINO.y + idY * Img.m_Step.y);
//			maxBox = minBox + make_float2(Img.m_Step.x, Img.m_Step.y);
//			curImg = minBox + 0.5f * make_float2(Img.m_Step.x, Img.m_Step.y);
//			ray.d = normalize(curImg - ray.o);
//
//			tnear = 0;
//			tfar = 0;
//			intersectBox(ray, minBox, maxBox, &tnear, &tfar); //�ཻ����;
//
//			weg = (tfar - tnear); //no weighting;
//			if (!IS_ZERO(flrID - ceilID))
//			{
//				summ += ((donp[angIdx * FanGeo.m_DetN + flrID] * (ceilID - detID) + donp[angIdx * FanGeo.m_DetN + ceilID] * (detID - flrID)) * weg);
//			}
//			else
//			{
//				summ += donp[angIdx * FanGeo.m_DetN + flrID] * weg;
//			}
//		}
//		dimg[idY * Img.m_Reso.x + idX] = (summ * mask);
//	}
//}
//void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker << <gid, blk >> >(donp, dimg, dmsk, FanGeo, Img);
//}
//
//
//__global__ void bakproj_PIXEL_Ker(float* donp, float* dimg, const FanEAGeo FanGeo, const Image Img, cuint sliceNum)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	if (idX < Img.m_Reso.x && idY < Img.m_Reso.y)
//	{
//		float2 MINO = make_float2(-Img.m_Size.x / 2.0f + Img.m_Bias.x, -Img.m_Size.y / 2.0f + Img.m_Bias.y);
//		//cur rotation angle;
//		float curAng;// = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT;//
//		float sinT;
//		//current source position, assume the initial position is on the positive Y axis;
//		Ray2D ray;
//
//		float2 curImg, initDet;
//		float biasAng, tnear, tfar, weg, detID;
//		int flrID, ceilID;
//		unsigned int angIdx = 0;
//		unsigned int sliceID;
//
//
//		float2 minBox = make_float2(
//			MINO.x + idX * Img.m_Step.x,
//			MINO.y + idY * Img.m_Step.y);
//		float2 maxBox = minBox + make_float2(Img.m_Step.x, Img.m_Step.y);
//		curImg = minBox + 0.5f * make_float2(Img.m_Step.x, Img.m_Step.y);
//		float upd[DEMO10_SLICESNUM];
//		for (sliceID = 0; sliceID != sliceNum; ++sliceID)
//		{
//			upd[sliceID] = 0;
//		}
//		//�������нǶ�;
//		for (angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
//		{
//			curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//
//			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//			ray.d = normalize(curImg - ray.o);
//
//			//current detector element Cartesian coordinates;
//
//			//original detector position;
//			initDet = rotation(ray.o + ray.d * FanGeo.m_S2D, cosT, -sinT);
//
//			//Bias angle;
//			biasAng = atan(initDet.x / (FanGeo.m_S2O - initDet.y));
//
//			//judging the postive or negative angle direction;
//			detID = biasAng / FanGeo.m_DetStp + FanGeo.m_DetCntIdx - 0.5f;
//			if (detID < 0)
//			{
//				dimg[(sliceID * Img.m_Reso.y + idY) * Img.m_Reso.x + idX] += 0;
//				continue;
//			}
//			if (detID > FanGeo.m_DetN)
//			{
//				dimg[(sliceID * Img.m_Reso.y + idY) * Img.m_Reso.x + idX] += 0;
//				continue;
//			}
//			tnear = 0;
//			tfar = 0;
//			flrID = detID;
//			ceilID = ceilf(detID);
//			intersectBox(ray, minBox, maxBox, &tnear, &tfar);
//
//			weg = (tfar - tnear); //no weighting;
//
//			//ͬʱ����64��
//			for (sliceID = 0; sliceID != sliceNum; ++sliceID)
//			{
//				if (!IS_ZERO(flrID - ceilID))
//				{
//					upd[sliceID] +=
//						((donp[(sliceID * FanGeo.m_ViwN + angIdx) * FanGeo.m_DetN + flrID] * (ceilID - detID)
//						+ donp[(sliceID * FanGeo.m_ViwN + angIdx) * FanGeo.m_DetN + ceilID] * (detID - flrID)) * weg);
//				}
//				else
//				{
//					upd[sliceID] += donp[(sliceID * FanGeo.m_ViwN + angIdx) * FanGeo.m_DetN + flrID] * weg;
//				}
//			}//End 64 slices;
//		}//End all angles;
//
//		for (sliceID = 0; sliceID != sliceNum; ++sliceID)
//		{
//			dimg[(sliceID * Img.m_Reso.y + idY) * Img.m_Reso.x + idX] += upd[sliceID];
//		}
//
//	}
//}
////��ϵ�к����ʾͬ�������configuration���ڲ�ͬ����ͬʱ��ͶӰ�ĺ���Ӧ����patientRec �Ǹ�64����ع�
//void bakproj_PIXEL(float* donp, float* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid, cuint sliceNum)
//{
//	bakproj_PIXEL_Ker << <gid, blk >> >(donp, dimg, FanGeo, Img, sliceNum);
//}
//
//
//
//__global__ void bakproj_PIXEL_Ker(float* donp, float* dimg, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const FanEAGeo FanGeo, const Image Img, cuint sliceNum)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	if (idX < Img.m_Reso.x && idY < Img.m_Reso.y)
//	{
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//		//cur rotation angle;
//		float curAng;// = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT;//
//		float sinT;
//		//current source position, assume the initial position is on the positive Y axis;
//		Ray2D ray;
//
//		float2 minBox, maxBox, curImg, curDet, initDet, sincosAng;
//		float biasAng, tnear, tfar, weg, detID;
//		int flrID, ceilID;
//		unsigned int angIdx = 0;
//		unsigned int sliceID;
//		unsigned int prjIdx = 0;
//
//		minBox = make_float2(
//			MINO.x + idX * Img.m_Step.x,
//			MINO.y + idY * Img.m_Step.y);
//		maxBox = minBox + make_float2(Img.m_Step.x, Img.m_Step.y);
//		curImg = minBox + 0.5f * make_float2(Img.m_Step.x, Img.m_Step.y);
//
//
//		float upd[DEMO11_SLICESNUM];  //We have to change here if necessary;
//
//		for (sliceID = 0; sliceID != sliceNum; ++sliceID)
//		{
//			upd[sliceID] = 0;
//		}
//		//�������нǶ�;
//		for (prjIdx = 0; prjIdx != numPerSubSet; ++prjIdx)
//		{
//			angIdx = prjIdx * subSetNum + curSubSetIdx;
//			curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//
//			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//			ray.d = normalize(curImg - ray.o);
//
//			//current detector element Cartesian coordinates;
//			curDet = ray.o + ray.d * FanGeo.m_S2D;
//			//original detector position;
//			initDet = rotation(curDet, cosT, -sinT);
//
//			//Bias angle;
//			//sincosAng;
//			sincosAng.x = initDet.x / FanGeo.m_S2D; //sin(ang)
//			sincosAng.y = (initDet.y - FanGeo.m_S2O) / (-FanGeo.m_S2D); //cos(ang)
//			biasAng = atan(sincosAng.x / sincosAng.y);
//
//			//judging the postive or negative angle direction;
//			detID = biasAng / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			if (detID < 0)
//			{
//				dimg[(sliceID * Img.m_Reso.y + idY) * Img.m_Reso.x + idX] += 0;
//				return;
//			}
//			if (detID > FanGeo.m_DetN)
//			{
//				dimg[(sliceID * Img.m_Reso.y + idY) * Img.m_Reso.x + idX] += 0;
//				return;
//			}
//			tnear = 0;
//			tfar = 0;
//			flrID = detID;
//			ceilID = ceilf(detID);
//			intersectBox(ray, minBox, maxBox, &tnear, &tfar);
//
//			weg = (tfar - tnear); //no weighting;
//
//			//ͬʱ����22��
//			for (sliceID = 0; sliceID != sliceNum; ++sliceID)
//			{
//				if (!IS_ZERO(flrID - ceilID))
//				{
//					upd[sliceID] +=
//						((donp[(sliceID * numPerSubSet + prjIdx) * FanGeo.m_DetN + flrID] * (ceilID - detID)
//						+ donp[(sliceID * numPerSubSet + prjIdx) * FanGeo.m_DetN + ceilID] * (detID - flrID)) * weg);
//				}
//				else
//				{
//					upd[sliceID] += donp[(sliceID * numPerSubSet + prjIdx) * FanGeo.m_DetN + flrID] * weg;
//				}
//			}//End 22 slices;
//		}//End all angles;
//
//		for (sliceID = 0; sliceID != sliceNum; ++sliceID)
//		{
//			dimg[(sliceID * Img.m_Reso.y + idY) * Img.m_Reso.x + idX] += upd[sliceID];
//		}
//
//	}
//}
//void bakproj_PIXEL(float* donp, float* dimg, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubIdx, const dim3& blk, const dim3& gid, cuint sliceNum, const cudaStream_t& streams)
//{
//	//float* donp, float* dimg, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const FanEAGeo FanGeo, const Image Img,cuint sliceNum
//	bakproj_PIXEL_Ker << <gid, blk, 0, streams >> >(donp, dimg, numPerSubSet, subSetNum, curSubIdx, FanGeo, Img, sliceNum);
//}
//
//
//
//template<typename T>
//__global__ void _bakproj_AIM_ker(T* dprj, T* dimg, const FanEAGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, cuint sliceNum)
//{
//	cuint xi = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint yi = threadIdx.y + blockIdx.y * blockDim.y;
//
//	if (xi < Img.m_Reso.x&&yi < Img.m_Reso.y)
//	{
//		const T cntImgX = static_cast<T>(Img.m_Reso.x - 1.0) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
//		const T cntImgY = static_cast<T>(Img.m_Reso.y - 1.0) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
//		const T area = Img.m_Step.x * Img.m_Step.y;
//		T sour[2], cosT, sinT;
//		T grid[4][3];
//		grid[0][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
//		grid[0][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
//
//		grid[1][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
//		grid[1][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
//
//		grid[2][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
//		grid[2][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
//
//		grid[3][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
//		grid[3][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
//		T det[4][3];
//		T tmp[2];
//		T beta[4];
//		int minDetIdx, maxDetIdx;
//		T summ[22];
//		T SVA[3], SVB[3];
//		unsigned int angIdx;
//		unsigned int levelIdx;
//		unsigned int prjIdx(0);
//		T ttpp[22]; // used in 22 slices
//		for (levelIdx = 0; levelIdx != sliceNum; ++levelIdx)
//		{
//			summ[levelIdx] = 0;
//		}
//		for (prjIdx = 0; prjIdx < numPerSubSet; ++angIdx)
//		{
//			angIdx = prjIdx * subSetNum + curSubSetIdx;
//			cosT = cos(FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp);
//			sinT = sin(FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp);
//
//			sour[0] = -FanGeo.m_S2O * sinT;
//			sour[1] = FanGeo.m_S2O * cosT;
//
//			det[0][0] = grid[0][0] - sour[0]; //used as the direction from source to the corner of the pixel
//			det[0][1] = grid[0][1] - sour[1];
//			det[1][0] = grid[1][0] - sour[0];
//			det[1][1] = grid[1][1] - sour[1];
//			det[2][0] = grid[2][0] - sour[0];
//			det[2][1] = grid[2][1] - sour[1];
//			det[3][0] = grid[3][0] - sour[0];
//			det[3][1] = grid[3][1] - sour[1];
//
//
//			tmp[1] = sqrt(det[0][0] * det[0][0] + det[0][1] * det[0][1]);
//			det[0][0] /= tmp[1];
//			det[0][1] /= tmp[1];
//			tmp[1] = sqrt(det[1][0] * det[1][0] + det[1][1] * det[1][1]);
//			det[1][0] /= tmp[1];
//			det[1][1] /= tmp[1];
//			tmp[1] = sqrt(det[2][0] * det[2][0] + det[2][1] * det[2][1]);
//			det[2][0] /= tmp[1];
//			det[2][1] /= tmp[1];
//			tmp[1] = sqrt(det[3][0] * det[3][0] + det[3][1] * det[3][1]);
//			det[3][0] /= tmp[1];
//			det[3][1] /= tmp[1];
//
//			det[0][0] = sour[0] + det[0][0] * FanGeo.m_S2D;  //det used as the detector coordinates
//			det[0][1] = sour[1] + det[0][1] * FanGeo.m_S2D;
//			det[1][0] = sour[0] + det[1][0] * FanGeo.m_S2D;
//			det[1][1] = sour[1] + det[1][1] * FanGeo.m_S2D;
//			det[2][0] = sour[0] + det[2][0] * FanGeo.m_S2D;
//			det[2][1] = sour[1] + det[2][1] * FanGeo.m_S2D;
//			det[3][0] = sour[0] + det[3][0] * FanGeo.m_S2D;
//			det[3][1] = sour[1] + det[3][1] * FanGeo.m_S2D;
//
//			tmp[0] = det[0][0] * cosT + det[0][1] * sinT;
//			tmp[1] = -det[0][0] * sinT + det[0][1] * cosT;
//			det[0][0] = tmp[0];
//			det[0][1] = tmp[1];
//			tmp[0] = det[1][0] * cosT + det[1][1] * sinT;
//			tmp[1] = -det[1][0] * sinT + det[1][1] * cosT;
//			det[1][0] = tmp[0];
//			det[1][1] = tmp[1];
//			tmp[0] = det[2][0] * cosT + det[2][1] * sinT;
//			tmp[1] = -det[2][0] * sinT + det[2][1] * cosT;
//			det[2][0] = tmp[0];
//			det[2][1] = tmp[1];
//			tmp[0] = det[3][0] * cosT + det[3][1] * sinT;
//			tmp[1] = -det[3][0] * sinT + det[3][1] * cosT;
//			det[3][0] = tmp[0];
//			det[3][1] = tmp[1];
//
//
//			beta[0] = atan(-det[0][0] / (det[0][1] - FanGeo.m_S2O));
//			beta[1] = atan(-det[1][0] / (det[1][1] - FanGeo.m_S2O));
//			beta[2] = atan(-det[2][0] / (det[2][1] - FanGeo.m_S2O));
//			beta[3] = atan(-det[3][0] / (det[3][1] - FanGeo.m_S2O));
//
//			grid[0][2] = beta[0] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			grid[1][2] = beta[1] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			grid[2][2] = beta[2] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			grid[3][2] = beta[3] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//			SortProjection<T>(grid);
//			minDetIdx = static_cast<int>(grid[0][2]);
//			maxDetIdx = (int(ceil(grid[3][2])));
//			minDetIdx = (minDetIdx < 0) ? 0 : minDetIdx;
//			maxDetIdx = (maxDetIdx > static_cast<int>(FanGeo.m_DetN)) ? static_cast<int>(FanGeo.m_DetN) : maxDetIdx;
//
//			tmp[1] = (hypot((xi - cntImgX)*Img.m_Step.x - sour[0], (yi - cntImgY)*Img.m_Step.y - sour[1]));
//			//T summ = 0;
//			//T pangle;
//			//T direct[2], legth, SVA[3], SVB[3], coef;
//
//
//			//�����߱߼���;
//			beta[3] = (minDetIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp; //���С�߱�ƫ�ƽ� // Beta[3] used as the pangle
//			beta[0] = sin(beta[3]) * FanGeo.m_S2D; // beta used as the initial X,Y and current X,Y with rotation
//			beta[1] = -cos(beta[3]) * FanGeo.m_S2D + FanGeo.m_S2O;
//
//			beta[2] = beta[0] * cosT - beta[1] * sinT;
//			beta[3] = beta[0] * sinT + beta[1] * cosT;
//			beta[0] = beta[2] - sour[0]; //Beta used as the directions
//			beta[1] = beta[3] - sour[1];
//			beta[2] = sqrt(beta[0] * beta[0] + beta[1] * beta[1]);//Beta[3] used as the length of the image
//			SVA[0] = beta[0] / beta[2];
//			SVA[1] = beta[1] / beta[2];
//			SVA[2] = static_cast<T>(minDetIdx);
//			tmp[0] = 0;
//			for (levelIdx = 0; levelIdx < sliceNum; levelIdx++)
//			{
//				ttpp[levelIdx] = 0;
//			}
//			for (; minDetIdx < maxDetIdx; ++minDetIdx)
//			{
//				//�����߱߼���;
//
//				beta[3] = (minDetIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp + FanGeo.m_DetStp; //��߱�ƫ�ƽ�
//				beta[0] = sin(beta[3]) * FanGeo.m_S2D;
//				beta[1] = -cos(beta[3]) * FanGeo.m_S2D + FanGeo.m_S2O;
//				beta[2] = beta[0] * cosT - beta[1] * sinT;
//				beta[3] = beta[0] * sinT + beta[1] * cosT;
//				beta[0] = beta[2] - sour[0];
//				beta[1] = beta[3] - sour[1];
//				beta[2] = sqrt(beta[0] * beta[0] + beta[1] * beta[1]);
//				SVB[0] = beta[0] / beta[2];
//				SVB[1] = beta[1] / beta[2];
//				SVB[2] = static_cast<T>(minDetIdx)+1;
//
//				//Compute the weighting coefficient for a special projection data
//				beta[3] = ComputeCoefficient<T>(grid, SVA, SVB, sour, area); // beta[3] used as the coef
//				beta[3] /= (tmp[1] * fabs(FanGeo.m_DetStp));
//				//tmp[0] += dprj[prjIdx* FanGeo.m_DetN + minDetIdx] * beta[3];
//				for (levelIdx = 0; levelIdx < sliceNum; levelIdx++)
//				{
//					ttpp[levelIdx] += dprj[(levelIdx * numPerSubSet + prjIdx) * FanGeo.m_DetN + minDetIdx] * beta[3];
//				}
//				SVA[0] = SVB[0]; //��߱߱��С�߱�
//				SVA[1] = SVB[1];
//				SVA[2] = SVB[2];
//			}
//			for (levelIdx = 0; levelIdx != sliceNum; ++levelIdx)
//			{
//				summ[levelIdx] += ttpp[levelIdx];
//			}
//		}
//		for (levelIdx = 0; levelIdx != sliceNum; ++levelIdx)
//		{
//			dimg[(levelIdx * Img.m_Reso.y + yi) * Img.m_Reso.x + xi] = summ[levelIdx];
//		}
//	}
//}
//template<typename T>
//void bakproj_AIM_temp(T* dprj, T* dimg, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubIdx, const dim3& blk, const dim3& gid, cuint sliceNum, const cudaStream_t& streams)
//{
//	_bakproj_AIM_ker<T> << <gid, blk, 0, streams >> >(dprj, dimg, FanGeo, Img, numPerSubSet, subSetNum, curSubIdx, sliceNum);
//}
//
//void bakproj_AIM(float* dprj, float* dimg, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubIdx, const dim3& blk, const dim3& gid, cuint sliceNum, const cudaStream_t& streams)
//{
//	bakproj_AIM_temp<float>(dprj, dimg, FanGeo, Img, numPerSubSet, subSetNum, curSubIdx, blk, gid, sliceNum, streams);
//}
//void bakproj_AIM(double* dprj, double* dimg, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubIdx, const dim3& blk, const dim3& gid, cuint sliceNum, const cudaStream_t& streams)
//{
//	bakproj_AIM_temp<double>(dprj, dimg, FanGeo, Img, numPerSubSet, subSetNum, curSubIdx, blk, gid, sliceNum, streams);
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
//
//
//
//
//
//__global__ void _bakproj_BOXED_ker(float* correction, float* volume, float* dmask, const unsigned int angIdx, const ConeEDGeo ConeGeo, const Volume Vol)
//{
//	/*uint3 idx;
//	idx.x = threadIdx.x + blockIdx.x * blockDim.x;
//	idx.y = threadIdx.y + blockIdx.y * blockDim.y;
//	idx.z = threadIdx.z + blockIdx.z * blockDim.z;*/
//	//cuint index = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint index = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
//	if (index >= Vol.m_Reso.x * Vol.m_Reso.y * Vol.m_Reso.z) return; //if (idx.x < Vol.m_Reso.x && idx.y < Vol.m_Reso.y && idx.z < Vol.m_Reso.z)
//
//	const float cosT = cosf(ConeGeo.m_ViwBeg + (float)angIdx * ConeGeo.m_ViwStp);
//	const float sinT = sinf(ConeGeo.m_ViwBeg + (float)angIdx * ConeGeo.m_ViwStp);
//	float2 MIND = make_float2(
//		-ConeGeo.m_DetCntIdx.x * ConeGeo.m_DetStp.x,
//		-ConeGeo.m_DetCntIdx.y * ConeGeo.m_DetStp.y);
//
//	float M[10];
//	M[0] = (-(ConeGeo.m_S2D) * cosT / ConeGeo.m_DetStp.x + MIND.x * sinT / ConeGeo.m_DetStp.x);
//	M[1] = (-(ConeGeo.m_S2D) * sinT / ConeGeo.m_DetStp.x - MIND.x * cosT / ConeGeo.m_DetStp.x);
//	M[2] = MIND.x * ConeGeo.m_S2O / ConeGeo.m_DetStp.x;
//	M[3] = MIND.y / ConeGeo.m_DetStp.y * sinT;
//	M[4] = (-MIND.y / ConeGeo.m_DetStp.y) * cosT;
//	M[5] = (-(ConeGeo.m_S2D) / ConeGeo.m_DetStp.y);
//	M[6] = MIND.y * ConeGeo.m_S2O / ConeGeo.m_DetStp.y;
//	M[7] = -sinT;
//	M[8] = cosT;
//	M[9] = -ConeGeo.m_S2O;
//
//	/*__shared__ float M[10];
//	M[1] = (-(__S2O__ + __O2D__) * sinT / __DETSTPL__ - __MINDL__ * cosT / __DETSTPL__);
//	M[2] = __MINDL__ * __S2O__ / __DETSTPL__;
//	M[3] = __MINDH__ / __DETSTPH__ * sinT;
//	M[4] = (-__MINDH__ / __DETSTPH__) * cosT;
//	M[5] = (-(__S2O__ + __O2D__) / __DETSTPH__);
//	M[6] = __MINDH__ * __S2O__ / __DETSTPH__;
//	M[7] = -sinT;
//	M[8] = cosT;
//	M[9] = -__S2O__;*/
//	//__syncthreads();
//
//
//	uint3 idx;
//	idx.z = index / (static_cast<int>(Vol.m_Reso.x) * static_cast<int>(Vol.m_Reso.y));
//	idx.y = (index - idx.z * (static_cast<int>(Vol.m_Reso.x) * static_cast<int>(Vol.m_Reso.y))) / static_cast<int>(Vol.m_Reso.x);
//	idx.x = index - idx.z * (static_cast<int>(Vol.m_Reso.x) * static_cast<int>(Vol.m_Reso.y)) - idx.y * static_cast<int>(Vol.m_Reso.x);
//
//	float minX = 10000, minZ = 10000, maxX = -10000, maxZ = -10000;
//
//	float3 projection, origPos, p, curP;
//	//int i(0), j(0), k(0);
//	//float tt;
//	Ray eyeRay;
//	eyeRay.o = make_float3(-ConeGeo.m_S2O * sinT, ConeGeo.m_S2O * cosT, 0);
//
//	float tnear = 0.0f;
//	float tfar = 0.0f;
//	float3 boxMin, boxMax;
//	float3 MINO = make_float3(
//		-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//		-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//		-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//
//	boxMin.x = MINO.x + idx.x * Vol.m_Step.x;
//	boxMin.y = MINO.y + idx.y * Vol.m_Step.y;
//	boxMin.z = MINO.z + idx.z * Vol.m_Step.z;
//	boxMax.x = boxMin.x + Vol.m_Step.x;
//	boxMax.y = boxMin.y + Vol.m_Step.y;
//	boxMax.z = boxMin.z + Vol.m_Step.z;
//
//
//	origPos.x = MINO.x + idx.x * Vol.m_Step.x;
//	origPos.y = MINO.y + idx.y * Vol.m_Step.y;
//	origPos.z = MINO.z + idx.z * Vol.m_Step.z;
//	//Trans from the world coordinate to detector index coordinate;
//	projection.x = M[0] * origPos.x + M[1] * origPos.y + M[2];
//	projection.y = M[3] * origPos.x + M[4] * origPos.y + M[5] * origPos.z + M[6];
//	projection.z = M[7] * origPos.x + M[8] * origPos.y + M[9];
//	projection.x = projection.x / projection.z;
//	projection.y = projection.y / projection.z;
//	projection.z = projection.z / projection.z;
//
//	//Obtain the bounding box in the projection plane;
//	minX = MY_MIN<float>(minX, projection.x);
//	maxX = MY_MAX<float>(maxX, projection.x);
//	minZ = MY_MIN<float>(minZ, projection.y);
//	maxZ = MY_MAX<float>(maxZ, projection.y);
//
//
//	origPos.x += Vol.m_Step.x;
//	//Trans from the world coordinate to detector index coordinate;
//	projection.x = M[0] * origPos.x + M[1] * origPos.y + M[2];
//	projection.y = M[3] * origPos.x + M[4] * origPos.y + M[5] * origPos.z + M[6];
//	projection.z = M[7] * origPos.x + M[8] * origPos.y + M[9];
//	projection.x = projection.x / projection.z;
//	projection.y = projection.y / projection.z;
//	projection.z = projection.z / projection.z;
//
//	//Obtain the bounding box in the projection plane;
//	minX = MY_MIN<float>(minX, projection.x);
//	maxX = MY_MAX<float>(maxX, projection.x);
//	minZ = MY_MIN<float>(minZ, projection.y);
//	maxZ = MY_MAX<float>(maxZ, projection.y);
//
//
//
//
//	origPos.x -= Vol.m_Step.x;
//	origPos.y += Vol.m_Step.y;
//	//Trans from the world coordinate to detector index coordinate;
//	projection.x = M[0] * origPos.x + M[1] * origPos.y + M[2];
//	projection.y = M[3] * origPos.x + M[4] * origPos.y + M[5] * origPos.z + M[6];
//	projection.z = M[7] * origPos.x + M[8] * origPos.y + M[9];
//	projection.x = projection.x / projection.z;
//	projection.y = projection.y / projection.z;
//	projection.z = projection.z / projection.z;
//
//	//Obtain the bounding box in the projection plane;
//	minX = MY_MIN<float>(minX, projection.x);
//	maxX = MY_MAX<float>(maxX, projection.x);
//	minZ = MY_MIN<float>(minZ, projection.y);
//	maxZ = MY_MAX<float>(maxZ, projection.y);
//
//
//
//	origPos.x += Vol.m_Step.x;
//	//Trans from the world coordinate to detector index coordinate;
//	projection.x = M[0] * origPos.x + M[1] * origPos.y + M[2];
//	projection.y = M[3] * origPos.x + M[4] * origPos.y + M[5] * origPos.z + M[6];
//	projection.z = M[7] * origPos.x + M[8] * origPos.y + M[9];
//	projection.x = projection.x / projection.z;
//	projection.y = projection.y / projection.z;
//	projection.z = projection.z / projection.z;
//
//	//Obtain the bounding box in the projection plane;
//	minX = MY_MIN<float>(minX, projection.x);
//	maxX = MY_MAX<float>(maxX, projection.x);
//	minZ = MY_MIN<float>(minZ, projection.y);
//	maxZ = MY_MAX<float>(maxZ, projection.y);
//
//
//
//	origPos.x -= Vol.m_Step.x;
//	origPos.y -= Vol.m_Step.y;
//	origPos.z += Vol.m_Step.z;
//	//Trans from the world coordinate to detector index coordinate;
//	projection.x = M[0] * origPos.x + M[1] * origPos.y + M[2];
//	projection.y = M[3] * origPos.x + M[4] * origPos.y + M[5] * origPos.z + M[6];
//	projection.z = M[7] * origPos.x + M[8] * origPos.y + M[9];
//	projection.x = projection.x / projection.z;
//	projection.y = projection.y / projection.z;
//	projection.z = projection.z / projection.z;
//
//	//Obtain the bounding box in the projection plane;
//	minX = MY_MIN<float>(minX, projection.x);
//	maxX = MY_MAX<float>(maxX, projection.x);
//	minZ = MY_MIN<float>(minZ, projection.y);
//	maxZ = MY_MAX<float>(maxZ, projection.y);
//
//
//
//	origPos.x += Vol.m_Step.x;
//	//Trans from the world coordinate to detector index coordinate;
//	projection.x = M[0] * origPos.x + M[1] * origPos.y + M[2];
//	projection.y = M[3] * origPos.x + M[4] * origPos.y + M[5] * origPos.z + M[6];
//	projection.z = M[7] * origPos.x + M[8] * origPos.y + M[9];
//	projection.x = projection.x / projection.z;
//	projection.y = projection.y / projection.z;
//	projection.z = projection.z / projection.z;
//
//	//Obtain the bounding box in the projection plane;
//	minX = MY_MIN<float>(minX, projection.x);
//	maxX = MY_MAX<float>(maxX, projection.x);
//	minZ = MY_MIN<float>(minZ, projection.y);
//	maxZ = MY_MAX<float>(maxZ, projection.y);
//
//
//
//
//	origPos.x -= Vol.m_Step.x;
//	origPos.y += Vol.m_Step.y;
//	//Trans from the world coordinate to detector index coordinate;
//	projection.x = M[0] * origPos.x + M[1] * origPos.y + M[2];
//	projection.y = M[3] * origPos.x + M[4] * origPos.y + M[5] * origPos.z + M[6];
//	projection.z = M[7] * origPos.x + M[8] * origPos.y + M[9];
//	projection.x = projection.x / projection.z;
//	projection.y = projection.y / projection.z;
//	projection.z = projection.z / projection.z;
//
//	//Obtain the bounding box in the projection plane;
//	minX = MY_MIN<float>(minX, projection.x);
//	maxX = MY_MAX<float>(maxX, projection.x);
//	minZ = MY_MIN<float>(minZ, projection.y);
//	maxZ = MY_MAX<float>(maxZ, projection.y);
//
//
//
//	origPos.x += Vol.m_Step.x;
//	//Trans from the world coordinate to detector index coordinate;
//	projection.x = M[0] * origPos.x + M[1] * origPos.y + M[2];
//	projection.y = M[3] * origPos.x + M[4] * origPos.y + M[5] * origPos.z + M[6];
//	projection.z = M[7] * origPos.x + M[8] * origPos.y + M[9];
//	projection.x = projection.x / projection.z;
//	projection.y = projection.y / projection.z;
//	projection.z = projection.z / projection.z;
//
//	//Obtain the bounding box in the projection plane;
//	minX = MY_MIN<float>(minX, projection.x);
//	maxX = MY_MAX<float>(maxX, projection.x);
//	minZ = MY_MIN<float>(minZ, projection.y);
//	maxZ = MY_MAX<float>(maxZ, projection.y);
//
//
//
//	if (minX != 10000)
//		minX = (minX < 0) ? 0 : minX;
//	if (maxX != -10000)
//		maxX = (maxX > ConeGeo.m_DetN.x - 1) ? ConeGeo.m_DetN.x - 1 : maxX;
//	if (minZ != 10000)
//		minZ = (minZ < 0) ? 0 : minZ;
//	if (maxZ != -10000)
//		maxZ = (maxZ > ConeGeo.m_DetN.y - 1) ? ConeGeo.m_DetN.y - 1 : maxZ;
//
//
//
//	int x(0), y(0);
//	int hit(0);
//	float sumWeight(0), sumError(0);
//	for (x = minX; x <= maxX; ++x)
//	{
//		for (y = minZ; y <= maxZ; ++y)
//		{
//			p = make_float3(
//				MIND.x + (x + 0.5f) * ConeGeo.m_DetStp.x,
//				-ConeGeo.m_O2D,
//				MIND.y + (y + 0.5f) * ConeGeo.m_DetStp.y);
//			curP = rotation(p, cosT, sinT);
//
//			eyeRay.d = curP - eyeRay.o;
//			eyeRay.d = normalize(eyeRay.d);
//
//			hit = intersectBox(eyeRay, boxMin, boxMax, &tnear, &tfar);
//			if (hit)
//			{
//				sumWeight += (tfar - tnear);
//				sumError += (tfar - tnear) * correction[y * ConeGeo.m_DetN.x + x];
//			}
//		}
//	}
//	if (!IS_ZERO(sumWeight))
//		volume[index] += (sumError / sumWeight * dmask[idx.y * Vol.m_Reso.x + idx.x]);
//}
//void bakproj_BOXED(float* donp, float* dvol, float* dmsk, cuint angIdx, const ConeEDGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid)
//{
//	_bakproj_BOXED_ker << <gid, blk >> >(donp, dvol, dmsk, angIdx, ConeGeo, Vol);
//}
//
//
//__global__ void _bakproj_BOXED_ker(float* correction, float* volume, const unsigned int angIdx, const ConeEDGeo ConeGeo, const Volume Vol)
//{
//	/*uint3 idx;
//	idx.x = threadIdx.x + blockIdx.x * blockDim.x;
//	idx.y = threadIdx.y + blockIdx.y * blockDim.y;
//	idx.z = threadIdx.z + blockIdx.z * blockDim.z;*/
//	//cuint index = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint index = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
//	if (index >= Vol.m_Reso.x * Vol.m_Reso.y * Vol.m_Reso.z) return; //if (idx.x < Vol.m_Reso.x && idx.y < Vol.m_Reso.y && idx.z < Vol.m_Reso.z)
//
//	const float cosT = cosf(ConeGeo.m_ViwBeg + (float)angIdx * ConeGeo.m_ViwStp);
//	const float sinT = sinf(ConeGeo.m_ViwBeg + (float)angIdx * ConeGeo.m_ViwStp);
//	float2 MIND = make_float2(
//		-ConeGeo.m_DetCntIdx.x * ConeGeo.m_DetStp.x,
//		-ConeGeo.m_DetCntIdx.y * ConeGeo.m_DetStp.y);
//
//	float M[10];
//	M[0] = (-(ConeGeo.m_S2D) * cosT / ConeGeo.m_DetStp.x + MIND.x * sinT / ConeGeo.m_DetStp.x);
//	M[1] = (-(ConeGeo.m_S2D) * sinT / ConeGeo.m_DetStp.x - MIND.x * cosT / ConeGeo.m_DetStp.x);
//	M[2] = MIND.x * ConeGeo.m_S2O / ConeGeo.m_DetStp.x;
//	M[3] = MIND.y / ConeGeo.m_DetStp.y * sinT;
//	M[4] = (-MIND.y / ConeGeo.m_DetStp.y) * cosT;
//	M[5] = (-(ConeGeo.m_S2D) / ConeGeo.m_DetStp.y);
//	M[6] = MIND.y * ConeGeo.m_S2O / ConeGeo.m_DetStp.y;
//	M[7] = -sinT;
//	M[8] = cosT;
//	M[9] = -ConeGeo.m_S2O;
//
//	/*__shared__ float M[10];
//	M[1] = (-(__S2O__ + __O2D__) * sinT / __DETSTPL__ - __MINDL__ * cosT / __DETSTPL__);
//	M[2] = __MINDL__ * __S2O__ / __DETSTPL__;
//	M[3] = __MINDH__ / __DETSTPH__ * sinT;
//	M[4] = (-__MINDH__ / __DETSTPH__) * cosT;
//	M[5] = (-(__S2O__ + __O2D__) / __DETSTPH__);
//	M[6] = __MINDH__ * __S2O__ / __DETSTPH__;
//	M[7] = -sinT;
//	M[8] = cosT;
//	M[9] = -__S2O__;*/
//	//__syncthreads();
//
//
//	uint3 idx;
//	idx.z = index / (static_cast<int>(Vol.m_Reso.x) * static_cast<int>(Vol.m_Reso.y));
//	idx.y = (index - idx.z * (static_cast<int>(Vol.m_Reso.x) * static_cast<int>(Vol.m_Reso.y))) / static_cast<int>(Vol.m_Reso.x);
//	idx.x = index - idx.z * (static_cast<int>(Vol.m_Reso.x) * static_cast<int>(Vol.m_Reso.y)) - idx.y * static_cast<int>(Vol.m_Reso.x);
//
//	float minX = 10000, minZ = 10000, maxX = -10000, maxZ = -10000;
//
//	float3 projection, origPos, p, curP;
//	//int i(0), j(0), k(0);
//	//float tt;
//	Ray eyeRay;
//	eyeRay.o = make_float3(-ConeGeo.m_S2O * sinT, ConeGeo.m_S2O * cosT, 0);
//
//	float tnear = 0.0f;
//	float tfar = 0.0f;
//	float3 boxMin, boxMax;
//	float3 MINO = make_float3(
//		-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//		-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//		-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//
//	boxMin.x = MINO.x + idx.x * Vol.m_Step.x;
//	boxMin.y = MINO.y + idx.y * Vol.m_Step.y;
//	boxMin.z = MINO.z + idx.z * Vol.m_Step.z;
//	boxMax.x = boxMin.x + Vol.m_Step.x;
//	boxMax.y = boxMin.y + Vol.m_Step.y;
//	boxMax.z = boxMin.z + Vol.m_Step.z;
//
//
//	origPos.x = MINO.x + idx.x * Vol.m_Step.x;
//	origPos.y = MINO.y + idx.y * Vol.m_Step.y;
//	origPos.z = MINO.z + idx.z * Vol.m_Step.z;
//	//Trans from the world coordinate to detector index coordinate;
//	projection.x = M[0] * origPos.x + M[1] * origPos.y + M[2];
//	projection.y = M[3] * origPos.x + M[4] * origPos.y + M[5] * origPos.z + M[6];
//	projection.z = M[7] * origPos.x + M[8] * origPos.y + M[9];
//	projection.x = projection.x / projection.z;
//	projection.y = projection.y / projection.z;
//	projection.z = projection.z / projection.z;
//
//	//Obtain the bounding box in the projection plane;
//	minX = MY_MIN<float>(minX, projection.x);
//	maxX = MY_MAX<float>(maxX, projection.x);
//	minZ = MY_MIN<float>(minZ, projection.y);
//	maxZ = MY_MAX<float>(maxZ, projection.y);
//
//
//	origPos.x += Vol.m_Step.x;
//	//Trans from the world coordinate to detector index coordinate;
//	projection.x = M[0] * origPos.x + M[1] * origPos.y + M[2];
//	projection.y = M[3] * origPos.x + M[4] * origPos.y + M[5] * origPos.z + M[6];
//	projection.z = M[7] * origPos.x + M[8] * origPos.y + M[9];
//	projection.x = projection.x / projection.z;
//	projection.y = projection.y / projection.z;
//	projection.z = projection.z / projection.z;
//
//	//Obtain the bounding box in the projection plane;
//	minX = MY_MIN<float>(minX, projection.x);
//	maxX = MY_MAX<float>(maxX, projection.x);
//	minZ = MY_MIN<float>(minZ, projection.y);
//	maxZ = MY_MAX<float>(maxZ, projection.y);
//
//
//
//
//	origPos.x -= Vol.m_Step.x;
//	origPos.y += Vol.m_Step.y;
//	//Trans from the world coordinate to detector index coordinate;
//	projection.x = M[0] * origPos.x + M[1] * origPos.y + M[2];
//	projection.y = M[3] * origPos.x + M[4] * origPos.y + M[5] * origPos.z + M[6];
//	projection.z = M[7] * origPos.x + M[8] * origPos.y + M[9];
//	projection.x = projection.x / projection.z;
//	projection.y = projection.y / projection.z;
//	projection.z = projection.z / projection.z;
//
//	//Obtain the bounding box in the projection plane;
//	minX = MY_MIN<float>(minX, projection.x);
//	maxX = MY_MAX<float>(maxX, projection.x);
//	minZ = MY_MIN<float>(minZ, projection.y);
//	maxZ = MY_MAX<float>(maxZ, projection.y);
//
//
//
//	origPos.x += Vol.m_Step.x;
//	//Trans from the world coordinate to detector index coordinate;
//	projection.x = M[0] * origPos.x + M[1] * origPos.y + M[2];
//	projection.y = M[3] * origPos.x + M[4] * origPos.y + M[5] * origPos.z + M[6];
//	projection.z = M[7] * origPos.x + M[8] * origPos.y + M[9];
//	projection.x = projection.x / projection.z;
//	projection.y = projection.y / projection.z;
//	projection.z = projection.z / projection.z;
//
//	//Obtain the bounding box in the projection plane;
//	minX = MY_MIN<float>(minX, projection.x);
//	maxX = MY_MAX<float>(maxX, projection.x);
//	minZ = MY_MIN<float>(minZ, projection.y);
//	maxZ = MY_MAX<float>(maxZ, projection.y);
//
//
//
//	origPos.x -= Vol.m_Step.x;
//	origPos.y -= Vol.m_Step.y;
//	origPos.z += Vol.m_Step.z;
//	//Trans from the world coordinate to detector index coordinate;
//	projection.x = M[0] * origPos.x + M[1] * origPos.y + M[2];
//	projection.y = M[3] * origPos.x + M[4] * origPos.y + M[5] * origPos.z + M[6];
//	projection.z = M[7] * origPos.x + M[8] * origPos.y + M[9];
//	projection.x = projection.x / projection.z;
//	projection.y = projection.y / projection.z;
//	projection.z = projection.z / projection.z;
//
//	//Obtain the bounding box in the projection plane;
//	minX = MY_MIN<float>(minX, projection.x);
//	maxX = MY_MAX<float>(maxX, projection.x);
//	minZ = MY_MIN<float>(minZ, projection.y);
//	maxZ = MY_MAX<float>(maxZ, projection.y);
//
//
//
//	origPos.x += Vol.m_Step.x;
//	//Trans from the world coordinate to detector index coordinate;
//	projection.x = M[0] * origPos.x + M[1] * origPos.y + M[2];
//	projection.y = M[3] * origPos.x + M[4] * origPos.y + M[5] * origPos.z + M[6];
//	projection.z = M[7] * origPos.x + M[8] * origPos.y + M[9];
//	projection.x = projection.x / projection.z;
//	projection.y = projection.y / projection.z;
//	projection.z = projection.z / projection.z;
//
//	//Obtain the bounding box in the projection plane;
//	minX = MY_MIN<float>(minX, projection.x);
//	maxX = MY_MAX<float>(maxX, projection.x);
//	minZ = MY_MIN<float>(minZ, projection.y);
//	maxZ = MY_MAX<float>(maxZ, projection.y);
//
//
//
//
//	origPos.x -= Vol.m_Step.x;
//	origPos.y += Vol.m_Step.y;
//	//Trans from the world coordinate to detector index coordinate;
//	projection.x = M[0] * origPos.x + M[1] * origPos.y + M[2];
//	projection.y = M[3] * origPos.x + M[4] * origPos.y + M[5] * origPos.z + M[6];
//	projection.z = M[7] * origPos.x + M[8] * origPos.y + M[9];
//	projection.x = projection.x / projection.z;
//	projection.y = projection.y / projection.z;
//	projection.z = projection.z / projection.z;
//
//	//Obtain the bounding box in the projection plane;
//	minX = MY_MIN<float>(minX, projection.x);
//	maxX = MY_MAX<float>(maxX, projection.x);
//	minZ = MY_MIN<float>(minZ, projection.y);
//	maxZ = MY_MAX<float>(maxZ, projection.y);
//
//
//
//	origPos.x += Vol.m_Step.x;
//	//Trans from the world coordinate to detector index coordinate;
//	projection.x = M[0] * origPos.x + M[1] * origPos.y + M[2];
//	projection.y = M[3] * origPos.x + M[4] * origPos.y + M[5] * origPos.z + M[6];
//	projection.z = M[7] * origPos.x + M[8] * origPos.y + M[9];
//	projection.x = projection.x / projection.z;
//	projection.y = projection.y / projection.z;
//	projection.z = projection.z / projection.z;
//
//	//Obtain the bounding box in the projection plane;
//	minX = MY_MIN<float>(minX, projection.x);
//	maxX = MY_MAX<float>(maxX, projection.x);
//	minZ = MY_MIN<float>(minZ, projection.y);
//	maxZ = MY_MAX<float>(maxZ, projection.y);
//
//
//
//	if (minX != 10000)
//		minX = (minX < 0) ? 0 : minX;
//	if (maxX != -10000)
//		maxX = (maxX > ConeGeo.m_DetN.x - 1) ? ConeGeo.m_DetN.x - 1 : maxX;
//	if (minZ != 10000)
//		minZ = (minZ < 0) ? 0 : minZ;
//	if (maxZ != -10000)
//		maxZ = (maxZ > ConeGeo.m_DetN.y - 1) ? ConeGeo.m_DetN.y - 1 : maxZ;
//
//
//
//	int x(0), y(0);
//	int hit(0);
//	float sumWeight(0), sumError(0);
//	for (x = minX; x <= maxX; ++x)
//	{
//		for (y = minZ; y <= maxZ; ++y)
//		{
//			p = make_float3(
//				MIND.x + (x + 0.5f) * ConeGeo.m_DetStp.x,
//				-ConeGeo.m_O2D,
//				MIND.y + (y + 0.5f) * ConeGeo.m_DetStp.y);
//			curP = rotation(p, cosT, sinT);
//
//			eyeRay.d = curP - eyeRay.o;
//			eyeRay.d = normalize(eyeRay.d);
//
//			hit = intersectBox(eyeRay, boxMin, boxMax, &tnear, &tfar);
//			if (hit)
//			{
//				sumWeight += (tfar - tnear);
//				sumError += (tfar - tnear) * correction[(angIdx * ConeGeo.m_DetN.y + y) * ConeGeo.m_DetN.x + x];
//			}
//		}
//	}
//	if (!IS_ZERO(sumWeight))
//		volume[index] += (sumError / sumWeight);
//}
//void bakproj_BOXED(float* donp, float* dvol, cuint angIdx, const ConeEDGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid)
//{
//	_bakproj_BOXED_ker << <gid, blk >> >(donp, dvol, angIdx, ConeGeo, Vol);
//}
//
//
//
//template<typename T>
//inline __host__ __device__ void GenM(T* M, const ConeEDGeo& ConeGeo, const float2& MIND, const T& cosT, const T& sinT)
//{
//	M[0] = (-(ConeGeo.m_S2D) * cosT / ConeGeo.m_DetStp.x + MIND.x * sinT / ConeGeo.m_DetStp.x);
//	M[1] = (-(ConeGeo.m_S2D) * sinT / ConeGeo.m_DetStp.x - MIND.x * cosT / ConeGeo.m_DetStp.x);
//	M[2] = MIND.x * ConeGeo.m_S2O / ConeGeo.m_DetStp.x;
//	M[3] = MIND.y / ConeGeo.m_DetStp.y * sinT;
//	M[4] = (-MIND.y / ConeGeo.m_DetStp.y) * cosT;
//	M[5] = (-(ConeGeo.m_S2D) / ConeGeo.m_DetStp.y);
//	M[6] = MIND.y * ConeGeo.m_S2O / ConeGeo.m_DetStp.y;
//	M[7] = -sinT;
//	M[8] = cosT;
//	M[9] = -ConeGeo.m_S2O;
//}
//
//
//__global__ void _bakproj_BOXED_ker(float* correction, float* volume, const ConeEDGeo ConeGeo, const Volume Vol)
//{
//	cuint index = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
//	if (index >= Vol.m_Reso.x * Vol.m_Reso.y * Vol.m_Reso.z) return;
//	unsigned int angIdx = 0;
//	float2 MIND = make_float2(
//		-ConeGeo.m_DetCntIdx.x * ConeGeo.m_DetStp.x,
//		-ConeGeo.m_DetCntIdx.y * ConeGeo.m_DetStp.y);
//
//	uint3 idx;
//	idx.z = index / (static_cast<int>(Vol.m_Reso.x) * static_cast<int>(Vol.m_Reso.y));
//	idx.y = (index - idx.z * (static_cast<int>(Vol.m_Reso.x) * static_cast<int>(Vol.m_Reso.y))) / static_cast<int>(Vol.m_Reso.x);
//	idx.x = index - idx.z * (static_cast<int>(Vol.m_Reso.x) * static_cast<int>(Vol.m_Reso.y)) - idx.y * static_cast<int>(Vol.m_Reso.x);
//
//	Ray eyeRay;
//
//
//	float3 projection, origPos, curP;
//	float3 MINO = make_float3(
//		-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//		-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//		-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//
//	float3 boxMin(make_float3(MINO.x + idx.x * Vol.m_Step.x, MINO.y + idx.y * Vol.m_Step.y, MINO.z + idx.z * Vol.m_Step.z));
//	float3 boxMax(make_float3(boxMin.x + Vol.m_Step.x, boxMin.y + Vol.m_Step.y, boxMin.z + Vol.m_Step.z));
//
//	int kk(0), jj(0), ii(0), x(0), y(0), hit(0);
//
//	float cosT, sinT;
//	float M[10];
//	float minX = 10000, minZ = 10000, maxX = -10000, maxZ = -10000;
//	float sumWeight(0.0f), sumError(0.0f), summ(0.0f), tnear(0.0f), tfar(0.0f);
//
//
//	for (angIdx = 0; angIdx != ConeGeo.m_ViwN; ++angIdx)
//	{
//		cosT = cosf(ConeGeo.m_ViwBeg + (float)angIdx * ConeGeo.m_ViwStp);
//		sinT = sinf(ConeGeo.m_ViwBeg + (float)angIdx * ConeGeo.m_ViwStp);
//
//		GenM<float>(M, ConeGeo, MIND, cosT, sinT);
//
//
//		minX = 10000;
//		minZ = 10000;
//		maxX = -10000;
//		maxZ = -10000;
//
//		eyeRay.o = make_float3(-ConeGeo.m_S2O * sinT, ConeGeo.m_S2O * cosT, 0);
//
//		for (kk = 0; kk != 2; ++kk)
//		{
//			for (jj = 0; jj != 2; ++jj)
//			{
//				for (ii = 0; ii != 2; ++ii)
//				{
//					origPos.x = MINO.x + (idx.x + ii) * Vol.m_Step.x;
//					origPos.y = MINO.y + (idx.y + jj) * Vol.m_Step.y;
//					origPos.z = MINO.z + (idx.z + kk) * Vol.m_Step.z;
//					//Trans from the world coordinate to detector index coordinate;
//					projection.x = M[0] * origPos.x + M[1] * origPos.y + M[2];
//					projection.y = M[3] * origPos.x + M[4] * origPos.y + M[5] * origPos.z + M[6];
//					projection.z = M[7] * origPos.x + M[8] * origPos.y + M[9];
//					projection.x = projection.x / projection.z;
//					projection.y = projection.y / projection.z;
//					projection.z = projection.z / projection.z;
//
//					//Obtain the bounding box in the projection plane;
//					minX = MY_MIN<float>(minX, projection.x);
//					maxX = MY_MAX<float>(maxX, projection.x);
//					minZ = MY_MIN<float>(minZ, projection.y);
//					maxZ = MY_MAX<float>(maxZ, projection.y);
//				}
//			}
//		}
//
//		if (minX != 10000)
//			minX = (minX < 0) ? 0 : minX;
//		if (maxX != -10000)
//			maxX = (maxX > ConeGeo.m_DetN.x - 1) ? ConeGeo.m_DetN.x - 1 : maxX;
//		if (minZ != 10000)
//			minZ = (minZ < 0) ? 0 : minZ;
//		if (maxZ != -10000)
//			maxZ = (maxZ > ConeGeo.m_DetN.y - 1) ? ConeGeo.m_DetN.y - 1 : maxZ;
//
//		hit = 0;
//		sumWeight = 0;
//		sumError = 0;
//		for (x = minX; x <= maxX; ++x)
//		{
//			for (y = minZ; y <= maxZ; ++y)
//			{
//				curP = rotation(make_float3(
//					MIND.x + (x + 0.5f) * ConeGeo.m_DetStp.x,
//					-ConeGeo.m_O2D,
//					MIND.y + (y + 0.5f) * ConeGeo.m_DetStp.y), cosT, sinT);
//
//				eyeRay.d = normalize(curP - eyeRay.o);
//				tnear = 0.0f; tfar = 0.0f;
//				hit = intersectBox(eyeRay, boxMin, boxMax, &tnear, &tfar);
//				if (hit)
//				{
//					sumWeight += (tfar - tnear);
//					sumError += (tfar - tnear) * correction[(angIdx * ConeGeo.m_DetN.y + y) * ConeGeo.m_DetN.x + x];
//				}
//			}
//		}
//		if (!IS_ZERO(sumWeight))
//			summ += (sumError / sumWeight);
//	}//End for each angle
//	volume[index] = summ;
//}
//void bakproj_BOXED(float* donp, float* dvol, const ConeEDGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid)
//{
//	_bakproj_BOXED_ker << <gid, blk >> >(donp, dvol, ConeGeo, Vol);
//}
//
//
//__global__ void _bakproj_BOXED_ker(float* dprj, float* dimg, const FanEAGeo FanGeo, const Image Img)
//{
//	cuint i = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint j = threadIdx.y + blockIdx.y * blockDim.y;
//	if (i < Img.m_Reso.x && j < Img.m_Reso.y)
//	{
//		float2 minBox = make_float2(
//			-Img.m_Size.x * 0.5f + Img.m_Bias.x + i * Img.m_Step.x,
//			-Img.m_Size.y * 0.5f + Img.m_Bias.y + j * Img.m_Step.y);
//		float2 maxBox = make_float2(
//			-Img.m_Size.x * 0.5f + Img.m_Bias.x + (i + 1) * Img.m_Step.x,
//			-Img.m_Size.y * 0.5f + Img.m_Bias.y + (j + 1) * Img.m_Step.y);
//
//		float2 initD; //original detector position
//
//		unsigned int angIdx = 0;
//		Ray2D ray;
//		int minID = 10000;
//		int maxID = -10000;
//		int curID(0);
//		float summ(0.0f), tnear(0.0f), tfar(0.0f);
//		float cosT(0.0f);// = cosf(FanGeo.m_ViwBeg + (float) * )
//		float sinT(0.0f);// sinf
//		float detID(0.0f);
//		bool hit = false;
//		float curAng(0.0f);
//		//float2 curDet;
//		for (angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
//		{
//			minID = 10000;
//			maxID = -10000;
//			curAng = FanGeo.m_ViwBeg + (float)angIdx * FanGeo.m_ViwStp;
//
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT); //current ray position;
//
//			ray.d = normalize(minBox - ray.o);
//			initD = rotation(ray.o + FanGeo.m_S2D * ray.d, cosT, -sinT);
//			detID = atanf(initD.x / (FanGeo.m_S2O - initD.y)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			minID = (detID < minID) ? detID : minID;
//			maxID = (detID > maxID) ? ceilf(detID) : maxID;
//			ray.d = normalize(minBox + make_float2(Img.m_Step.x, 0) - ray.o);
//			initD = rotation(ray.o + FanGeo.m_S2D * ray.d, cosT, -sinT);
//			detID = atanf(initD.x / (FanGeo.m_S2O - initD.y)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			minID = (detID < minID) ? detID : minID;
//			maxID = (detID > maxID) ? ceilf(detID) : maxID;
//			ray.d = normalize(minBox + make_float2(0, Img.m_Step.y) - ray.o);
//			initD = rotation(ray.o + FanGeo.m_S2D * ray.d, cosT, -sinT);
//			detID = atanf(initD.x / (FanGeo.m_S2O - initD.y)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			minID = (detID < minID) ? detID : minID;
//			maxID = (detID > maxID) ? ceilf(detID) : maxID;
//			ray.d = normalize(maxBox - ray.o);
//			initD = rotation(ray.o + FanGeo.m_S2D * ray.d, cosT, -sinT);
//			detID = atanf(initD.x / (FanGeo.m_S2O - initD.y)) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			minID = (detID < minID) ? detID : minID;
//			maxID = (detID > maxID) ? ceilf(detID) : maxID;
//
//			if (minID != 10000){ minID = (minID < 0) ? 0 : minID; }
//			if (maxID != -10000){ maxID = (maxID > FanGeo.m_DetN - 1) ? FanGeo.m_DetN - 1 : maxID; }
//
//			for (curID = minID; curID <= maxID; ++curID)
//			{
//				curAng = (curID - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp;
//				initD = rotation(make_float2(sinf(curAng) * FanGeo.m_S2D,
//					-cosf(curAng)*FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//				ray.d = normalize(initD - ray.o);
//				tnear = 0.0f; tfar = 0.0f;
//				hit = intersectBox(ray, minBox, maxBox, &tnear, &tfar);
//				if (hit)
//				{
//					summ += (tfar - tnear) * dprj[angIdx * FanGeo.m_DetN + curID];
//				}
//			}// end curID
//		}//end angIdx
//		dimg[j * Img.m_Reso.x + i] = summ;
//	}
//}
//void bakproj_BOXED(float* dprj, float* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_bakproj_BOXED_ker << <gid, blk >> >(dprj, dimg, FanGeo, Img);
//}
//
//
//
//
//template<typename T>
//__global__ void _bakproj_AIM_ker(T* dimg, T* dprj, const FanEAGeo FanGeo, const Image Img)
//{
//	cuint xi = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint yi = threadIdx.y + blockIdx.y * blockDim.y;
//
//	if (xi < Img.m_Reso.x&&yi < Img.m_Reso.y)
//	{
//		const T cntImgX = static_cast<T>(Img.m_Reso.x - 1.0) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
//		const T cntImgY = static_cast<T>(Img.m_Reso.y - 1.0) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
//		const T area = Img.m_Step.x * Img.m_Step.y;
//		T sour[2], cosT, sinT;
//		T grid[4][3];
//		grid[0][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
//		grid[0][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
//
//		grid[1][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
//		grid[1][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
//
//		grid[2][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
//		grid[2][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
//
//		grid[3][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
//		grid[3][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
//		T det[4][3];
//		T tmp[2];
//		T beta[4];
//		int minDetIdx, maxDetIdx;
//		T summ = 0;
//		T SVA[3], SVB[3];
//		for (size_t angIdx = 0; angIdx < FanGeo.m_ViwN; ++angIdx)
//		{
//			cosT = cos(FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp);
//			sinT = sin(FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp);
//
//			sour[0] = -FanGeo.m_S2O * sinT;
//			sour[1] = FanGeo.m_S2O * cosT;
//
//			det[0][0] = grid[0][0] - sour[0]; //used as the direction from source to the corner of the pixel
//			det[0][1] = grid[0][1] - sour[1];
//			det[1][0] = grid[1][0] - sour[0];
//			det[1][1] = grid[1][1] - sour[1];
//			det[2][0] = grid[2][0] - sour[0];
//			det[2][1] = grid[2][1] - sour[1];
//			det[3][0] = grid[3][0] - sour[0];
//			det[3][1] = grid[3][1] - sour[1];
//
//
//			tmp[1] = sqrt(det[0][0] * det[0][0] + det[0][1] * det[0][1]);
//			det[0][0] /= tmp[1];
//			det[0][1] /= tmp[1];
//			tmp[1] = sqrt(det[1][0] * det[1][0] + det[1][1] * det[1][1]);
//			det[1][0] /= tmp[1];
//			det[1][1] /= tmp[1];
//			tmp[1] = sqrt(det[2][0] * det[2][0] + det[2][1] * det[2][1]);
//			det[2][0] /= tmp[1];
//			det[2][1] /= tmp[1];
//			tmp[1] = sqrt(det[3][0] * det[3][0] + det[3][1] * det[3][1]);
//			det[3][0] /= tmp[1];
//			det[3][1] /= tmp[1];
//
//			det[0][0] = sour[0] + det[0][0] * FanGeo.m_S2D;  //det used as the detector coordinates
//			det[0][1] = sour[1] + det[0][1] * FanGeo.m_S2D;
//			det[1][0] = sour[0] + det[1][0] * FanGeo.m_S2D;
//			det[1][1] = sour[1] + det[1][1] * FanGeo.m_S2D;
//			det[2][0] = sour[0] + det[2][0] * FanGeo.m_S2D;
//			det[2][1] = sour[1] + det[2][1] * FanGeo.m_S2D;
//			det[3][0] = sour[0] + det[3][0] * FanGeo.m_S2D;
//			det[3][1] = sour[1] + det[3][1] * FanGeo.m_S2D;
//
//			tmp[0] = det[0][0] * cosT + det[0][1] * sinT;
//			tmp[1] = -det[0][0] * sinT + det[0][1] * cosT;
//			det[0][0] = tmp[0];
//			det[0][1] = tmp[1];
//			tmp[0] = det[1][0] * cosT + det[1][1] * sinT;
//			tmp[1] = -det[1][0] * sinT + det[1][1] * cosT;
//			det[1][0] = tmp[0];
//			det[1][1] = tmp[1];
//			tmp[0] = det[2][0] * cosT + det[2][1] * sinT;
//			tmp[1] = -det[2][0] * sinT + det[2][1] * cosT;
//			det[2][0] = tmp[0];
//			det[2][1] = tmp[1];
//			tmp[0] = det[3][0] * cosT + det[3][1] * sinT;
//			tmp[1] = -det[3][0] * sinT + det[3][1] * cosT;
//			det[3][0] = tmp[0];
//			det[3][1] = tmp[1];
//
//
//			beta[0] = atan(-det[0][0] / (det[0][1] - FanGeo.m_S2O));
//			beta[1] = atan(-det[1][0] / (det[1][1] - FanGeo.m_S2O));
//			beta[2] = atan(-det[2][0] / (det[2][1] - FanGeo.m_S2O));
//			beta[3] = atan(-det[3][0] / (det[3][1] - FanGeo.m_S2O));
//
//			grid[0][2] = beta[0] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			grid[1][2] = beta[1] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			grid[2][2] = beta[2] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			grid[3][2] = beta[3] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//			SortProjection<T>(grid);
//			minDetIdx = static_cast<int>(grid[0][2]);
//			maxDetIdx = (int(ceil(grid[3][2])));
//			minDetIdx = (minDetIdx < 0) ? 0 : minDetIdx;
//			maxDetIdx = (maxDetIdx > static_cast<int>(FanGeo.m_DetN)) ? static_cast<int>(FanGeo.m_DetN) : maxDetIdx;
//
//			tmp[1] = (hypot((xi - cntImgX)*Img.m_Step.x - sour[0], (yi - cntImgY)*Img.m_Step.y - sour[1]));
//			//T summ = 0;
//			//T pangle;
//			//T direct[2], legth, SVA[3], SVB[3], coef;
//
//
//			//�����߱߼���;
//			beta[3] = (minDetIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp; //���С�߱�ƫ�ƽ� // Beta[3] used as the pangle
//			beta[0] = sin(beta[3]) * FanGeo.m_S2D; // beta used as the initial X,Y and current X,Y with rotation
//			beta[1] = -cos(beta[3]) * FanGeo.m_S2D + FanGeo.m_S2O;
//
//			beta[2] = beta[0] * cosT - beta[1] * sinT;
//			beta[3] = beta[0] * sinT + beta[1] * cosT;
//			beta[0] = beta[2] - sour[0]; //Beta used as the directions
//			beta[1] = beta[3] - sour[1];
//			beta[2] = sqrt(beta[0] * beta[0] + beta[1] * beta[1]);//Beta[3] used as the length of the image
//			SVA[0] = beta[0] / beta[2];
//			SVA[1] = beta[1] / beta[2];
//			SVA[2] = static_cast<T>(minDetIdx);
//			tmp[0] = 0;
//			for (; minDetIdx < maxDetIdx; ++minDetIdx)
//			{
//				//�����߱߼���;
//
//				beta[3] = (minDetIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp + FanGeo.m_DetStp; //��߱�ƫ�ƽ�
//				beta[0] = sin(beta[3]) * FanGeo.m_S2D;
//				beta[1] = -cos(beta[3]) * FanGeo.m_S2D + FanGeo.m_S2O;
//				beta[2] = beta[0] * cosT - beta[1] * sinT;
//				beta[3] = beta[0] * sinT + beta[1] * cosT;
//				beta[0] = beta[2] - sour[0];
//				beta[1] = beta[3] - sour[1];
//				beta[2] = sqrt(beta[0] * beta[0] + beta[1] * beta[1]);
//				SVB[0] = beta[0] / beta[2];
//				SVB[1] = beta[1] / beta[2];
//				SVB[2] = static_cast<T>(minDetIdx)+1;
//
//				//Compute the weighting coefficient for a special projection data
//				beta[3] = ComputeCoefficient<T>(grid, SVA, SVB, sour, area); // beta[3] used as the coef
//				beta[3] /= (tmp[1] * fabs(FanGeo.m_DetStp));
//				tmp[0] += dprj[angIdx* FanGeo.m_DetN + minDetIdx] * beta[3];
//
//				SVA[0] = SVB[0]; //��߱߱��С�߱�
//				SVA[1] = SVB[1];
//				SVA[2] = SVB[2];
//			}
//			summ += tmp[0];
//		}
//		dimg[yi* Img.m_Reso.x + xi] = summ;
//	}
//}
//
//
//
//template<typename T>
//void bakproj_AIM_temp(T* dimg, T* dprj, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_bakproj_AIM_ker<T> << <gid, blk >> >(dimg, dprj, FanGeo, Img);
//}
//void bakproj_AIM(float* dimg, float* dprj, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	bakproj_AIM_temp<float>(dimg, dprj, FanGeo, Img, blk, gid);
//}
//
//void bakproj_AIM(double* dimg, double* dprj, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	bakproj_AIM_temp<double>(dimg, dprj, FanGeo, Img, blk, gid);
//}
//
//
//
//
//
//template<typename T>
//__global__ void bakproj_AIM_ker(T*dimg, T* dprj, const unsigned int angIdx, const FanEAGeo FanGeo, const Image Img)
//{
//	cuint xi = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint yi = threadIdx.y + blockIdx.y * blockDim.y;
//	if (xi < Img.m_Reso.x && yi < Img.m_Reso.y)
//	{
//		const T cosT = cos(FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp);
//		const T sinT = sin(FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp);
//		const T cntImgX = static_cast<T>(Img.m_Reso.x - 1.0) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
//		const T cntImgY = static_cast<T>(Img.m_Reso.y - 1.0) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
//		const T area = Img.m_Step.x * Img.m_Step.y;
//
//
//		T sour[2];
//		sour[0] = -FanGeo.m_S2O * sinT;
//		sour[1] = FanGeo.m_S2O * cosT;
//
//
//		T grid[4][3];
//		grid[0][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
//		grid[0][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
//
//		grid[1][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
//		grid[1][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
//
//		grid[2][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
//		grid[2][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
//
//		grid[3][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
//		grid[3][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
//
//		T det[4][3];
//		//T dir[4][3];
//		//T lengths = 0;
//
//		det[0][0] = grid[0][0] - sour[0]; //used as the direction from source to the corner of the pixel
//		det[0][1] = grid[0][1] - sour[1];
//		det[1][0] = grid[1][0] - sour[0];
//		det[1][1] = grid[1][1] - sour[1];
//		det[2][0] = grid[2][0] - sour[0];
//		det[2][1] = grid[2][1] - sour[1];
//		det[3][0] = grid[3][0] - sour[0];
//		det[3][1] = grid[3][1] - sour[1];
//
//		T tmp[2];
//		tmp[1] = sqrt(det[0][0] * det[0][0] + det[0][1] * det[0][1]);
//		det[0][0] /= tmp[1];
//		det[0][1] /= tmp[1];
//		tmp[1] = sqrt(det[1][0] * det[1][0] + det[1][1] * det[1][1]);
//		det[1][0] /= tmp[1];
//		det[1][1] /= tmp[1];
//		tmp[1] = sqrt(det[2][0] * det[2][0] + det[2][1] * det[2][1]);
//		det[2][0] /= tmp[1];
//		det[2][1] /= tmp[1];
//		tmp[1] = sqrt(det[3][0] * det[3][0] + det[3][1] * det[3][1]);
//		det[3][0] /= tmp[1];
//		det[3][1] /= tmp[1];
//
//
//		//T initDet[4][3];
//
//		////���㵱ǰ���Ӧ���ĸ�̽����λ�����
//		det[0][0] = sour[0] + det[0][0] * FanGeo.m_S2D;  //det used as the detector coordinates
//		det[0][1] = sour[1] + det[0][1] * FanGeo.m_S2D;
//		det[1][0] = sour[0] + det[1][0] * FanGeo.m_S2D;
//		det[1][1] = sour[1] + det[1][1] * FanGeo.m_S2D;
//		det[2][0] = sour[0] + det[2][0] * FanGeo.m_S2D;
//		det[2][1] = sour[1] + det[2][1] * FanGeo.m_S2D;
//		det[3][0] = sour[0] + det[3][0] * FanGeo.m_S2D;
//		det[3][1] = sour[1] + det[3][1] * FanGeo.m_S2D;
//
//		tmp[0] = det[0][0] * cosT + det[0][1] * sinT;
//		tmp[1] = -det[0][0] * sinT + det[0][1] * cosT;
//		det[0][0] = tmp[0];
//		det[0][1] = tmp[1];
//		tmp[0] = det[1][0] * cosT + det[1][1] * sinT;
//		tmp[1] = -det[1][0] * sinT + det[1][1] * cosT;
//		det[1][0] = tmp[0];
//		det[1][1] = tmp[1];
//		tmp[0] = det[2][0] * cosT + det[2][1] * sinT;
//		tmp[1] = -det[2][0] * sinT + det[2][1] * cosT;
//		det[2][0] = tmp[0];
//		det[2][1] = tmp[1];
//		tmp[0] = det[3][0] * cosT + det[3][1] * sinT;
//		tmp[1] = -det[3][0] * sinT + det[3][1] * cosT;
//		det[3][0] = tmp[0];
//		det[3][1] = tmp[1];
//
//		T beta[4];
//		beta[0] = atan(-det[0][0] / (det[0][1] - FanGeo.m_S2O));
//		beta[1] = atan(-det[1][0] / (det[1][1] - FanGeo.m_S2O));
//		beta[2] = atan(-det[2][0] / (det[2][1] - FanGeo.m_S2O));
//		beta[3] = atan(-det[3][0] / (det[3][1] - FanGeo.m_S2O));
//
//		grid[0][2] = beta[0] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//		grid[1][2] = beta[1] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//		grid[2][2] = beta[2] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//		grid[3][2] = beta[3] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//		SortProjection<T>(grid);
//
//		int minDetIdx = static_cast<int>(grid[0][2]);
//		int maxDetIdx(int(ceil(grid[3][2])));
//		minDetIdx = (minDetIdx < 0) ? 0 : minDetIdx;
//		maxDetIdx = (maxDetIdx > static_cast<int>(FanGeo.m_DetN)) ? static_cast<int>(FanGeo.m_DetN) : maxDetIdx;
//
//		tmp[1] = (hypot((xi - cntImgX)*Img.m_Step.x - sour[0], (yi - cntImgY)*Img.m_Step.y - sour[1]));
//		//T summ = 0;
//		//T pangle;
//		//T direct[2], legth, SVA[3], SVB[3], coef;
//		T SVA[3], SVB[3];
//
//		//�����߱߼���;
//		beta[3] = (minDetIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp; //���С�߱�ƫ�ƽ� // Beta[3] used as the pangle
//		beta[0] = sin(beta[3]) * FanGeo.m_S2D; // beta used as the initial X,Y and current X,Y with rotation
//		beta[1] = -cos(beta[3]) * FanGeo.m_S2D + FanGeo.m_S2O;
//
//		beta[2] = beta[0] * cosT - beta[1] * sinT;
//		beta[3] = beta[0] * sinT + beta[1] * cosT;
//		beta[0] = beta[2] - sour[0]; //Beta used as the directions
//		beta[1] = beta[3] - sour[1];
//		beta[2] = sqrt(beta[0] * beta[0] + beta[1] * beta[1]);//Beta[3] used as the length of the image
//		SVA[0] = beta[0] / beta[2];
//		SVA[1] = beta[1] / beta[2];
//		SVA[2] = static_cast<T>(minDetIdx);
//		tmp[0] = 0;
//		for (; minDetIdx < maxDetIdx; ++minDetIdx)
//		{
//			//�����߱߼���;
//
//			beta[3] = (minDetIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp + FanGeo.m_DetStp; //��߱�ƫ�ƽ�
//			beta[0] = sin(beta[3]) * FanGeo.m_S2D;
//			beta[1] = -cos(beta[3]) * FanGeo.m_S2D + FanGeo.m_S2O;
//			beta[2] = beta[0] * cosT - beta[1] * sinT;
//			beta[3] = beta[0] * sinT + beta[1] * cosT;
//			beta[0] = beta[2] - sour[0];
//			beta[1] = beta[3] - sour[1];
//			beta[2] = sqrt(beta[0] * beta[0] + beta[1] * beta[1]);
//			SVB[0] = beta[0] / beta[2];
//			SVB[1] = beta[1] / beta[2];
//			SVB[2] = static_cast<T>(minDetIdx)+1;
//
//			//Compute the weighting coefficient for a special projection data
//			beta[3] = ComputeCoefficient<T>(grid, SVA, SVB, sour, area); // beta[3] used as the coef
//			beta[3] /= (tmp[1] * fabs(FanGeo.m_DetStp));
//			tmp[0] += dprj[angIdx* FanGeo.m_DetN + minDetIdx] * beta[3];
//
//			SVA[0] = SVB[0]; //��߱߱��С�߱�
//			SVA[1] = SVB[1];
//			SVA[2] = SVB[2];
//		}
//		dimg[yi * Img.m_Reso.x + xi] = tmp[0];
//	}
//}
//
//
//template<typename T>
//void bakproj_AIM_temp(T* dimg, T *dprj, const unsigned int angIdx, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3&gid)
//{
//	bakproj_AIM_ker<T> << < gid, blk >> >(dimg, dprj, angIdx, FanGeo, Img);
//}
//
//void bakproj_AIM(float* dimg, float* dprj, const unsigned int angIdx, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	bakproj_AIM_temp<float>(dimg, dprj, angIdx, FanGeo, Img, blk, gid);
//}
//
//void bakproj_AIM(double* dimg, double* dprj, const unsigned int angIdx, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	bakproj_AIM_temp<double>(dimg, dprj, angIdx, FanGeo, Img, blk, gid);
//}
//
//
//
//
//
//
//
//template<typename T>
//__global__ void _bakproj_AIM_ker(T* dprj, T* dimg, const FanEAGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint xi = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint yi = threadIdx.y + blockIdx.y * blockDim.y;
//
//	if (xi < Img.m_Reso.x&&yi < Img.m_Reso.y)
//	{
//		const T cntImgX = static_cast<T>(Img.m_Reso.x - 1.0) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
//		const T cntImgY = static_cast<T>(Img.m_Reso.y - 1.0) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
//		const T area = Img.m_Step.x * Img.m_Step.y;
//		T sour[2], cosT, sinT;
//		T grid[4][3];
//		grid[0][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
//		grid[0][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
//
//		grid[1][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
//		grid[1][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
//
//		grid[2][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
//		grid[2][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
//
//		grid[3][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
//		grid[3][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
//		T det[4][3];
//		T tmp[2];
//		T beta[4];
//		int minDetIdx, maxDetIdx;
//		T summ = 0;
//		T SVA[3], SVB[3];
//		size_t angIdx;
//		for (size_t prjIdx = 0; prjIdx < numPerSubSet; ++angIdx)
//		{
//			angIdx = prjIdx * subSetNum + curSubSetIdx;
//			cosT = cos(FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp);
//			sinT = sin(FanGeo.m_ViwBeg + angIdx * FanGeo.m_ViwStp);
//
//			sour[0] = -FanGeo.m_S2O * sinT;
//			sour[1] = FanGeo.m_S2O * cosT;
//
//			det[0][0] = grid[0][0] - sour[0]; //used as the direction from source to the corner of the pixel
//			det[0][1] = grid[0][1] - sour[1];
//			det[1][0] = grid[1][0] - sour[0];
//			det[1][1] = grid[1][1] - sour[1];
//			det[2][0] = grid[2][0] - sour[0];
//			det[2][1] = grid[2][1] - sour[1];
//			det[3][0] = grid[3][0] - sour[0];
//			det[3][1] = grid[3][1] - sour[1];
//
//
//			tmp[1] = sqrt(det[0][0] * det[0][0] + det[0][1] * det[0][1]);
//			det[0][0] /= tmp[1];
//			det[0][1] /= tmp[1];
//			tmp[1] = sqrt(det[1][0] * det[1][0] + det[1][1] * det[1][1]);
//			det[1][0] /= tmp[1];
//			det[1][1] /= tmp[1];
//			tmp[1] = sqrt(det[2][0] * det[2][0] + det[2][1] * det[2][1]);
//			det[2][0] /= tmp[1];
//			det[2][1] /= tmp[1];
//			tmp[1] = sqrt(det[3][0] * det[3][0] + det[3][1] * det[3][1]);
//			det[3][0] /= tmp[1];
//			det[3][1] /= tmp[1];
//
//			det[0][0] = sour[0] + det[0][0] * FanGeo.m_S2D;  //det used as the detector coordinates
//			det[0][1] = sour[1] + det[0][1] * FanGeo.m_S2D;
//			det[1][0] = sour[0] + det[1][0] * FanGeo.m_S2D;
//			det[1][1] = sour[1] + det[1][1] * FanGeo.m_S2D;
//			det[2][0] = sour[0] + det[2][0] * FanGeo.m_S2D;
//			det[2][1] = sour[1] + det[2][1] * FanGeo.m_S2D;
//			det[3][0] = sour[0] + det[3][0] * FanGeo.m_S2D;
//			det[3][1] = sour[1] + det[3][1] * FanGeo.m_S2D;
//
//			tmp[0] = det[0][0] * cosT + det[0][1] * sinT;
//			tmp[1] = -det[0][0] * sinT + det[0][1] * cosT;
//			det[0][0] = tmp[0];
//			det[0][1] = tmp[1];
//			tmp[0] = det[1][0] * cosT + det[1][1] * sinT;
//			tmp[1] = -det[1][0] * sinT + det[1][1] * cosT;
//			det[1][0] = tmp[0];
//			det[1][1] = tmp[1];
//			tmp[0] = det[2][0] * cosT + det[2][1] * sinT;
//			tmp[1] = -det[2][0] * sinT + det[2][1] * cosT;
//			det[2][0] = tmp[0];
//			det[2][1] = tmp[1];
//			tmp[0] = det[3][0] * cosT + det[3][1] * sinT;
//			tmp[1] = -det[3][0] * sinT + det[3][1] * cosT;
//			det[3][0] = tmp[0];
//			det[3][1] = tmp[1];
//
//
//			beta[0] = atan(-det[0][0] / (det[0][1] - FanGeo.m_S2O));
//			beta[1] = atan(-det[1][0] / (det[1][1] - FanGeo.m_S2O));
//			beta[2] = atan(-det[2][0] / (det[2][1] - FanGeo.m_S2O));
//			beta[3] = atan(-det[3][0] / (det[3][1] - FanGeo.m_S2O));
//
//			grid[0][2] = beta[0] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			grid[1][2] = beta[1] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			grid[2][2] = beta[2] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			grid[3][2] = beta[3] / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//			SortProjection<T>(grid);
//			minDetIdx = static_cast<int>(grid[0][2]);
//			maxDetIdx = (int(ceil(grid[3][2])));
//			minDetIdx = (minDetIdx < 0) ? 0 : minDetIdx;
//			maxDetIdx = (maxDetIdx > static_cast<int>(FanGeo.m_DetN)) ? static_cast<int>(FanGeo.m_DetN) : maxDetIdx;
//
//			tmp[1] = (hypot((xi - cntImgX)*Img.m_Step.x - sour[0], (yi - cntImgY)*Img.m_Step.y - sour[1]));
//			//T summ = 0;
//			//T pangle;
//			//T direct[2], legth, SVA[3], SVB[3], coef;
//
//
//			//�����߱߼���;
//			beta[3] = (minDetIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp; //���С�߱�ƫ�ƽ� // Beta[3] used as the pangle
//			beta[0] = sin(beta[3]) * FanGeo.m_S2D; // beta used as the initial X,Y and current X,Y with rotation
//			beta[1] = -cos(beta[3]) * FanGeo.m_S2D + FanGeo.m_S2O;
//
//			beta[2] = beta[0] * cosT - beta[1] * sinT;
//			beta[3] = beta[0] * sinT + beta[1] * cosT;
//			beta[0] = beta[2] - sour[0]; //Beta used as the directions
//			beta[1] = beta[3] - sour[1];
//			beta[2] = sqrt(beta[0] * beta[0] + beta[1] * beta[1]);//Beta[3] used as the length of the image
//			SVA[0] = beta[0] / beta[2];
//			SVA[1] = beta[1] / beta[2];
//			SVA[2] = static_cast<T>(minDetIdx);
//			tmp[0] = 0;
//			for (; minDetIdx < maxDetIdx; ++minDetIdx)
//			{
//				//�����߱߼���;
//
//				beta[3] = (minDetIdx - FanGeo.m_DetCntIdx) * FanGeo.m_DetStp + FanGeo.m_DetStp; //��߱�ƫ�ƽ�
//				beta[0] = sin(beta[3]) * FanGeo.m_S2D;
//				beta[1] = -cos(beta[3]) * FanGeo.m_S2D + FanGeo.m_S2O;
//				beta[2] = beta[0] * cosT - beta[1] * sinT;
//				beta[3] = beta[0] * sinT + beta[1] * cosT;
//				beta[0] = beta[2] - sour[0];
//				beta[1] = beta[3] - sour[1];
//				beta[2] = sqrt(beta[0] * beta[0] + beta[1] * beta[1]);
//				SVB[0] = beta[0] / beta[2];
//				SVB[1] = beta[1] / beta[2];
//				SVB[2] = static_cast<T>(minDetIdx)+1;
//
//				//Compute the weighting coefficient for a special projection data
//				beta[3] = ComputeCoefficient<T>(grid, SVA, SVB, sour, area); // beta[3] used as the coef
//				beta[3] /= (tmp[1] * fabs(FanGeo.m_DetStp));
//				tmp[0] += dprj[prjIdx* FanGeo.m_DetN + minDetIdx] * beta[3];
//
//				SVA[0] = SVB[0]; //��߱߱��С�߱�
//				SVA[1] = SVB[1];
//				SVA[2] = SVB[2];
//			}
//			summ += tmp[0];
//		}
//		dimg[yi* Img.m_Reso.x + xi] = summ;
//	}
//}
//
//
//
//
//template<typename T>
//void bakproj_AIM_temp(T* donp, T* dimg, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3&gid)
//{
//	_bakproj_AIM_ker<T> << < gid, blk >> >(donp, dimg, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx);
//}
//void bakproj_AIM(float* donp, float* dimg, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	bakproj_AIM_temp<float>(donp, dimg, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx, blk, gid);
//}
//void bakproj_AIM(double* donp, double* dimg, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	bakproj_AIM_temp<double>(donp, dimg, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx, blk, gid);
//}
//
//
//
//
//template<typename T>
//__global__ void _bakproj_AIM_ker(T* dprj, T* dimg, const FanEDGeo FanGeo, const Image Img)
//{
//	cuint xi = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint yi = threadIdx.y + blockDim.y * blockIdx.y;
//	if (xi < Img.m_Reso.x && yi < Img.m_Reso.y)
//	{
//		T cosT;
//		T sinT;
//		const T cntImgX = static_cast<T>(Img.m_Reso.x - 1.0) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
//		const T cntImgY = static_cast<T>(Img.m_Reso.y - 1.0) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
//		const T area = Img.m_Step.x * Img.m_Step.y;
//
//		T sour[2];
//		T grid[4][3];
//
//		int minDetIdx, maxDetIdx;
//		T pdist;
//		T SVA[3], SVB[3];
//		T vv;
//		T coef;
//		T summ = 0;
//		int angidx = 0;
//		for (angidx = 0; angidx != FanGeo.m_ViwN; ++angidx)
//		{
//			cosT = cos(FanGeo.m_ViwBeg + angidx * FanGeo.m_ViwStp);
//			sinT = sin(FanGeo.m_ViwBeg + angidx * FanGeo.m_ViwStp);
//
//			sour[0] = -FanGeo.m_S2O * sinT;
//			sour[1] = FanGeo.m_S2O * cosT;
//
//			grid[0][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
//			grid[0][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
//			pdist = grid[0][0] * cosT + grid[0][1] * sinT;
//			vv = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
//			grid[0][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//			grid[1][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
//			grid[1][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
//			pdist = grid[1][0] * cosT + grid[1][1] * sinT;
//			vv = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
//			grid[1][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//			grid[2][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
//			grid[2][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
//			pdist = grid[2][0] * cosT + grid[2][1] * sinT;
//			vv = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
//			grid[2][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//			grid[3][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
//			grid[3][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
//			pdist = grid[3][0] * cosT + grid[3][1] * sinT;
//			vv = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
//			grid[3][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			//����
//			SortProjection<T>(grid);
//
//			minDetIdx = static_cast<int>(grid[0][2]);
//			maxDetIdx = (int(ceil(grid[3][2])));
//			minDetIdx = (minDetIdx < 0) ? 0 : minDetIdx;
//			maxDetIdx = (maxDetIdx > static_cast<int>(FanGeo.m_DetN)) ? static_cast<int>(FanGeo.m_DetN) : maxDetIdx;
//			pdist = (hypot((xi - cntImgX)*Img.m_Step.x - sour[0], (yi - cntImgY)*Img.m_Step.y - sour[1]));
//
//			for (; minDetIdx < maxDetIdx; ++minDetIdx)
//			{
//				calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, minDetIdx);
//				vv = acos(abs(SVA[0] * SVB[0] + SVA[1] * SVB[1]));
//				coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//				coef = coef / (pdist * vv);
//				summ += dprj[angidx * FanGeo.m_DetN + minDetIdx] * coef;
//			}
//		}
//		dimg[yi* Img.m_Reso.x + xi] = summ;
//	}
//}
//
//template<typename T>
//void bakproj_AIM_temp(T* dprj, T* dimg, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_bakproj_AIM_ker<T> << <gid, blk >> >(dprj, dimg, FanGeo, Img);
//}
//
//void bakproj_AIM(float* dprj, float* dimg, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	bakproj_AIM_temp<float>(dprj, dimg, FanGeo, Img, blk, gid);
//}
//
//void bakproj_AIM(double* dprj, double* dimg, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	bakproj_AIM_temp<double>(dprj, dimg, FanGeo, Img, blk, gid);
//}
//
//
//
//
//
//
//template<typename T>
//__global__ void _bakproj_AIM_ker(T* dprj, T* dimg, cuint angidx, const FanEDGeo FanGeo, const Image Img)
//{
//	cuint xi = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint yi = threadIdx.y + blockDim.y * blockIdx.y;
//	if (xi < Img.m_Reso.x && yi < Img.m_Reso.y)
//	{
//		T cosT;
//		T sinT;
//		const T cntImgX = static_cast<T>(Img.m_Reso.x - 1.0) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
//		const T cntImgY = static_cast<T>(Img.m_Reso.y - 1.0) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
//		const T area = Img.m_Step.x * Img.m_Step.y;
//
//		T sour[2];
//		T grid[4][3];
//
//		int minDetIdx, maxDetIdx;
//		T pdist;
//		T SVA[3], SVB[3];
//		T vv;
//		T coef;
//		T summ = 0;
//
//		cosT = cos(FanGeo.m_ViwBeg + angidx * FanGeo.m_ViwStp);
//		sinT = sin(FanGeo.m_ViwBeg + angidx * FanGeo.m_ViwStp);
//
//		sour[0] = -FanGeo.m_S2O * sinT;
//		sour[1] = FanGeo.m_S2O * cosT;
//
//		grid[0][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
//		grid[0][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
//		pdist = grid[0][0] * cosT + grid[0][1] * sinT;
//		vv = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
//		grid[0][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//		grid[1][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
//		grid[1][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
//		pdist = grid[1][0] * cosT + grid[1][1] * sinT;
//		vv = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
//		grid[1][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//		grid[2][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
//		grid[2][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
//		pdist = grid[2][0] * cosT + grid[2][1] * sinT;
//		vv = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
//		grid[2][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//		grid[3][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
//		grid[3][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
//		pdist = grid[3][0] * cosT + grid[3][1] * sinT;
//		vv = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
//		grid[3][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//		//����
//		SortProjection<T>(grid);
//
//		minDetIdx = static_cast<int>(grid[0][2]);
//		maxDetIdx = (int(ceil(grid[3][2])));
//		minDetIdx = (minDetIdx < 0) ? 0 : minDetIdx;
//		maxDetIdx = (maxDetIdx > static_cast<int>(FanGeo.m_DetN)) ? static_cast<int>(FanGeo.m_DetN) : maxDetIdx;
//		pdist = (hypot((xi - cntImgX)*Img.m_Step.x - sour[0], (yi - cntImgY)*Img.m_Step.y - sour[1]));
//
//		for (; minDetIdx < maxDetIdx; ++minDetIdx)
//		{
//			calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, minDetIdx);
//			vv = acos(abs(SVA[0] * SVB[0] + SVA[1] * SVB[1]));
//			coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//			coef = coef / (pdist * vv);
//			summ += dprj[angidx * FanGeo.m_DetN + minDetIdx] * coef;
//		}
//
//		dimg[yi* Img.m_Reso.x + xi] = summ;
//	}
//}
//
//template<typename T>
//void bakproj_AIM_temp(T* dprj, T* dimg, cuint angidx, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_bakproj_AIM_ker<T> << <gid, blk >> >(dprj, dimg, angidx, FanGeo, Img);
//}
//void bakproj_AIM(float* dimg, float* dprj, const unsigned int angIdx, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	bakproj_AIM_temp(dprj, dimg, angIdx, FanGeo, Img, blk, gid);
//}
//void bakproj_AIM(double* dimg, double* dprj, const unsigned int angIdx, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	bakproj_AIM_temp(dprj, dimg, angIdx, FanGeo, Img, blk, gid);
//}
//
//
//template<typename T>
//__global__ void _bakproj_AIM_ker(T* dprj, T* dimg, const FanEDGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint xi = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint yi = threadIdx.y + blockDim.y * blockIdx.y;
//	if (xi < Img.m_Reso.x && yi < Img.m_Reso.y)
//	{
//		T cosT;
//		T sinT;
//		const T cntImgX = static_cast<T>(Img.m_Reso.x - 1.0) * 0.5 + (Img.m_Bias.x / Img.m_Step.x);
//		const T cntImgY = static_cast<T>(Img.m_Reso.y - 1.0) * 0.5 + (Img.m_Bias.y / Img.m_Step.y);
//		const T area = Img.m_Step.x * Img.m_Step.y;
//
//		T sour[2];
//		T grid[4][3];
//
//		int minDetIdx, maxDetIdx;
//		T pdist;
//		T SVA[3], SVB[3];
//		T vv;
//		T coef;
//		T summ = 0;
//		int angidx = 0;
//		int prjidx = 0;
//		for (prjidx = 0; prjidx < numPerSubSet; ++prjidx)
//		{
//			angidx = prjidx * subSetNum + curSubSetIdx;
//			cosT = cos(FanGeo.m_ViwBeg + angidx * FanGeo.m_ViwStp);
//			sinT = sin(FanGeo.m_ViwBeg + angidx * FanGeo.m_ViwStp);
//
//			sour[0] = -FanGeo.m_S2O * sinT;
//			sour[1] = FanGeo.m_S2O * cosT;
//
//			grid[0][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
//			grid[0][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
//			pdist = grid[0][0] * cosT + grid[0][1] * sinT;
//			vv = -grid[0][0] * sinT + grid[0][1] * cosT - FanGeo.m_S2O;
//			grid[0][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//			grid[1][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
//			grid[1][1] = (yi - cntImgY - 0.5f) * Img.m_Step.y;
//			pdist = grid[1][0] * cosT + grid[1][1] * sinT;
//			vv = -grid[1][0] * sinT + grid[1][1] * cosT - FanGeo.m_S2O;
//			grid[1][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//			grid[2][0] = (xi - cntImgX + 0.5f) * Img.m_Step.x;
//			grid[2][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
//			pdist = grid[2][0] * cosT + grid[2][1] * sinT;
//			vv = -grid[2][0] * sinT + grid[2][1] * cosT - FanGeo.m_S2O;
//			grid[2][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//
//			grid[3][0] = (xi - cntImgX - 0.5f) * Img.m_Step.x;
//			grid[3][1] = (yi - cntImgY + 0.5f) * Img.m_Step.y;
//			pdist = grid[3][0] * cosT + grid[3][1] * sinT;
//			vv = -grid[3][0] * sinT + grid[3][1] * cosT - FanGeo.m_S2O;
//			grid[3][2] = (-(pdist * FanGeo.m_S2D) / vv) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			//����
//			SortProjection<T>(grid);
//
//			minDetIdx = static_cast<int>(grid[0][2]);
//			maxDetIdx = (int(ceil(grid[3][2])));
//			minDetIdx = (minDetIdx < 0) ? 0 : minDetIdx;
//			maxDetIdx = (maxDetIdx > static_cast<int>(FanGeo.m_DetN)) ? static_cast<int>(FanGeo.m_DetN) : maxDetIdx;
//			pdist = (hypot((xi - cntImgX)*Img.m_Step.x - sour[0], (yi - cntImgY)*Img.m_Step.y - sour[1]));
//
//			for (; minDetIdx < maxDetIdx; ++minDetIdx)
//			{
//				calSVASVB<T>(SVA, SVB, sour, cosT, sinT, FanGeo, Img, minDetIdx);
//				vv = acos(abs(SVA[0] * SVB[0] + SVA[1] * SVB[1]));
//				coef = ComputeCoefficient<T>(grid, SVA, SVB, sour, area);
//				coef = coef / (pdist * vv);
//				summ += dprj[prjidx * FanGeo.m_DetN + minDetIdx] * coef;
//			}
//		}
//		dimg[yi* Img.m_Reso.x + xi] = summ;
//	}
//}
//template<typename T>
//void bakproj_AIM_temp(T* donp, T*dimg, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_bakproj_AIM_ker<T> << <gid, blk >> >(donp, dimg, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx);
//}
//
//void bakproj_AIM(float* donp, float* dimg, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	bakproj_AIM_temp<float>(donp, dimg, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx, blk, gid);
//}
//void bakproj_AIM(double* donp, double* dimg, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	bakproj_AIM_temp<double>(donp, dimg, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx, blk, gid);
//}
//
//
//
//__global__ void _proj_Ker_DEMO17(float* dvol, float* dprj, float* draw, const ConeEDGeo ConeGeo, const Volume Vol, const int angIdx)
//{
//	cuint detIdx = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint detIdz = threadIdx.y + blockIdx.y * blockDim.y;
//	//cuint angIdx = threadIdx.z + blockIdx.z * blockDim.z;
//	if (detIdx < ConeGeo.m_DetN.x && detIdz < ConeGeo.m_DetN.y)
//	{
//		float3 MINO = make_float3(
//			-Vol.m_Size.x * 0.5f + Vol.m_Bias.x,
//			-Vol.m_Size.y * 0.5f + Vol.m_Bias.y,
//			-Vol.m_Size.z * 0.5f + Vol.m_Bias.z);
//
//		//current rotation
//		float cosT = cos(ConeGeo.m_ViwBeg + (float)angIdx * ConeGeo.m_ViwStp);
//		float sinT = sin(ConeGeo.m_ViwBeg + (float)angIdx * ConeGeo.m_ViwStp);
//
//		Ray ray;
//		ray.o = rotation(make_float3(0, ConeGeo.m_S2O, 0), cosT, sinT);
//		float3 curDet = rotation(make_float3(
//			((float)detIdx - ConeGeo.m_DetCntIdx.x + 0.5f) * ConeGeo.m_DetStp.x,
//			-ConeGeo.m_O2D,
//			((float)detIdz - ConeGeo.m_DetCntIdx.y + 0.5f) * ConeGeo.m_DetStp.y), cosT, sinT);
//		ray.d = normalize(curDet - ray.o);
//		float totWeight(0.0f);
//		float res = draw[detIdz * ConeGeo.m_DetN.x + detIdx] - //һ��ֻloadһ���Ƕ�;
//			calSiddonOneRayKer(
//			ray.o.x, ray.o.y, ray.o.z,
//			curDet.x, curDet.y, curDet.z,
//			MINO.x, MINO.y, MINO.z,
//			Vol.m_Step.x, Vol.m_Step.y, Vol.m_Step.z,
//			Vol.m_Reso.x, Vol.m_Reso.y, Vol.m_Reso.z, dvol, &totWeight);
//		if (!IS_ZERO(totWeight))
//		{
//			dprj[detIdz * ConeGeo.m_DetN.x + detIdx] = res / totWeight;
//		}
//		else
//		{
//			dprj[detIdz * ConeGeo.m_DetN.x + detIdx] = 0;
//		}
//	}
//}
//void proj_DEMO17(float* dvol, float* dprj, float* draw, const ConeEDGeo& ConeGeo, const Volume& Vol, const int angIdx, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker_DEMO17 << <gid, blk >> >(dvol, dprj, draw, ConeGeo, Vol, angIdx);
//}
//
//__global__ void backProj_Ker_DEMO17(float* correction, float* volume, const ConeEDGeo ConeGeo, const Volume Vol, const int angIdx)
//{
//	//cuint index = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
//
//	int3 idx;
//	idx.x = threadIdx.x + blockIdx.x * blockDim.x;
//	idx.y = threadIdx.y + blockIdx.y * blockDim.y;
//	idx.z = threadIdx.z + blockIdx.z * blockDim.z;
//	if (idx.x < Vol.m_Reso.x && idx.y < Vol.m_Reso.y && idx.z < Vol.m_Reso.z)
//	{
//		int index = (idx.z * Vol.m_Reso.y + idx.y) * Vol.m_Reso.x + idx.x;
//
//		//cuint volReso = Vol.m_Reso.x * Vol.m_Reso.y * Vol.m_Reso.z;
//		const float cosT = cos(ConeGeo.m_ViwBeg + angIdx * ConeGeo.m_ViwStp);
//		const float sinT = sin(ConeGeo.m_ViwBeg + angIdx * ConeGeo.m_ViwStp);
//
//		float minX = 10000, minZ = 10000, maxX = -10000, maxZ = -10000;
//
//		float3 projection, origPos, p, curP;
//		int i(0), j(0), k(0);
//		//float tt;
//		Ray eyeRay;
//
//		eyeRay.o = make_float3(-ConeGeo.m_S2O * sinT, ConeGeo.m_S2O * cosT, 0);
//
//		float tnear = 0.0;
//		float tfar = 0.0;
//		float3 boxMin, boxMax;
//
//		boxMin.x = -Vol.m_Size.x * 0.5 + idx.x * Vol.m_Step.x;
//		boxMin.y = -Vol.m_Size.y * 0.5 + idx.y * Vol.m_Step.y;
//		boxMin.z = -Vol.m_Size.z * 0.5 + idx.z * Vol.m_Step.z;
//		boxMax.x = boxMin.x + Vol.m_Step.x;
//		boxMax.y = boxMin.y + Vol.m_Step.y;
//		boxMax.z = boxMin.z + Vol.m_Step.z;
//		float __MINDL__ = -ConeGeo.m_DetCntIdx.x * ConeGeo.m_DetStp.x;
//		float __MINDH__ = -ConeGeo.m_DetCntIdx.y * ConeGeo.m_DetStp.y;
//		for (i = 0; i != 2; ++i)
//		{
//			for (j = 0; j != 2; ++j)
//			{
//				for (k = 0; k != 2; ++k)
//				{
//
//					origPos.x = -Vol.m_Size.x * 0.5 + (idx.x + i) * Vol.m_Step.x;
//					origPos.y = -Vol.m_Size.y * 0.5 + (idx.y + j) * Vol.m_Step.y;
//					origPos.z = -Vol.m_Size.z * 0.5 + (idx.z + k) * Vol.m_Step.z;
//
//					//Trans from the world coordinate to detector index coordinate;
//					projection.x = (-ConeGeo.m_S2D * cosT / ConeGeo.m_DetStp.x + __MINDL__ * sinT / ConeGeo.m_DetStp.x) * origPos.x +
//						(-ConeGeo.m_S2D * sinT / ConeGeo.m_DetStp.x - __MINDL__ * cosT / ConeGeo.m_DetStp.x) * origPos.y +
//						0 * origPos.z + __MINDL__ * ConeGeo.m_S2O / ConeGeo.m_DetStp.x;
//					projection.y = __MINDH__ / ConeGeo.m_DetStp.y * sinT * origPos.x + (-__MINDH__ / ConeGeo.m_DetStp.y) * cosT * origPos.y +
//						(-ConeGeo.m_S2D / ConeGeo.m_DetStp.y) * origPos.z + __MINDH__ * ConeGeo.m_S2O / ConeGeo.m_DetStp.y;
//					projection.z = -sinT * origPos.x + cosT * origPos.y - ConeGeo.m_S2O;
//
//					projection.x = projection.x / projection.z;
//					projection.y = projection.y / projection.z;
//					//projection.z = projection.z / projection.z;
//
//					//Obtain the bounding box in the projection plane;
//
//					minX = min(minX, projection.x);
//					maxX = max(maxX, projection.x);
//
//					minZ = min(minZ, projection.y);
//					maxZ = max(maxZ, projection.y);
//				}
//			}
//		}
//
//		if (minX != 10000)
//			minX = (minX < 0) ? 0 : minX;
//		if (maxX != -10000)
//			maxX = (maxX >= ConeGeo.m_DetN.x) ? ConeGeo.m_DetN.x : maxX;
//
//		if (minZ != 10000)
//			minZ = (minZ < 0) ? 0 : minZ;
//		if (maxZ != -10000)
//			maxZ = (maxZ >= ConeGeo.m_DetN.y) ? ConeGeo.m_DetN.y : maxZ;
//
//		int2 mbrMin, mbrMax;
//
//		mbrMin.x = minX;
//		mbrMin.y = minZ;
//		mbrMax.x = maxX;
//		mbrMax.y = maxZ;
//
//
//		int x(0), y(0);
//		bool hit(0);
//		float sumWeight(0), sumError(0);
//		for (x = mbrMin.x; x <= mbrMax.x; ++x)
//		{
//			for (y = mbrMin.y; y <= mbrMax.y; ++y)
//			{
//				p = make_float3(
//					__MINDL__ + (x + 0.5) * ConeGeo.m_DetStp.x,
//					-ConeGeo.m_O2D,
//					__MINDH__ + (y + 0.5) * ConeGeo.m_DetStp.y);
//				curP = rotation(p, cosT, sinT);
//				//curP = rotWithConstantSinCos(p, angIdx);
//
//				eyeRay.d = curP - eyeRay.o;
//				eyeRay.d = normalize(eyeRay.d);
//
//				hit = intersectBox(eyeRay, boxMin, boxMax, &tnear, &tfar);
//				if (hit)
//				{
//					sumWeight += (tfar - tnear);
//					sumError += (tfar - tnear) * correction[y * ConeGeo.m_DetN.x + x];
//				}
//
//			}
//		}
//		//volume[index] = 1.0;
//		if (!IS_ZERO(sumWeight))
//			volume[index] += sumError / sumWeight;
//	}
//}
//void back_proj_DEMO17(float* correction, float* volume, const ConeEDGeo& ConeGeo, const Volume& Vol, const int angIdx, const dim3& blockSize, const dim3& gridSize)
//{
//	backProj_Ker_DEMO17 << < gridSize, blockSize >> >(correction, volume, ConeGeo, Vol, angIdx);
//}
//
//
//
//
//__global__ void _proj_Ker(float* dimg, float* donp, float* draw, float* dpho, const float maxY, const FanEDGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint prjIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detId < FanGeo.m_DetN && prjIdx < numPerSubSet)
//	{
//
//		cuint angIdx = prjIdx * subSetNum + curSubSetIdx;
//
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//
//		//Current rotation angle;
//		float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//Current source position, assuming initial position is on the positive Y axis
//		Ray2D ray;
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//		float ang(0); //bias angle from the center of the Fan Beams
//
//		float2 curDetPos; //the current detector element position;
//		float totWeight(0);
//
//
//		// judging two indices addressing method, from large to small or small to large
//
//		ang = ((float)detId - FanGeo.m_DetCntIdx + 0.5f) * FanGeo.m_DetStp;
//		// current detector element position
//		//curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//		curDetPos = rotation(make_float2(ang, -FanGeo.m_O2D), cosT, sinT);
//		// x-ray direction;
//		ray.d = normalize(curDetPos - ray.o);
//
//		totWeight = 0;
//		float difff = draw[angIdx * FanGeo.m_DetN + detId] - calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
//			MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg, &totWeight);
//		if (!IS_ZERO(totWeight))
//		{
//			donp[prjIdx * FanGeo.m_DetN + detId] = difff / ((totWeight * maxY) / dpho[angIdx * FanGeo.m_DetN + detId]);
//		}
//		else
//		{
//			donp[prjIdx * FanGeo.m_DetN + detId] = 0;
//		}
//	}
//}
//void proj(float* dimg, float* donp, float* draw, float* dpho, const float maxY, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker << <gid, blk >> >(dimg, donp, draw, dpho, maxY, FanGeo, Img, numPerSubSet, subSetNum, curSubSetIdx);
//}
//
//
//
//__global__ void _proj_Ker_DEMO18v4_2D(float* dimg, float* draw, const FanEDGeo FanGeo, const Image Img)
//{
//	cuint detId = threadIdx.x + blockIdx.x * blockDim.x;
//	cuint prjIdx = threadIdx.y + blockIdx.y * blockDim.y;
//	if (detId < FanGeo.m_DetN && prjIdx < FanGeo.m_ViwN)
//	{
//		cuint angIdx = prjIdx;
//
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//
//		//Current rotation angle;
//		float curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//		float cosT = cosf(curAng);
//		float sinT = sinf(curAng);
//
//		//Current source position, assuming initial position is on the positive Y axis
//		Ray2D ray;
//		ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//
//		float ang(0); //bias angle from the center of the Fan Beams
//
//		float2 curDetPos; //the current detector element position;
//		float totWeight(0);
//
//
//		// judging two indices addressing method, from large to small or small to large
//
//		ang = ((float)detId - FanGeo.m_DetCntIdx + 0.5f) * FanGeo.m_DetStp;
//		// current detector element position
//		//curDetPos = rotation(make_float2(sinf(ang) * FanGeo.m_S2D, -cosf(ang) * FanGeo.m_S2D + FanGeo.m_S2O), cosT, sinT);
//		curDetPos = rotation(make_float2(ang, -FanGeo.m_O2D), cosT, sinT);
//		// x-ray direction;
//		ray.d = normalize(curDetPos - ray.o);
//
//		totWeight = 0;
//		draw[prjIdx * FanGeo.m_DetN + detId] = calSiddonOneRayKer2D(ray.o.x, ray.o.y, curDetPos.x, curDetPos.y,
//			MINO.x, MINO.y, Img.m_Step.x, Img.m_Step.y, Img.m_Reso.x, Img.m_Reso.y, dimg, &totWeight);
//	}
//}
//void proj_DEMO18v4_2D(float* dimg, float* draw, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_proj_Ker_DEMO18v4_2D << <gid, blk >> >(dimg, draw, FanGeo, Img);
//}
//
//
//__global__ void _bakproj_PIXEL_Ker_DEMO18v4_2D(float* donp, float* dimg, float* dmsk, const FanEDGeo FanGeo, const Image Img)//;(float* donp, float* dimg, const FanEDGeo FanGeo, const Image Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx)
//{
//	cuint idX = threadIdx.x + blockDim.x * blockIdx.x;
//	cuint idY = threadIdx.y + blockDim.y * blockIdx.y;
//	if (idX < Img.m_Reso.x && idY < Img.m_Reso.y)
//	{
//		float2 MINO = make_float2(-Img.m_Size.x * 0.5f + Img.m_Bias.x, -Img.m_Size.y * 0.5f + Img.m_Bias.y);
//		//cur rotation angle;
//		//unsigned int angIdx = 0;
//		float curAng(0), cosT(0), sinT(0), summ(0), mask(dmsk[idY * Img.m_Reso.x + idX]);
//		Ray2D ray;
//		float2 minBox, maxBox, curImg, curPix;
//		float tnear(0), tfar(0), weg(0), detID(0);
//		int flrID, ceilID;
//
//		for (unsigned int angIdx = 0; angIdx != FanGeo.m_ViwN; ++angIdx)
//		{
//			//angIdx = prjIdx * subSetNum + curSubSetIdx;
//			curAng = FanGeo.m_ViwBeg + FanGeo.m_ViwStp * angIdx;
//			cosT = cosf(curAng);
//			sinT = sinf(curAng);
//
//			//current pixel position
//			curPix = rotation(make_float2(
//				MINO.x + Img.m_Step.x * (idX + 0.5f),
//				MINO.y + Img.m_Step.y * (idY + 0.5f)), cosT, -sinT); //��ת���λ��;
//			//���λ����
//
//			//�����Ӧ��detidx;			//���ﲻ���� angDir; ���Ǵ�С�����;
//			detID = curPix.x / (FanGeo.m_S2O - curPix.y) * FanGeo.m_S2D / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;//FanGeo.m_S2D * curPix.x / (curPix.y - FanGeo.m_S2O) / FanGeo.m_DetStp + FanGeo.m_DetCntIdx;
//			if (detID < 0 || detID > FanGeo.m_DetN - 1)
//			{
//				dimg[idY * Img.m_Reso.x + idX] += 0.0f;
//				return;
//			}
//			flrID = detID;
//			ceilID = ceilf(detID);
//
//			ray.o = rotation(make_float2(0, FanGeo.m_S2O), cosT, sinT);
//			minBox = make_float2(
//				MINO.x + idX * Img.m_Step.x,
//				MINO.y + idY * Img.m_Step.y);
//			maxBox = minBox + make_float2(Img.m_Step.x, Img.m_Step.y);
//			curImg = minBox + 0.5f * make_float2(Img.m_Step.x, Img.m_Step.y);
//			ray.d = normalize(curImg - ray.o);
//
//			tnear = 0;
//			tfar = 0;
//			intersectBox(ray, minBox, maxBox, &tnear, &tfar); //�ཻ����;
//
//			weg = (tfar - tnear); //no weighting;
//			if (!IS_ZERO(flrID - ceilID))
//			{
//				//dimg[idY * Img.m_Reso.x + idX] += 1.0f;
//				summ += ((donp[angIdx * FanGeo.m_DetN + flrID] * (ceilID - detID) + donp[angIdx * FanGeo.m_DetN + ceilID] * (detID - flrID)) * weg);
//			}
//			else
//			{
//				//dimg[idY * Img.m_Reso.x + idX] += 100.0f;
//				summ += donp[angIdx * FanGeo.m_DetN + flrID] * weg;
//			}
//		}
//		dimg[idY * Img.m_Reso.x + idX] = (summ * mask);
//	}
//}
//void bakproj_PIXEL_DEMO18v4_2D(float* donp, float* dimg, float* dmsk, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid)
//{
//	_bakproj_PIXEL_Ker_DEMO18v4_2D << <gid, blk >> >(donp, dimg, dmsk, FanGeo, Img);
//}
