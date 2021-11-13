#include "Geometry.h"

#include <cmath>
#include <algorithm>
#ifndef M_PI
#define M_PI (3.14159265358979323846264)
#endif

Geometry::Geometry() :mDetShape(ARC)
{
}

Geometry::~Geometry()
{
}

Geometry::Geometry(const Geometry& rhs):
	mSourToObj(rhs.getSourToObj()),
	mSourToDet(rhs.getSourToDet()),
	mStartView(rhs.getStartView()),
	mViewPerRot(rhs.getViewPerRot()),
	mViewNumber(rhs.getViewNumber()),
	mViewStep(rhs.getViewStep()),
	mPitch(rhs.getPitch()),
	mZStep(rhs.getZStep()),
	mStartZ(rhs.getStartZ()),
	mVoxelSize(make_float3(rhs.getVoxelSizeX(), rhs.getVoxelSizeY(), rhs.getVoxelSizeZ())),
	mObjDim(make_int3(rhs.getObjDimX(), rhs.getObjDimY(), rhs.getObjDimZ())),
	mObjCtrCoord(make_float3(rhs.getObjCtrCoordX(), rhs.getObjCtrCoordY(), rhs.getObjCtrCoordZ())),
	mObjCtrIdx(make_float3(rhs.getObjCtrIdxX(), rhs.getObjCtrIdxY(), rhs.getObjCtrIdxZ())),
	mDetCell(make_float2(rhs.getDetCellWidth(),rhs.getDetCellHeight())),
	mDetNum(make_int2(rhs.getDetNumWidth(), rhs.getDetNumHeight())),
	mDetCtrIdx(make_float2(rhs.getDetCtrIdxWidth(), rhs.getDetCtrIdxHeight())),
	mXds(rhs.getXds()),
	mYds(rhs.getYds()),
	mZds(rhs.getZds()),
	mSour(rhs.getSour())
{

}

void Geometry::correctObjCtrIdx()
{
	mObjCtrIdx.x = (static_cast<float>(mObjDim.x) - 1.0) * 0.5 - mObjCtrCoord.x / mVoxelSize.x;
	mObjCtrIdx.y = (static_cast<float>(mObjDim.y) - 1.0) * 0.5 - mObjCtrCoord.y / mVoxelSize.y;
	mObjCtrIdx.z = (static_cast<float>(mObjDim.z) - 1.0) * 0.5 - mObjCtrCoord.z / mVoxelSize.z;
}



void Geometry::setSourToObj(const float s)
{
	mSourToObj = s;
}
void Geometry::setSourToDet(const float s)
{
	mSourToDet = s;
}
void Geometry::setStartView(const float s)
{
	mStartView = s;
}
void Geometry::setViewPerRot(const int s)
{
	mViewPerRot = s;
	mViewStep = M_PI * 2.0 / static_cast<float>(mViewPerRot);
}
void Geometry::setViewNumber(const int s)
{
	mViewNumber = s;
}
void Geometry::setPitch(const float s)
{
	mPitch = s;
	mZStep = mPitch / static_cast<float>(mViewPerRot);
}
void Geometry::setStartZ(const float s)
{
	mStartZ = s;
}
void Geometry::setVoxelSize(const float3 s)
{
	mVoxelSize = s;
	correctObjCtrIdx();
}
void Geometry::setObjDim(const int3 s)
{
	mObjDim = s;
	correctObjCtrIdx();
}
void Geometry::setObjCtrCoord(const float3 s)
{
	mObjCtrCoord = s;
	correctObjCtrIdx();
}

void Geometry::setDetCell(const float2 s)
{
	mDetCell = s;
}
void Geometry::setDetNum(const int2 s)
{
	mDetNum = s;
}
void Geometry::setDetCtrIdx(const float2 s)
{
	mDetCtrIdx = s;
}

void Geometry::setDetShape(int)
{

}

void Geometry::generateDetCoords()
{
	mXds.clear();
	mYds.clear();
	mZds.clear();
	mXds.resize(mDetNum.x);
	mYds.resize(mDetNum.x);
	mZds.resize(mDetNum.y);

	if (mDetShape == ARC)
	{
		const float stepTheta = atan((mDetCell.x * 0.5) / mSourToDet) * 2.0;
		for (int i = 0; i != mDetNum.x; ++i)
		{
			mXds[i] = sin((static_cast<float>(i) - mDetCtrIdx.x) * stepTheta) * mSourToDet;
			mYds[i] = mSourToObj - cos((static_cast<float>(i) - mDetCtrIdx.x) * stepTheta) * mSourToDet;
		}
		for (int i = 0; i != mDetNum.y; ++i)
		{
			mZds[i] = (static_cast<float>(i) - mDetCtrIdx.y) * mDetCell.y;
		}
	}
	else
	{

	}
}

void Geometry::generateSourCoords()
{
	mSour.clear();
	mSour.resize(mViewNumber);

	for (int i = 0; i != mViewNumber; ++i)
	{
		mSour[i].x = -sin(mStartView + static_cast<float>(i) * mViewStep) * mSourToObj;
		mSour[i].y = cos(mStartView + static_cast<float>(i) * mViewStep) * mSourToObj;
		mSour[i].z = mStartZ + static_cast<float>(i) * mZStep;
	}
}

std::vector<float> Geometry::getAllAngles() const
{
	std::vector<float> toReturn(mViewNumber, 0);
	size_t i = 0;
	std::generate(toReturn.begin(), toReturn.end(), [&]() {
		const float v = mStartView + mViewStep * static_cast<float>(i);
		++i;
		return v; });
	return toReturn;
}
std::vector<float> Geometry::getAllZPoses() const
{
	std::vector<float> toReturn(mViewNumber, 0);
	size_t i = 0;
	std::generate(toReturn.begin(), toReturn.end(), [&]() {
		const float z = mStartZ + mZStep * static_cast<float>(i);
		++i;
		return z;
	});
	return toReturn;
}
