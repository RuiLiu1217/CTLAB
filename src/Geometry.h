#pragma once
#include <cuda_runtime.h>
#include <vector>

class Geometry
{
public:
	enum DetShape { ARC, PANEL };
	Geometry();
	virtual ~Geometry();
	Geometry(const Geometry&);
private:
	float mSourToObj;
	float mSourToDet;
	float mStartView;
	int mViewPerRot;
	int mViewNumber;
	float mViewStep;

	float mPitch;
	float mZStep;

	float mStartZ;

	float3 mVoxelSize;
	int3 mObjDim;
	float3 mObjCtrCoord;
	float3 mObjCtrIdx;


	float2 mDetCell;
	int2 mDetNum;
	float2 mDetCtrIdx;

	std::vector<float> mXds;
	std::vector<float> mYds;
	std::vector<float> mZds;

	std::vector<float3> mSour;
	void correctObjCtrIdx();
	DetShape mDetShape;
public:
	float getSourToObj() const { return mSourToObj; }
	float getSourToDet() const { return mSourToDet; }
	float getStartView() const { return mStartView; }
	int getViewPerRot() const { return mViewPerRot; }
	int getViewNumber() const { return mViewNumber; }
	float getViewStep() const { return mViewStep; }
	float getPitch() const { return mPitch; }
	float getStartZ() const { return mStartZ; }
	float getZStep() const { return mZStep; }

	float getVoxelSizeX() const { return mVoxelSize.x; }
	float getVoxelSizeY() const { return mVoxelSize.y; }
	float getVoxelSizeZ() const { return mVoxelSize.z; }

	int getObjDimX() const { return mObjDim.x; }
	int getObjDimY() const { return mObjDim.y; }
	int getObjDimZ() const { return mObjDim.z; }

	float getObjCtrCoordX() const { return mObjCtrCoord.x; }
	float getObjCtrCoordY() const { return mObjCtrCoord.y; }
	float getObjCtrCoordZ() const { return mObjCtrCoord.z; }

	float getObjCtrIdxX() const { return mObjCtrIdx.x; }
	float getObjCtrIdxY() const { return mObjCtrIdx.y; }
	float getObjCtrIdxZ() const { return mObjCtrIdx.z; }

	float getDetCellWidth() const { return mDetCell.x; }
	float getDetCellHeight() const { return mDetCell.y; }

	int getDetNumWidth() const { return mDetNum.x; }
	int getDetNumHeight() const { return mDetNum.y; }

	float getDetCtrIdxWidth() const { return mDetCtrIdx.x; }
	float getDetCtrIdxHeight() const { return mDetCtrIdx.y; }

	std::vector<float>& getXds() { return mXds; }
	std::vector<float>& getYds() { return mYds; }
	std::vector<float>& getZds() { return mZds; }

	const std::vector<float>& getXds() const { return mXds; }
	const std::vector<float>& getYds() const { return mYds; }
	const std::vector<float>& getZds() const { return mZds; }

	std::vector<float3>& getSour() { return mSour; }
	const std::vector<float3>& getSour() const { return mSour; }

	void setSourToObj(const float);
	void setSourToDet(const float);
	void setStartView(const float);
	void setViewPerRot(const int);
	void setViewNumber(const int);
	void setPitch(const float);
	void setStartZ(const float);
	void setVoxelSize(const float3);
	void setObjDim(const int3);
	void setObjCtrCoord(const float3);

	void setDetCell(const float2);
	void setDetNum(const int2);
	void setDetCtrIdx(const float2);

	void setDetShape(int); // not implemented

	void generateDetCoords(); // generate the coordinates of the detector cells
	void generateSourCoords(); // generate the source coordinates

	std::vector<float> getAllAngles() const;
	std::vector<float> getAllZPoses() const;

	using byte = unsigned char;
};
