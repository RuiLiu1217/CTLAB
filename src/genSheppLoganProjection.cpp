#include "genSheppLoganProjection.h"



std::vector<float> calculateProjection(CTLAB::SheppLoganPhantom& phantom,
	const float sid, const float sdd, const int DNU, std::vector<float> angs, const float col_size, const float col_offset)
{
	const int PN = angs.size();
	std::vector<float> proj(DNU * PN, 0);
	CTLAB::float2 sour(0.0f, sid);

	std::vector<CTLAB::float2> ends(DNU);
	float stepTheta = atanf((col_size * 0.5) / sdd) * 2.0;
	float curBeta = 0;
	for (int ii = 0; ii != DNU; ++ii)
	{
		curBeta = (ii - (DNU - 1.0) * 0.5 + col_offset) * stepTheta;
		ends[ii].x = sinf(curBeta) * sdd;
		ends[ii].y = sid - cosf(curBeta) * sdd;
	}


	CTLAB::float2 startPoint;
	CTLAB::float2 endPoint;
	for (int angIdx = 0; angIdx != PN; ++angIdx)
	{
		startPoint = sour.rot(angs[angIdx]);
		for (int detIdx = 0; detIdx != DNU; ++detIdx)
		{
			endPoint = ends[detIdx].rot(angs[angIdx]);
			for (int elpIdx = 0; elpIdx != 10; ++elpIdx)
			{
				proj[angIdx * DNU + detIdx] += phantom.ele[elpIdx].intersectionLength(startPoint, endPoint) * phantom.ele[elpIdx].getIntensity();
			}
		}
	}

	return proj;
}
