#include "AnalyticalRecon.h"
#include "AnalyticalRecon_CBFBP.h"
#include "Geometry.h"
#include <vector>
// Parallel Rebinning CB curve detector;
static void ParallelRebinningCBCurve_GPU(std::vector<float>& outputProj,
	std::vector<float>& Proj, const Geometry& geo, const int threadX,
	const int threadY, const int threadZ)
{

}

void AnalyticalRecon_CBFBP::rebinning_GPU()
{
	if (getUseGPU())
	{
		ParallelRebinningCBCurve_GPU(std::vector<float>& outputProj,
			std::vector<float>& Proj, const Geometry& geo, const int threadX,
			const int threadY, const int threadZ)
	}
	else
	{

	}
}