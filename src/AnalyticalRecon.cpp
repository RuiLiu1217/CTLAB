#include "AnalyticalRecon.h"
#include "Geometry.h"

AnalyticalRecon::AnalyticalRecon(const Geometry& g, bool ug):geo(g),useGPU(ug)
{
}


AnalyticalRecon::~AnalyticalRecon()
{
}

const bool AnalyticalRecon::getUseGPU() const
{
	return useGPU;
}
bool AnalyticalRecon::getUseGPU()
{
	return useGPU;
}
const Geometry& AnalyticalRecon::getGeometry() const { return geo; }
Geometry& AnalyticalRecon::getGeometry() { return geo; }