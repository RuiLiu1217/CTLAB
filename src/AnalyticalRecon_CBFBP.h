#pragma once
#include "AnalyticalRecon.h"
class Geometry;
// This method is based on the paper 
// "A three-dimentional-weighted cone beam filtered backprojection (CB-FBP) algorithm for
// image reconstruction in volumetric CT-helical scanning", Phys Med Biol
class AnalyticalRecon_CBFBP :
	public AnalyticalRecon
{
public:
	AnalyticalRecon_CBFBP(const Geometry&, bool useGPU);
	virtual ~AnalyticalRecon_CBFBP();
	virtual void rebinning(
		std::vector<float>&,
		const std::vector<float>&, float deltaT); //Rebinning the projection data
	virtual void filtering(
		std::vector<float>&,
		const std::vector<float>&);
	virtual void backprojection(
		std::vector<float>&,
		const std::vector<float>&);
private:
	
	void rebinning_CPU(std::vector<float>&, const std::vector<float>&, float deltaT);
	void rebinning_GPU(std::vector<float>&, const std::vector<float>&, float deltaT);
	void filtering_CPU(std::vector<float>&, const std::vector<float>&); // we do not provide this function
	void filtering_GPU(std::vector<float>&, const std::vector<float>&);
	void backprojection_CPU(std::vector<float>& vol, std::vector<float>& proj);
	void backprojection_GPU(std::vector<float>& vol, std::vector<float>& proj);
};

