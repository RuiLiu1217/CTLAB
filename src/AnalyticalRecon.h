#pragma once
class AnalyticalRecon
{
private:
	bool useGPU;
	Geometry geo; // Reconstruction geometry
public:
	AnalyticalRecon(const Geometry& geo, bool useGPU);
	virtual ~AnalyticalRecon();

	const bool getUseGPU() const; 
	const Geometry& getGeometry() const;

	bool getUseGPU();
	Geometry& getGeometry();
};

