#pragma once
#include <vector>
#include <string>

class Geometry;

class CT
{
public:
	static void Proj(std::vector<float>& hvol, std::vector<float>& hprj, Geometry geo, const std::string& projModel);
	static void Back(std::vector<float>& hvol, std::vector<float>& hprj, Geometry geo, const std::string& backModel);
};