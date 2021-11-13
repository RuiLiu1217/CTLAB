#include "BackprojectionModelMap.h"


std::map<std::string, int> BackprojectionModelMap::createMap()
{
	std::map<std::string, int> m;
	m["Branchless DD"] = 0;
	m["Pseudo DD"] = 1;
	m["Double Precision DD"] = 2;
	m["Z Line"] = 3;
	m["Distance Driven"] = 4;
	return m;
}
BackprojectionModelMap::BackprojectionModelMap() :m(createMap()) {}

std::map<std::string, int>& BackprojectionModelMap::get()
{
	return m;
}