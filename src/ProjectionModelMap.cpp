#include "ProjectionModelMap.h"

ProjectionModelMap::ProjectionModelMap():m(createMap())
{

}
std::map<const std::string, int> ProjectionModelMap::createMap()
{
	std::map<const std::string, int> m;
	m["Branchless DD"] = 0;
	m["Volume Rendering"] = 1;
	m["Double Precision DD"] = 2;
	m["Pseudo DD"] = 3;
	m["Siddon's Algorithm"] = 4;
	m["Distance Driven"] = 5;
	return m;
}
std::map<const std::string, int>& ProjectionModelMap::get() 
{
	return m;
}