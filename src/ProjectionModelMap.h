#pragma once
#ifndef PROJECTIONMODELMAP
#define PROJECTIONMODELMAP
#include <map>
#include <string>
class ProjectionModelMap
{
public:
	ProjectionModelMap();
	std::map<const std::string, int>& get();
private:
	std::map<const std::string, int> m;
	std::map<const std::string, int> createMap();
};
#endif