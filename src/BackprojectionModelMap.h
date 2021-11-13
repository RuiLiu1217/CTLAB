#pragma once
#ifndef BACKPROJECTIONMODELMAP
#define BACKPROJECTIONMODELMAP
#include <map>
#include <string>
class BackprojectionModelMap
{
public:
	BackprojectionModelMap();// : m(createMap()) {}
	std::map<std::string, int>& get();
private:
	std::map<std::string, int> m;
	std::map<std::string, int> createMap();
};

#endif