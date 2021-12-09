/*!
* \file readwrite.hpp
* \brief This file mainly includes read and write functions called by users. The author strongly suggest the users do NOT modify any codes but you can add your own functions and classes.
* \author Rui Liu
* \version 1.0
* \date 08/24/14
*/
#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
template<typename T>
std::vector<T> readData(const std::string& fileName, const size_t fileSize) {
	std::vector<T> res(fileSize, 0);
	std::ifstream fin(fileName.c_str(), std::ios::binary);
	if (!fin.is_open()) {
		std::cout << "Cannot open File " << fileName << std::endl;
		exit(-1);
	}
	fin.read((char*)res, sizeof(T) * fileSize);
	fin.close();
	return res;
}
template<typename T>
bool writeData(const std::vector<T>& source, const std::string& fileName) {
	std::ofstream fout(fileName.c_str(), std::ios::binary);
	if (!fout.is_open()) {
		std::cout << "Cannot open file " << fileName << " to write\n";
		return -1;
	}
	const size_t fileSize = source.size();
	fout.write((char*)(&source[0]), sizeof(T) * fileSize);
	fout.close();
	return 0;
}

template<typename T>
bool writeData(T* source, const size_t fileSize, const std::string& fileName) {
	std::ofstream fout(fileName.c_str(), std::ios::binary);
	if (!fout.is_open()) {
		std::cout << "Cannot open file " << fileName << " to write\n";
		return -1;
	}
	
	fout.write((char*)(source), sizeof(T) * fileSize);
	fout.close();
	return 0;
}