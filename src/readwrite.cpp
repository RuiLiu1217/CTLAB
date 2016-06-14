#include "readwrite.hpp"
typedef unsigned int uint;
#include <fstream>
template<typename T>
void readData_host_template(const std::string& FileName, thrust::host_vector<T>& inputData, cuint FileSize)
{
	inputData.clear();
	inputData.resize(FileSize);
	std::ifstream fin(FileName.c_str(),std::ios::binary);
	if (!fin.is_open())
	{
		std::cout<<"Cannot open File "<<FileName<<std::endl;
		exit(-1);
	}
	fin.read((char*)(&inputData[0]),sizeof(T) * FileSize);
	fin.close();

}

template<typename T>
void readData_host_template(const std::string& FileName, std::vector<T>& inputData, cuint FileSize)
{
	inputData.clear();
	inputData.resize(FileSize);
	std::ifstream fin(FileName.c_str(),std::ios::binary);
	if (!fin.is_open())
	{
		std::cout<<"Cannot open File "<<FileName<<std::endl;
		exit(-1);
	}
	fin.read((char*)(&inputData[0]),sizeof(T) * FileSize);
	fin.close();

}

template<typename T>
void readData_host_template(const std::string& FileName, T* inputData, cuint FileSize)
{
	if (inputData != NULL)
	{
		delete[] inputData;
		inputData = NULL;
	}
	inputData = new T[FileSize];
	std::ifstream fin(FileName.c_str(),std::ios::binary);
	if (!fin.is_open())
	{
		std::cout<<"Cannot open File "<<FileName<<std::endl;
		exit(-1);
	}
	fin.read((char*)inputData,sizeof(T) * FileSize);
	fin.close();

}

void readData_host(const std::string& FileName, thrust::host_vector<float>& inputData, cuint FileSize)
{
	readData_host_template<float>(FileName,inputData,FileSize);
}
void readData_host(const std::string& FileName, thrust::host_vector<double>& inputData, cuint FileSize)
{
	readData_host_template<double>(FileName,inputData,FileSize);
}

void readData_host(const std::string& FileName, std::vector<float>& inputData, cuint FileSize)
{
	readData_host_template<float>(FileName,inputData,FileSize);
}
void readData_host(const std::string& FileName, std::vector<double>& inputData, cuint FileSize)
{
	readData_host_template<double>(FileName,inputData,FileSize);
}

void readData_host(const std::string& FileName, float* inputData, cuint FileSize)
{
	readData_host_template<float>(FileName,inputData,FileSize);
}
void readData_host(const std::string& FileName, double* inputData, cuint FileSize)
{
	readData_host_template<double>(FileName,inputData,FileSize);
}



template<typename T>
void writeData_host_template(const std::string& FileName,const thrust::host_vector<T>& outputData,cuint FileSize)
{
	std::ofstream fout(FileName.c_str(),std::ios::binary);
	if (!fout.is_open())
	{
		std::cout<<"Cannot open file "<<FileName<<" to write\n";
		exit(-1);
	}
	fout.write((char*)(&outputData[0]),sizeof(T) * FileSize);
	fout.close();
	std::cout<<"Save File "<<FileName<<" finished\n";
}


template<typename T>
void writeData_device_template(const std::string& FileName,const thrust::device_vector<T>& outData,cuint FileSize)
{
	thrust::host_vector<T> outputData = outData;
	std::ofstream fout(FileName.c_str(),std::ios::binary);
	if (!fout.is_open())
	{
		std::cout<<"Cannot open file "<<FileName<<" to write\n";
		exit(-1);
	}
	fout.write((char*)(&outputData[0]),sizeof(T) * FileSize);
	fout.close();
	std::cout<<"Save File "<<FileName<<" finished\n";
}


template<typename T>
void writeData_host_template(const std::string& FileName,const std::vector<T>& outputData,cuint FileSize)
{
	std::ofstream fout(FileName.c_str(),std::ios::binary);
	if (!fout.is_open())
	{
		std::cout<<"Cannot open file "<<FileName<<" to write\n";
		exit(-1);
	}
	fout.write((char*)(&outputData[0]),sizeof(T) * FileSize);
	fout.close();
	std::cout<<"Save File "<<FileName<<" finished\n";
}

template<typename T>
void writeData_host_template(const std::string& FileName, T* outputData,cuint FileSize)
{
	std::ofstream fout(FileName.c_str(),std::ios::binary);
	if (!fout.is_open())
	{
		std::cout<<"Cannot open file "<<FileName<<" to write\n";
		exit(-1);
	}
	fout.write((char*)outputData,sizeof(T) * FileSize);
	fout.close();
	std::cout<<"Save File "<<FileName<<" finished\n";
}





void writeData_host(const std::string& FileName,const thrust::host_vector<float>& inputData, cuint FileSize)
{
	writeData_host_template<float>(FileName, inputData, FileSize);
}
void writeData_host(const std::string& FileName,const thrust::host_vector<double>& inputData, cuint FileSize)
{
	writeData_host_template<double>(FileName, inputData, FileSize);
}
void writeData_device(const std::string& FileName,const thrust::device_vector<float>& inputData, cuint FileSize)
{
	writeData_device_template<float>(FileName, inputData, FileSize);
}
void writeData_device(const std::string& FileName,const thrust::device_vector<double>& inputData, cuint FileSize)
{
	writeData_device_template<double>(FileName, inputData, FileSize);
}

void writeData_host(const std::string& FileName,const std::vector<float>& inputData, cuint FileSize)
{
	writeData_host_template<float>(FileName, inputData, FileSize);
}
void writeData_host(const std::string& FileName,const std::vector<double>& inputData, cuint FileSize)
{
	writeData_host_template<double>(FileName, inputData, FileSize);
}
void writeData_host(const std::string& FileName, float* inputData, cuint FileSize)
{
	writeData_host_template<float>(FileName, inputData, FileSize);
}
void writeData_host(const std::string& FileName, double* inputData, cuint FileSize)
{
	writeData_host_template<double>(FileName, inputData, FileSize);
}




template<typename T>
void writeFileFromDeviceToDisk_template(const std::string& FileName, T* data, cuint bytes)
{
	T* hdata = new T[bytes];
	cudaMemcpy(hdata,data,sizeof(T) * bytes,cudaMemcpyDeviceToHost);
	std::ofstream fid(FileName.c_str(),std::ios::binary);
	fid.write((char*)hdata,sizeof(T)* bytes);
	fid.close();
	//delete[] hdata;
	//hdata = nullptr;
}
void writeFileFromDeviceToDisk(const std::string& FileName, double* data, cuint bytes)
{
	writeFileFromDeviceToDisk_template<double>(FileName, data, bytes);
}


void writeFileFromDeviceToDisk(const std::string& FileName, float* data, cuint bytes)
{
	writeFileFromDeviceToDisk_template<float>(FileName, data, bytes);
}
