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
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
typedef const unsigned int cuint;

/// \brief Read data in host memory
/// \param FileName File name to be stored
/// \param inputData the data to be stored
/// \param FileSize elements of number
void readData_host(const std::string& FileName, thrust::host_vector<float>& inputData, cuint FileSize);
/// \brief Read data in host memory
/// \param FileName File name to be stored
/// \param inputData the data to be stored
/// \param FileSize elements of number
void readData_host(const std::string& FileName, thrust::host_vector<double>& inputData, cuint FileSize);

/// \brief Read data in device memory
/// \param FileName File name to be stored
/// \param inputData the data to be stored
/// \param FileSize elements of number
void readData_host(const std::string& FileName, thrust::device_vector<float>& inputData, cuint FileSize);
/// \brief Read data in device memory
/// \param FileName File name to be stored
/// \param inputData the data to be stored
/// \param FileSize elements of number
void readData_host(const std::string& FileName, thrust::device_vector<double>& inputData, cuint FileSize);


/// \brief Read data in host memory
/// \param FileName File name to be stored
/// \param inputData the data to be stored
/// \param FileSize elements of number
void readData_host(const std::string& FileName, std::vector<float>& inputData, cuint FileSize);
/// \brief Read data in host memory
/// \param FileName File name to be stored
/// \param inputData the data to be stored
/// \param FileSize elements of number
void readData_host(const std::string& FileName, std::vector<double>& inputData, cuint FileSize);
/// \brief Read data in host memory
/// \param FileName File name to be stored
/// \param inputData the data to be stored
/// \param FileSize elements of number
void readData_host(const std::string& FileName, float* inputData, cuint FileSize);
/// \brief Read data in host memory
/// \param FileName File name to be stored
/// \param inputData the data to be stored
/// \param FileSize elements of number
void readData_host(const std::string& FileName, double* inputData, cuint FileSize);

/// \brief Write data from host memory
/// \param FileName File name to be stored
/// \param inputData the data to be stored
/// \param FileSize elements of number
void writeData_host(const std::string& FileName,const thrust::host_vector<float>& inputData, cuint FileSize);

/// \brief Write data from host memory
/// \param FileName File name to be stored
/// \param inputData the data to be stored
/// \param FileSize elements of number
void writeData_host(const std::string& FileName,const thrust::host_vector<double>& inputData, cuint FileSize);

/// \brief Write data from device memory
/// \param FileName File name to be stored
/// \param inputData the data to be stored
/// \param FileSize elements of number
void writeData_host(const std::string& FileName,const thrust::device_vector<float>& inputData, cuint FileSize);

/// \brief Write data from device memory
/// \param FileName File name to be stored
/// \param inputData the data to be stored
/// \param FileSize elements of number
void writeData_host(const std::string& FileName,const thrust::device_vector<double>& inputData, cuint FileSize);


/// \brief Write data from host memory
/// \param FileName File name to be stored
/// \param inputData the data to be stored
/// \param FileSize elements of number
void writeData_host(const std::string& FileName,const std::vector<float>& inputData, cuint FileSize);
/// \brief Write data from host memory
/// \param FileName File name to be stored
/// \param inputData the data to be stored
/// \param FileSize elements of number
void writeData_host(const std::string& FileName,const std::vector<double>& inputData, cuint FileSize);
/// \brief Write data from host memory
/// \param FileName File name to be stored
/// \param inputData the data to be stored
/// \param FileSize elements of number
void writeData_host(const std::string& FileName, float* inputData, cuint FileSize);
/// \brief Write data from host memory
/// \param FileName File name to be stored
/// \param inputData the data to be stored
/// \param FileSize elements of number
void writeData_host(const std::string& FileName, double* inputData, cuint FileSize);
/// \brief Write data from device memory
/// \param FileName File name to be stored
/// \param data the data to be stored
/// \param bytes elements of number
void writeFileFromDeviceToDisk(const std::string& FileName, double* data, cuint bytes);
/// \brief Write data from device memory
/// \param FileName File name to be stored
/// \param data the data to be stored
/// \param bytes elements of number
void writeFileFromDeviceToDisk(const std::string& FileName, float* data, cuint bytes);
