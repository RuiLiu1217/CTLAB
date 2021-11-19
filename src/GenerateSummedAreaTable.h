#ifndef GENERATE_SUMMED_AREA_TABLE_H
#define GENERATE_SUMMED_AREA_TABLE_H
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/pair.h>
void genSAT_fof_Volume(float* hvol,	thrust::device_vector<float>& ZXY, thrust::device_vector<float>& ZYX, int XN, int YN, int ZN);
void genSAT_fof_Volume_alreadyinGPU(const thrust::device_vector<float>& vol, thrust::device_vector<float>& ZXY, thrust::device_vector<float>& ZYX, int XN, int YN, int ZN);
thrust::pair<thrust::device_vector<float>, thrust::device_vector<float>> generateSumAreaTableForVolume(float* hvol,
	int XN, int YN, int ZN);
void genSAT_fof_Volume_alreadyinGPU(const thrust::device_vector<float>& vol, thrust::device_vector<float>& ZXY, thrust::device_vector<float>& ZYX, int XN, int YN, int ZN);
#endif