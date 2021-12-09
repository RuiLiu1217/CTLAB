#include "Image.h"
#include "readwrite.hpp"
//
//
//
//
//Image::Image(void)
//{
//	m_Bias.x = 0.0f;
//	m_Bias.y = 0.0f;  //���ƫ�Ƶĵ�λ����ʵ���?λ;
//	m_Reso.x = m_Reso.y = 512;
//	m_Size.x = m_Size.y = 400.0;// static_cast<float>(4.484740011196460e+02);
//	m_Step.x = m_Size.x / m_Reso.x;
//	m_Step.y = m_Size.y / m_Reso.y;
//}
//
//Image::Image(const Image& rhs)
//{
//	m_Bias = rhs.m_Bias;
//	m_Reso = rhs.m_Reso;
//	m_Size = rhs.m_Size;
//	m_Step = rhs.m_Step;
//}
//
//
//
//Image::Image(
//	const unsigned int resoL,
//	const unsigned int resoW,
//	const float sizeL,
//	const float sizeW,
//	const float BiasL,
//	const float BiasW) :m_Reso(make_uint2(resoL, resoW)),
//	m_Size(make_float2(sizeL, sizeW)),
//	m_Step(make_float2(sizeL / resoL, sizeW / resoW)),
//	m_Bias(make_float2(BiasL, BiasW)){}


Image::Image() {
	resoX = 512;
	resoY = 512;
	resoZ = 512;

	centX = 0.0f;
	centY = 0.0f;
	centZ = 0.0f;

	dx = 200.0f / resoX;
	dy = 200.0f / resoY;
	dz = dx;

	data.resize(resoX * resoY * resoZ);
	for (int i = 0; i < 512 * 512 * 512; ++i) {
		data[i] = 1.0f;
	}
}

Image::Image(std::string fileName) {}

int Image::getResoX() const {
	return resoX;
}
int Image::getResoY() const {
	return resoY;
}
int Image::getResoZ() const {
	return resoZ;
}

float Image::getDx() const {
	return dx;
}
float Image::getDy() const {
	return dy;
}
float Image::getDz() const {
	return dz;
}

float Image::getCentX() const {
	return centX;
}
float Image::getCentY() const {
	return centY;
}
float Image::getCentZ() const {
	return centZ;
}

float* Image::getDataPtr() {
	return &(data[0]);
}

void Image::saveRaw(std::string fileName) {
	writeData<float>(data, fileName);
}