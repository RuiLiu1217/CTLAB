#include "Image.h"




Image::Image(void)
{
	m_Bias.x = 0.0f;
	m_Bias.y = 0.0f;  //���ƫ�Ƶĵ�λ����ʵ���?λ;
	m_Reso.x = m_Reso.y = 512;
	m_Size.x = m_Size.y = 400.0;// static_cast<float>(4.484740011196460e+02);
	m_Step.x = m_Size.x / m_Reso.x;
	m_Step.y = m_Size.y / m_Reso.y;
}

Image::Image(const Image& rhs)
{
	m_Bias = rhs.m_Bias;
	m_Reso = rhs.m_Reso;
	m_Size = rhs.m_Size;
	m_Step = rhs.m_Step;
}



Image::Image(
	const unsigned int resoL,
	const unsigned int resoW,
	const float sizeL,
	const float sizeW,
	const float BiasL,
	const float BiasW) :m_Reso(make_uint2(resoL, resoW)),
	m_Size(make_float2(sizeL, sizeW)),
	m_Step(make_float2(sizeL / resoL, sizeW / resoW)),
	m_Bias(make_float2(BiasL, BiasW)){}
