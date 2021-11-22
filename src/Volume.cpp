#include "Volume.h"


Volume::Volume(void)
{
	m_Reso.x = m_Reso.y = 512;
	m_Reso.z = 256;
	m_Size.x = m_Size.y = 20.0f;
	m_Size.z = 10.0f;
	m_Step.x = m_Size.x / m_Reso.x;
	m_Step.y = m_Size.y / m_Reso.y;
	m_Step.z = m_Size.z / m_Reso.z;
	m_Bias = make_float3(0.0f, 0.0f, 0.0f);
}

Volume::Volume(const Volume& rhs)
{
	m_Bias = rhs.m_Bias;
	m_Reso = rhs.m_Reso;
	m_Size = rhs.m_Size;
	m_Step = rhs.m_Step;
}

Volume::Volume(
	const unsigned int resoL,
	const unsigned int resoW,
	const unsigned int resoH,
	const float sizeL,
	const float sizeW,
	const float sizeH,
	const float biasL,
	const float biasW,
	const float biasH) :m_Reso(make_uint3(resoL, resoW, resoH)),
	m_Size(make_float3(sizeL, sizeW, sizeH)),
	m_Step(make_float3(sizeL / resoL, sizeW / resoW, sizeH / resoH)),
	m_Bias(make_float3(biasL, biasW, biasH)){}