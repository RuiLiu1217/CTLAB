#ifndef VOLUME_H_
#define VOLUME_H_
#include <cuda_runtime.h>
/// \brief Volume configuration class
class Volume {
public:
	uint3 m_Reso; ///< Volume resolution
	float3 m_Size; ///< Image size
	float3 m_Step; ///< Image step size
	float3 m_Bias; ///< Bias of the image
public:
	/// \brief constructor
	Volume(void);
	/// \brief destructor
	~Volume(void) {}
	/// \brief copy constructor
	Volume(const Volume& rhs);
	/// \brief
	Volume(
		const unsigned int resoL, ///< object length resolution
		const unsigned int resoW,///<object width resolution
		const unsigned int resoH, ///< object height resolution
		const float sizeL, ///< object size on length direction
		const float sizeW,///< object size on width direction
		const float sizeH,///< object size on height direction
		const float biasL,///< object bias on length direction
		const float biasW,///< object bias on width direction
		const float biasH///< object bias on height direction
	);
};

#endif