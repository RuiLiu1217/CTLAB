#ifndef IMAGE_H_
#define IMAGE_H_
#include <cuda_runtime.h>
/// \brief Image configuration class
class Image {
public:
	uint2 m_Reso; ///< Image resolution
	float2 m_Size;///< Image size
	float2 m_Step; ///< Image Step
	float2 m_Bias; ///< The bias of the image
public:
	/// \brief constructor
	Image(void);
	/// \brief destructor
	~Image(void) {};
	/// \brief copy constructor
	Image(const Image& rhs);
	/// \brief constructor
	Image(
		const unsigned int resoL,///< resolution on length direction
		const unsigned int resoW,///< resolution on width direction
		const float sizeL, ///< length size of the image
		const float sizeW,///< width size of the image
		const float BiasL, ///< bias on length direction
		const float BiasW ///<bias on width direction
	);
};
#endif