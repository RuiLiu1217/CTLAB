#ifndef IMAGE_H_
#define IMAGE_H_
#include <string>
#include <vector>
///// \brief Image configuration class
//class Image {
//public:
//	uint2 m_Reso; ///< Image resolution
//	float2 m_Size;///< Image size
//	float2 m_Step; ///< Image Step
//	float2 m_Bias; ///< The bias of the image
//public:
//	/// \brief constructor
//	Image(void);
//	/// \brief destructor
//	~Image(void) {};
//	/// \brief copy constructor
//	Image(const Image& rhs);
//	/// \brief constructor
//	Image(
//		const unsigned int resoL,///< resolution on length direction
//		const unsigned int resoW,///< resolution on width direction
//		const float sizeL, ///< length size of the image
//		const float sizeW,///< width size of the image
//		const float BiasL, ///< bias on length direction
//		const float BiasW ///<bias on width direction
//	);
//};

class Image {
public:
	Image();
	Image(std::string fileName);
	int getResoX() const;
	int getResoY() const;
	int getResoZ() const;

	float getDx() const;
	float getDy() const;
	float getDz() const;

	float getCentX() const;
	float getCentY() const;
	float getCentZ() const;

	float* getDataPtr();

	void save(std::string fileName) {
		// Todo: not implemented yet
	}

	void saveRaw(std::string fileName);
private:
	int resoX;
	int resoY;
	int resoZ;

	float dx;
	float dy;
	float dz;

	float centX;
	float centY;
	float centZ;

	std::vector<float> data;
};

#endif