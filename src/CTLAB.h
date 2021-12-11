#pragma once
#include <vector>
#include <string>

/*!
 * \brief Total header files
 * \author Rui Liu
 * \date Nov. 22, 2021
 * \version 1.1
 * \email liurui1217@gmail.com
*/

extern "C" { // Tested
	

	/// \brief Generate the modified Shepp-Logan phantom with given size and saved to the disk with a given name
	/// \param FileName Name of the file to be stored for 3d Shepp-Logan phantom
	/// \param lenReso length of the phantom resolution
	/// \param widReso width of the phantom resolution
	/// \param heiReso height of the phantom resolution
	void generateModiPhantomFloat(
		const std::string& FileName, // File Name
		const unsigned int& lenReso = 256, // Length of the phantom
		const unsigned int& widReso = 256, // Width of the phantom
		const unsigned int& heiReso = 256);// Height of the phantom
	void generateModiPhantomDouble(
		const std::string& FileName,
		const unsigned int& lenReso = 256,
		const unsigned int& widReso = 256,
		const unsigned int& heiReso = 256);
}

/// \brief Generate the analytical projection of a Shepp-Logan Phantom. Assume the detector shape is arc.
/// \param sid source to iso-center distance
/// \param sdd source to detector distance
/// \param DNU number of detector cell along channel direction
/// \param angs all views
/// \param col_size detector size 
/// \param col_offset offset of the detector along column direction
/// Example: calculateProjection(3.0f, 5.0f, 888, angs, 0.003f, 0.0f); // The phantom size is 1.0x1.0
std::vector<float> calculateProjection(
	const float sid, const float sdd, const int DNU, std::vector<float> angs, const float col_size, const float col_offset);

// class Geometry;
//
//class CT
//{
//public:
//	static void Proj(std::vector<float>& hvol, std::vector<float>& hprj, Geometry geo, const std::string& projModel);
//	static void Back(std::vector<float>& hvol, std::vector<float>& hprj, Geometry geo, const std::string& backModel);
//};