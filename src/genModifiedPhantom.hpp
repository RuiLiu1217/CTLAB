/*!
 * \brief GPU based generating the modified Shepp-Logan phantom routine and save it to the disk with a given name.
 * \author Rui Liu
 * \date Sep. 25, 2014
 * \version 1.0
 * \email liurui1217@gmail.com
*/
#include <string>
/// \brief Generate the modified Shepp-Logan phantom with given size and saved to the disk with a given name single floating type
/// \param FileName Name of the file to be stored for 3d Shepp-Logan phantom
/// \param lenReso length of the phantom resolution
/// \param widReso width of the phantom resolution
/// \param heiReso height of the phantom resolution
void generateModiPhantomFloat(
			 const std::string& FileName, // File Name
			 const unsigned int& lenReso = 256, // Length of the phantom
			 const unsigned int& widReso = 256, // Width of the phantom
			 const unsigned int& heiReso = 256);// Height of the phantom
/// \brief Generate the modified Shepp-Logan phantom with given size and saved to the disk with a given name double floating type
/// \param FileName Name of the file to be stored for 3d Shepp-Logan phantom
/// \param lenReso length of the phantom resolution
/// \param widReso width of the phantom resolution
/// \param heiReso height of the phantom resolution
void generateModiPhantomDouble(
			 const std::string& FileName,
			 const unsigned int& lenReso = 256,
			 const unsigned int& widReso = 256,
			 const unsigned int& heiReso = 256);
