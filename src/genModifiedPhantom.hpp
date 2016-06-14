/*
 * GPU based generating the modified Shepp-Logan phantom 
 * routine and save it to the disk with a given name.
 * Author: Rui Liu
 * Date: 09/24/14
 * Version: 1.0
 * Email: liurui1217@gmail.com
*/

#include <string>
// Generate the modified Shepp-Logan phantom with given size and saved to the disk with a given name
// Single Floating type
void generateModiPhantomFloat(
			 const std::string& FileName, // File Name
			 const unsigned int& lenReso = 256, // Length of the phantom
			 const unsigned int& widReso = 256, // Width of the phantom
			 const unsigned int& heiReso = 256);// Height of the phantom
// Double floating type
void generateModiPhantomDouble(
			 const std::string& FileName,
			 const unsigned int& lenReso = 256,
			 const unsigned int& widReso = 256,
			 const unsigned int& heiReso = 256);
