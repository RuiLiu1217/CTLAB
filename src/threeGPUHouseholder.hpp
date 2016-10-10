/*! 
 * \file threeGPUHouseholder.hpp
 * \brief Use Three GPUs to implement the Householder transform
 * \author Rui Liu
 * \date Sep.1, 2015
 * \version 1.0
 * \email liurui1217@gmail.com
 */
#include <thrust/device_vector.h>
namespace SVD{
  /// \brief Householder transform 
  /// \param FileName Read the file into the memory
  /// \param dataType define float or double to be used
  /// \param m number of rows of the matrix
  /// \param n number of columns of the matrix
  void threeGPUHouseHolder(const std::string& FileName,const std::string& dataType, const int m, const int n);
  
//  void HouseHolder_v2(thrust::device_vector<float>& B, const int m, const int n);
  void HouseHolder_v2(thrust::device_vector<float>& B, const int m, const int n);
  void HouseHolder_v2(thrust::device_vector<double>& B, const int m, const int n);
  /// \brief test the correctness of the Householder transform
  int testHouseHolder();
}
