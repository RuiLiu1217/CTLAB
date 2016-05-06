#include <thrust/device_vector.h>
namespace SVD{
  //FileName : the matrix in raw form that you want to transform
  //datatype : use "float" or "double"
  void threeGPUHouseHolder(const std::string& FileName,const std::string& dataType, const int m, const int n);
//  void HouseHolder_v2(thrust::device_vector<float>& B, const int m, const int n);
  void HouseHolder_v2(thrust::device_vector<float>& B, const int m, const int n);
  void HouseHolder_v2(thrust::device_vector<double>& B, const int m, const int n);
  int testHouseHolder();
}
