/*!
* \file projbackproj.hpp
* \brief This file mainly includes projection processes and backprojection processes with different configuration and different projection and backprojection models.
* \author Rui Liu
* \version 1.0
* \date 08/24/14
*/
#include "utilities.hpp"


/// \brief Define the brieviation of the const unsigned int;
typedef const unsigned int cuint;


//////////////////////////////////////////////////////////////////////////
// This function is the projection difference with weighting;
// dimg: image to be projected
// donp: one angle projection data difference with SART weighting
// draw: raw projection data
// angIdx: angle index
// angDir: detector from positive to negative or verse
// FanGeo: Fan Beam geometry
// Img: Image configuration;
// blk: it is better to be 1024, 1
// gid: correspondingly configurated
//////////////////////////////////////////////////////////////////////////
/// \brief Projection difference with weighting
/// \param dimg device image
/// \param donp result the project difference with weighting
/// \param draw original projection
/// \param angIdx the index of the projection angle
/// \param angDir the indexing direction from large to smaller or inversely
/// \param FanGeo Fan beam geometry with equal angle detector
/// \param Img image configuration
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void proj(float* dimg, float* donp, float* draw, cuint angIdx, const bool angDir, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);
/// \brief Projection difference with weighting
/// \param dimg device image
/// \param donp result the project difference with weighting
/// \param draw original projection
/// \param angIdx the index of the projection angle
/// \param angDir the indexing direction from large to smaller or inversely
/// \param FanGeo Fan beam geometry with equal distance detector
/// \param Img image configuration
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void proj(float* dimg, float* donp, float* draw, cuint angIdx, const bool angDir, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid); //Not tested
/// \brief Projection difference with weighting
/// \param dvol device image
/// \param donp result the project difference with weighting
/// \param draw original projection
/// \param angIdx the index of the projection angle
/// \param angDir the indexing direction from large to smaller or inversely
/// \param ConeGeo Cone beam geometry with equal angle detector
/// \param Vol volume configuration
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void proj(float* dvol, float* donp, float* draw, cuint angIdx, const bool angDir, const ConeEAGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid); //Not tested
/// \brief Projection difference with weighting
/// \param dvol device image
/// \param donp result the project difference with weighting
/// \param draw original projection
/// \param angIdx the index of the projection angle
/// \param angDir the indexing direction from large to smaller or inversely
/// \param ConeGeo Cone beam geometry with equal angle detector
/// \param Vol volume configuration
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void proj(float* dvol, float* donp, float* draw, cuint angIdx, const bool angDir, const ConeEDGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid); //Not tested



//////////////////////////////////////////////////////////////////////////
// This function is the projection difference with weighting;
// �����������subset �е�ͶӰ�в�ļ�Ȩ;
// dimg: image to be projected
// donp: projection data difference with weighting inside the subset
// draw: raw projection data
// angDir: detector from positive to negative or verse
// FanGeo: Fan Beam geometry
// Img: Image configuration;
// numPerSubSet: how many angles are included in the subset
// subSetNum: How many subsets are divided from the original projection;
// curSubSetIdx: current subset index
// blk: it is better to be 1024, 1
// gid: correspondingly configurated
//////////////////////////////////////////////////////////////////////////




/// \brief Projection difference with weighting in subsets
/// \param dimg the image for projection
/// \param donp the projection difference with weighting
/// \param draw the raw projection
/// \param angDir the indexing order from large to small or from small to large, usually we use "true"
/// \param FanGeo The fan beam geometry with equal angle detector
/// \param Img the image configuration
/// \param numPerSubSet number of projections in one subset
/// \param subSetNum how many subsets are there
/// \param curSubSetIdx current subset index
/// \param blk block size of CUDA
/// \param gid grid size for CUDA
void proj(float* dimg, float* donp, float* draw, const bool angDir, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);
/// \brief Projection difference with weighting in subsets
/// \param dimg the image for projection
/// \param donp the projection difference with weighting
/// \param draw the raw projection
/// \param angDir the indexing order from large to small or from small to large, usually we use "true"
/// \param FanGeo The fan beam geometry with equal distance detector
/// \param Img the image configuration
/// \param numPerSubSet number of projections in one subset
/// \param subSetNum how many subsets are there
/// \param curSubSetIdx current subset index
/// \param blk block size of CUDA
/// \param gid grid size for CUDA
void proj(float* dimg, float* donp, float* draw, const bool angDir, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);
/// \brief Projection difference with weighting in subsets
/// \param dvol the image for projection
/// \param donp the projection difference with weighting
/// \param draw the raw projection
/// \param angDir the indexing order from large to small or from small to large, usually we use "true"
/// \param ConeGeo The fan beam geometry with equal angle detector
/// \param Vol the image configuration
/// \param numPerSubSet number of projections in one subset
/// \param subSetNum how many subsets are there
/// \param curSubSetIdx current subset index
/// \param blk block size of CUDA
/// \param gid grid size for CUDA
void proj(float* dvol, float* donp, float* draw, const bool angDir, const ConeEAGeo& ConeGeo, const Volume& Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);
/////////////// \brief Projection difference with weighting in subsets
/////////////// \param dvol the image for projection
/////////////// \param donp the projection difference with weighting
/////////////// \param draw the raw projection
/////////////// \param angDir the indexing order from large to small or from small to large, usually we use "true"
/////////////// \param ConeGeo The fan beam geometry with equal distance detector
/////////////// \param Vol the image configuration
/////////////// \param numPerSubSet number of projections in one subset
/////////////// \param subSetNum how many subsets are there
/////////////// \param curSubSetIdx current subset index
/////////////// \param blk block size of CUDA
/////////////// \param gid grid size for CUDA
////////////void proj(float* dvol, float* donp, float* draw, const bool angDir, const ConeEDGeo& ConeGeo, const Volume& Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);




//////////////////////////////////////////////////////////////////////////
// This function is the projection difference with weighting;
// �����������subset �е�ͶӰ�в�ļ�Ȩ; ����������,ֻ��Ϊ������;
// dimg: image to be projected
// donp: projection data difference with weighting inside the subset
// draw: raw projection data
// FanGeo: Fan Beam geometry
// Img: Image configuration;
// numPerSubSet: how many angles are included in the subset
// subSetNum: How many subsets are divided from the original projection;
// curSubSetIdx: current subset index
// blk: it is better to be 1024, 1
// gid: correspondingly configurated
//////////////////////////////////////////////////////////////////////////

/// \brief Projection difference with weighting in subsets
/// \param dimg the image for projection
/// \param donp the projection difference with weighting
/// \param draw the raw projection
/// \param FanGeo The fan beam geometry with equal angle detector
/// \param Img the image configuration
/// \param numPerSubSet number of projections in one subset
/// \param subSetNum how many subsets are there
/// \param curSubSetIdx current subset index
/// \param blk block size of CUDA
/// \param gid grid size for CUDA
void proj(float* dimg, float* donp, float* draw, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);
/// \brief Projection difference with weighting in subsets
/// \param dimg the image for projection
/// \param donp the projection difference with weighting
/// \param draw the raw projection
/// \param FanGeo The fan beam geometry with equal distance detector
/// \param Img the image configuration
/// \param numPerSubSet number of projections in one subset
/// \param subSetNum how many subsets are there
/// \param curSubSetIdx current subset index
/// \param blk block size of CUDA
/// \param gid grid size for CUDA
void proj(float* dimg, float* donp, float* draw, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);

/// \brief Projection difference with weighting in subsets
/// \param dimg the image for projection
/// \param donp the projection difference with weighting
/// \param draw the raw projection
/// \param ConeGeo The fan beam geometry with equal angle detector
/// \param Vol the image configuration
/// \param numPerSubSet number of projections in one subset
/// \param subSetNum how many subsets are there
/// \param curSubSetIdx current subset index
/// \param blk block size of CUDA
/// \param gid grid size for CUDA
void proj(float* dimg, float* donp, float* draw, const ConeEAGeo& ConeGeo, const Volume& Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);
//void proj(float* dimg, float* donp, float* draw, const ConeEDGeo& FanGeo, const Volume& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid); 




//////////////////////////////////////////////////////////////////////////
// This function is the projection difference with weighting;
// ��������������нǶȵ�ͶӰ;
// dimg: image to be projected
// draw: projection data
// angDir: detector from positive to negative or verse
// FanGeo: Fan Beam geometry
// Img: Image configuration;
// blk: it is better to be 1024, 1
// gid: correspondingly configurated
//////////////////////////////////////////////////////////////////////////

/// \brief Projection in all angles
/// \param dimg the image for projection
/// \param draw the raw projection
/// \param angDir indexing direction from small to large or from large to small
/// \param FanGeo The fan beam geometry with equal angle detector
/// \param Img the image configuration
/// \param blk block size of CUDA
/// \param gid grid size for CUDA
void proj(float* dimg, float* draw, const bool angDir, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);



void proj_CPU(float*dimg, float* dprj, const FanEAGeo& FanGeo, const Image& Img);
void proj_CPU_OPENMP(float*dimg, float* dprj, const FanEAGeo& FanGeo, const Image& Img);



/// \brief Projection in all angles
/// \param dimg the image for projection
/// \param draw the raw projection
/// \param angDir indexing direction from small to large or from large to small
/// \param FanGeo The fan beam geometry with equal distance detector
/// \param Img the image configuration
/// \param blk block size of CUDA
/// \param gid grid size for CUDA
void proj(float* dimg, float* draw, const bool angDir, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);

/// \brief Projection in all angles
/// \param dvol the image for projection
/// \param draw the raw projection
/// \param angDir indexing direction from small to large or from large to small
/// \param ConeGeo The Cone beam geometry with equal angle detector
/// \param Vol the image configuration
/// \param blk block size of CUDA
/// \param gid grid size for CUDA
void proj(float* dvol, float* draw, const bool angDir, const ConeEAGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid);

/// \brief Projection in all angles
/// \param dvol the image for projection
/// \param draw the raw projection
/// \param angDir indexing direction from small to large or from large to small
/// \param ConeGeo The Cone beam geometry with equal angle detector
/// \param Vol the image configuration
/// \param blk block size of CUDA
/// \param gid grid size for CUDA
void proj(float* dvol, float* draw, const bool angDir, const ConeEDGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid);



/// \brief Projection in all angles
/// \param dimg the image for projection
/// \param draw the raw projection
/// \param FanGeo The fan beam geometry with equal angle detector
/// \param Img the image configuration
/// \param blk block size of CUDA
/// \param gid grid size for CUDA
void proj(float* dimg, float* draw, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);

/// \brief Projection in all angles
/// \param dimg the image for projection
/// \param draw the raw projection
/// \param FanGeo The fan beam geometry with equal distance detector
/// \param Img the image configuration
/// \param blk block size of CUDA
/// \param gid grid size for CUDA
void proj(float* dimg, float* draw, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);


/// \brief Projection in all angles
/// \param dvol the image for projection
/// \param draw the raw projection
/// \param ConeGeo The Cone beam geometry with equal angle detector
/// \param Vol the image configuration
/// \param blk block size of CUDA
/// \param gid grid size for CUDA
void proj(float* dvol, float* draw, const ConeEAGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid); //OK

/// \brief Projection in all angles
/// \param dvol the image for projection
/// \param draw the raw projection
/// \param ConeGeo The Cone beam geometry with equal angle detector
/// \param Vol the image configuration
/// \param blk block size of CUDA
/// \param gid grid size for CUDA
void proj(float* dvol, float* draw, const ConeEDGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid);


/// \brief Projection in all angles
/// \param dvol the image for projection
/// \param draw the result that projected
/// \param ConeGeo Cone Beam geometry with equal angle detector
/// \param Vol Volume configuration
void proj_CPU(float* dvol, float* draw, const ConeEAGeo ConeGeo, const Volume Vol); //OK

/// \brief Projection in all angles with OpenMP
/// \param dvol the image for projection
/// \param draw the result that projected
/// \param ConeGeo Cone Beam geometry with equal angle detector
/// \param Vol Volume configuration
void proj_CPU_OPENMP(float* dvol, float* draw, const ConeEAGeo ConeGeo, const Volume Vol);

/// \brief Projection differences with weighting in all angles
/// \param dimg the image for projection
/// \param dcor the result 
/// \param draw the raw projection data
/// \param FanGeo The Cone beam geometry with equal angle detector
/// \param Img the image configuration
/// \param blk block size of CUDA
/// \param gid grid size for CUDA
void proj(float* dimg, float* dcor, float* draw, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);

/// \brief Projection differences with weighting in all angles
/// \param dimg the image for projection
/// \param dcor the result 
/// \param draw the raw projection data
/// \param FanGeo The Cone beam geometry with equal distance detector
/// \param Img the image configuration
/// \param blk block size of CUDA
/// \param gid grid size for CUDA
void proj(float* dimg, float* dcor, float* draw, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);

void proj_CPU(float* dimg, float* ddiff, float* ddraw, const FanEDGeo& FanGeo, const Image& Img);
void proj_CPU_OPENMP(float* dimg, float* ddif, float* draw, const FanEDGeo& FanGeo, const Image& Img);
//void proj_CPU(double* dimg, double* ddiff, double* ddraw, const FanEDGeo& FanGeo, const Image& Img);

//��ϵ�к����ʾͬ�������configuration���ڲ�ͬ����ͬʱ��ͶӰ�ĺ���Ӧ����patientRec �Ǹ�64����ع�
/// \brief This is the function used for 64 slices Fan Beam reconstruction simultaneously
/// \param dimg the reconstructed result that is 64 slices 
/// \param dcor the difference of the projection with weighting
/// \param draw the raw projection
/// \param FanGeo The Fan beam geometry with equal angle detector
/// \param Img the image configuration
/// \param blk The block size for CUDA
/// \param gid the grid size for CUDA
/// \param sliceNum the slice number for the volume
void proj(float* dimg, float* dcor, float* draw, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid, cuint sliceNum);

/// \brief This is the function used for 66 slices Fan Beam reconstruction simultaneously
/// \param dimg the reconstructed result that is 64 slices 
/// \param donp the difference of the projection with weighting
/// \param draw the raw projection
/// \param FanGeo The Fan beam geometry with equal angle detector
/// \param Img the image configuration
/// \param numPerSubSet number of projections per subset
/// \param subSetNum subset number
/// \param curSubSetIdx current subset Index
/// \param blk The block size for CUDA
/// \param gid the grid size for CUDA
/// \param sliceNum the slice number for the volume
/// \param streams The stream for multi GPU
void proj(float* dimg, float* donp, float* draw, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid, cuint sliceNum, const cudaStream_t& streams);



/// \brief This is the function used for 66 slices Fan Beam Equal Angular Detector Based reconstruction simultaneously with AIM 
/// \param dimg the reconstructed result that is 64 slices
/// \param donp the difference of the projection with weighting
/// \param draw the raw projection data
/// \param FanGeo FanEAGeo geometry with equal angle detector
/// \param Img the image configuration
/// \param numPerSubSet number of projections per subset
/// \param subSetNum subset number
/// \param curSubSetIdx current subset Index
/// \param blk The block size for CUDA
/// \param gid the grid size for CUDA
/// \param sliceNum the slice number for the volume
/// \param streams The stream for multi GPU
void proj_AIM(float* dimg, float* donp, float* draw, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid, cuint sliceNum, const cudaStream_t& streams);

/// \brief This is the function used for 66 slices Fan Beam Equal Angular Detector Based reconstruction simultaneously with AIM 
/// \param dimg the reconstructed result that is 64 slices
/// \param donp the difference of the projection with weighting
/// \param draw the raw projection data
/// \param FanGeo FanEAGeo geometry with equal angle detector
/// \param Img the image configuration
/// \param numPerSubSet number of projections per subset
/// \param subSetNum subset number
/// \param curSubSetIdx current subset Index
/// \param blk The block size for CUDA
/// \param gid the grid size for CUDA
/// \param sliceNum the slice number for the volume
/// \param streams The stream for multi GPU
void proj_AIM(double* dimg, double* donp, double* draw, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid, cuint sliceNum, const cudaStream_t& streams);






/// \brief This is the Golden standard Area Integral Based projection method in 2D, which is actually slow in GPU implementation
/// Each threads corresponds to one detector cell in one angle. For calculating this, it traversals the whole image to determine 
/// the bijection fan intersects with the current pixel or not, therefore it is slow.
/// \param dprj The projection result
/// \param dimg The image to be projected
/// \param FanGeo The Fan Geometry
/// \param Img The Image configuration
/// \param blk The block configuration for CUDA
/// \param gid The grid configuration for CUDA
void proj_AIMslow(float* dprj, float* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);

/// \brief This is the Golden standard Area Integral Based projection method in 2D, which is actually slow in GPU implementation
/// Each threads corresponds to one detector cell in one angle. For calculating this, it traversals the whole image to determine 
/// the bijection fan intersects with the current pixel or not, therefore it is slow.
/// \param dprj The projection result
/// \param dimg The image to be projected
/// \param FanGeo The Fan Geometry
/// \param Img The Image configuration
/// \param blk The block configuration for CUDA
/// \param gid The grid configuration for CUDA
void proj_AIMslow(double* dprj, double* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);




/// \brief This is the Area Integral Based projection method in 2D
/// \param dprj The projection result
/// \param dimg The Image to be projected
/// \param FanGeo Fan Beam Geometry with Equal Angular Detector
/// \param Img Image Configuration
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void proj_AIM(float* dprj, float* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);
/// \brief This is the Area Integral Based projection method in 2D
/// \param dprj The projection result
/// \param dimg The Image to be projected
/// \param FanGeo Fan Beam Geometry with Equal Angular Detector
/// \param Img Image Configuration
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void proj_AIM(double* dprj, double*dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);

/// \brief This is the Area Integral Based projection method in 2D, we get the projection difference with weight
/// \param dprj The projection result
/// \param draw The raw projection data
/// \param dimg The Image to be projected
/// \param FanGeo Fan Beam Geometry with Equal Angular Detector
/// \param Img Image Configuration
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void proj_AIM(float* dprj, float* draw, float* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);
/// \brief This is the Area Integral Based projection method in 2D, we get the projection difference with weight
/// \param dprj The projection result
/// \param draw The raw projection data
/// \param dimg The Image to be projected
/// \param FanGeo Fan Beam Geometry with Equal Angular Detector
/// \param Img Image Configuration
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void proj_AIM(double* dprj, double* draw, double* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);


/// \brief This is the Area Integral Based projection method in 2D, we get the projection difference with weight inside the subset
/// \param dimg The Image to be projected
/// \param donp The projection result
/// \param draw The raw projection data
/// \param FanGeo Fan Beam Geometry with Equal Angular Detector
/// \param Img Image Configuration
/// \param numPerSubSet How many projections are there in one subset
/// \param subSetNum How many subsets
/// \param curSubSetIdx current subset index
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void proj_AIM(float* dimg, float* donp, float* draw, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);
/// \brief This is the Area Integral Based projection method in 2D, we get the projection difference with weight inside the subset
/// \param dimg The Image to be projected
/// \param donp The projection result
/// \param draw The raw projection data
/// \param FanGeo Fan Beam Geometry with Equal Angular Detector
/// \param Img Image Configuration
/// \param numPerSubSet How many projections are there in one subset
/// \param subSetNum How many subsets
/// \param curSubSetIdx current subset index
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void proj_AIM(double* dimg, double* donp, double* draw, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);







/// \brief This is the Area Integral Based projection method in 2D, we get the projection with SHARED memory. The maximum detector resolution is 1023, the maximum project angle is 65535
/// \param dprj The projection result
/// \param dimg The raw projection data
/// \param FanGeo Fan Beam Geometry with Equal Angular Detector
/// \param Img Image Configuration
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void proj_AIM2(float* dprj, float* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);
/// \brief This is the Area Integral Based projection method in 2D, we get the projection with SHARED memory.
/// \param dprj The projection result
/// \param dimg The raw projection data
/// \param FanGeo Fan Beam Geometry with Equal Angular Detector
/// \param Img Image Configuration
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void proj_AIM2(double* dprj, double* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);



/// \brief This is the Area Integral Based projection method in 2D in CPU
/// \param prj The projection result
/// \param img The raw projection data
/// \param FanGeo FanEAGeo configuration with Equal Angular Detector
/// \param Img Image Configuration
void proj_AIM_CPU(float* prj, float* img, const FanEAGeo& FanGeo, const Image& Img);
/// \brief This is the Area Integral Based projection method in 2D in CPU
/// \param prj The projection result
/// \param img The raw projection data
/// \param FanGeo FanEAGeo configuration with Equal Angular Detector
/// \param Img Image Configuration
void proj_AIM_CPU(double* prj, double* img, const FanEAGeo& FanGeo, const Image& Img);

/// \brief This is the Area Integral Based projection method in 2D in CPU with OpenMP
/// \param prj The projection result
/// \param img The raw projection data
/// \param FanGeo FanEAGeo configuration with Equal Angular Detector
/// \param Img Image Configuration
void proj_AIM_CPU_OPENMP(float* prj, float* img, const FanEAGeo& FanGeo, const Image& Img);
/// \brief This is the Area Integral Based projection method in 2D in CPU with OpenMP
/// \param prj The projection result
/// \param img The raw projection data
/// \param FanGeo FanEAGeo configuration with Equal Angular Detector
/// \param Img Image Configuration
void proj_AIM_CPU_OPENMP(double* prj, double* img, const FanEAGeo& FanGeo, const Image& Img);





/// \brief This is the Area Integral Based projection method in 2D in CPU. Maybe the float precision is not enough
/// \param prj The projection result
/// \param img The raw projection data
/// \param FanGeo FanEDGeo configuration with Equal Angular Detector
/// \param Img Image Configuration
void proj_AIM_CPU(float* dprj, float* dimg, const FanEDGeo& FanGeo, const Image& Img);
/// \brief This is the Area Integral Based projection method in 2D in CPU
/// \param prj The projection result
/// \param img The raw projection data
/// \param FanGeo FanEDGeo configuration with Equal Angular Detector
/// \param Img Image Configuration
void proj_AIM_CPU(double* dprj, double* dimg, const FanEDGeo& FanGeo, const Image& Img);
/// \brief This is the Area Integral Based projection method in 2D in CPU with OpenMP. Maybe the float precision is not enough
/// \param prj The projection result
/// \param img The raw projection data
/// \param FanGeo FanEDGeo configuration with Equal Angular Detector
/// \param Img Image Configuration
void proj_AIM_CPU_OPENMP(float* dprj, float* dimg, const FanEDGeo& FanGeo, const Image& Img);
/// \brief This is the Area Integral Based projection method in 2D in CPU with OpenMP
/// \param prj The projection result
/// \param img The raw projection data
/// \param FanGeo FanEDGeo configuration with Equal Angular Detector
/// \param Img Image Configuration
void proj_AIM_CPU_OPENMP(double* dprj, double* dimg, const FanEDGeo& FanGeo, const Image& Img);
/// \brief This is the Area Integral Based projection method in 2D in GPU. Maybe the float precision is not enought
/// \param dprj The projection result
/// \param dimg The raw projection data
/// \param FanGeo FanEDGeo configuration with Equal Angular Detector
/// \param Img Image Configuration
/// \param blk Block configuration in CUDA
/// \param gid Grid configuration in CUDA
void proj_AIM(float* dprj, float* dimg, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);


/// \brief This is the Area Integral Based projection method in 2D in 3 GPUs
void proj_AIM_3GPUs(float* dprj, float* dimg, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);




/// \brief This is the Area Integral Based projection method in 2D in GPU
/// \param dprj The projection result
/// \param dimg The raw projection data
/// \param FanGeo FanEDGeo configuration with Equal Angular Detector
/// \param Img Image Configuration
/// \param blk Block configuration in CUDA
/// \param gid Grid configuration in CUDA
void proj_AIM(double* dprj, double* dimg, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);



/// \brief This is the Area Integral Based projection method in 2D in GPU, we get the differences with weighting . Maybe the float precision is not enough
/// \param dprj The projection result difference with weighting
/// \param dimg The raw projection data
/// \param FanGeo FanEDGeo configuration with Equal Angular Detector
/// \param Img Image Configuration
/// \param blk Block configuration in CUDA
/// \param gid Grid configuration in CUDA
void proj_AIM(float* dprj, float* draw, float* dimg, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);




/// \brief This is the Area Integral Based projection method in 2D in GPU, we get the differences with weighting . Maybe the float precision is not enough
/// \param dprj The projection result difference with weighting
/// \param dimg The raw projection data
/// \param FanGeo FanEDGeo configuration with Equal Angular Detector
/// \param Img Image Configuration
/// \param blk Block configuration in CUDA
/// \param gid Grid configuration in CUDA
void proj_AIM(double* dprj, double* draw, double* dimg, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);




/// \brief Area Integral Based projection method in 2D, we get the projection difference with weight inside the subset
/// \param dprj The difference with weight
/// \param draw The raw projection
/// \param dimg The image to be projected
/// \param FanGeo Fan Beam Geometry with Equal Angular Detector
/// \param Img Image Configuration
/// \param numPerSubSet How many projections are there in one subset
/// \param subSetNum How many subsets
/// \param curSubSetIdx current subset index
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void proj_AIM(float* dprj, float* draw, float* dimg, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);
/// \brief Area Integral Based projection method in 2D, we get the projection difference with weight inside the subset
/// \param dprj The difference with weight
/// \param draw The raw projection
/// \param dimg The image to be projected
/// \param FanGeo Fan Beam Geometry with Equal Angular Detector
/// \param Img Image Configuration
/// \param numPerSubSet How many projections are there in one subset
/// \param subSetNum How many subsets
/// \param curSubSetIdx current subset index
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void proj_AIM(double* dprj, double* draw, double* dimg, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);


/// \brief projection in one angle with weighting in Cone ED geometry, with GPU
/// \param dvol the volume
/// \param dprj the projection data to be projected
/// \param draw the raw projection data
/// \param ConeGeo The ED cone geometry
/// \param Vol The volume configuration
/// \param ang The angle index
/// \param blk The block configuration
/// \param gid The grid configuration

void proj_oneang(float* dvol, float* dprj, float* draw, const ConeEDGeo& ConeGeo, const Volume& Vol, const int angIdx, const dim3& blk, const dim3& gid);

void proj_V2(float* dvol, float* donp, float* draw, cuint angIdx, const bool angDir, const ConeEDGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid);
/// \brief backprojection in one angle
/// \param correction The correction projection data with weight
/// \param volume the Volume to be reconstructed
/// \param ConeGeo The cone ED geometry
/// \param Vol The Volume configuration
/// \param angIdx The angle index for backprojection
void back_proj(float* correction, float* volume, const ConeEDGeo ConeGeo, const Volume Vol, const unsigned int angIdx);




//NOTE: ���е�FOV����ֻ�ǰ���;�������κ�Ч��;
//NOTE: ���е�FOV����ֻ�ǰ���;�������κ�Ч��;
//NOTE: ���е�FOV����ֻ�ǰ���;�������κ�Ч��;
//////////////////////////////////////////////////////////////////////////
// Backprojection from one angle projection data
// donp: projection data from one angle;
// dimg: image;
// angIdx: angle index;
// angDir: positive or negative direction;
// FOV: do we only calculate inside the FOV;
// FanGeo: Fan Beam Geometry
// Img: Image configuration
// blk: 32,32
// gid: configured
//////////////////////////////////////////////////////////////////////////



/// \brief Back projection from one angle projection data
/// \param donp the projection in one angle to be backprojected
/// \param dimg the result
/// \param angIdx current angle index
/// \param angDir the direction for indexing
/// \param FOV reconstructed inside of FOV or not
/// \param FanGeo The Fan Beam Geometry Equal Angle detector
/// \param Img the image configuration
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dimg, cuint angIdx, const bool angDir, const bool FOV, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);
/// \brief Back projection from one angle projection data
/// \param donp the projection in one angle to be backprojected
/// \param dimg the result
/// \param angIdx current angle index
/// \param angDir the direction for indexing
/// \param FOV reconstructed inside of FOV or not
/// \param FanGeo The Fan Beam Geometry Equal Distance detector
/// \param Img the image configuration
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dimg, cuint angIdx, const bool angDir, const bool FOV, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);
/// \brief Back projection from one angle projection data
/// \param donp the projection in one angle to be backprojected
/// \param dvol the result
/// \param angIdx current angle index
/// \param angDir the direction for indexing
/// \param FOV reconstructed inside of FOV or not
/// \param ConeGeo The Cone Beam Geometry Equal Angle detector
/// \param Vol the volume configuration
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dvol, cuint angIdx, const bool angDir, const bool FOV, const ConeEAGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid);





//////////////////////////////////////////////////////////////////////////
// Back projection inside the subset with given subset index
// donp: the projection data inside the subset 
// dimg: image 
// angDir: we use positive or negative index direction
// FOV: do we only calculate inside the FOV
// FanGeo: Fan Beam Geometry
// Img: Image configuration
// numPerSubSet: how many angles are there inside the subset 
// subSetNum: how many subsets 
// curSubSetIdx: current subset index;
// blk: 32,32
// gid: configured
//////////////////////////////////////////////////////////////////////////


/// \brief Back projection from subsets projection data
/// \param donp the projection in subsets angle to be back-projected
/// \param dimg the result
/// \param angDir the direction for indexing
/// \param FOV reconstructed inside of FOV or not
/// \param FanGeo The Fan Beam Geometry Equal Angle detector
/// \param Img the image configuration
/// \param numPerSubSet number of projections in one subset
/// \param subSetNum how many subsets are divided in projection
/// \param curSubSetIdx current subset index
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dimg, const bool angDir, const bool FOV, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);

/// \brief Back projection from subsets projection data
/// \param donp the projection in subsets angle to be back-projected
/// \param dimg the result
/// \param angDir the direction for indexing
/// \param FOV reconstructed inside of FOV or not
/// \param FanGeo The Fan Beam Geometry Equal Distance detector
/// \param Img the image configuration
/// \param numPerSubSet number of projections in one subset
/// \param subSetNum how many subsets are divided in projection
/// \param curSubSetIdx current subset index
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dimg, const bool angDir, const bool FOV, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);

/// \brief Back projection from subsets projection data
/// \param donp the projection in subsets angle to be back-projected
/// \param dvol the result
/// \param angDir the direction for indexing
/// \param FOV reconstructed inside of FOV or not
/// \param ConeGeo The Cone Beam Geometry Equal Distance detector
/// \param Vol the volume configuration
/// \param numPerSubSet number of projections in one subset
/// \param subSetNum how many subsets are divided in projection
/// \param curSubSetIdx current subset index
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dvol, const bool angDir, const bool FOV, const ConeEAGeo& ConeGeo, const Volume& Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);














//////////////////////////////////////////////////////////////////////////
// Backprojection from the subset with given subset index
// donp: the projection data from subset
// dimg: image
// FanGeo: Fan Beam Geometry
// Img: Image configuration
// numPerSubSet: how many angles are there inside the subset 
// subSetNum: how many subsets 
// curSubSetIdx: current subset index;
// blk: 32,32
// gid: configured
//////////////////////////////////////////////////////////////////////////



/// \brief Backprojection from the subset with given subset index
/// \param donp the projection data from subset
/// \param dimg the result
/// \param FanGeo Fan Beam geometry with equal angle detector
/// \param Img Image configuration
/// \param numPerSubSet number of projections per subset
/// \param subSetNum How many subsets are there
/// \param curSubSetIdx current subset index
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dimg, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);

/// \brief Backprojection from the subset with given subset index
/// \param donp the projection data from subset
/// \param dimg the result
/// \param FanGeo Fan Beam geometry with equal distance detector
/// \param Img Image configuration
/// \param numPerSubSet number of projections per subset
/// \param subSetNum How many subsets are there
/// \param curSubSetIdx current subset index
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dimg, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);

/// \brief Backprojection from the subset with given subset index
/// \param donp the projection data from subset
/// \param dvol the result
/// \param ConeGeo Cone Beam geometry with equal distance detector
/// \param Vol Volume configuration
/// \param numPerSubSet number of projections per subset
/// \param subSetNum How many subsets are there
/// \param curSubSetIdx current subset index
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dvol, const ConeEAGeo& ConeGeo, const Volume& Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);




//////////////////////////////////////////////////////////////////////////
// Backprojection from all the angles
// donp: projection data
// dimg: image;
// FanGeo: Fan Beam Geometry
// Img: Image configuration
// blk: 32,32
// gid: configured
//////////////////////////////////////////////////////////////////////////


/// \brief Backprojection from all angles
/// \param donp the raw projection data 
/// \param dimg the result
/// \param FanGeo Fan Beam geometry with equal angle detector
/// \param Img Image configuration
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);

/// \brief Backprojection from all angles
/// \param donp the raw projection data 
/// \param dimg the result
/// \param FanGeo Fan Beam geometry with equal distance detector
/// \param Img Image configuration
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dimg, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);

/// \brief Backprojection from all angles
/// \param donp the raw projection data 
/// \param dvol the result
/// \param CoeGeo Cone Beam geometry with equal angle detector
/// \param Vol Volume configuration
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dvol, const ConeEAGeo& CoeGeo, const Volume& Vol, const dim3& blk, const dim3& gid);
//////////////////////////////////////////////////////////////////////////
// Backprojection from one angle 
// donp: projection data from one angle;
// dimg: image;
// dmsk: mask;
// angIdx: current angle index;
// angDir: index from positive or negative direction;
// FanGeo: Fan Beam Geometry
// Img: Image configuration
// blk: 32,32
// gid: configured
//////////////////////////////////////////////////////////////////////////

/// \brief Backprojection from one angle
/// \param donp the projection data in one angle;
/// \param dimg the image
/// \param dmsk device mask
/// \param angIdx angle index
/// \param angDir indexing direction from small to large or large to small
/// \param FanGeo Fan Beam geometry with equal angle detector
/// \param Img Image configuration
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, cuint angIdx, const bool angDir, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);

/// \brief Backprojection from one angle
/// \param donp the projection data in one angle;
/// \param dimg the image
/// \param dmsk device mask
/// \param angIdx angle index
/// \param angDir indexing direction from small to large or large to small
/// \param FanGeo Fan Beam geometry with equal distance detector
/// \param Img Image configuration
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, cuint angIdx, const bool angDir, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);
/// \brief Backprojection from one angle
/// \param donp the projection data in one angle;
/// \param dvol the volume
/// \param dmsk device mask
/// \param angIdx angle index
/// \param angDir indexing direction from small to large or large to small
/// \param ConeGeo Fan Beam geometry with equal distance detector
/// \param Vol Volume configuration
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dvol, float* dmsk, cuint angIdx, const bool angDir, const ConeEAGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid);



//////////////////////////////////////////////////////////////////////////
// Backprojection from subset angles
// donp: projection inside the subsets;
// dimg: image;
// dmsk: mask;
// angDir: index direction from positive or negative;
// FanGeo: Fan Beam Geometry
// Img: Image configuration
// numPerSubSet: How many projections are inside the subset;
// subSetNum: How many subsets are distributed;
// curSubSetIdx: current subset Index
// blk: 32,32
// gid: configured
//////////////////////////////////////////////////////////////////////////
/// \brief Backprojection from subset
/// \param donp the projection data in the subset
/// \param dimg the volume
/// \param dmsk device mask
/// \param angDir indexing direction from small to large or large to small
/// \param FanGeo Fan Beam geometry with equal angle detector
/// \param Img Image configuration
/// \param numPerSubSet number of projections per subset
/// \param subSetNum How many subsets are in projections
/// \param curSubSetIdx current subset 
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, const bool angDir, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);
/// \brief Backprojection from subset
/// \param donp the projection data in the subset
/// \param dimg the volume
/// \param dmsk device mask
/// \param angDir indexing direction from small to large or large to small
/// \param FanGeo Fan Beam geometry with equal distance detector
/// \param Img Image configuration
/// \param numPerSubSet number of projections per subset
/// \param subSetNum How many subsets are in projections
/// \param curSubSetIdx current subset 
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, const bool angDir, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);
/// \brief Backprojection from subset
/// \param donp the projection data in the subset
/// \param dvol the volume
/// \param dmsk device mask
/// \param angDir indexing direction from small to large or large to small
/// \param ConeGeo Cone Beam geometry with equal angle detector
/// \param Vol Image configuration
/// \param numPerSubSet number of projections per subset
/// \param subSetNum How many subsets are in projections
/// \param curSubSetIdx current subset 
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dvol, float* dmsk, const bool angDir, const ConeEAGeo& ConeGeo, const Volume& Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);
//////////////////////////////////////////////////////////////////////////
// Backprojection from subset angles
// donp: projection inside the subsets;
// dimg: image;
// dmsk: mask;
// FanGeo: Fan Beam Geometry
// Img: Image configuration
// numPerSubSet: How many projections are inside the subset;
// subSetNum: How many subsets are distributed;
// curSubSetIdx: current subset Index
// blk: 32,32
// gid: configured
//////////////////////////////////////////////////////////////////////////


/// \brief Backprojection from subset
/// \param donp the projection data in the subset
/// \param dimg the volume
/// \param dmsk device mask
/// \param FanGeo Fan Beam geometry with equal angle detector
/// \param Img Image configuration
/// \param numPerSubSet number of projections per subset
/// \param subSetNum How many subsets are in projections
/// \param curSubSetIdx current subset 
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);
/// \brief Backprojection from subset
/// \param donp the projection data in the subset
/// \param dimg the volume
/// \param dmsk device mask
/// \param FanGeo Fan Beam geometry with equal distance detector
/// \param Img Image configuration
/// \param numPerSubSet number of projections per subset
/// \param subSetNum How many subsets are in projections
/// \param curSubSetIdx current subset 
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);

/// \brief Backprojection from subset
/// \param donp the projection data in the subset
/// \param dimg the volume
/// \param dmsk device mask
/// \param ConeGeo Cone Beam geometry with equal angle detector
/// \param Vol Volume configuration
/// \param numPerSubSet number of projections per subset
/// \param subSetNum How many subsets are in projections
/// \param curSubSetIdx current subset 
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, const ConeEAGeo& ConeGeo, const Volume& Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);

/// \brief Backprojection from subset
/// \param donp the projection data in the subset
/// \param dimg the volume
/// \param dmsk device mask
/// \param ConeGeo Cone Beam geometry with equal distance detector
/// \param Vol Volume configuration
/// \param numPerSubSet number of projections per subset
/// \param subSetNum How many subsets are in projections
/// \param curSubSetIdx current subset 
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, const ConeEDGeo& ConeGeo, const Volume& Vol, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);
//////////////////////////////////////////////////////////////////////////
// Backprojection from all angles
// donp: projection inside the subsets;
// dimg: image;
// dmsk: mask;
// FanGeo: Fan Beam Geometry
// Img: Image configuration
// blk: 32,32
// gid: configured
//////////////////////////////////////////////////////////////////////////

/// \brief Backprojection from all angles
/// \param donp the projection data in the subset
/// \param dimg the volume
/// \param dmsk device mask
/// \param FanGeo Fan Beam geometry with equal angle detector
/// \param Img Image configuration
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);
void bakproj_PIXEL_CPU(float* donp, float* dimg, float* dmsk, const FanEAGeo& FanGeo, const Image& Img);
void bakproj_PIXEL_CPU_OPENMP(float* donp, float* dimg, float* dmsk, const FanEAGeo& FanGeo, const Image& Img);



/// \brief Backprojection from all angles
/// \param donp the projection data in the subset
/// \param dimg the volume
/// \param dmsk device mask
/// \param FanGeo Fan Beam geometry with equal distance detector
/// \param Img Image configuration
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);
/// \brief Backprojection from all angles
/// \param donp the projection data in the subset
/// \param dimg the volume
/// \param dmsk device mask
/// \param ConeGeo Cone Beam geometry with equal angle detector
/// \param Vol Volume configuration
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_PIXEL(float* donp, float* dimg, float* dmsk, const ConeEAGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid);










//��ϵ�к����ʾͬ�������configuration���ڲ�ͬ����ͬʱ��ͶӰ�ĺ���Ӧ����patientRec �Ǹ�64����ع�
/// \brief Backprojection for GPU with 64 slices
/// \param donp the projection data for back projection
/// \param dimg the image result
/// \param FanGeo Fan Beam Geometry for equal angle geometry
/// \param Img Image configuration
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
/// \param sliceNum slices number
void bakproj_PIXEL(float* donp, float* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid, cuint sliceNum);

/// \brief Backprojection for 3 GPUs for 66 slices
/// \param donp the projection data for back projection
/// \param dimg the image result
/// \param FanGeo Fan Beam Geometry for equal angle geometry
/// \param Img Image configuration
/// \param numPerSubSet number of projections per subset
/// \param subSetNum How many subsets are there
/// \param curSubIdx current subset index
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
/// \param sliceNum slices number
/// \param streams the streams for different GPUs
void bakproj_PIXEL(float* donp, float* dimg, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubIdx, const dim3& blk, const dim3& gid, cuint sliceNum, const cudaStream_t& streams);



//////////////////////////////////////////////////////////////////////////
// Backprojection from one angle
// This is the boxed based model

/// \brief Boxed based back-projection from one angle
/// \param donp the projection data for back-projection
/// \param dvol the volume result
/// \param dmsk the mask for the image
/// \param angIdx the angle index
/// \param ConeGeo Cone Beam Geometry for equal angle geometry
/// \param Vol Volume configuration
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_BOXED(float* donp, float* dvol, float* dmsk, cuint angIdx, const ConeEDGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid);

/// \brief Boxed based back-projection from one angle
/// \param donp the projection data for back-projection
/// \param dvol the volume result
/// \param angIdx the angle index
/// \param ConeGeo Cone Beam Geometry for equal angle geometry
/// \param Vol Volume configuration
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_BOXED(float* donp, float* dvol, cuint angIdx, const ConeEDGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid);

/// \brief Boxed based back-projection from all angles
/// \param donp the projection data for back-projection
/// \param dvol the volume result
/// \param ConeGeo Cone Beam Geometry for equal angle geometry
/// \param Vol Volume configuration
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_BOXED(float* donp, float* dvol, const ConeEDGeo& ConeGeo, const Volume& Vol, const dim3& blk, const dim3& gid);

// Backprojection from all angle
// This is the boxed based model

/// \brief NOT FINISHED YET!! Boxed based back-projection from all angles
/// \param dprj the projection data for back-projection
/// \param dimg the image result
/// \param FanGeo Fan Beam Geometry for equal angle geometry
/// \param Img Image configuration
/// \param blk block size for CUDA
/// \param gid grid size for CUDA
void bakproj_BOXED(float* dprj, float* dimg, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);


/// \brief Boxed based back-projection from all angles in CPU
/// \param prj the projection data for back-projection
/// \param img the image result
/// \param FanGeo Fan Beam Geometry for equal angle geometry
/// \param Img Image configuration
void bakproj_BOXED_CPU(float* prj, float* img, const FanEAGeo& FanGeo, const Image& Img);

/// \brief Boxed based back-projection from all angles in CPU with OpenMP
/// \param prj the projection data for back-projection
/// \param img the image result
/// \param FanGeo Fan Beam Geometry for equal angle geometry
/// \param Img Image configuration
void bakproj_BOXED_CPU_OPENMP(float* dprj, float* dimg, const FanEAGeo& FanGeo, const Image& Img);




/// \brief The Area Integral based backprojection from certain angle in CPU
/// \param dprj the data for backprojection
/// \param dimg the image for backprojection, it is the result
/// \param angIdx the angle index of the projection data
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
void bakproj_AIM_CPU(float* dprj, float* dimg, cuint angIdx, const FanEAGeo& FanGeo, const Image& Img);
/// \brief The Area Integral based backprojection from certain angle in CPU
/// \param dprj the data for backprojection
/// \param dimg the image for backprojection, it is the result
/// \param angIdx the angle index of the projection data
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
void bakproj_AIM_CPU(double* dprj, double* dimg, cuint angIdx, const FanEAGeo& FanGeo, const Image& Img);
/// \brief The Area Integral based backprojection from all angles in CPU with OpenMP
/// \param dprj the data for backprojection
/// \param dimg the image for backprojection, it is the result
/// \param angIdx the angle index of the projection data
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
void bakproj_AIM_CPU(float* dprj, float* dimg, const FanEAGeo& FanGeo, const Image& Img);
/// \brief The Area Integral based backprojection from all angles in CPU with OpenMP
/// \param dprj the data for backprojection
/// \param dimg the image for backprojection, it is the result
/// \param angIdx the angle index of the projection data
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
void bakproj_AIM_CPU(double* dprj, double* dimg, const FanEAGeo& FanGeo, const Image& Img);
/// \brief The Area Integral based backprojection from certain angle in CPU with OpenMP
/// \param dprj the data for backprojection
/// \param dimg the image for backprojection, it is the result
/// \param angIdx the angle index of the projection data
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
void bakproj_AIM_CPU_OPENMP(float* dprj, float* dimg, cuint angIdx, const FanEAGeo& FanGeo, const Image& Img);
/// \brief The Area Integral based backprojection from certain angle in CPU with OpenMP
/// \param dprj the data for backprojection
/// \param dimg the image for backprojection, it is the result
/// \param angIdx the angle index of the projection data
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
void bakproj_AIM_CPU_OPENMP(double* dprj, double* dimg, cuint angIdx, const FanEAGeo& FanGeo, const Image& Img);
/// \brief The Area Integral based backprojection from all angles in CPU with OpenMP
/// \param dprj the data for backprojection
/// \param dimg the image for backprojection, it is the result
/// \param angIdx the angle index of the projection data
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
void bakproj_AIM_CPU_OPENMP(float* dprj, float* dimg, const FanEAGeo& FanGeo, const Image& Img);
/// \brief The Area Integral based backprojection from all angles in CPU with OpenMP
/// \param dprj the data for backprojection
/// \param dimg the image for backprojection, it is the result
/// \param angIdx the angle index of the projection data
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
void bakproj_AIM_CPU_OPENMP(double* dprj, double* dimg, const FanEAGeo& FanGeo, const Image& Img);



/// \brief The Area Integral based backprojection from all angle in GPU 
/// \param dprj the data in device for backprojection 
/// \param dimg the image in device for backprojection, it is the result
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void bakproj_AIM(float* dimg, float* dprj, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);

/// \brief The Area Integral based backprojection from all angle in GPU 
/// \param dprj the data in device for backprojection 
/// \param dimg the image in device for backprojection, it is the result
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void bakproj_AIM(double* dimg, double* dprj, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);

/// \brief The Area Integral based backprojection from certain angle in GPU 
/// \param dimg the image in device for backprojection, it is the result
/// \param dprj the projection data for backprojection in device
/// \param angIdx the angle index of the projection data
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void bakproj_AIM(float* dimg, float* dprj, const unsigned int angIdx, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);
/// \brief The Area Integral based backprojection from certain angle in GPU 
/// \param dimg the image in device for backprojection, it is the result
/// \param dprj the projection data for backprojection in device
/// \param angIdx the angle index of the projection data
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void bakproj_AIM(double* dimg, double* dprj, const unsigned int angIdx, const FanEAGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);

/// \brief The Area Integral based backprojection from a subset angles in GPU 
/// \param donp the projection data for backprojection in device
/// \param dimg the image in device for backprojection, it is the result
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
/// \param numPerSubSet number of projection in one subset
/// \param subSetNum How many subsets
/// \param curSubSetIdx current subset index
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void bakproj_AIM(float* donp, float* dimg, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);
/// \brief The Area Integral based backprojection from a subset angles in GPU 
/// \param donp the projection data for backprojection in device
/// \param dimg the image in device for backprojection, it is the result
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
/// \param numPerSubSet number of projection in one subset
/// \param subSetNum How many subsets
/// \param curSubSetIdx current subset index
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void bakproj_AIM(double* donp, double* dimg, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);










/// \brief The Area Integral based backprojection from certain angle in CPU
/// \param dprj the data for backprojection
/// \param dimg the image for backprojection, it is the result
/// \param angIdx the angle index of the projection data
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
void bakproj_AIM_CPU(float* dprj, float* dimg, cuint angIdx, const FanEDGeo& FanGeo, const Image& Img);
/// \brief The Area Integral based backprojection from certain angle in CPU
/// \param dprj the data for backprojection
/// \param dimg the image for backprojection, it is the result
/// \param angIdx the angle index of the projection data
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
void bakproj_AIM_CPU(double* dprj, double* dimg, cuint angIdx, const FanEDGeo& FanGeo, const Image& Img);
/// \brief The Area Integral based backprojection from all angles in CPU with OpenMP
/// \param dprj the data for backprojection
/// \param dimg the image for backprojection, it is the result
/// \param angIdx the angle index of the projection data
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
void bakproj_AIM_CPU(float* dprj, float* dimg, const FanEDGeo& FanGeo, const Image& Img);
/// \brief The Area Integral based backprojection from all angles in CPU with OpenMP
/// \param dprj the data for backprojection
/// \param dimg the image for backprojection, it is the result
/// \param angIdx the angle index of the projection data
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
void bakproj_AIM_CPU(double* dprj, double* dimg, const FanEDGeo& FanGeo, const Image& Img);
/// \brief The Area Integral based backprojection from certain angle in CPU with OpenMP
/// \param dprj the data for backprojection
/// \param dimg the image for backprojection, it is the result
/// \param angIdx the angle index of the projection data
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
void bakproj_AIM_CPU_OPENMP(float* dprj, float* dimg, cuint angIdx, const FanEDGeo& FanGeo, const Image& Img);
/// \brief The Area Integral based backprojection from certain angle in CPU with OpenMP
/// \param dprj the data for backprojection
/// \param dimg the image for backprojection, it is the result
/// \param angIdx the angle index of the projection data
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
void bakproj_AIM_CPU_OPENMP(double* dprj, double* dimg, cuint angIdx, const FanEDGeo& FanGeo, const Image& Img);
/// \brief The Area Integral based backprojection from all angles in CPU with OpenMP
/// \param dprj the data for backprojection
/// \param dimg the image for backprojection, it is the result
/// \param angIdx the angle index of the projection data
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
void bakproj_AIM_CPU_OPENMP(float* dprj, float* dimg, const FanEDGeo& FanGeo, const Image& Img);
/// \brief The Area Integral based backprojection from all angles in CPU with OpenMP
/// \param dprj the data for backprojection
/// \param dimg the image for backprojection, it is the result
/// \param angIdx the angle index of the projection data
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
void bakproj_AIM_CPU_OPENMP(double* dprj, double* dimg, const FanEDGeo& FanGeo, const Image& Img);
/// \brief The Area Integral based backprojection from all angle in GPU 
/// \param dprj the data in device for backprojection 
/// \param dimg the image in device for backprojection, it is the result
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void bakproj_AIM(float* dimg, float* dprj, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);
/// \brief The Area Integral based backprojection from all angle in GPU 
/// \param dprj the data in device for backprojection 
/// \param dimg the image in device for backprojection, it is the result
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void bakproj_AIM(double* dimg, double* dprj, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);
/// \brief The Area Integral based backprojection from certain angle in GPU 
/// \param dimg the image in device for backprojection, it is the result
/// \param dprj the projection data for backprojection in device
/// \param angIdx the angle index of the projection data
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void bakproj_AIM(float* dimg, float* dprj, const unsigned int angIdx, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);
/// \brief The Area Integral based backprojection from certain angle in GPU 
/// \param dimg the image in device for backprojection, it is the result
/// \param dprj the projection data for backprojection in device
/// \param angIdx the angle index of the projection data
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void bakproj_AIM(double* dimg, double* dprj, const unsigned int angIdx, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);
/// \brief The Area Integral based backprojection from a subset angles in GPU 
/// \param donp the projection data for backprojection in device
/// \param dimg the image in device for backprojection, it is the result
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
/// \param numPerSubSet number of projection in one subset
/// \param subSetNum How many subsets
/// \param curSubSetIdx current subset index
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void bakproj_AIM(float* donp, float* dimg, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);
/// \brief The Area Integral based backprojection from a subset angles in GPU 
/// \param donp the projection data for backprojection in device
/// \param dimg the image in device for backprojection, it is the result
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
/// \param numPerSubSet number of projection in one subset
/// \param subSetNum How many subsets
/// \param curSubSetIdx current subset index
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
void bakproj_AIM(double* donp, double* dimg, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);

/// \brief The Area Integral based backprojection from a subset angles in GPU 
/// \param donp the projection data for backprojection in device
/// \param dimg the image in device for backprojection, it is the result
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
/// \param numPerSubSet number of projection in one subset
/// \param subSetNum How many subsets
/// \param curSubSetIdx current subset index
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
/// \param sliceNum is defined as 22 slices
/// \param streams which stream is applied
void bakproj_AIM(float* dprj, float* dimg, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubIdx, const dim3& blk, const dim3& gid, cuint sliceNum, const cudaStream_t& streams);
/// \brief The Area Integral based backprojection from a subset angles in GPU 
/// \param donp the projection data for backprojection in device
/// \param dimg the image in device for backprojection, it is the result
/// \param FanGeo Fan Beam Geometry with equal angle detector
/// \param Img Image configuration
/// \param numPerSubSet number of projection in one subset
/// \param subSetNum How many subsets
/// \param curSubSetIdx current subset index
/// \param blk Block Size for CUDA
/// \param gid Grid Size for CUDA
/// \param sliceNum is defined as 22 slices
/// \param streams which stream is applied
void bakproj_AIM(double* dprj, double* dimg, const FanEAGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubIdx, const dim3& blk, const dim3& gid, cuint sliceNum, const cudaStream_t& streams);



/// \brief The projection function applied in 2D SIR reconstruction in GPU
/// \param dimg The image to be reconstructed
/// \param donp The weighted projection difference
/// \param draw The raw projection data
/// \param dpho The photon distribution
/// \param maxY The maximum number for the photon on the detector
/// \param FanGeo FanBeam geometry in FanEDGeo configuration
/// \param Img The Image configuration
/// \param numPerSubSet The number of projections in one subset
/// \param subSetNum The subset number for reconstruction;
/// \param curSubSetIdx The current subset index
/// \param blk projection CUDA block configuration
/// \param gid projection CUDA grid configuration 
void proj(float* dimg, float* donp, float* draw, float* dpho, const float maxY, const FanEDGeo& FanGeo, const Image& Img, cuint numPerSubSet, cuint subSetNum, cuint curSubSetIdx, const dim3& blk, const dim3& gid);




void proj_DEMO18v4_2D(float* dimg, float* draw, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);
void bakproj_PIXEL_DEMO18v4_2D(float* donp, float* dimg, float* dmsk, const FanEDGeo& FanGeo, const Image& Img, const dim3& blk, const dim3& gid);


void proj_DEMO17(float* dvol, float* dprj, float* draw, const ConeEDGeo& ConeGeo, const Volume& Vol, const int angIdx, const dim3& blk, const dim3& gid);
void back_proj_DEMO17(float* correction, float* volume, const ConeEDGeo& ConeGeo, const Volume& Vol, const int angIdx, const dim3& blockSize, const dim3& gridSize);
void back_proj_DEMO18v3(float* donp, float* dimg, const ConeEDGeo& ConeGeo, const Volume& Vol, const int angIdx, const dim3& blockSize, const dim3& gridSize);


