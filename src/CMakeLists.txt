cmake_minimum_required(VERSION 3.22)
STRING(REPLACE "-O2" "-Od" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
FIND_PACKAGE(CUDA)
MESSAGE( STATUS "CTLAB - project is being generated" )

#For OpenMP
FIND_PACKAGE(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_FOUND)

##For MKL
#find_package(MKL REQUIRED)
#include_directories(${MKL_INCLUDE_DIRS})
#link_directories(${MKL_LIBRARIES})


option(CUDA_DEVICE_DEBUGGING
  "On to enable debugging device code."
  OFF
  )  

if (CUDA_DEVICE_DEBUGGING)
  set(DEVICE_DEBUGGING -G0)
endif (CUDA_DEVICE_DEBUGGING)



SET(CTLAB_HEADERS
	utilities.cuh    	# The head files for DD library
	utilities.hpp    	# The head files for DD library
	DD3_GPU_Proj.h   	# The Branchless DD projection
	DD3_GPU_Back.h   	# The Branchless DD backprojection
	DD3_GPU_demo.h		# Demo of using the DD3 lib
	FastIterativeShrinkageAlgorithm.h # Fast Iterative Shrinkage Algorithm
	SARTWeighting.h         # SART weighting algorithm
	DD3_GPU_recon.h         # GPU based reconstruction methods
	DDM_proj.h       	# The conventional DD projection/backprojeciton head file
	SF_proj.hpp	        # The head file interface for the Separable footprint implementation
	differentialFuncs.hpp 	# The interface of the differential operation
	useful.hpp		# head files useful functions for CT library
	projbackproj.hpp 	# The interface of conventional projection and backprojection 
	threeGPUHouseholder.hpp # The interface of GPU Based Householder transform, it requires three GPUs
	spectralCTscan.hpp      # The interface of spectral CT projection simulation
	#DD3_demo.hpp	# The interface of the demo of the Branchless DD
	readwrite.hpp           # Read and write the raw data (in GPU or CPU) from/to disk
	FastMatrixVectorMultiplication.hpp # Fast Matrix-Vector multiplication function in GPU
	genModifiedPhantom.hpp  # Generate modified Shepp-Logan phantom in GPU
	Ellipsoid.h             # kernel function for generation Shepp-Logan projection
	genSheppLoganProjection.h # Generate Shepp Logan phantom based projection
	genSysMatrix.hpp        # Generate system matrix with SIDDON/AIM model
	multiSlices_ker.cuh     # multi slices based projection and backprojection
	ProjectionModelMap.h
	BackprojectionModelMap.h
	Geometry.h
	
	
)

SET(CTLAB_SRC 
    DD3_GPU_Proj.cu 	# implementation of branchless DD projection
    DD3_GPU_Back.cu 	# implementation of branchless DD backprojection
	DD3_GPU_demo.cu   	# Demo of using the DD3 lib
	DD3_GPU_recon.cu         # GPU based reconstruction methods
	DDM_proj.cpp    	# conventional DD projection/backprojection in CPU
	DDM_proj.cu		# conventional DD projectoin/backprojection in GPU
	SF_proj.cpp		# The CPU implementation of the SF projection and backprojection
	differentialFuncs.cpp   # The CPU implementation of the differential operations
	differentialFuncs.cu    # The GPU implementation of the differential operations
	useful.cpp		# useful functions for CT library
	useful.cu		# useful functions for CT library
	projbackproj.cpp 	# The conventional projection and backprojection in CPU
	projbackproj.cu 	# The conventional projection and backprojection in GPU
	threeGPUHouseholder.cu  # The interface of GPU Based Householder transform, it requires three GPUs
	spectralCTscan.cpp 	# The spectral CT projection simulation in CPU
	spectralCTscan.cu 	# The spectral CT projection simulation in GPU
#	DD3_demo.cu	# The interface of the demo of the Branchless DD
	FastMatrixVectorMultiplication.cu # Fast Matrix-Vector multiplication function in GPU
	readwrite.cpp           # Read and write the raw data (in GPU or CPU) from/to disk
	genModifiedPhantom.cu   # Generate modified Shepp-Logan phantom in GPU
	genSheppLoganProjection.cpp # Generate Shepp Logan phantom based projection
	genSysMatrix.cpp        # Generate system matrix with SIDDON/AIM model
	multiSlices_ker.cu      # multi slices based projection and backprojection
	ProjectionModelMap.cpp
	BackprojectionModelMap.cpp
	Geometry.cpp
)



if (WIN32)
  if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    set( LIB_PATH $ENV{CUDA_SDK_ROOT_DIR}/common/lib/x64/ )
    
    set( GLEW_NAME glew64 )
  else( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    set( LIB_PATH $ENV{CUDA_SDK_ROOT_DIR}/common/lib/win32/ )
    
    set( GLEW_NAME glew32 )
  endif( CMAKE_SIZEOF_VOID_P EQUAL 8 )
else (WIN32)
  if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    set( LIB_PATH $ENV{CUDA_SDK_ROOT_DIR}/common/lib/linux/x86_64/ )
    
    set( GLEW_NAME GLEW_x86_64 )
  else( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    set( LIB_PATH $ENV{CUDA_SDK_ROOT_DIR}/common/lib/linux/i686/ )
    
    set( GLEW_NAME GLEW )  
  endif( CMAKE_SIZEOF_VOID_P EQUAL 8 )
endif (WIN32)

include_directories($ENV{CUDA_SDK_ROOT_DIR}/common/inc)

link_directories("C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.5/bin")

set(GENCODE -O3 --use_fast_math -Xcompiler)

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

CUDA_ADD_LIBRARY(CTLAB ${CTLAB_SRC} ${CTLAB_HEADERS} SHARED OPTIONS ${GENCODE} ${DEVICE_DEBUGGING})

target_link_libraries(CTLAB "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.5/lib/x64/cublas.lib")
target_link_libraries(CTLAB "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.5/lib/x64/cusparse.lib")
