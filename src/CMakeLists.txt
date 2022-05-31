# we need at least cmake version 3.1
cmake_minimum_required(VERSION 3.1)

project(qmc_worm CXX)

set(CMAKE_BUILD_TYPE Release)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the project relies on the ALPSCore package. If ALPSCore is not found
# automatically, specify its location using:
# export ALPSCore_DIR=/location/to/ALPSCORE/
find_package(ALPSCore 2.0  REQUIRED mc accumulators params)

# check if Eigen is present:
if (ALPSCore_HAS_EIGEN_VERSION)
  message(STATUS "ALPSCore provides Eigen, version ${ALPSCore_HAS_EIGEN_VERSION}")
else()
  message("WARNING: ALPSCore does not provide Eigen, the code won't compile!")
endif()

add_executable(${PROJECT_NAME} worm.run.cpp worm.cpp worm.update.cpp worm.output.cpp)


# Use ALPSCore_LIBRARIES variable to link to ALPSCore 
target_link_libraries(${PROJECT_NAME} ${ALPSCore_LIBRARIES})
