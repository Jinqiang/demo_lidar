cmake_minimum_required(VERSION 2.4.6)
include($ENV{ROS_ROOT}/core/rosbuild/rosbuild.cmake)

# Set the build type.  Options are:
#  Coverage       : w/ debug symbols, w/o optimization, w/ code-coverage
#  Debug          : w/ debug symbols, w/o optimization
#  Release        : w/o debug symbols, w/ optimization
#  RelWithDebInfo : w/ debug symbols, w/ optimization
#  MinSizeRel     : w/o debug symbols, w/ optimization, stripped binaries
#set(ROS_BUILD_TYPE RelWithDebInfo)

rosbuild_init()

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

#uncomment if you have defined messages
#rosbuild_genmsg()
#uncomment if you have defined services
#rosbuild_gensrv()

#common commands for building c++ executables and libraries
#rosbuild_add_library(${PROJECT_NAME} src/example.cpp)
#target_link_libraries(${PROJECT_NAME} another_library)
#rosbuild_add_boost_directories()
#rosbuild_link_boost(${PROJECT_NAME} thread)
#rosbuild_add_executable(example examples/example.cpp)
#target_link_libraries(example ${PROJECT_NAME})

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories("${PROJECT_SOURCE_DIR}/include")

add_subdirectory(isamlib)
link_libraries(isamlib)

rosbuild_add_executable(featureTracking src/featureTracking.cpp)
rosbuild_add_executable(visualOdometry src/visualOdometry.cpp)
rosbuild_add_executable(bundleAdjust src/bundleAdjust.cpp)
rosbuild_add_executable(processDepthmap src/processDepthmap.cpp)
rosbuild_add_executable(stackDepthPoint src/stackDepthPoint.cpp)
rosbuild_add_executable(transformMaintenance src/transformMaintenance.cpp)
rosbuild_add_executable(registerPointCloud src/registerPointCloud.cpp)

target_link_libraries(bundleAdjust ${EIGEN3_INCLUDE_DIR})
