cmake_minimum_required(VERSION 3.16)
project(super_lio)


#DEBUG CONFIGURATION
# set(CMAKE_BUILD_TYPE "Debug")
# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_FLAGS_DEBUG "-O0 -Wall -g -pthread -fsanitize=address")

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -pthread -fexceptions -flto=auto")
add_definitions(-DROOT=\"${CMAKE_CURRENT_SOURCE_DIR}/\")

# Find packages
find_package(Eigen3 REQUIRED)
find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  basic
  pcl_ros
  std_msgs
  geometry_msgs
  sensor_msgs
  visualization_msgs
  nav_msgs
  message_generation
)
find_package(PythonLibs REQUIRED)

add_message_files(
  FILES
  CloudPose.msg
  CloudPose2.msg
)

generate_messages(
  DEPENDENCIES
  std_msgs
  sensor_msgs
  geometry_msgs
)

# Include directories for libraries
include_directories(
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
  include
  3rdparty
)

# Package setup for catkin
catkin_package(
  CATKIN_DEPENDS geometry_msgs nav_msgs roscpp rospy std_msgs message_runtime
  DEPENDS EIGEN3 PCL
  INCLUDE_DIRS include
)

# LIO library source files
file(GLOB_RECURSE LIO_LIB_SRC CONFIGURE_DEPENDS "src/lio/*.cpp" "src/ros/*.cpp")

file(GLOB_RECURSE LIO_LIB_HEADERS "include/lio/*.hpp" "include/lio/*.h")

# library build
add_library(lio ${LIO_LIB_SRC} ${LIO_LIB_HEADERS})

# Target specific includes
target_include_directories(lio PUBLIC
  include
  3rdparty
  ${EIGEN3_INCLUDE_DIRS}
  ${catkin_INCLUDE_DIRS}
  ${PCL_INCLUDE_DIRS}
)

# Linking libraries
target_link_libraries(lio
  ${EIGEN3_LIBRARIES}
  ${catkin_LIBRARIES}
  ${PCL_LIBRARIES}
  ${PYTHON_LIBRARIES}
  tbb 
  gflags 
  glog
)


add_executable(super_lio_node src/apps/super_lio_node.cpp)
target_link_libraries(super_lio_node
  lio
  ${catkin_LIBRARIES}
)


add_executable(relocation_node src/apps/relocation_node.cpp)
target_link_libraries(relocation_node
  lio
  ${catkin_LIBRARIES}
)

if(CATKIN_ENABLE_TESTING)
  add_executable(map_io_test tests/map_io_test.cpp)
  target_link_libraries(map_io_test lio)
  add_test(NAME map_io COMMAND map_io_test)
  add_executable(deep_logic_test tests/deep_logic_test.cpp)
  target_link_libraries(deep_logic_test lio)
  add_test(NAME deep_logic COMMAND deep_logic_test)
  add_executable(logic_review_test tests/logic_review_test.cpp)
  target_link_libraries(logic_review_test lio)
  add_test(NAME logic_review COMMAND logic_review_test)
  add_executable(motion_uncertainty_test tests/motion_uncertainty_test.cpp)
  target_link_libraries(motion_uncertainty_test lio)
  add_test(NAME motion_uncertainty COMMAND motion_uncertainty_test)
  add_executable(input_review_test tests/input_review_test.cpp)
  target_link_libraries(input_review_test lio)
  add_test(NAME input_review COMMAND input_review_test)
  add_executable(spectral_reliability_test tests/spectral_reliability_test.cpp src/lio/geometry.cpp)
  add_test(NAME spectral_reliability COMMAND spectral_reliability_test)
  add_executable(lio_numerics_test tests/lio_numerics_test.cpp)
  target_link_libraries(lio_numerics_test lio)
  add_test(NAME lio_numerics COMMAND lio_numerics_test)
  add_executable(geometry_test tests/geometry_test.cpp src/lio/geometry.cpp)
  add_test(NAME geometry_core COMMAND geometry_test)
  add_executable(geometry_filter_test tests/geometry_filter_test.cpp)
  target_link_libraries(geometry_filter_test lio)
  add_test(NAME geometry_filter COMMAND geometry_filter_test)
  add_executable(geometry_map_test tests/geometry_map_test.cpp)
  target_link_libraries(geometry_map_test lio)
  add_test(NAME geometry_map COMMAND geometry_map_test)
  add_executable(geometry_sampling_test tests/geometry_sampling_test.cpp)
  target_link_libraries(geometry_sampling_test lio)
  add_test(NAME geometry_sampling COMMAND geometry_sampling_test)
  add_executable(geometry_fallback_test tests/geometry_fallback_test.cpp)
  target_link_libraries(geometry_fallback_test lio)
  add_test(NAME geometry_fallback COMMAND geometry_fallback_test)
endif()
