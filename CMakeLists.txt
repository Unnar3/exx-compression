cmake_minimum_required(VERSION 2.8.3)
project(exx_compression)

set(CMAKE_CXX_FLAGS "-O4 -g -fPIC -std=c++0x -fpermissive")

find_package(catkin REQUIRED COMPONENTS
	roscpp
	pcl_conversions
  	pcl_ros
  	sensor_msgs
  	std_msgs
  	message_generation
)

add_message_files(FILES
	planes.msg
)

generate_messages(DEPENDENCIES
	std_msgs
	sensor_msgs
)

catkin_package(
	DEPENDS
	CATKIN_DEPENDS roscpp message_runtime
	INCLUDE_DIRS include
	LIBRARIES compression features
)
include_directories(include ${catkin_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
add_library(compression src/compression.cpp)
add_library(features src/plane_features.cpp)
#add_executable(plane_compression src/plane_compression.cpp)
target_link_libraries(features ${catkin_LIBRARIES})
target_link_libraries(compression ${catkin_LIBRARIES} libvtkCommon.so libvtkFiltering.so libvtkGraphics.so)