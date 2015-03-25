cmake_minimum_required(VERSION 2.8.3)
project(exx_compression)

find_package(catkin REQUIRED COMPONENTS
	roscpp
	pcl_conversions
  	pcl_ros
)

#add_message_files(FILES
#	distPose.msg
#)

#generate_messages(DEPENDENCIES
#	std_msgs
#)

catkin_package(
	DEPENDS
	CATKIN_DEPENDS roscpp message_runtime
	INCLUDE_DIRS include
	LIBRARIES compression
)
include_directories(include ${catkin_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
add_library(compression src/compression.cpp)
#add_executable(plane_compression src/plane_compression.cpp)
target_link_libraries(compression ${catkin_LIBRARIES} libvtkCommon.so libvtkFiltering.so libvtkGraphics.so)