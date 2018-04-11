#project cmakelist
cmake_minimum_required(VERSION 2.8.3)

project(VINS_SW)
# compile and build options
# settle gcc old version newversion conflict
ADD_DEFINITIONS(-fpermissive)
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -w  -O3 -march=native ")
set(CMAKE_CXX_FLAGS "-std=c++11")
#-DEIGEN_USE_MKL_ALL")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -w -g -march=native")

# solve libraries dependency
find_package(Eigen3 REQUIRED)
# what exactly does find ceres do?
# where to find this find_ceres.cmake or what?
find_package(Ceres REQUIRED)
find_package(OpenCV REQUIRED)
include_directories(
	${OpenCV_INCLUDE_DIRS}
	${CERES_INCLUDE_DIRS}
	${EIGEN3_INCLUDE_DIR}
)
message(STATUS "EIGEN3 include dir: " ${EIGEN3_INCLUDE_DIR}
)
MESSAGE(STATUS "OPENCV INCLUDE DIR: " ${OpenCV_INCLUDE_DIRS})
MESSAGE(STATUS "CERES INCLUDE DIR: " ${CERES_INCLUDE_DIRS}
)
add_subdirectory(./src/)