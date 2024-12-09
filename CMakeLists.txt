cmake_minimum_required(VERSION 3.5)
project (libSFS)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(SOURCE_FILES src/generateSFS.cpp)
add_executable (generateSFS ${SOURCE_FILES})
target_link_libraries (generateSFS Eigen3::Eigen)

set(SOURCE_FILES src/singleCountSFS.cpp)
add_executable (singleCountSFS ${SOURCE_FILES})
target_link_libraries (singleCountSFS Eigen3::Eigen)
