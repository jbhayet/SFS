cmake_minimum_required (VERSION 3.0)

project (libSFS)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(SOURCE_FILES src/main.cpp)
add_executable (testSFS ${SOURCE_FILES})
target_link_libraries (testSFS Eigen3::Eigen)
