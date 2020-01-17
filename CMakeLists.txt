cmake_minimum_required (VERSION 3.0)

project (libSFS)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

set(SOURCE_FILES src/main.cpp)
add_executable (testSFS ${SOURCE_FILES})
target_link_libraries (testSFS Eigen3::Eigen)
