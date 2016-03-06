project(bbn)

cmake_minimum_required(VERSION 2.8)

set(BBN_VERSION "0.1.0-devel") # cannot use project before required cmake 3.0

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

include(CheckCPP11Flags)

include(FindOpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_FOUND)


if(APPLE)
  set(CMAKE_MACOSX_RPATH ON)
endif(APPLE)

include_directories(include)

enable_testing()

add_subdirectory(include)
add_subdirectory(lib)
add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(doc)
