project(bbn)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

include(CheckCPP11Flags)

if(APPLE)
  set(CMAKE_MACOSX_RPATH ON)
endif(APPLE)

include_directories(include)

enable_testing()

add_subdirectory(include)
add_subdirectory(lib)
add_subdirectory(src)
add_subdirectory(tests)

