project(bbn)

cmake_minimum_required(VERSION 2.8)

enable_testing()
set(CMAKE_MACOSX_RPATH ON)

include_directories(include)

add_subdirectory(include)
add_subdirectory(lib)
add_subdirectory(src)
add_subdirectory(tests)
