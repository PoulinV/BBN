project(BBNcode)
cmake_minimum_required(VERSION 2.8)
file(GLOB include "${CMAKE_CURRENT_DIR}/*.h")

install(FILES ${include} DESTINATION include/BBNcode)
