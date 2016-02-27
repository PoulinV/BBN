
file(GLOB cbbn_sources *.cpp)
add_executable(cbbn_binary ${cbbn_sources})
set_target_properties(cbbn_binary PROPERTIES OUTPUT_NAME cbbn)
target_link_libraries(cbbn_binary cbbn_library)
