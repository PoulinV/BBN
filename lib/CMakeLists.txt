
file(GLOB libcbbn_sources *.cpp)
add_library(cbbn_library SHARED ${libcbbn_sources})
set_target_properties(cbbn_library PROPERTIES OUTPUT_NAME cbbn)
