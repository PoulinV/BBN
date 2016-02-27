
add_executable(cbbn_test test.cpp)
target_link_libraries(cbbn_test cbbn_library)

add_test(tests cbbn_test)
