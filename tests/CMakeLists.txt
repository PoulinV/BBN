
add_executable(cbbn_test test.cpp)
target_link_libraries(cbbn_test cbbn_library)
add_test(test_uno cbbn_test)

add_executable(test_ini ini.cpp)
target_link_libraries(test_ini cbbn_library)
add_test(test_ini test_ini)
