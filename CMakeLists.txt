project(BBNcode)
cmake_minimum_required(VERSION 2.8)
add_executable(cbbn ../source/BBN_constraints.cpp	../source/EM_cascade.cpp ../source/Injected_spectrum.cpp ../source/tools.cpp	../source/test_functions.cpp ../source/main.cpp)
add_executable(cbbn_test ../source/BBN_constraints.cpp	../source/EM_cascade.cpp ../source/Injected_spectrum.cpp ../source/tools.cpp ../source/test_functions.cpp	../source/test.cpp)
