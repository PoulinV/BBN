
# BBN

## Introduction

...TO BE CONTINUED...

## Dependencies

You will need the following programs on your computer to build the code and generate the documentation :

  * [CMake](http://cmake.org/)
  * C++ Compiler ([clang](http://clang.llvm.org/) / [g++](http://gcc.gnu.org/))
  * [Doxygen](http://www.stack.nl/~dimitri/doxygen/)

## Build

### Basic instruction

    mkdir build
    cd build
    cmake ..
    make

### OpenMP

On MacOSx, the current version of clang not support OpenMP. If you wish to use BBN with multicore support, you have to use gcc instead of clang.

    cmake -DCMAKE_CXX_COMPILER=/usr/local/bin/g++ -DCMAKE_C_COMPILER=/usr/local/bin/gcc ..

### Documentation

To build the documentation, run the following command :

    cd build
    make doc

The documentation will be generated in "build/doc".

## Run

### Add executable in your path

You can add the executable to your path in order to avoid switching directories while working. Here is the command.

    export PATH=$PATH:$(pwd)/build/src/

### Choose the number of core

To choose the number of cores with which cBBNfast is running you can use the following command before launching the executable.

    OMP_NUM_THREADS=4 ./cbbn param.ini
