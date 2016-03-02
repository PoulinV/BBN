# INSTALLATION OF THE CODE

To use the code, here are the steps you have to follow.

      mkdir build
      cd build/
      cmake .. ***WARNING FOR MAC USERS : to compile with openMP, check next instruction line***      
      make (-jN)


# FOR MAC USERS :

To compile with openMP, use gcc instead of clang. Here is an example.

    cmake -DCMAKE_CXX_COMPILER=/usr/local/bin/g++ -DCMAKE_C_COMPILER=/usr/local/bin/gcc ..




# EXPORT PATH

You can add the executable to you path in order to avoid switching directories while working. Here is the command.

    export PATH=$PATH:$(pwd)/build/src/

# MAKE FROM THE MAIN REPOSITORY

You can do the following command in order to launch 'make' directly from the BBN directory.

    cd build ; make clean ; make -j ; cd ..

# CHOOSE THE NUMBER OF CORES

To choose the number of cores with which cBBNfast is running you can use the following command before launching the executable.

    OMP_NUM_THREADS=4 ./cbbn param.ini
