# FOR MAC USERS :

To compile with openMP, use gcc instead of clang. Here is an example.

    cmake -DCMAKE_CXX_COMPILER=/usr/local/bin/g++ -DCMAKE_C_COMPILER=/usr/local/bin/gcc ..




# EXPORT PATH

You can add the executable to you path in order to avoid switching diretories while working. Here is the command.

    export PATH=$PATH:$(pwd)/build/src/

# CHOOSE THE NUMBER OF CORES

To choose the number of cores with which cBBNfast is running you can use the following command before launching the executable.

    OMP_NUM_THREADS=4 ./cbbn param.ini
