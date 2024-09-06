# Readme
**Disclaimer:**&nbsp; 
_I have not written the GNUplot library. This is a copy from GNUplot which is used to display the resulting graphs_

You can find the whole library [here](https://github.com/dstahlke/gnuplot-iostream/tree/master)

## Before running the program

For plotting, this code uses GNUplot. Hence, make sure it installed.

Also, for compability reasons, the program runs on a C++ 17 compiler so also please check that you have installed that version.

## Running the program

The GNUplot library uses Boost which we need to add to the program parameters when running it. Thus, you have to figure out where you have saved boost and add the path to the command below to compile the program.

> g++ -std=c++17 black_scholes_heston.cpp -o black_scholes_heston -I/opt/homebrew/opt/boost/include -L/opt/homebrew/opt/boost/lib -lboost_iostreams -lboost_system -lboost_filesystem

After you run this, you can continue as usually with any other C++ program and just run it using

> ./black_scholes_heston
