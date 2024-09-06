# Readme
**Disclaimer:**&nbsp; 
_I have not written the GNUplot library. This is a copy from GNUplot which is used to display the resulting graphs_

You can find the whole library [here](https://github.com/dstahlke/gnuplot-iostream/tree/master)


# About the program
This Code is the addition for a thesis I wrote during my university course Trading and Sales. The aim of the thesis was to compare both of those models on runtime complexity and how far the results for option pricing drift away to see if there is a chance to use a faster model for arbitraging options or to find out if they are beeing misspriced at the current time. The programm calculates the implied volatility and uses it for pricing a stock that can be entered and also taking the timevalue into consideration. The results are plotted with GNUplot and the sensitivity on reactions to volatility and time of those models are beeing compared. The models I am comparing in this example are the Black-Scholes model and the Merton Jump Diffusion Model.

## Before running the program

For plotting, this code uses GNUplot. Hence, make sure it installed.

Also, for compability reasons, the program runs on a C++ 17 compiler so also please check that you have installed that version. Theoretically C++ 14 should work as well but hasn't been tested.

## Running the program

The GNUplot library uses Boost which we need to add to the program parameters when running it. Thus, you have to figure out where you have saved boost and add the path to the command below to compile the program.

> g++ -std=c++17 black_scholes_heston.cpp -o black_scholes_heston -I/opt/homebrew/opt/boost/include -L/opt/homebrew/opt/boost/lib -lboost_iostreams -lboost_system -lboost_filesystem

After you run this, you can continue as usually with any other C++ program and just run it using

> ./black_scholes_heston
