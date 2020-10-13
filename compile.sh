g++ -std=c++11 -O3 -c -o ising.o ising.cpp
g++ -std=c++11 -O3 -c -o main.o main.cpp
g++ -std=c++11 -O3 -o isingsim ising.o main.o
