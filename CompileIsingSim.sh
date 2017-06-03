g++ -std=c++11 -O3 -c -o Ising.o Ising.cpp
g++ -std=c++11 -O3 -c -o main.o main.cpp
g++ -std=c++11 -O3 -o IsingSim Ising.o main.o
