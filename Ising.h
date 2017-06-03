/*
Ising.h
Ising class header file
2D Ising Model MC Simulation
Alan Morningstar
November 2016
*/

#ifndef __ISING_H_INCLUDED__
#define __ISING_H_INCLUDED__

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

class Ising {
  private:
    bool coldStart;

    // square lattice details
    int d; // dimension
    int L; // length of the lattice in each dimension
    double J; // ferromagnetic spin-spin coupling
    int N; // number of spins in total
    double doubleN; // double version of N

    // lattice spins
    int **s; // array of spin variables

    // external parameters
    double T; // temperature
    double beta; // inverse temperature

    // measured quantities
    double magnetization;
    double susceptability;
    double energyPerSpin;
    double specificHeat;
    double mcTime;

    // lists for plotting
    std::vector<double> magnetizationData;
    std::vector<double> energyPerSpinData;
    std::vector<double> susceptabilityData;
    std::vector<double> specificHeatData;
    std::vector<double> mcTimeData;

    // running sums for observables
    double energyPerSpinRunningSum;
    double energyPerSpinSquaredRunningSum;
    double magnetizationRunningSum;
    double magnetizationSquaredRunningSum;
    int numSamples;

    // Wolff algorithm variables
    double addProb; // probability of adding a neighbor to the cluster
    bool **cluster; // true if the spin is already in the cluster
    int clusterSize;
    int netSpin; // keep track of total spin

    // output file for states
    std::ofstream statesFile;

  public:

    // functions
    Ising(int,double,bool);
    ~Ising(void);

    void growCluster(int,int);
    void attemptAdd(int,int,int);
    void WolffUpdate(void);
    void MetropolisUpdate(void);
    void updateMeasurements(int);
    void updateData(void);
    void printMeasurements(void);
    void printLattice(void);
    void writeData(std::string);
    void writeMeasurements(std::string,bool);
    void writeState(std::string);
};

#endif
