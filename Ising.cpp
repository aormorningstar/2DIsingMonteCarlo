/*
ising.cpp
source for Ising class
2D Ising Model MC Simulation
Alan Morningstar
November 2016
*/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <string>
#include "ising.h"
#include "isingutils.h"

using namespace std;
using namespace IsingUtils;

// constructor
Ising::Ising(int length,double temperature,bool coldStart = true) {
  // ** lattice properties
  d = 2;
  J = 1.;
  L = length;
  N = length*length;
  doubleN = double(length*length);
  T = temperature;
  beta = 1./temperature;

  // ** running sums and counters for observables
  energyPerSpinRunningSum = 0;
  energyPerSpinSquaredRunningSum = 0;
  magnetizationRunningSum = 0;
  magnetizationSquaredRunningSum = 0;
  numSamples = 0;
  mcTime = 0;
  netSpin = 0;

  // ** Wolff algorithm variables
  // probability of adding a spin to the cluster
  addProb = 1.-exp(-2.*J/temperature);
  clusterSize = 0;

  // allocate spin and cluster array
  s = new int* [length];
  cluster = new bool* [length];
  for (int i=0;i<L;i++) {
    s[i] = new int [length];
    cluster[i] = new bool [length];
  }

  // initial conditions, empty cluster
  for (int i=0;i<L;i++) {
    for (int j=0;j<L;j++) {
      // cold start
      if (coldStart==true) {
        s[i][j] = 1;
      }
      // hot start
      else{
        if (unifRand()<.5) {
          s[i][j] = 1;
        }
        else {
          s[i][j] = -1;
        }
      }
      cluster[i][j] = false;
      netSpin += s[i][j];
    }
  }
}

// destructor
Ising::~Ising(void) {
  delete s;
  delete cluster;
  // if file is open, close it
  if(statesFile.is_open()){
    statesFile.close();
  }
}

// attempt to add a spin to the cluster
// component of growCluster(i,j)
void Ising::attemptAdd(int i,int j,int spin) {
  // if neighbor is aligned with current spin, add spin to cluster with probability addProb
  if (s[i][j] == spin && unifRand() < addProb) {
    // add spin to cluster with probability addProb
    growCluster (i,j);
  }
}

// grow the cluster
// component of WolffUpdate()
void Ising::growCluster(int i,int j) {
  // add spin to cluster and flip it, keep track of up spins and cluster size
  int spin = s[i][j];
  cluster[i][j] = true;
  clusterSize += 1;
  netSpin += -2*s[i][j];
  s[i][j] = -s[i][j];

  // locations of neighbors, impose periodic boundaries
  int iPrev = i == 0 ? L-1 : i-1;
  int iNext = i == L-1? 0 : i+1;
  int jPrev = j == 0 ? L-1 : j-1;
  int jNext = j == L-1? 0 : j+1;

  // try to add neighbors to the cluster
  if (!cluster[iPrev][j]){
    attemptAdd(iPrev,j,spin);
  }
  if (!cluster[iNext][j]){
    attemptAdd(iNext,j,spin);
  }
  if (!cluster[i][jPrev]){
    attemptAdd(i,jPrev,spin);
  }
  if (!cluster[i][jNext]){
    attemptAdd(i,jNext,spin);
  }
}

// perform one Wolff update step
void Ising::WolffUpdate() {
  // no spins in the cluster yet
  for (int i=0;i<L;i++) {
    for (int j=0;j<L;j++) {
    cluster[i][j] = false;
    }
  }
  clusterSize = 0;

  // choose a random seed spin from the lattice
  int iSeed = rand() % L;
  int jSeed = rand() % L;

  // grow the cluster
  growCluster(iSeed,jSeed);
}

// perform one Metropolis update step
void Ising::MetropolisUpdate() {
  // choose a random seed spin from the lattice
  int iSeed = rand() % L;
  int jSeed = rand() % L;
  int seedSpin = s[iSeed][jSeed];

  // locations of neighbors, impose periodic boundaries
  int iPrev = iSeed == 0 ? L-1 : iSeed-1;
  int iNext = iSeed == L-1? 0 : iSeed+1;
  int jPrev = jSeed == 0 ? L-1 : jSeed-1;
  int jNext = jSeed == L-1? 0 : jSeed+1;
  int neighbors[4][2] = {{iPrev,jSeed},{iNext,jSeed},{iSeed,jPrev},{iSeed,jNext}};

  // calculate energy change due to flip
  double energyDiff = 0.;
  for (int nbr = 0; nbr < 4; nbr++) {
    int iLoc = neighbors[nbr][0];
    int jLoc = neighbors[nbr][1];
    if (seedSpin == s[iLoc][jLoc]) {
      energyDiff += 2.*J;
    }
    else {
      energyDiff -= 2.*J;
    }
  }

  // attempt to flip spin
  if (energyDiff <= 0.) {
    s[iSeed][jSeed] = -s[iSeed][jSeed];
    netSpin += 2*s[iSeed][jSeed];
  }
  else if (unifRand() < exp(-beta * energyDiff)) {
    s[iSeed][jSeed] = -s[iSeed][jSeed];
    netSpin += 2*s[iSeed][jSeed];
  }
  else {
    // don't flip the spin
  }
}

// print a depiction of the current lattice state
void Ising::printLattice() {
  for (int i=0;i<L;i++) {
    for (int j=0;j<L;j++) {
    if (s[i][j]==1) {cout << "o ";}
    else {cout << "x ";}
    }
    cout << endl;
  }
}

// print current observable values
void Ising::printMeasurements() {
  cout << "monte carlo time: " << mcTime << endl;
  cout << "magnetization: " << magnetization << endl;
  cout << "energy per spin: " << energyPerSpin << endl;
  cout << "susceptability: " << susceptability << endl;
  cout << "specific heat: " << specificHeat << endl;
}

// update measured values
void Ising::updateMeasurements(int steps) {
  // update energy per spin and magnetization
  int sum = 0;
  int iPrev,iNext,jPrev,jNext;
  for (int i=0;i<L;i++) {
    for (int j=0;j<L;j++) {
    int iPrev = i == 0 ? L-1 : i-1;
    int iNext = i == L-1? 0 : i+1;
    int jPrev = j == 0 ? L-1 : j-1;
    int jNext = j == L-1? 0 : j+1;
    sum += s[i][j]*(s[iNext][j] + s[iPrev][j] + s[i][jPrev] + s[i][jNext]);
    }
  }

  // current samples
  double energyPerSpinSample = -0.5*J*double(sum)/doubleN;
  double magnetizationSample = abs(double(netSpin))/doubleN;

  // update running sums
  energyPerSpinRunningSum += energyPerSpinSample;
  energyPerSpinSquaredRunningSum += energyPerSpinSample * energyPerSpinSample;
  magnetizationRunningSum += magnetizationSample;
  magnetizationSquaredRunningSum += magnetizationSample * magnetizationSample;
  numSamples += 1;

  // update observables
  energyPerSpin = energyPerSpinRunningSum / double(numSamples);
  magnetization = magnetizationRunningSum / double(numSamples);
  specificHeat = doubleN*(energyPerSpinSquaredRunningSum/double(numSamples) - energyPerSpin*energyPerSpin)/(T*T);
  susceptability = doubleN*(magnetizationSquaredRunningSum/double(numSamples) - magnetization*magnetization)/T;
  mcTime = double(steps)/doubleN;
}

// append current measurements to data vectors if keeping
// a history of how the measurements evolve in the simulation
void Ising::updateData(){
  // push observables into data holders
  energyPerSpinData.push_back(energyPerSpin);
  magnetizationData.push_back(magnetization);
  specificHeatData.push_back(specificHeat);
  susceptabilityData.push_back(susceptability);
  mcTimeData.push_back(mcTime);
}

// write data vectors to file
void Ising::writeData(string csvFileName) {
  // create output file for data
  ofstream outputFile(csvFileName);

  // header
  outputFile << "mcTime" << "," << "magnetization" << "," << "energyPerSpin" << "," << "susceptability" << "," << "specificHeat" << endl;

  // loop over and write data
  int dataLength = mcTimeData.size();
  for (int i=0;i<dataLength;i++) {
    outputFile << mcTimeData[i] << "," << magnetizationData[i] << "," << energyPerSpinData[i] << "," << susceptabilityData[i] << "," << specificHeatData[i] << endl;
  }

  // close file
  outputFile.close();
}

// write measured quantities to file
void Ising::writeMeasurements(string csvFileName,bool header = false) {
  // open output file for measurements
  ofstream outputFile(csvFileName,ios::app);

  // add header and measurements
  if (header) {
    outputFile << "T,energyPerSpin,magnetization,specificHeat,susceptability" << endl;
  }
  outputFile << T << "," << energyPerSpin << "," << magnetization << "," << specificHeat << "," << susceptability << endl;

  // close file
  outputFile.close();
}

// write state to file as a single csv line, let the spin variables take values {0,1}
void Ising::writeState(string fileName) {
  if (!statesFile.is_open()) {
    statesFile.open(fileName);
  }

  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      if (i==0 && j==0) {
        statesFile << (s[i][j]+1)/2;
      }
      else{
        statesFile << "," << (s[i][j]+1)/2;
      }
    }
  }
  statesFile << endl;
}
