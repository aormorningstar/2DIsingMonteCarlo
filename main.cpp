/*
main.cpp
2D Ising Model MC Simulation
Alan Morningstar
November 2016
note: compile with c++11
*/

#include <iostream>
#include "ising.h"
#include <string>

int main(){

  int latticeLength = 8;
  int N = latticeLength*latticeLength;
  std::string stringN = std::to_string(N);

  // simulation parameters
  int intervalSteps = N*N;
  int maxSteps = (10^5)*intervalSteps;

  // temperature range
  double maxTemp = 3.5;
  double minTemp = 1.0;

  // number of simulations, set array of temperatures, include critical temp
  int numSims = 26;
  double temperature[numSims];
  for(int i = 0; i < numSims; i++) {
    temperature[i] = minTemp + i*(maxTemp-minTemp)/(numSims-1.);
  }

  // file name for writing observables data
  std::string measurementFileName = "IsingMeasurements_N=" + stringN + ".csv";

  // run simulations
  for(int i = 0; i < numSims; i++){
    std::cout << " " << std::endl;
    std::cout << "Temperature: " << temperature[i] << std::endl;
    std::cout << " " << std::endl;

    // set initial conditions
    bool coldStart = false;

    // file name for state data
    std::string stateFileName = "states_N=" + stringN + "_T=" + std::to_string(temperature[i]) + ".csv";

    // initialize simulations
    Ising IsingSim(latticeLength,temperature[i],coldStart);

    // run monte carlo steps
    for (int WolffSteps = 0; WolffSteps < maxSteps; WolffSteps++) {
      // a cluster update
      IsingSim.WolffUpdate();
      // a sweep of Metopolis updates
      for (int MetropSteps = 0; MetropSteps < N; MetropSteps++) {
        IsingSim.MetropolisUpdate();
      }

      // update measurements and write state
      if ((WolffSteps > 1000*N) && (WolffSteps % intervalSteps == 0)) {
        // measure the lattice and update running observables
        IsingSim.updateMeasurements(WolffSteps);
        /* IsingSim.updateData(); // if tracking history of observables */
        // save a snapshot of the spins
        IsingSim.writeState(stateFileName);
      }

    }

    /* IsingSim.writeData("IsingSimData.csv"); // if tracking history of obs. */

    std::cout << "Completed simulation " << i+1 << " of " << numSims << std::endl;
    std::cout << "------------------------" << std::endl;

    if (i==0) {
      // put a header in
      IsingSim.writeMeasurements(measurementFileName,true);
    }
    else{
      // no header needed
      IsingSim.writeMeasurements(measurementFileName,false);
    }
    // print final results
    IsingSim.printMeasurements();

  }

  return 0;
}
