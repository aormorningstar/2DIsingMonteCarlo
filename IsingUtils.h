/*
utils.h
source code for utilities used in Ising
Alan Morningstar
January 2017
*/

#ifndef __UTILS_H_INCLUDED__
#define __UTILS_H_INCLUDED__

#include <random>

using namespace std;

namespace IsingUtils {

  // SET UP RANDOM NUMBER GENERATOR
  // random_device rd;
  // mt19937 e2(rd());
  // uniform_real_distribution<double> dist(0, 1);
  // call dist(e2) for random float between 0 and 1

  double unifRand() {
      return rand() / double(RAND_MAX);
  };

}

#endif
