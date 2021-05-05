#include "cities.hh"
#include "deme.hh"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <numeric>

// Randomized search, my own code from part 1
Cities::permutation_t random_search(Cities& cities, unsigned iters) {
  Cities::permutation_t bestroute;
  double bestdist;
  for (int i=0; i<iters; i++) {
    const Cities::permutation_t& nroute = cities.random_permutation();  // generate random permutation
    double ndist = cities.total_path_distance(nroute);
    if (ndist < bestdist || !i) {  // update bestroute if better ordering of route or first iteration
      bestroute = nroute;
      bestdist = ndist;
      std::cout << i << " " << bestdist << std::endl;
    }
  }
  return bestroute;
}

// Gentic algorithm search!
// Create a random population of permutations (Deme), with a given population size and mutation rate, then repeatedly evolves for higher fitness (ie lower distance).
Cities::permutation_t ga_search(Cities& cities, unsigned iters, unsigned population, double mutation) {
  Deme deme(&cities, population, mutation);
  Cities::permutation_t bestroute;
  double bestdist;

//  double avg;  // FOR TESTING

  for (int i=0; i<(iters/population); i++) {
    deme.compute_next_generation();  // generate new generation
    Cities::permutation_t nroute = deme.get_best() -> get_ordering();  // get best out of the generation
    double ndist = cities.total_path_distance(nroute);

    if (ndist < bestdist || !i) {  // update bestroute if better ordering of route or first iteration

//      std::cout << "Average improvement: " << avg << std::endl;  // FOR TESTING

      bestroute = nroute;
      bestdist = ndist;
      std::cout << i << " " << bestdist << std::endl;
    }
//    avg = (avg + bestdist - ndist)/2;  // FOR TESTING
  }
  return bestroute;
}

// heavily reworked from the version from moodle, to work with how I have Cities set up
int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cerr << "Required arguments: filename for cities, population size, and mutation rate\n";
    return -1;
  }

  srand(time(NULL));

  char* filename = argv[1];  // if the user gives multiple file arguments, just use the first
  std::ifstream ifile(filename);  // file stream to read
  Cities cities;
  ifile >> cities;  // read file map to cities object

  unsigned popsize = atoi(argv[2]);
  double mutrate = atof(argv[3]);

  constexpr unsigned itercount = 100000;
  Cities::permutation_t bestroute = ga_search(cities, itercount, popsize, mutrate);
//  Cities::permutation_t bestroute = random_search(cities, itercount);

  Cities bestcities = cities.reorder(bestroute);
  std::ofstream ofile("shortest.tsv");  // file stream to write
//  std::ofstream ofile("randomized.tsv");
  ofile << bestcities;
  return 0;
}
