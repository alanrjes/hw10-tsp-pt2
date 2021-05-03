#include "cities.hh"
#include "deme.hh"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <numeric>

Cities::permutation_t random_search(Cities& cities, unsigned iters) {  // my own code from part 1
  Cities::permutation_t bestroute;
  double bestdist;
  for (int i=0; i<iters; i++) {
    Cities::permutation_t nroute = cities.random_permutation();
    double ndist = cities.total_path_distance(nroute);
    if (ndist < bestdist || !i) {  // update bestroute if better ordering of route or first iteration
      bestroute = nroute;
      bestdist = ndist;
      std::cout << i << " " << bestdist << std::endl;
    }
  }
  return bestroute;
}

// ga_search uses a genetic algorithm to solve the traveling salesperson problem for a given list of cities.
// This function then creates a randomly generated population of permutations for traveling to those cities.
// The function also requires a population size and mutation rate, to indicate how aggressively the population's individuals mutate.
// The function then repeatedly evolves the population to generate increasingly better (i.e. shorter total distances) city permutations.
// The best cities permutation is returned.
Cities::permutation_t ga_search(Cities& cities, unsigned iters, unsigned pop_size, double mutation_rate) {
  auto best_dist = 1e100;
  auto best_ordering = Cities::permutation_t(cities.size());

  Deme deme(&cities, pop_size, mutation_rate);

  // Evolve the population to make it fitter and keep track of
  // the shortest distance generated
  for (long i = 1; i <= iters/pop_size; ++i) {
    deme.compute_next_generation();    // generate next generation

    // Find best individual in this population
    const auto ordering = deme.get_best()->get_ordering();
    if (is_improved(cities, ordering, best_dist, i * pop_size)) {
      best_ordering = ordering;
    }
  }
  return best_ordering;
}


int main(int argc, char* argv[]) {  // heavily reworked the version from moodle, to work with how I have Cities set up
  if (argc != 4) {
    std::cerr << "Required arguments: filename for cities, population size, and mutation rate\n";
    return -1;
  }

  char* filename = argv[1];  // if the user gives multiple file arguments, just use the first
  std::ifstream ifile(filename);  // file stream to read
  Cities cities;
  ifile >> cities;  // read file map to cities object

  auto popsize = atoi(argv[2]);
  auto mutrate = atof(argv[3]);

  constexpr unsigned itercount = 100000;
//  Cities::permutation_t bestroute = ga_search(cities, itercount, popsize, mutrate);
  Cities::permutation_t bestroute = random_search(cities, itercount);

  Cities bestcities = cities.reorder(bestroute);
  std::ofstream ofile("shortest.tsv");  // file stream to write
//  std::ofstream ofile("randomized.tsv");
  ofile << bestcities;
  return 0;
}
