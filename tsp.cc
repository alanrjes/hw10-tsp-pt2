#include "cities.hh"
#include <fstream>
#include <cmath>  // std::pow

int main(int argc, char* argv[]) {
  char* filename = argv[1];  // if the user gives multiple file arguments, just use the first
  std::ifstream ifile(filename);  // file stream to read
  Cities cities;
  ifile >> cities;  // read file map to cities object

  int itercount = std::pow(10, 7);
  Cities::permutation_t bestroute;
  double bestdist;
  for (int i=0; i<itercount; i++) {
    Cities::permutation_t nroute = cities.random_permutation();
    double ndist = cities.total_path_distance(nroute);
    if (ndist < bestdist || !i) {  // update bestroute if better ordering of route or first iteration
      bestroute = nroute;
      bestdist = ndist;
      std::cout << i << " " << bestdist << std::endl;
    }
  }
  Cities bestcities = cities.reorder(bestroute);
  std::ofstream ofile("shortest.tsv");  // file stream to write
  ofile << bestcities;  // write new ordered cities map to shortest file
  return 0;
}
