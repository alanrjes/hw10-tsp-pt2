#include "cities.hh"
#include <cmath>  // std::sqrt
#include <numeric>  // std::iota
#include <algorithm>  // std::random_shuffle
#include <random>
#include <ctime>
#include <iostream>

std::istream& operator>> (std::istream& is, Cities& cities) {
  Cities::coord_t coords;
  int i = 0;
  while (is>> coords.first>> coords.second) {
    cities.map_.push_back(coords);
    i++;
  }
  return is;
}

std::ostream& operator<< (std::ostream& os, Cities& cities) {
  int len = cities.map_.size();
  for (int i=0; i<len; i++) {
    Cities::coord_t coords = cities.map_[i];
    os<< coords.first<< " "<< coords.second<< "\n";
  }
  return os;
}

double Cities::total_path_distance(const permutation_t& ordering) const {
  int len = ordering.size();
  coord_t firstcoords = map_[ordering[0]];
  coord_t prevcoords;
  double dist;  // sum of distance

  for (int i=1; i<=len; i++) {
    coord_t coords;
    if (i == len) {  // if it's the last iteration
      coords = firstcoords;  // find returning distance to the first city
    }
    else {
      coords = map_[ordering[i]];  // coords of next city in permutation, adjusted by ordering index
    }
    double subdist = std::sqrt(std::pow((coords.first - prevcoords.first), 2) + std::pow((coords.second - prevcoords.second), 2));
    dist += subdist;
    prevcoords = coords;
  }
  return dist;
}

Cities Cities::reorder(const permutation_t& ordering) const {
  int len = ordering.size();
  std::vector<coord_t> nmap(len);
  for (int i=0; i<len; i++) {
    coord_t coords = map_[ordering[i]];
    nmap[i] = coords;
  }
  Cities ncities;
  ncities.map_ = nmap;
  return ncities;
}

Cities::permutation_t Cities::random_permutation() {
  int len = map_.size();
  permutation_t ordering(len);
  std::iota(std::begin(ordering), std::end(ordering), 0);  // populate ordering vector {0, 1, 2...len}

  std::srand(unsigned(std::time(0)));  // seed randomization
  std::random_shuffle(std::begin(ordering), std::end(ordering));  // randomly shuffle vector
  return ordering;
}
