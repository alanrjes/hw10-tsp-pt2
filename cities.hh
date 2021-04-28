#pragma once

#include <iostream>
#include <vector>
#include <utility>

// Representation of an ordering of cities
class Cities {
public:
  // A pair of integral coordinates for each city
  using coord_t = std::pair<int, int>;
  // An ordering of cities. Each value represents a unique index into the current city ordering in some container.
  using permutation_t = std::vector<unsigned int>;

  // Given a permutation, return a new Cities object with the city coords indexed in the order given by the permutation.
  Cities reorder(const permutation_t& ordering) const;
  // Given a permutation in a Cities object, compute the distance through all the cities in the given order, returning to the first city.
  double total_path_distance(const permutation_t& ordering) const;
  // Generates and returns a new permutation of all the numbers from 0 to len-1 (with no skipped or repeated indices).
  permutation_t random_permutation();

  friend std::istream& operator>> (std::istream& is, Cities& cities);
  friend std::ostream& operator<< (std::ostream& os, Cities& cities);

private:
  std::vector<coord_t> map_;
};
