#include <algorithm>
#include <cassert>
#include <random>
#include <cmath>

#include "chromosome.hh"

// Constructor: generate a completely random permutation from a list of cities
Chromosome::Chromosome(Cities* cities_ptr) {
  cities_ptr_ = cities_ptr;
  order_ = cities_ptr -> random_permutation();
  assert(is_valid());
}

// Destructor: clean up as necessary
Chromosome::~Chromosome() {
  assert(is_valid());
}

// A chromosome is valid if it has no repeated values in its permutation, as well as no indices above the range (length) of the chromosome.
// We implement this check with a sort, which is a bit inefficient, but simple
bool Chromosome::is_valid() const {
  Cities::permutation_t sorted = order_;
  std::sort(sorted.begin(), sorted.end());
  int len = sorted.size();

  for (int i=1; i<len; i++) {
    if (sorted[i-1] == sorted[i]) {  // check for duplicates
      return false;
    }
  }

  if (sorted[len - 1] >= len) {  // index out of range of chromosome length
    return false;
  }

  else {
    return true;
  }
}

// Return a positive fitness value, with higher numbers representing fitter solutions (shorter total-city traversal path).
double Chromosome::get_fitness() const {
  double dist = cities_ptr_ -> total_path_distance(order_);
  return std::pow(1/dist, 2);
//  return 1/dist;
}

// Perform a single mutation on this chromosome by swapping two values
void Chromosome::mutate() {
  int len = order_.size() - 1;
  int i1 = std::rand() %len;
  int i2 = std::rand() %len;
  while (i1 == i2) {
    i2 = std::rand() %len;
  }
  int val1 = order_[i1];
  order_[i1] = order_[i2];
  order_[i2] = val1;
  assert(is_valid());
}

// Return a pair of offsprings by recombining with another chromosome. Note: this method allocates memory for the new offsprings
std::pair<Chromosome*, Chromosome*> Chromosome::recombine(const Chromosome* other) {
  std::pair<Chromosome*, Chromosome*> children;
  for (int i=0; i<2; i++) {
    assert(is_valid());
    assert(other->is_valid());

    int len = order_.size();
    int lower = std::rand() %len;
    int upper = std::rand() %len;
    while (lower==upper) {  // generate different lower and upper bounds of a range
      upper = std::rand() %len;
    }
    if (lower>upper) {  // if lower and upper bounds are swapped, switch them
      int lowholder = lower;
      lower = upper;
      upper = lowholder;
    }

    Chromosome* child = create_crossover_child(this, other, lower, upper);
    if (!i) {
      children.first = child;
    }
    else {
      children.second = child;
    }
  }
  return children;
}

// Find whether a certain value appears in a given range of the chromosome.
// Returns true if value is within the specified the range specified [begin, end) and false otherwise.
bool Chromosome::is_in_range(unsigned value, unsigned begin, unsigned end) const {
  for (int i=begin; i<end; i++) {
    if (value == order_[i]) {
      return true;
    }
  }
  return false;
}

// For an ordered set of parents, return a child using the ordered crossover.
// The child will have the same values as p1 in the range [b,e), and all the other values in the same order as in p2.
Chromosome* Chromosome::create_crossover_child(const Chromosome* p1, const Chromosome* p2, unsigned b, unsigned e) const {
  Chromosome* child = p1->clone();

  // We iterate over both parents separately, copying from parent1 if the
  // value is within [b,e) and from parent2 otherwise
  unsigned i = 0, j = 0;

  for ( ; i < p1->order_.size() && j < p2->order_.size(); ++i) {
    if (i >= b and i < e) {
      child->order_[i] = p1->order_[i];
    }
    else { // Increment j as long as its value is in the [b,e) range of p1
      while (p1->is_in_range(p2->order_[j], b, e)) {
        ++j;
      }
      assert(j < p2->order_.size());
      child->order_[i] = p2->order_[j];
      j++;
    }
  }

  assert(child->is_valid());
  return child;
}
