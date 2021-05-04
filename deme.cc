#include <vector>
#include "chromosome.hh"
#include "deme.hh"

// Generate a Deme of the specified size with all-random chromosomes.
// Also receives a mutation rate in the range [0-1].
Deme::Deme(Cities* cities_ptr, unsigned pop_size, double mut_rate) {
  mut_rate_ = mut_rate;
  for (int i=0; i<pop_size; i++) {
    Chromosome* chromosome = new Chromosome(cities_ptr);
    pop_.push_back(chromosome);
  }
}

// Clean up as necessary
Deme::~Deme() {
}

// Return a copy of the chromosome with the highest fitness.
const Chromosome* Deme::get_best() const {
  double bestfit = 0;
  Chromosome* best;
  for (int i=0; i<pop_.size(); i++) {
    Chromosome* thischromo = pop_[i];
    double fitness = thischromo -> get_fitness();
    if (fitness>bestfit) {
      bestfit = fitness;
      best = thischromo;
    }
  }
  return best;
}

// Randomly select a chromosome in the population based on fitness and return a pointer to that chromosome.
Chromosome* Deme::select_parent() {
  double fitsum = 0;
  std::vector<double> fitscores;  // store so that we only have to run get_fitness on each chromosome once
  for (int i=0; i<pop_.size(); i++) {
    Chromosome* thischromo = pop_[i];
    double fitness = thischromo -> get_fitness();
    fitsum += fitness;
    fitscores.push_back(fitness);
  }
  // get randomized double for "roulette wheel" between 0 and the sum of fitness scores
  std::uniform_real_distribution<double> unif(0, fitsum);
  std::default_random_engine re;
  double target = unif(re);

  int i = -1;  // -1 because i increments after the loop check
  fitsum = 0;
  while (fitsum < target) {
    i++;
    fitsum += fitscores[i];
  }
  return pop_[i];
}

// Evolve a single generation of new chromosomes, as follows: We select pop_size/2 pairs of chromosomes (using the select() method below).
// Each chromosome in the pair can be randomly selected for mutation, with probability mut_rate, in which case it calls the chromosome mutate() method.
// Then, the pair is recombined once (using the recombine() method) to generate a new pair of chromosomes, which are stored in the Deme.
// After we've generated pop_size new chromosomes, we delete all the old ones.
void Deme::compute_next_generation() {
  std::vector<Chromosome*> newpop;
  for (int i=0; i<(pop_.size()/2); i++) {
    Chromosome* p1 = this -> mutate_parent();
    Chromosome* p2 = this -> mutate_parent();
    std::pair<Chromosome*, Chromosome*> children = p1 -> recombine(p2);
    newpop.push_back(children.first);
    newpop.push_back(children.second);
  }
  pop_ = newpop;
}

Chromosome* Deme::mutate_parent() {
  Chromosome* parent = this -> select_parent();
  int mutroll = rand() % 100;
  if (mutroll < mut_rate_*100) {
    parent -> mutate();
  }
  return parent;
}
