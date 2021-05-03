**TSP Project**

The traveling salesperson problem, with randomized search and genetic algorithm to find a good path.

**Description**

Uh oh! You waited too long to sell last term's textbooks back to the Reed bookstore, and now they won't take them anymore! There are several students across campus who will buy these books, but you have limited time left on campus to deliver them and make the transaction. How are you going to figure out a path for delivering all the books as quickly as possible?

**HW 9 : Part 1**
1. *Cities operator>> ()* : Take a stream and empty *Cities* object and read from the stream into the *Cities* object map.
2. *Cities operator<< ()* : Take a stream and *Cities* object and write the city object map to the stream.
3. *Cities total_path_distance()* : Return the total distance to travel through all cities in the permutation of the *Cities* object, returning to the first city.
4. *Cities reorder()* : Take a permutation and return a new *Cities* object with the city coordinates indexed in the order of the permutation.
5. *Cities random_permutation()* : Return a "very random" new permutation of a cities object.
6. *tsp.cc random_search()* : Take a map and iterate over random permutations of the cities. On each iteration, compare the *total_path_distance* of each, and print out the iteration number and distance result. After a certain number of iterations, save the best permutation to a *shortest.tsv* file.
7. *speed.gif* and *shortest.gif* : Visualize results using gnuplot. Sample graphs are generated from 10 million iterations, for moderately better consistency with number of cities in *challenge.tsv* than the shorter runtime of a million iterations.

**HW 10 : Part 2**

*Chromosome* class
1. *is_valid()* : Check object's chromosome for duplicate indices, or any indices out of range. If either is found, return false.
2. *get_fitness()* : Return a value representing the object's "fitness", higher for better and lower if worse.
3. Constructor/destructor : No changes.
4. *mutate()* : Create variation by swapping two values in the object's chromosome permutation.
5. *recombine()* : Take an other parent chromosome, and use *create_crossover_child()* to combine its object parent and the other parent via ordered crossover.

*Deme* class
1. Constructor : Take a *Cities* object, population size, and mutation rate, and generate a random population of *Chromosome* objects from the *Cities* object.
2. *get_best()* : Return the chromosome with the highest fitness in the population.
3. *select_parent()* : Return a chromosome to become a parent, selected from the population via a fitness proportionate randomized selection.
4. *compute_next_generation()* : Pick two parents at random, potentially mutate based on population object's mutation rate, and recombine. Add both children to a new generation population, repeat for half the population size number of iterations, and when done, replace the object's population with the new generation population.
3. *tsp.cc ga_search()* : Create a random population of permutations (a Deme object), with a given population size and mutation rate. Then, repeatedly evolve for higher fitness (ie lower distance), comparing the distances of the best permutation in each generation of the population on each iteration.

**Authors**

Alan Jessup, certified hermit. [Github](https://github.com/alwritescode).
