**HW 9: TSP Part 1**

The traveling salesperson problem part 1: randomized search.

**Description**

Uh oh! You waited too long to sell last term's textbooks back to the Reed bookstore, and now they won't take them anymore! There are several students across campus who will buy these books, but you have limited time left on campus to deliver them and make the transaction. How are you going to figure out a path for delivering all the books as quickly as possible?

**Parts**
1. *Cities operator>> ()* : Take a stream and empty *Cities* object and read from the stream into the *Cities* object map.
2. *Cities operator<< ()* : Take a stream and *Cities* object and write the city object map to the stream.
3. *Cities total_path_distance()* : Return the total distance to travel through all cities in the permutation of the *Cities* object, returning to the first city.
4. *Cities reorder()* : Take a permutation and return a new *Cities* object with the city coordinates indexed in the order of the permutation.
5. *Cities random_permutation()* : Return a "very random" new permutation of a cities object.
6. *tsp.cc main()* : Take a map and iterate over random permutations of the cities. On each iteration, compare the *total_path_distance* of each, and print out the iteration number and distance result. After a certain number of iterations, save the best permutation to a *shortest.tsv* file.
7. *speed.gif* and *shortest.gif* : Visualize results using gnuplot. Sample graphs are generated from 10 million iterations, for moderately better consistency with number of cities in *challenge.tsv* than the shorter runtime of a million iterations.

**Authors**

Alan Jessup, certified hermit. [Github](https://github.com/alwritescode).
