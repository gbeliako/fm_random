# fm_random
Random generation of fuzzy measures
Two functions generate random fuzzy measures (general or k-interactive) by the MiimalsPlus and topologial sort methods
followed by Markov chain.

They are based on generating linear extensions of partial orders. Each linear extension corresponds to a simplex in the polytope of fuzzy measures.
Then selecting a random point within a simplex results in a random point within the order polytope (uniformly distributed)

In addition, the functions may keep track of the linear extensions generated and then can print their distributions and also 
calculates the distance from the uniform. 

The method is based on randomisation and topological sorting as the initialisation step, followed by a few Markov chain steps.


More information is in the header file.

Gleb Beliakov, April 2020

