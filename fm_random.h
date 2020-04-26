/*

These two functions generate random fuzzy measures (general or k-interactive) by the MiimalsPlus and topologial sort methods
followed by Markov chain.

They are based on generating linear extensions of partial orders. Each linear extension corresponds to a simplex in the polytope of fuzzy measures.
Then selecting a random point within a simplex results in a random point within the order polytope (uniformly distributed)

In addition, the functions may keep track of the linear extensions generated and then can print their distributions and also 
calculates the distance from the uniform. 

The method is based on randomisation and topological sorting as the initialisation step, followed by a few Markov chain steps.


Gleb Beliakov, April 2020


The functions take parameters:  (num m  kint Markovsteps option  K   vv)

num is the number of linear extensions to generate
m is the number of criteria so that 2^m is the size of the fuzzy measure (or size of the binary lattice Bn)

kint - (<= m) is the number k in k-interactive fuzzy measures. For example kint=1 means the measure is (almost) additive 
(only singletons are free parameters, and there are m! linear extensions).

markovsteps - number of Markov chain steps to make. Start from 100 then 1000, 5000, ... the more is steps the beter quality 
is the distribution but the more time the algorithm works.

option - not used for minimalsplus, but for topological sort , if 1 (but it is overwritten to 0 if the number of entries to sort * num is too large)
   then the Markov chain is also followed by a selection process which looks at the history of generated linear extensions (stored in a map with the help of hash keys),
   and choosing the orders which have not been used a large number of times. It uses reservoir sampling and favours less frequently used extensions to compensate
   for potential nonuniformity. For large m, kint (m>5) topological sort results in non-repeated extensions anyway, so this selection is not needed. hence set option =0
   Only for small m it may be set to 1.

K - parameter for the k-interactive fuzzy measures from (0,1]

Output: vv matrix of suze fmsize x num, with num rows, where the reandom fuzzy measure values are stored. 
	    The type myfloat needs to be defined as float or double in this file
		The fmsize is found by calling function fm_arraysize( m, 2^m,  kint)
		The fuzzy measures are arranged in cardinality based ordering with the last m-kint values corresponding to the cardinalities kint+1,...m
		the first element is always 0 (v(emptyset))

The function generate_fm_tsort is based on the topological sort of randomised entries of the partial order graph, followed by Markov chain and optionally
reservoir sampling to discourage repeating linear extensions.

The method is outlined n the paper by Beliakov, Cabrerizo, Herrera-Viedma and Wu "Random generation of k-interactive capacities"
currently under review.


The generate_fm_minplus is based on the minimalsPlus method, which
 is based on the program supplied by Elias Combarro and Susana Irene Diaz Rodribuez in October 2019,
which is also published in their paper

"Minimals plus: An improved algorithm for the random generation of linear extensions of partially ordered sets." Information
Sciences, 501:50–67, 2019.

The code can be optimised by not using the matrix representation of the poset, relying on the relation preceeds.Not attempted here due to lack of time.

The function generate_fm_tsort is faster than generate_fm_minplus based on our experiments. It requires less markov steps for convergence to uniform

Compiling:
	while the files binarylattice.cpp and minimalsplus.cpp are automatically included, the library fuzzymeasuretools.cpp has to be compiled separately and linked


Usage:

		m=5; kint=3;
		int_64 n ;

		Preparations_FM(m, &n); // to initialise some global variables. the probram will crash otherwise
		int arraysize = fm_arraysize(m, n, kint);

		myfloat* VV = new myfloat[arraysize*total];

		int option =1;

		myfloat K = 0.1;
		// use either one of the. tsort generally is faster
//		generate_fm_minplus(total, m, kint, Markov, option, K, VV);
		generate_fm_tsort(total,  m,  kint,  Markov,  option,  K, VV);

		// print the resulting fuzzy measures
		for (int i = 0;i < total;i++) {
			for (int j = 0;j < arraysize;j++)
				cout << VV[i*arraysize + j] << " ";
			cout << endl;
		}


Email me if there are any questions. It is a freeware.

Gleb Beliakov, April 2020
gleb@deakin.edu.au

*/

// These definitions can be changed to support other data types
//#include "fuzzymeasuretools.h"


//typedef int_64 myint;
//typedef unsigned int myint;
typedef uint16_t myint;
// any of the above, for m<16 use uint_16

typedef unsigned int uint;

// float or double
typedef float  myfloat;

unsigned int bitweight(int_64 i);

int fm_arraysize(int m, int_64 n, int kint);
// calculates the size of the array to store one k-interctive fuzzy measure



// generate fuzzy measures randomly using topological sort
int generate_fm_tsort(int_64 num, int m, int kint, int markov, int option, myfloat K, myfloat * vv);


// generate fuzzy measures randomly using MinimalsPlus method
int generate_fm_minplus(int_64 num, int m, int kint, int markov, int option, myfloat K, myfloat * vv);


