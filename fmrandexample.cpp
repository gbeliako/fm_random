// fmrandlib.cpp : This file contains the 'main' function. Program execution begins and ends there.
//


#include <iostream>
#include "fuzzymeasuretools.h"
#include "fm_random.h"

using namespace std;

int main(int argc, char *argv[])
	{

		int m = 8;    //Number of variables of the fuzzy measure
		long long int total = 100000; //Number of measures to generate
		int kint = 2;
		int Markov = 100;
		long int NumLinExt = 100;

		if (argc < 4) {
			cerr << "\n\nCorrect Usage: test n iters k Markov > outfile.txt" << endl;
			//		goto L1;
			return 1;
		}
		for (int i1 = 1; i1 < argc; i1++) {
			m = atoi(argv[i1]);
			i1++;
			total = atol(argv[i1]);
			i1++;
			kint = atoi(argv[i1]);
			i1++;
			Markov = atoi(argv[i1]);
			i1++;
		}

		int_64 n = (int_64)1 << m; //Number of coefficients of the fuzzy measure

		Preparations_FM(m, &n);
		int arraysize = fm_arraysize(m, n, kint);

		myfloat* VV = new myfloat[arraysize*total];

		int option =1;

		cout << arraysize << " " << total << endl;

		myfloat K = 0.1;
//		generate_fm_minplus(total, m, kint, Markov, option, K, VV);
		generate_fm_tsort(total,  m,  kint,  Markov,  option,  K, VV);

		if(1)
		for (int i = 0;i < total;i++) {
			for (int j = 0;j < arraysize;j++)
				cout << VV[i*arraysize + j] << " ";
			cout << endl;
		}
		return 0;
	
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
