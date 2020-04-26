/********************* Fuzzy measure toolkit ******************************************

This is a set of useful routines for manipulations with fuzzy measures 
(and other set functions). They include binary encoding of a discrete set 
(as  integers (up to 32 elements in a set)), simple set operations: 
intersection, union, inclusion, difference, etc. various representations of 
fuzzy measures (standard, Moebius), orderings of their values, conversions, 
calculations of Shapley, Banzhaf and other interaction indices, orness, entropy, etc.
Calculation of Choquet and Sugeno integrals for a given input x.

--------------------------------------------------------------------------------------
 *
 *      begin                : May 10 2007
 *		end					 : June 3 2018
 *              version                          : 3.0 
 *              copyright            : (C) 2007-2018 by Gleb Beliakov
 *              email                : gleb@deakin.edu.au
 *
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the Lesser GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
**************************************************************************************/

#define tolerance 0.00001

double ElapsedTime();
void   ResetTime();

struct valindex {
        double v;
        int i;
};
struct Less_than {  
  bool operator()(const valindex& a, const valindex& b) { return a.v < b.v;  }
  };
struct Greater_than {  
  bool operator()(const valindex& a, const valindex& b) {  return a.v > b.v;  }
  };

//#define int_64 unsigned long long
typedef unsigned long long int_64;
extern valindex tempxi[100];

extern double   *m_factorials;  // keeps factorials  n! up to n
extern int              *card;                  // array to keep set cardinalities in binary ordering
extern int *cardpos;   // array to store the indices of elements of different cardinalities in the cardinality ordering

extern int_64 *bit2card; // arrays to transform from one ordering to another
extern int_64 *card2bit;

extern int *cardposm;   // array to store the indices of elements of different cardinalities in the cardinality ordering

extern int_64 *bit2cardm; // arrays to transform from one ordering to another
extern int_64 *card2bitm;

// this routine should be called first to prepare all the arrays
void Preparations_FM(int n, int_64 *m);

void Preparations_FM_marginal(int n, int_64 *m, int Kinter);

void ConvertCard2Bit(double* dest, double* src,  int_64 m);

// this routine should be called last to clean all the arrays
void Cleanup_FM();


/* useful routines */
extern double   minf(double a, double b); 
double  maxf(double a, double b);
int             sign(int i);   // sign of i
int             IsOdd(int i);  // is i odd ?




unsigned int cardf(int_64 i); // count how many bits in i are set
double xlogx(double t);         // x log x (but takes care of x close to 0, using tolerance parameter


/* Set manipulations. A set is represented by an unsigned int (32 bits) */
void    RemoveFromSet(int_64* A, int i);  // remove a from set i
void    AddToSet(int_64* A, int i);               // add a to set i
int             IsInSet(int_64 A, int i);                 // does a belong to set i?
int             IsSubset(int_64 A, int_64 B);       // is j subset of i ?
int_64  Setunion(int_64 A, int_64 B); // returns the  union of sets i and j
int_64  Setintersection(int_64 A, int_64 B); // returns the  intersection of sets i and j
int_64  Setdiff(int_64 A, int_64 B);                  // returns set difference  i \ j
double min_subset(double* x, int n, int_64 S); // returns minimum of x_i, such that i belongs to set S
int_64 UniversalSet(int n);
int Removei_th_bitFromSet(int_64* A, int i);


unsigned int ShowValue(int_64 s); // shows the elements of a subset as a decimal string (up to 10 elements)

double Choquet(double*x, double* v, int n, int_64 m);
/* Calculates the value of a descrete Choquet integral of x, wrt fuzzy measure v. 
   Parameters: x array[n] ,v array[m], n, m=2^n 
   This proceduce requires sorting the array of pairs (v[i],i) in non-decreasing order
   (perfored by standard STL sort function, and RemoveFromSet() function, to remove
   an indicated bit from a set (in its binary representation). 
*/

double ChoquetMob(double*x, double* Mob, int n, int_64 m);
/* This is an alternative calculation of the Choquet integral from the Moebius transform v. 
   It is not as efficient as Choquet(). Provided for testing purposes.
*/

double Sugeno(double*x, double* v, int n, int_64 m);
/* Calculates the value of a descrete Sugeno integral of x, wrt fuzzy measure v 
Parameters: x array[n] ,v array[m], n, m=2^n 
This proceduce requires sorting the array of pairs (v[i],i) in non-decreasing order
(perfored by standard STL sort function, and RemoveFromSet() function, to remove
an indicated bit from a set (in its binary representation)
Also requires maxf and minf functions.
*/

double OWA(double*x, double* v, int n );
/* Calculates the value of OWA */

double WAM(double*x, double* v, int n );
/* Calculates the value of WAM */

void ConstructLambdaMeasure(double *singletons, double *lambda, double *v, int n, int_64 m);
/* Given the values of the fuzzy measure at singletons, finds the appropriate
lambda, and constructs the rest of the fuzzy measure. Returns lambda and v at the output
*/

/* ---------------Operations on fuzzy measures -------------------------*/

double Orness(double* Mob, int n, int_64 m); // calculates orness value of a fuzzy measure

double Entropy(double* v, int n, int_64 m);// calculates Entropy of a fuzzy measure

double OrnessOWA(double* w, int n); // calculates orness value of a fuzzy measure

void Mobius(double* v, double* Mob, int n, int_64 m); // calculates Moebius representation of v
/* the output array w should have the same size 2^n=m as v */

void Zeta(double* Mob, double* v, int n, int_64 m);// calculates inverse Moebius transform

void Shapley(double* v, double* x, int n, int_64 m); // calculates the array x of Shapley values
void Banzhaf(double* v, double* x, int n, int_64 m);// calculates the array x of Banzhaf indices
void Interaction(double* Mob, double* w, int_64 m); // calculates all 2^n interaction indices (returned in w)
void InteractionB(double* Mob, double* w, int_64 m); // calculates all 2^n Banzhaf interaction indices (returned in w)

void NonadditivityIndex(double* v, double* w, int n, int_64 m); // calculates  all 2^n nonadditivity indices (returned in w)
void NonadditivityIndexMob(double* Mob, double* w, int n, int_64 m); // calculates  all 2^n nonadditivity indices (returned in w) using Mobius transform
void BipartitionShapleyIndex(double* v, double* w, int n, int_64 m); // calculates  all 2^n bipartition Shapley indices (returned in w)
void BipartitionBanzhafIndex(double* v, double* w, int n, int_64 m); // calculates  all 2^n bipartition Banzhaf indices (returned in w)

void dualm(double* v, double* w, int n, int_64 m); // calculates the dual fuzzy measure, returns it in w

// Various queries about a fuzzy measure. All performed with a given tolerance. 
int IsMeasureAdditive(double* v, int n, int_64 m);
int IsMeasureBalanced(double* v, int_64 m);
int IsMeasureSelfdual(double* v, int_64 m);
int IsMeasureSubadditive(double* v, int_64 m);
int IsMeasureSubmodular(double* v, int_64 m);
int IsMeasureSuperadditive(double* v, int_64 m);
int IsMeasureSupermodular(double* v, int_64 m);
int IsMeasureSymmetric(double* v, int n, int_64 m);
int IsMeasureKMaxitive(double* v, int n, int_64 m);
