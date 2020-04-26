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




#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <ctime>

#include <numeric>


clock_t clockS, clockF;
double TotalTime;
void   ResetTime() { TotalTime = 0; clockS = clock(); }
double ElapsedTime()
{
	clockF = clock();
	double duration = (double)(clockF - clockS) / CLOCKS_PER_SEC;
	TotalTime += duration;
	clockS = clockF;
	return TotalTime;
}




//#include <R.h>
#include "fuzzymeasuretools.h"
using namespace std;

  Less_than less_than;              /* declare a comparison function object, to */
  Greater_than  greater_than ;      /*  pass to sort and search algorithms */

 valindex tempxi[100];
 double         *m_factorials;  // keeps factorials  n! up to n
 int            *card;                  // array to keep set cardinalities in binary ordering
 int            *cardpos;   // array to store the indices of elements of different cardinalities in the cardinality ordering

 int_64 *bit2card;        // arrays to transform from one ordering to another
 int_64 *card2bit;


 // marginal representation
 int            *cardposm;   // array to store the indices of elements of different cardinalities in the cardinality ordering
 int_64 *bit2cardm;        // arrays to transform from one ordering to another
 int_64 *card2bitm;

int sign(int i) {if(i<0) return -1; else return 1;}
int signd(double i) {if(i<0) return -1; else return 1;}
typedef double ( *USER_FUNCTION)(double );
double bisection(double a, double b, USER_FUNCTION f, int nmax)
{
        double u,v,c,w;
        int i;
        u=f(a); v=f(b);
        if(signd(u)==signd(v)) { return -10e10;} // no solution
        i=nmax;
        while(i>0) {
                i--;
                c=(a+b)/2.0;
                w=f(c);
                if( (b-a) < 1.0e-10 ) break;
                if(signd(u)==signd(w)) {
                        u=w; a=c;
                } else {
                        v=w;b=c;
                }
        }
        return (a+b)/2.0;
}


double minf(double a, double b) {if (a<b) return a; else return b;}
double maxf(double a, double b) {if (a>b) return a; else return b;}

int IsOdd(int i) {return ((i & 0x1)?1:0); }
unsigned int cardf(int_64 A) // count how many bits in i are set
{
        int s=0;
        int_64 t=A;
                while(t>0) {
                        if(t & 0x1) s++;
                        t=(t>>1);
                }
        return s;
}


int_64 UnivSetTable[20] = { 0, 1, 3, 7, 15, 31, 63, 127, 255,511,1023,2047,4095,8191,16383,32767,65535,131071,262143 };
//int_64 UnivSetTable[7] = { 0, 1, 3, 7, 15, 31, 63};


double xlogx(double t) { if(t<tolerance) return 0; else return t*log(t);} 

void RemoveFromSet(int_64* A, int i) { *A &= (~(int_64(1) << i)); }
void AddToSet(int_64* A, int i) { *A |= (int_64(1) << i); }
int  IsInSet(int_64 A, int i) { return int((A >> i) & 0x1); }
int  IsSubset(int_64 A, int_64 B) { return ((A & B) == B); }
int_64 Setunion(int_64 A, int_64 B) { return (A | B); }
int_64 Setintersection(int_64 A, int_64 B) { return (A & B); }
int_64 Setdiff(int_64 A, int_64 B) { return (A & ~(A & B)); }

int Removei_th_bitFromSet(int_64* A, int i) { int_64 B = (*A) & (~(int_64(1) << i)); 
                          if (B == *A) return 1; else *A = B; return 0; 
}

int_64 UniversalSet(int n){
	int_64 A = UnivSetTable[min(n,19)];
    while (n > 19)  AddToSet(&A, --n); return A; 
};


unsigned int     ShowValue(int_64 s) {
                int i,k;
                k=0;
                
                for(i=0;i<9;i++) {
                        if(IsInSet(s,i)) {
                                k *= 10;
                                k += (i+1);
                        }
                }
				for (i = 9; i<10; i++) {
					if (IsInSet(s, i)) {
						k *= 10;
						//k += (i + 1);
					}
				}
                return k;
}

double Choquet(double*x, double* v, int n, int_64 m)
/* Calculates the value of a discrete Choquet integral of x, wrt fuzzy measure v 
Parameters: x array[n] ,v array[m], n, m=2^n 
This proceduce requires sorting the array of pairs (v[i],i) in non-decreasing order
(perfored by standard STL sort function, and RemoveFromSet() function, to remove
an indicated bit from a set (in its binary representation) */
{       double s=0;
        int i;
        for(i=0;i<n;i++) { (tempxi[i]).v=x[i]; (tempxi[i]).i=i;}
        sort(&(tempxi[0]),&(tempxi[n]),less_than); // sorted in increasing order

        int_64 id = m-1; // full set N (11111... in binary)

        s=tempxi[0].v*v[id];
        RemoveFromSet(&id, tempxi[0].i);
        for(i=1;i<n;i++) {
                s+=(tempxi[i].v - tempxi[i-1].v)* v[id];
                RemoveFromSet(&id, tempxi[i].i);
        }
        return s;
}




double ChoquetMob(double*x, double* Mob, int n, int_64 m)
/* This is an alternative calculation of the Choquet integral from the Moebius transform. 
   It is not as efficient as Choquet(). Provided for testing purposes. 
*/
{
    double s=0;
    int_64 A;
    for(A=1; A < m; A++)
        s += Mob[A] * min_subset( x, n, A);
    return s;
}



double Sugeno(double*x, double* v, int n, int_64 m)
/* Calculates the value of a discrete Sugeno integral of x, wrt fuzzy measure v 
   Parameters: x array[n] ,v array[m], n, m=2^n 
   This proceduce requires sorting the array of pairs (v[i],i) in non-decreasing order
   (perfored by standard STL sort function, and RemoveFromSet() function, to remove
   an indicated bit from a set (in its binary representation)
   Also requires maxf and minf functions. 
*/
{
    double s=0;
    int i;
    for(i=0;i<n;i++) { (tempxi[i]).v=x[i]; (tempxi[i]).i=i;}
    sort(&(tempxi[0]),&(tempxi[n]),less_than); // sorted in decreasing order

    int_64 id = m-1; // full set N (11111... in binary)

    s=0;
    for(i=0;i<n;i++) {
        s =maxf(s, minf(tempxi[i].v , v[id]));
        RemoveFromSet(&id, tempxi[i].i);
    }
    return s;
}


double OWA(double*x, double* v, int n )
/* Calculates the value of an OWA 
Parameters: x array[n] ,v array[n], n, 
This proceduce requires sorting the array of pairs (v[i],i) in non-decreasing order
 */
{       double s=0;
        int i;
        for(i=0;i<n;i++) { (tempxi[i]).v=x[i]; (tempxi[i]).i=i;}
        sort(&(tempxi[0]),&(tempxi[n]),less_than); // sorted in increasing order

        for(i=0;i<n;i++) {
                s+=  tempxi[n-i-1].v * v[i];
        }
        return s;
}
double WAM(double*x, double* v, int n )
/* Calculates the value of a WAM 
Parameters: x array[n] ,v array[n], n, 
 */
{       double s=0;
        int i;
        for(i=0;i<n;i++) {
                s+=  x[i] * v[i];
        }
        return s;
}


double auxarray[100];
int auxN;
double auxfun(double lam)
{
        int i;
        double s=1;
        for(i=0;i<auxN;i++) s*= (1 + lam* auxarray[i]);
        s -= (lam+1);
        return s;
}
void ConstructLambdaMeasure(double *singletons, double *lambda, double *v, int n, int_64 m)
/* Given the values of the fuzzy measure at singletons, finds the appropriate
lambda, and constructs the rest of the fuzzy measure. Returns lambda and v at the output
*/
{
        double tol=1.0e-8;
        int i;
        auxN=n;
	 double cond=0;
        double a,b,c;             
		int_64 j;
        double s;
   
	for(i=0;i<n;i++) {cond +=singletons[i]; auxarray[i]=singletons[i];}
        
	 if(fabs(cond-1)<tol) // additive
	{
		*lambda=0;
		c=0;
		goto E1;
	}

  
        a=-1+tol;
        b=0-tol;
        c=bisection(a,b,auxfun,10000);
        if(c<-1) { //means we have to use another interval
                a=tol;
            b=10000;
                c=bisection(a,b,auxfun,100000);
        }
        // so lambda is c now
	 tol *=tol;
	
	if(fabs(c)<tol) goto E1;

        v[0]=0;
        for(j=1;j<m;j++) {
                s=1;
                for(i=0;i<n;i++) if(IsInSet(j, i)) 
                        s *= (1+ c* auxarray[i]);
                s = (s-1)/c;
                v[j]=s;
        }
        *lambda=c;
	return;

// special calse lambda=0

E1:
	*lambda=0;
        v[0]=0;
        for(j=1;j<m;j++) {
                s=0;
                for(i=0;i<n;i++) if(IsInSet(j, i)) 
                        s += auxarray[i];
                
                v[j]=s;
        }

	        
}

double Orness(double* Mob, int n, int_64 m)
{
		int_64 i;
        double s;
        s=0;
        for(i=1;i<m;i++) { 
                s += Mob[i] * (n-card[i])/(card[i]+1.);
        }
        return s/(n-1);
}

double OrnessOWA(double* w,  int n)
{
        int i;
        double s;
        s=0;
        for(i=1;i<=n;i++) { 
                s += w[i-1] * (n-i+0.0)/(n-1.);
        }
        return s;
}

double Entropy(double* v, int n, int_64 m)
{
        int i;
		int_64 id, tempid;
        double s=0;
        double nfac=m_factorials[n];
        for(i=0;i<n;i++) {
                tempid=0;       AddToSet(&tempid,i);
                for(id=0;id<m;id++) if(!IsInSet(id,i)) { 
                        s += -xlogx(v[ Setunion(id,tempid) ] - v[id]) * m_factorials[n-card[id]-1]*m_factorials[card[id]]/nfac;
                }
        }
        return s;
}

void Mobius(double* v, double* Mob, int n, int_64 m)
{
	int_64 i;
	int_64 id;
        double s;
        for(i=0;i<m;i++) {
                s=0;
                for(id=0;id <= i;id++) if(IsSubset(i,id)) {
                        s += v[id] * (IsOdd( card[Setdiff(i,id) ]) ? -1:1); ;
                }
                Mob[i]=s;
        }
}
void Zeta(double* Mob, double* v, int n, int_64 m)
//inverse Moebius transform
{
		int_64 i;
		int_64 id;
        double s;
        for(i=0;i<m;i++) {
                s=0;
                for(id=0;id <= i;id++) if(IsSubset(i,id)) {
                        s += Mob[id] ;
                }
                v[i]=s;
        }
}


void Shapley(double* v, double* x, int n, int_64 m)
{
		int j;
		int_64 i;
		int_64 id;

        for(j=0;j<n;j++) {
                id=0; AddToSet(&id, j); 
                x[j]=0;
                for(i=0;i<m;i++) 
                        if(!IsInSet(i,j)) {
                                x[j] += (m_factorials[n-card[i]-1]*m_factorials[card[i]])/m_factorials[n] *
                                        ( v[ Setunion(i,id) ] - v[i]);
                        }
        }
}

void Banzhaf(double* v, double* x, int n, int_64 m)
{
		int j;
		int_64 i;
		int_64 id;

        for(j=0;j<n;j++) {
                id=0; AddToSet(&id, j); 
                x[j]=0;
                for(i=0;i<m;i++) 
                        if(!IsInSet(i,j)) {
                                x[j] +=  (v[ Setunion(i,id) ] - v[i]);
                        }
                x[j] /= (1<<(n-1));
        }
}


void Interaction(double* Mob, double* w, int_64 m)
{
	int_64 j, i;
	int_64 id;
        for(i=0;i<m;i++) {
                w[i]=0;
                j=card[i];
                for(id=i;id<m;id++) // supersets only
                        if(IsSubset(id,i))
                                w[i]+=Mob[id]/( card[id]- j +1);
        }
}


void InteractionB(double* Mob, double* w, int_64 m)
{
	int_64 j, i;
	int_64 id;
        for(i=0;i<m;i++) {
                w[i]=0;
                j=card[i];
                for(id=i;id<m;id++) // supersets only
                        if(IsSubset(id,i))
                                w[i]+=Mob[id]/ (1<<( card[id]- j ));
        }
}



void NonadditivityIndex(double* v, double* w, int n, int_64 m) // calculates  all 2^n nonadditivity indices (returned in w)
{
	int_64 j, i;
	double r;
	int_64 id;
	w[0] = 0;
	for (i = 1; i<m; i++) {
		w[i] = 0;
		j = card[i];

		if (j == 1){ w[i] = v[i]; }
		else{

			r = (j>1) ? 1. / ((int(1) << (j - 1)) - 1.) : 1;
			for (id = 1; id < i; id++) // subsets only
			if (IsSubset(i, id))
				w[i] += v[id];
			w[i] *= r;
			w[i] = v[i] - w[i];
		}
	}
}


void NonadditivityIndexMob(double* Mob, double* w, int n, int_64 m) // calculates  all 2^n nonadditivity indices (returned in w) using Mobius transform
{
	unsigned int j; int_64 i;
	double r;
	int_64 id;
	w[0] = 0;
	for (i = 1; i<m; i++) {
		w[i] = 0;
		j = card[i];
		
		for (id = 0; id<i; id++) // subsets only
		if (IsSubset(i, id))
		{
			r = (j>1) ? (  (int(1) << (j - 1)) - (int(1) << (j-card[id]))  ) / ((int(1) << (j - 1)) - 1.)  : 1;
			w[i] += Mob[id]*r;
		}
		w[i] += Mob[i] ;
	}
}

void CalculateDeltaHat(double* v, double* w, int_64 A, int_64 B, int_64 m)
{
	unsigned int j;
	int_64 id;

	*w = 0;	j = card[A];
	if (j == 0) return;
	if (j == 1) { *w = v[Setunion(A, B)] - v[B];  return; }

	for (id = 1; id < A; id++)
	if (IsSubset(A, id)){
		*w += v[ Setunion(B, id) ];
	}

	*w *= 1. / ( (int(1)<<(j-1)) -1.);
	id = Setunion(A, B);
	*w = v[id] + v[B] - *w;
}
double ChooseInverse(int n, int a, int b){ return m_factorials[b]*m_factorials[n-a-b] / m_factorials[n-a];  }

void BipartitionShapleyIndex(double* v, double* w, int n, int_64 m) // calculates  all 2^n bipartition Shapley indices (returned in w)
{
	unsigned int j,k;
	int_64 A, B;
	double r,d,ch;
	w[0] = 0;
	for (A = 1; A < m; A++) {
		w[A] = 0;
		j = card[A];
		r = 1. / (n - j + 1.);
		for (B = 0; B < m; B++)
		if (Setintersection(A, B) == 0){
			CalculateDeltaHat(v, &d, A, B, m);
			k = card[B];
			ch = ChooseInverse(n, j, k);
			w[A] += r* ch* d;
		}
	}
}
void BipartitionBanzhafIndex(double* v, double* w, int n, int_64 m) // calculates  all 2^n bipartition Shapley indices (returned in w)
{
	unsigned int j;
	int_64 A, B;
	double r, d;
	w[0] = 0;
	for (A = 1; A < m; A++) {
		w[A] = 0;
		j = card[A];
		r = 1. / (int(1)<<(n-j));
		for (B = 0; B < m; B++)
		if (Setintersection(A, B) == 0){
			CalculateDeltaHat(v, &d, A, B, m);
			w[A] += r* d;
		}
	}
}


void dualm(double* v, double* w, int n, int_64 m)
{
        unsigned int i;
        for(i=0;i<m;i++)
        {
                w[ (~i) & (m-1) ] = 1-v[i];
        }
}


int IsMeasureAdditive(double* v, int n, int_64 m)
{       
	int j;
	int_64 i;
        double s;
        for(i=3;i<m;i++) {
                if(card[i]>1) {
                        s=0;
                        for(j=0;j<n;j++)
                         if(IsInSet(i,j)) s+=v[(int_64) 1<<j ];
                        if(fabs(s-v[i])>tolerance) return 0;
                }
        }
        return 1;
}

int IsMeasureKMaxitive(double* v, int n, int_64 m)
{
	int j;
	int_64 i;
	double s;
	int K = 1;
	for (i = 1; i<m; i++) {
		if (card[i]>1) {
			s = 0;
			for (j = 0; j<n; j++)
			if (IsInSet(i, j)) s = maxf(v[Setdiff(i,1ULL <<j)], s);
			if (fabs(s - v[i])>tolerance) K=max(K,card[i]); // fails for cardinality K
		}
	}
	return K;
}

int IsMeasureBalanced(double* v, int_64 m)
{
	int_64 i, j;
        for(i=0;i<m;i++) {
                for(j=i;j<m;j++) {
                        if((card[i]<card[j]) && (v[i]>v[j])) return 0;
                        if((card[i]>card[j]) && (v[i]<v[j])) return 0;
                }
        }
        return 1;
}


int IsMeasureSelfdual(double* v, int_64 m)
{       
	int_64 i;
        for(i=0;i<m;i++)
        {
                if(fabs(v[ (~i) & (m-1) ] + v[i]-1) > tolerance) return 0;
        }
        return 1;
}


int IsMeasureSubadditive(double* v, int_64 m)
{
	int_64 i, j;
        for(i=0;i<m;i++) {
                for(j=i+1;j<m;j++) if(Setintersection(i,j)==0) {
                        if(v[i]+v[j] - v[Setunion(i,j)] < -tolerance) return 0;
                }
        }
        return 1;
}


int IsMeasureSubmodular(double* v, int_64 m)
{
	int_64 i, j;
        for(i=0;i<m;i++) {
                for(j=i+1;j<m;j++) if(Setintersection(i,j)==0) {
                        if(v[i]+v[j] - v[Setunion(i,j)]- v[Setintersection(i,j)] < -tolerance) return 0;
                }
        }
        return 1;

}


int IsMeasureSuperadditive(double* v, int_64 m)
{
	int_64 i, j;
        for(i=0;i<m;i++) {
                for(j=i+1;j<m;j++) if(Setintersection(i,j)==0) {
                        if(v[i]+v[j] - v[Setunion(i,j)] > tolerance) return 0;
                }
        }
        return 1;
}


int IsMeasureSupermodular(double* v, int_64 m)
{
	int_64 i, j;
        for(i=0;i<m;i++) {
                for(j=i+1;j<m;j++) {
                        if(v[i]+v[j] - v[Setunion(i,j)] - v[Setintersection(i,j)] > tolerance) return 0;
                }
        }
        return 1;
}


int IsMeasureSymmetric(double* v, int n, int_64 m)
{
	int_64 i, j;
        double *w=new double[n+1];
        for(i=0;i<(unsigned int)n+1;i++) w[i]=-1;

        for(i=0;i<m;i++) {
                j=card[i];
                if(w[j]<0) w[j]=v[i]; else
                        if(fabs(w[j]-v[i])>tolerance) {delete[] w; return 0; }
        }
        delete[] w;
        return 1;
}


double min_subset(double* x, int n, int_64 S)
{ // returns min x_i when i \in S, or 0 if S is empty
        int i;
        double r=10e10;
        for(i=0;i<n;i++)
                if( IsInSet(S,i)) r=minf(r,x[i]);
        if(r>1) r=0;
        return r;
}



void ConvertCard2Bit(double* dest, double* src, int_64 m)
{
	for (int_64 i = 0; i < m; i++)
		dest[card2bit[i]] = src[i];
}

// this is a recursive procedure which helps build all subsets of a given cardinality, and 
// set up conversion arrays
void recursive_card(int_64* k, unsigned int level, unsigned int maxlevel,
                                        unsigned int start, unsigned int finish,
										int_64* b2c, int_64* c2b, int_64 *s, int n)
{
        unsigned int i1;
        for(i1=start; i1 <= finish; i1++) { AddToSet(s,i1);
                if(level == maxlevel) {
                        b2c[*s]=*k;
                        c2b[*k]=(*s);
                        (*k)++;
                } else {
                        recursive_card(k,level+1,maxlevel,i1+1,finish+1,b2c,c2b,s,n);
                }
                RemoveFromSet(s,i1);
        }
}


void main_card(int_64* k, unsigned int level, int_64* b2c, int_64* c2b, int n)
{
        // we recursively construct all subsets of cardinality "level"
		int_64 s = 0;
        recursive_card(k,1,level,0, n-level, b2c,c2b, &s,n);
}



// this is a recursive procedure which helps build all subsets of a given cardinality, and 
// set up conversion arrays
void recursive_card_marginal(int_64* k, unsigned int level, unsigned int maxlevel,
	unsigned int start, unsigned int finish,
	int_64* b2c, int_64* c2b, int_64 *s, int n)
{
	unsigned int i1;
	for (i1 = start; i1 <= finish; i1++) {
		AddToSet(s, i1);
		if (level == maxlevel) {
			b2c[*s] = *k;
			c2b[*k] = (*s);
			for (unsigned int r = 1; r < level; r++)
			{
				c2b[++(*k)] = (*s);
			}
			(*k)++;
		}
		else {
			recursive_card_marginal(k, level + 1, maxlevel, i1 + 1, finish + 1, b2c, c2b, s, n);
		}
		RemoveFromSet(s, i1);
	}
}


void main_card_marginal(int_64* k, unsigned int level, int_64* b2c, int_64* c2b, int n, int Kinter)
{
	// we recursively construct all subsets of cardinality "level"
	int_64 s = 0;
	recursive_card_marginal(k, 1, level, 0, n - level, b2c, c2b, &s, Kinter);
}

void Preparations_FM(int n, int_64 *m)
{
        int i;
        int_64 j;
        *m= (int_64)1<<(n);

    // calculate the array containing factorials of i! (faster than calculating them every time)
    m_factorials=new double[n+1];
        m_factorials[0]=1;
        for(i=1;i<=n;i++) m_factorials[i] = m_factorials[i-1]*i;

    // this array will contains cardinailities of subsets (coded as binaries), i.e. the number of bits in i.
        card=new int[(int) *m];
        cardpos=new int[n+1];
        card[0]=0;
        for(j=1;j<*m;j++) card[j] = cardf(j);

// these two arrays are used to pass from binary to cardinality ordering
// they are precomputed 
// in binary ordering the subsets are ordered as
// 0 1 2 12 3 13 23 123 4 14 24 124 34 134 234 1234,...
// (which corresponds to the order 0,1,2,3,... in binary form)
// in cardinality ordering they are ordered as
// 0 1 2 3 4 5 6 12 13 14 15 16 23 24 25 26 34 35 36 45 46 56 123 124,...
// (empty, singletons, pairs,triples, etc.)
// for a given subset s in cardinality ordering, to find its binary code use  card2bit[s]
// and vice versa
// cardpos[i] is the index at which subsets with cardinality i+1 start in the cardinality ordering
// i.e. cardpos[0]..cardpos[1]-1 - singletons, cardpos[1]..cardpos[2]-1 - pairs, etc.

		bit2card = (int_64*) new int_64[*m];
		card2bit = (int_64*) new int_64[*m];

		int_64 k; int l;
        bit2card[0]=card2bit[0]=0;

        cardpos[0]=1; // positions where singletons start, the 0th element is empyset

        k=1;
        for(l=1;l<=n-1;l++) {
                main_card(&k, l, bit2card, card2bit,  n);
                cardpos[l]=(int)k;
        }
        cardpos[n]=cardpos[n-1]+1;
        
        bit2card[*m-1]=card2bit[*m-1]=*m-1;

		card2bitm=NULL;
		bit2cardm=NULL;
		cardposm=NULL;
}


void Preparations_FM_marginal(int n, int_64 *m, int Kinter)
{
	int i;
//	int_64 j;
	int sz=n;
	*m = (int_64)1 << (n);

	// do standard job, factorials and card
	Preparations_FM(n, m);

	//*m *= 2;

	// calculate the array containing factorials of i! (faster than calculating them every time)
//	m_factorials = new double[n + 1];
//	m_factorials[0] = 1;
//	for (i = 1; i <= n; i++) m_factorials[i] = m_factorials[i - 1] * i;

	// this array will contains cardinailities of subsets (coded as binaries), i.e. the number of bits in i.
//	card = new int[(int)*m];

//	card[0] = 0;
//	for (j = 1; j<*m; j++) card[j] = cardf(j);

	// these two arrays are used to pass from binary to cardinality ordering
	// they are precomputed 
	// in binary ordering the subsets are ordered as
	// 0 1 2 12 3 13 23 123 4 14 24 124 34 134 234 1234,...
	// (which corresponds to the order 0,1,2,3,... in binary form)
	// in cardinality ordering they are ordered as
	// 0 1 2 3 4 5 6 12 13 14 15 16 23 24 25 26 34 35 36 45 46 56 123 124,...
	// (empty, singletons, pairs,triples, etc.)
	// for a given subset s in cardinality ordering, to find its binary code use  card2bit[s]
	// and vice versa
	// cardpos[i] is the index at which subsets with cardinality i+1 start in the cardinality ordering
	// i.e. cardpos[0]..cardpos[1]-1 - singletons, cardpos[1]..cardpos[2]-1 - pairs, etc.

	cardposm = new int[n + 1];
	for (i = 2; i <= Kinter; i++) sz += (int)(m_factorials[n] / m_factorials[i - 1] / m_factorials[n - i]);


	bit2cardm = new int_64[*m];
	card2bitm = new int_64[sz*2];

	int_64 k; int l;
	bit2cardm[0] = card2bitm[0] = 0;

	cardposm[0] = 1; // positions where singletons start, the 0th element is empyset

	k = 1;
	for (l = 1; l <= Kinter; l++) {
		main_card_marginal(&k, l, bit2cardm, card2bitm, n, Kinter);
		cardposm[l] = (int)k;
	}
	if(Kinter<n)	cardposm[n] = cardposm[n - 1] + 1;

//	bit2card[*m - 1] = card2bit[*m - 1] = *m - 1;
}


void Preparations_FM(int n, int_64 *m, int Kinteractive)
{
	// todo: complete this procedure, and provide a set of tools for Choque, conversions and others to work with K-interactive measures
	int i;
	int_64 j;

	if (Kinteractive > n)Kinteractive = n;
	if (Kinteractive < 1)Kinteractive = 1;



	// calculate the array containing factorials of i! (faster than calculating them every time)
	m_factorials = new double[n + 1];
	m_factorials[0] = 1;
	for (i = 1; i <= n; i++) m_factorials[i] = m_factorials[i - 1] * i;

	*m = 1;
	for (i = 1; i <= Kinteractive; i++) *m += (int)(m_factorials[n] / m_factorials[i] / m_factorials[n - i]);
	*m +=  n-Kinteractive;



	// this array will contains cardinailities of subsets (coded as binaries), i.e. the number of bits in i.
	card = new int[(int)*m];
	cardpos = new int[n + 1];
	card[0] = 0;

	for (j = 1; j<(*m - (n - Kinteractive)); j++) card[j] = cardf(j);

	for (j = Kinteractive + 1; j <= (int_64)n; j++) card[j + (*m - (n - Kinteractive)) - Kinteractive - 1] = (int)j;

	// these two arrays are used to pass from binary to cardinality ordering
	// they are precomputed 
	// in binary ordering the subsets are ordered as
	// 0 1 2 12 3 13 23 123 4 14 24 124 34 134 234 1234,...
	// (which corresponds to the order 0,1,2,3,... in binary form)
	// in cardinality ordering they are ordered as
	// 0 1 2 3 4 5 6 12 13 14 15 16 23 24 25 26 34 35 36 45 46 56 123 124,...
	// (empty, singletons, pairs,triples, etc.)
	// for a given subset s in cardinality ordering, to find its binary code use  card2bit[s]
	// and vice versa
	// cardpos[i] is the index at which subsets with cardinality i+1 start in the cardinality ordering
	// i.e. cardpos[0]..cardpos[1]-1 - singletons, cardpos[1]..cardpos[2]-1 - pairs, etc.

	bit2card = new int_64[*m];
	card2bit = new int_64[*m];

	int_64 k; int l;
	bit2card[0] = card2bit[0] = 0;

	cardpos[0] = 1; // positions where singletons start, the 0th element is empyset

	k = 1;
	for (l = 1; l <= n - 1; l++) {
		main_card(&k, l, bit2card, card2bit, n);
		cardpos[l] = (int)k;
	}
	cardpos[n] = cardpos[n - 1] + 1;

	bit2card[*m - 1] = card2bit[*m - 1] = *m - 1;
}




void Cleanup_FM()
{
        delete [] card2bit;
        delete [] bit2card;

        delete [] m_factorials;
        delete [] card;
        delete [] cardpos;

		if (card2bitm!=NULL) delete[] card2bitm;
		if (bit2cardm != NULL)delete[] bit2cardm;
		if (cardposm != NULL)delete[] cardposm;
}
