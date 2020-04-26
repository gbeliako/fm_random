/*
/*
	See .h for description

*/



#include<vector>
#include <map>
#include<random>
#include<iostream>
#include <fstream>
#include<algorithm>
#include <string>
#include <unordered_map>
#include "fuzzymeasuretools.h"

#include "fm_random.h"


using namespace std;

random_device rd;
mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)

std::uniform_real_distribution<double> distribution(0.0,1.0);


// These definitions can be changed to support other data types. They sin in the .h file

//typedef uint16_t myint;
//typedef unsigned int myint;
//typedef int_64 myint;
//typedef unsigned int uint;
//typedef float  myfloat;


//int __builtin_popcount (unsigned int x);
#ifdef __GNUC__
unsigned int bitweight(int_64 i) {
	return __builtin_popcountl(i);
}

#elif _MSC_VER
#  include <intrin.h>

#ifdef  _WIN64
#  define __builtin_popcountl  __popcnt64  //_mm_popcnt_u64
inline unsigned int bitweight(int_64 i) {
	return __builtin_popcountl(i);
}
#else
#  define __builtin_popcountl  __popcnt  //_mm_popcnt_u64
inline unsigned int bitweight(int_64 i) {
	return __builtin_popcountl((uint32_t)(i >> 32)) + __builtin_popcountl((uint32_t)i);
}
#endif

#else 
uint bitweight( int_64 v) {
	v = v - ((v >> 1) & (int_64)~(int_64)0 / 3);                           // temp
	v = (v & (int_64)~(int_64)0 / 15 * 3) + ((v >> 2) & (int_64)~(int_64)0 / 15 * 3);      // temp
	v = (v + (v >> 4)) & (int_64)~(int_64)0 / 255 * 15;                      // temp
	unsigned int c = (int_64)(v * ((int_64)~(int_64)0 / 255)) >> (sizeof(int_64) - 1) * CHAR_BIT; // count
	return (unsigned int)c;
}
//#endif
/*
unsigned int bitweight(unsigned int i)
{
	 i = i - ((i >> 1) & 0x55555555);
	 i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
	 return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}*/
#endif
	


// the functions below are for encoding an array into strings to be used as keys in unordered_map structure (hash keys)
typedef unsigned char BYTE;

std::string base64_encode(BYTE const* buf, unsigned int bufLen);
std::vector<BYTE> base64_decode(std::string const&);

static const std::string base64_chars = 
             "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
             "abcdefghijklmnopqrstuvwxyz"
             "0123456789+/";


static inline bool is_base64(BYTE c) {
  return (isalnum(c) || (c == '+') || (c == '/'));
}

std::string base64_encode(BYTE const* buf, unsigned int bufLen) {
  std::string ret;
  int i = 0;
  int j = 0;
  BYTE char_array_3[3];
  BYTE char_array_4[4];

  while (bufLen--) {
    char_array_3[i++] = *(buf++);
    if (i == 3) {
      char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
      char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
      char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
      char_array_4[3] = char_array_3[2] & 0x3f;

      for(i = 0; (i <4) ; i++)
        ret += base64_chars[char_array_4[i]];
      i = 0;
    }
  }

  if (i)
  {
    for(j = i; j < 3; j++)
      char_array_3[j] = '\0';

    char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
    char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
    char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
    char_array_4[3] = char_array_3[2] & 0x3f;

    for (j = 0; (j < i + 1); j++)
      ret += base64_chars[char_array_4[j]];

    while((i++ < 3))
      ret += '=';
  }

  return ret;
}

std::vector<BYTE> base64_decode(std::string const& encoded_string) {
  int in_len = encoded_string.size();
  int i = 0;
  int j = 0;
  int in_ = 0;
  BYTE char_array_4[4], char_array_3[3];
  std::vector<BYTE> ret;

  while (in_len-- && ( encoded_string[in_] != '=') && is_base64(encoded_string[in_])) {
    char_array_4[i++] = encoded_string[in_]; in_++;
    if (i ==4) {
      for (i = 0; i <4; i++)
        char_array_4[i] = base64_chars.find(char_array_4[i]);

      char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
      char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
      char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

      for (i = 0; (i < 3); i++)
          ret.push_back(char_array_3[i]);
      i = 0;
    }
  }

  if (i) {
    for (j = i; j <4; j++)
      char_array_4[j] = 0;

    for (j = 0; j <4; j++)
      char_array_4[j] = base64_chars.find(char_array_4[j]);

    char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
    char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
    char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

    for (j = 0; (j < i - 1); j++) ret.push_back(char_array_3[j]);
  }

  return ret;
}



int_64 remove1bit(int_64 a, int i)
// counting fom 0
{
	int_64 b=1;
	b=b<<i;
	b=~b;
	b=a & b;
	return b;
}


// various versions of the preceeds relation in the partial order

int preceeds(myint i, myint j)
{
	myint b= i & j;
	if(b==i) // i subset of j
	{
		if(bitweight( int_64(~i & j))==1) return 1;
	}
	if(b==j)
	{
		if(bitweight(int_64(~j & i))==1) return -1;
	}
	return 0;
}

int preceedsP(myint i, myint j, myint r)
{
	if(j==r) return 1;
	if(i==r) return -1;
	return preceeds(i,j);
}

int preceedsa(myint i, myint j)
{
	myint b= i & j;
	if(b==i) // i subset of j
	{
		return 1;
	}
	if(b==j)
	{
	    return -1;
	}
	return 0;
}

int preceedsaP(myint i, myint j, myint r)
{
	if(j==r) return 1;
	if(i==r) return -1;
	return preceedsa(i,j);
}

//does not seem to work
inline bool preceedsB(myint i, myint j)
{
	// return 1 if i subset of j
	return ((i & j) == i);
}
inline double sqr(double a) {return a*a;}




// Data structure to store graph edges
struct Edge {
	myint src, dest;
};

// Class to represent a graph object
class Graph
{
public:	
	// construct a vector of vectors to represent an adjacency list
	vector<vector<myint>> adjList;

	// Graph Constructor
	Graph(vector<Edge> const &edges, int N)
	{
		// resize the vector to N elements of type vector<int>
		adjList.resize(N);

		// add edges to the Directed graph
		for (auto &edge: edges)
			adjList[edge.src].push_back(edge.dest);
	}
	void clear()
	{
		for (auto gi = adjList.begin();gi < adjList.end();gi++)
			gi->clear();
		
		adjList.clear();
	}	
};

// Perform DFS on graph and set departure time of all
// vertices of the graph
void DFS(Graph const &graph, int v, vector<bool> &discovered, vector<int> &departure, int& time)
{
	// mark current node as discovered
	discovered[v] = true;

	// set arrival time
	time++;

	// do for every edge (v -> u)
	for (int u : graph.adjList[v])
	{
		// u is not discovered
		if (!discovered[u])
			DFS(graph, u, discovered, departure, time);
	}
	
	// ready to backtrack
	// set departure time of vertex v
	departure[time] = v;
	time++;
}

// performs Topological Sort on a given DAG
void doTopologicalSort(Graph const& graph, int N, vector<myint>& v,vector<myint>& v1)
{
	// departure[] stores the vertex number using departure time as index
	vector<int> departure(2*N, -1);
	
	// Note if we had done the other way around i.e. fill the
	// array with departure time by using vertex number
	// as index, we would need to sort the array later

	// stores vertex is discovered or not
	vector<bool> discovered(N);
	int time = 0;

	// perform DFS on all undiscovered vertices
	for (int i = 0; i < N; i++)
		if (!discovered[i])
			DFS(graph, i, discovered, departure, time);
	
	// Print the vertices in order of their decreasing
	// departure time in DFS i.e. in topological order
	for (int i = 2*N - 1; i >= 0; i--) {
		if (departure[i] != -1)
		{
			//cout << departure[i] << " ";
			v.push_back(v1[departure[i]]);
		}
	}
}

void DoMarkovChain(vector<myint>& v, int k, myint r)
{
	uniform_int_distribution<int> uni(0,v.size()-2); // guaranteed unbiased
//	uniform_int_distribution<int> coin(0,1); // guaranteed unbiased

	
	for(int j = 0;j < k; j++)
	{
		//if(coin(rng))
		{
			int pos = uni(rng);	
			if(preceedsaP(v[pos] , v[pos+1], r)==0)
			//if(((v[pos] & v[pos+1])!= v[pos]))
				std::swap(v[pos] , v[pos+1]);
		}
	} 
	
	
}

/*  construct binary lattrices and other posets for k-interactive capacities */

#include "binarylattice.cpp"



/* Combarro, Diaz method included here
*/

#include "minimalsplus.cpp"


void random_coefficients(int n, vector<myfloat> & c)
//Generates a vector of n random real numbers 1=X1>=X1>=...>=Xn=0
{
	uniform_real_distribution<> dis(0.0, 1.0);
	
	//c[0] = 1.0;
	//c[1] = 0.0;
	for (int i = 0; i < n; i++)
		c[i] = (myfloat) dis(rng);

	sort(c.begin(), c.end(), greater<myfloat>());
}



int fm_arraysize(int m, int_64 n, int kint)
{
	// calculates the number of parameers needed in cardinal representation for kinteractive capacity

	int extra = m - kint;
	if (kint >= m) extra = 1;
	// count the number of items as the sum of Cin

	int r = 1;
	for (int i = 1;i <= kint; i++)
		r += (int)(choose(i, m));

	r += extra;  // for emptyset
	return r;
}

myfloat fm_delta(int m, int kint, myfloat K)
{  // delta is the fixed marginal contribution in the k-interactive fuzzy measures
	if (m <= kint + 1) return 0;
	return (myfloat)(1.0 - K) / (m - kint - 1);
}


int generate_fm_minplus(int_64 num, int m, int kint, int markov, int option, myfloat K, myfloat * vv)
{
	/* generates num random vectors representing kinteractive fuzzy measures of m arguments in cardinality ordering
	uses markov iterations of markov chain
	kint and K are parameters for k-interactivity
	option reserved, unused
	vv output, needs to be allocated by the calling routine, which also needs to call Preparations_FM(m, &n);
	method: uses minimalsPlus algorithm by Combarro et al followed by Markov chain
	*/
	unordered_map<string, myint> mymap;

	int_64 n = (int_64)1 << m; //Number of coefficients of the fuzzy measure
	int r;

	r = (int) n;// initially for compatibility

//	Preparations_FM(m, &n); in the calling routine

	int arraysize = fm_arraysize(m, n, kint);
	myfloat delta = fm_delta(m, kint, K);

	vector<bool> P = booleanlatticerestricted(m, kint, r);

//	vector<bool> PNR = booleanlatticerestrictednonredundant(m, kint, r);

	/* for minplus */
/*
	vector<float> A((int)n*m*r);
	vector<float> b(m*r);
	vector<int> dir(m*r);

	int numcon = convertintomatrix(PNR, A, b, dir, r);

	cout << numcon << endl;
	*/

	int j1, j;

	int length = n - 2;
	length = r - 1;// no need 0 and 1, but the  top one is needed. Although ay be not, but just leave it for now.

	vector<myint> w = losw(P, r);

	// for counting linear xtensions
	string s5, s6;
	std::vector<BYTE> decodedData1;

	vector<myint> v1;

	for (auto i = 0;i < length;i++) v1.push_back(card2bit[i + 1]);

	vector<myfloat> coef(length);

	for (auto i = 0; i < num; i++)
	{
		vector<myint> le = minimals_w(P, w, r);
	
		vector<myint> new_le = markovKKclassic(P, r, le, markov);

		/*  this is for counting different linear extensions
		s5 = base64_encode((BYTE *)&new_le[0], length * sizeof(myint));
		unordered_map<string, myint>::iterator it = mymap.find(s5);
		if (it != mymap.end())
			it->second++; else
			mymap[s5] = 1;
		*/
/**/
		 random_coefficients(length, coef);

		// generate on simplex
		vv[i* arraysize + 0] = 0; //emptyset
		for ( j = 0;j < length;j++) vv[i* arraysize +  new_le[j]  ] = coef[j] *K;
		for ( j = arraysize - 1,  j1 = m; j1 >=  kint + 1; j--, j1--) vv[i*arraysize + j] = (myfloat) (1.0 - delta * (m - j1));
		/**/

	}

	return 0;
}


int generate_fm_tsort(int_64 num, int m, int kint, int markov, int option, myfloat K, myfloat * vv)
{
	/* generates num random vectors representing kinteractive fuzzy measures of m arguments in cardinality ordering
	uses markov iterations of markov chain
	kint and K are parameters for k-interactivity
	option 1 if using rejection (memory hungry) 0 otherwise
	vv output, needs to be allocated by the calling routine, which also needs to call Preparations_FM(m, &n);

	Method: randomised topological sotring of the preorder to get a linear extension, followed by Markov chain and optionally by rejection mthod,
	which records all generated linear extensions and attempts to reject the ones that have been already seen. Uses reservoir sampling for that.
	*/
	unordered_map<string, myint> mymap;

	int_64 n = (int_64)1 << m; //Number of coefficients of the fuzzy measure
	int r;

	r = (int)n;// initially for compatibility

	//Preparations_FM(m, &n);  done before in the calling routine

	int dorejection = 0;

	int arraysize = fm_arraysize(m, n, kint);
	myfloat delta = fm_delta(m, kint, K);


//	vector<bool> P = booleanlatticerestricted(m, kint, r);
//	vector<bool> PNR = booleanlatticerestrictednonredundant(m, kint, r);
// We do not need matrix representation of adjacency, just count the number of vertices r
	sizeindependent(m, kint, r);


	vector<myint>  v;

	int length = n - 2;
	length = (int) r - 2;// no need 0 and 1, but the  top one is needed. Although ay be not, but just leave it for now.

//	cout << length << endl;

	if (option == 1 && length*num < 20*200000) dorejection = 1;

	vector<myint>  v1, v2, v0(length);

	string s5, s6;

	std::vector<BYTE> decodedData1;


	int NN = 0;
	myint i,j1,j;
	double Wei, WeiS, p, u;
	Edge E;
	vector<Edge> edges;

	for ( i = 0;i < length;i++) v1.push_back(card2bit[i + 1]);

	vector<myfloat> coef(length);

	for (auto j2 = 0;j2 < num;j2++) {

		NN = 0;
		WeiS = 0;
		edges.clear();

		std::shuffle(std::begin(v1), std::end(v1), rng);
	
		for ( i = 0; i < length;i++) {
			for ( j = i + 1; j < length;j++)
				if (i != j) switch (preceedsP(v1[i], v1[j], card2bit[length])) {
				case 1:
					E.src = i; E.dest = j; edges.push_back(E);
					break;
				case -1:
					E.src = j; E.dest = i; edges.push_back(E);
					break;
					//case 1: on1<<i<<" "<<j<<endl; break;	
					//case -1: on1<<j<<" "<<i<<endl; break;	
				}

		}

		// topological sort here
		Graph graph(edges, length);

		doTopologicalSort(graph, length, v2, v1);
		graph.clear();

		DoMarkovChain(v2, markov, card2bit[length]);

		if (dorejection) {
			s5 = base64_encode((BYTE *)&v2[0], length * sizeof(myint));
			s6 = s5;
			unordered_map<string, myint>::iterator it = mymap.find(s5);
			if (it != mymap.end()) Wei = 1. / sqr(sqr(it->second + 0.4)); else Wei = 1;
			WeiS = Wei;
			//			it->second++; else
			//		mymap[s5]=1;

			for ( i = 1;i < length;i++) {
				if (preceedsaP(v2[i - 1], v2[i], card2bit[length]) == 0) {
					NN++;
					std::swap(v2[i - 1], v2[i]);
					s5 = base64_encode((BYTE *)&v2[0], length * sizeof(myint));

					/* reservoir sampling */
					it = mymap.find(s5);
					if (it != mymap.end()) Wei = 1. / sqr(sqr(it->second + 0.4)); else Wei = 1;
					WeiS += Wei;

					p = Wei / WeiS;
					u = distribution(rng);

					if (u <= p) s6 = s5;

					std::swap(v2[i - 1], v2[i]);
				}

			}// for i
			it = mymap.find(s6);
			if (it != mymap.end()) it->second++;
			else
				mymap[s6] = 1;

			//decode s6
			decodedData1 = base64_decode(s6);
			for( i=0;i<length;i++) v2[i]= *((myint*) (&decodedData1[i*sizeof(myint)]));	
		}  // do rejection

		random_coefficients(length,coef);
		// use v2
		vv[j2* arraysize + 0] = 0; //emptyset
//		for (i = 0; i < length; i++) cout<<" "<<bit2card[v2[i]]<<" ";
//		cout << endl;
		for ( i = 0; i < length; i++) vv[j2* arraysize + bit2card[v2[i]]] = coef[length-i-1] * K;
		for ( i = arraysize - 1,  j1=m;  j1 >= kint+1;    i--, j1--) vv[j2*arraysize + i] = myfloat( 1.0 - delta * (m-j1));

		v2.clear();
	}
	//cout<<mymap.size()<<" "<<num<< endl;

	return 0;
}


