/* 
   to be included into fm_random
   a few methods to get a matrix representation of the prtial order, needed by the minimalspls method

*/


inline int_64 choose(int i, int n)
{
	return (int_64)(m_factorials[n]/m_factorials[i]/m_factorials[n-i]);
}


void sizeindependent(int n, int k, int &r)
{
	// count the number of items as the sum of Cin
//	long int
	r = 1;
	for (int i = 1;i <= k; i++)
		r += choose(i, n);

	r += 1;
}

vector<bool> booleanlatticerestricted(int n, int k, int &r)
{
//Devuelve la matriz del retículo booleano Bn
	long int s = (int_64)1<<n;
	// count the number of items as the sum of Cin

	sizeindependent(n,k,r);

	vector<bool> B(r*r,false);
	
	for(int i=1;i<r-1;i++)
	{
		int_64 a = card2bit[i];
		for(int j=0;j<i;j++){
			int_64 b = card2bit[j];
			
			if(IsSubset(a,b)) B[r*i+j]=true;
		}
		
	}
	for(int i=0;i<r;i++)B[r*i+i]=true;  // reflexive
	
	for(int i=0;i<r;i++)B[r*(r-1)+i]=true;  // top element
	

	return B;
}	


int preceedsby1(int_64 i, int_64 j)
{
	int_64 b= i & j;
	if(b==i) // i subset of j
	{
		if(bitweight(int_64(~i & j))==1) return 1;
	}
	if(b==j)
	{
		if(bitweight(int_64(~j & i))==1) return -1;
	}
	return 0;
}


vector<bool> booleanlatticerestrictednonredundant(int n, int k, int &r)
{
//Devuelve la matriz del retículo booleano Bn
//	int_64 s = (int_64)1<<n;
	// count the number of items as the sum of Cin
	sizeindependent(n, k, r);
	
	vector<bool> B(r*r,false);
	
	for(int i=1;i<r;i++)
	{
		int_64 a = card2bit[i];
		for(int j=0;j<i;j++){
			int_64 b = card2bit[j];
			
			if(IsSubset(a,b) && (preceedsby1(a,b)==-1)) B[r*i+j]=true;
		}
		
	}
//	for(int i=0;i<r;i++)B[r*i+i]=true;  // reflexive
	
//	for(int i=0;i<r;i++)B[r*(r-1)+i]=true;  // top element
	

	return B;
}	




int convertintomatrix(vector<bool> & B, vector<float> & A, vector<float> &b, vector<int>& dir, int r)
{
	int t=0;
	for(int i=0;i<r;i++){
		for(int j=0; j<r;j++)
			if(B[r*i+j]){
				b[t]=0;
				dir[t]=2; // <=
				A[t*r+i]=-1;
				A[t*r+j]=1;
				t++;
			}
	}
	return t;
	
}
