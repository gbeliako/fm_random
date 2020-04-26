/*

This code is based on the program supplied by Elias Combarro and Susana Irene Diaz Rodribuez in October 2019, 
which is also published in their paper

"Minimals plus: An improved algorithm for the random generation of linear extensions of partially ordered sets." Information
Sciences, 501:50–67, 2019. 


The program generates random linear extensions of a partial order, in particular Boolean lattices and its sublattices.
In addition, the program keeps track of the linear extensions generated and then can prin their distribution, and also calculates the distance from the uniform.

The code can be optimised by not using the matrix representation of the poset, relying on the relation preceeds. Not attempted here due to lack of time.

Gleb Beliakov, December 2019

*/


vector<bool> booleanlattice(int n)
{
//Devuelve la matriz del retículo booleano Bn
	long int s = 1<<n;
	vector<bool> B(s*s,false);
	if(n==0)
	{	
		B[0]=true;
	}
	else
	{
		long int t = 1<<(n-1);
		vector<bool> A=booleanlattice(n-1);
		for(auto i = 0; i < t; i++)
		{	
			for(auto j = 0; j < t; j++)
			{
				//upperleft
				B[i*s+j]=A[i*t+j];
				//downleft
				B[(i+t)*s+j]=A[i*t+j];
				//downright
				B[(i+t)*s+j+t]=A[i*t+j];
			}
		}
	}
	return B;
}	

vector<vector<myint> > levels(const vector<bool> &M, int n)
{
//Calcula los conjuntos P(k) del poset P de matriz M (empezando a contar en P(0) en lugar de en P(1))
	vector<vector<myint> > L;
	vector<bool> Q(n,true);
	int t = n;
	int h = 0; 
	while(t!=0)
	{
		h++;
		vector<myint> P;
		for(int a = 0; a < n; a++) //for a in Q
		{
			//cout << "a" << endl;
			if(Q[a])
			{
				bool maximal = true;
				for(int b = 0; b < n; b++) //for b in Q
				{
					if(M[a*n+b]&&Q[b]&&a!=b)
					{
						maximal = false;
						break;
					}
				}	
				if(maximal)
				{
					P.push_back(a);
					t--;
					//cout << "t " << t << endl;
				}
			}
		}
		for(unsigned int i = 0; i < P.size(); i++)
			Q[P[i]]=false;
		L.push_back(P);
	}
	return L;
}

vector<vector<myint> > bar(const vector<bool> &M, int n)
{
//Devuelve \bar{a} (el conjunto de elementos que lo cubren) para cada elemento a de P
	vector<vector<myint> > L = levels(M,n);
	uint h = L.size();
	vector<vector<myint> > C(n,vector<myint>(0));
	for(uint k = 0; k < h; k++)
	{
		vector<myint> & Pk = L[k];
		for(uint i = 0; i < Pk.size(); i++)
		{
			myint a = Pk[i];
			vector<myint> A;
			for(int l = k-1; l >=0; l--)
			{
				vector<myint> & Pl = L[l];
				for(int j = 0; j < Pl.size(); j++)
				{
					myint b = Pl[j];
					if(M[a*n+b])
					{
						bool covers = true;
						for(uint m = 0; m < A.size(); m++)
						{
							myint c = A[m];
							if(M[c*n+b])
							{
								covers = false;
								break;
							}
						}
						if(covers)
							A.push_back(b);
					}
				}
			}
			C[a]=A;
		}
	}
	return C;
}	

vector<myint> losw(const vector<bool> & M, int n)
// Devuelve los pesos w_a 
{
	vector<myint> w(n,0);
	vector<vector<myint> > L=levels(M,n);
	uint h = L.size();
	vector<vector<myint> > C=bar(M,n);
	
	for(uint i =0; i < L[0].size(); i++)
	{
		myint a = L[0][i];
		w[a]=1;
	}
	for(uint k = 1; k < h; k++)
	{
		//cout << "k " << k << endl;
		vector<myint> & Pk=L[k];
		//vector<int> & Pkmu=L[k-1];
		vector<myint> Ak;
		vector<myint> Bk;
		for(uint i = 0; i < Pk.size(); i++)
		{
			myint a = Pk[i];	
			bool contenido = true;
			for(uint l = 0; l < k; l++)
			{
				vector<myint> &Pl = L[l];
				for(uint j=0; j < Pl.size(); j++)
				{
					myint b = Pl[j];
					if(!M[a*n+b]) //Hay un elemento de P(k-1) que no están en up(a)
					{
						contenido=false;
						break;
					}
				}
				if(!contenido)
					break;
			}
			if(contenido)
				Ak.push_back(a);
			else
				Bk.push_back(a);
		}
		if(Bk.size()==0)
		{
			for(uint i = 0; i < Pk.size(); i++)
			{
				myint a = Pk[i];
				w[a]=1;
				/*for(int j = 0; j < C[a].size(); j++)
				{
					int x = C[a][j];
					w[a]+=w[x];
				}*/

			}
		}
		else
		{
			for(uint i = 0; i < Bk.size(); i++)
			{
				myint a = Bk[i];
				bool encontrado = false;
				myint ele;
				for(uint l = 0; l < k; l++)
				{
					vector<myint> & Pl = L[l];
					for(uint j = 0; j < Pl.size(); j++)
					{
						myint b = Pl[j];
						//cout << "Probando " << a << " " << b << endl;
						if(!M[a*n+b]) //a y b no comparables
						{
							//Comprobar	si bar(b) <= up(a)
							//cout << "Probando " << a << " " << b << endl;
							bool contenido = true;
							for(uint m = 0; m < C[b].size(); m++)
							{
								myint c = C[b][m];
								if(!M[a*n+c])
								{
									contenido = false;
									break;
								}
							}
							if(contenido)
							{
								encontrado = true;
								ele = b;
								break;
							}
						}
						if(encontrado)
							break;
					}
					if(encontrado)
						break;
				}
				if(!encontrado)
					cout << "ERROR!!!" << a <<endl;
				else
				{
					myint b = ele;
					w[a]+=w[b];
					for(uint j = 0; j < C[a].size(); j++)
					{
						int x = C[a][j];
						if(!M[b*n+x])
							w[a]+=w[x];
					}
				}
			}
			for(uint i = 0; i < Ak.size(); i++)
			{
				myint a = Ak[i];
				bool encontrado = false;
				myint ele;
				for(uint l = 0; l < k; l++)
				{
					vector<myint> & Pl = L[l];
					for(uint j = 0; j < Pl.size(); j++)
					{
						myint b = Pl[j];
						if(!M[a*n+b]) //a y b no comparables
						{
							//Comprobar	si bar(b) <= up(a)
							bool contenido = true;
							for(uint m = 0; m < C[b].size(); m++)
							{
								myint c = C[b][m];
								if(!M[a*n+c])
								{
									contenido = false;
									break;
								}
							}
							if(contenido)
							{
								encontrado = true;
								ele = b;
								break;
							}
						}
						if(encontrado)
							break;
					}
					if(encontrado)
						break;
				}
				if(!encontrado) //Nos vale el primer elemento de Bk				
				{
					ele = Bk[0];
					//cout << "b" << ele << endl;
				}
						
				myint b = ele;
				//cout << "a " << a << " b " << b << endl;
				w[a]+=w[b];
				for(uint j = 0; j < C[a].size(); j++)
				{
					myint x = C[a][j];
					if(!M[b*n+x])
					{
						//cout << "Sumando w" << x << endl; 
						w[a]+=w[x];
					}
				}
			}
		}
	}
	return w;	
}	


vector<myint> markovKKclassic(const vector<bool> & p, int n, const vector<myint> & inicial, int k)
{

	uniform_int_distribution<int> uni(0,n-2); // guaranteed unbiased
	uniform_int_distribution<int> coin(0,1); // guaranteed unbiased
		
	vector<myint> ext = inicial;
	
	for(int i = 0; i < k; i++)
	{
		//Decidimos si movernos o no
		
	//	if(coin(rng))   Not needed
		{
			int pos = uni(rng);	
			myint a = ext[pos];
			myint b = ext[pos+1];
			if(!p[a*n+b])
			{
				ext[pos]=b;
				ext[pos+1]=a;
			}
		}
	} 
	
	return ext;
}


vector<myint> minimals_w(const vector<bool> & p, const vector<myint> & w,int n)
//Devuelve una extension lineal aleatoria usando el algoritmo minimals
//De acuerdo a los pesos w
{

	vector<myint> le;
	vector<bool> q = p;
	vector<myint> down(n,0);
	vector<myint> up(n,0);
	//Calculamos las cotas superiores e inferiores
	for(auto i = 0; i < n; i++)
	{
		for(auto j = 0; j < n; j++)
		{
			if(q[j*n+i])
			{
				down[i]++;
				//up[j++];
			}
			//if(q[i*n+j])
				//up[i]++;
		}

	}
	int faltan = n;
	vector<bool> usado (n,false);
	while(faltan>0)
	{
		vector<myint> boletos;
		for(auto i = 0; i < n; i++)
		{	
			if(down[i]==1&&!usado[i]) //Es un minimal
			{
				for(myint j = 0; j < w[i]; j++)
					boletos.push_back(i);
			}		
		}
		
		uniform_int_distribution<int> uni(0,boletos.size()-1); // guaranteed unbiased	
		myint chosen = boletos[uni(rng)];	


		//eliminamos el elemento
		usado[chosen]=true;
		for(myint j=0; j < n; j++)		
		{
			if(q[chosen*n+j])
				down[j]--;
		}
		//Lo añadimos
		le.push_back(chosen);
		faltan--;
		
	}	
	return le;
}


	







