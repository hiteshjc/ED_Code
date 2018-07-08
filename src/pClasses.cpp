#include "pClasses.h"
// Vector for storing Quantum numbers and energy eigenvalues
/*
class pVector{
private:
	struct nQuantum
	{
		size_t id;
		double En, px, py, Sz;  
	};
	vector<nQuantum> V;
	size_t n;

	void merge(vector<nQuantum> &E, vector<nQuantum> E1, size_t n1, vector<nQuantum> E2, size_t n2);
	void merge_sort_En(vector<nQuantum> &E, size_t m);
	void sort_by_Sz(vector<nQuantum> &E);
	void swapVij(vector<nQuantum> &E, size_t i, size_t j);


public:
	pVector() : n(0) {}
	pVector(size_t n_) { V.resize(n_);	n = n_;	}

	void resize(size_t n_)
	{ n = n_; V.resize(n);}

	size_t size()	{return n;}

	void push_back(size_t i, size_t idd, double tEn, double tpx, double tpy, double tSz){
		if(i < n){
			V[i].id = idd;	V[i].En = tEn;	V[i].px = tpx;	V[i].py = tpy;	V[i].Sz = tSz;
		}
		else
			cout << "\nError in pVector push_back \n";
	}

	void push_back(size_t idd, double tEn, double tpx, double tpy, double tSz){

		nQuantum ntmp;
		ntmp.id = idd; ntmp.En = tEn; ntmp.px = tpx; ntmp.py = tpy; ntmp.Sz = tSz;
		V.push_back(ntmp);
		n++;
	}

	void update_id(size_t i, size_t ni){
		if(i < n)
			V[i].id = ni;
		else
			cout << "\n Error in update_id in pVector \n";
	}

	void update_p(size_t id, double px, double py){
		if(id < n){
			V[id].px = px;	V[id].py = py;
		}
		else
			cout << "\n Error in update_p in pVector \n";
	}

	void update_Sz(size_t id, double Sz){
		if(id < n)
			V[id].Sz = Sz;
		else
			cout << "\n Error in update_Sz in pVector \n";
	}

	void update_En(size_t id, double En){
		if(id < n)
			V[id].En = En;    
		else
			cout << "\n Error in update_En in pVector \n";
	}

	void state_p(size_t id, double &px, double &py){
		if(id < n){
			px = V[id].px;	py = V[id].py;
		}
		else
			cout << "\n Error in state_p in pVector \n";
	}

	double state_Sz(size_t id){	return V[id].Sz;}
	double state_En(size_t id){	return V[id].En;}
	size_t state_id(size_t i){	return V[i].id;	}

	void swap(size_t i, size_t j){
		nQuantum tmpqn = V[i];
		V[i] = V[j];
		V[j] = tmpqn;
	}

	void outfile(ofstream &outfile){
		outfile << "\n\nID\tEigenvalues\t\tSz\t(px, py) \n\n";
		for(size_t i = 0; i < n; i++)
			outfile << V[i].id << "\t" << fixed << setprecision(15) << V[i].En << setprecision(3)
		<< "\t" << V[i].Sz << "\t(" << V[i].px << "," << V[i].py << ")" << endl;
	}

	void outfile(	ofstream &outfile,
			int &n_spins,
			double &Jz,
			double &J2)
	{
		outfile << "\n\nID\tEnergies\tEnergies/site\tSz\tpx\tpy\tJz\tJ2\tEnergies-E0\n\n";
		for(size_t i = 0; i < n; i++)
			outfile << V[i].id << "\t" << fixed << setprecision(10) << V[i].En << "\t" << V[i].En/double(n_spins) << setprecision(3)
		<< "\t" << V[i].Sz << "\t" << V[i].px << "\t" << V[i].py << "\t" << Jz << "\t" << J2 << "\t" << setprecision(10) << V[i].En-V[0].En<< endl;
	}

	
};

void pVector::sort_by_En()
  {	pVector::merge_sort_En(V, n);	pVector::sort_by_Sz(V);	}
*/
void pVector::swapVij(	vector<nQuantum>& E,
			int64_t i,
			int64_t j)
{
	nQuantum tmp = E[i];
	E[i] = E[j];
	E[j] = tmp;
}
 
void pVector::merge(	vector<nQuantum>& E,
			vector<nQuantum> E1,
			int64_t n1,
			vector<nQuantum> E2,
			int64_t n2)
{
  for (int i = 0, i1 = 0, i2 = 0; i < n1 + n2; i++)
  {
    if (i1 < n1 && i2 < n2)
    {
      if (E1[i1].En > E2[i2].En)
	E[i] = E2[i2++];
      else if(E1[i1].En == E2[i2].En)
      {
	if(E1[i1].Sz > E2[i2].Sz)
	  E[i] = E2[i2++];
	else
	  E[i] = E1[i1++];
      }
      else
	E[i] = E1[i1++];
    }
    else if (i1 < n1)
      E[i] = E1[i1++];
    else if (i2 < n2)
      E[i] = E2[i2++];	
    
  }
}

void pVector::merge_sort_En(	vector<nQuantum> &E,
				int64_t m)
{
  int64_t n1 = m/2, n2;
  if (m % 2 == 0)
    n2 = m/2;
  else
    n2 = m/2 + 1;
    
  vector<nQuantum> E1(n1), E2(n2);
  
  for (int64_t i = 0; i < n1; i++)
    E1[i] = E[i];
  for(int64_t i = 0; i < n2; i++)
    E2[i] = E[i + n1]; 

  if (n1 > 2)
    merge_sort_En(E1, n1);
  else if(n1 == 2)		// Base case with just two elements
    if (E1[0].En > E1[1].En)
      swapVij(E1,0,1);
    
  if (n2 > 2)
    merge_sort_En(E2, n2);
  else if (n2 == 2)		// Base case with just two elements
    if (E2[0].En > E2[1].En)
      swapVij(E2, 0, 1);
    
  merge(E, E1, n1, E2, n2);
    
  E1.clear();
  E2.clear();
}


void pVector::sort_by_Sz(	vector< pVector::nQuantum >& E)
{ 
  for(int64_t i = 0; i < E.size(); i++)
  {
    double tE = E[i].En;
    double tSz = E[i].Sz;
    int64_t j = i, k = i;
    while(abs(E[j].En-tE) < pres)
    {
      if(tSz > E[j].Sz)
      {	
	k = j;
	tSz = E[j].Sz;
      }
      j++;
    }    
    swapVij(E, k, i);    
  }
}


