#ifndef GUARD_PCLASSES_H
#define GUARD_PCLASSES_H

#include "global.h"
// Vector for storing Quantum numbers and energy eigenvalues

class pVector{
private:
	struct nQuantum
	{
		int64_t id;
		double En, px, py, Sz;  
	};
	vector<nQuantum> V;
	int64_t n;

	void merge(vector<nQuantum> &E, vector<nQuantum> E1, int64_t n1, vector<nQuantum> E2, int64_t n2);
	void merge_sort_En(vector<nQuantum> &E, int64_t m);
	void sort_by_Sz(vector<nQuantum> &E);
	void swapVij(vector<nQuantum> &E, int64_t i, int64_t j);


public:
	pVector() : n(0) {}
	pVector(int64_t n_) { V.resize(n_);	n = n_;	}

	void resize(int64_t n_)
	{ n = n_; V.resize(n);}

	int64_t size()	{return n;}

	void push_back(int64_t i, int64_t idd, double tEn, double tpx, double tpy, double tSz){
		if(i < n){
			V[i].id = idd;	V[i].En = tEn;	V[i].px = tpx;	V[i].py = tpy;	V[i].Sz = tSz;
		}
		else
			cout << "\nError in pVector push_back \n";
	}

	void push_back(int64_t idd, double tEn, double tpx, double tpy, double tSz){

		nQuantum ntmp;
		ntmp.id = idd; ntmp.En = tEn; ntmp.px = tpx; ntmp.py = tpy; ntmp.Sz = tSz;
		V.push_back(ntmp);
		n++;
	}

	void update_id(int64_t i, int64_t ni){
		if(i < n)
			V[i].id = ni;
		else
			cout << "\n Error in update_id in pVector \n";
	}

	void update_p(int64_t id, double px, double py){
		if(id < n){
			V[id].px = px;	V[id].py = py;
		}
		else
			cout << "\n Error in update_p in pVector \n";
	}

	void update_Sz(int64_t id, double Sz){
		if(id < n)
			V[id].Sz = Sz;
		else
			cout << "\n Error in update_Sz in pVector \n";
	}

	void update_En(int64_t id, double En){
		if(id < n)
			V[id].En = En;    
		else
			cout << "\n Error in update_En in pVector \n";
	}

	void state_p(int64_t id, double &px, double &py){
		if(id < n){
			px = V[id].px;	py = V[id].py;
		}
		else
			cout << "\n Error in state_p in pVector \n";
	}

	double state_Sz(int64_t id){	return V[id].Sz;}
	double state_En(int64_t id){	return V[id].En;}
	int64_t state_id(int64_t i){	return V[i].id;	}

	void swap(int64_t i, int64_t j){
		nQuantum tmpqn = V[i];
		V[i] = V[j];
		V[j] = tmpqn;
	}

	void outfile(ofstream &outfile){
		outfile << "\n\nID\tEigenvalues\t\tSz\t(px, py) \n\n";
		for(int64_t i = 0; i < n; i++)
			outfile << V[i].id << "\t" << fixed << setprecision(15) << V[i].En << setprecision(3)
		<< "\t" << V[i].Sz << "\t(" << V[i].px << "," << V[i].py << ")" << endl;
	}

	void outfile(	ofstream &outfile,
			int &n_spins,
			double &Jz,
			double &J2)
	{
		outfile << "\n\nID\tEnergies\tEnergies/site\tSz\tpx\tpy\tJz\tJ2\tEnergies-E0\n\n";
		for(int64_t i = 0; i < n; i++)
			outfile << V[i].id << "\t" << fixed << setprecision(10) << V[i].En << "\t" << V[i].En/double(n_spins) << setprecision(3)
		<< "\t" << V[i].Sz << "\t" << V[i].px << "\t" << V[i].py << "\t" << Jz << "\t" << J2 << "\t" << setprecision(10) << V[i].En-V[0].En<< endl;
	}

	void sort_by_En()  {	merge_sort_En(V, n);	sort_by_Sz(V);	}
};
/*
void pVector::swapVij(	vector<nQuantum>& E,
			int64_t i,
			int64_t j);


void pVector::merge(	vector<nQuantum>& E,
			vector<nQuantum> E1,
			int64_t n1,
			vector<nQuantum> E2,
			int64_t n2);

void pVector::merge_sort_En(	vector<nQuantum> &E,
				int64_t m);

void pVector::sort_by_Sz(	vector< pVector::nQuantum >& E);
*/
#endif
