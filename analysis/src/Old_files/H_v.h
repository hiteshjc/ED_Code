#ifndef GUARD_H_V_H
#define GUARD_H_V_H

#include "global.h"
#include "Heisenberg_Hamiltonian.h"

using namespace std;

// H*v for ARPACK calls (with OMP parallelization) with sparse Hamiltonian
void sHv(	sMatrix &pBij, 
		complex< double > *v,
		complex< double > *w);

// H*v for ARPACK calls (Works for OMP parallelization) on the fly
void sHv(	complex< double > *v,
		complex< double > *w,
		vector<int64_t> &pblock_states,
		vector<double> &Norm,
		double &kx,
		double &ky,
		double &lambda,
		vector<coordinates> &adj_list,
		vector<int> &T1,
		vector<int> &T2,
		ofstream &outfile);


// H*v DONT USE FOR OMP
/*
void sHv(	complex< double > *v,
		complex< double > *w,
		vector<int64_t> &pblock_states,
		vector<double> &Norm,
		double &kx,
		double &ky,
		double &lambda,
		vector<coordinates> &adj_list,
		vector<int> &T1,
		vector<int> &T2,
		ofstream &outfile);

*/

#endif
