#ifndef GUARD_H_V_H
#define GUARD_H_V_H

#include "global.h"
#include "Heisenberg_Hamiltonian.h"

using namespace std;

void makeH(	const vector<int64_t> &pblock_states,
		const vector<double> &norm,
		const vector<double> &norm_inv,
		const vector<int64_t> &rep_loc,
	        const vector< complex<double> > &state_ph,
		const double &lambda,
		const vector<coordinates> &adj_list,
		const double &Ch,
		const vector<triangles> &ch_list,
		const double &ChCh,
		const vector<bowties> &chch_list,
		const int &n_sites,
		const vector<int64_t> &Ia,
		const vector<int64_t> &Ib,
		const int &bits_right,
		const vector<int> &T1,
		const vector<int> &T2,
		const double &kx,
		const double &ky,
		ofstream &outfile, zMatrix &H);

// H*v for ARPACK calls (with OMP parallelization) with sparse Hamiltonian
void sHv(	sMatrix &pBij, 
		complex< double > *v,
		complex< double > *w);

// H*v for ARPACK calls (Works for OMP parallelization) on the fly
void sHv(	complex< double > *v,
		complex< double > *w,
		const vector<int64_t> &pblock_states,
		const vector<double> &norm,
		const vector<double> &norm_inv,
		const vector<int64_t> &rep_loc,
	        const vector< complex<double> > &state_ph,
		const double &lambda,
		const vector<coordinates> &adj_list,
		const double &Ch,
		const vector<triangles> &ch_list,
		const double &ChCh,
		const vector<bowties> &chch_list,
		const int &n_sites,
		const vector<int64_t> &Ia,
		const vector<int64_t> &Ib,
		const int &bits_right,
		const vector<int> &T1,
		const vector<int> &T2,
		const double &kx,
		const double &ky,
		ofstream &outfile);

void sHv(	complex< double > *v,
		complex< double > *w,
		const vector<int64_t> &pblock_states,
		const vector<double> &norm,
		const vector<double> &norm_inv,
		const vector<int64_t> &rep_loc,
	        const vector< complex<double> > &state_ph,
		const double &lambda,
		const vector<coordinates> &adj_list,
		const double &Ch,
		const vector<triangles> &ch_list,
		const double &ChCh,
		const vector<bowties> &chch_list,
		const int &n_sites,
		const vector<int64_t> &Ia,
		const vector<int64_t> &Ib,
		const int &bits_right,
		const vector<int> &T1,
		const vector<int> &T2,
		const vector<int> &T3,
		const double &kx,
		const double &ky,
		const double &kz,
		ofstream &outfile);
#endif
