
#ifndef GUARD_krishna_lanczos_h
#define GUARD_krishna_lanczos_h

#include "global.h"
#include "H_v.h"

void lanczos_evecs(	vector<int64_t> &pblock_states,
			vector<double> &Norm,
			vector<double> &Norm_inv,
			vector<int64_t> &rep_loc, 
			vector< complex<double> > &state_ph,
			vector< complex<double> > &pblock_eigs,
			vector< vector< complex<double> > > &pblock_eigvecs,
			bool find_ev,
			double &lambda,
			int &n_spins,
			double &px, 
			double &py,
			double &pz,
			int &n_ones,
			int64_t &n_p,
			vector<coordinates> &adj_list,
			double &Ch,
			vector<triangles> &ch_list,
			double &ChCh,
			vector<bowties> &chch_list,
			vector<int> &T1,
			vector<int> &T2,
			vector<int> &T3,
			int n_ev,
			double tol,
			int iterations,
			int ncycles,
			vector<int64_t> &Ia,
			vector<int64_t> &Ib,
			int bits_right,
			ofstream &outfile);


void lanczos_evecs(	vector<int64_t> &pblock_states,
			vector<double> &Norm,
			vector<double> &Norm_inv,
			vector<int64_t> &rep_loc, 
			vector< complex<double> > &state_ph,
			vector< complex<double> > &pblock_eigs,
			vector< vector< complex<double> > > &pblock_eigvecs,
			bool find_ev,
			double &lambda,
			int &n_spins,
			double &px, 
			double &py,
			int &n_ones,
			int64_t &n_p,
			vector<coordinates> &adj_list,
			double &Ch,
			vector<triangles> &ch_list,
			double &ChCh,
			vector<bowties> &chch_list,
			vector<int> &T1,
			vector<int> &T2,
			int n_ev,
			double tol,
			int iterations,
			int ncycles,
			vector<int64_t> &Ia,
			vector<int64_t> &Ib,
			int bits_right,
			ofstream &outfile);


#endif
