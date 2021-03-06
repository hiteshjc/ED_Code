#ifndef GUARD_ARPACK_SECTOR_H
#define GUARD_ARPACK_SECTOR_H

#include "global.h"
#include "ARPACK_Functions.h"
#include "LAPACK_matrix_functions.h"
#include "H_v.h"  
#include "pClasses.h"
#include "Lin_Tables.h"


void ARPACK_Iteration_EV(	vector<int64_t> &pblock_states,
				vector<double> &norm,
				vector<double> &norm_inv,
				vector<int64_t> &rep_loc, 
				vector< complex<double> > &state_ph,
				int &n_spins,
				int &n_it,
		      		double &px,
				double &py,
		      		vector < complex<double> > &Evals,
				zMatrix &Evecs,
		      		bool find_ev,
				double lambda,
				vector<coordinates> &adj_list,
				double &Ch,
				vector<triangles> &ch_list,
				double &ChCh,
				vector<bowties> &chch_list,
		      		double tol,
		      		int maxiter,
				vector<int> &T1,
				vector<int> &T2, 
				vector<int64_t> &Ia,
				vector<int64_t> &Ib,
				int &bits_right,
				ofstream &outfile);

void ARPACK_sector(	vector<int64_t> &pblock_states,
			vector<double> &norm,
			vector<double> &norm_inv,
			vector<int64_t> &rep_loc,
			vector< complex<double> > &state_ph,
			vector<complex<double> > &pblock_eigs,
			zMatrix &pblock_eigvecs,
			bool find_ev,
			double &lambda,
			int &n_spins,
			double &px, double &py,
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
			int maxiter,
			vector<int64_t> &Ia,
			vector<int64_t> &Ib,
			int &bits_right,
			ofstream &outfile);

void sector_eigs_eigvecs(	int &kx,
				int &ky,
				double &cSz,
				int &n_spins,
				double &lambda,
				double &J2,
				vector<coordinates> &adj_list,
				double &Ch,
				vector<triangles> &ch_list,
				double &ChCh,
				vector<bowties> &chch_list,
				vector<int> &T1,
				vector<int> &T2,
				int &n_sec_ev,
				bool &find_ev,
				bool &unwrap_peigvec,
				string &filepath_EV,
				double &tol,
				int &maxiter,
				ofstream &outfile);
#endif
