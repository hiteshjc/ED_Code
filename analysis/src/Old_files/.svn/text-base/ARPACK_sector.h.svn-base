#ifndef GUARD_ARPACK_SECTOR_H
#define GUARD_ARPACK_SECTOR_H

#include "global.h"
#include "ARPACK_Functions.h"
#include "LAPACK_matrix_functions.h"
#include "H_v.h"  
#include "pClasses.h"

void ARPACK_Iteration_EV(	vector<int64_t> &pblock_states,
				vector<double> &norm,
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
		      		double tol,
		      		int maxiter,
				vector<int> &T1,
				vector<int> &T2, 
		      		ofstream &outfile);

void ARPACK_sector(	vector<int64_t> &pblock_states,
			vector<double> &Norm,
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
			vector<int> &T1,
			vector<int> &T2,
			int n_ev,
			double tol,
			int maxiter,
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
