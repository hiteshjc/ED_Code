#ifndef GUARD_ARPACK_ITERATION_H
#define GUARD_ARPACK_ITERATION_H

#include "global.h"
#include "matrix.h"
#include "pClasses.h"
#include "ARPACK_Functions.h"
#include "LAPACK_matrix_functions.h"
#include "Momentum_Basis.h"
#include "H_v.h"

void ARPACK_Iteration(	vector<int64_t> &pblock_states,
			vector<double> &norm,
			int &n_spins,
			int &n_it,
	//		int64_t &n,
			double &px,
			double &py,
			vector < complex<double> > &Evals,
		      	double lambda,
			vector<coordinates> &adj_list,
			double &Ch,
			vector<triangles> &ch_list,
			double tol,
		      	int maxiter,
			vector<int> &T1,
			vector<int> &T2, 
			ofstream &outfile);

void ARPACK(	pVector &LanczosEij,
		double &lambda,
		int &n_spins,
		double &sector_Sz,
		vector<coordinates> &adj_list,
		double &Ch,
		vector<triangles> &ch_list,
		vector<int> &T1,
		vector<int> &T2,
		int n_ev,
		double tol,
		int maxiter,
		ofstream &outfile);
#endif





