#ifndef GUARD_unwrap_pstate_h
#define GUARD_unwrap_pstate_h

#include "Momentum_Basis.h"
#include "read_inputfile.h"

/////////////////////////////////////////////////////////////
//
// Modified updated version that parallelizes unwrapping ////
void Unwrap_pState(	vector< complex <double> > &v_coeff,
			vector<int64_t> &pBasis,
			vector<double> &Norm,
			double &px,
			double &py,
			vector<int> &T1,
			vector<int> &T2,
			vector<ids> &psi,
			int n_site,
			int n_ones,
			ofstream &outfile);


void read_n_unwrap(	std::string filein,
			std::string filepEV,
			vector<int> &T1,
			vector<int> &T2,
			vector<ids> &psi,
			ofstream &outfile);

#endif


