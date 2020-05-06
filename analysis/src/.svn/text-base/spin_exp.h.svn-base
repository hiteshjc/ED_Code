#ifndef GUARD_spin_exp_h
#define GUARD_spin_exp_h

#include "global.h"
#include "idsClass.h"
#include "Basis_States.h"
#include "matrix.h"
#include "Momentum_Basis.h"


int64_t findstate(	int64_t &rep_id,
	       		vector<ids> &Basis);

vector< complex<double> > spin_spin_exp(	int &i, 
						int &j, 
						int64_t &n1, 
						vector<ids> &psi,
						ofstream &outfile);



vector< complex<double> > bond_bond_exp(	links &i,
						links &j,
						int64_t &n1,
						vector<ids> &psi,
						ofstream &out);


#endif
