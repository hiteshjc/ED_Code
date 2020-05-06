#ifndef GUARD_READ_VECTOR_H
#define GUARD_READ_VECTOR_H

#include "global.h"
#include "idsClass.h"

void read_pEV(	int64_t &n_p,
		double &px,
		double &py,
		vector<int64_t> &pblock_states,
		vector<double> &norm,
		vector< complex <double> > &v_coeff, 
		string &filepEV);

void readEV(	vector<ids> &psi,
		string &fileEV);



#endif
