#ifndef GUARD_CHIRALITY_EXP_H
#define GUARD_CHIRALITY_EXP_H


#include "global.h"
#include "idsClass.h"
#include "Basis_States.h"
#include "matrix.h"
#include "Momentum_Basis.h"
complex<double> sxsxsysy_op(	int &i,
				int &j,
				int64_t &n1,
				vector<ids> &psi1,
				vector<ids> &psi2,
				ofstream &out);
complex<double> szsz_op(	int &i,
				int &j,
				int64_t &n1,
				vector<ids> &psi1,
				vector<ids> &psi2,
				ofstream &out);

void act_chirality(int i, int j, int k, int64_t &id, vector<int64_t> &new_ids, vector<complex<double> > &hints);

complex<double> ch_ch_op(	int &i,
				int &j,
				int &k,
				int &p,
				int &q,
				int &r,
				int64_t &n1,
				vector<ids> &psi1,
				vector<ids> &psi2,
				ofstream &out);

complex<double> SipSjm(	int &i,
				int &j,
				int64_t &n1,
				vector<ids> &psi1,
				vector<ids> &psi2,
				ofstream &out);

complex<double> ch_op(	int &i,
			int &j,
			int &k,
			int64_t &n1,
			vector<ids> &psi1,
			vector<ids> &psi2,
			ofstream &out);

complex<double> curr_op(	int &i,
				int &j,
				int64_t &n1,
				vector<ids> &psi1,
				vector<ids> &psi2,
				ofstream &out);

#endif
