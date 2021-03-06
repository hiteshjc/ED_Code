#ifndef GUARD_Momentum_Basis_h
#define GUARD_Momentum_Basis_h

#include "global.h"
#include "Matrix_Functions.h"
#include "Basis_States.h"
#include "sort.h"
#include "Lin_Tables.h"

using namespace std;

// Function to translates Basis state
void translateT(	int64_t &BasisID,
			const vector<int> &T);

// Function that determines the maximum translations along two directions L1 & L2
void LxLy(	int &L1, int &L2,
	  	vector<int> &T1,
		vector<int> &T2);

void LxLyLz(	int &L1, int &L2, int &L3,
	  	vector<int> &T1,
		vector<int> &T2,
		vector<int> &T3);


// Funtion that checks if a basis state is the representative
bool ifnotrepresentative(	int64_t &v_id,
				int &L1, int &L2,
			 	vector<int> &T1,
				vector<int> &T2);
// Funtion that checks if particular Basis state was already visited
bool visited(	vector<int64_t> &visit,
		int64_t &tmp_id);


void Construct_States(	vector<int64_t> &v_sz,
			vector<double> &norm,
		      	vector<int> &T1,
			vector<int> &T2,
		      	int &num_ones,
			double &px,
			double &py,
			int64_t &n_states,
			ofstream &outfile);

// Function that constructs the momentum basis states
// Function also stores all the representative ids

void Construct_States(	vector<int64_t> &v_sz,
			vector<double> &norm,
			vector<double> &norm_inv,
			vector<int> &T1,
			vector<int> &T2,
		      	int &num_ones,
			double &px,
			double &py,
			int64_t &n_states,
			vector<int64_t> &rep_loc,
			vector< complex<double> > &state_ph,
			vector<int64_t> &Ia,
			vector<int64_t> &Ib,
			int &bits_right,
			ofstream &outfile);



void pConstruct_States(	vector<int64_t> &v_sz,
			vector<double> &norm,
			vector<double> &norm_inv,
			vector<int> &T1,
			vector<int> &T2,
		      	int &num_ones,
			double &px,
			double &py,
			int64_t &n_states,
			vector<int64_t> &rep_loc,
			vector< complex<double> > &state_ph,
			vector<int64_t> &Ia,
			vector<int64_t> &Ib,
			int &bits_right,
			ofstream &outfile);

void pConstruct_States(	vector<int64_t> &v_sz,
			vector<double> &norm,
		      	vector<double> &norm_inv,
			vector<int> &T1,
			vector<int> &T2,
			vector<int> &T3,
		      	int &num_ones,
			double &px,
			double &py,
			double &pz,
			int64_t &n_states,
			vector<int64_t> &rep_loc, 
		        vector< complex<double> > &state_ph,
			vector<int64_t> &Ia,
			vector<int64_t> &Ib,
			int &bits_right,
			ofstream &outfile);


// Function to Unwrap a pState back to its original Basis states
vector< ids > Unwrap_pState(	vector< complex <double> > &v_coeff,
				vector<int64_t> &pBasis,
				vector<double> &Norm,
				double &px,
				double &py,
				vector<int> &T1,
				vector<int> &T2);

vector< ids > Unwrap_pState(	vector< complex <double> > &v_coeff,
				vector<int64_t> &pBasis,
				vector<double> &Norm,
				double &px,
				double &py,
				double &pz,
				vector<int> &T1,
				vector<int> &T2,
				vector<int> &T3);


#endif
