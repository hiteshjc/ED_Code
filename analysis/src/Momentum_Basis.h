#ifndef GUARD_Momentum_Basis_h
#define GUARD_Momentum_Basis_h

#include "global.h"
#include "Matrix_Functions.h"
#include "Basis_States.h"
#include "sort.h"


using namespace std;

// Function to translates Basis state
void translateT(	int64_t &BasisID,
			vector<int> &T);

// Function that determines the maximum translations along two directions L1 & L2
void LxLy(	int &L1, int &L2,
	  	vector<int> &T1,
		vector<int> &T2);

// Funtion that checks if a basis state is the representative
bool ifnotrepresentative(	int64_t &v_id,
				int &L1, int &L2,
			 	vector<int> &T1,
				vector<int> &T2);

// Funtion that checks if particular Basis state was already visited
bool visited(	vector<int64_t> &visit,
		int64_t &tmp_id);

int64_t findstate_basis(	int64_t &rep_id,
				vector<int64_t> &basis);

// Function that constructs the momentum basis states
// Function outputs the state's representative id and its normalization

void Construct_States(	vector<int64_t> &v_sz,
			vector<double> &norm,
		      	vector<int> &T1,
			vector<int> &T2,
		      	int &num_ones,
			double &px,
			double &py,
			int64_t &n_states,
			ofstream &outfile);

#endif
