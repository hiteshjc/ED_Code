#ifndef GUARD_Heisenberg_Hamiltonian_h
#define GUARD_Heisenberg_Hamiltonian_h

// Function to compute the matrix elements of the XXZ Heisenberg model from the Adjacency list

#include "global.h"
#include "Momentum_Basis.h"

// Function to find the representative given a Basis stateID
// Function also returns the number of hops (lx, ly) need to get to he representative state
int64_t representative(	int &lx,
			int &ly,
			int &L1,
			int &L2,
			int64_t &v_id,
			vector<int> &T1,
			vector<int> &T2);

// Locates the representative pBasis state using binary search
int64_t findstate(	int64_t &rep_id,
	       		vector<int64_t> &pBasis);

void findLxLy(	int &Lx,
		int &Ly,
		vector<int> &T1,
	      	vector<int> &T2);

void zXXZ_Heisenberg(	zMatrix &pBij,
			int &n_spins, 
			vector<int64_t> &pblock_states,
			vector<double> &norm,
		   	double &kx,
			double &ky,
			double &lambda,
			vector<coordinates> &adj_list,
			double &Ch,
			vector<triangles> &ch_list,
			vector<int> &T1,
			vector<int> &T2);

void sXXZ_Heisenberg(	sMatrix &pBij,
			int &n_spins, 
			vector<int64_t> &pblock_states,
			vector<double> &norm,
			double &kx,
			double &ky,
			double &lambda,
			vector<coordinates> &adj_list,
			double &Ch,
			vector<triangles> &ch_list,
			vector<int> &T1,
			vector<int> &T2,
			ofstream &outfile);
#endif

