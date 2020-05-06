#ifndef GUARD_IDSCLASS_H
#define GUARD_IDSCLASS_H

#include "global.h"


// Structure that holds the State Basis label and its coefficient
struct ids{
	int64_t stateID;
	complex< double > coeff;  
};


// Function to add an element to vector<ids>

void add2w(	vector<ids> &w,
	   	int64_t j,
		complex< double > coeff);

// Function that adds an element to vector<ids>
// Used to unwrap a momentum state back to its basis elements
void add2ids(	vector<ids> &w,
	   	int64_t j,
		complex< double > coeff);


// Function to add two vector<ids>

vector<ids> add_ids(	vector<ids> &psi1,
			vector<ids> &psi2,
			complex<double> &c1,
			complex<double> &c2);

// Dot product of two vector<ids>

complex<double> dot_product_ids(	vector<ids> &psi1,
					vector<ids> &psi2);

/*
void swap_ids(	vector<ids> &psi,
		int64_t i,
		int64_t j);

void merge_ids(	vector<ids> &psi,
		int64_t st,
		int64_t mid,
		int64_t end);


void msort_ids(	vector<ids> &psi, 
		int64_t &st,
		int64_t &end);
*/
int64_t find_ids_state(	vector<ids> &psi,
			int64_t &id);


#endif
