#ifndef Basis_States_h
#define Basis_States_h

// Function to increment a string of bools that represent Spin 1/2 Basis States

#include "global.h"
#include "Simple_Math_Func.h"
#include <bitset>

using namespace std;


// btest64 tells you whether a bit is 0 or 1  at position pos

bool btest64(int64_t a,int pos);

// Takes an integer and sets a bit to 1 at location pos

int64_t ibset64(int64_t a,int pos);

// Takes an integer and sets a bit to 0 at location pos

int64_t ibclr64(int64_t a,int pos);


int64_t flip(int64_t v, int i, int j,
	     bool ti, bool tj);

void constrained_dets(int num_sites,int num_ones,std::vector<int64_t> &dets);

#endif
