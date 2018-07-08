#ifndef GUARD_MATRIX_FUNCTIONS_H
#define GUARD_MATRIX_FUNCTIONS_H

#include "global.h"
#include "LAPACK_matrix_functions.h"
#include "matrix.h"

void Matrix_Diagonalize(	Matrix &a,
				vector< complex<double> > &eigs,
				Matrix &eigenvecs);

void zMatrix_Diagonalize(	zMatrix &zM,
				vector< complex<double> > &W,
				zMatrix &VR);

void zMatrix_Diagonalize_Hermitian(	zMatrix &zM,
				vector< double > &W,
				zMatrix &VR);

complex< double > zDot_Product(	complex<double> alpha, 
				vector< complex< double > > &A,
		    		vector< complex< double > > &B);

void zMatrix_Vector_Product(	zMatrix M,
			    	vector< complex< double > > v,
			    	vector< complex< double > > &w);

void zVector_Sum(	complex< double > a,
			vector< complex< double > > X,
			vector< complex< double > > &Y);

#endif
