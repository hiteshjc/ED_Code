#ifndef GUARD_LAPACK_matrix_functions_h
#define GUARD_LAPACK_matrix_functions_h

#include<complex>

////////////////////////////////////////////////////////////////
// LAPACK function for matrix multiplication: C = alpha*A(m,k)*B(k,n) + beta*C(m,n) 

extern "C" void dgemm_(const char *transa,	// 'N'
		       const char *transb,	// 'N'
		       const int *m,		// m = A.NRows(); or C.NRows();
		       const int *n,		// n = B.NCols();
		       const int *k,		// k = A.NCols(); or B.NRows();
		       const double *alpha,
		       double *A,		
		       const int *lda,		// m = Dim of A 
		       double *B,
		       const int *ldb,		// k = Dim of B
		       const double *beta,
		       double *C,		// m = Dim of C
		       const int *ldc); 	

inline void dgemm(const char transa,	// 'N'
		       const char transb,	// 'N'
		       const int m,		// m = A.NRows(); or C.NRows();
		       const int n,		// n = B.NCols();
		       const int k,		// k = A.NCols(); or B.NRows();
		       const double alpha,
		       double *A,		
		       const int lda,		// m = Dim of A 
		       double *B,
		       const int ldb,		// k = Dim of B
		       const double beta,
		       double *C,		// m = Dim of C
		       const int ldc)	
{
  dgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

extern "C" void zgemm_(const char *transa,
		      const char *transb,
		      const int *m,
		      const int *n,
		      const int *k,
		      complex<double> *alpha,
		      complex<double> *A,
		      const int *lda,
		      complex<double> *B,
		      const int *ldb, 
		      complex< double > *beta,
		      complex< double > *C,
		      const int *ldc);


inline void zgemm(const char transa,		// 'N'
		  const char transb,		// 'N'
		  const int m,			// m = # rows in A/C
		  const int n,			// n = # cols in B
		  const int k,			// k = # cols in A = # rows in B
		  complex<double> alpha,	// scalar multiplier
		  complex<double> *A,		// complex array of dim (lda, k)
		  const int lda,		// dim of A = m
		  complex<double> *B,		// complex array of dim (ldb, k)
		  const int ldb,		// dim of k = n
		  complex< double > beta,	// additive const
		  complex< double > *C,		// final vector
		  const int ldc)		// dim of C
{
  zgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}


///////////////////////////////////////////////////////////////

// LAPACK function to perform Y = alpha*A*X + beta*Y
// Scalars: alpha, beta
// Vectors: X, Y
// Matrix: A

extern "C" void dgemv_(const char *TRANS,
		   const int *M,
		   const int *N,
		   const double *ALPHA,
		   double *A,
		   const int *LDA,
		   double *X,
		   const int *INCX,
		   const double *BETA,
		   double *Y,
		   const int *INCY);

inline void dgemv(const char TRANS,
		   const int M,
		   const int N,
		   const double ALPHA,
		   double *A,
		   const int LDA,
		   double *X,
		   const int INCX,
		   const double BETA,
		   double *Y,
		   const int INCY)	
{
    dgemv_(&TRANS, &M, &N, &ALPHA, A, &LDA, X, &INCX, &BETA, Y, &INCY);
}

// INCX & INCY are both 1 and correspond the the way c++ vs fortran labels vectors a(0)
//corresponds to a(1) in the other.

/////////////////////////////////////////////////////////////

// LAPACK function for adding vectors: Y = a*X + Y;

extern "C" void daxpy_(const int *n,
		     const double *a,
		     double *X,
		     const int *incx,
		     double *Y,
		     const int *incy);

inline void dapxy(const int n,		// Size of vector
		 const double a,
		 double *X,
		 const int incx,
		 double *Y,
		 const int incy)
{
  daxpy_(&n, &a, X, &incx, Y, &incy);
}

extern "C" void zaxpy_(const int *n,
		       complex< double > *a,
		       complex< double > *X,
		       const int *incx,
		       complex< double > *Y,
		       const int *incy); 	

inline void zaxpy(const int n,
		  complex< double > a,
		  complex< double > *X,
		  const int incx,
		  complex< double > *Y,
		  const int incy)
{
  zaxpy_(&n, &a, X, &incx, Y, &incy);
}


///////////////////////////////////////////////////////////////

// LAPACK function for diagonalizing and finding out eigenvalues and eigenvectors

extern "C" void dgeev_(const char *jobvl,
		       const char *jobvr,
		       const int *n,
		       double *A,
		       const int *lda,
		       double *WR,
		       double *WL,
		       double *VL,
		       const int *ldvl,
		       double *VR,
		       const int *ldvr,
		       double *WORK,
		       const int *lwork,
		       const int *info);

inline void dgeev(const char jobvl,	// 'N': Don't compute left eigenvectors		
		  const char jobvr,	// 'V': Compute right eigenvectors
		  const int n, 		// Order of matrix A
		  double *A, 		// A(lda, n), A is overwritten
		  const int lda,	// n
		  double *WR, 		// n
		  double *WI,		// n
		  double *VL,		// Left eigenvector...not used
		  const int ldvl,	// n
		  double *VR, 		// Compute right eigenvectors
		  const int ldvr,	// n
		  double *WORK,		// array of size 4*n
		  const int lwork,	// >= 4*n
		  const int info)	// if info = 0 implies successful exit
{
    dgeev_(&jobvl, &jobvr, &n, A, &lda, WR, WI, VL, &ldvl, VR, &ldvr, WORK, &lwork, &info);
}

////////////////////////////////////////////////////////////

// LAPACK function for diagonalizing complex matrices
extern "C" void zheev_(char *, char *, int *, complex<double> *, int *, double *, complex<double> *, int *,double *, int *);

inline void zheev(char jobz, char uplo, int n, complex<double> *a, int lda,
                double * w, complex<double> * work, int lwork, double *rwork, int & info)
{
    zheev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &info);
}


extern "C" void cgeev_(const char *jobvl,	// 'N': Don't compute left eigenvectors
		       const char *jobvr,	// 'V' : Compute right eigenvectors
		       const int *n,		// order of matrix A
		       complex<double> *A,	// Complex Array A(lda, n)
		       const int *lda,		// Dim of array A
		       complex<double> *W,	// Contains computed eigenvalues
		       complex<double> *VL,	// left eigenvectors...not used
		       const int *ldvl,		// dimension of VL
		       complex<double> *VR,	// Right complex array (ldvr, n)
		       const int *ldvr,		// Dimension of VR
		       complex<double> *WORK,	// Complex array of dim lwork
		       const int *lwork,		// lwork must be atleast 2*n
		       double *RWORK,
		       const int *info);		// 0 for successful exit

inline void cgeev(const char jobvl,
		       const char jobvr,
		       const int n,
		       complex<double> *A,
		       const int lda,
		       complex<double> *W,
		       complex<double> *VL,
		       const int ldvl,
		       complex<double> *VR,
		       const int ldvr,
		       complex<double> *WORK,
		       const int lwork,
		       double *RWORK,
		       const int info)
{
	cgeev_(&jobvl, &jobvr, &n, A, &lda, W, VL, &ldvl, VR, &ldvr, WORK, &lwork, RWORK, &info);
}

extern "C" void zgeev_(const char *jobvl,
		       const char *jobvr,
		       const int *n,
		       complex< double > *A,
		       const int *lda,
		       complex< double > *W,
		       complex< double > *VL,
		       const int *ldvl,
		       complex< double > *VR,
		       const int *ldvr,
		       complex< double > *WORK,
		       const int *lwork,
		       double *RWORK,
		       const int *info);

inline void zgeev(const char jobvl,		// 'N'
		  const char jobvr,		// 'V'
		  const int n,			// order of A
		  complex< double > *A,		// A complex array dim(lda, n)
		  const int lda,		// n
		  complex< double > *W,		// contains computed eigenvalues
		  complex< double > *VL,	// left eigenvectors
		  const int ldvl,		// n
		  complex< double > *VR,	// right eigenvectors
		  const int ldvr,		// n
		  complex< double > *WORK,	// work array of dim(1, lwork)
		  const int lwork,		// lwork of dimension 2*n
		  double *RWORK,		// double precision array of size 2*n
		  const int info)		// successful exit info = 0;
{
	zgeev_(&jobvl, &jobvr, &n, A, &lda, W, VL, &ldvl, VR, &ldvr, WORK, &lwork, RWORK, &info);
}
#endif
