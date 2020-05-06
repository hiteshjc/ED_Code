#ifndef GUARD_ARPACK_FUNCTIONS_H
#define GUARD_ARPACK_FUNCTIONS_H

#include "global.h"

using namespace std;

extern "C" void znaupd_(int *ido, char *bmat, int64_t *n, char *which,
			int *nev, double *tol, complex<double> *resid,
			int *ncv, complex<double> *v, int *ldv,
			int *iparam, int *ipntr, complex<double> *workd,
			complex<double> *workl, int *lworkl,
			double *rwork, int *info);

extern "C" void zneupd_(int *rvec, char *All, int *select,
			complex<double> *d, complex<double> *z, int *ldz,
			double *sigma, complex<double> *workev, char *bmat,
			int64_t *n, char *which, int *nev, double *tol,
			complex<double> *resid, int *ncv,
			complex<double> *v, int *ldv, int *iparam,
			int *ipntr, complex<double> *workd,
			complex<double> *workl, int *lworkl,
			double *rwork, int *ierr);



#endif
