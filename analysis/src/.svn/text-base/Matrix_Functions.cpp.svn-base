#include "Matrix_Functions.h"


void Matrix_Diagonalize(	Matrix &a,
				vector< complex<double> > &eigs,
				Matrix &eigenvecs)
{
	int n = int(a.NRows());
	int info = 0;

	vector<double> reigs(n), imeigs(n);
	vector<double> work(4*n);
	Matrix vl(n,n), b(n,n);

	b = a;
	eigs.clear();
	eigenvecs.resize(n,n);
	work.assign(4*n, 0.);

	dgeev(	'N', 'V', n, &*b.begin(), n, &*reigs.begin(), &*imeigs.begin(),
		&*vl.begin(), n, &*eigenvecs.begin(), n,
		&*work.begin(), int(work.size()),info);

	for(int i=0; i<n; i++)
		eigs.push_back(complex<double>(reigs[i],imeigs[i]));
}

void zMatrix_Diagonalize(	zMatrix &zM,
				vector< complex<double> > &W,
				zMatrix &VR)
{
	int n = int(zM.NRows());
	int info = 0;

	vector< complex<double> > WORK(2*n);
	zMatrix VL(n,n), A = zM;
	vector<double> RWORK(2*n);

	zgeev(	'N', 'V', n, &*A.begin(), n, &*W.begin(), &*VL.begin(),
		 n, &*VR.begin(), n, &*WORK.begin(), 2*n, &*RWORK.begin(), info);

	if(info != 0)
		cout << "\nerrors in complex matrix diag\n";
}

void zMatrix_Diagonalize_Hermitian(	zMatrix &zM,
					vector< double > &W,
					zMatrix &VR)
{
	int n = int(zM.NRows());
	int info = 0;

	vector< complex<double> > WORK(4*n);
	VR=zM;
	vector<double> RWORK(4*n);

         zheev('V','U',n,&*VR.begin(),n,
           &*W.begin(),&*WORK.begin(),int(WORK.size()),&*RWORK.begin(),info);

	if(info != 0)
		cout << "\nerrors in complex matrix diag\n";
}

complex< double > zDot_Product(	complex<double> alpha, 
				vector< complex< double > > &A,
		    		vector< complex< double > > &B)
{
	vector< complex< double > > zA(A.size());

	for(int64_t l = 0; l < A.size(); l++)
		zA[l] = conj(A[l]);

	complex< double > prod;
	if(zA.size() == B.size() )
	{
		int64_t m = zA.size();
		zgemm('N', 'N', 1, 1, m, alpha, &*zA.begin(), 1, &*B.begin(), m, 0.0, &prod, 1);
	}
	else
		cout << "\n\n Error in zInner_Product sizes \n\n";  

	return prod;
}

void zMatrix_Vector_Product(	zMatrix M,
			    	vector< complex< double > > v,
			    	vector< complex< double > > &w)
{
	if(M.NCols() == v.size())
	{
		int64_t m = M.NRows();
		int64_t k = M.NCols();

		zgemm('N', 'N', m, 1, k, 1.0, &*M.begin(), m, &*v.begin(), k, 0.0, &*w.begin(), m);
	}
	else
		cout << "\n\n Error in zMatrix_Vector_Product sizes \n\n";
} 

void zVector_Sum(	complex< double > a,
			vector< complex< double > > X,
			vector< complex< double > > &Y)
{
	int64_t n = X.size();
	zaxpy(n, a, &*X.begin(), 1, &*Y.begin(), 1);
}


