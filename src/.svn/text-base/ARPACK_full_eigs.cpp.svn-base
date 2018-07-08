
#include "ARPACK_full_eigs.h"

void ARPACK_Iteration(	vector<int64_t> &pblock_states,
			vector<double> &norm,
			int &n_spins,
			int &n_it,
	//		int64_t &n,
			double &px,
			double &py,
			vector < complex<double> > &Evals,
		      	double lambda,
			vector<coordinates> &adj_list,
			double &Ch,
			vector<triangles> &ch_list,
			double tol,
		      	int maxiter,
			vector<int> &T1,
			vector<int> &T2, 
			ofstream &outfile)
{
 	int64_t n = pblock_states.size();
	int nev = n_it;

	if(nev >= n - 1 && n > 0)
	nev = n - 2;

	int ido = 0;
	char bmat[2] = "I";
	char which[3] = "SR";

	complex< double > *resid;
	resid = new complex< double >[n];

	int ncv = 4*nev; //4*nev;
	if (ncv > n)
	ncv = n;

	int ldv = n;
	//  zMatrix v(ldv, ncv);
	complex<double> *v;
	v = new complex< double >[ldv*ncv];

	int *iparam; iparam = new int[11];
	iparam[0] = 1;
	iparam[2] = maxiter;	// # of Max interations
	iparam[6] = 1;

	int *ipntr;
	ipntr = new int[14];
	//  vector<int> ipntr(14);

	complex<double> *workd;
	workd = new complex<double>[3*n];

	int lworkl = 3*ncv*ncv + 5*ncv;
	complex< double > *workl;
	workl = new complex< double >[lworkl];

	double *rwork;
	rwork = new double[ncv];

	int info = 0;
	int rvec = 0;  // Changed from above

	int *select; 
	select = new int[ncv];

	complex< double > *d; 
	d = new complex< double >[ncv];
	double sigma;

	complex< double > *workev;
	workev = new complex< double >[3*ncv];

	int ierr;


	clock_t t_Ham;
	clock_t t_Hv;
	t_Hv = clock();
	double time_Hv_max = 0.0;
	double time_arpacki = 0.0;
	double wc_Hv = 0.0; 
	///// Stores the Full Sparse Hamiltonian each time
	sMatrix pBij(n, n);
	t_Ham = clock();
	sXXZ_Heisenberg(pBij, n_spins, pblock_states, norm, px, py,
			   lambda, adj_list, Ch, ch_list, T1, T2, outfile);
	t_Ham = clock() - t_Ham;
	
	pblock_states.clear();
	norm.clear();	
  
 	outfile << "Time to construct H\t\t" << double(t_Ham)/double(CLOCKS_PER_SEC)/60 << " minutes" << endl;
	int itit = 0, itzit = 0;
	double total_zna = 0.0;
	double total_Hv = 0.0;
	do {
		clock_t t_zna; t_zna = clock();
		znaupd_(&ido, bmat, &n, which, &nev, &tol, resid, 
			    &ncv, v, &ldv,
			    iparam, ipntr, workd, 
			    workl,  &lworkl, rwork, &info);
		t_zna = clock() - t_zna;
		total_zna = total_zna + t_zna;
		itzit++;
		if ((ido==1)||(ido==-1))
		{
			itit++;
			t_Hv = clock();	
			double wc_start = wall_clock();
			sHv(pBij, workd + ipntr[0]-1, workd + ipntr[1]-1);
			t_Hv = clock() - t_Hv;
			double wc_end = wall_clock();
			total_Hv = total_Hv + t_Hv;
			wc_Hv = wc_Hv + (wc_end-wc_start);
		}
	} while ((ido==1)||(ido==-1));

	outfile << "\nARPACK iterations\n";
	outfile << "# iterations (# of H*v calls)\t" << itit << endl;

	outfile << "Total time for H*v\t\t" << double(total_Hv)/double(CLOCKS_PER_SEC)/60 << " minutes" << endl;
	outfile << "Wall Clock time H*v\t\t" << wc_Hv/60.0 << " minutes\n\n";
	outfile << "Max iter for Arnoldi\t\t" << iparam[2] << endl;
	outfile << "Total time for znaupd \t\t" << total_zna/double(CLOCKS_PER_SEC)/60 << " minutes" << endl;
	if (info<0) {
		outfile << "Error with znaupd, info = " << info << "\n";
		outfile << "Check documentation in dsaupd\n\n";
	}
	else
	{
		char howmany = 'A';
		clock_t t_z;
		t_z = clock();
		zneupd_(&rvec, &howmany, select, d,
		    v, &ldv, &sigma, workev,
		    bmat, &n, which, &nev, &tol, resid, &ncv,
		    v, &ldv, iparam, ipntr,
		    workd, workl, &lworkl, rwork, &ierr);
		t_z = clock() - t_z;
		outfile << "Time for zneupd \t\t" << double(t_z)/double(CLOCKS_PER_SEC)/60 << " mintues " << endl;
		if (ierr!=0) {
			outfile << "Error with zneupd, info = " << ierr << "\n";
			outfile << "Check the documentation of zneupd.\n\n";
		} else if (info==1) {
			outfile << "Maximum number of iterations reached.\n\n";
		} else if (info==3) {
			outfile << "No shifts could be applied during implicit\n";
			outfile << "Arnoldi update, try increasing NCV.\n\n";
		}

		for (int i=0; i<nev; i++){
			Evals[i] = d[i];
		}

		delete resid;
		delete v;
		delete d;
		delete iparam;
		delete ipntr;
		delete workd;
		delete workl;
		delete select;
	}
}

void ARPACK(	pVector &LanczosEij,
		double &lambda,
		int &n_spins,
		double &sector_Sz,
		vector<coordinates> &adj_list,
		double &Ch,
		vector<triangles> &ch_list,
		vector<int> &T1,
		vector<int> &T2,
		int n_ev,
		double tol,
		int maxiter,
		ofstream &outfile)
{
	int L1 = 0, L2 = 0;
	LxLy(L1, L2, T1, T2);
	int global_n_it = 0;
	int Sz_min;
	int Sz_max = n_spins;

	if(n_spins % 2 == 0)
		Sz_min = n_spins/2;
	else
		Sz_min = (n_spins+1)/2;

	if(sector_Sz > -1.0)
	{
		Sz_min += int(sector_Sz);
		Sz_max = Sz_min;
	}	

  	for(int n_ones = Sz_min; n_ones <= Sz_max; n_ones++) 
	{
		for(int kx = 0; kx < L1; kx++)
		{
			for(int ky = 0; ky < L2; ky++)
			{
				double px = 2*kx*M_PI/L1, py = 2*ky*M_PI/L2;
				vector<int64_t> pblock_states;
				vector<double> Norm;
				int64_t n_p = 0;

				double Sz = 0.5*n_ones - 0.5*(n_spins - n_ones);

				outfile << "\n--------------------------------------------------------------------\n"; 
				clock_t t_states;
				t_states = clock();
				Construct_States(pblock_states, Norm, T1, T2, n_ones, px, py, n_p, outfile);
				t_states = clock() - t_states;

				outfile << "\n SECTOR : " <<  Sz << "\t(" << px << "," << py << ")\t" << n_p << endl;	
				outfile << "Time to construct pStates\t " << double(t_states)/double(CLOCKS_PER_SEC)/60 << " minutes\n";

				if(n_p == 0)
				  continue;

				int it = 0;
				if (n_ev > n_p - 1)
				  it = n_p;
				else
				  it = n_ev;
				 
				vector< complex<double> > pblock_eigs(it);
				 
				if (n_p > n_ev + 1)
				{
					clock_t t_arpack;
					t_arpack = clock();
					double	start = wall_clock();		
						  
					ARPACK_Iteration(pblock_states, Norm, n_spins, it, px, py,
						   pblock_eigs, // pblock_eigvecs,
						   lambda, adj_list, Ch, ch_list, tol, maxiter, T1, T2, outfile);
					t_arpack = clock() - t_arpack;

					double	end = wall_clock();

					outfile << "\nTime required by ARPACK\t" << double(t_arpack)/double(CLOCKS_PER_SEC)/60 << " minutes " << endl;
					outfile << "\n Wall Clock for ARPACK\t" << (end - start)/60 << " minutes\n";

				}
				else
				{
					it = n_p;
					pblock_eigs.clear();
					pblock_eigs.resize(n_p);
					zMatrix pBij(n_p, n_p), peigvecs(n_p, n_p);

					zXXZ_Heisenberg(pBij, n_spins, pblock_states, Norm, px, py, lambda, adj_list, Ch, ch_list, T1, T2);
					zMatrix_Diagonalize(pBij, pblock_eigs, peigvecs);
					
				}
				
				pblock_states.clear();
				Norm.clear();

				outfile << "\nEigenvalues\n";

				for(int k = global_n_it, kb = 0; kb < it; k++, kb++)
				{	
					outfile << setprecision(10) << pblock_eigs[kb] << "\t";
					LanczosEij.push_back(k, real(pblock_eigs[kb]), px, py, Sz);
				}

				global_n_it = global_n_it + it;
			}
		}
	}
}
    

