
#include "ARPACK_sector.h"

void ARPACK_Iteration_EV(	vector<int64_t> &pblock_states,
				vector<double> &norm,
				vector<double> &norm_inv,
				vector<int64_t> &rep_loc,
				vector< complex<double> > &state_ph,
				int &n_spins,
				int &n_it,
		      		double &px,
				double &py,
		      		vector < complex<double> > &Evals,
				zMatrix &Evecs,
		      		bool find_ev,
				double lambda,
				vector<coordinates> &adj_list,
				double &Ch,
				vector<triangles> &ch_list,
				double &ChCh,
				vector<bowties> &chch_list,
		      		double tol,
		      		int maxiter,
				vector<int> &T1,
				vector<int> &T2, 
				vector<int64_t> &Ia,
				vector<int64_t> &Ib,
				int &bits_right,
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

	int ncv = 4*nev; //4*nev;	// Parameter that sets the number of Arnoldi vectors stored in ARPACK
	if (ncv > n)
		ncv = n;

	int ldv = n;
	complex<double> *v;
	v = new complex< double >[ldv*ncv];

	outfile << " ARPACK PARAMETERS " << endl;
	outfile << " ncv = " << ncv << endl;
	outfile << " nev = " << nev << endl;
	outfile << " n = " << n << endl;
	outfile << " ldv = " << ldv << endl;
	outfile << " n_it = " << n_it << endl;
	outfile << "----------------------\n" << endl;

	int *iparam; iparam = new int[11];
	iparam[0] = 1;
	iparam[2] = maxiter;	// # of Max interations
	iparam[6] = 1;

	int *ipntr;
	ipntr = new int[14];

	complex<double> *workd;
	workd = new complex<double>[3*n];

	int lworkl = 3*ncv*ncv + 5*ncv;
	complex< double > *workl;
	workl = new complex< double >[lworkl];

	double *rwork;
	rwork = new double[ncv];

	int info = 0;
	int rvec = find_ev;  // Changed from above

	int *select; 
	select = new int[ncv];

	complex< double > *d; 
	d = new complex< double >[ncv];
	double sigma;

	complex< double > *workev;
	workev = new complex< double >[3*ncv];

	int ierr;

	double wc_Hv = 0.0; 
	double wc_zna = 0.0;
	double wc_zne = 0.0;
	  
	int itit = 0, itzit = 0;

	do {
		clock_t t_zna; t_zna = clock();
		
		double tmp_wc_zna = wall_clock();
		znaupd_(&ido, bmat, &n, which, &nev, &tol, resid, 
			    &ncv, v, &ldv,
			    iparam, ipntr, workd, 
			    workl,  &lworkl, rwork, &info);

		tmp_wc_zna = wall_clock() - tmp_wc_zna;
		wc_zna += tmp_wc_zna;
		outfile << " Wall clock for each iteration of znaupd " << tmp_wc_zna/60 << " minutes" << endl;

		itzit++;
		if ((ido==1)||(ido==-1))
		{
			itit++;
			double wc_start = wall_clock();
		
// H_v routine	
			sHv(	workd + ipntr[0]-1, workd + ipntr[1]-1,
				pblock_states, norm, norm_inv, rep_loc, state_ph,
				lambda, adj_list, Ch, ch_list, ChCh, chch_list, n_spins, Ia, Ib, bits_right, T1, T2, px, py, outfile);
	
			double wc_end = wall_clock();
			outfile << "Wall clock for Hv step " << (wc_end - wc_start)/60.0 << endl;
			wc_Hv += wc_end-wc_start;
		}
	} while ((ido==1)||(ido==-1));

	outfile << "\nARPACK iterations\n";
	outfile << "# iterations (# of H*v calls)\t" << itit << endl;

	outfile << "Wall Clock time H*v\t\t" << wc_Hv/60.0 << " minutes\n\n";
	outfile << "Max iter for Arnoldi\t\t" << iparam[2] << endl;
	outfile << "Wall Clock time znapud\t\t" << wc_zna/60.0 << " minutes\n\n";
	

	if (info<0) {
		outfile << "Error with znaupd, info = " << info << "\n";
		outfile << "Check documentation in dsaupd\n\n";
	}
	else
	{
		char howmany = 'A';

		zneupd_(&rvec, &howmany, select, d,
		    	v, &ldv, &sigma, workev,
		    	bmat, &n, which, &nev, &tol, resid, &ncv,
		   	v, &ldv, iparam, ipntr,
			workd, workl, &lworkl, rwork, &ierr);

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
//		delete v;
		delete d;
		delete iparam;
		delete ipntr;
		delete workd;
		delete workl;
		delete select;
		
		if(find_ev)
			for (int i=0; i<nev; i++)
				for (int64_t j=0; j<n; j++)
				{
					Evecs.resize(n_it, n);	
					Evecs(i,j) = v[i*n+j];
				}

		delete v;
	}
}

void ARPACK_sector(	vector<int64_t> &pblock_states,
			vector<double> &norm,
			vector<double> &norm_inv,
			vector<int64_t> &rep_loc, 
			vector< complex<double> > &state_ph,
			vector<complex<double> > &pblock_eigs,
			zMatrix &pblock_eigvecs,
			bool find_ev,
			double &lambda,
			int &n_spins,
			double &px, double &py,
			int &n_ones,
			int64_t &n_p,
			vector<coordinates> &adj_list,
			double &Ch,
			vector<triangles> &ch_list,
			double &ChCh,
			vector<bowties> &chch_list,
			vector<int> &T1,
			vector<int> &T2,
			int n_ev,
			double tol,
			int maxiter,
			vector<int64_t> &Ia,
			vector<int64_t> &Ib,
			int &bits_right,
			ofstream &outfile)
{
	int it = 0;
	if (n_ev > n_p - 1)
	  it = n_p;
	else
	  it = n_ev;
	 
	pblock_eigs.resize(it);

	if (n_p > n_ev + 1)
	{
		double	start = wall_clock();		
			  
		ARPACK_Iteration_EV(	pblock_states, norm, norm_inv, rep_loc, state_ph, n_spins, it, px, py,
			  	 	pblock_eigs, pblock_eigvecs, find_ev,
			   		lambda, adj_list, Ch, ch_list, ChCh, chch_list, tol, maxiter, T1, T2, Ia, Ib, bits_right, outfile);

		double	end = wall_clock();
		outfile << "\n Wall Clock for ARPACK\t" << (end - start)/60 << " minutes\n";

	}
	else
	{
		outfile << " Here \n\n";
		it = n_p;
		pblock_eigs.clear();
		pblock_eigs.resize(n_p);
		zMatrix pBij(n_p, n_p); pblock_eigvecs.resize(n_p, n_p);

		zXXZ_Heisenberg(pBij, n_spins, pblock_states, norm, px, py, lambda, adj_list, Ch, ch_list, T1, T2);
		zMatrix_Diagonalize(pBij, pblock_eigs, pblock_eigvecs);
		
		
		for(int64_t ii = 0; ii < pBij.NRows(); ii++)
		{
			for(int64_t jj = 0; jj < pBij.NCols(); jj++)
				outfile << fixed << setprecision(5) << pBij(ii,jj) << "\t\t";
			outfile << endl;
		}		
	}

	outfile << "\nEigenvalues\n";
	
	for(int kb = 0; kb < it; kb++)
	{outfile << setprecision(10) << pblock_eigs[kb] << "\t";}

}


void sector_eigs_eigvecs(	int &kx,
				int &ky,
				double &cSz,
				int &n_spins,
				double &lambda,
				double &J2,
				vector<coordinates> &adj_list,
				double &Ch,
				vector<triangles> &ch_list,
				double &ChCh,
				vector<bowties> &chch_list,
				vector<int> &T1,
				vector<int> &T2,
				int &n_sec_ev,
				bool &find_ev,
				bool &unwrap_peigvec,
				string &filepath_EV,
				double &tol,
				int &maxiter,
				ofstream &outfile)
{

	int n_ones = n_spins/2;
	
	for(double i = 0.0; i < cSz; i += 1)
	{	n_ones++;} 

	int L1 = 0, L2 = 0;
	LxLy(L1, L2, T1, T2);

	if(kx >= L1 || ky >= L2 || kx < 0 || ky < 0)
	{
		outfile << "\nInvalid kx or ky\n";
		exit(1);
	}

	double px = 2*kx*M_PI/L1, py = 2*ky*M_PI/L2;
	vector<int64_t> pblock_states;
	vector<double> norm, norm_inv;
	vector<int64_t> rep_loc;
	vector< complex<double> > state_ph;
	int64_t n_p = 0;
	vector<int64_t> Ia, Ib;
	int bits_right;

	double Sz = 0.5*n_ones - 0.5*(n_spins - n_ones);

	outfile << "\n--------------------------------------------------------------------\n"; 
	double t_states = wall_clock();

	pConstruct_States(	pblock_states, norm, norm_inv, T1, T2, n_ones, px, py, n_p, 
				rep_loc, state_ph, Ia, Ib, bits_right, outfile);

	t_states = wall_clock() - t_states;

	outfile << endl << "Total size of all states/rep loc vector " << rep_loc.size() << endl;
	outfile << " SECTOR : " <<  Sz << "\t(" << px << "," << py << ")\t" << n_p << endl;	
	outfile << "Time to construct pStates\t " << t_states/60.0 << " minutes" << endl;

	if(n_p == 0)
	{
		outfile << "\nNo states in psector\n";
		exit(1);
	}
	outfile << "-------------------------------------------------------\n" << endl;

	vector< complex<double> > pblock_eigs;
	zMatrix pblock_eigvecs;

	ARPACK_sector(	pblock_states, norm, norm_inv, rep_loc, state_ph,
			pblock_eigs, pblock_eigvecs, find_ev,
			lambda, n_spins, px, py, n_ones, n_p,
			adj_list, Ch, ch_list, ChCh, chch_list, T1, T2, n_sec_ev, tol, maxiter, Ia, Ib, bits_right, outfile);


	outfile  << "\n\nEnergies\t\tEnergies/site\t\tSz\tkx\tky\tJz\tJ2\n\n";
	
	for(int i = pblock_eigs.size() - 1; i >= 0; i--)	
		outfile << fixed << setprecision(15) << real(pblock_eigs[i]) << "\t" << real(pblock_eigs[i])/double(n_spins) <<"\t" << setprecision (3) << Sz << "\t" << px << "\t" << py << "\t" << lambda << "\t" << J2 << "\n";  	

	
	if(find_ev)
	{	
		// double ts_unwrap = wall_clock();
		
		for(int iib = 0; iib < n_sec_ev; iib++)
		{
			int t_tmp_int = iib + 1 + (int)'0';
			char t_tmp = (char)t_tmp_int;

			string filepath_tmp;

			if(unwrap_peigvec)
				filepath_tmp = filepath_EV + "_EV" + t_tmp + ".txt";		
			else
				filepath_tmp = filepath_EV + "_pEV" + t_tmp + ".txt";
	
			ofstream outfile_vecs(filepath_tmp.c_str());
			outfile_vecs << "\n Eigenvector for state \n";
			outfile_vecs << "Sz\t(kx,ky)\tJz/lambda\n";
			outfile_vecs << Sz <<"\t(" << kx << "," << ky << ")\t" << lambda << endl;
		
			outfile_vecs << "\n\n-------------------------------------------";
			outfile_vecs << "\nState Energy " << pblock_eigs[iib] << endl;

			if(unwrap_peigvec)
			{
				vector< complex <double> > v_coeff(n_p);
				
				for(int64_t jb = 0; jb < n_p; jb++)
					v_coeff[jb] = pblock_eigvecs(iib, jb);

				vector<ids> psi = Unwrap_pState(v_coeff, pblock_states, norm, px, py, T1, T2);	
			
				outfile_vecs << "\nEIGENVECTOR\n";
				for(int64_t jb = 0; jb < psi.size(); jb++)
				{
					outfile_vecs << psi[jb].stateID << "\t" << fixed << setprecision(15) << psi[jb].coeff << endl;
				}
			
			}
			else
			{
				outfile_vecs << "\nMomentum EIGENVECTOR\n";
			
				outfile_vecs << "\nnp \t\t px \t\tpy\n";
				outfile_vecs << n_p << "\t\t" << px << "\t\t" << py << endl;

				outfile_vecs << "\n\npblock_states \t\t Norm \t\t\t v_coeff\n\n";
				
				for(int64_t jb = 0; jb < n_p; jb++)
					outfile_vecs << pblock_states[jb] << "\t\t\t" << norm[jb] << "\t\t\t" << fixed << setprecision(15) << pblock_eigvecs(iib, jb) << endl;

			}
			outfile_vecs.close();
		}
		// double te_unwrap = wall_clock();
		// outfile << "\nTime to unwrap state = " << (te_unwrap - ts_unwrap)/60 << " minutes\n";			
	}
}

