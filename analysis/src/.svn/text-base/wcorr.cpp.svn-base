//#include "corr_readfile.h"
#include "read_inputfile.h"
#include "spin_exp.h"
#include "unwrap_pstate.h"
#include "read_vector.h"


int main(int argc, char *argv[])
{
	string filein, filepEV, fileout;
 
	double t_sp = wall_clock();

	if (argc <= 1)
	{
		cout << "Usage: " << argv[0] << " <Filename>" << endl;
		exit(1);
	}

	filein=argv[1];
	fileout=argv[2];

	// Read triangles and eigenvector file
	vector<int> list;
	vector<links> bonds;
	bool BB = 0, SS = 0;
	bool Szcorr = 1;
	vector<int> T1, T2;

	read_SSBB(filein, SS, BB);
	read_list(filein, list);
	read_bonds(filein, bonds);
	read_T1T2(filein, T1, T2);
	read_pEV(filein, filepEV);

//	readCorr_List(filein, SS, list, BB, bonds, T1, T2, filepEV1);

	vector<ids> psi;
	int64_t n_p;
	double kx, ky;	
 	vector<int64_t> pBasis;
	vector<double> norm;
	vector< complex <double> > v_coeff;

	ofstream outfile(fileout.c_str());
	outfile << "\n pEV_file_path = " << filepEV;
	outfile << "\n Compute <SS> = " << SS;	
	outfile << "\n Compute <BB> = " << BB;

	outfile << "\nT1 = ";
	for(int i = 0; i < T1.size(); i++)
		outfile << T1[i] << " ";

	outfile << "\nT2 = ";	
	for(int i = 0; i < T2.size(); i++)
		outfile << T2[i] << " ";
	outfile << endl << endl;	


	double t_read = wall_clock();
	read_pEV(n_p, kx, ky, pBasis, norm, v_coeff, filepEV);
	t_read = wall_clock() - t_read; 	
	outfile << "\nTime to read momentum eigenvector = " << (t_read)/60.0 << " minutes";

//	double t_unwrap = wall_clock();
//	Unwrap_pState(v_coeff, pblock_states, norm, kx, ky, T1, T2, psi);	
//	t_unwrap = wall_clock() - t_unwrap;
//	outfile << "\nTime to unwrap eigenvector = " << (t_unwrap)/60.0 << " minutes";

	outfile << "\nSize of pket: " << pBasis.size() << endl;

	outfile << "\nList of site. Correlations are measured w.r.t first site\n";
	for(int m=0; m < list.size(); m++)
		outfile << list[m] << " "; 
	outfile << endl;

	t_sp = wall_clock();	

	bool com_ch = 0;
	int nlist = list.size();
	int L1, L2;
	LxLy(L1, L2, T1, T2);

	complex<double> cz = complex<double> (0.0, 0.0);
	vector< complex<double> > Szexp(nlist, cz);
	complex< double > argx = complex< double >(0.0, -kx);
	complex< double > argy = complex< double >(0.0, -ky);
			
	for(int m = 0; m < list.size(); m++)
	{
		int site_i = list[m];
	

		double Szreal = 0.0;
		double Szimag = 0.0;				

		#pragma omp parallel for default(shared) reduction(+:Szreal, Szimag)
		for(int64_t n1 = 0; n1 < v_coeff.size(); n1++)
		{	

			int64_t p_id = pBasis[n1];
			complex<double> tmp_coeff = v_coeff[n1]/norm[n1];
			vector<ids> tmp_ids;
			complex<double> ph = complex<double> (1.0, 0.0);

			for(int x = 0; x < L1; x++)
			{
				if(x > 0)
				{
					translateT(p_id, T1);	
					ph = ph*exp(argx);
					complex<double> tmpx_coeff = ph*tmp_coeff;
					add2ids(tmp_ids, p_id, tmpx_coeff); 

				}

				for(int y = 0; y < L2; y++)
				{
					if(y > 0)
					{
						translateT(p_id, T2);
						ph = ph*exp(argy);
						complex<double> tmpy_coeff = ph*tmp_coeff;	
						add2ids(tmp_ids, p_id, tmpy_coeff); 
					}
				} 
			}

			complex<double> tmpexp = complex<double> (0.0, 0.0);

			for(int l = 0; l < tmp_ids.size(); l++)
			{
		 		bool loci = btest64(tmp_ids[l].stateID, site_i);

				complex<double> st_coeff = tmp_ids[l].coeff;

				if (loci == 0)
				{	tmpexp -= 0.5*conj(st_coeff)*st_coeff;}
				else
				{	tmpexp += 0.5*conj(st_coeff)*st_coeff;}		
			}

			Szreal += real(tmpexp);
			Szimag += imag(tmpexp);
		}

		Szexp[m] = complex<double> (Szreal, Szimag);
		if(abs(Szimag) > 1e-6)
			com_ch = 1;
	}

	t_sp = wall_clock() - t_sp;

	outfile << "\n Total time for Sz expectations = " << (t_sp)/60.0 << " minutes\n\n";
	if(com_ch)
		outfile << "\n\nComplex expectations-Change OUTPUT TO COMPLEX no.\n\n";
	outfile << "\nSz EXPECTATION VALUES\n";

	for(int m=0; m < nlist; m++)
	{
		outfile << fixed << setprecision(8) << real(Szexp[m]) << "\t";
	}		

/*
	if(SS)
	{
		t_sp = wall_clock();
	
		bool com_ch = 0;
		int nlist = list.size();
		zMatrix SSmat(nlist, nlist, 0.0);
		
		zMatrix SSzzmat(nlist, nlist, 0.0);

		for(int m = 0; m < list.size(); m++)
		{
			int site_i = list[m];
		
			for(int n = m; n < list.size(); n++)
			{
				//outfile << " n = " << n << endl;
				int site_j = list[n];

				double sexpr = 0.0, sexpi = 0.0, szzr = 0.0, szzi = 0.0;	
//				complex<double> sexp = complex<double> (0.0, 0.0);				
//				complex<double> szzexp = complex<double> (0.0, 0.0);

				#pragma omp parallel for default(shared) reduction(+:sexpr, sexpi, szzr, szzi)
				for(int64_t n1 = 0; n1 < psi.size(); n1++)
				{	
					vector< complex<double> > SSij;
	
					SSij = spin_spin_exp(site_i, site_j, n1, psi, outfile);
					sexpr += real(SSij[0]);
					sexpi += imag(SSij[0]);
					szzr += real(SSij[1]);
					szzi += imag(SSij[1]);
				}
	
				SSmat(m, n) = complex<double> (sexpr, sexpi);
				SSzzmat(m,n) = complex<double> (szzr, szzi);

				if(abs(sexpi) > 1e-6)
					com_ch = 1;
				if (m != n)
				{					
					SSmat(n, m) = complex<double> (sexpr, -sexpi);
					SSzzmat(n, m) = complex<double> (szzr, -szzi);
				}
			}
		}

		t_sp = wall_clock() - t_sp;

		outfile << "\n\n\nTotal time for spin-spin correlations = " << (t_sp)/60.0 << " minutes";
		if(com_ch)
			outfile << "\n\nComplex expectations-Change OUTPUT TO COMPLEX no.\n\n";
		outfile << "\n\nSPIN SPIN EXPECTATION VALUES\n";

		for(int m=0; m < list.size(); m++)
		{
			for(int n=0; n < list.size(); n++)
				outfile << fixed << setprecision(8) << real(SSmat(m,n)) << "\t";
			outfile << endl;
		}		

		outfile << "\n\nSz Sz EXPECTATION VALUES\n";

		for(int m=0; m < list.size(); m++)
		{
			for(int n=0; n < list.size(); n++)
				outfile << fixed << setprecision(8) << real(SSzzmat(m,n)) << "\t";
			outfile << endl;
		}		

	}

	if(BB)
	{
		t_sp = wall_clock();
	
		bool com_ch = 0;
		int nbonds = bonds.size();
		zMatrix BBmat(nbonds, nbonds, 0.0);		
		zMatrix BBzzmat(nbonds, nbonds, 0.0);

		for(int m = 0; m < nbonds; m++)
		{
			links bond_i = bonds[m];

			for(int n=m; n < nbonds; n++)
			{
				links bond_j = bonds[n];

				double bexpr = 0.0, bexpi = 0.0, bzzr = 0.0, bzzi = 0.0;

				#pragma omp parallel for default(shared) reduction(+:bexpr, bexpi, bzzr, bzzi)
				for(int64_t n1 = 0; n1 < psi.size(); n1++)
				{
					vector< complex<double> > BBij;
					
					BBij = bond_bond_exp(bond_i, bond_j, n1, psi, outfile);

					bexpr += real(BBij[0]);
					bexpi += imag(BBij[0]);
					bzzr += real(BBij[1]);
					bzzi += imag(BBij[1]);
				}
		
			
				BBmat(m, n) = complex<double> (bexpr, bexpi);				
				BBzzmat(m, n) = complex<double> (bzzr, bzzi);
			
				if(abs(bexpi) > 1e-6)
					com_ch = 1;
				if (m != n)		
				{			
					BBmat(n, m) = complex<double> (bexpr, -bexpi);
					BBzzmat(n, m) = complex<double> (bzzr, bzzi);
				}
			}		
		}

		
		t_sp = wall_clock() - t_sp;
		outfile << "\n\n\nTotal time for bond-bond correlations = " << (t_sp)/60.0 << " minutes";
		if(com_ch)
			outfile << "\n\nComplex expectations-Change OUTPUT TO COMPLEX no.\n\n";
		
		outfile << "\n\nBOND BOND EXPECTATION VALUES\n";
		for(int m=0; m < nbonds; m++)
		{
			for(int n=0; n < nbonds; n++)
				outfile << fixed << setprecision(8) << real(BBmat(m,n)) << "\t";
			outfile << endl;
		}

		outfile << "\n\nBOND BOND ZZ EXPECTATION VALUES\n";
		for(int m=0; m < nbonds; m++)
		{
			for(int n=0; n < nbonds; n++)
				outfile << fixed << setprecision(8) << real(BBzzmat(m,n)) << "\t";
			outfile << endl;
		}
	}
*/
	outfile.close();

	return 0;
}
