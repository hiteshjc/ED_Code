#include "global.h"
#include "idsClass.h"
#include "pmes_readfile.h"
#include "Entropy.h"
#include "read_inputfile.h"
#include "unwrap_pstate.h"

int main(int argc, char *argv[])
{
	string filein, filepEV1, filepEV2, filepEV3, filepEV4, fileout;

	if (argc <= 1)
	{
		cout << "Usage: " << argv[0] << " <Filename>" << endl;
		exit(1);
	}

	filein=argv[1];
	fileout=argv[2];

	double c_min = 0.0, interval = 1.0;	// default values
	double delta = 0.01;

	// Read two Cuts along non-trivial directions from file
	vector<int> Cut, T1, T2;
	int alpha=2;
	read_mes_chiral(filein, Cut, filepEV1, filepEV2, filepEV3, filepEV4, T1, T2, alpha, c_min, interval, delta);

	int N; double Sz;
	read_N_Sz(filein, N, Sz);
	int nones = int(2*Sz);
	nones = (nones+N)/2;

	int ntot = 1/delta;
	double c_max = c_min + interval;

	ofstream outfile(fileout.c_str());
	
	outfile << " Vector file path 1: " << filepEV1 << endl;
	outfile << " Vector file path 2: " << filepEV2 << endl << endl;
	outfile << " Vector file path 3: " << filepEV3 << endl;
	outfile << " Vector file path 4: " << filepEV4 << endl << endl;

	outfile << "\nT1=";
	for(int i = 0; i < T1.size(); i++)
		outfile << T1[i] << " " ;
	outfile << "\nT2=";
	for(int i = 0; i < T2.size(); i++)
		outfile << T2[i] << " " ;
	outfile << endl;

	double t_unwrap = wall_clock();
	
	vector<ids> psi1, psi12, psi2, psi22;
	read_n_unwrap(filein, filepEV1, T1, T2, psi1, outfile);
	read_n_unwrap(filein, filepEV2, T1, T2, psi12, outfile);
	read_n_unwrap(filein, filepEV3, T1, T2, psi2, outfile);
	read_n_unwrap(filein, filepEV4, T1, T2, psi22, outfile);
	
	t_unwrap = wall_clock() - t_unwrap;
	outfile << "Time to unwrap p-states = " << t_unwrap/(60.0) << " minutes" << endl; 
	outfile << "\n Size of unwrapped EV1: " << psi1.size() << endl;
	outfile << " Size of unwrapped EV2: " << psi12.size() << endl << endl;
	outfile << " Size of unwrapped EV3: " << psi2.size() << endl << endl;
	outfile << " Size of unwrapped EV4: " << psi22.size() << endl << endl;

	complex<double> c1=complex<double>(0.70710678118,0.0); 
	complex<double> c2=complex<double>(0.0,0.70710678118); 

	for(int i = 0; i< psi1.size(); i++) psi1[i].coeff = c1*psi1[i].coeff + c2*psi12[i].coeff; 	
	for(int i = 0; i< psi2.size(); i++) psi2[i].coeff = c1*psi2[i].coeff + c2*psi22[i].coeff; 	
	
	outfile << "\n Got both psi's of the same (positive) chirality " << endl;
	// Sizes of reduced density matrices along two different cuts
	int64_t n_dm = int64_t_power(2, Cut.size());
	
	outfile << "\n Size of density matrix = " << n_dm << " x " << n_dm << endl;

	// Initialize and compute the density matrices for the two vectors
	double t_dm = wall_clock();
	zMatrix	dm11(n_dm, n_dm, 0.0);
	zMatrix dm12(n_dm, n_dm, 0.0);
	zMatrix	dm21(n_dm, n_dm, 0.0);
	zMatrix dm22(n_dm, n_dm, 0.0);
	
	density_matrix(dm11, psi1, Cut);
	density_matrix(dm22, psi2, Cut);
	density_matrixij(dm12, psi1, psi2, Cut);
	density_matrixij(dm21, psi2, psi1, Cut);

	t_dm = wall_clock() - t_dm;
	outfile << "\n Cut \t";
	for(int i = 0; i < Cut.size(); i++)
		outfile << Cut[i] << " ";
	outfile << "\n Total time to construct density matrices  " << (t_dm)/60.0 << " minutes\n\n\n";
	/////////////////////////////////////////////////////////////////////

	double ts_EE = wall_clock();
	outfile << fixed << setprecision(6) << "c1|psi1>" << "\t\t" << "c2|psi2>" << "\t\t" << "E.Entropy \t\t c \t\t ph \t\t E.Entropy" << endl;

	int n_iter = interval/delta;
	double t_tot_add = 0.0;
	double t_tot_diag = 0.0;	

	for(int i = 0; i < n_iter; i++)
	{
		for(int j = 0; j <= ntot; j++)
		{
			double c = double(i)*delta+c_min;
			complex<double> ph = complex<double>(0.0, 2*M_PI*j/double(ntot));
				
			complex<double> c1 = c, c2 = sqrt(1-c*c)*exp(ph);
			zMatrix dmsum(n_dm, n_dm, 0.0);
			// Adding all the four density matrices together with the co-efficients

			double t_add = wall_clock();
			add_density_matrices(c1, c2, dm11, dm12, dm21, dm22, dmsum);
			t_add = wall_clock() - t_add;
			t_tot_add += t_add;

			double t_diag = wall_clock();
			complex<double> S_tmp;
			if(alpha == 1)
				S_tmp = Von_Neumann_Entropy(dmsum);
			else if (alpha == 2)
				S_tmp = Renyi_Entropy(dmsum);
			t_diag = wall_clock() - t_diag;
			t_tot_diag += t_diag;

			outfile << fixed << setprecision(6) << c1 << "\t" << c2 << "\t" << S_tmp << "\t" << c << "\t" << imag(ph) << "\t" << real(S_tmp) << endl;
		}
	}
	ts_EE = wall_clock() - ts_EE;

	outfile << "\n Total time to add density matrices " << (t_tot_add)/60.0 << " minutes \n"; 		
	outfile << "\n Total time to diagonalize density matrix " << (t_tot_diag)/60.0 << " minutes \n"; 		
	outfile << "\n Total time for EE calc  " << (ts_EE)/60.0 << " minutes\n";
	outfile.close();

	return 0;
}
