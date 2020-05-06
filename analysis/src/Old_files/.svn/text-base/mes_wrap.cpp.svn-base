
#include "global.h"
#include "idsClass.h"
#include "mes_readfile.h"
#include "Entropy.h"
#include "read_vector.h"

void add_density_matrices(	complex<double> &c1,
				complex<double> &c2,
				zMatrix &dm11,
				zMatrix &dm12,
				zMatrix &dm21,
				zMatrix &dm22,
				zMatrix &dmsum)
{
	#pragma omp parallel for
	for(int64_t i = 0; i < dm11.NRows(); i++)
		for(int64_t j = 0; j < dm11.NRows(); j++)
			dmsum(i,j) = conj(c1)*c1*dm11(i,j) + conj(c1)*c2*dm12(i,j) + conj(c2)*c1*dm21(i,j) + conj(c2)*c2*dm22(i,j);	
}


int main(int argc, char *argv[])
{
	string filein, filepEV1, filepEV2, fileout;

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
	read_mes(filein, Cut, filepEV1, filepEV2, T1, T2, alpha, c_min, interval, delta);
		
	int ntot = 1/delta;
	double c_max = c_min + interval;

	ofstream outfile(fileout.c_str());
	
	outfile << " Vector file path 1: " << filepEV1 << endl;
	outfile << " Vector file path 2: " << filepEV2 << endl << endl;

	int64_t n_p1, n_p2;
	double px1, px2, py1, py2;
	vector<int64_t> pBasis1, pBasis2;
	vector<double> norm1, norm2;
	vector< complex <double> > v_coeff1, v_coeff2;

	outfile << "\nT1=";
	for(int i = 0; i < T1.size(); i++)
		outfile << T1[i] << " " ;
	outfile << "\nT2=";
	for(int i = 0; i < T2.size(); i++)
		outfile << T2[i] << " " ;
	outfile << endl;

	double t_read = wall_clock();
	read_pEV(n_p1, px1, py1, pBasis1, norm1, v_coeff1, filepEV1);
	read_pEV(n_p2, px2, py2, pBasis2, norm2, v_coeff2, filepEV2);
	t_read = wall_clock() - t_read;
	outfile << "\nTime to read p-states = " << t_read/(60.0) << " minutes" << endl; 
/*
	vector<ids> psi1, psi2;
	
	double t_unwrap = wall_clock();
//	Unwrap_pState(v_coeff1, pBasis1, norm1, px1, py1, T1, T2, psi1);
//	Unwrap_pState(v_coeff2, pBasis2, norm2, px2, py2, T1, T2, psi2);
	t_unwrap = wall_clock() - t_unwrap;
	outfile << "Time to unwrap p-states = " << t_unwrap/(60.0) << " minutes" << endl; 
	outfile << " \nSize of unwrapped EV1: " << psi1.size() << endl;
	outfile << " Size of unwrappedn EV2: " << psi2.size() << endl << endl;
*/

	// Sizes of reduced density matrices along two different cuts
	int64_t n_dm = int64_t_power(2, Cut.size());
	
	outfile << "\n Size of density matrix = " << n_dm << " x " << n_dm << endl;

	// Initialize and compute the density matrices for the two vectors
	double t_dm = wall_clock();
	zMatrix	dm11(n_dm, n_dm, 0.0);
	zMatrix dm12(n_dm, n_dm, 0.0);
	zMatrix	dm21(n_dm, n_dm, 0.0);
	zMatrix dm22(n_dm, n_dm, 0.0);

	pdensity_matrix(dm11, v_coeff1, pBasis1, norm1, px1, py1, v_coeff1, pBasis1, norm1, px1, py1, T1, T2, Cut);
	pdensity_matrix(dm12, v_coeff1, pBasis1, norm1, px1, py1, v_coeff2, pBasis2, norm2, px2, py2, T1, T2, Cut);
	pdensity_matrix(dm21, v_coeff2, pBasis2, norm2, px2, py2, v_coeff1, pBasis1, norm1, px1, py1, T1, T2, Cut);
	pdensity_matrix(dm22, v_coeff2, pBasis2, norm2, px2, py2, v_coeff2, pBasis2, norm2, px2, py2, T1, T2, Cut);
		
	outfile << "\n EE dm11 = " << Renyi_Entropy(dm11);
	outfile << "\n EE dm22 = " << Renyi_Entropy(dm22);

	t_dm = wall_clock() - t_dm;
	outfile << "\n Cut \t";
	for(int i = 0; i < Cut.size(); i++)
		outfile << Cut[i] << " ";
	outfile << "\n Total time to construct density matrices  " << (t_dm)/60.0 << " minutes\n";
	/////////////////////////////////////////////////////////////////////

	// Look for the minimally entangled linear combination of eigenstates
//	int n_tot = 100;	


	double ts_EE = wall_clock();
//	outfile << fixed << setprecision(6) << "c1|psi1>" << "\t\t" << "c2|psi2>" << "\t\t" << "E.Entropy" << endl;
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

	outfile << "\n Total time to add density matrices " << (t_tot_add)/60.0 << " minutes \n"; 	outfile << "\n Total time to diagonalize density matrix " << (t_tot_diag)/60.0 << " minutes \n"; 		
	outfile << "\n Total time for EE calc  " << (ts_EE)/60.0 << " minutes\n";
	outfile.close();

	return 0;
}
