
#include "global.h"
#include "idsClass.h"
#include "dm_readfile.h"
#include "Momentum_Basis.h"
#include "Entropy.h"
#include "read_vector.h"

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

	vector<int> Cut, T1, T2;

	readpCuts(filein, Cut, T1, T2, filepEV1, filepEV2);

	int64_t n_p1, n_p2;
	double kx1, ky1, kx2, ky2;
	vector<int64_t> pBasis1, pBasis2;
	vector<double> norm1, norm2;
	vector< complex<double> > v_coeff1, v_coeff2;

	ofstream outfile(fileout.c_str());
	outfile << "\n pEV1_file path = " << filepEV1;
	outfile << "\n pEV2_file path = " << filepEV2;
	outfile << endl;

	outfile << "T1=";
	for(int i = 0; i < T1.size();i++)
		outfile << T1[i] << " ";
	outfile << endl;
	outfile << "T2=";
	for(int i = 0; i < T2.size();i++)
		outfile << T2[i] << " ";
	outfile << endl;

	outfile << "\nCUT\t";	
	for(int i = 0; i < Cut.size(); i++)
		outfile << Cut[i] << " ";
	outfile << endl;


	double t_read = wall_clock();
	read_pEV(n_p1, kx1, ky1, pBasis1, norm1, v_coeff1, filepEV1);
	read_pEV(n_p2, kx2, ky2, pBasis2, norm2, v_coeff2, filepEV2);
	t_read = wall_clock() - t_read; 	
	outfile << "\nTime to read momentum eigenvectors = " << (t_read)/60.0 << " minutes";
	outfile << endl;

	int64_t n_dm = int64_t_power(2, Cut.size());

	double ts_dm = wall_clock();
	zMatrix denm(n_dm, n_dm, 0.0);

	pdensity_matrixij(denm, v_coeff1, pBasis1, norm1, kx1, ky1, v_coeff2, pBasis2, norm2, kx2, ky2, T1, T2, Cut);

	ts_dm = wall_clock()-ts_dm;

	outfile << "\n\n Time to construct density matrices : " << (ts_dm)/60.0 << " minutes \n"; 

	double ts_diag = wall_clock();

	int count = 0;	
	for(int i = 0; i < n_dm; i++)
	{
		for(int j = 0; j < n_dm; j++)
			if ( abs(real(denm(i, j) ) ) > 0.0 || abs(imag(denm(i, j) ) > 0.0) )		
				count++;
	}		

	outfile << "\n Number of non-zero elements = " << count << endl;

//	complex<double> EEdm = Von_Neumann_Entropy(denm);
	complex<double> EEdm = Renyi_Entropy(denm);	

	outfile << "\nEntanglement entropy " << EEdm;
	ts_diag = wall_clock()-ts_diag;
	
	outfile << "\n\nTime of matrix diagonalization: " << (ts_diag)/60.0 << " minutes \n"; 
	
	outfile << "\nDENSITY_MATRIX\n\n";
	
	for(int i = 0; i < n_dm; i++)
	{
		for(int j = 0; j < n_dm; j++)
			outfile << fixed << setprecision(8) << denm(i,j) << "\t";
		outfile << endl;
	}		

	outfile.close();

	return 0;
}
