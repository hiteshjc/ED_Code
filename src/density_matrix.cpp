
#include "global.h"
#include "idsClass.h"
#include "density_matrix_readfile.h"
#include "Entropy.h"

void add_density_matrices(	complex<double> &c1,
				complex<double> &c2,
				zMatrix &dm11,
				zMatrix &dm12,
				zMatrix &dm21,
				zMatrix &dm22,
				zMatrix &dmsum)
{
	for(int64_t i = 0; i < dm11.NRows(); i++)
		for(int64_t j = 0; j < dm11.NRows(); j++)
			dmsum(i,j) = conj(c1)*c1*dm11(i,j) + conj(c1)*c2*dm12(i,j) + conj(c2)*c1*dm21(i,j) + conj(c2)*c2*dm22(i,j);	
}


int main(int argc, char *argv[])
{
	string filein, fileEV1, fileEV2, fileout;

	if (argc <= 1)
	{
		cout << "Usage: " << argv[0] << " <Filename>" << endl;
		exit(1);
	}

	filein=argv[1];
	fileout=argv[2];

	vector<int> Cut;
	readCuts(filein, Cut, fileEV1, fileEV2);

	vector<ids> psi[2]; // = add_ids(psi1, psi2, c, ph);

	readEV(psi[0], fileEV1);
	readEV(psi[1], fileEV2);

	ofstream outfile(fileout.c_str());

	outfile << "\nCUT\t";	
	for(int i = 0; i < Cut.size(); i++)
		outfile << Cut[i] << " ";
	outfile << endl;
	
	outfile << fileEV1 << endl;
	outfile << fileEV2 << endl << endl;

	for(int i = 0; i < 2; i++)
		outfile << " Size Vector psi = " << psi[i].size() << endl;
		
	outfile << "\n\n Dot Product for vectors " << dot_product_ids(psi[0],psi[1]);
	
	int64_t n_dm = int64_t_power(2, Cut.size());

	double ts_dm = wall_clock();

	zMatrix denm(n_dm, n_dm, 0.0);

	density_matrixij(denm, psi[1], psi[0], Cut);

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
