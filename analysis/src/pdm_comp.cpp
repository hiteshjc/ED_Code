
#include "global.h"
#include "idsClass.h"
#include "Entropy.h"
#include "read_inputfile.h"
#include "unwrap_pstate.h"


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
	read_Cut(filein, Cut);
	read_T1T2(filein, T1, T2);
	read_pev1ev2(filein, filepEV1, filepEV2);

	ofstream outfile(fileout.c_str());

	outfile << filepEV1 << endl;
	outfile << filepEV2 << endl << endl;

	outfile << "\nCUT\t";	
	for(int i = 0; i < Cut.size(); i++)
		outfile << Cut[i] << " ";
	outfile << endl;

	vector <ids> psi1, psi2;

	double t_unwrap = wall_clock();
	read_n_unwrap( filein, filepEV1, T1, T2, psi1, outfile);
	read_n_unwrap( filein, filepEV2, T1, T2, psi2, outfile);
	t_unwrap = wall_clock() - t_unwrap;
	outfile << "\nTime to unwrap pEV's = " << t_unwrap/(60.0) << " minutes" << endl;

	outfile << "\n\n Dot Product for vectors " << dot_product_ids(psi1,psi2);
	
	int64_t n_dm = int64_t_power(2, Cut.size());

	double t_dm = wall_clock();

	zMatrix denm(n_dm, n_dm, 0.0);
	density_matrixij(denm, psi1, psi2, Cut);

	t_dm = wall_clock()-t_dm;
	outfile << "\n\n Time to construct density matrices : " << (t_dm)/60.0 << " minutes \n"; 

	double t_diag = wall_clock();

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
	t_diag = wall_clock()-t_diag;
	
	outfile << "\n\nTime of matrix diagonalization: " << (t_diag)/60.0 << " minutes \n"; 
	
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
