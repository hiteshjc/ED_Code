
#include "global.h"
#include "idsClass.h"
#include "Entropy.h"
#include "read_inputfile.h"
#include "Momentum_Basis.h"

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
	int n_ones;
	read_Cut(filein, Cut);
	read_ev1ev2(filein, fileEV1, fileEV2);
	read_dm_nones(filein, n_ones);
	
	ofstream outfile(fileout.c_str());

	outfile << " File path EV1 = " << fileEV1 << endl;
	outfile << " File path EV2 = " << fileEV2 << endl << endl;

	outfile << "\nCUT\t";	
	for(int i = 0; i < Cut.size(); i++)
		outfile << Cut[i] << " ";
	outfile << endl;

	vector <ids> psi1, psi2;
	double t_readev = wall_clock();
	readEV(psi1, fileEV1);
	readEV(psi2, fileEV2);
	t_readev = wall_clock() - t_readev;
	outfile << "\nTime to read unwrapped EV's = " << t_readev/(60.0) << " minutes" << endl;

	outfile << "\n\n Dot Product for vectors " << dot_product_ids(psi1,psi2) << endl << endl;
	
	double t_dm = wall_clock();
	int n_sites = Cut.size();

//	for(int n_ones = 0; n_ones <= n_sites; n_ones++)
//	{
		int64_t n_dm = n_choose_k(n_sites, n_ones);
		outfile << "Num states = " << n_dm << endl;
		outfile << "Num_ones = " << n_ones << endl << endl;
		
		zMatrix denm(n_dm, n_dm, 0.0);
		density_matrixij_Sz(denm, psi1, psi2, Cut, n_sites, n_ones, outfile);
		t_dm = wall_clock()-t_dm;
		outfile << "Time to construct density matrices : " << (t_dm)/60.0 << " minutes" << endl << endl; 


	
		outfile << "DENSITY_MATRIX" << endl;
		
		for(int i = 0; i < n_dm; i++)
		{
			for(int j = 0; j < n_dm; j++)
				outfile << fixed << setprecision(24) << denm(i,j) << "\t";
			outfile << endl;
		}	

//	}
	
	outfile.close();

	return 0;
}
