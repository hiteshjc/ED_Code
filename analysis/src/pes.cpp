
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

	double c_min = 0.0, interval = 1.0;	// default values
	double delta = 0.01;

	// Read two Cuts along non-trivial directions from file
	vector<int> Cut, T1, T2;
	int alpha=2;
	read_Cut(filein, Cut);
	read_pev1(filein, filepEV1);
	read_T1T2(filein, T1, T2);
	int N; double Sz;
	read_N_Sz(filein, N, Sz);
	int nones = int(2*Sz);
	nones = (nones+N)/2;

	int ntot = 1/delta;
	double c_max = c_min + interval;

	ofstream outfile(fileout.c_str());
	
	outfile << " Vector file path 1: " << filepEV1 << endl;
	outfile << "\nT1=";
	for(int i = 0; i < T1.size(); i++)
		outfile << T1[i] << " " ;
	outfile << "\nT2=";
	for(int i = 0; i < T2.size(); i++)
		outfile << T2[i] << " " ;
	outfile << endl;

	double t_unwrap = wall_clock();
	
	vector<ids> psi1;
	read_n_unwrap(filein, filepEV1, T1, T2, psi1, outfile);
	t_unwrap = wall_clock() - t_unwrap;
	outfile << "Time to unwrap p-states = " << t_unwrap/(60.0) << " minutes" << endl; 
	outfile << "\n Size of unwrapped EV1: " << psi1.size() << endl;
	// Sizes of reduced density matrices along two different cuts
	int64_t n_dm = int64_t_power(2, Cut.size());
	outfile << "\n Size of density matrix = " << n_dm << " x " << n_dm << endl;
	// Initialize and compute the density matrices for the two vectors
	double t_dm = wall_clock();
	zMatrix	dm11(n_dm, n_dm, 0.0);
	density_matrix(dm11, psi1, Cut);
	zMatrix	VR(n_dm, n_dm, 0.0);
	vector< complex<double> > lambdai(n_dm);
	zMatrix_Diagonalize(dm11, lambdai, VR);
	for (int i=0;i<n_dm;i++)
	{
		outfile<<real(lambdai[i])<<endl;
	}
	outfile.close();

	return 0;
}
