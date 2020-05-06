
#include "global.h"
#include "idsClass.h"
#include "Entropy.h"
#include "read_inputfile.h"
#include "unwrap_pstate.h"

template <class T>
std::string to_string (T const &t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}
template std::string to_string<int>(int const &);
template std::string  to_string<double>(double const &);


int main(int argc, char *argv[])
{
	string filein, filepEV, fileout;

	if (argc <= 1)
	{
		cout << "Usage: " << argv[0] << " <Filename>" << endl;
		exit(1);
	}

	filein=argv[1];
	fileout=argv[2];

	// Read two Cuts along non-trivial directions from file
	vector<int> Cut, T1, T2;
	read_pevfamily(filein, filepEV);
	read_T1T2(filein, T1, T2);
	read_Cut(filein, Cut);
	int N; double Sz;
	read_N_Sz(filein, N, Sz);
	int nones = int(2*Sz);
	nones = (nones+N)/2;
	ofstream outfile(fileout.c_str());
	
	outfile << " Vector file path: " << filepEV << endl;

	outfile << "\nT1=";
	for(int i = 0; i < T1.size(); i++)
		outfile << T1[i] << " " ;
	outfile << "\nT2=";
	for(int i = 0; i < T2.size(); i++)
		outfile << T2[i] << " " ;
	outfile << endl;

	//double t_unwrap = wall_clock();
	
	outfile << "\n Cut \t";
	for(int i = 0; i < Cut.size(); i++)
		outfile << Cut[i] << " ";

	outfile << endl;
	//return 0;

	int hilbert_space=10000;

	for (int i=0;i<hilbert_space;i++)
	{
		vector<ids> psi;
		read_n_unwrap(filein, filepEV+to_string(i)+".txt", T1, T2, psi, outfile);
	
		//t_unwrap = wall_clock() - t_unwrap;
		//outfile << "Time to unwrap p-states = " << t_unwrap/(60.0) << " minutes" << endl; 
		outfile << "\n Size of unwrapped EV: " << psi.size() << endl;

		// Sizes of reduced density matrices along two different cuts
		int64_t n_dm = int64_t_power(2, Cut.size());
	
		outfile << "\n Size of density matrix = " << n_dm << " x " << n_dm << endl;

		// Initialize and compute the density matrices for the two vectors
		zMatrix	dm(n_dm, n_dm, 0.0);
		density_matrix(dm, psi, Cut);
		complex<double> S_von = Von_Neumann_Entropy_basee(dm);
		complex<double> S_ren = Renyi_Entropy_basee(dm);

		outfile << "Von and Renyi Entropies (base e) of state "<< i <<" = "<< real(S_von) << " , " <<real(S_ren)<< endl;
	}
	return 0;
}
