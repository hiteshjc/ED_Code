
#include "global.h"
#include "idsClass.h"
#include "read_inputfile.h"
#include "unwrap_pstate.h"
#include "Momentum_Basis.h"

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

	// Read two Cuts along non-trivial directions from file
	int N; double Sz;

	ofstream outfile(fileout.c_str());

	read_pev1ev2(filein, filepEV1, filepEV2);
	read_N_Sz(filein, N, Sz);

	int nones = int(2*Sz);
	nones = (nones + N)/2;

	complex< double > alpha = complex<double> (1.0, 0.0);
	complex< double > beta = complex<double> (1.0, 0.0);

// Reads vectors 1 and 2 unwraps
	vector<ids> psi1, psi2;
	readEV(psi1, filepEV1);
	readEV(psi2, filepEV2);

	double arg1 = atan( imag(psi1[0].coeff) / real(psi1[0].coeff) );
	double arg2 = atan( imag(psi2[0].coeff) / real(psi2[0].coeff) );
	complex<double> rho = complex<double> (0.0, -arg1 - arg2);
	

// Forms linear combination and outputs
	vector<ids> psi(psi1.size() );
	for(int64_t jb = 0; jb < psi1.size(); jb++)
	{
		psi[jb].coeff = alpha*psi1[jb].coeff + exp(rho)*beta*psi2[jb].coeff;
		psi[jb].stateID = psi1[jb].stateID;

	}

	outfile << "\nEIGENVECTOR\n";
	for(int64_t jb = 0; jb < psi.size(); jb++)
	{
		outfile << psi[jb].stateID << "\t" << fixed << setprecision(15) << psi[jb].coeff << endl;
	}
	
	outfile.close();

	return 0;
}
