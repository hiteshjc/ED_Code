
#include "global.h"
#include "idsClass.h"
#include "read_inputfile.h"
#include "Entropy.h"
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

	double c_min, interval, delta;
	int N;
	double Sz;

	vector<int> T1, T2;
	read_T1T2(filein, T1, T2);
	read_pev1ev2(filein, filepEV1, filepEV2);
	read_N_Sz(filein, N, Sz);

	int nones = int(2*Sz);
	nones = (nones + N)/2;

	ofstream outfile(fileout.c_str());

	outfile << " Num sites = " << N << endl;
	outfile << " Sz sector = " << Sz << endl;
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

	vector<ids> psi1, psi2;
	
	double t_unwrap = wall_clock();
	read_pEV(n_p1, px1, py1, pBasis1, norm1, v_coeff1, filepEV1);

	outfile << " size of wrapped states " << n_p1 << "\t" << px1 <<"\t" << py1 << "\t" << pBasis1.size() << "\t" << norm1.size() << "\t" << v_coeff1.size() << endl;

	Unwrap_pState(v_coeff1, pBasis1, norm1, px1, py1, T1, T2, psi1, N, nones, outfile);
	t_unwrap = wall_clock() - t_unwrap;
	outfile << "Time to read and unwrap p-state1 = " << t_unwrap/(60.0) << " minutes" << endl; 
	outfile << "\n Size of unwrapped EV1: " << psi1.size() << endl;
	v_coeff1.clear(); pBasis1.clear(); norm1.clear();

	read_pEV(n_p2, px2, py2, pBasis2, norm2, v_coeff2, filepEV2);
	outfile << " size of wrapped states " << n_p2 << "\t" << px2 <<"\t" << py2 << "\t" << pBasis2.size() << "\t" << norm2.size() << "\t" << v_coeff2.size() << endl;

	Unwrap_pState(v_coeff2, pBasis2, norm2, px2, py2, T1, T2, psi2, N, nones, outfile);
	t_unwrap = wall_clock() - t_unwrap;
	outfile << "Time to read and unwrap p-state2 = " << t_unwrap/(60.0) << " minutes" << endl; 
	outfile << " Size of unwrapped EV2: " << psi2.size() << endl << endl;
	v_coeff2.clear(); pBasis2.clear(); norm2.clear();

	double t_dp = wall_clock();
	outfile << "\n Dot product of psi1 and psi1 " << dot_product_ids(psi1, psi1) << endl;
	outfile << "\n Dot product of psi1 and psi2 " << dot_product_ids(psi1, psi2) << endl;
	outfile << "\n Dot product of psi2 and psi2 " << dot_product_ids(psi2, psi2) << endl;
	t_dp = wall_clock() - t_dp;
	outfile << "\n Total time to compute dot products " << (t_dp)/60.0 << " minutes\n\n\n";
	outfile.close();

	return 0;
}
