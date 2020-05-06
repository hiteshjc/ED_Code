#include "read_inputfile.h"
#include "Momentum_Basis.h"
#include "unwrap_pstate.h"
#include "chirality_exp.h"

int main(int argc, char *argv[])
{
	string filein, filepEV1, filepEV2, fileout;

	double ts_ch = wall_clock();

	if (argc <= 1)
	{
		cout << "Usage: " << argv[0] << " <Filename>" << endl;
		exit(1);
	}

	filein=argv[1];
	fileout=argv[2];
	
// Read input file parametrs

	vector<int> T1, T2;
	int N; double Sz;
	vector<triangles> tri; // Read Triangles for computing expectations

	read_pev1ev2(filein, filepEV1, filepEV2);
	read_T1T2(filein, T1, T2);	
	read_N_Sz(filein, N, Sz);
	read_Triangles(filein, tri);

	int nones = int(2*Sz);
	nones = (nones + N)/2;

	vector<ids> psi;

	ofstream outfile(fileout.c_str());
	outfile << "pEV1_file path = " << filepEV1 << endl;
	outfile << "pEV2_file path = " << filepEV2 << endl << endl;

	outfile << "\nChirality of triangles\n";
	for(int m=0; m < tri.size(); m++)
		outfile << tri[m].i << " " << tri[m].j << " " << tri[m].k << endl;

// Read p-states from vector file
	int64_t np1, np2;
	double px1, px2, py1, py2;
	vector<int64_t> pblock1, pblock2;
	vector<double> norm1, norm2;
	vector< complex<double> > vcoeff1, vcoeff2;

	double t_read = wall_clock();
	
	read_pEV(np1, px1, py1, pblock1, norm1, vcoeff1, filepEV1);
	read_pEV(np2, px2, py2, pblock2, norm2, vcoeff2, filepEV2);

	t_read = wall_clock() - t_read;
	outfile << "Time to read eigenvector = " << (t_read)/60.0 << " minutes" << endl << endl;

	outfile << "Size of pEV1: " << pblock1.size() << endl;
	outfile << "Size of pEV2: " << pblock2.size() << endl << endl;

// Unwrap p-states
	double t_unwrap = wall_clock();
	
	vector<ids> psi1, psi2;
	Unwrap_pState(vcoeff1, pblock1, norm1, px1, py1, T1, T2, psi1, N, nones, outfile);	
	Unwrap_pState(vcoeff2, pblock2, norm2, px2, py2, T1, T2, psi2, N, nones, outfile);

	t_unwrap = wall_clock() - t_unwrap;
	outfile << " Total time for unwrapping p-vector  " << (t_unwrap)/60.0 << " minutes" << endl;

	outfile << " Size of psi1 = " << psi1.size() << endl;
	outfile << " Size of psi2 = " << psi2.size() << endl << endl;

	complex<double> ch_exp = complex<double> (0.0, 0.0);
	complex<double> curr_exp = complex<double> (0.0, 0.0);
	double ch_real = 0.0, ch_imag = 0.0;
	double curr_real = 0.0, curr_imag = 0.0;

	for(int n = 0; n < tri.size(); n++)
	{
//		#pragma omp parallel for
//	#pragma omp parallel for default(shared) reduction(+:ch_real, ch_imag, curr_real, curr_imag)
		for(int64_t n1 = 0; n1 < psi1.size(); n1++)
		{
			complex<double> tmp_ch, tmp_curr;

			tmp_ch = ch_op(tri[n].i, tri[n].j, tri[n].k, n1, psi1, psi2, outfile)
				+ ch_op(tri[n].j, tri[n].k, tri[n].i, n1, psi1, psi2, outfile)
				+ ch_op(tri[n].k, tri[n].i, tri[n].j, n1, psi1, psi2, outfile); 

			tmp_curr = curr_op(tri[n].i, tri[n].j, n1, psi1, psi2, outfile)
			 	 + curr_op(tri[n].j, tri[n].k, n1, psi1, psi2, outfile)
				 + curr_op(tri[n].k, tri[n].i, n1, psi1, psi2, outfile); 

			ch_real += real(tmp_ch);
			ch_imag += imag(tmp_ch);
			curr_real += real(tmp_curr);
			curr_imag += imag(tmp_curr);
		}
	}
	
	// Need to multiply by an additional factor of i
	ch_exp = complex<double>(-ch_imag, ch_real);
	curr_exp = complex<double>(-curr_imag, curr_real);

	outfile << "\n\n Expectation value of chirality operator is = " << ch_exp  << endl;	
	outfile << "\n\n Expectation value for current operator is = " << curr_exp << endl;	

	ts_ch = wall_clock() - ts_ch;

	outfile << "\n Total time for chirality operator = " << (ts_ch)/60.0 << " minutes\n";
	outfile.close();

	return 0;
}
