
#include "global.h"
#include "idsClass.h"
#include "Entropy.h"
#include "read_inputfile.h"
#include "Matrix_Functions.h"
#include "Momentum_Basis.h"

double find_phase(	vector< complex<double> > &vec,
			vector<int> &T1,
			int L1,
			vector<int64_t> &dets)
{
	int N = vec.size();
	vector< complex<double> > tmp_v(N);

	for(int i = 0; i < N; i++)
	{
		int64_t new_id = dets[i];
		translateT(new_id, T1);

		int64_t loc = findstate_basis(new_id, dets);

		tmp_v[loc] = vec[i];
	}

	complex<double> sum = 0.0;

	for(int i = 0; i < N; i++)
		sum += conj(tmp_v[i])*vec[i];

/*
	if(vec[0] == complex<double>(0.0, 0.0) )
		return 0;

	complex<double> trans_eig = vec[loc]/vec[0]; 

	if(trans_eig != 0.0)
		trans_eig = log(trans_eig);
*/
//	double ph = imag(trans_eig);

	sum = log(sum);

	double ph = imag(sum);
	double kx = ph*L1/2/M_PI;

	if(kx < 0.0)
		kx += L1;

	if( abs(kx - L1) < 1e-2)
		kx = 0.0;

	return abs(kx);
}


int main(int argc, char *argv[])
{
	string filein, fileout, file_dm;

	if (argc <= 1)
	{
		cout << "Usage: " << argv[0] << " <Filename>" << endl;
		exit(1);
	}

	filein=argv[1];
	fileout=argv[2];

	// Read two Cuts along non-trivial directions from file
	vector<int> Cut, T1, T2;
	read_EVfile(filein, file_dm);
	read_T1T2(filein, T1, T2);
	read_Cut(filein, Cut);

	ofstream outfile(fileout.c_str());

	outfile << " Cut = ";
	for(int64_t i = 0; i < Cut.size(); i++)
		outfile << Cut[i] << " ";

	outfile << endl;
	
	outfile << " T1 = ";
	for(int64_t i = 0; i < T1.size(); i++)
		outfile << T1[i] << " ";
	outfile << endl;

	outfile << " T2 = ";
	for(int64_t i = 0; i < T2.size(); i++)
		outfile << T2[i] << " ";
	outfile << endl;
	
	outfile << file_dm << endl;

	double ts_dm = wall_clock();

	int Nvals = 50;

	int L1, L2;
	LxLy(L1, L2, T1, T2);
	
	for(int ones = 0; ones <= Cut.size(); ones++)
	{
		int64_t n_dm = n_choose_k(Cut.size(),ones);
		outfile << "\n Size of density matrix = " << n_dm << " x " << n_dm << endl;

		string sones;
		ostringstream convert; convert << ones; 
		sones = convert.str(); //= string(tmp_ones);

		string tmp_file_dm = file_dm + sones + ".txt";
		outfile << "Num of ones " << ones << endl;

		zMatrix dm(n_dm, n_dm);	

		int count = 0;
		read_density_matrix(dm, tmp_file_dm, count);	
		outfile << "Num of non-zero elements read for dm = " << count << endl;
		
		outfile << endl;
		/////////////////////////////////////////////////////////////////////

		outfile << tmp_file_dm << endl;
		
		double t_diag = wall_clock();
		
		vector<double> Eigs(n_dm);
		zMatrix Evecs(n_dm, n_dm);
		zMatrix_Diagonalize_Hermitian(dm, Eigs, Evecs);
		t_diag = wall_clock() - t_diag;
		outfile << "\n Time to diagonalize matrix = " << t_diag/(60.0) << endl;

		int Nmax = n_dm;

		if(Nmax > Nvals)
			Nmax = Nvals;

		outfile << "Eigenvalues\tSz\t-log(Eig)\tkx\tcol\n\n";

		vector<int64_t> dets;
		constrained_dets(Cut.size(), ones, dets);

		int dSz = 2*ones - Cut.size();
		double Sz = dSz/2.0; 

		std::vector<double> kx_vals;
		double kx0 = 0;
		double Ent0 = 10e10;

		for(int i = n_dm-1 ; i >= n_dm - Nmax; i--)
		{
			vector< complex<double> > tmp_vec(n_dm);
			
			for(int j = 0; j < n_dm; j++)
				tmp_vec[j] = Evecs(j, i);

			double kx = find_phase(tmp_vec, T1, L1, dets);
			if( (-log(Eigs[i])) < Ent0)
			{
				kx0 = kx;
				Ent0 = -log(Eigs[i]);
			}
		}

		for(int i = n_dm-1 ; i >= n_dm - Nmax; i--)
		{
			vector< complex<double> > tmp_vec(n_dm);
			
			for(int j = 0; j < n_dm; j++)
				tmp_vec[j] = Evecs(j, i);

			double kx = find_phase(tmp_vec, T1, L1, dets) - kx0;

			if( kx > 0.01)
				kx = kx - L1;
		
			int col = 0;

			outfile << fixed << setprecision(6) << Eigs[i] << "\t" << setprecision(2) << Sz << "\t" << setprecision(6) << -log(Eigs[i]) << "\t" << setprecision(1) << kx << "\t" << setprecision(1) << int(0) << std::endl;
		}
	
	}

	outfile.close();

	return 0;
}
