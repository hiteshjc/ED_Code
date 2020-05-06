
#include "global.h"
#include "idsClass.h"
#include "Entropy.h"
#include "Momentum_Basis.h"
#include "read_inputfile.h"
#include "unwrap_pstate.h"


int main(int argc, char *argv[])
{
	string filein, filepEV1, filepEV2, filepEV3, filepEV4, fileout;

	if (argc <= 1)
	{
		cout << "Usage: " << argv[0] << " <Filename>" << endl;
		exit(1);
	}

	filein=argv[1];
	fileout=argv[2];

	vector<int> Cut, T1, T2;
	int N; double Sz; 
	int n_ones;
	int ev1,ev2;

	read_Cut( filein, Cut);
	read_dm_nones(filein, n_ones);
	read_pev1ev2ev3ev4( filein, filepEV1, filepEV2, filepEV3, filepEV4);
	read_ev_op(filein, ev1, ev2);
	read_T1T2(filein, T1, T2);
	read_N_Sz(filein, N, Sz);

	int unwrap_nones = int(2*Sz);
	unwrap_nones = (unwrap_nones + N)/2;


	ofstream outfile(fileout.c_str());

	outfile << " File path pEV1 = " << filepEV1 << endl;
	outfile << " File path pEV2 = " << filepEV2 << endl;
	outfile << " File path pEV3 = " << filepEV3 << endl;
	outfile << " File path pEV4 = " << filepEV4 << endl << endl;

	outfile << " Sites N = " << N << endl;
	outfile << " Sector Sz = " << Sz << endl;
	outfile << " Unwrap_nones = " << unwrap_nones << endl;

	outfile << "\nCUT\t";	
	for(int i = 0; i < Cut.size(); i++)
		outfile << Cut[i] << " ";
	outfile << endl;

	outfile << "\nT1 = ";
	for(int i = 0; i < T1.size(); i++)
		outfile << T1[i] << " ";

	outfile << "\nT2 = ";	
	for(int i = 0; i < T2.size(); i++)
		outfile << T2[i] << " ";
	outfile << endl << endl;	
	
	outfile <<"DM of ev1 = "<<ev1<<" ev2 = "<<ev2<<" being calculated now...."<<endl;
	outfile.flush();

	double t_unwrap = wall_clock();
	
	vector<ids> psi12, psi22;
	vector<ids> psi1, psi2;
	
	read_n_unwrap(filein, filepEV1, T1, T2, psi1, outfile);
	read_n_unwrap(filein, filepEV2, T1, T2, psi12, outfile);
	read_n_unwrap(filein, filepEV3, T1, T2, psi2, outfile);
	read_n_unwrap(filein, filepEV4, T1, T2, psi22, outfile);
	
	t_unwrap = wall_clock() - t_unwrap;

	outfile << "Time to unwrap p-states = " << t_unwrap/(60.0) << " minutes" << endl; 
	outfile << " Size of unwrapped EV1: " << psi1.size() << endl;
	outfile << " Size of unwrapped EV2: " << psi12.size() << endl;
	outfile << " Size of unwrapped EV3: " << psi2.size() << endl;
	outfile << " Size of unwrapped EV4: " << psi22.size() << endl;
	outfile.flush();

//Right vectors
//( -0.184861 +0.000000 j) ( -0.682515 -0.000000 j) ( +0.682515 +0.000000 j) ( -0.184861 +0.000000 j)
//( -0.000008 -0.184869 j) ( +0.000002 +0.682513 j) ( -0.000002 +0.682513 j) ( -0.000008 +0.184869 j)
//( -0.108233 +0.673877 j) ( +0.029315 +0.182528 j) ( -0.029315 +0.182528 j) ( -0.108233 -0.673877 j)
//( -0.673878 -0.108230 j) ( +0.182522 -0.029324 j) ( -0.182522 -0.029324 j) ( -0.673878 +0.108230 j)


// Wrong vectors I mixed up bra-ket when diagonalizing chirality - Eigenvectors of chirality in the 1,2,4,5 space
//( -0.184861 +0.000000 j) ( -0.682515 -0.000000 j) ( +0.682515 +0.000000 j) ( -0.184861 +0.000000 j)
//( -0.000008 +0.184869 j) ( +0.000002 -0.682513 j) ( -0.000002 -0.682513 j) ( -0.000008 -0.184869 j)
//( -0.108233 -0.673877 j) ( +0.029315 -0.182528 j) ( -0.029315 -0.182528 j) ( -0.108233 +0.673877 j)
//( -0.673878 +0.108230 j) ( +0.182522 +0.029324 j) ( -0.182522 +0.029324 j) ( -0.673878 -0.108230 j)

// Eigenvalues of chirality
//[-0.01214185 -0.00997014  0.00997014  0.01214185]

// Basically - my bra-ket messup just changes sign of chirality - nothing else

	complex<double> c1=complex<double>(+0.682515,+0.000000); 
	complex<double> c2=complex<double>(-0.000002,-0.682513); 
	complex<double> c3=complex<double>(-0.029315,-0.182528); 
	complex<double> c4=complex<double>(-0.182522,+0.029324); 
	complex<double> c5=complex<double>(-0.184861,+0.000000); 
	complex<double> c6=complex<double>(-0.000008,-0.184869); 
	complex<double> c7=complex<double>(-0.108233,+0.673877); 
	complex<double> c8=complex<double>(-0.673878,-0.108230); 

	complex<double> c5p=c5/c1;
	complex<double> c6p=c6 - (c2*c5/c1);
	complex<double> c7p=c7 - (c3*c5/c1);
	complex<double> c8p=c8 - (c4*c5/c1);

	for(int i = 0; i< psi1.size(); i++) psi1[i].coeff = c1*psi1[i].coeff  + c2*psi12[i].coeff  + c3*psi2[i].coeff  + c4*psi22[i].coeff; 	
	for(int i = 0; i< psi2.size(); i++) psi2[i].coeff = c5p*psi1[i].coeff + c6p*psi12[i].coeff + c7p*psi2[i].coeff + c8p*psi22[i].coeff; 	
	
	outfile << "\n Got both psi's of the same (positive) chirality " << endl;

	if(psi1.size() != psi2.size())
	outfile << "\n ERROR IN UNWRAPPING VECTORS...SIZES DON'T MATCH \n" << endl;

	double t_readev = wall_clock();
	outfile << "\n\n Dot Product for vectors 1 and 1 (both positive chirality)" << dot_product_ids(psi1,psi1) << endl << endl;
	outfile << "\n\n Dot Product for vectors 1 and 2 (both positive chirality)" << dot_product_ids(psi1,psi2) << endl << endl;
	outfile << "\n\n Dot Product for vectors 2 and 2 (both positive chirality)" << dot_product_ids(psi2,psi2) << endl << endl;
	t_readev = wall_clock() - t_readev;

	outfile << "\nTime to take dot product of unwrapped vectors = " << t_readev/(60.0) << " minutes" << endl;
	
	outfile <<"DM of ev1 = "<<ev1<<" ev2 = "<<ev2<<" being calculated now...."<<endl;

	double t_dm = wall_clock();
	int n_sites = Cut.size();

//	for(int n_ones = 0; n_ones <= n_sites; n_ones++)
//	{
		int64_t n_dm = n_choose_k(n_sites, n_ones);
		outfile << "Num states = " << n_dm << endl;
		outfile << "Num_ones = " << n_ones << endl << endl;
		
		zMatrix denm(n_dm, n_dm, 0.0);
		if (ev1==1 and ev2==1) density_matrixij_Sz(denm, psi1, psi1, Cut, n_sites, n_ones, outfile);
		if (ev1==1 and ev2==2) density_matrixij_Sz(denm, psi1, psi2, Cut, n_sites, n_ones, outfile);
		if (ev1==2 and ev2==1) density_matrixij_Sz(denm, psi2, psi1, Cut, n_sites, n_ones, outfile);
		if (ev1==2 and ev2==2) density_matrixij_Sz(denm, psi2, psi2, Cut, n_sites, n_ones, outfile);
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
