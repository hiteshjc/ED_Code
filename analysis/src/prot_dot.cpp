
#include "read_inputfile.h"
#include "unwrap_pstate.h"
#include "vector_operators.h"

int main(int argc, char *argv[])
{
	string filein, filepEV1, filepEV2, fileout;
 
	double t_sp = wall_clock();

	if (argc <= 1)
	{
		cout << "Usage: " << argv[0] << " <Filename>" << endl;
		exit(1);
	}

	filein=argv[1];
	fileout=argv[2];

	
// Read input file parametrs

	vector<int> Rot;
	vector<int> T1, T2;
	int N; double Sz;

	read_pev1ev2(filein, filepEV1, filepEV2);
	read_Rot(filein, Rot);
	read_T1T2(filein, T1, T2);	
	read_N_Sz(filein, N, Sz);

	int nones = int(2*Sz);
	nones = (nones + N)/2;

	vector<ids> psi;

	ofstream outfile(fileout.c_str());
	outfile << "pEV1_file path = " << filepEV1 << endl;
	outfile << "pEV2_file path = " << filepEV2 << endl << endl;

	outfile << "Rot = ";
	for(int i = 0; i < Rot.size(); i++)
		outfile << Rot[i] << " ";
	outfile << endl << endl;	



// Read p-states and unwrap from vector file
	double t_unwrap = wall_clock();
	
	vector<ids> psi1, psi2;
	read_n_unwrap(filein, filepEV1, T1, T2, psi1, outfile);
	read_n_unwrap(filein, filepEV2, T1, T2, psi2, outfile);

	t_unwrap = wall_clock() - t_unwrap;
	outfile << " Total time for unwrapping p-vector  " << (t_unwrap)/60.0 << " minutes" << endl << endl;

	outfile << " Size of psi1 = " << psi1.size() << endl;
	outfile << " Size of psi2 = " << psi2.size() << endl << endl;

// Rotate vectors
	double t_rot = wall_clock();

	vector<ids> rpsi1, rpsi2;
	rotate_vector(psi1, Rot, rpsi1);
	rotate_vector(psi2, Rot, rpsi2);

	t_rot = wall_clock() - t_rot;
	outfile << "Time for rotating vector  " << (t_rot)/60.0 << " minutes" << endl << endl ;


// Compute dot products
	double t_dp = wall_clock();
	
	outfile << "--- Psi overlaps ---\n";
	outfile << " Dot product of psi1 and psi1 " << dot_product_ids(psi1, psi1) << endl;
	outfile << " Dot product of psi1 and psi2 " << dot_product_ids(psi1, psi2) << endl;
	outfile << " Dot product of psi2 and psi2 " << dot_product_ids(psi2, psi2) << endl << endl;

	outfile << "--- Rotated Psi overlaps ---\n";
	outfile << " Dot product of rpsi1 and rpsi1 " << dot_product_ids(rpsi1, rpsi1) << endl;
	outfile << " Dot product of rpsi1 and rpsi2 " << dot_product_ids(rpsi1, rpsi2) << endl;
	outfile << " Dot product of rpsi2 and rpsi2 " << dot_product_ids(rpsi2, rpsi2) << endl << endl;

	outfile << "--- Psi and Rotated Psi overlaps ---\n";
	outfile << " Dot product of psi1 and rpsi1 " << dot_product_ids(psi1, rpsi1) << endl;
	outfile << " Dot product of psi1 and rpsi2 " << dot_product_ids(psi1, rpsi2) << endl;
	outfile << " Dot product of psi2 and rpsi1 " << dot_product_ids(psi2, rpsi1) << endl;
	outfile << " Dot product of psi2 and rpsi2 " << dot_product_ids(psi2, rpsi2) << endl << endl;

	t_dp = wall_clock() - t_dp;
	outfile << " Total time to compute dot products " << (t_dp)/60.0 << " minutes";

	outfile.close();

	return 0;
}
