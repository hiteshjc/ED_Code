
#include "read_inputfile.h"
#include "vector_operators.h"
#include "unwrap_pstate.h"

int main(int argc, char *argv[])
{
	string filein, fileEV1, fileEV2, fileout;
 
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
	read_ev1ev2(filein, fileEV1, fileEV2);
	read_Rot(filein, Rot);

	vector<ids> psi;

	ofstream outfile(fileout.c_str());
	outfile << "EV1_file path = " << fileEV1 << endl;
	outfile << "EV2_file path = " << fileEV2 << endl << endl;

	outfile << "Rot = ";
	for(int i = 0; i < Rot.size(); i++)
		outfile << Rot[i] << " ";
	outfile << endl << endl;	


// Read states from vector fiel
	vector<ids> psi1, psi2;
	double t_readev = wall_clock();
	readEV(psi1, fileEV1);
	readEV(psi2, fileEV2);

	t_readev = wall_clock() - t_readev;
	outfile << " Total time for reading vectors  " << (t_readev)/60.0 << " minutes" << endl << endl;

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
