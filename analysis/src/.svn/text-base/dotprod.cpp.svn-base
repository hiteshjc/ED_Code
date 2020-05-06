
#include "global.h"
#include "idsClass.h"
#include "mes_readfile.h"
#include "Entropy.h"
#include "read_inputfile.h"

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

	double c_min, interval, delta;

	// Read two Cuts along non-trivial directions from file
	vector<int> Cut;
	int alpha=2;
	read_mes(filein, Cut, fileEV1, fileEV2, alpha, c_min, interval, delta);
		
	ofstream outfile(fileout.c_str());
	
	outfile << " Vector file path 1: " << fileEV1 << endl;
	outfile << " Vector file path 2: " << fileEV2 << endl << endl;

	vector<ids> psi1, psi2;
	double t_read = wall_clock();
	readEV(psi1, fileEV1);
	readEV(psi2, fileEV2);
	t_read = wall_clock() - t_read;
	outfile << "\nTime to read p-states = " << t_read/(60.0) << " minutes" << endl; 

	outfile << "\n Size of psi1 = " << psi1.size();
	outfile << "\n Size of psi2 = " << psi2.size() << endl;

	double t_dp = wall_clock();
	outfile << "\n Dot product of psi1 and psi1 " << dot_product_ids(psi1, psi1) << endl;
	outfile << "\n Dot product of psi1 and psi2 " << dot_product_ids(psi1, psi2) << endl;
	outfile << "\n Dot product of psi2 and psi2 " << dot_product_ids(psi2, psi2) << endl;
	t_dp = wall_clock() - t_dp;
	outfile << "\n Total time to compute dot products " << (t_dp)/60.0 << " minutes\n\n\n";
	outfile.close();

	return 0;
}
