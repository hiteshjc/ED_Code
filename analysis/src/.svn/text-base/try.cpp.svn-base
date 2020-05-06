
#include "read_inputfile.h"
#include "Momentum_Basis.h"
#include "read_vector.h"


int main(int argc, char *argv[])
{
	vector<ids> psi(5);

	psi[0].stateID = 12; psi[0].coeff = 0.04;
	psi[1].stateID = 5; psi[1].coeff = 0.01;
	psi[2].stateID = 10; psi[2].coeff = 0.05;
	psi[3].stateID = 0; psi[3].coeff = 0.06;
	psi[4].stateID = -1; psi[4].coeff = 0.08;

	string fileout = "ttry.txt";

	ofstream outfile(fileout.c_str());

	merge_sort_ids(psi,0,4);

	for (int i = 0; i < psi.size(); i++)
		outfile << psi[i].stateID << "\t" << psi[i].coeff << endl;
	
	outfile.close();

	return 0;
}
