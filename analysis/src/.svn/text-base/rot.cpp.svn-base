
#include "read_inputfile.h"
#include "Momentum_Basis.h"


int main(int argc, char *argv[])
{
	string filein, fileEV, fileout;
 
	double t_sp = wall_clock();

	if (argc <= 1)
	{
		cout << "Usage: " << argv[0] << " <Filename>" << endl;
		exit(1);
	}

	filein=argv[1];
	fileout=argv[2];

	
	// Read EV file location
	vector<int> Rot;
	
	read_EVfile(filein, fileEV);
	read_Rot(filein, Rot);

//	vector<int> T1, T2;
//	read_T1T2(filein, T1, T2);	

	vector<ids> psi;

	ofstream outfile(fileout.c_str());
	outfile << "\n EV_file path = " << fileEV;

	outfile << "\n\nRot = ";
	for(int i = 0; i < Rot.size(); i++)
		outfile << Rot[i] << " ";

	outfile << endl << endl;	


	double t_read = wall_clock();
	readEV(psi, fileEV);
	t_read = wall_clock() - t_read; 	
	outfile << "\n\nTime to read eigenvector = " << (t_read)/60.0 << " minutes";

	int64_t nstates = psi.size();
	outfile << "\n\nSize of ket: " << nstates << endl;

	vector<ids> rpsi;
	rpsi.resize(nstates);

	double t_rot = wall_clock();

	#pragma omp parallel for
	for(int64_t i = 0; i < nstates; i++)
	{
		int64_t tmpid = psi[i].stateID;

		int64_t newid = 0;

		for(int64_t j = 0; j < Rot.size(); j++)
		{
			bool bval = btest64(tmpid, Rot[j]);
			if(bval == 1)
				newid = ibset64(newid, j);

		}

		int64_t loc = find_ids_state(psi, newid);

		if (loc >= 0)
		{
			rpsi[loc].stateID = newid;
			rpsi[loc].coeff = psi[i].coeff;
		}
	}

	t_rot = wall_clock() - t_rot;
	outfile << "\n\n Time for rotating vector  " << (t_rot)/60.0 << " minutes\n";
	

	outfile << "\nEIGENVECTOR\n";
	for(int64_t j = 0; j < rpsi.size(); j++)
	{
		outfile << rpsi[j].stateID << "\t" << fixed << setprecision(15) << rpsi[j].coeff << endl;
	}
	
	outfile.close();

	return 0;
}
