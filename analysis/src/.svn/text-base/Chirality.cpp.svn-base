#include "chirality_exp.h"
#include "read_inputfile.h"

int main(int argc, char *argv[])
{
	string filein, fileEV1, fileEV2, fileout;

	double ts_ch = wall_clock();

	if (argc <= 1)
	{
		cout << "Usage: " << argv[0] << " <Filename>" << endl;
		exit(1);
	}

	filein=argv[1];
	fileout=argv[2];

	// Read triangles and eigenvector file
	vector<triangles> tri;
	read_Triangles(filein, tri);
	read_ev1ev2(filein, fileEV1, fileEV2);

	vector<ids> psi1, psi2;
	readEV(psi1, fileEV1);
	readEV(psi2, fileEV2);
	
	ofstream outfile(fileout.c_str());

	outfile << "\n EV1_file path = " << fileEV1;
	outfile << "\n EV2_file path = " << fileEV2;
	outfile << endl;

	outfile << "\nChirality of triangles\n";
	for(int m=0; m < tri.size(); m++)
		outfile << tri[m].i << " " << tri[m].j << " " << tri[m].k << endl;

	outfile << "\nSize of ket: " << psi1.size() << endl;
	outfile << "\nSize of ket: " << psi1.size() << endl;
	
	complex<double> ch_exp = complex<double> (0.0, 0.0);
	complex<double> curr = complex<double> (0.0, 0.0);

	for(int n = 0; n < tri.size(); n++)
	{
//		#pragma omp parallel for
		for(int64_t n1 = 0; n1 < psi1.size(); n1++)
		{
			ch_exp += ch_op(tri[n].i, tri[n].j, tri[n].k, n1, psi1, psi2, outfile)
				+ ch_op(tri[n].j, tri[n].k, tri[n].i, n1, psi1, psi2, outfile)
				+ ch_op(tri[n].k, tri[n].i, tri[n].j, n1, psi1, psi2, outfile); 

			curr += curr_op(tri[n].i, tri[n].j, n1, psi1, psi2, outfile)
				+ curr_op(tri[n].j, tri[n].k, n1, psi1, psi2, outfile)
				+ curr_op(tri[n].k, tri[n].i, n1, psi1, psi2, outfile); 

		}
	}

	outfile << "\n\n Expectation value of chirality operator is = i" << ch_exp << endl;	
	outfile << "\n\n Expectation value for current operator is = i" << curr << endl;	

	ts_ch = wall_clock() - ts_ch;

	outfile << "\n Total time for chirality operator = " << (ts_ch)/60.0 << " minutes\n";
	outfile.close();

	return 0;
}
