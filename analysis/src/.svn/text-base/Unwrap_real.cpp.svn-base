
#include "global.h"
#include "idsClass.h"
#include "read_inputfile.h"
#include "unwrap_pstate.h"

void reads_pvectors(	string filepEV,
			vector<int64_t> &pblock_states,
			vector<double> &Norm,
			vector< complex <double> > &v_coeff)
{
	int64_t jb = 0;

	string phead("pblock_states 		 Norm 			 v_coeff");

	ifstream infile(filepEV.c_str());

	string line;
	while(getline(infile, line))
	{
		int64_t tmp_state; double tmp_norm; complex<double> tmp_coeff;
 
		if(line.compare(phead) == 0)
		{
			while(infile >> pblock_states[jb] >> Norm[jb] >> tmp_coeff)
			{	
				v_coeff[jb++] = tmp_coeff;
			}
		}
	}
	infile.close();

}



int main(int argc, char *argv[])
{
	double ts_uw = wall_clock();
	
	string filein, filepEV, fileout;

	if (argc <= 1)
	{
		cout << "Usage: " << argv[0] << " <Filename>" << endl;
		exit(1);
	}

	filein=argv[1];
	fileout=argv[2];

	// Read two Cuts along non-trivial directions from file
	vector<int> T1, T2;
	int N; double Sz;

	ofstream outfile(fileout.c_str());

	read_T1T2(filein, T1, T2);
	read_pEV(filein, filepEV);
	read_N_Sz(filein, N, Sz);

	int nones = int(2*Sz);
	nones = (nones + N)/2;

	ifstream infile(filepEV.c_str());
	string line;
	string p1p2("np 		 px 		py");	
	double px, py;
	int64_t n_p;

	while(getline(infile, line))
	{
		if(line.compare(p1p2) == 0)
		{
			infile >> n_p >> px >> py;
			break;
		}
	}	
	infile.close();

// Reads vectors 1 and 2 unwraps
	outfile << px << "\t" << py << "\t" << n_p << endl;	
	vector<int64_t> pblock_states1(n_p);
	vector<double> Norm1(n_p);
	vector< complex <double> > v_coeff1(n_p);
	
	reads_pvectors(filepEV, pblock_states1, Norm1, v_coeff1);
	vector<ids> psi1;
	Unwrap_pState(v_coeff1, pblock_states1, Norm1, px, py, T1, T2, psi1, N, nones, outfile);	
	v_coeff1.clear(); pblock_states1.clear(); Norm1.clear();

// Forms linear combination and outputs
	vector<ids> psi(psi1.size() );

	double rat = imag( psi1[0].coeff)/real( psi1[0].coeff);
	outfile << "\n ratio = " << rat << endl;
	double ph1 = atan( rat );
	outfile << " phase = " << ph1 << endl;
	complex<double> alpha = complex<double> (0.0, -ph1);
	
	for(int64_t jb = 0; jb < psi1.size(); jb++)
	{
		psi[jb].coeff = exp(alpha)*psi1[jb].coeff;
		psi[jb].stateID = psi1[jb].stateID;
	}


	ts_uw = wall_clock() - ts_uw;
	outfile << "\n\n Total time for unwrapping p-vector  " << (ts_uw)/60.0 << " minutes\n";
		
	outfile << "\nEIGENVECTOR\n";
	for(int64_t jb = 0; jb < psi.size(); jb++)
	{
		outfile << psi[jb].stateID << "\t" << fixed << setprecision(15) << psi[jb].coeff << endl;
	}
	
	outfile.close();

	return 0;
}
