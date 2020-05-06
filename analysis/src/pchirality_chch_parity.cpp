#include "read_inputfile.h"
#include "Momentum_Basis.h"
#include "unwrap_pstate.h"
#include "chirality_exp.h"

int main(int argc, char *argv[])
{
	string filein, filepEV1, filepEV2, filepEV3, fileout;
	int ev1,ev2;

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

	read_pev1ev2ev3(filein, filepEV1, filepEV2, filepEV3);
	read_T1T2(filein, T1, T2);	
	read_N_Sz(filein, N, Sz);
	read_Triangles(filein, tri);
	read_ev_op(filein,ev1,ev2);

	int nones = int(2*Sz);
	nones = (nones + N)/2;

	ofstream outfile(fileout.c_str());
	outfile << "pEV1_file path = " << filepEV1 << endl;
	outfile << "pEV2_file path = " << filepEV2 << endl;
	outfile << "pEV3_file path = " << filepEV3 << endl;

	outfile<<"Note EV notation 1,2 = parity eigenstates made out of HJC real vectors and 3 = R0"<<endl;
 
	outfile << "\nChirality of triangles\n";
	for(int m=0; m < tri.size(); m++)
		outfile << "Triangle number "<<m<<": "<<tri[m].i << " " << tri[m].j << " " << tri[m].k << endl;

// Read p-states from vector file
	int64_t np1, np2, np3;
	double px1, px2, px3, py1, py2, py3;
	vector<int64_t> pblock1, pblock2, pblock3;
	vector<double> norm1, norm2, norm3;
	vector< complex<double> > vcoeff1, vcoeff2, vcoeff3;

	double t_read = wall_clock();
	
	read_pEV(np1, px1, py1, pblock1, norm1, vcoeff1, filepEV1);
	read_pEV(np2, px2, py2, pblock2, norm2, vcoeff2, filepEV2);
	read_pEV(np3, px3, py3, pblock3, norm3, vcoeff3, filepEV3);

	t_read = wall_clock() - t_read;
	outfile << "Time to read eigenvector = " << (t_read)/60.0 << " minutes" << endl << endl;

	outfile << "Size of pEV1: " << pblock1.size() << endl;
	outfile << "Size of pEV2: " << pblock2.size() << endl;
	outfile << "Size of pEV3: " << pblock3.size() << endl << endl;

// Unwrap p-states
	double t_unwrap = wall_clock();
	
	vector< vector<ids> > psis; 
	vector<ids>  psi; 
	psis.push_back(psi);
	psis.push_back(psi);
	psis.push_back(psi);
	Unwrap_pState(vcoeff1, pblock1, norm1, px1, py1, T1, T2, psis[0], N, nones, outfile);	
	Unwrap_pState(vcoeff2, pblock2, norm2, px2, py2, T1, T2, psis[1], N, nones, outfile);
	Unwrap_pState(vcoeff3, pblock3, norm3, px3, py3, T1, T2, psis[2], N, nones, outfile);

	t_unwrap = wall_clock() - t_unwrap;
	outfile << " Total time for unwrapping p-vector  " << (t_unwrap)/60.0 << " minutes" << endl;

	outfile << " Size of psi1 = " << psis[0].size() << endl;
	outfile << " Size of psi2 = " << psis[1].size() << endl;
	outfile << " Size of psi3 = " << psis[2].size() << endl << endl;

	outfile<<"Making rotation eigenstate and saving it in psi1"<<endl;
	complex<double> c1=complex<double>(+0.99210124,0.0); 
	complex<double> c2=complex<double>(-0.12543978,0.0);
	complex<double> c3=complex<double>(-0.12543978,0.0);
	complex<double> c4=complex<double>(-0.99210124,0.0); 
	
	complex<double> c3p=c3/c1;
	complex<double> c4p=c4-(c3*c2/c1);

	for(int i = 0; i< psis[0].size(); i++) psis[0][i].coeff = c1*psis[0][i].coeff + c2*psis[1][i].coeff; 	
	for(int i = 0; i< psis[1].size(); i++) psis[1][i].coeff = c3p*psis[0][i].coeff+ c4p*psis[1][i].coeff; 	
	outfile<<"FINISHED Making both parity eigenstates (1=+,2=-)"<<endl;
	outfile.flush();
	
	outfile<<"Calculating observables for ev_op ev1 (ket) = "<<ev1<<" and ev2 (bra) = "<<ev2<<endl;
	outfile.flush();
	for (int n=0;n<1;n++)
	{
		// Check with respect to triangle 0
		for (int m=0;m<tri.size();m++)
		{
			complex<double> ch_ch_exp=complex<double> (0.0,0.0);
			double sr=0.0,si=0.0;
			#pragma omp parallel for default(shared) reduction(+:sr, si)
			for (int64_t n1=0;n1<psis[0].size();n1++)
			{
			complex<double> tmp_ch_ch;
			tmp_ch_ch = ch_ch_op(tri[n].i, tri[n].j, tri[n].k, tri[m].i, tri[m].j, tri[m].k, n1, psis[ev1-1], psis[ev2-1], outfile);
			sr+=real(tmp_ch_ch);
			si+=imag(tmp_ch_ch);
			}
			ch_ch_exp=complex<double>(sr,si);
			outfile<<"Ch ("<<n<<") Ch ("<<m<<") = "<<ch_ch_exp<<endl;
			outfile.flush();
		}
	}
	
	ts_ch = wall_clock() - ts_ch;

	outfile << "\n Total time for Ch-Ch operator = " << (ts_ch)/60.0 << " minutes\n";
	outfile.flush();
	outfile.close();

	return 0;
}
