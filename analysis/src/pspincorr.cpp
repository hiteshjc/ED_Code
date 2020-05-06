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

	read_pev1ev2(filein, filepEV1, filepEV2);
	read_T1T2(filein, T1, T2);	
	read_N_Sz(filein, N, Sz);
	read_Triangles(filein, tri);
	read_ev_op(filein,ev1,ev2);

	int nones = int(2*Sz);
	nones = (nones + N)/2;

	ofstream outfile(fileout.c_str());
	outfile << "pEV1_file path = " << filepEV1 << endl;
	outfile << "pEV2_file path = " << filepEV2 << endl;

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
	t_read = wall_clock() - t_read;
	outfile << "Time to read eigenvector = " << (t_read)/60.0 << " minutes" << endl << endl;
	outfile << "Size of pEV1: " << pblock1.size() << endl;
	outfile << "Size of pEV2: " << pblock2.size() << endl;

// Unwrap p-states
	double t_unwrap = wall_clock();
	
	vector< vector<ids> > psis; 
	vector<ids>  psi; 
	psis.push_back(psi);
	psis.push_back(psi);
	Unwrap_pState(vcoeff1, pblock1, norm1, px1, py1, T1, T2, psis[0], N, nones, outfile);	
	Unwrap_pState(vcoeff2, pblock2, norm2, px2, py2, T1, T2, psis[1], N, nones, outfile);

	t_unwrap = wall_clock() - t_unwrap;
	outfile << " Total time for unwrapping p-vector  " << (t_unwrap)/60.0 << " minutes" << endl;

	outfile << " Size of psi1 = " << psis[0].size() << endl;
	outfile << " Size of psi2 = " << psis[1].size() << endl;
	outfile.flush();
	
	outfile<<"Calculating observables for ev_op ev1 (ket) = "<<ev1<<" and ev2 (bra) = "<<ev2<<endl;
	outfile.flush();

	//Sz Sz operator
	outfile<<" SzSz "<<endl;
	for (int n=0;n<tri.size();n++)
	{
			complex<double> sz_sz_exp1=complex<double> (0.0,0.0);
			complex<double> sz_sz_exp2=complex<double> (0.0,0.0);
			complex<double> sz_sz_exp3=complex<double> (0.0,0.0);
			double sr1=0.0,si1=0.0;
			double sr2=0.0,si2=0.0;
			double sr3=0.0,si3=0.0;
			#pragma omp parallel for default(shared) reduction(+:sr1, si1,sr2,si2,sr3,si3)
			for (int64_t n1=0;n1<psis[0].size();n1++)
			{
			complex<double> tmp_szsz1,tmp_szsz2,tmp_szsz3;
			tmp_szsz1 = szsz_op(tri[n].i,tri[n].j, n1, psis[ev1-1], psis[ev2-1], outfile);
			tmp_szsz2 = szsz_op(tri[n].j,tri[n].k, n1, psis[ev1-1], psis[ev2-1], outfile);
			tmp_szsz3 = szsz_op(tri[n].i,tri[n].k, n1, psis[ev1-1], psis[ev2-1], outfile);
			sr1+=real(tmp_szsz1);
			si1+=imag(tmp_szsz1);
			sr2+=real(tmp_szsz2);
			si2+=imag(tmp_szsz2);
			sr3+=real(tmp_szsz3);
			si3+=imag(tmp_szsz3);
			}
			sz_sz_exp1=complex<double>(sr1,si1);
			sz_sz_exp2=complex<double>(sr2,si2);
			sz_sz_exp3=complex<double>(sr3,si3);
			outfile<<"Sz ("<<tri[n].i<<") Sz ("<<tri[n].j<<") = "<<sz_sz_exp1<<endl;
			outfile<<"Sz ("<<tri[n].j<<") Sz ("<<tri[n].k<<") = "<<sz_sz_exp2<<endl;
			outfile<<"Sz ("<<tri[n].i<<") Sz ("<<tri[n].k<<") = "<<sz_sz_exp3<<endl;
			outfile.flush();
	}
	
	//SxSx+SySy = (S+ S- + S-S+)/2 operator
	outfile<<" SxSx+SySy "<<endl;
	for (int n=0;n<tri.size();n++)
	{
			complex<double> sz_sz_exp1=complex<double> (0.0,0.0);
			complex<double> sz_sz_exp2=complex<double> (0.0,0.0);
			complex<double> sz_sz_exp3=complex<double> (0.0,0.0);
			double sr1=0.0,si1=0.0;
			double sr2=0.0,si2=0.0;
			double sr3=0.0,si3=0.0;
			#pragma omp parallel for default(shared) reduction(+:sr1, si1,sr2,si2,sr3,si3)
			for (int64_t n1=0;n1<psis[0].size();n1++)
			{
			complex<double> tmp_szsz1,tmp_szsz2,tmp_szsz3;
			tmp_szsz1 = sxsxsysy_op(tri[n].i,tri[n].j, n1, psis[ev1-1], psis[ev2-1], outfile);
			tmp_szsz2 = sxsxsysy_op(tri[n].j,tri[n].k, n1, psis[ev1-1], psis[ev2-1], outfile);
			tmp_szsz3 = sxsxsysy_op(tri[n].i,tri[n].k, n1, psis[ev1-1], psis[ev2-1], outfile);
			sr1+=real(tmp_szsz1);
			si1+=imag(tmp_szsz1);
			sr2+=real(tmp_szsz2);
			si2+=imag(tmp_szsz2);
			sr3+=real(tmp_szsz3);
			si3+=imag(tmp_szsz3);
			}
			sz_sz_exp1=complex<double>(sr1,si1);
			sz_sz_exp2=complex<double>(sr2,si2);
			sz_sz_exp3=complex<double>(sr3,si3);
			outfile<<"Sp ("<<tri[n].i<<") Sm ("<<tri[n].j<<") = "<<sz_sz_exp1<<endl;
			outfile<<"Sp ("<<tri[n].j<<") Sm ("<<tri[n].k<<") = "<<sz_sz_exp2<<endl;
			outfile<<"Sp ("<<tri[n].i<<") Sm ("<<tri[n].k<<") = "<<sz_sz_exp3<<endl;
			outfile.flush();
	}
	
	
	ts_ch = wall_clock() - ts_ch;

	outfile << "\n Total time for Sz-Sz operator = " << (ts_ch)/60.0 << " minutes\n";
	outfile.flush();


	outfile.close();

	return 0;
}
