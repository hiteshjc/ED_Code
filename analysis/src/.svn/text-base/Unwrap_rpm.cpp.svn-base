
#include "global.h"
#include "idsClass.h"
#include "read_inputfile.h"
#include "unwrap_pstate.h"
#include "Momentum_Basis.h"

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


void rotate_vector(	vector<ids> &psi,
			vector<int> &Rot,
			vector<ids> &rpsi)
{

	int64_t nstates = psi.size();
	rpsi.resize(nstates);

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

}




int main(int argc, char *argv[])
{
	double ts_uw = wall_clock();
	
	string filein, filepEV1, filepEV2, fileout;

	if (argc <= 1)
	{
		cout << "Usage: " << argv[0] << " <Filename>" << endl;
		exit(1);
	}

	filein=argv[1];
	fileout=argv[2];

	// Read two Cuts along non-trivial directions from file
	vector<int> Rot;
	vector<int> T1, T2;
	int N; double Sz;

	string ofile = fileout + ".txt";
	ofstream outfile(ofile.c_str());

	read_Rot(filein, Rot);
//	read_T1T2(filein, T1, T2);
	read_pev1ev2(filein, filepEV1, filepEV2);
	read_N_Sz(filein, N, Sz);

	int nones = int(2*Sz);
	nones = (nones + N)/2;

	complex< double > alpha = complex<double> (1.0, 0.0);
	complex< double > beta = complex<double> (1.0, 0.0);

	ifstream infile(filepEV1.c_str());
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
//
	double t_un = wall_clock();
	outfile << px << "\t" << py << "\t" << n_p << endl;	
/*	vector<int64_t> pblock_states1(n_p);
	vector<double> Norm1(n_p);
	vector< complex <double> > v_coeff1(n_p);
	
	reads_pvectors(filepEV1, pblock_states1, Norm1, v_coeff1);
*/

	vector<ids> psi1, psi2;
	read_n_unwrap(filein, filepEV1, T1, T2, psi1, outfile);
	read_n_unwrap(filein, filepEV2, T1, T2, psi2, outfile);


//	Unwrap_pState(v_coeff1, pblock_states1, Norm1, px, py, T1, T2, psi1, N, nones, outfile);	
//	v_coeff1.clear(); pblock_states1.clear(); Norm1.clear();


/*	vector<int64_t> pblock_states2(n_p);
	vector<double> Norm2(n_p);
	vector< complex <double> > v_coeff2(n_p);

	reads_pvectors(filepEV2, pblock_states2, Norm2, v_coeff2);
	vector<ids> psi2;
	Unwrap_pState(v_coeff2, pblock_states2, Norm2, px, py, T1, T2, psi2, N, nones, outfile);	
	v_coeff2.clear(); pblock_states2.clear(); Norm2.clear();

	t_un = wall_clock() - t_un;
	outfile << "Time for unwrapping vectors  " << (t_un)/60.0 << " minutes" << endl << endl ;
*/

// Rotate vectors and compute dot-products
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


	zMatrix dprot(2,2);
	dprot(0,0) = dot_product_ids(psi1, rpsi1);	dprot(0,1) = dot_product_ids(psi1, rpsi2);
	dprot(1,0) = dot_product_ids(psi2, rpsi1);	dprot(1,1) = dot_product_ids(psi2, rpsi2);


	outfile << "--- Psi and Rotated Psi overlaps ---\n";
	outfile << " Dot product of psi1 and rpsi1 " << dprot(0,0) << endl;
	outfile << " Dot product of psi1 and rpsi2 " << dprot(0,1) << endl;
	outfile << " Dot product of psi2 and rpsi1 " << dprot(1,0) << endl; 
	outfile << " Dot product of psi2 and rpsi2 " << dprot(1,1) << endl;

	t_dp = wall_clock() - t_dp;
	outfile << "\n Total time to compute dot products " << (t_dp)/60.0 << " minutes";

// Diagonalize rotation matrix
//
	zMatrix VR(2,2);
	vector<complex <double> > W(2);
	zMatrix_Diagonalize(dprot, W, VR);

	zMatrix VRT(2,2);

	outfile << "\n Rotation matrix components \n";
	for(int ii = 0; ii < 2; ii++)
	{
		for(int jj = 0; jj < 2; jj++)
		{
			if(ii == jj)
				VRT(ii, jj) = VR(ii, jj);
			else	VRT(ii, jj) = VR(jj, ii);
			outfile << VRT(ii, jj) << "\t";
		}	
		outfile << endl;
	}

	double arg1 = atan( imag(psi1[0].coeff) / real(psi1[0].coeff) );
	double arg2 = atan( imag(psi2[0].coeff) / real(psi2[0].coeff) );
	complex<double> rho = complex<double> (0.0, -arg1 - arg2);
	

// Form rotation basis states

	double ts_rb = wall_clock();

	int64_t nstate = psi1.size();
	vector<ids> psi_rp(nstate );
	vector<ids> psi_rm(nstate);
	
	outfile << " \n Vectors sizes " << psi1.size() << "\t" << psi2.size() << endl;

	for(int64_t jb = 0; jb < psi1.size(); jb++)
	{
		psi_rp[jb].coeff = VRT(0,0)*psi1[jb].coeff + VRT(0,1)*psi2[jb].coeff;
		psi_rp[jb].stateID = psi1[jb].stateID;
		psi_rm[jb].coeff = VRT(1,0)*psi1[jb].coeff + VRT(1,1)*psi2[jb].coeff;
		psi_rm[jb].stateID = psi1[jb].stateID;
	}


	outfile << " \n Rotation basis vectors sizes " << psi_rp.size() << "\t" << psi_rm.size() << endl;


	ts_rb = wall_clock() - ts_rb;
	outfile << "\n\n Total time for computing rotated plus minus vectors  " << (ts_rb)/60.0 << " minutes\n";
	outfile.close();


// Output rplus vector
	string splus = fileout + "_plus.txt";
	ofstream oplus(splus.c_str());
	
	oplus << px << "\t" << py << "\t" << n_p << endl;	
		
	oplus << "\nEIGENVECTOR\n";
	for(int64_t jb = 0; jb < psi_rp.size(); jb++)
	{
		oplus << psi_rp[jb].stateID << "\t" << fixed << setprecision(15) << psi_rp[jb].coeff << endl;
	}
	
	oplus.close();


// Output rminus vector
	string sminus = fileout + "_minus.txt";
	ofstream ominus(sminus.c_str());
	
	ominus << px << "\t" << py << "\t" << n_p << endl;	
		
	ominus << "\nEIGENVECTOR\n";
	for(int64_t jb = 0; jb < psi_rm.size(); jb++)
	{
		ominus << psi_rm[jb].stateID << "\t" << fixed << setprecision(15) << psi_rm[jb].coeff << endl;
	}
	
	ominus.close();


	return 0;
}
