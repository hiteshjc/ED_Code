
#include "Lanczos_sector.h"

template <class T>
std::string to_string (T const &t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}
template std::string to_string<int>(int const &);
template std::string  to_string<double>(double const &);


void sector_lanczos(	int &kx,
			int &ky,
			double &cSz,
			int &n_spins,
			double &lambda,
			double &J2,
			vector<coordinates> &adj_list,
			double &Ch,
			vector<triangles> &ch_list,
			double &ChCh,
			vector<bowties> &chch_list,
			vector<int> &T1,
			vector<int> &T2,
			int &n_sec_ev,
			bool &find_ev,
			bool &unwrap_peigvec,
			string &filepath_EV,
			double &tol,
			int &iterations,
			int &ncycles,
			ofstream &outfile)
{

	/*int n_ones = n_spins/2;
	
	for(double i = 0.0; i < cSz; i += 1)
	{	n_ones++;} */

	int n_ones;

        if (cSz>=0) n_ones=int(cSz + n_spins/2+1.0e-6);
        if (cSz<0) n_ones=int(cSz + n_spins/2+1.0e-6);

	int L1 = 0, L2 = 0;
	LxLy(L1, L2, T1, T2);

	if(kx >= L1 || ky >= L2 || kx < 0 || ky < 0)
	{
		outfile << "\nInvalid kx or ky\n";
		exit(1);
	}

	double px = 2*kx*M_PI/L1, py = 2*ky*M_PI/L2;
	vector<int64_t> pblock_states;
	vector<double> norm, norm_inv;
	vector<int64_t> rep_loc;
	vector< complex<double> > state_ph;
	int64_t n_p = 0;
	vector<int64_t> Ia, Ib;
	int bits_right;

	double Sz = 0.5*n_ones - 0.5*(n_spins - n_ones);

	outfile << "\n--------------------------------------------------------------------\n"; 
	double t_states = wall_clock();

        outfile.flush();
	pConstruct_States(	pblock_states, norm, norm_inv, T1, T2, n_ones, px, py, n_p, 
				rep_loc, state_ph, Ia, Ib, bits_right, outfile);
        outfile.flush();
/*
	outfile << endl << endl;
	for(int ii = 0; ii < pblock_states.size(); ii++)
		outfile << pblock_states[ii] << "\t" << norm[ii] << "\t" << state_ph[ii] << endl;
	outfile << endl;

	outfile << endl << endl;
	for(int ii = 0; ii < rep_loc.size(); ii++)
		outfile << ii << "\t" << Ia[ii] << "\t" << Ib[ii] << "\t" << rep_loc[ii] << endl;
	outfile << endl;
*/

	t_states = wall_clock() - t_states;

	outfile << endl << "Total size of all states/rep loc vector " << rep_loc.size() << endl;
	outfile << " SECTOR : " <<  Sz << "\t(" << px << "," << py << ")\t" << n_p << endl;	
	outfile << "Time to construct pStates\t " << t_states/60.0 << " minutes" << endl;
	
	if(n_p == 0)
	{
		outfile << "\nNo states in psector\n";
		exit(1);
	}
	outfile << "------------------------------------------------------------------" << endl;
	
	vector< complex<double> > pblock_eigs;
	vector< vector< complex<double> > > pblock_eigvecs;

	double lanczos_time = wall_clock();

        if (n_p>10000)
	{
		lanczos_evecs(	pblock_states, norm, norm_inv, rep_loc, state_ph,
			pblock_eigs, pblock_eigvecs, find_ev,
			lambda, n_spins, px, py, n_ones, n_p,
			adj_list, Ch, ch_list, ChCh, chch_list, T1, T2, n_sec_ev, tol, iterations, ncycles,
			Ia, Ib, bits_right, outfile);
	}
        else
	{
		exact_evecs(	pblock_states, norm, norm_inv, rep_loc, state_ph,
			pblock_eigs, pblock_eigvecs, find_ev,
			lambda, n_spins, px, py, n_ones, n_p,
			adj_list, Ch, ch_list, ChCh, chch_list, T1, T2, n_sec_ev, tol, iterations, ncycles,
			Ia, Ib, bits_right, outfile);
	}
	
	outfile << "-----------------------------------------------------------" << endl;
	outfile << "-----------------------------------------------------------" << endl;
	outfile << " Wall clock for Lanczos " << (wall_clock() - lanczos_time)/60.0 << endl;


	outfile<<"pblock_eigs.size() = "<<pblock_eigs.size()<<endl;
	outfile  << "\n\nEnergies\t\tEnergies/site\t\tSz\tkx\tky\tJz\tJ2\n\n";
	outfile.flush();
		
	for(int i = 0; i < pblock_eigs.size(); i++)	
		outfile << fixed << setprecision(15) << real(pblock_eigs[i]) << "\t" << real(pblock_eigs[i])/double(n_spins) <<"\t" << setprecision (3) << Sz << "\t" << px << "\t" << py << "\t" << lambda << "\t" << J2 << "\n";  	

	outfile.flush();
	
	if(find_ev)
	{	
		// double ts_unwrap = wall_clock();
		
		for(int iib = 0; iib < min(n_sec_ev,int(pblock_eigs.size())); iib++)
		{
			/*int t_tmp_int = iib + 1 + (int)'0';
			char t_tmp = (char)t_tmp_int;

			string filepath_tmp;

			if(unwrap_peigvec)
				filepath_tmp = filepath_EV + "_EV" + t_tmp + ".txt";		
			else
				filepath_tmp = filepath_EV + "_pEV" + t_tmp + ".txt";*/
	
			string filepath_tmp;

			if(unwrap_peigvec)
				filepath_tmp = filepath_EV + "_EV" + to_string(iib) + ".txt";		
			else
				filepath_tmp = filepath_EV + "_pEV" + to_string(iib) + ".txt";
	
			ofstream outfile_vecs(filepath_tmp.c_str());
			outfile_vecs << "\n Eigenvector for state \n";
			outfile_vecs << "Sz\t(kx,ky)\tJz/lambda\n";
			outfile_vecs << Sz <<"\t(" << kx << "," << ky << ")\t" << lambda << endl;
		
			outfile_vecs << "\n\n-------------------------------------------";
			outfile_vecs << "\nState Energy " << pblock_eigs[iib] << endl;

			if(unwrap_peigvec)
			{
				vector< complex <double> > v_coeff = pblock_eigvecs[iib];

				vector<ids> psi = Unwrap_pState(v_coeff, pblock_states, norm, px, py, T1, T2);	
			
				outfile_vecs << "\nEIGENVECTOR\n";
				for(int64_t jb = 0; jb < psi.size(); jb++)
				{
					outfile_vecs << psi[jb].stateID << "\t" << fixed << setprecision(15) << psi[jb].coeff << endl;
				}
			
			}
			else
			{
				outfile_vecs << "\nMomentum EIGENVECTOR\n";
			
				outfile_vecs << "\nnp \t\t px \t\tpy\n";
				outfile_vecs << n_p << "\t\t" << px << "\t\t" << py << endl;

				outfile_vecs << "\n\npblock_states \t\t Norm \t\t\t v_coeff\n\n";
				
				for(int64_t jb = 0; jb < n_p; jb++)
					outfile_vecs << pblock_states[jb] << "\t\t\t" << norm[jb] << "\t\t\t" << fixed << setprecision(15) << pblock_eigvecs[iib][jb]<< endl;

			}
			outfile_vecs.close();
		}
		// double te_unwrap = wall_clock();
		// outfile << "\nTime to unwrap state = " << (te_unwrap - ts_unwrap)/60 << " minutes\n";			
	}
}


void sector_lanczos(	int &kx,
			int &ky,
			int &kz,
			double &cSz,
			int &n_spins,
			double &lambda,
			double &J2,
			vector<coordinates> &adj_list,
			double &Ch,
			vector<triangles> &ch_list,
			double &ChCh,
			vector<bowties> &chch_list,
			vector<int> &T1,
			vector<int> &T2,
			vector<int> &T3,
			int &n_sec_ev,
			bool &find_ev,
			bool &unwrap_peigvec,
			string &filepath_EV,
			double &tol,
			int &iterations,
			int &ncycles,
			ofstream &outfile)
{

	int n_ones = n_spins/2;
	
	for(double i = 0.0; i < cSz; i += 1)
	{	n_ones++;} 

	int L1 = 0, L2 = 0, L3 = 0 ;
	LxLyLz(L1, L2, L3, T1, T2, T3);

	if(kx >= L1 || ky >= L2 || kz>=L3 || kx < 0 || ky < 0 || kz<0 )
	{
		outfile << "\nInvalid kx or ky\n";
		exit(1);
	}

	double px = 2*kx*M_PI/L1, py = 2*ky*M_PI/L2 , pz=2*kz*M_PI/L3 ;
	vector<int64_t> pblock_states;
	vector<double> norm, norm_inv;
	vector<int64_t> rep_loc;
	vector< complex<double> > state_ph;
	int64_t n_p = 0;
	vector<int64_t> Ia, Ib;
	int bits_right;

	double Sz = 0.5*n_ones - 0.5*(n_spins - n_ones);

	outfile << "\n--------------------------------------------------------------------\n"; 
	double t_states = wall_clock();

	pConstruct_States(	pblock_states, norm, norm_inv, T1, T2, T3, n_ones, px, py, pz, n_p, 
				rep_loc, state_ph, Ia, Ib, bits_right, outfile);
/*
	outfile << endl << endl;
	for(int ii = 0; ii < pblock_states.size(); ii++)
		outfile << pblock_states[ii] << "\t" << norm[ii] << "\t" << state_ph[ii] << endl;
	outfile << endl;

	outfile << endl << endl;
	for(int ii = 0; ii < rep_loc.size(); ii++)
		outfile << ii << "\t" << Ia[ii] << "\t" << Ib[ii] << "\t" << rep_loc[ii] << endl;
	outfile << endl;
*/

	t_states = wall_clock() - t_states;

	outfile << endl << "Total size of all states/rep loc vector " << rep_loc.size() << endl;
	outfile << " SECTOR : " <<  Sz << "\t(" << px << "," << py << "," <<pz <<")\t" << n_p << endl;	
	outfile << "Time to construct pStates\t " << t_states/60.0 << " minutes" << endl;
	
	if(n_p == 0)
	{
		outfile << "\nNo states in psector\n";
		exit(1);
	}
	outfile << "------------------------------------------------------------------" << endl;
	
	vector< complex<double> > pblock_eigs;
	vector< vector< complex<double> > > pblock_eigvecs;

	double lanczos_time = wall_clock();

	lanczos_evecs(	pblock_states, norm, norm_inv, rep_loc, state_ph,
			pblock_eigs, pblock_eigvecs, find_ev,
			lambda, n_spins, px, py, pz, n_ones, n_p,
			adj_list, Ch, ch_list, ChCh, chch_list, T1, T2, T3, n_sec_ev, tol, iterations, ncycles,
			Ia, Ib, bits_right, outfile);
	
	outfile << "-----------------------------------------------------------" << endl;
	outfile << "-----------------------------------------------------------" << endl;
	outfile << " Wall clock for Lanczos " << (wall_clock() - lanczos_time)/60.0 << endl;


	outfile<<"pblock_eigs.size() = "<<pblock_eigs.size()<<endl;
	outfile  << "\n\nEnergies\t\tEnergies/site\t\tSz\tkx\tky\tkz\tJz\tJ2\n\n";
	outfile.flush();
		
	for(int i = 0; i < pblock_eigs.size(); i++)	
		outfile << fixed << setprecision(15) << real(pblock_eigs[i]) << "\t" << real(pblock_eigs[i])/double(n_spins) <<"\t" << setprecision (3) << Sz << "\t" << px << "\t" << py << "\t" << pz << "\t" << lambda << "\t" << J2 << "\n";  	

	outfile.flush();
	outfile<<" find_ev = "<< find_ev << endl;
	outfile.flush();
	
	if(find_ev)
	{	
		// double ts_unwrap = wall_clock();
		
		for(int iib = 0; iib < min(n_sec_ev,int(pblock_eigs.size())); iib++)
		{
			int t_tmp_int = iib + 1 + (int)'0';
			char t_tmp = (char)t_tmp_int;

			string filepath_tmp;

			/*if(unwrap_peigvec)
				filepath_tmp = filepath_EV + "_EV" + t_tmp + ".txt";		
			else
				filepath_tmp = filepath_EV + "_pEV" + t_tmp + ".txt";*/
			
			if(unwrap_peigvec)
				filepath_tmp = filepath_EV + "_EV" + to_string(iib) + ".txt";		
			else
				filepath_tmp = filepath_EV + "_pEV" + to_string(iib) + ".txt";
	
	
			ofstream outfile_vecs(filepath_tmp.c_str());
			outfile_vecs << "\n Eigenvector for state \n";
			outfile_vecs << "Sz\t(kx,ky,kz)\tJz/lambda\n";
			outfile_vecs << Sz <<"\t(" << kx << "," << ky << ","<< kz <<")\t" << lambda << endl;
		
			outfile_vecs << "\n\n-------------------------------------------";
			outfile_vecs << "\nState Energy " << pblock_eigs[iib] << endl;

			if(unwrap_peigvec)
			{
				vector< complex <double> > v_coeff = pblock_eigvecs[iib];

				vector<ids> psi = Unwrap_pState(v_coeff, pblock_states, norm, px, py, pz, T1, T2, T3);	
			
				outfile_vecs << "\nEIGENVECTOR\n";
				for(int64_t jb = 0; jb < psi.size(); jb++)
				{
					outfile_vecs << psi[jb].stateID << "\t" << fixed << setprecision(15) << psi[jb].coeff << endl;
				}
			
			}
			else
			{
				outfile_vecs << "\nMomentum EIGENVECTOR\n";
			
				outfile_vecs << "\nnp \t\t px \t\tpy \t\tpz\n";
				outfile_vecs << n_p << "\t\t" << px << "\t\t" << py << "\t\t" << pz <<endl;

				outfile_vecs << "\n\npblock_states \t\t Norm \t\t\t v_coeff\n\n";
				
				for(int64_t jb = 0; jb < n_p; jb++)
					outfile_vecs << pblock_states[jb] << "\t\t\t" << norm[jb] << "\t\t\t" << fixed << setprecision(15) << pblock_eigvecs[iib][jb]<< endl;

			}
			outfile_vecs.close();
		}
		// double te_unwrap = wall_clock();
		// outfile << "\nTime to unwrap state = " << (te_unwrap - ts_unwrap)/60 << " minutes\n";			
	}
}

