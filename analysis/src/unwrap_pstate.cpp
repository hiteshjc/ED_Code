
#include "unwrap_pstate.h"

double vpres = 1e-15;

// Function to Unwrap a pState back to its original Basis states
void Unwrap_pState(	vector< complex <double> > &v_coeff,
			vector<int64_t> &pBasis,
			vector<double> &Norm,
			double &px,
			double &py,
			vector<int> &T1,
			vector<int> &T2,
			vector<ids> &psi,
			int n_sites,
			int n_ones,
			ofstream &outfile)
{

	int64_t n_states = n_choose_k(n_sites, n_ones);

	vector<int64_t> basis;
	constrained_dets(n_sites, n_ones, basis);

//	outfile << basis.size() << "\tn_states = " << n_states << endl;

	psi.resize(n_states);
	for(int64_t i = 0; i < n_states; i++)
	{
		psi[i].coeff = 0.0;
		psi[i].stateID = basis[i];
	}

	int64_t n = v_coeff.size();

//	outfile << "n_states " << n_states << "\t psi size = " << psi.size() << endl;

	int L1, L2;
	LxLy(L1,L2, T1, T2);

	#pragma omp parallel for
	for(int64_t i = 0; i < n; i++)
	{
		int64_t v_id = pBasis[i];
		complex<double> tmp_coeff = v_coeff[i]/Norm[i];
		vector<ids> tmp_ids;
			
		add2ids(tmp_ids, v_id, tmp_coeff);		
		complex< double > ph = complex< double >(1.0,0.0);
	
		for(int x = 0; x < L1; x++)
		{
			if(x > 0)
			{
				translateT(v_id, T1);	
				complex< double > argx = complex< double >(0.0, -px);
				ph = ph*exp(argx);
				
				complex<double> tmpx_coeff = v_coeff[i]*ph/Norm[i];
				add2ids(tmp_ids, v_id, tmpx_coeff); 

			}

			for(int y = 0; y < L2; y++)
			{
				if(y > 0)
				{
					translateT(v_id, T2);
					complex< double > argy = complex< double >(0.0, -py);
					ph = ph*exp(argy);
					complex<double> tmpy_coeff = v_coeff[i]*ph/Norm[i];				
					add2ids(tmp_ids, v_id, tmpy_coeff); 
				}
			} 
		}

		for(int l = 0; l < tmp_ids.size(); l++)
		{

			int64_t tmp_stateID = tmp_ids[l].stateID;
			int64_t basis_j = findstate_basis(tmp_stateID, basis);
			if(basis_j >=0)
			{	
				psi[basis_j] = tmp_ids[l];
			}
			
		} 
	}	

}




void read_n_unwrap(	std::string filein,
			std::string filepEV,
			vector<int> &T1,
			vector<int> &T2,
			vector<ids> &psi,
			ofstream &outfile)
{
	int N; double Sz;
	read_N_Sz(filein, N, Sz);

	int nones = int(2*Sz);
	nones = (nones + N)/2;

	int64_t n_p;
	double px, py;
	vector<int64_t> pBasis;
	vector<double> norm;
	vector< complex <double> > v_coeff;

	read_pEV(n_p, px, py, pBasis, norm, v_coeff, filepEV);
	
	Unwrap_pState( v_coeff, pBasis, norm, px, py, T1, T2, psi, N, nones, outfile);
}
