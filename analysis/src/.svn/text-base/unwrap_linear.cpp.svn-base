
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
			vector<ids> &Basis)
{
	int64_t n = v_coeff.size();

	int L1, L2;
	LxLy(L1,L2, T1, T2);

//	#pragma omp parallel for
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
			if( abs( real(tmp_ids[l].coeff) ) > vpres || abs( imag(tmp_ids[l].coeff) ) > vpres )
			{	Basis.push_back(tmp_ids[l]); }
		} 
	}	

	merge_sort_ids(Basis, 0, Basis.size() - 1);

}

int64_t findstate_basis(	int64_t &rep_id,
				vector<int64_t> &basis)
{
	int64_t j = -1;
	int64_t n = basis.size();

	int64_t b_min = 0, b_max = n - 1;
	do{
		int64_t b = b_min + (b_max - b_min)/2;
		if(rep_id < basis[b] )
			b_max = b - 1;
		else if (rep_id > basis[b])
			b_min = b + 1;
		else
		{
			j = b;
			break;
		}
	}while(b_max >= b_min);

	return j;

}



void read_n_unwrap(	vector<ids> &psi,
			vector<int> &T1,
			vector<int> &T2,
			string &filepEV)
{
	
	int L1, L2;
	LxLy(L1,L2, T1, T2);

	ifstream infile(filepEV.c_str());
	string line;
	string p1p2("np 		 px 		py");	

	int64_t n_p;
	double px, py; 
	while(getline(infile, line))
	{
		if(line.compare(p1p2) == 0)
		{
			infile >> n_p >> px >> py;
			break;
		}
	}	
	infile.close();
	
	string phead("pblock_states 		 Norm 			 v_coeff");
	ifstream infile1(filepEV.c_str());

	int64_t npsi = 0;

	while(getline(infile1, line))
	{
		complex<double> v_coeff;
		double norm;
		int64_t pstateID;

		if(line.compare(phead) == 0)
		{
			while(infile1 >> pstateID >> norm >> v_coeff)
			{

				int64_t v_id = pstateID;
				complex<double> tmp_coeff = v_coeff/norm;
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
						
						complex<double> tmpx_coeff = tmp_coeff*ph;
						add2ids(tmp_ids, v_id, tmpx_coeff); 

					}

					for(int y = 0; y < L2; y++)
					{
						if(y > 0)
						{
							translateT(v_id, T2);
							complex< double > argy = complex< double >(0.0, -py);
							ph = ph*exp(argy);
							complex<double> tmpy_coeff = tmp_coeff*ph;				
							add2ids(tmp_ids, v_id, tmpy_coeff); 
						}
					} 
				}

				for(int l = 0; l < tmp_ids.size(); l++)
				{
					if( abs( real(tmp_ids[l].coeff) ) > vpres || abs( imag(tmp_ids[l].coeff) ) > vpres )
					{	npsi++;//psi.push_back(tmp_ids[l]);
					}
				}
			}
		}
	}	
	infile1.close();
	ifstream infile2(filepEV.c_str());

	int64_t ipsi=0;
	psi.resize(npsi);

	while(getline(infile2, line))
	{
		complex<double> v_coeff;
		double norm;
		int64_t pstateID;

		if(line.compare(phead) == 0)
		{
			while(infile2 >> pstateID >> norm >> v_coeff)
			{

				int64_t v_id = pstateID;
				complex<double> tmp_coeff = v_coeff/norm;
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
						
						complex<double> tmpx_coeff = tmp_coeff*ph;
						add2ids(tmp_ids, v_id, tmpx_coeff); 

					}

					for(int y = 0; y < L2; y++)
					{
						if(y > 0)
						{
							translateT(v_id, T2);
							complex< double > argy = complex< double >(0.0, -py);
							ph = ph*exp(argy);
							complex<double> tmpy_coeff = tmp_coeff*ph;				
							add2ids(tmp_ids, v_id, tmpy_coeff); 
						}
					} 
				}

				for(int l = 0; l < tmp_ids.size(); l++)
				{
					if( abs( real(tmp_ids[l].coeff) ) > vpres || abs( imag(tmp_ids[l].coeff) ) > vpres )
					{	psi[ipsi++] = tmp_ids[l]; }
				}
			}
		}
	}	
	infile2.close();
	merge_sort_ids(psi, 0, psi.size() - 1);

}





void read_n_unwrap(	vector<ids> &psi,
			vector<int> &T1,
			vector<int> &T2,
			string &filepEV,
			int n_sites,
			int n_ones)
{
	
	int L1, L2;
	LxLy(L1,L2, T1, T2);

	int64_t n_states = n_choose_k(n_sites, n_ones);

	vector<int64_t> basis;
	constrained_dets(n_sites, n_ones, basis);

	psi.resize(n_states);
	for(int64_t i = 0; i < n_states; i++)
	{
		psi[i].coeff = 0.0;
		psi[i].stateID = basis[i];
	}


	ifstream infile(filepEV.c_str());
	string line;
	string p1p2("np 		 px 		py");	

	int64_t n_p;
	double px, py; 
	while(getline(infile, line))
	{
		if(line.compare(p1p2) == 0)
		{
			infile >> n_p >> px >> py;
			break;
		}
	}	
	infile.close();
	
	string phead("pblock_states 		 Norm 			 v_coeff");
	ifstream infile2(filepEV.c_str());


	while(getline(infile2, line))
	{
		complex<double> v_coeff;
		double norm;
		int64_t pstateID;

		if(line.compare(phead) == 0)
		{
			while(infile2 >> pstateID >> norm >> v_coeff)
			{

				int64_t v_id = pstateID;
				complex<double> tmp_coeff = v_coeff/norm;
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
						
						complex<double> tmpx_coeff = tmp_coeff*ph;
						add2ids(tmp_ids, v_id, tmpx_coeff); 

					}

					for(int y = 0; y < L2; y++)
					{
						if(y > 0)
						{
							translateT(v_id, T2);
							complex< double > argy = complex< double >(0.0, -py);
							ph = ph*exp(argy);
							complex<double> tmpy_coeff = tmp_coeff*ph;				
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
	}	
	infile2.close();

}




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
			int n_ones)
{

	int64_t n_states = n_choose_k(n_sites, n_ones);

	vector<int64_t> basis;
	constrained_dets(n_sites, n_ones, basis);

	psi.resize(n_states);
	for(int64_t i = 0; i < n_states; i++)
	{
		psi[i].coeff = 0.0;
		psi[i].stateID = basis[i];
	}

	int64_t n = v_coeff.size();

	int L1, L2;
	LxLy(L1,L2, T1, T2);

//	#pragma omp parallel for
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

