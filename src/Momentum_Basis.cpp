#include "Momentum_Basis.h"

// Function to translates Basis state
void translateT(	int64_t &BasisID,
			const vector<int> &T)
{
	int64_t tmp_id = 0;

	for(int i = 0; i < T.size(); i++)
	{
		bool t = btest64(BasisID, i);
		if (t == 1)
			tmp_id = ibset64(tmp_id, T[i]);
	} 
	BasisID = tmp_id;
}

// Function that determines the maximum translations along two directions L1 & L2
void LxLy(	int &L1, int &L2,
	  	vector<int> &T1,
		vector<int> &T2)
{
	L1 = 0, L2 = 0;
	int64_t tmp_id = 1;
	do{
		translateT(tmp_id, T1);
		L1++;
	}while(tmp_id != 1);

	tmp_id = 1;

	do{
		translateT(tmp_id, T2); 
		L2++;
	}while(tmp_id != 1);  
}

void LxLyLz(	int &L1, int &L2, int &L3,
	  	vector<int> &T1,
		vector<int> &T2,
		vector<int> &T3)
{
	L1 = 0, L2 = 0, L3 = 0;
	int64_t tmp_id = 1;
	do{
		translateT(tmp_id, T1);
		L1++;
	}while(tmp_id != 1);

	tmp_id = 1;

	do{
		translateT(tmp_id, T2); 
		L2++;
	}while(tmp_id != 1);  
	
	tmp_id = 1;

	do{
		translateT(tmp_id, T3); 
		L3++;
	}while(tmp_id != 1);  
}


// Funtion that checks if a basis state is the representative
bool ifnotrepresentative(	int64_t &v_id,
				int &L1, int &L2,
			 	vector<int> &T1,
				vector<int> &T2)
{
	bool b = 0;
	int64_t tmp_id = v_id, rep_id = v_id;

	for(int x = 0; x < L1; x++)
	{
		if(x > 0)
		{
			translateT(tmp_id, T1); 
			if(rep_id > tmp_id)
			{
				rep_id = tmp_id;
				b = 1;
				x = L1;// y = L2;
				break;
			}
		}

		for(int y = 0; y < L2 && x < L1; y++)
		{
			if(y > 0)
			{
				translateT(tmp_id, T2);
				if(rep_id > tmp_id)
				{
					rep_id = tmp_id;
					b = 1;
					x = L1; y = L2;
					break;
				}
			}
		}
	}

	return b;
}

// Funtion that checks if particular Basis state was already visited
bool visited(	vector<int64_t> &visit,
		int64_t &tmp_id)
{
	bool b = 0;

	for(int64_t i = 0; i < visit.size(); i++)
	{
		if(visit[i] == tmp_id)
		{
			b = 1; i = visit.size();     
			break;      
		}    
	} 
	return b;
}

// Function that constructs the momentum basis states
// Function outputs the state's representative id and its normalization

void Construct_States(	vector<int64_t> &v_sz,
			vector<double> &norm,
		      	vector<int> &T1,
			vector<int> &T2,
		      	int &num_ones,
			double &px,
			double &py,
			int64_t &n_states,
			ofstream &outfile)
{
	v_sz.clear(); norm.clear();
	int n_sites = T1.size();
	n_states = 0;
	int L1, L2;
	// Routine to determine the total length/dimensions of the system in both directions
	LxLy(L1, L2, T1, T2);
	////////////////////////////////////////////////////

	// Routine for generating basis states 
	int64_t temp,temp1,temp2,temp3,temp4,temp5,temp6;
	int64_t i,num_configs;

	if (n_sites<num_ones)
	{
		outfile<<"ERROR: Num sites < Num_ones"<<endl;
		exit(1);  
	}

	num_configs = n_choose_k(n_sites,num_ones);
	//  outfile << "\nNUM of CONfig " << num_configs;
	int64_t tmp_id; 

	for( i = 0; i < num_configs; i++)
	{ 
		if(i == 0)
		tmp_id = pow(2,num_ones) - 1;
		else
		{
			temp  = (tmp_id | tmp_id - 1) + 1;
			temp2 = (temp) & (-temp);
			temp3=  (tmp_id) & (-tmp_id);
			temp4=temp2/temp3;
			temp5=temp4>>1;
			temp6=temp5-1;
			tmp_id = temp|temp6;
		}  

		//    //////// Check to see if state was already explored)
		bool rep = ifnotrepresentative(tmp_id, L1, L2, T1, T2);
		if(rep) 
		continue;

		complex< double > ph = complex< double >(1.0,0.0);
		complex< double > id_coeff = complex<double> (1.0,0.0);
		int64_t v_id = tmp_id;
		int lid = 1; 
		for(int x = 0; x < L1; x++)
		{
			if(x > 0)
			{
				translateT(v_id, T1);	
				complex< double > argx = complex< double >(0.0, -px);
				ph = ph*exp(argx);
				if(v_id == tmp_id)
				{  id_coeff = id_coeff + ph;     lid++;}
			}

			for(int y = 0; y < L2; y++)
			{
				if(y > 0)
				{
					translateT(v_id, T2);
					complex< double > argy = complex< double >(0.0, -py);
					ph = ph*exp(argy);
					if(v_id == tmp_id)
					{   id_coeff = id_coeff + ph;	lid++;}
				}
			} 
		}

		if(abs(real(id_coeff)) > 1e-8 || abs(imag(id_coeff)) > 1e-8)
		{
			v_sz.push_back(tmp_id);
			int64_t N = L1*L2/lid; 
			norm.push_back(sqrt(N));
			n_states++;
		}   
	}

}

// Locates the representative pBasis state using binary search
int64_t findstate(	int64_t &rep_id,
	       		vector<int64_t> &pBasis,
			int64_t st,
			int64_t end)
{
	int64_t j = -1;

	int64_t b_min = st, b_max = end;
	do{
		int64_t b = b_min + (b_max - b_min)/2;
		if(rep_id < pBasis[b] )
			b_max = b - 1;
		else if (rep_id > pBasis[b] )
			b_min = b + 1;
		else
		{	j = b;	break;}
	}while(b_max >= b_min);

	return j;

}

// Funtion that checks if a basis state is the representative and returns rep

int64_t representative2(	int &lx,
				int &ly,
				int &L1,
				int &L2,
				int64_t &v_id,
				vector<int> &T1,
				vector<int> &T2)
{
	lx = 0, ly = 0;
	int64_t rep_id = v_id, tmp_id = v_id;
	int tx = 0, ty = 0;			  

	for(int x = 0; x < L1; x++)
	{
		if(x > 0)
		{
			translateT(tmp_id, T1); 
			tx++;

			if(rep_id > tmp_id)
			{
				rep_id = tmp_id;
				lx = L1 - tx; ly = L2 - ty;
				break;
			}
		}
		for(int y = 0; y < L2; y++)
		{
			if(y>0)
			{
				translateT(tmp_id, T2);
				ty++;

				if(rep_id > tmp_id)
				{
					rep_id = tmp_id;
					lx = L1 - tx; ly = L2 - ty;
					break;
				}
			}
		}
	}
	
	return rep_id;
}




// Function to construct p-states and store representative ids

void Construct_States(	vector<int64_t> &v_sz,
			vector<double> &norm,
		      	vector<double> &norm_inv,
			vector<int> &T1,
			vector<int> &T2,
		      	int &num_ones,
			double &px,
			double &py,
			int64_t &n_states,
			vector<int64_t> &rep_loc, 
		        vector< complex<double> > &state_ph,
			vector<int64_t> &Ia,
			vector<int64_t> &Ib,
			int &bits_right,
			ofstream &outfile)
{

	v_sz.clear(); norm.clear();
	int n_sites = T1.size();
	int half_sites=n_sites/2;
	if(n_sites%2 == 1) ++half_sites;

	/*int half_sitesa;
        if (n_sites%2==0)
        {
         half_sites= n_sites/2;
	 half_sitesa=half_sites;
        }
	else
	{
	half_sites=(n_sites+1)/2;
	half_sitesa=(n_sites-1)/2;
	}*/
	n_states = 0;
	int L1, L2;
	LxLy(L1, L2, T1, T2);
	////////////////////////////////////////////////////

	// Routine for generating basis states 
	int64_t temp,temp1,temp2,temp3,temp4,temp5,temp6;
	int64_t i,num_configs;

	if (n_sites<num_ones)
	{
		outfile<<"ERROR: Num sites < Num_ones"<<endl;
		exit(1);  
	}

	num_configs = n_choose_k(n_sites,num_ones);

	int64_t tmp_id;
	vector<int64_t> rep_ids, all_states;
	rep_ids.resize(num_configs);
	rep_loc.resize(num_configs);	all_states.resize(num_configs);

	vector<int> state_lx, state_ly;
	state_lx.resize(num_configs);	state_ly.resize(num_configs);
	state_ph.resize(num_configs);

// Set up Lin Tables
// Function sets up default indices for Ia
// Labels for Ib will be detrmined below
	Lin_Tables(Ia, Ib, half_sites, outfile);
	bits_right = pow(2, half_sites) - 1;

	for( i = 0; i < num_configs; i++)
	{ 
		if(i == 0)
		tmp_id = pow(2,num_ones) - 1;
		else
		{
			temp  = (tmp_id | tmp_id - 1) + 1;
			temp2 = (temp) & (-temp);
			temp3=  (tmp_id) & (-tmp_id);
			temp4=temp2/temp3;
			temp5=temp4>>1;
			temp6=temp5-1;
			tmp_id = temp|temp6;
		}  

		all_states[i] = tmp_id;

///////////////////////////////////////////////////
// Complete Lin Table generation /////
		int left, right;
		split_id(tmp_id, left, right, bits_right, n_sites);

		if(Ib[left] == -1)
			Ib[left] = i - Ia[right] + 1;
		else if (i != (Ib[left] + Ia[right] - 1) )
			outfile << "\n\n ERROR in LIN TABLE GENERATION " << endl << endl;
///////////////////////////////////////////////////
//
		int lx, ly;
		int64_t tmp_rep_id = representative2(lx, ly, L1, L2, tmp_id, T1, T2);
	
		if(tmp_rep_id != tmp_id)
		{
			int64_t j = findstate(tmp_rep_id, all_states, 0, i);
			rep_ids[i] = rep_ids[j];
			if(j != i)
			lx = lx + state_lx[j] - L1; ly = ly + state_ly[j] - L2;

			rep_loc[i] = findstate(rep_ids[i], v_sz, 0, v_sz.size() - 1);
		}	
		else
		{
			rep_ids[i] = tmp_id;
			rep_loc[i] = v_sz.size();
		}
	
		lx = lx % L1;
		ly = ly % L2;

		if(lx < 0)
			lx = lx + L1;
		if(ly < 0)
			ly = ly + L2;

		state_lx[i] = lx; state_ly[i] = ly;
		state_ph[i] = exp( complex<double>(0.0,-px*lx - py*ly) );

		if(rep_ids[i] == tmp_id)
		{
			complex<double> ph = complex<double>(1.0, 0.0);
			complex< double > id_coeff = complex<double> (1.0,0.0);
			int lid = 1;	

			int64_t v_id = tmp_id;
			for(int x = 0; x < L1; x++)
			{
				if(x > 0)
				{
					translateT(v_id, T1);	
					complex< double > argx = complex< double >(0.0, -px);
					ph = ph*exp(argx);
					if(v_id == tmp_id)
					{	
						id_coeff = id_coeff + ph; 
						lid++;
					}
				}

				for(int y = 0; y < L2; y++)
				{
					if(y > 0)
					{
						translateT(v_id, T2);
						complex< double > argy = complex< double >(0.0, -py);
						ph = ph*exp(argy);
						if(v_id == tmp_id)
						{
							id_coeff = id_coeff + ph;
							lid++;
						}
					}
				} 
			}

			if(abs(real(id_coeff)) > 1e-8 || abs(imag(id_coeff)) > 1e-8)
			{
				v_sz.push_back(tmp_id);
				int64_t N = L1*L2/lid; 
				norm.push_back(sqrt(N));
				norm_inv.push_back(1/sqrt(N));
				n_states++;
			}
			else
				rep_loc[i] = -1;
		}		
	}
}



// Function to construct p-states and store representative ids
/////////////////////////////////////////////////////////////////////
int64_t representative_p(	int &lx,
				int &ly,
				int &L1,
				int &L2,
				int64_t &v_id,
				vector<int> &T1,
				vector<int> &T2,
				double &px,
				double &py,
				bool &p_id,
				int &lid)
{
	lx = 0, ly = 0;
	int64_t rep_id = v_id, tmp_id = v_id;
	int tx = 0, ty = 0;			  

	complex<double> ph = complex<double>(1.0, 0.0);
	complex< double > id_coeff = complex<double> (1.0,0.0);
	lid = 1;	


	for(int x = 0; x < L1; x++)
	{
		if(x > 0)
		{
			translateT(tmp_id, T1); 
			tx++;
			complex< double > argx = complex< double >(0.0, -px);
			ph = ph*exp(argx);

			if(v_id == tmp_id)
			{	
				id_coeff = id_coeff + ph; 
				lid++;
			}
					
			if(rep_id > tmp_id)
			{
				rep_id = tmp_id;
				lx = L1 - tx; ly = L2 - ty;
			}
		}
		for(int y = 0; y < L2; y++)
		{
			if(y>0)
			{
				translateT(tmp_id, T2);
				ty++;
				complex< double > argy = complex< double >(0.0, -py);
				ph = ph*exp(argy);

				if(v_id == tmp_id)
				{
					id_coeff = id_coeff + ph;
					lid++;
				}

				if(rep_id > tmp_id)
				{
					rep_id = tmp_id;
					lx = L1 - tx; ly = L2 - ty;
				}
			}
		}
	}
	
	if(abs(real(id_coeff)) > 1e-8 || abs(imag(id_coeff)) > 1e-8)
		p_id = 1;
	else 	p_id = 0;

	return rep_id;
}
//////////////// 3 dims //////////////////////////////////////////////////////////
int64_t representative_p(	int &lx,
				int &ly,
				int &lz,
				int &L1,
				int &L2,
				int &L3,
				int64_t &v_id,
				vector<int> &T1,
				vector<int> &T2,
				vector<int> &T3,
				double &px,
				double &py,
				double &pz,
				bool &p_id,
				int &lid)
{
	lx = 0, ly = 0; lz = 0;
	int64_t rep_id = v_id, tmp_id = v_id;
	int tx = 0, ty = 0, tz = 0;			  

	complex<double> ph = complex<double>(1.0, 0.0);
	complex< double > id_coeff = complex<double> (1.0,0.0);
	lid = 1;	


	for(int x = 0; x < L1; x++)
	{
		if(x > 0)
		{
			translateT(tmp_id, T1); 
			tx++;
			complex< double > argx = complex< double >(0.0, -px);
			ph = ph*exp(argx);

			if(v_id == tmp_id)
			{	
				id_coeff = id_coeff + ph; 
				lid++;
			}
					
			if(rep_id > tmp_id)
			{
				rep_id = tmp_id;
				lx = L1 - tx; ly = L2 - ty; lz = L3 - tz;
			}
		}
		for(int y = 0; y < L2; y++)
		{
			if(y>0)
			{
				translateT(tmp_id, T2);
				ty++;
				complex< double > argy = complex< double >(0.0, -py);
				ph = ph*exp(argy);

				if(v_id == tmp_id)
				{
					id_coeff = id_coeff + ph;
					lid++;
				}

				if(rep_id > tmp_id)
				{
					rep_id = tmp_id;
					lx = L1 - tx; ly = L2 - ty; lz = L3 - tz;
				}
			}
			for(int z = 0; z < L3; z++)
			{
				if(z>0)
				{
					translateT(tmp_id, T3);
					tz++;
					complex< double > argz = complex< double >(0.0, -pz);
					ph = ph*exp(argz);

					if(v_id == tmp_id)
					{
						id_coeff = id_coeff + ph;
						lid++;
					}

					if(rep_id > tmp_id)
					{
						rep_id = tmp_id;
						lx = L1 - tx; ly = L2 - ty; lz = L3 - tz;
					}
				}
			}
		}
	}
	
	if(abs(real(id_coeff)) > 1e-8 || abs(imag(id_coeff)) > 1e-8)
		p_id = 1;
	else 	p_id = 0;

	return rep_id;
}

/////////////////////////////////////////////////////////////////////////////////
void pConstruct_States(	vector<int64_t> &v_sz,
			vector<double> &norm,
		      	vector<double> &norm_inv,
			vector<int> &T1,
			vector<int> &T2,
		      	int &num_ones,
			double &px,
			double &py,
			int64_t &n_states,
			vector<int64_t> &rep_loc, 
		        vector< complex<double> > &state_ph,
			vector<int64_t> &Ia,
			vector<int64_t> &Ib,
			int &bits_right,
			ofstream &outfile)
{

	v_sz.clear(); norm.clear();
	int n_sites = T1.size();
	//int half_sites;
	/*int half_sitesa;
        if (n_sites%2==0)
        {
         half_sites= n_sites/2;
	 half_sitesa=half_sites;
        }
	else
	{
	half_sites=(n_sites+1)/2;
	half_sitesa=(n_sites-1)/2;
	}*/
	int half_sites = n_sites/2;
	if(n_sites%2 == 1) ++half_sites;
        n_states = 0;
	int L1, L2;
	LxLy(L1, L2, T1, T2);
	////////////////////////////////////////////////////

	// Routine for generating basis states 
	int64_t temp,temp1,temp2,temp3,temp4,temp5,temp6;
	int64_t i,num_configs;

	if (n_sites<num_ones)
	{
		outfile<<"ERROR: Num sites < Num_ones"<<endl;
		exit(1);  
	}

	num_configs = n_choose_k(n_sites,num_ones);

	vector<int64_t> rep_ids(num_configs), all_states(num_configs);
	vector<bool> p_ids(num_configs); vector<int> l_ids(num_configs);
	rep_loc.resize(num_configs);


//	vector<int> state_lx(num_configs), state_ly(num_configs);
//	state_lx.resize(num_configs);	state_ly.resize(num_configs);
	state_ph.resize(num_configs);

// Set up Lin Tables
// Function sets up default indices for Ia
// Labels for Ib will be detrmined below in the for loop
        outfile.flush();
        outfile<<"Before Lin"<<endl;	
        outfile.flush();
	Lin_Tables(Ia, Ib, half_sites, outfile);
        outfile.flush();
        outfile<<"After Lin"<<endl;	
        outfile.flush();
// Mask for extracting bits from one half of int64_t integer
	bits_right = pow(2, half_sites) - 1;


        outfile<<"Before constrained dets"<<endl;	
        outfile.flush();
	constrained_dets(n_sites, num_ones, all_states);
        outfile<<"After constrained dets"<<endl;	
        outfile.flush();

	#pragma omp parallel for
	for( i = 0; i < num_configs; i++)
	{ 
		int64_t tmp_id = all_states[i];

///////////////////////////////////////////////////
// Complete Lin Table generation /////
		int left, right;
		split_id(tmp_id, left, right, bits_right, n_sites);

		if(Ib[left] == -1)
			Ib[left] = i - Ia[right] + 1;
		else if (i != (Ib[left] + Ia[right] - 1) )
			outfile << "\n\n ERROR in LIN TABLE GENERATION " << endl << endl;
///////////////////////////////////////////////////
//
		int lx, ly;	bool tmp_p_id;	int tmp_l_id;
		rep_ids[i] = representative_p(lx, ly, L1, L2, tmp_id, T1, T2, px, py, tmp_p_id, tmp_l_id);
		p_ids[i] = tmp_p_id;
		l_ids[i] = tmp_l_id;
	
		lx = lx % L1;
		ly = ly % L2;

		if(lx < 0)
			lx = lx + L1;
		if(ly < 0)
			ly = ly + L2;

		state_ph[i] = exp( complex<double>(0.0,-px*lx - py*ly) );
	}

	for(int64_t i = 0; i < num_configs; i++)
	{
		int64_t tmp_id = all_states[i];
		if(rep_ids[i] == tmp_id)
		{
			if(p_ids[i])//abs(real(id_coeff)) > 1e-8 || abs(imag(id_coeff)) > 1e-8)
			{
				v_sz.push_back(tmp_id);
				int64_t N = L1*L2/l_ids[i]; 
				norm.push_back(sqrt(N));
				norm_inv.push_back(1/sqrt(N));
				n_states++;
			}
		}
	}

	#pragma omp parallel for
	for(int64_t i = 0; i < num_configs; i++)
	{
		int64_t tmp_id = rep_ids[i];
		rep_loc[i] = findstate(tmp_id, v_sz, 0, v_sz.size() - 1);
	}
}

/////////////////////////////////////////////////////////////////////////////////
void pConstruct_States(	vector<int64_t> &v_sz,
			vector<double> &norm,
		      	vector<double> &norm_inv,
			vector<int> &T1,
			vector<int> &T2,
			vector<int> &T3,
		      	int &num_ones,
			double &px,
			double &py,
			double &pz,
			int64_t &n_states,
			vector<int64_t> &rep_loc, 
		        vector< complex<double> > &state_ph,
			vector<int64_t> &Ia,
			vector<int64_t> &Ib,
			int &bits_right,
			ofstream &outfile)
{

	v_sz.clear(); norm.clear();
	int n_sites = T1.size();
	int half_sites = n_sites/2;
	if(n_sites%2 == 1) ++half_sites;
        n_states = 0;
	int L1, L2, L3;
	LxLyLz(L1, L2, L3, T1, T2, T3);
	////////////////////////////////////////////////////

	// Routine for generating basis states 
	int64_t temp,temp1,temp2,temp3,temp4,temp5,temp6;
	int64_t i,num_configs;

	if (n_sites<num_ones)
	{
		outfile<<"ERROR: Num sites < Num_ones"<<endl;
		exit(1);  
	}

	num_configs = n_choose_k(n_sites,num_ones);

	vector<int64_t> rep_ids(num_configs), all_states(num_configs);
	vector<bool> p_ids(num_configs); vector<int> l_ids(num_configs);
	rep_loc.resize(num_configs);


//	vector<int> state_lx(num_configs), state_ly(num_configs);
//	state_lx.resize(num_configs);	state_ly.resize(num_configs);
	state_ph.resize(num_configs);

// Set up Lin Tables
// Function sets up default indices for Ia
// Labels for Ib will be detrmined below in the for loop
	Lin_Tables(Ia, Ib, half_sites, outfile);
// Mask for extracting bits from one half of int64_t integer
	bits_right = pow(2, half_sites) - 1;


	constrained_dets(n_sites, num_ones, all_states);

	#pragma omp parallel for
	for( i = 0; i < num_configs; i++)
	{ 
		int64_t tmp_id = all_states[i];

///////////////////////////////////////////////////
// Complete Lin Table generation /////
		int left, right;
		split_id(tmp_id, left, right, bits_right, n_sites);

		if(Ib[left] == -1)
			Ib[left] = i - Ia[right] + 1;
		else if (i != (Ib[left] + Ia[right] - 1) )
			outfile << "\n\n ERROR in LIN TABLE GENERATION " << endl << endl;
///////////////////////////////////////////////////
//
		int lx, ly, lz;	bool tmp_p_id;	int tmp_l_id;
		rep_ids[i] = representative_p(lx, ly, lz, L1, L2, L3, tmp_id, T1, T2, T3, px, py, pz, tmp_p_id, tmp_l_id);
		p_ids[i] = tmp_p_id;
		l_ids[i] = tmp_l_id;
	
		lx = lx % L1;
		ly = ly % L2;
		lz = lz % L3;

		if(lx < 0)
			lx = lx + L1;
		if(ly < 0)
			ly = ly + L2;
		if(lz < 0)
			lz = lz + L3;

		state_ph[i] = exp( complex<double>(0.0,-px*lx - py*ly - pz*lz) );
	}

	for(int64_t i = 0; i < num_configs; i++)
	{
		int64_t tmp_id = all_states[i];
		if(rep_ids[i] == tmp_id)
		{
			if(p_ids[i])//abs(real(id_coeff)) > 1e-8 || abs(imag(id_coeff)) > 1e-8)
			{
				v_sz.push_back(tmp_id);
				int64_t N = L1*L2*L3/l_ids[i]; 
				norm.push_back(sqrt(N));
				norm_inv.push_back(1/sqrt(N));
				n_states++;
			}
		}
	}

	#pragma omp parallel for
	for(int64_t i = 0; i < num_configs; i++)
	{
		int64_t tmp_id = rep_ids[i];
		rep_loc[i] = findstate(tmp_id, v_sz, 0, v_sz.size() - 1);
	}
}


// Function to Unwrap a pState back to its original Basis states
vector< ids > Unwrap_pState(	vector< complex <double> > &v_coeff,
				vector<int64_t> &pBasis,
				vector<double> &Norm,
				double &px,
				double &py,
				double &pz,
				vector<int> &T1,
				vector<int> &T2,
				vector<int> &T3)
{
	vector<ids> Basis;
	int64_t n = v_coeff.size();

	int L1, L2, L3;
	LxLyLz(L1,L2, L3, T1, T2, T3);

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

				for(int z = 0; z < L3; z++)
				{
					if(z > 0)
					{
						translateT(v_id, T3);
						complex< double > argz = complex< double >(0.0, -pz);
						ph = ph*exp(argz);
						complex<double> tmpy_coeff = v_coeff[i]*ph/Norm[i];				
						add2ids(tmp_ids, v_id, tmpy_coeff); 
					}
				} 

			} 
		}

		for(int l = 0; l < tmp_ids.size(); l++)
		{
			if( abs( real(tmp_ids[l].coeff) ) > 1e-8 || abs( imag(tmp_ids[l].coeff) ) > 1e-8 )
			{	Basis.push_back(tmp_ids[l]); }
		} 
	}	

	merge_sort_ids(Basis, 0, Basis.size() - 1);
	
	return Basis;
}

// Function to Unwrap a pState back to its original Basis states
vector< ids > Unwrap_pState(	vector< complex <double> > &v_coeff,
				vector<int64_t> &pBasis,
				vector<double> &Norm,
				double &px,
				double &py,
				vector<int> &T1,
				vector<int> &T2)
{
	vector<ids> Basis;
	int64_t n = v_coeff.size();

	int L1, L2;
	LxLy(L1,L2, T1, T2);

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
			if( abs( real(tmp_ids[l].coeff) ) > 1e-8 || abs( imag(tmp_ids[l].coeff) ) > 1e-8 )
			{	Basis.push_back(tmp_ids[l]); }
		} 
	}	

	merge_sort_ids(Basis, 0, Basis.size() - 1);
	
	return Basis;
}
