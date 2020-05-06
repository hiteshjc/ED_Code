#include "Momentum_Basis.h"

// Function to translates Basis state
void translateT(	int64_t &BasisID,
			vector<int> &T)
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


