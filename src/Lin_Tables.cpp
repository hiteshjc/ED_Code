#include "Lin_Tables.h"

// Code given by Hitesh
//
void gen_Sz_states(	vector<int64_t> &v_sz,
			int &n_ones,
			int &n_sites,
			ofstream &outfile)
{
	// Routine for generating basis states 
	int64_t temp,temp1,temp2,temp3,temp4,temp5,temp6;

	if (n_sites < n_ones)
	{
		outfile<<"ERROR: Num sites < Num_ones"<<endl;
		exit(1);  
	}

	int64_t n_configs = n_choose_k(n_sites,n_ones);
	int64_t tmp_id; 
	v_sz.resize(n_configs);

	for(int64_t i = 0; i < n_configs; i++)
	{ 
		if(i == 0)
			tmp_id = pow(2,n_ones) - 1;
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
	
		v_sz[i] = tmp_id;
	}
}

void split_id( 	const int64_t &id,
		int &left,
		int &right,
		const int &bits_right,
		const int &n_sites)
{
	int half_sites;
	half_sites=n_sites/2;
	//if (n_sites%2==0) half_sites=n_sites/2;
	//else half_sites=(n_sites-1)/2;
	right = (id & bits_right);
	left = id >> half_sites;
}


void Lin_Tables(	vector<int64_t> &Ia,
			vector<int64_t> &Ib,
			int half_sites,
			ofstream &outfile)
{
	int64_t n_config = pow(2, half_sites);
	Ia.resize(n_config);
	Ib.resize(n_config);

	for(int i = 0; i <= half_sites; i++)
	{
		vector<int64_t> v_sz;
		gen_Sz_states(v_sz, i, half_sites, outfile);

		int64_t ia_tmp = 1;
		for(int64_t j = 0; j < v_sz.size(); j++)
		{
			Ia[v_sz[j]] = ia_tmp++;
			Ib[v_sz[j]] = -1;
		}
	}
}


