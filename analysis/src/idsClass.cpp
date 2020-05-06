
#include "global.h"
#include "idsClass.h"


// Function to add an element to vector<ids>

void add2w(	vector<ids> &w,
	   	int64_t j,
		complex< double > coeff)
{
	bool addj = 0;
	for(int64_t l = 0; l < w.size(); l++)
	{
		if(w[l].stateID == j)
		{
			addj = 1;
			w[l].coeff = w[l].coeff + coeff;
			break;
		}  
	}

	if(!addj)
	{
		ids tmp; 
		tmp.stateID = j; tmp.coeff = coeff; 
		w.push_back(tmp);
	}  
}

// Function that adds an element to vector<ids>
// Used to unwrap a momentum state back to its basis elements
void add2ids(	vector<ids> &w,
	   	int64_t j,
		complex< double > coeff)
{
	bool addj = 0;
	for(int64_t l = 0; l < w.size(); l++)
	{
		if(w[l].stateID == j)
		{
			addj = 1;
			break;
		}  
	}

	if(!addj)
	{
		ids tmp; 
		tmp.stateID = j; tmp.coeff = coeff; 
		w.push_back(tmp);
	}  
}


// Function to add two vector<ids>
vector<ids> add_ids(	vector<ids> &psi1,
			vector<ids> &psi2,
			complex<double> &c1,
			complex<double> &c2)
{
	vector<ids> psi;

	for(int i = 0, i1 = 0, i2 = 0; i1 < psi1.size() || i2 < psi2.size(); i++)
	{
		if(psi1[i1].stateID == psi2[i2].stateID)
		{
			ids tmp;
			tmp.stateID = psi1[i1].stateID;
			tmp.coeff = c1*psi1[i1].coeff + c2*psi2[i2].coeff; 	
			psi.push_back(tmp);
			i1++; i2++;
		}
		else if(psi1[i1].stateID < psi2[i2].stateID)
		{
			ids tmp;
			tmp.stateID = psi1[i1].stateID;
			tmp.coeff = c1*psi1[i1].coeff; 	
			psi.push_back(tmp);
			i1++; 

		}
		else
		{
			ids tmp;
			tmp.stateID = psi2[i2].stateID;
			tmp.coeff = c2*psi2[i2].coeff; 	
			psi.push_back(tmp);
			i2++;
		}
	} 
	
	return psi;

}

// ////////////////////////////// //
// Dot product of two vector<ids> //
// ////////////////////////////// //

complex<double> dot_product_ids(	vector<ids> &psi1,
					vector<ids> &psi2)
{
	complex<double> sum = complex<double>(0.0,0.0);

	for(int64_t i = 0; i < psi1.size(); i++)
		sum += conj(psi1[i].coeff)*psi2[i].coeff;

	return sum;
}

// ///////////// //
// Binary search //
// ///////////// //


int64_t find_ids_state( vector<ids> &psi,
			int64_t &id)
{
	int64_t j = -1;

	int64_t b_min = 0, b_max = psi.size() - 1;
	do{
		int64_t b = b_min + (b_max - b_min)/2;
		if(id < psi[b].stateID )
			b_max = b - 1;
		else if (id > psi[b].stateID )
			b_min = b + 1;
		else
		{	j = b;	break;}
	}while(b_max >= b_min);

	return j;

}











