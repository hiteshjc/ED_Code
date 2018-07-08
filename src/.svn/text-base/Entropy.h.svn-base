#ifndef GUARD_ENTROPY_H
#define GUARD_ENTROPY_H

#include "global.h"
#include "matrix.h"
#include "Basis_States.h"
#include "Matrix_Functions.h"



void set_bits_zero_Cut(	int64_t &a,
			vector<int> &Cut)
{
	for(int i = 0; i < Cut.size(); i++)
	{	a = ibclr64(a, Cut[i]);}
}

int64_t set_bits_Cut(	int64_t big_string, 
			int64_t &small_string,
			vector<int> &Cut)
{
	for(int i = 0; i < Cut.size(); i++)
	{
		bool loc = btest64(small_string, i);
		if(loc == 1)
			big_string = ibset64(big_string, Cut[i]);
		else
			big_string = ibclr64(big_string, Cut[i]);
	}

	return big_string;
}



int64_t findstate_psi(	int64_t &rep_id,
			vector<ids> &psi)
{
	int64_t j = -1;
	int64_t n = psi.size();

	int64_t b_min = 0, b_max = n - 1;
	do{
		int64_t b = b_min + (b_max - b_min)/2;
		if(rep_id < psi[b].stateID )
			b_max = b - 1;
		else if (rep_id > psi[b].stateID )
			b_min = b + 1;
		else
		{
			j = b;
			break;
		}
	}while(b_max >= b_min);

	return j;

}

void density_matrix(	zMatrix &dm,
			vector< ids > &psi,
			vector<int> &Cut)
{
	int64_t n_states = psi.size();
	int64_t n_dm_states = int64_t_power(2, Cut.size());

	for(int64_t i = 0; i < n_states; i++)
	{
		int64_t big_string = psi[i].stateID;	
		int64_t statej, statek;

		int64_t big_stringj = big_string; //set_bits_Cut(big_string, j, Cut);	
		int64_t j = 0;
		for(int ii = 0; ii < Cut.size(); ii++)
		{
			bool bii = btest64(big_stringj, Cut[ii]);
			if (bii == 1)
				j = ibset64(j, ii);
		}			

		statej = findstate_psi(big_stringj, psi);
			
		for(int64_t k = j; k < n_dm_states; k++)
		{
			int64_t big_stringk = set_bits_Cut(big_string, k, Cut);
			statek = findstate_psi(big_stringk, psi);
				
			if(statej >=0 && statek >=0 ) 
			{
				dm(j, k) += conj(psi[statej].coeff)*psi[statek].coeff;
				if(j != k)
					dm(k,j) += conj(psi[statek].coeff)*psi[statej].coeff;
			}
		}
	}

}

void density_matrixij(	zMatrix &dm,
			vector< ids > &psi1,
			vector< ids > &psi2,
			vector<int> Cut)
{
	int64_t n1 = psi1.size(), n2 = psi2.size();
	int64_t n_dm = int64_t_power(2, Cut.size());

	#pragma omp parallel for
	for(int64_t i = 0; i < n1; i++)
	{
		int64_t big_string = psi1[i].stateID;	
		int64_t statej, statek;

		int64_t big_stringj = big_string; //set_bits_Cut(big_string, j, Cut);	
		int64_t j = 0;
		for(int l = 0; l < Cut.size(); l++)
		{
			bool tmp_b = btest64(big_stringj, Cut[l]);
			if (tmp_b == 1)
				j = ibset64(j, l);
		}			

//		statej = findstate_psi(big_stringj, psi1);
		statej = i;
			
		for(int64_t k = 0; k < n_dm; k++)
		{
			int64_t big_stringk = set_bits_Cut(big_string, k, Cut);
			statek = findstate_psi(big_stringk, psi2);
				
			if(statej >=0 && statek >=0 ) 
			{
				dm(j, k) += conj(psi1[statej].coeff)*psi2[statek].coeff;
	//			if(j != k)
	//				dm(k,j) += conj(psi2[statek].coeff)*psi1[statej].coeff;
			}
		}
	}

}

	
complex<double> Von_Neumann_Entropy(	zMatrix &dm)
{
	int64_t n_dm_states = dm.NRows();

	complex<double> sum = complex<double>(0.0,0.0);

	vector< complex<double> > lambdai(n_dm_states);
	zMatrix VR(n_dm_states,n_dm_states);
	
	zMatrix_Diagonalize(dm, lambdai, VR);

	for(int64_t i = 0; i < n_dm_states; i++)
	{		
		if(abs(real(lambdai[i])) > 1e-8 || abs(imag(lambdai[i])) > 1e-8 )
		{	sum = sum + lambdai[i]*log2(real(lambdai[i]));}

	}

	return -sum;

}

complex<double> Renyi_Entropy(	zMatrix &dm)
{
	complex<double> sum = complex<double>(0.0,0.0);

	for(int64_t i = 0; i < dm.NRows(); i++)
		for(int64_t j = 0; j < dm.NRows(); j++)
			sum = sum + dm(i, j)*dm(j,i);
	
	return -log2(real(sum));

}

void entanglement_entropies(	vector<ids> &psi,
				vector<int> &Cut,
				ofstream &outfile)
{
		double ts_dm = wall_clock();
		int64_t n_dm_states = int64_t_power(2, Cut.size());
		zMatrix	dm(n_dm_states, n_dm_states, 0.0);
		density_matrix(dm, psi, Cut);
		double te_dm = wall_clock();
		outfile << "\n\n Time for Density Matrix Construction " << (te_dm - ts_dm)/60 << " minutes\n";
	
		double ts_Rentropy = wall_clock();
		outfile << "\n Renyi Entropy = " << Renyi_Entropy(dm);
		double te_Rentropy = wall_clock();
		outfile << "\n Time for Renyi Entropy " << (te_Rentropy - ts_Rentropy)/60 << " minutes\n\n";

		
		double ts_VNentropy = wall_clock();
		outfile << "\n Von-Neumann Entropy = " << Von_Neumann_Entropy(dm);
		double te_VNentropy = wall_clock();
		outfile << "\n Time for VN Entropy " << (te_VNentropy - ts_VNentropy)/60 << " minutes\n\n";
	
}

complex<double> VN_EE(	vector<ids> &psi,
			vector<int> &Cut,
			ofstream &outfile)
{
		int64_t n_dm_states = int64_t_power(2, Cut.size());
		zMatrix	dm(n_dm_states, n_dm_states, 0.0);
		density_matrix(dm, psi, Cut);
		return Von_Neumann_Entropy(dm);
	
}


#endif


