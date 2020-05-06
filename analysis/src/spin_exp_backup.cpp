#include "spin_exp.h"


int64_t findstate(	int64_t &rep_id,
	       		vector<ids> &Basis)
{
	int64_t j = -1;

	int64_t b_min = 0, b_max = Basis.size() - 1;
	do{
		int64_t b = b_min + (b_max - b_min)/2;
		if(rep_id < Basis[b].stateID )
			b_max = b - 1;
		else if (rep_id > Basis[b].stateID )
			b_min = b + 1;
		else
		{	j = b;	break;}
	}while(b_max >= b_min);

	return j;

}

vector< complex<double> > spin_spin_exp(	int &i, 
						int &j, 
						int64_t &n1, 
						vector<ids> &psi,
						ofstream &outfile)
{
	int64_t tmp_id1 = psi[n1].stateID;
			
	bool loci = btest64(tmp_id1, i);			
	bool locj = btest64(tmp_id1, j);
	
	complex<double> cz = complex<double> (0.0, 0.0);
	vector< complex< double > > sp_exp(2, cz);

	if (loci == locj)
	{
		sp_exp[0] = sp_exp[0] + 0.25*conj(psi[n1].coeff)*psi[n1].coeff;		
		sp_exp[1] = sp_exp[1] + 0.25*conj(psi[n1].coeff)*psi[n1].coeff;
	}
	else
	{
		sp_exp[0] = sp_exp[0] - 0.25*conj(psi[n1].coeff)*psi[n1].coeff;				
		sp_exp[1] = sp_exp[1] - 0.25*conj(psi[n1].coeff)*psi[n1].coeff;

		int64_t n2 = -1;
		int64_t tmp_id2 = tmp_id1;
		if(locj	== 0)
		{
			tmp_id2 = ibset64(tmp_id2, j);
			tmp_id2 = ibclr64(tmp_id2, i);
		}
		else
		{
			tmp_id2 = ibset64(tmp_id2, i);
			tmp_id2 = ibclr64(tmp_id2, j);
		}
					
		n2 = findstate(tmp_id2, psi);
		if(n2 >= 0)
		{	sp_exp[0] = sp_exp[0] + 0.5*conj(psi[n2].coeff)*psi[n1].coeff;}		
	
	}

	// Additional contribution for diagonal terms due to (S_x)^2 + (S_y)^2
	if (i == j)
	{
		int64_t n2 = -1;
		int64_t tmp_id2 = tmp_id1;
		n2 = findstate(tmp_id2, psi);
		if(n2 >= 0)
		{	sp_exp[0] = sp_exp[0] + 0.5*conj(psi[n2].coeff)*psi[n1].coeff;}
	}
	
	return sp_exp;
}


vector< complex<double> > bond_bond_exp(	links &i,
						links &j,
						int64_t &n1,
						vector<ids> &psi,
						ofstream &out)
{
	int64_t tmp_id1 = psi[n1].stateID;
			
	bool loc1 = btest64(tmp_id1, i.x);			
	bool loc2 = btest64(tmp_id1, i.y);
	bool loc3 = btest64(tmp_id1, j.x);
	bool loc4 = btest64(tmp_id1, j.y);

	complex<double> cz = complex<double> (0.0, 0.0);
	vector< complex< double > > bb_exp(2, cz);

	if (loc1 == loc2 && loc3 == loc4)
	{
//		bb_exp[0] = bb_exp[0] + 0.5*0.5*0.5*0.5*conj(psi[n1].coeff)*psi[n1].coeff;
		bb_exp[0] = bb_exp[0] + (0.25+0.5*0.25+0.5*0.25+0.25*0.25)*conj(psi[n1].coeff)*psi[n1].coeff;
		bb_exp[1] = bb_exp[1] + 0.25*0.25*conj(psi[n1].coeff)*psi[n1].coeff;
	}
	else if (loc1 == loc2 && loc3 != loc4)
	{
		// First add zz components
		bb_exp[0] = bb_exp[0] - 0.25*0.25*conj(psi[n1].coeff)*psi[n1].coeff;
		bb_exp[1] = bb_exp[1] - 0.25*0.25*conj(psi[n1].coeff)*psi[n1].coeff;

		int64_t n2 = -1;
		int64_t tmp_id2 = tmp_id1;
		if(loc4	== 0)
		{
			tmp_id2 = ibset64(tmp_id2, j.y);
			tmp_id2 = ibclr64(tmp_id2, j.x);
		}
		else
		{
			tmp_id2 = ibset64(tmp_id2, j.x);
			tmp_id2 = ibclr64(tmp_id2, j.y);
		}
			
		n2 = findstate(tmp_id2, psi);
		if(n2 >= 0)
		{
//			bb_exp[0] = bb_exp[0] + 0.125*conj(psi[n2].coeff)*psi[n1].coeff;
			bb_exp[0] = bb_exp[0] + 0.25*conj(psi[n2].coeff)*psi[n1].coeff;		
		}	
	}
	else if (loc1 != loc2 && loc3 == loc4)
	{
		bb_exp[0] = bb_exp[0] - 0.25*0.25*conj(psi[n1].coeff)*psi[n1].coeff;
		bb_exp[1] = bb_exp[1] - 0.25*0.25*conj(psi[n1].coeff)*psi[n1].coeff;
	
		int64_t n2 = -1;
		int64_t tmp_id2 = tmp_id1;
		if(loc2	== 0)
		{
			tmp_id2 = ibset64(tmp_id2, i.y);
			tmp_id2 = ibclr64(tmp_id2, i.x);
		}
		else
		{
			tmp_id2 = ibset64(tmp_id2, i.x);
			tmp_id2 = ibclr64(tmp_id2, i.y);
		}
			
		n2 = findstate(tmp_id2, psi);
		if(n2 >= 0)
		{
		//	bb_exp[0] = bb_exp[0] + 0.125*conj(psi[n2].coeff)*psi[n1].coeff;
			bb_exp[0] = bb_exp[0] + 0.25*conj(psi[n2].coeff)*psi[n1].coeff;
		}	
	}
	else
	{

	// SzSz terms
		bb_exp[0] = bb_exp[0] + 0.5*0.5*0.5*0.5*conj(psi[n1].coeff)*psi[n1].coeff;
		bb_exp[1] = bb_exp[1] + 0.5*0.5*0.5*0.5*conj(psi[n1].coeff)*psi[n1].coeff;

	// Flip bond i
		int64_t n2 = -1;
		int64_t tmp_id2 = tmp_id1;
		if(loc2	== 0)
		{	tmp_id2 = ibset64(tmp_id2, i.y);	tmp_id2 = ibclr64(tmp_id2, i.x);}
		else
		{	tmp_id2 = ibset64(tmp_id2, i.x);	tmp_id2 = ibclr64(tmp_id2, i.y);}
			
		n2 = findstate(tmp_id2, psi);
		if(n2 >= 0)
		{	bb_exp[0] = bb_exp[0] - 0.125*conj(psi[n2].coeff)*psi[n1].coeff;}	


	// Flip bond j		
		n2 = -1;
		tmp_id2 = tmp_id1;
		if(loc4	== 0)
		{	tmp_id2 = ibset64(tmp_id2, j.y);	tmp_id2 = ibclr64(tmp_id2, j.x);}
		else
		{	tmp_id2 = ibset64(tmp_id2, j.x);	tmp_id2 = ibclr64(tmp_id2, j.y);}
			
		n2 = findstate(tmp_id2, psi);
		if(n2 >= 0)
		{	bb_exp[0] = bb_exp[0] - 0.125*conj(psi[n2].coeff)*psi[n1].coeff;}


	// Flip both bonds i and j
		n2 = -1;
		tmp_id2 = tmp_id1;
		if(loc1 == 1 && loc2 == 0 && loc3 == 1 && loc4 == 0)
		{
			tmp_id2 = ibset64(tmp_id2, i.y);	tmp_id2 = ibclr64(tmp_id2, i.x);
			tmp_id2 = ibset64(tmp_id2, j.y);	tmp_id2 = ibclr64(tmp_id2, j.x);		
		}		
		else if (loc1 == 1 && loc2 == 0 && loc3 == 0 && loc4 == 1)
		{
			tmp_id2 = ibset64(tmp_id2, i.y);	tmp_id2 = ibclr64(tmp_id2, i.x);
			tmp_id2 = ibset64(tmp_id2, j.x);	tmp_id2 = ibclr64(tmp_id2, j.y);
		}
		else if (loc1 == 0 && loc2 == 1 && loc3 == 1 && loc4 == 0)
		{
			tmp_id2 = ibset64(tmp_id2, i.x);	tmp_id2 = ibclr64(tmp_id2, i.y);
			tmp_id2 = ibset64(tmp_id2, j.y);	tmp_id2 = ibclr64(tmp_id2, j.x);
		}
		else if (loc1 == 0 && loc2 == 1 && loc3 == 0 && loc4 == 1)
		{
			tmp_id2 = ibset64(tmp_id2, i.x);	tmp_id2 = ibclr64(tmp_id2, i.y);
			tmp_id2 = ibset64(tmp_id2, j.x);	tmp_id2 = ibclr64(tmp_id2, j.y);
		}
		n2 = findstate(tmp_id2, psi);
		if(n2 >= 0)
		{	bb_exp[0] = bb_exp[0] + 0.25*conj(psi[n2].coeff)*psi[n1].coeff;}		
	

	// If bond i and j are the same
		if (i.x == j.x && i.y == j.y)
			bb_exp[0] = bb_exp[0] + 0.25*conj(psi[n1].coeff)*psi[n1].coeff;
		
		if (i.x == j.y && i.y == j.x)
			bb_exp[0] = bb_exp[0] + 0.25*conj(psi[n1].coeff)*psi[n1].coeff;
			
		
	// If bond i and j share only one common site
		n2 = -1;
		tmp_id2 = tmp_id1;				
		if (i.x == j.x && i.y != j.y)
		{	
			if (loc2 == 0 && loc4 == 1)
			{	tmp_id2 = ibset64(tmp_id2, i.y);	tmp_id2 = ibclr64(tmp_id2, j.y);}
			else if (loc2 == 1 && loc4 == 0)
			{	tmp_id2 = ibclr64(tmp_id2, i.y);	tmp_id2 = ibset64(tmp_id2, j.y);}
		}
		else if (i.x == j.y && i.y != j.x)
		{	
			if (loc2 == 0 && loc3 == 1)
			{	tmp_id2 = ibset64(tmp_id2, i.y);	tmp_id2 = ibclr64(tmp_id2, j.x);}
			else if (loc2 == 1 && loc3 == 0)
			{	tmp_id2 = ibclr64(tmp_id2, i.y);	tmp_id2 = ibset64(tmp_id2, j.x);}
		}
		else if (i.x != j.x && i.y == j.y)
		{	
			if (loc1 == 0 && loc3 == 1)
			{	tmp_id2 = ibset64(tmp_id2, i.x);	tmp_id2 = ibclr64(tmp_id2, j.x);}
			else if (loc1 == 1 && loc3 == 0)
			{	tmp_id2 = ibclr64(tmp_id2, i.x);	tmp_id2 = ibset64(tmp_id2, j.x);}
		}
		else if (i.x != j.y && i.y == j.x)
		{	
			if (loc1 == 0 && loc4 == 1)
			{	tmp_id2 = ibset64(tmp_id2, i.x);	tmp_id2 = ibclr64(tmp_id2, j.y);}
			else if (loc1 == 1 && loc4 == 0)
			{	tmp_id2 = ibclr64(tmp_id2, i.x);	tmp_id2 = ibset64(tmp_id2, j.y);}
		}

		n2 = findstate(tmp_id2, psi);
		if (n2 >= 0)		
			bb_exp[0] = bb_exp[0] + 0.25*conj(psi[n2].coeff)*psi[n1].coeff;
	

	}
	return bb_exp;
}

