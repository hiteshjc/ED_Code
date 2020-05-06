#include "chirality_exp.h"

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

complex<double> ch_op(	int &i,
			int &j,
			int &k,
			int64_t &n1,
			vector<ids> &psi1,
			vector<ids> &psi2,
			ofstream &out)
{
	int64_t tmp_id1 = psi1[n1].stateID;
			
	bool loci = btest64(tmp_id1, i);			
	bool locj = btest64(tmp_id1, j);
	bool lock = btest64(tmp_id1, k);

//	out << i << "\t" << j << "\t" << k << endl;
//	out << loci << "\t" << locj << "\t" << lock << endl;

	complex< double > ch_exp = complex<double> (0.0, 0.0);

	if (loci != locj)
	{
		int sgn = 1;
		if (lock == 0) sgn = -1;
		int64_t n2 = -1;			

		int64_t tmp_id2 = tmp_id1;
		if(locj	== 0)
		{
			tmp_id2 = ibset64(tmp_id2, j);
			tmp_id2 = ibclr64(tmp_id2, i);
			
			n2 = findstate(tmp_id2, psi2);
			if(n2 >= 0)
			{	ch_exp = ch_exp - 0.25*sgn*conj(psi2[n2].coeff)*psi1[n1].coeff;}
		}
		else
		{
			tmp_id2 = ibset64(tmp_id2, i);
			tmp_id2 = ibclr64(tmp_id2, j);
						
			n2 = findstate(tmp_id2, psi2);
			if(n2 >= 0)
			{	ch_exp = ch_exp + 0.25*sgn*conj(psi2[n2].coeff)*psi1[n1].coeff;}

		}			
	
	}

	return ch_exp;
}


void act_chirality(int i, int j, int k, int64_t &id, vector<int64_t> &new_ids, vector<complex<double> > &hints)
{
	hints.clear();new_ids.clear();
	complex<double> im=complex<double>(0,1);	
	bool loci = btest64(id, i);			
	bool locj = btest64(id, j);
	bool lock = btest64(id, k);

	complex< double > ch_exp = complex<double> (0.0, 0.0);

	if (loci != locj)
	{
		int sgn = 1; if (lock == 0) sgn = -1;
		int64_t tmp_id2 = id;
		if(locj	== 0)
		{
			tmp_id2 = ibset64(tmp_id2, j);
			tmp_id2 = ibclr64(tmp_id2, i);
			new_ids.push_back(tmp_id2);
			hints.push_back(0.25*sgn/im);
		}
		else
		{
			tmp_id2 = ibset64(tmp_id2, i);
			tmp_id2 = ibclr64(tmp_id2, j);
			new_ids.push_back(tmp_id2);
			hints.push_back(-0.25*sgn/im);
		}			
	}
	
	if (lock != locj)
	{
		int sgn = 1; if (loci == 0) sgn = -1;
		int64_t tmp_id2 = id;
		if(lock	== 0)
		{
			tmp_id2 = ibset64(tmp_id2, k);
			tmp_id2 = ibclr64(tmp_id2, j);
			new_ids.push_back(tmp_id2);
			hints.push_back(0.25*sgn/im);
		}
		else
		{
			tmp_id2 = ibset64(tmp_id2, j);
			tmp_id2 = ibclr64(tmp_id2, k);
			new_ids.push_back(tmp_id2);
			hints.push_back(-0.25*sgn/im);
		}			
	}

	if (lock != loci)
	{
		int sgn = 1; if (locj == 0) sgn = -1;
		int64_t tmp_id2 = id;
		if(loci	== 0)
		{
			tmp_id2 = ibset64(tmp_id2, i);
			tmp_id2 = ibclr64(tmp_id2, k);
			new_ids.push_back(tmp_id2);
			hints.push_back(0.25*sgn/im);
		}
		else
		{
			tmp_id2 = ibset64(tmp_id2, k);
			tmp_id2 = ibclr64(tmp_id2, i);
			new_ids.push_back(tmp_id2);
			hints.push_back(-0.25*sgn/im);
		}			
	}

}

complex<double> ch_ch_op(	int &i,
				int &j,
				int &k,
				int &p,
				int &q,
				int &r,
				int64_t &n1,
				vector<ids> &psi1,
				vector<ids> &psi2,
				ofstream &out)
{
	// Unlike Krishna's ch_op code I sum over all terms here
	// Or I spit out all the n2's with the matrix element 
	
	int64_t tmp_id1 = psi1[n1].stateID;
	complex< double > ch_ch_exp = complex<double> (0.0, 0.0);
	vector<complex<double> > hints; 
	vector<int64_t> tmp_ids; 
	
	act_chirality(i,j,k,tmp_id1,tmp_ids,hints);

	for (int num_id=0;num_id<tmp_ids.size();num_id++)
	{
	 if (tmp_ids[num_id]!=-1) 
	 {
		vector<complex<double> > new_hints; 
		vector<int64_t> new_tmp_ids;
 
		act_chirality(p,q,r,tmp_ids[num_id],new_tmp_ids,new_hints);
		for (int n=0;n<new_tmp_ids.size();n++)
		{
			if (new_tmp_ids[n]!=-1)
			{
				int64_t n2=findstate(new_tmp_ids[n],psi2);
				if (n2>=0) ch_ch_exp+=hints[num_id]*new_hints[n]*psi1[n1].coeff*conj(psi2[n2].coeff);
			}
		}
	 }
	}
	return ch_ch_exp;
}

complex<double> curr_op(	int &i,
				int &j,
				int64_t &n1,
				vector<ids> &psi1,
				vector<ids> &psi2,
				ofstream &out)
{
	int64_t tmp_id1 = psi1[n1].stateID;
			
	bool loci = btest64(tmp_id1, i);			
	bool locj = btest64(tmp_id1, j);

	complex< double > curr = complex<double> (0.0, 0.0);

	if (loci != locj)
	{
		int64_t n2 = -1;			

		int64_t tmp_id2 = tmp_id1;
		if(locj	== 0)
		{
			tmp_id2 = ibset64(tmp_id2, j);
			tmp_id2 = ibclr64(tmp_id2, i);
			
			n2 = findstate(tmp_id2, psi2);
			if(n2 >= 0)
			{	curr = curr - 0.25*conj(psi2[n2].coeff)*psi1[n1].coeff;}
		}
		else
		{
			tmp_id2 = ibset64(tmp_id2, i);
			tmp_id2 = ibclr64(tmp_id2, j);
						
			n2 = findstate(tmp_id2, psi2);
			if(n2 >= 0)
			{	curr = curr + 0.25*conj(psi2[n2].coeff)*psi1[n1].coeff;}

		}			
	
	}

	return curr;
}

complex<double> szsz_op(	int &i,
				int &j,
				int64_t &n1,
				vector<ids> &psi1,
				vector<ids> &psi2,
				ofstream &out)
{
	int64_t tmp_id1 = psi1[n1].stateID;
			
	bool loci = btest64(tmp_id1, i);			
	bool locj = btest64(tmp_id1, j);

	complex< double > val = complex<double> (0.0, 0.0);
	val=val+(double(loci-0.5)*double(locj-0.5)*conj(psi2[n1].coeff)*psi1[n1].coeff);
	return val;
}

complex<double> sxsxsysy_op(	int &i,
				int &j,
				int64_t &n1,
				vector<ids> &psi1,
				vector<ids> &psi2,
				ofstream &out)
{
	int64_t tmp_id1 = psi1[n1].stateID;
			
	bool loci = btest64(tmp_id1, i);			
	bool locj = btest64(tmp_id1, j);

	complex< double > val = complex<double> (0.0, 0.0);

	if (loci != locj)
	{
		int64_t n2 = -1;			

		int64_t tmp_id2 = tmp_id1;
		if(locj	== 1)
		{
			tmp_id2 = ibclr64(tmp_id2, j);
			tmp_id2 = ibset64(tmp_id2, i);
			
			n2 = findstate(tmp_id2, psi2);
			if(n2 >= 0)
			{	val = val+conj(psi2[n2].coeff)*psi1[n1].coeff;}
		}
		if(locj	== 0)
		{
			tmp_id2 = ibset64(tmp_id2, j);
			tmp_id2 = ibclr64(tmp_id2, i);
			
			n2 = findstate(tmp_id2, psi2);
			if(n2 >= 0)
			{	val = val+conj(psi2[n2].coeff)*psi1[n1].coeff;}
		}
	}

	return 0.5*val;
}


complex<double> SipSjm(	int &i,
				int &j,
				int64_t &n1,
				vector<ids> &psi1,
				vector<ids> &psi2,
				ofstream &out)
{
	int64_t tmp_id1 = psi1[n1].stateID;
			
	bool loci = btest64(tmp_id1, i);			
	bool locj = btest64(tmp_id1, j);

	complex< double > val = complex<double> (0.0, 0.0);

	if (loci != locj)
	{
		int64_t n2 = -1;			

		int64_t tmp_id2 = tmp_id1;
		if(locj	== 1)
		{
			tmp_id2 = ibclr64(tmp_id2, j);
			tmp_id2 = ibset64(tmp_id2, i);
			
			n2 = findstate(tmp_id2, psi2);
			if(n2 >= 0)
			{	val = val+conj(psi2[n2].coeff)*psi1[n1].coeff;}
		}
	}

	return val;
}


