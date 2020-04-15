#include "H_v.h"

void sHv(	sMatrix &pBij, 
		complex< double > *v,
		complex< double > *w)
{
	int64_t n_p = pBij.NRows();
		
	#pragma omp parallel for
	for(int64_t i = 0; i < n_p; i++)
	{
		w[i] = 0;
		for(int64_t j = 0; j < pBij.RowSize(i); j++)
		{
			w[i] += pBij.me(i,j).coeff*v[pBij.me(i,j).stateID];      
		}
		
	} 
}

void makeH(	const vector<int64_t> &pblock_states,
		const vector<double> &norm,
		const vector<double> &norm_inv,
		const vector<int64_t> &rep_loc,
	        const vector< complex<double> > &state_ph,
		const double &lambda,
		const vector<coordinates> &adj_list,
		const double &Ch,
		const vector<triangles> &ch_list,
		const double &ChCh,
		const vector<bowties> &chch_list,
		const int &n_sites,
		const vector<int64_t> &Ia,
		const vector<int64_t> &Ib,
		const int &bits_right,
		const vector<int> &T1,
		const vector<int> &T2,
		const double &kx,
		const double &ky,
		ofstream &outfile, zMatrix &H)
{
	int64_t n_p = pblock_states.size();
	double Ezzfac = 0.25*lambda;

	H.resize(n_p,n_p);
	
        for(int64_t i = 0; i < n_p*n_p; i++) H[i]=0.0;

	//#pragma omp parallel for
	for(int64_t i = 0; i < n_p; i++)
	{
		int64_t idi = pblock_states[i];

		double Ezz = 0.0;
		for(int k = 0; k < adj_list.size(); k++)
		{
			int x = adj_list[k].x;
			int y = adj_list[k].y;
			double J = adj_list[k].J;
			complex<double> expphase=adj_list[k].expphase;
	
			if(abs(J) < 1e-16)
				continue;

			bool tx = btest64(idi, x);
			bool ty = btest64(idi, y);  
			
			if(tx == ty)
				Ezz += (Ezzfac*J);
			else
			{
				if (tx==1) expphase=conj(expphase);
				Ezz += - (Ezzfac*J); 
			
				int64_t tmp = flip(idi, x, y, tx, ty);

				int left, right;
				split_id(tmp, left, right, bits_right, n_sites);
				int64_t tmp_id_loc = Ia[right] + Ib[left] - 1;

				int64_t j = rep_loc[tmp_id_loc];

				// Need if statement to check if states rep_id actually exists
				// rep_id might have cancelled out due to addition of phases
				if(j >= 0)
					H(i,j) += expphase*J*0.5*state_ph[tmp_id_loc]*norm[i]*norm_inv[j];
					//w[i] += J*0.5*v[j]*state_ph[tmp_id_loc]*norm[i]*norm_inv[j];
			}
		}

		H(i,i) += Ezz;
		//w[i] += v[i]*Ezz;

		if(abs(Ch) > 1e-16)
		{
			for(int k = 0; k < ch_list.size(); k++)
			{
				int x = ch_list[k].i;
				int y = ch_list[k].j;			
				int z = ch_list[k].k;
				complex<double> ep =complex<double> (1.0,0.0);
				ep =ch_list[k].expphase;
				bool tx = btest64(idi, x);
				bool ty = btest64(idi, y);			
				
				if(tx != ty)
				{
					int64_t tmp = flip(idi, x, y, tx, ty);

					int left, right;
					split_id(tmp, left, right, bits_right, n_sites);
					int64_t tmp_id_loc = Ia[right] + Ib[left] - 1;

					int64_t j = rep_loc[tmp_id_loc];

					if(j >= 0)
					{
						complex< double > ch_sgn = complex< double > (0, 0);
						if (ty == 1)
							ch_sgn = complex< double> (0, 1) * ep;
						else 
							ch_sgn = complex< double> (0, -1) * conj(ep);
		
						bool tz = btest64(idi, z);
						if (tz == 0)
							ch_sgn = ch_sgn*complex <double> (-1, 0);			
					
						complex <double> ele_ij = Ch*0.25*ch_sgn*state_ph[tmp_id_loc]*norm[i]*norm_inv[j];

						H(i,j) += ele_ij;
						//w[i] += ele_ij*v[j];
					}	
				}
			}
		}
		
		if(abs(ChCh) > 1e-16)
		{
			for(int k = 0; k < chch_list.size(); k++)
			{
				int x = chch_list[k].x;
				int y = chch_list[k].y;			
				int z = chch_list[k].z;
				int p = chch_list[k].p;
				int q = chch_list[k].q;
				int r = chch_list[k].r;
	
				bool tx = btest64(idi, x);
				bool ty = btest64(idi, y);			
				complex< double > ch_sgn1 = complex< double > (0, 0);
				if (ty == 1) ch_sgn1 = complex< double> (0, 1);
				else 	     ch_sgn1 = complex< double> (0, -1);
				
				if(tx != ty)
				{
					int64_t tmp = flip(idi, x, y, tx, ty);
					bool tz = btest64(tmp, z);	
					bool tp = btest64(tmp, p);
					bool tq = btest64(tmp, q);			
					if (tz == 0) ch_sgn1 = ch_sgn1*complex <double> (-1, 0);			
					complex< double > ch_sgn2 = complex< double > (0, 0);
					if (tq == 1) ch_sgn2 = complex< double> (0, 1);
					else 	     ch_sgn2 = complex< double> (0, -1);
				
					if (tp !=tq)
					{
						int64_t tmp2 = flip(tmp, p, q, tp, tq);
						bool tr = btest64(tmp2, r);

						int left, right;
						split_id(tmp2, left, right, bits_right, n_sites);
						int64_t tmp_id_loc = Ia[right] + Ib[left] - 1;

						int64_t j = rep_loc[tmp_id_loc];

						if(j >= 0)
						{
							if (tr == 0) ch_sgn2 = ch_sgn2*complex <double> (-1, 0);			
							complex <double> ele_ij = ChCh*0.25*0.25*ch_sgn1*ch_sgn2*state_ph[tmp_id_loc]*norm[i]*norm_inv[j];
							H(i,j) += ele_ij;
							//w[i] += ele_ij*v[j];
						}
					}	
				}
			}
		}

	}	

}

void sHv(	complex< double > *v,
		complex< double > *w,
		const vector<int64_t> &pblock_states,
		const vector<double> &norm,
		const vector<double> &norm_inv,
		const vector<int64_t> &rep_loc,
	        const vector< complex<double> > &state_ph,
		const double &lambda,
		const vector<coordinates> &adj_list,
		const double &Ch,
		const vector<triangles> &ch_list,
		const double &ChCh,
		const vector<bowties> &chch_list,
		const int &n_sites,
		const vector<int64_t> &Ia,
		const vector<int64_t> &Ib,
		const int &bits_right,
		const vector<int> &T1,
		const vector<int> &T2,
		const double &kx,
		const double &ky,
		ofstream &outfile)
{
	int64_t n_p = pblock_states.size();
	double Ezzfac = 0.25*lambda;

	#pragma omp parallel for
	for(int64_t i = 0; i < n_p; i++)
	{
		w[i] = 0.0;
		int64_t idi = pblock_states[i];

		double Ezz = 0.0;
		for(int k = 0; k < adj_list.size(); k++)
		{
			int x = adj_list[k].x;
			int y = adj_list[k].y;
			double J = adj_list[k].J;
			complex<double> expphase=adj_list[k].expphase;
	
			if(abs(J) < 1e-16)
				continue;

			bool tx = btest64(idi, x);
			bool ty = btest64(idi, y);  
			
			if(tx == ty)
				Ezz += (Ezzfac*J);
			else
			{
				if (tx==1) expphase=conj(expphase);
				Ezz += - (Ezzfac*J); 
			
				int64_t tmp = flip(idi, x, y, tx, ty);

				int left, right;
				split_id(tmp, left, right, bits_right, n_sites);
				int64_t tmp_id_loc = Ia[right] + Ib[left] - 1;

				int64_t j = rep_loc[tmp_id_loc];

				// Need if statement to check if states rep_id actually exists
				// rep_id might have cancelled out due to addition of phases
				if(j >= 0)
					w[i] += expphase*J*0.5*v[j]*state_ph[tmp_id_loc]*norm[i]*norm_inv[j];
					//w[i] += J*0.5*v[j]*state_ph[tmp_id_loc]*norm[i]*norm_inv[j];
			}
		}

		w[i] += v[i]*Ezz;

		if(abs(Ch) > 1e-16)
		{
			for(int k = 0; k < ch_list.size(); k++)
			{
				int x = ch_list[k].i;
				int y = ch_list[k].j;			
				int z = ch_list[k].k;
				complex<double> ep =complex<double> (1.0,0.0);
				ep =ch_list[k].expphase;
				bool tx = btest64(idi, x);
				bool ty = btest64(idi, y);			
				
				if(tx != ty)
				{
					int64_t tmp = flip(idi, x, y, tx, ty);

					int left, right;
					split_id(tmp, left, right, bits_right, n_sites);
					int64_t tmp_id_loc = Ia[right] + Ib[left] - 1;

					int64_t j = rep_loc[tmp_id_loc];

					if(j >= 0)
					{
						complex< double > ch_sgn = complex< double > (0, 0);
						if (ty == 1)
							ch_sgn = complex< double> (0, 1) * ep;
						else 
							ch_sgn = complex< double> (0, -1) * conj(ep);
		
						bool tz = btest64(idi, z);
						if (tz == 0)
							ch_sgn = ch_sgn*complex <double> (-1, 0);			
					
						complex <double> ele_ij = Ch*0.25*ch_sgn*state_ph[tmp_id_loc]*norm[i]*norm_inv[j];

						w[i] += ele_ij*v[j];
					}	
				}
			}
		}
		
		if(abs(ChCh) > 1e-16)
		{
			for(int k = 0; k < chch_list.size(); k++)
			{
				int x = chch_list[k].x;
				int y = chch_list[k].y;			
				int z = chch_list[k].z;
				int p = chch_list[k].p;
				int q = chch_list[k].q;
				int r = chch_list[k].r;
	
				bool tx = btest64(idi, x);
				bool ty = btest64(idi, y);			
				complex< double > ch_sgn1 = complex< double > (0, 0);
				if (ty == 1) ch_sgn1 = complex< double> (0, 1);
				else 	     ch_sgn1 = complex< double> (0, -1);
				
				if(tx != ty)
				{
					int64_t tmp = flip(idi, x, y, tx, ty);
					bool tz = btest64(tmp, z);	
					bool tp = btest64(tmp, p);
					bool tq = btest64(tmp, q);			
					if (tz == 0) ch_sgn1 = ch_sgn1*complex <double> (-1, 0);			
					complex< double > ch_sgn2 = complex< double > (0, 0);
					if (tq == 1) ch_sgn2 = complex< double> (0, 1);
					else 	     ch_sgn2 = complex< double> (0, -1);
				
					if (tp !=tq)
					{
						int64_t tmp2 = flip(tmp, p, q, tp, tq);
						bool tr = btest64(tmp2, r);

						int left, right;
						split_id(tmp2, left, right, bits_right, n_sites);
						int64_t tmp_id_loc = Ia[right] + Ib[left] - 1;

						int64_t j = rep_loc[tmp_id_loc];

						if(j >= 0)
						{
							if (tr == 0) ch_sgn2 = ch_sgn2*complex <double> (-1, 0);			
							complex <double> ele_ij = ChCh*0.25*0.25*ch_sgn1*ch_sgn2*state_ph[tmp_id_loc]*norm[i]*norm_inv[j];
							w[i] += ele_ij*v[j];
						}
					}	
				}
			}
		}

	}	

}


void sHv(	complex< double > *v,
		complex< double > *w,
		const vector<int64_t> &pblock_states,
		const vector<double> &norm,
		const vector<double> &norm_inv,
		const vector<int64_t> &rep_loc,
	        const vector< complex<double> > &state_ph,
		const double &lambda,
		const vector<coordinates> &adj_list,
		const double &Ch,
		const vector<triangles> &ch_list,
		const double &ChCh,
		const vector<bowties> &chch_list,
		const int &n_sites,
		const vector<int64_t> &Ia,
		const vector<int64_t> &Ib,
		const int &bits_right,
		const vector<int> &T1,
		const vector<int> &T2,
		const vector<int> &T3,
		const double &kx,
		const double &ky,
		const double &kz,
		ofstream &outfile)
{
	int64_t n_p = pblock_states.size();
	double Ezzfac = 0.25*lambda;

	#pragma omp parallel for
	for(int64_t i = 0; i < n_p; i++)
	{
		w[i] = 0.0;
		int64_t idi = pblock_states[i];

		double Ezz = 0.0;
		for(int k = 0; k < adj_list.size(); k++)
		{
			int x = adj_list[k].x;
			int y = adj_list[k].y;
			double J = adj_list[k].J;
			complex<double> expphase=adj_list[k].expphase;
	
			if(abs(J) < 1e-16)
				continue;

			bool tx = btest64(idi, x);
			bool ty = btest64(idi, y);  
			
			if(tx == ty)
				Ezz += (Ezzfac*J);
			else
			{
				if (tx==1) expphase=conj(expphase);
				Ezz += - (Ezzfac*J); 
			
				int64_t tmp = flip(idi, x, y, tx, ty);

				int left, right;
				split_id(tmp, left, right, bits_right, n_sites);
				int64_t tmp_id_loc = Ia[right] + Ib[left] - 1;

				int64_t j = rep_loc[tmp_id_loc];

				// Need if statement to check if states rep_id actually exists
				// rep_id might have cancelled out due to addition of phases
				if(j >= 0)
					w[i] += expphase*J*0.5*v[j]*state_ph[tmp_id_loc]*norm[i]*norm_inv[j];
					//w[i] += J*0.5*v[j]*state_ph[tmp_id_loc]*norm[i]*norm_inv[j];
			}
		}

		w[i] += v[i]*Ezz;

		if(abs(Ch) > 1e-16)
		{
			for(int k = 0; k < ch_list.size(); k++)
			{
				int x = ch_list[k].i;
				int y = ch_list[k].j;			
				int z = ch_list[k].k;
				complex<double> ep =complex<double> (1.0,0.0);
				ep =ch_list[k].expphase;
				bool tx = btest64(idi, x);
				bool ty = btest64(idi, y);			
				
				if(tx != ty)
				{
					int64_t tmp = flip(idi, x, y, tx, ty);

					int left, right;
					split_id(tmp, left, right, bits_right, n_sites);
					int64_t tmp_id_loc = Ia[right] + Ib[left] - 1;

					int64_t j = rep_loc[tmp_id_loc];

					if(j >= 0)
					{
						complex< double > ch_sgn = complex< double > (0, 0);
						if (ty == 1)
							ch_sgn = complex< double> (0, 1) * ep;
						else 
							ch_sgn = complex< double> (0, -1) * conj(ep);
		
						bool tz = btest64(idi, z);
						if (tz == 0)
							ch_sgn = ch_sgn*complex <double> (-1, 0);			
					
						complex <double> ele_ij = Ch*0.25*ch_sgn*state_ph[tmp_id_loc]*norm[i]*norm_inv[j];

						w[i] += ele_ij*v[j];
					}	
				}
			}
		}
		
		if(abs(ChCh) > 1e-16)
		{
			for(int k = 0; k < chch_list.size(); k++)
			{
				int x = chch_list[k].x;
				int y = chch_list[k].y;			
				int z = chch_list[k].z;
				int p = chch_list[k].p;
				int q = chch_list[k].q;
				int r = chch_list[k].r;
	
				bool tx = btest64(idi, x);
				bool ty = btest64(idi, y);			
				complex< double > ch_sgn1 = complex< double > (0, 0);
				if (ty == 1) ch_sgn1 = complex< double> (0, 1);
				else 	     ch_sgn1 = complex< double> (0, -1);
				
				if(tx != ty)
				{
					int64_t tmp = flip(idi, x, y, tx, ty);
					bool tz = btest64(tmp, z);	
					bool tp = btest64(tmp, p);
					bool tq = btest64(tmp, q);			
					if (tz == 0) ch_sgn1 = ch_sgn1*complex <double> (-1, 0);			
					complex< double > ch_sgn2 = complex< double > (0, 0);
					if (tq == 1) ch_sgn2 = complex< double> (0, 1);
					else 	     ch_sgn2 = complex< double> (0, -1);
				
					if (tp !=tq)
					{
						int64_t tmp2 = flip(tmp, p, q, tp, tq);
						bool tr = btest64(tmp2, r);

						int left, right;
						split_id(tmp2, left, right, bits_right, n_sites);
						int64_t tmp_id_loc = Ia[right] + Ib[left] - 1;

						int64_t j = rep_loc[tmp_id_loc];

						if(j >= 0)
						{
							if (tr == 0) ch_sgn2 = ch_sgn2*complex <double> (-1, 0);			
							complex <double> ele_ij = ChCh*0.25*0.25*ch_sgn1*ch_sgn2*state_ph[tmp_id_loc]*norm[i]*norm_inv[j];
							w[i] += ele_ij*v[j];
						}
					}	
				}
			}
		}

	}	

}

