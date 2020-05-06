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


void sHv(	complex< double > *v,
		complex< double > *w,
		vector<int64_t> &pblock_states,
		vector<double> &Norm,
		double &kx,
		double &ky,
		double &lambda,
		vector<coordinates> &adj_list,
		vector<int> &T1,
		vector<int> &T2,
		ofstream &outfile)
{
	int Lx = 0, Ly = 0;
	findLxLy(Lx, Ly, T1, T2);
	int64_t n_p = pblock_states.size();

	for(int64_t i = 0; i < n_p; i++)
		w[i] = 0.0;

	#pragma omp parallel for
	for(int64_t i = 0; i < n_p; i++)
	{
		int64_t idi = pblock_states[i];
		double Ni = Norm[i]; //findNormalization(idi, T1, T2, Lx, Ly, pblock_states);
		
		double Ezz = 0.0;

		for(int k = 0; k < adj_list.size(); k++)
		{
			if(adj_list[k].J == 0.0)
				continue;
	
			int x = adj_list[k].x;
			int y = adj_list[k].y;
			double J = adj_list[k].J;

			bool tx = btest64(idi, x);
			bool ty = btest64(idi, y);  


			if(tx == ty)
				Ezz = Ezz + J*0.25*lambda;
			else
			{
				Ezz = Ezz - J*0.25*lambda; 

				int64_t tmp = flip(idi, x, y, tx, ty);

				int lx, ly;
				int64_t rep_id = representative(lx, ly, Lx, Ly, tmp, T1, T2);
				lx = lx % Lx;
				ly = ly % Ly;

				if(lx < 0)
					lx = lx + Lx;
				if(ly < 0)
					ly = ly + Ly;

				int64_t j = findstate(rep_id, pblock_states);
			
				if(j >= 0)
				{
					complex<double> argx = complex< double >(0,-kx*lx), argy = complex< double >(0, -ky*ly);
					complex<double> phx = exp(argx), phy = exp(argy);

					double Nj = Norm[j]; //findNormalization(idj, T1, T2, Lx, Ly, pblock_states);
					w[i] += J*0.5*v[j]*phx*phy*Ni/Nj;
				}	
			}
		}
		w[i] += v[i]*Ezz;
	}	
}



/*
void sHv(	complex< double > *v,
		complex< double > *w,
		vector<int64_t> &pblock_states,
		vector<double> &Norm,
		double &kx,
		double &ky,
		double &lambda,
		vector<coordinates> &adj_list,
		vector<int> &T1,
		vector<int> &T2,
		ofstream &outfile)
{
	int Lx = 0, Ly = 0;
	findLxLy(Lx, Ly, T1, T2);
	int64_t n_p = pblock_states.size();

	for(int64_t i = 0; i < n_p; i++)
		w[i] = 0.0;

	#pragma omp parallel for
	for(int64_t i = 0; i < n_p; i++)
	{
		int64_t idi = pblock_states[i];
		double Ni = Norm[i]; //findNormalization(idi, T1, T2, Lx, Ly, pblock_states);
		complex< double > v_coeff = v[i];

		double Ezz = 0.0;

		for(int k = 0; k < adj_list.size(); k++)
		{

			int x = adj_list[k].x;
			int y = adj_list[k].y;
			double J = adj_list[k].J;

			bool tx = btest64(idi, x);
			bool ty = btest64(idi, y);  


			if(tx == ty)
				Ezz = Ezz + J*0.25*lambda;
			else
			{
				Ezz = Ezz - J*0.25*lambda; 

				int64_t tmp = flip(idi, x, y, tx, ty);

				int lx, ly;
				int64_t rep_id = representative(lx, ly, Lx, Ly, tmp, T1, T2);
				lx = lx % Lx;
				ly = ly % Ly;

				if(lx < 0)
					lx = lx + Lx;
				if(ly < 0)
					ly = ly + Ly;

				int64_t j = findstate(rep_id, pblock_states);

				if(j >= 0)
				{
					complex<double> argx = complex< double >(0,-kx*lx), argy = complex< double >(0, -ky*ly);
					complex<double> phx = exp(argx), phy = exp(argy);

					double Nj = Norm[j]; //findNormalization(idj, T1, T2, Lx, Ly, pblock_states);
					w[j] += v_coeff*J*0.5*phx*phy*Nj/Ni;
				}	
			}
		}
		w[i] += v_coeff*Ezz;
	}	
}

*/

