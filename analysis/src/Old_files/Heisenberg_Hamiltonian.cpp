#include "Heisenberg_Hamiltonian.h"

// Function to find the representative given a Basis stateID
// Function also returns the number of hops (lx, ly) need to get to he representative state
int64_t representative(	int &lx,
			int &ly,
			int &L1,
			int &L2,
			int64_t &v_id,
			vector<int> &T1,
			vector<int> &T2)
{
	lx = 0, ly = 0;
	int64_t rep_id = v_id;
	int64_t tmp_id = v_id;
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
			}
			}
	    	}
	}
	
	return rep_id;
}

// Locates the representative pBasis state using binary search
int64_t findstate(	int64_t &rep_id,
	       		vector<int64_t> &pBasis)
{
	int64_t j = -1;

	int64_t b_min = 0, b_max = pBasis.size() - 1;
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

void findLxLy(	int &Lx,
		int &Ly,
		vector<int> &T1,
	      	vector<int> &T2)
{
	int64_t tmp = 1;

	do{
		translateT(tmp, T1);
		Lx++;
	}while(tmp!= 1);

	tmp = 1;

	do{
		translateT(tmp, T2);
		Ly++;
	}while(tmp!= 1);
}


void zXXZ_Heisenberg(	zMatrix &pBij,
			int &n_spins, 
			vector<int64_t> &pblock_states,
			vector<double> &norm,
		   	double &kx,
			double &ky,
			double &lambda,
			vector<coordinates> &adj_list,
			double &Ch,
			vector<triangles> &ch_list,
			vector<int> &T1,
			vector<int> &T2)
{
	int64_t idi, idj;
	int Lx = 0, Ly = 0;
	findLxLy(Lx, Ly, T1, T2);

	for(int64_t i = 0; i < pblock_states.size(); i++)
	{
		idi = pblock_states[i];
		double Ni = norm[i]; //= findNormalization(idi, T1, T2, Lx, Ly, pblock_states);

		// Compute the Heisenberg Hamiltonian	
		vector<int64_t> id;
		vector< complex <double> > coeff;
		double Ezz = 0.0;

		for(int k = 0; k < adj_list.size(); k++)
		{
			int x = adj_list[k].x;
			int y = adj_list[k].y;
			double J = adj_list[k].J;

			bool tx = btest64(idi, x);
			bool ty = btest64(idi, y);

			if(tx == ty)
				Ezz = Ezz + 0.25*lambda;
			else
			{
				Ezz = Ezz - 0.25*lambda;   

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

				if(j >= i && j >= 0)
				{
					complex<double> argx = complex< double >(0,-kx*lx), argy = complex< double >(0, -ky*ly);
					complex<double> phx = exp(argx), phy = exp(argy);

					double Nj = norm[j]; // findNormalization(rep_id, T1, T2, Lx, Ly, pblock_states);
					pBij(i,j) = pBij(i,j) + J*0.5*phx*phy*Ni/Nj;//sqrt(Ni/Nj);
					if(j != i)
						pBij(j,i) = pBij(j,i) + J*0.5*conj(phx*phy)*Ni/Nj;//sqrt(Ni/Nj);
				}	
			}
		}

		pBij(i,i) = pBij(i,i) + Ezz;
	
		// Compute the Chirality Hamiltonian
		if(Ch != 0.0)
		for(int k = 0; k < ch_list.size(); k++)
		{
			int x = ch_list[k].i;
			int y = ch_list[k].j;			
			int z = ch_list[k].k;

			bool tx = btest64(idi, x);
			bool ty = btest64(idi, y);			
			
			if(tx != ty)
			{
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

				if(j >= i && j >= 0)
				{
					complex< double > ch_sgn = complex< double > (0, 0);

					if (ty == 1)
						ch_sgn = complex< double> (0, 1);
					else 
						ch_sgn = complex< double> (0, -1);
	
					bool tz = btest64(idi, z);
		
					if (tz == 0)
						ch_sgn = ch_sgn*complex <double> (-1, 0);			
				
					complex<double> argx = complex< double >(0,-kx*lx), argy = complex< double >(0, -ky*ly);
					complex<double> phx = exp(argx), phy = exp(argy);

					double Nj = norm[j]; // findNormalization(rep_id, T1, T2, Lx, Ly, pblock_states);
					complex <double> ele_ij = Ch*0.25*ch_sgn*phx*phy*Ni/Nj;
 
					pBij(i,j) = pBij(i,j) + ele_ij;//sqrt(Ni/Nj);
					if(j != i)
						pBij(j,i) = pBij(j,i) + conj(ele_ij);//sqrt(Ni/Nj);
				}	
			}
		}
	}


}

void sXXZ_Heisenberg(	sMatrix &pBij,
			int &n_spins, 
			vector<int64_t> &pblock_states,
			vector<double> &norm,
			double &kx,
			double &ky,
			double &lambda,
			vector<coordinates> &adj_list,
			double &Ch,
			vector<triangles> &ch_list,
			vector<int> &T1,
			vector<int> &T2,
			ofstream &outfile)
{
	int64_t idi, idj;
	int Lx = 0, Ly = 0;
	findLxLy(Lx, Ly, T1, T2);
	int64_t n = pBij.NRows();

	double t_tot_findstate = 0.0;
	double t_tot_rep = 0.0;
	int64_t n_rep = 0;

	for(int64_t i = 0; i < n; i++)
	{
		idi = pblock_states[i];
		double Ni = norm[i]; //= findNormalization(idi, T1, T2, Lx, Ly, pblock_states);
		vector<ids> w;

		vector< complex <double> > coeff;

		double Ezz = 0.0;
		complex< double > eleij = complex< double >(0.0, 0.0);

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
				Ezz = Ezz + 0.25*lambda;
			else
			{
				Ezz = Ezz - 0.25*lambda;   

				int64_t tmp1 = flip(idi, x, y, tx, ty);

				int lx, ly;
				clock_t t_rep;
				t_rep = clock();
				int64_t rep_id = representative(lx, ly, Lx, Ly, tmp1, T1, T2);       
				t_rep = clock() - t_rep;
				t_tot_rep = t_tot_rep + t_rep;
				n_rep++;
				lx = lx % Lx;
				ly = ly % Ly;

				if(lx < 0)
					lx = lx + Lx;
				if(ly < 0)
					ly = ly + Ly;

				clock_t t_findstate;
				t_findstate = clock();
				int64_t j = findstate(rep_id, pblock_states);
				t_findstate = clock() - t_findstate;
				t_tot_findstate = t_tot_findstate + t_findstate;

				if(j >= i && j >= 0)
				{
					  complex<double> argx = complex< double >(0,-kx*lx), argy = complex< double >(0, -ky*ly);
					  complex<double> phx = exp(argx), phy = exp(argy);
					  double Nj = norm[j]; // findNormalization(rep_id, T1, T2, Lx, Ly, pblock_states);
					  complex<double> coeff = J*0.5*phx*phy*Ni/Nj;
					  
					  add2w(w, j, coeff);
				}
			}
		}

		if(Ezz != 0.0)
			add2w(w, i, Ezz);

		// Compute the Chirality Hamiltonian
		if(Ch != 0.0)
		for(int k = 0; k < ch_list.size(); k++)
		{
			int x = ch_list[k].i;
			int y = ch_list[k].j;			
			int z = ch_list[k].k;

			bool tx = btest64(idi, x);
			bool ty = btest64(idi, y);			
			
			if(tx != ty)
			{
				int64_t tmp1 = flip(idi, x, y, tx, ty);

				int lx, ly;
				clock_t t_rep;
				t_rep = clock();
				int64_t rep_id = representative(lx, ly, Lx, Ly, tmp1, T1, T2);       
				t_rep = clock() - t_rep;
				t_tot_rep = t_tot_rep + t_rep;
				n_rep++;
				lx = lx % Lx;
				ly = ly % Ly;


				if(lx < 0)
				  lx = lx + Lx;
				if(ly < 0)
				  ly = ly + Ly;

				
				clock_t t_findstate;
				t_findstate = clock();
				int64_t j = findstate(rep_id, pblock_states);		
				t_findstate = clock() - t_findstate;
				t_tot_findstate = t_tot_findstate + t_findstate;

				if(j >= i && j >= 0)
				{
					complex<double> ch_sgn = complex<double>(0, 0);

					if (ty == 1)
						ch_sgn = complex< double> (0, Ch*0.25);
					else 
						ch_sgn = complex< double> (0, -Ch*0.25);
	
					bool tz = btest64(idi, z);
		
					if (tz == 0)
						ch_sgn = -1.0*ch_sgn;			
				
					complex<double> argx = complex< double >(0,-kx*lx), argy = complex< double >(0, -ky*ly);
					complex<double> phx = exp(argx), phy = exp(argy);

					double Nj = norm[j]; // findNormalization(rep_id, T1, T2, Lx, Ly, pblock_states);
				//	complex <double> coeff = Ch*0.25*ch_sgn*phx*phy*Ni/Nj;
 					complex <double> coeff = ch_sgn*phx*phy*Ni/Nj;
					add2w(w, j, coeff);
				}	
			}
		}

		for(int64_t l = 0; l < w.size(); l++)
		{
			complex< double > tmp_coeff = w[l].coeff;
			int64_t tmp_id = w[l].stateID;

			if(abs(real(tmp_coeff)) > 1e-8 || abs(imag(tmp_coeff)) > 1e-8)
			{
				pBij.push_back(i, tmp_id, tmp_coeff);
				if(i != tmp_id)
				  pBij.push_back(tmp_id, i, conj(tmp_coeff));
			}
		}

	}  
/*
	for(int ii = 0; ii < pBij.NRows(); ii++)
	{
		vector<ids> tmp = pBij.return_sVector(ii);
		for(int jj = 0; jj < tmp.size(); jj++)
			outfile << tmp[jj].stateID << " " << tmp[jj].coeff << "\t";
		outfile << endl;	
	}
*/
	outfile << "\nHamiltonian Construction\n";
	outfile << "Total # of rep calls \t\t" << n_rep << endl;
	outfile << "Total T for representative \t" << double(t_tot_rep)/double(CLOCKS_PER_SEC)/60 << " minutes" << endl;

	outfile << "Total T for findstate \t\t" << double(t_tot_findstate)/double(CLOCKS_PER_SEC)/60 << " minutes" << endl;

}

