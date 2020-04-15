
#include "lanczos.h"
#include <boost/format.hpp>
#include <ctime>

using namespace std;

void zscalk(const int64_t size,
	    complex<double> a,
	    vector< complex< double> > &x)
{
	#pragma omp parallel for
	for(int64_t i = 0; i < size; ++i) {x[i]*=(a);}

}


void zaxpyk(const int64_t size,
	    complex<double> a,
	    vector< complex< double> > &x,
	    vector< complex< double> > &y)
{
	#pragma omp parallel for
	for(int64_t i = 0; i < size; ++i) {y[i]+=(a*x[i]);}

}


complex<double> zdotc(	const int64_t &size,
			const vector< complex< double> > &v1,
			const vector< complex< double> > &v2)
{
	double sumr = 0.0;
	double sumi = 0.0;

	#pragma omp parallel for default(shared) reduction (+ : sumr,sumi)
	for(int64_t i = 0; i < size; ++i)
	{
		complex<double> a=conj(v1[i])*v2[i];
		sumr += real(a);
		sumi += imag(a);
	}

	complex<double> sum=complex<double>(sumr,sumi);

	return sum;

}

void normalize(std::vector< complex<double> > & v) 
{
	int64_t size=v.size();
	complex<double> norminv=1.0/sqrt(real(zdotc(size, v, v)));
	zscalk(size,norminv,v);
}

void orth_wrt_previous_evecs(ofstream &outfile, 
		       std::vector< complex<double> > & v, 
		       std::vector< std::vector< complex<double> > > & previous_evecs)
{	
	int64_t size=v.size();
	std::vector< complex<double> > qs;

	for (int i=0; i < previous_evecs.size();i++)
	{
		complex<double> q=conj(zdotc(size, v, previous_evecs[i]));
		qs.push_back(q);
	}
	
	for (int i=0; i < previous_evecs.size();i++)
	{
		zaxpyk(size,-qs[i],previous_evecs[i],v);
	}
}

double unif_rnd()
{
	return rand()/double(RAND_MAX);
}

void exact_evecs(	vector<int64_t> &pblock_states,
			vector<double> &Norm,
			vector<double> &Norm_inv,
			vector<int64_t> &rep_loc, 
			vector< complex<double> > &state_ph,
			vector<complex<double> > &eigs,
			vector< vector< complex<double> > > &previous_evecs,
			bool find_ev,
			double &lambda,
			int &n_spins,
			double &px, 
			double &py,
			int &n_ones,
			int64_t &n_p,
			vector<coordinates> &adj_list,
			double &Ch,
			vector<triangles> &ch_list,
			double &ChCh,
			vector<bowties> &chch_list,
			vector<int> &T1,
			vector<int> &T2,
			int n_ev,
			double tol,
			int iterations,
			int ncycles,
			vector<int64_t> &Ia,
			vector<int64_t> &Ib,
			int bits_right,
			ofstream &outfile)
{
	bool ipr = true;
	int64_t size = pblock_states.size();
	vector< double > 			reigs(size);//v_o(size);
	vector< complex<double> > 		v_p(size);//v_o(size);
	zMatrix					H_mat,H_eigenvecs;
	makeH(pblock_states, Norm, Norm_inv, rep_loc, state_ph, lambda, adj_list, Ch, ch_list, ChCh, chch_list, n_spins, Ia, Ib, bits_right, T1, T2, px, py, outfile, H_mat);
	zMatrix_Diagonalize_Hermitian(H_mat, reigs, H_eigenvecs);
	eigs.clear(); for (int b=0;b<reigs.size();b++) eigs.push_back(reigs[b]);
		
	for (int i=0;i<size;i++)
	{		
		for (int k=0;k<size;k++)
		{
			v_p[k]=H_eigenvecs(k,i);
		}
		previous_evecs.push_back(v_p);
	}

	outfile << "\n\n--------------------------------------------------\n";
	outfile << "\n EXACT SUMMARY \n" << endl;
}

void lanczos_evecs(	vector<int64_t> &pblock_states,
			vector<double> &Norm,
			vector<double> &Norm_inv,
			vector<int64_t> &rep_loc, 
			vector< complex<double> > &state_ph,
			vector<complex<double> > &eigs,
			vector< vector< complex<double> > > &previous_evecs,
			bool find_ev,
			double &lambda,
			int &n_spins,
			double &px, 
			double &py,
			int &n_ones,
			int64_t &n_p,
			vector<coordinates> &adj_list,
			double &Ch,
			vector<triangles> &ch_list,
			double &ChCh,
			vector<bowties> &chch_list,
			vector<int> &T1,
			vector<int> &T2,
			int n_ev,
			double tol,
			int iterations,
			int ncycles,
			vector<int64_t> &Ia,
			vector<int64_t> &Ib,
			int bits_right,
			ofstream &outfile)
{
	int how_many_evecs = min(n_ev,int(pblock_states.size()));
	bool ipr = true;
	int base;

	int64_t size = pblock_states.size();
	double spin;
	bool orth_failed;

	double 					dif,tmp;
	complex< double	>			q,alpha,beta,norm;
	vector< complex<double> > 		alphas,betas;
	vector< complex<double> >		w(size);
	vector< complex<double> > 		v_p(size);//v_o(size);
	vector< vector< complex<double> > > 	vs;
	zMatrix					t_mat,t_eigenvecs;
	vector< complex<double> > 		all_eigs;
	complex<double>				previous_eig;
	int Hv_calls = 0;
	double time_Hv_total = 0.0;
	double time_ortho_total = 0.0;
	int 					wait_cycles=0;

	//iterations=min(iterations,int(size));
	for (int num = 0; num < how_many_evecs; num++)
	{
		iterations=min(iterations,int(size)-num);
		wait_cycles=0;
		//if (num==0)
		//{
		for (int64_t i = 0; i < size; i++) v_p[i]=complex<double>(2.0*unif_rnd()-1.0,2.0*unif_rnd()-1.0);
		//wait_cycles=30;
		//}
		previous_eig=1000.0;
		for (int cycle=0;cycle<ncycles;cycle++)
		{

		       orth_wrt_previous_evecs(outfile,v_p,previous_evecs);
	  	       normalize(v_p);
		       //zscalk(size,complex<double>(0.0,0.0),v_o);
	  	       zscalk(size,complex<double>(0.0,0.0),w);

			vs.clear();
			beta=0.0;
			betas.clear();
			betas.push_back(beta);
			alphas.clear();

			if (ipr)
			{
				outfile<<"iterations = "<<iterations<<endl;
				outfile<<"vs.size()  = "<<vs.size()<<endl;
			}
			outfile.flush();

			for (int it=0;it<iterations;it++)
			{
			       	if (ipr)
			       	{
					outfile<<"================================================================="<<endl;
					outfile<<"Doing Lanczos iteration,cycle,evec = "<<it<<" "<<cycle<<"  "<<num<<endl;
				}

			       	//vs.push_back(v_p);     
			       	
				/////////////////////////
				////////// Hv ///////////
				/////////////////////////
				
				double time_Hv = wall_clock(); Hv_calls++;
				sHv(&*v_p.begin(), &*w.begin(), pblock_states, Norm, Norm_inv, rep_loc, state_ph, lambda, adj_list, Ch, ch_list, ChCh, chch_list, n_spins, Ia, Ib, bits_right, T1, T2, px, py, outfile);
				outfile << "\n\t Wall clock for Hv step: " << (wall_clock() - time_Hv)/60.0 << endl;
				outfile.flush();
				time_Hv_total += wall_clock() - time_Hv;
			       	//if (it!=0) zaxpyk(size,-beta,v_o,w);
			       	if (it!=0) zaxpyk(size,-beta,vs[it-1],w);
				alpha=conj(zdotc(size, w, v_p));
			       	alphas.push_back(alpha);
			       	zaxpyk(size,-alpha,v_p,w);
			       	vs.push_back(v_p);     
			       	//v_o=v_p;
			       	beta=sqrt(real(zdotc(size, w, w)));
			       	#pragma omp parallel for 
				for (int64_t i=0;i<size;i++) v_p[i]=w[i];

			       	zscalk(size,1.0/beta,v_p);
			       	betas.push_back(beta);
			       	zscalk(size,0.0,w);

				//////////////////////////////////////////////////////////////////
			       	// Orthogonlize v_p with respect to all previous eigenvectors ////

				double time_ortho = wall_clock();
			       	orth_wrt_previous_evecs(outfile,v_p,previous_evecs);
				if (cycle==ncycles-1)
				{
					// Bad guess needs no orth
			       		orth_wrt_previous_evecs(outfile,v_p,vs);
			       	}
				normalize(v_p);

			       	if (ipr)
			       	{
					outfile <<"Time to perform Lanczos iteration "<<it<<" was "<<dif<<" seconds"<<endl;
					outfile <<"================================================================="<<endl;
			       	}

				if (ipr) {outfile <<"vs.size() = "<<vs.size()<<endl;}

				int tmp=alphas.size();
				t_mat.resize(tmp,tmp);
				t_eigenvecs.resize(tmp,tmp);
				for (int j=0;j<tmp*tmp;j++){t_mat[j]=0.0;t_eigenvecs[j]=0.0;}
				for (int j=0;j<tmp;j++)
				{
					t_mat(j,j)=alphas[j];
					if (j+1<tmp)
					{t_mat(j,j+1)=betas[j+1];t_mat(j+1,j)=betas[j+1];}
				}

//				if (ipr) {outfile <<"Time to build  T was "<<dif<<" seconds"<<endl;}
				
				std::vector<double> reigs(tmp);
				zMatrix_Diagonalize_Hermitian(t_mat, reigs, t_eigenvecs);
				eigs.clear(); for (int b=0;b<reigs.size();b++) eigs.push_back(reigs[b]);
				
//				if (ipr) {outfile<<"Time to diagonalize T was "<<dif<<" seconds"<<endl;}
				
				if (cycle!=0 and cycle!=ncycles-1 and it>1) {
					outfile <<"Eigs[0]      = "<<boost::format("%+.15f") %eigs[0]<<endl;
					outfile <<"Previous Eig = "<<boost::format("%+.15f") %previous_eig<<endl;
					if (abs(eigs[0]-previous_eig)<tol) {outfile <<"Convergence, doing last cycle with orth !!"<<endl;it=iterations;cycle=ncycles-2;}
					previous_eig=eigs[0];
				}
				else if (cycle==ncycles-1) {
					outfile <<"Final cycle, Eigs[0]      = "<<boost::format("%+.15f") %eigs[0]<<endl;
				}
				else{
					outfile <<"Eigs[0]      = "<<boost::format("%+.15f") %eigs[0]<<endl;
					outfile <<"Previous Eig = "<<boost::format("%+.15f") %previous_eig<<endl;
					previous_eig=eigs[0];
				}	

				outfile << "\n\t Wall clock for orthogonalize " << (wall_clock() - time_ortho)/60.0 << endl;
				time_ortho_total += wall_clock() - time_ortho;
				//////////////////////////////////////////////////////////////////////

			}

			if (ipr) {outfile <<"Making lowest Ritz eigenvector"<<endl;}
			#pragma omp parallel for
			for (int i=0;i<size;i++)
			{
				v_p[i]=0.0;
				for (int k=0;k<vs.size();k++)
				{v_p[i]+=vs[k][i]*t_eigenvecs(k,0);}
			}
		       orth_wrt_previous_evecs(outfile,v_p,previous_evecs);
		       normalize(v_p);
		}	
		previous_evecs.push_back(v_p);
		/*if (v_p.size()<30 and )
		{
			for (int i=0;i<v_p.size();i++)
			{
				
			}
		}*/
		/*// Use first excited state as start guess for next vector
		#pragma omp parallel for
		for (int i=0;i<size;i++)
		{
			v_p[i]=0.0;
			for (int k=0;k<vs.size();k++)
			{v_p[i]+=vs[k][i]*t_eigenvecs(k,1);}
		}*/
			
		outfile << "Eigenvalue number " << num << " = " << boost::format("%+.15f") %eigs[0] << endl;
		all_eigs.push_back(eigs[0]);
	}
	eigs=all_eigs;

	outfile << "\n\n--------------------------------------------------\n";
	outfile << "\n Lanczos SUMMARY \n" << endl;
	outfile << " Hv calls\t\t\t" << Hv_calls << endl;
	outfile << " Hv total time \t\t\t" << time_Hv_total/60.0 << " minutes"<< endl;
	outfile << " Orthogonalize total time \t" << time_ortho_total/60.0 << " minutes" << endl;
}

void lanczos_evecs(	vector<int64_t> &pblock_states,
			vector<double> &Norm,
			vector<double> &Norm_inv,
			vector<int64_t> &rep_loc, 
			vector< complex<double> > &state_ph,
			vector<complex<double> > &eigs,
			vector< vector< complex<double> > > &previous_evecs,
			bool find_ev,
			double &lambda,
			int &n_spins,
			double &px, 
			double &py,
			double &pz,
			int &n_ones,
			int64_t &n_p,
			vector<coordinates> &adj_list,
			double &Ch,
			vector<triangles> &ch_list,
			double &ChCh,
			vector<bowties> &chch_list,
			vector<int> &T1,
			vector<int> &T2,
			vector<int> &T3,
			int n_ev,
			double tol,
			int iterations,
			int ncycles,
			vector<int64_t> &Ia,
			vector<int64_t> &Ib,
			int bits_right,
			ofstream &outfile)
{
	int how_many_evecs = min(n_ev,int(pblock_states.size()));
	bool ipr = true;
	int base;

	int64_t size = pblock_states.size();
	double spin;
	bool orth_failed;

	double 					dif,tmp;
	complex< double	>			q,alpha,beta,norm;
	vector< complex<double> > 		alphas,betas;
	vector< complex<double> >		w(size);
	vector< complex<double> > 		v_p(size);//v_o(size);
	vector< vector< complex<double> > > 	vs;
	zMatrix					t_mat,t_eigenvecs;
	vector< complex<double> > 		all_eigs;
	complex<double>				previous_eig;
	int Hv_calls = 0;
	double time_Hv_total = 0.0;
	double time_ortho_total = 0.0;
	int 					wait_cycles=0;

	//iterations=min(iterations,int(size));
	for (int num = 0; num < how_many_evecs; num++)
	{
		iterations=min(iterations,int(size)-num);
		wait_cycles=0;
		//if (num==0)
		//{
		for (int64_t i = 0; i < size; i++) v_p[i]=complex<double>(2.0*unif_rnd()-1.0,2.0*unif_rnd()-1.0);
		//wait_cycles=30;
		//}
		previous_eig=1000.0;
		for (int cycle=0;cycle<ncycles;cycle++)
		{

		       orth_wrt_previous_evecs(outfile,v_p,previous_evecs);
	  	       normalize(v_p);
		       //zscalk(size,complex<double>(0.0,0.0),v_o);
	  	       zscalk(size,complex<double>(0.0,0.0),w);

			vs.clear();
			beta=0.0;
			betas.clear();
			betas.push_back(beta);
			alphas.clear();

			if (ipr)
			{
				outfile<<"iterations = "<<iterations<<endl;
				outfile<<"vs.size()  = "<<vs.size()<<endl;
			}
			outfile.flush();

			for (int it=0;it<iterations;it++)
			{
			       	if (ipr)
			       	{
					outfile<<"================================================================="<<endl;
					outfile<<"Doing Lanczos iteration,cycle,evec = "<<it<<" "<<cycle<<"  "<<num<<endl;
				}

			       	//vs.push_back(v_p);     
			       	
				/////////////////////////
				////////// Hv ///////////
				/////////////////////////
				
				double time_Hv = wall_clock(); Hv_calls++;
				sHv(&*v_p.begin(), &*w.begin(), pblock_states, Norm, Norm_inv, rep_loc, state_ph, lambda, adj_list, Ch, ch_list, ChCh, chch_list, n_spins, Ia, Ib, bits_right, T1, T2, T3, px, py, pz, outfile);
				outfile << "\n\t Wall clock for Hv step: " << (wall_clock() - time_Hv)/60.0 << endl;
				outfile.flush();
				time_Hv_total += wall_clock() - time_Hv;
			       	//if (it!=0) zaxpyk(size,-beta,v_o,w);
			       	if (it!=0) zaxpyk(size,-beta,vs[it-1],w);
				alpha=conj(zdotc(size, w, v_p));
			       	alphas.push_back(alpha);
			       	zaxpyk(size,-alpha,v_p,w);
			       	vs.push_back(v_p);     
			       	//v_o=v_p;
			       	beta=sqrt(real(zdotc(size, w, w)));
			       	#pragma omp parallel for 
				for (int64_t i=0;i<size;i++) v_p[i]=w[i];

			       	zscalk(size,1.0/beta,v_p);
			       	betas.push_back(beta);
			       	zscalk(size,0.0,w);

				//////////////////////////////////////////////////////////////////
			       	// Orthogonlize v_p with respect to all previous eigenvectors ////

				double time_ortho = wall_clock();
			       	orth_wrt_previous_evecs(outfile,v_p,previous_evecs);
				if (cycle==ncycles-1)
				{
					// Bad guess needs no orth
			       		orth_wrt_previous_evecs(outfile,v_p,vs);
			       	}
				normalize(v_p);

			       	if (ipr)
			       	{
					outfile <<"Time to perform Lanczos iteration "<<it<<" was "<<dif<<" seconds"<<endl;
					outfile <<"================================================================="<<endl;
			       	}

				if (ipr) {outfile <<"vs.size() = "<<vs.size()<<endl;}

				int tmp=alphas.size();
				t_mat.resize(tmp,tmp);
				t_eigenvecs.resize(tmp,tmp);
				for (int j=0;j<tmp*tmp;j++){t_mat[j]=0.0;t_eigenvecs[j]=0.0;}
				for (int j=0;j<tmp;j++)
				{
					t_mat(j,j)=alphas[j];
					if (j+1<tmp)
					{t_mat(j,j+1)=betas[j+1];t_mat(j+1,j)=betas[j+1];}
				}

//				if (ipr) {outfile <<"Time to build  T was "<<dif<<" seconds"<<endl;}
				
				std::vector<double> reigs(tmp);
				zMatrix_Diagonalize_Hermitian(t_mat, reigs, t_eigenvecs);
				eigs.clear(); for (int b=0;b<reigs.size();b++) eigs.push_back(reigs[b]);
				
//				if (ipr) {outfile<<"Time to diagonalize T was "<<dif<<" seconds"<<endl;}
				
				if (cycle!=0 and cycle!=ncycles-1 and it>1) {
					outfile <<"Eigs[0]      = "<<boost::format("%+.15f") %eigs[0]<<endl;
					outfile <<"Previous Eig = "<<boost::format("%+.15f") %previous_eig<<endl;
					if (abs(eigs[0]-previous_eig)<tol) {outfile <<"Convergence, doing last cycle with orth !!"<<endl;it=iterations;cycle=ncycles-2;}
					previous_eig=eigs[0];
				}
				else if (cycle==ncycles-1) {
					outfile <<"Final cycle, Eigs[0]      = "<<boost::format("%+.15f") %eigs[0]<<endl;
				}
				else{
					outfile <<"Eigs[0]      = "<<boost::format("%+.15f") %eigs[0]<<endl;
					outfile <<"Previous Eig = "<<boost::format("%+.15f") %previous_eig<<endl;
					previous_eig=eigs[0];
				}	

				outfile << "\n\t Wall clock for orthogonalize " << (wall_clock() - time_ortho)/60.0 << endl;
				time_ortho_total += wall_clock() - time_ortho;
				//////////////////////////////////////////////////////////////////////

			}

			if (ipr) {outfile <<"Making lowest Ritz eigenvector"<<endl;}
			#pragma omp parallel for
			for (int i=0;i<size;i++)
			{
				v_p[i]=0.0;
				for (int k=0;k<vs.size();k++)
				{v_p[i]+=vs[k][i]*t_eigenvecs(k,0);}
			}
		       orth_wrt_previous_evecs(outfile,v_p,previous_evecs);
		       normalize(v_p);
		}	
		previous_evecs.push_back(v_p);
		/*if (v_p.size()<30 and )
		{
			for (int i=0;i<v_p.size();i++)
			{
				
			}
		}*/
		/*// Use first excited state as start guess for next vector
		#pragma omp parallel for
		for (int i=0;i<size;i++)
		{
			v_p[i]=0.0;
			for (int k=0;k<vs.size();k++)
			{v_p[i]+=vs[k][i]*t_eigenvecs(k,1);}
		}*/
			
		outfile << "Eigenvalue number " << num << " = " << boost::format("%+.15f") %eigs[0] << endl;
		all_eigs.push_back(eigs[0]);
	}
	eigs=all_eigs;

	outfile << "\n\n--------------------------------------------------\n";
	outfile << "\n Lanczos SUMMARY \n" << endl;
	outfile << " Hv calls\t\t\t" << Hv_calls << endl;
	outfile << " Hv total time \t\t\t" << time_Hv_total/60.0 << " minutes"<< endl;
	outfile << " Orthogonalize total time \t" << time_ortho_total/60.0 << " minutes" << endl;
}
