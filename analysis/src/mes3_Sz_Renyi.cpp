
#include "global.h"
#include "idsClass.h"
#include "read_inputfile.h"
#include "Entropy.h"

complex<double> trace_matrix_prod(	zMatrix &m1,
					zMatrix &m2)
{
	complex<double> sum = complex<double>(0.0,0.0);
	double sr = 0.0;
	double si = 0.0;

	#pragma omp parallel for default(shared) reduction(+:sr, si)
	for(int64_t i = 0; i < m1.NRows(); i++)
		for(int64_t j = 0; j < m1.NRows(); j++)
		{
			complex<double> tmp = m1(i,j)*m2(j,i);
			sr += real(tmp);
			si += imag(tmp);
		}
	
	sum = complex<double> (sr, si);
	return sum;
//	return -log2(real(sum));

}

int main(int argc, char *argv[])
{
	string filein, fileout;
	int n = 9;
	vector <string> file_dm(n);


	if (argc <= 1)
	{
		cout << "Usage: " << argv[0] << " <Filename>" << endl;
		exit(1);
	}

	filein=argv[1];
	fileout=argv[2];

	double c_min = 0.0;	// default values
	double delta = 0.01;
	double interval = 1.0;
	// Read two Cuts along non-trivial directions from file
	vector<int> Cut;
	read_dmmes3(filein, Cut, file_dm[0], file_dm[1], file_dm[2], file_dm[3], file_dm[4], file_dm[5], file_dm[6], file_dm[7], file_dm[8], c_min, interval, delta);
	
	int n_tot = 1/delta;
	double c_max = c_min + interval;	

	ofstream outfile(fileout.c_str());

	outfile << " Cut = ";
	for(int64_t i = 0; i < Cut.size(); i++)
		outfile << Cut[i] << " ";

	outfile << endl;
	for(int i = 0; i < n; i++)
		outfile << file_dm[i] << endl;

	zMatrix t(n,n,0.0);

	double ts_dm = wall_clock();
		
	for(int ones = 0; ones <= Cut.size(); ones++)
	{
		int64_t n_dm = n_choose_k(Cut.size(),ones);
		outfile << "\n Size of density matrix = " << n_dm << " x " << n_dm << endl;

//		int tmp_ones = ones + (int)'0';
		string sones;
		ostringstream convert; convert << ones; 
		sones = convert.str(); //= string(tmp_ones);

		vector<string> tmp_file_dm(n);
		for(int i = 0; i < n; i++)
			tmp_file_dm[i] = file_dm[i] + sones + ".txt";


		// Initialize and read in the density matrix from file
		vector < zMatrix > dm(n);

		outfile << "Num of ones " << ones << endl;

		for(int i = 0; i < n; i++)
		{
			dm[i].resize(n_dm, n_dm);	

			int count = 0;
			read_density_matrix(dm[i], tmp_file_dm[i], count);	
			outfile << "Num of non-zero elements read for dm = " << count << endl;
		}
		outfile << endl;
		/////////////////////////////////////////////////////////////////////

//		cout << outfile << tmp_file_dm11 << endl;

		for(int i = 0 ; i < n; i++)
		{
			for(int j = 0; j < n; j++)
			{
				t(i, j) += trace_matrix_prod(dm[i],dm[j]);
			}
		}
	}

	ts_dm = wall_clock() - ts_dm;
	outfile << "\n\n Total time to read in density matrices and compute traces  " << (ts_dm)/60.0 << " minutes\n\n\n";
	
	double ts_EE = wall_clock();
	int n_iter = interval/delta;
	n_tot = 100;
	double ph_min = 0.0;

	
	Matrix4 srenyi(n_iter, n_iter, n_tot + 1, n_tot + 1, 0.0);

	double S_min = 1000;

	for(int i = 0; i < n_iter; i++)
	{

		double d1 = double(i)*delta+c_min;

		double c2_max = sqrt(1-d1*d1);

		for(int j = 0; j < n_iter ; j++)
		{
			double d2 = double(j)*delta+c_min;

			if(d2 <= c2_max)
			for(int k = 0; k <= n_tot; k++)
			{
				complex<double> ph2 = complex<double>(0.0, 2*M_PI*k/double(n_tot));
					
				for(int l=0; l <= n_tot; l++)
				{
					complex<double> c1 = d1;
					complex<double> c2 = d2*exp(ph2);
					complex<double> ph3 = complex<double>(0.0, 2*M_PI*l/double(n_tot));
					complex<double> c3 = sqrt(1-d1*d1-d2*d2)*exp(ph3);
			
					complex<double> c1s = conj(c1), c2s = conj(c2), c3s = conj(c3);

					vector < complex<double> > coeff(9);
					coeff[0] = c1*c1s; coeff[1] = c1*c2s; coeff[2] = c1*c3s;
					coeff[3] = c2*c1s; coeff[4] = c2*c2s; coeff[5] = c2*c3s;
					coeff[6] = c3*c1s; coeff[7] = c3*c2s; coeff[8] = c3*c3s;

					complex<double> S_tmp = complex<double> (0.0, 0.0);
					
					for(int s1 = 0; s1 < n; s1++)
						for(int s2 = 0; s2 < n; s2++)
							S_tmp += coeff[s1]*coeff[s2]*t(s1,s2);

		
					S_tmp = -log2(real(S_tmp));
		
					if(real(S_tmp) < S_min)
					{
						S_min = real(S_tmp);
//						outfile << fixed << setprecision(6) << c1 << "\t" << c2 << "\t" << c3 << "\t" << S_tmp << "\t" << d1 << "\t" << d2 << "\t" << imag(ph2) << "\t" << imag(ph3) << "\t" << S_tmp << endl;
					}	

					if (c1 == 1.0 && c2 == 0.0 && c3 == 0.0 && k == 0 && l == 0)
						outfile << fixed << setprecision(6) << "c1 =\t" << c1 << "\t" << S_tmp << endl; 

					if (c1 == 0.0 && c2 == 1.0 && c3 == 0.0 && k == 0 && l == 0)
						outfile << fixed << setprecision(6) << "c2 =\t" << c2 << "\t" << S_tmp << endl; 

					if (c1 == 0.0 && c2 == 0.0 && c3 == 1.0 && k == 0 && l == 0)
						outfile << fixed << setprecision(6) << "c3 =\t" << c3 << "\t" << S_tmp << endl; 


					srenyi(i,j,k,l) = real(S_tmp);
				}
			}

		}
	}

	outfile << endl << endl;
	outfile << "Minimum values of entropy = " << S_min << endl;
	outfile << " MIN VALUES " << endl;
	outfile << "\tc1\t\t\tc2\t\tphi2\t\t\tc3\t\t\tph3\t\t\t S_Renyi \n\n";
	outfile << endl;

	for (int i = 1; i < (n_iter-1); i++)
	{
		double d1 = double(i)*delta+c_min;
		double c2_max = sqrt(1-(d1+delta)*(d1+delta) );


		for (int j = 1; j < (n_iter-1); j++)
		{
			double d2 = double(j)*delta+c_min;
	
			if(d2 <= c2_max)
			for (int k = 1; k < n_tot; k++)
			{
				for (int l = 1; l < n_tot; l++)
				{
					if(srenyi(i,j,k,l) <= srenyi(i+1,j,k,l)  && srenyi(i,j,k,l) <= srenyi(i-1,j,k,l) )
					if(srenyi(i,j,k,l) <= srenyi(i,j+1,k,l)  && srenyi(i,j,k,l) <= srenyi(i,j-1,k,l) )
					if(srenyi(i,j,k,l) <= srenyi(i,j,k+1,l)  && srenyi(i,j,k,l) <= srenyi(i,j,k-1,l) )
					if(srenyi(i,j,k,l) <= srenyi(i,j,k,l+1)  && srenyi(i,j,k,l) <= srenyi(i,j,k,l-1) )
					{
							complex<double> c1 = d1;
							complex<double> ph2 = complex<double>(0.0, 2*M_PI*k/double(n_tot));
							complex<double> c2 = d2*exp(ph2);
							complex<double> ph3 = complex<double>(0.0, 2*M_PI*l/double(n_tot));
							complex<double> c3 = sqrt(1-d1*d1-d2*d2)*exp(ph3);

							outfile << fixed << setprecision(6) << c1 << "\t" << c2 << "\t" << ph2 << "\t" << c3 << "\t" << ph3 << "\t" <<  srenyi(i,j,k,l) << endl;
					}
				}
			}
		}
	}

	outfile << endl;
//	if( srenyi(0, 0, 0, 0) <= S_Renyi(1,0,0,0) )
//	{	double c1 = 0.0, ph = 0.0;	outfile << fixed << setprecision(6) << c << "\t" << c << "\t" << ph << "\t" << c << "\t" << ph << "\t" <<  srenyi(0,0,0,0) << endl;}
	if( srenyi(n_iter-1,0,0,0) < srenyi(n_iter-2,0,0,0) )
	{	
		complex<double> c1 = complex<double>(1.0,0.0) , c2 = complex<double>(0.0,0.0), ph = 0.0; 	
		outfile << fixed << setprecision(6) << c1 << "\t" << c2 << "\t" << ph << "\t" << c2 << "\t" << ph << "\t" <<  srenyi(n_iter-1,0,0,0) << endl;
	}
	if( srenyi(0,n_iter-1,0,0) < srenyi(0,n_iter-2,0,0) )
	{	
		complex<double> c1 = complex<double>(1.0,0.0) , c2 = complex<double>(0.0,0.0), ph = 0.0; 	
		outfile << fixed << setprecision(6) << c2 << "\t" << c1 << "\t" << ph << "\t" << c2 << "\t" << ph << "\t" <<  srenyi(0,n_iter-1,0,0) << endl;
	}
	if( srenyi(0,0,0,0) < srenyi(1,0,0,0) && srenyi(0,0,0,0) < srenyi(0,1,0,0) )
	{
		complex<double> c1 = complex<double>(1.0,0.0) , c2 = complex<double>(0.0,0.0), ph = 0.0; 	
		outfile << fixed << setprecision(6) << c2 << "\t" << c2 << "\t" << ph << "\t" << c1 << "\t" << ph << "\t" <<  srenyi(0,0,0,0) << endl;
	}






	outfile.close();

	return 0;
}
