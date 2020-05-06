
#include "global.h"
#include "idsClass.h"
#include "mes_dmmat_readfile.h"
#include "Entropy.h"


void read_density_matrix(	zMatrix &dmat,
				string &file_dm,
				int &count)
{
	ifstream infile(file_dm.c_str());
	string line;
	string sdm("DENSITY_MATRIX");

	count = 0;
	while(getline(infile, line))
	{
		if(line.compare(sdm) == 0)
		{
			for (int64_t i = 0; i < dmat.NRows(); i++)
				for(int64_t j = 0; j < dmat.NCols(); j++)
				{
					complex<double> tmp;
					infile >> tmp;
					if (abs(real( tmp )) > 0.0 || abs(imag( tmp )) > 0.0 )
							count++;
					 dmat(i,j) = tmp;
				}
		}
	
	}

	infile.close();
}

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
	int n = 4;
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
	read_dmmes(filein, Cut, file_dm[0], file_dm[1], file_dm[2], file_dm[3], c_min, interval, delta);

	int n_tot = 1/delta;
	double c_max = c_min + interval;	

	ofstream outfile(fileout.c_str());

	outfile << " Cut = ";
	for(int64_t i = 0; i < Cut.size(); i++)
		outfile << Cut[i] << " ";

	outfile << endl;
	for(int i = 0; i < n; i++)
		outfile << file_dm[i] << endl;

//	complex <double> t00 = 0.0, t11 = 0.0, t22 = 0.0, t33 = 0.0;
//	complex <double> t01 = 0.0, t02 = 0.0, t03 = 0.0;
//	complex <double> t12 = 0.0, t13 = 0.0, t23 = 0.0;

	zMatrix traces(n,n,0.0);

	double ts_dm = wall_clock();
		
	for(int ones = 0; ones <= Cut.size(); ones++)
	{
		int n = 4;

		int64_t n_dm = n_choose_k(Cut.size(),ones);
		outfile << "\n Size of density matrix = " << n_dm << " x " << n_dm << endl;

//		int tmp_ones = ones + (int)'0';
		string sones;
		ostringstream convert; convert << ones; 
		sones = convert.str(); //= string(tmp_ones);

		vector<string> tmp_file_dm(4);
		for(int i = 0; i < n; i++)
			tmp_file_dm[i] = file_dm[i] + sones + ".txt";
	
		// Initialize and read in the density matrix from file
		vector < zMatrix > dm(4);

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

	//	zMatrix traces(n, n, 0.0);

		for(int i = 0 ; i < n; i++)
		{
			for(int j = 0; j < n; j++)
			{
				traces(i, j) += trace_matrix_prod(dm[i],dm[j]);
			}
		}

//		t00 += traces(0,0);	t11 += traces(1,1);	t22 += traces(2,2);	t33 += traces(3,3);
		
//		t01 += traces(0,1) + traces(1,0);	t02 += traces(0,2) + traces(2,0);	t03 += traces(0,3) + traces(3,0);
//		t12 += traces(1,2) + traces(2,1);	t13 += traces(1,3) + traces(3,1);	t23 += traces(2,3) + traces(3,2);

	}

	ts_dm = wall_clock() - ts_dm;
	outfile << "\n\n Total time to read in density matrices and compute traces  " << (ts_dm)/60.0 << " minutes\n\n\n";
	
	double ts_EE = wall_clock();
	int n_iter = interval/delta;
//	n_tot = 100;
	double ph_min = 0.0;

	Matrix S_Renyi(n_iter, n_tot + 1, 0.0);

	for(int i = 0; i < n_iter; i++)
	{
		for(int j = 0; j <= n_tot; j++)
		{
			double c = double(i)*delta+c_min;
			complex<double> ph = complex<double>(0.0, -ph_min + 2*M_PI*j/double(n_tot));
				
			complex<double> c1 = c, c2 = sqrt(1-c*c)*exp(ph);
/*
//			complex<double> S_tmp =   c1*c1*conj(c1)*conj(c1)*t00
						+ c1*c1*conj(c2)*conj(c2)*t11
						+ conj(c1)*conj(c1)*c2*c2*t22
						+ c2*c2*conj(c2)*conj(c2)*t33
						+ c1*conj(c1)*c1*conj(c2)*t01
						+ c1*conj(c1)*conj(c1)*c2*t02
						+ c1*conj(c1)*c2*conj(c2)*t03
						+ c1*conj(c1)*c2*conj(c2)*t12
						+ c1*conj(c2)*c2*conj(c2)*t13
						+ conj(c1)*c2*c2*conj(c2)*t23;
*/
			vector< complex<double> > coeff(4);
			coeff[0] = c1*conj(c1);	coeff[1] = c1*conj(c2);
			coeff[2] = c2*conj(c1); coeff[3] = c2*conj(c2);


			complex<double> S_tmp = complex<double> (0.0,0.0);

			for(int s1 = 0; s1 < n; s1++)
				for(int s2 = 0; s2 < n; s2++)
					S_tmp += coeff[s1]*coeff[s2]*traces(s1,s2);

			S_tmp = -log2(real(S_tmp));

			S_Renyi(i, j) = real(S_tmp);


		}
	}

	ts_EE = wall_clock() - ts_EE;

	outfile << "\n Total time for EE calc  " << (ts_EE)/60.0 << " minutes\n";


	outfile << " MIN VALUES " << endl;
	outfile << "c\t\tphi\t\t S_Renyi \n\n";

	for (int i = 1; i < (n_iter - 1); i++)
	{
		int ip = i+1, im = i-1;
		for (int j = 1; j < n_tot; j++)
		{
			int jm = j-1, jp = j+1;

			if(S_Renyi(i,j) <= S_Renyi(im,j)  && S_Renyi(i,j) <= S_Renyi(ip,j) )
				if(S_Renyi(i,j) <= S_Renyi(i,jm)  && S_Renyi(i,j) <= S_Renyi(i,jp) )
				{
					double c = double(i)*delta+c_min;
					complex<double> ph = complex<double>(0.0, ph_min + 2*M_PI*j/double(n_tot));
				
					outfile << fixed << setprecision(6) << c << "\t" << imag(ph) << "\t" <<  S_Renyi(i,j) << endl;

				}
		}
	}

	if( S_Renyi(0, 0) <= S_Renyi(1,0) )
	{	double c = 0.0, ph = 0.0;	outfile << fixed << setprecision(6) << c << "\t" << ph << "\t" <<  S_Renyi(0,0) << endl;}
	if( S_Renyi(n_iter-1, 0) <= S_Renyi(n_iter-2,0) )
	{	double c = 1.0, ph = 0.0; 	outfile << fixed << setprecision(6) << c << "\t" << ph << "\t" <<  S_Renyi(n_iter-1,0) << endl;}




	outfile << endl << endl;

  	outfile << fixed << setprecision(6) << "c1|psi1>" << "\t\t" << "c2|psi2>" << "\t\t" << "E.Entropy \t\t c \t\t ph \t\t E.Entropy" << endl;
	for(int i = 0; i < n_iter; i++)
	{
		for(int j = 0; j <= n_tot; j++)
		{
			double c = double(i)*delta+c_min;
			complex<double> ph = complex<double>(0.0, 2*M_PI*j/double(n_tot));
				
			complex<double> c1 = c, c2 = sqrt(1-c*c)*exp(ph);

			complex<double> S_tmp = complex<double> (S_Renyi(i,j), 0);

			S_Renyi(i, j) = real(S_tmp);
			outfile << fixed << setprecision(6) << c1 << "\t" << c2 << "\t" << S_tmp << "\t" << c << "\t" << imag(ph) << "\t" << real(S_tmp) << endl;
		}
	}

	



	outfile.close();

	return 0;
}
