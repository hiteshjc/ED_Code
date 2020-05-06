
#include "global.h"
#include "idsClass.h"
#include "MES_mat_readfile.h"
#include "Entropy.h"


void read_density_matrix(	zMatrix &dmat,
				string &file_dm,
				int &count)
{
	ifstream infile(file_dm.c_str());
	string line;
	string sdm("DENSITY_MATRIX");

	count = 0;
//	int i = 0, j = 0;
//	bool start = 0;
//	check = 1;

	while(getline(infile, line))
	{
/*		if(line.compare(sdm) == 0)
		{
			getline(infile,line);		
			
			string lines;
			while(getline(infile, lines) )
			{
				j = 0;
				istringstream linestream(lines);
				complex<double> tmp;
	
				while(linestream >> tmp)
				{		
					dmat(i, j++) = tmp;
					if (abs(real( tmp )) > 0.0 || abs(imag( tmp )) > 0.0 )
							count++;
				}

				if (j != dmat.NCols() )
					check = 0;
				i++;
			}
		}
	}

	if (i != dmat.NRows() )
		check = 0;
*/	
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


void add_density_matrices(	complex<double> &c1,
				complex<double> &c2,
				zMatrix &dm11,
				zMatrix &dm12,
				zMatrix &dm21,
				zMatrix &dm22,
				zMatrix &dmsum)
{
	for(int64_t i = 0; i < dm11.NRows(); i++)
		for(int64_t j = 0; j < dm11.NRows(); j++)
			dmsum(i,j) = conj(c1)*c1*dm11(i,j) + conj(c1)*c2*dm12(i,j) + conj(c2)*c1*dm21(i,j) + conj(c2)*c2*dm22(i,j);	
}


int main(int argc, char *argv[])
{
	string filein, file_dm11, file_dm12, file_dm21, file_dm22, fileout;

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
	readCuts(filein, Cut, file_dm11, file_dm12, file_dm21, file_dm22, c_min, interval, delta);

	int n_tot = 1/delta;
	double c_max = c_min + interval;	

	ofstream outfile(fileout.c_str());

	outfile << " Cut = ";
	for(int64_t i = 0; i < Cut.size(); i++)
		outfile << Cut[i] << " ";
	outfile << endl;
	outfile << file_dm11 << endl;
	outfile << file_dm12 << endl;
	outfile << file_dm21 << endl;
	outfile << file_dm22 << endl;

	// Sizes of reduced density matrices along two different cuts
	int64_t n_dm = int64_t_power(2, Cut.size());
	
	outfile << "\n Size of density matrix = " << n_dm << " x " << n_dm << endl;

	// Initialize and read in the density matrix from file
	double ts_dm = wall_clock();
	zMatrix	dm11(n_dm, n_dm, 0.0);
	zMatrix dm12(n_dm, n_dm, 0.0);
	zMatrix	dm21(n_dm, n_dm, 0.0);
	zMatrix dm22(n_dm, n_dm, 0.0);

	int count = 0;

	read_density_matrix(dm11, file_dm11, count);	
	outfile << "\nNum of non-zero elements read for dm11 = " << count;
	outfile << "\n Entropy dm11 = " << Renyi_Entropy(dm11); //Von_Neumann_Entropy(dm11);
	
	read_density_matrix(dm12, file_dm12, count);
	outfile << "\nNum of non-zero elements read for dm12 = " << count;
	outfile << "\n Entropy dm12 = " << Renyi_Entropy(dm12); // Von_Neumann_Entropy(dm12);

	read_density_matrix(dm21, file_dm21, count);	
	outfile << "\nNum of non-zero elements read for dm21 = " << count;	
	outfile << "\n Entropy dm21 = " << Renyi_Entropy(dm21); // Von_Neumann_Entropy(dm21);

	read_density_matrix(dm22, file_dm22, count);
	outfile << "\nNum of non-zero elements read for dm22 = " << count;	
	outfile << "\n Entropy dm22 = " << Renyi_Entropy(dm22); // Von_Neumann_Entropy(dm22);

	ts_dm = wall_clock() - ts_dm;
	outfile << "\n\n Total time to read in density matrices  " << (ts_dm)/60.0 << " minutes\n\n\n";
	/////////////////////////////////////////////////////////////////////

	// Look for the minimally entangled linear combination of eigenstates
	n_tot = 50;	

	double ts_EE = wall_clock();
	outfile << fixed << setprecision(6) << "c1|psi1>" << "\t\t" << "c2|psi2>" << "\t\t" << "E.Entropy \t\t c \t\t ph \t\t E.Entropy" << endl;

	double t_tot_add = 0.0;
	double t_tot_diag = 0.0;	

	int n_iter = interval/delta;

	for(int i = 0; i < n_iter; i++)
	{
		for(int j = 0; j <= n_tot; j++)
		{
			double c = double(i)*delta+c_min;
			complex<double> ph = complex<double>(0.0, 2*M_PI*j/double(n_tot));
				
			complex<double> c1 = c, c2 = sqrt(1-c*c)*exp(ph);
			zMatrix dmsum(n_dm, n_dm, 0.0);
			// Adding all the four density matrices together with the co-efficients

			double t_add = wall_clock();
			add_density_matrices(c1, c2, dm11, dm12, dm21, dm22, dmsum);
			t_add = wall_clock() - t_add;
			t_tot_add += t_add;

			double t_diag = wall_clock();
			complex<double> S_tmp = Renyi_Entropy(dmsum); //Von_Neumann_Entropy(dmsum);
			t_diag = wall_clock() - t_diag;
			t_tot_diag += t_diag;

			outfile << fixed << setprecision(6) << c1 << "\t" << c2 << "\t" << S_tmp << "\t" << c << "\t" << imag(ph) << "\t" << real(S_tmp) << endl;
		}
	}

	ts_EE = wall_clock() - ts_EE;

	outfile << "\n Total time to add density matrices " << (t_tot_add)/60.0 << " minutes \n"; 		
	outfile << "\n Total time to diagonalize density matrix " << (t_tot_diag)/60.0 << " minutes \n"; 		
	outfile << "\n Total time for EE calc  " << (ts_EE)/60.0 << " minutes\n";
	outfile.close();

	return 0;
}
