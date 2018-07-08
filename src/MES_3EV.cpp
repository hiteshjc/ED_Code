
#include "global.h"
#include "idsClass.h"
#include "MES_readfile.h"
#include "Entropy.h"

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
	string filein, fileEV1, fileEV2, fileEV3, fileEV4, fileEV5, fileEV6, fileEV7, fileEV8, fileout;

	if (argc <= 1)
	{
		cout << "Usage: " << argv[0] << " <Filename>" << endl;
		exit(1);
	}

	filein=argv[1];
	fileout=argv[2];

	vector<int> Cut1, Cut2;
	readCuts(filein, Cut1, Cut2, fileEV1, fileEV2, fileEV3);

	vector<ids> psi[3];
	readEV(psi[0], fileEV1);
	readEV(psi[1], fileEV2);
	readEV(psi[2], fileEV3);

	ofstream outfile(fileout.c_str());

	for(int i = 0; i < 3; i++)
	{
		outfile << " Size Vector psi = " << psi[i].size() << endl;
	}
	int64_t n_dmI = int64_t_power(2, Cut1.size()), n_dmII = int64_t_power(2, Cut2.size());

	double ts_dm = wall_clock();
	int64_t n_dm = n_dmI;
	vector<int> Cut = Cut1;
	zMatrix	dm11(n_dm, n_dm, 0.0);
	zMatrix dm12(n_dm, n_dm, 0.0);
	zMatrix	dm13(n_dm, n_dm, 0.0);

	zMatrix dm21(n_dm, n_dm, 0.0);
	zMatrix	dm22(n_dm, n_dm, 0.0);
	zMatrix dm23(n_dm, n_dm, 0.0);

	zMatrix	dm31(n_dm, n_dm, 0.0);
	zMatrix dm32(n_dm, n_dm, 0.0);
	zMatrix	dm33(n_dm, n_dm, 0.0);

	density_matrix(dm11, psi[0], Cut);	
	density_matrixij(dm12, psi[0], psi[1], Cut);
	density_matrixij(dm13, psi[0], psi[2], Cut);

	density_matrix(dm22, psi[1], Cut);	
	density_matrixij(dm21, psi[1], psi[0], Cut);
	density_matrixij(dm23, psi[1], psi[2], Cut);

	density_matrixij(dm31, psi[2], psi[0], Cut);
	density_matrixij(dm32, psi[2], psi[1], Cut);
	density_matrix(dm33, psi[2], Cut);
	
	ts_dm = wall_clock() - ts_dm;
	outfile << "\n Cut \t";
	for(int i = 0; i < Cut.size(); i++)
		outfile << Cut[i] << " ";
	outfile << "\n Total time to construct density matrices  " << (ts_dm)/60.0 << " minutes\n";

	int n_tot = 100;	

	double ts_EE = wall_clock();
	outfile << fixed << setprecision(6) << "c1|psi1>" << "\t\t" << "c2|psi2>" << "\t\t" << "E.Entropy" << endl;

	for(int i = 0; i <= n_tot; i++)
	{
		for(int j = 0; j <= n_tot; j++)
		{
			double c = double(i)/double(n_tot);
			complex<double> ph = complex<double>(0.0, 2*M_PI*j/double(n_tot));
				
			complex<double> c1 = c, c2 = sqrt(1-c*c)*exp(ph);
			zMatrix dmsum(n_dm, n_dm, 0.0);
			add_density_matrices(c1, c2, dm11, dm12, dm21, dm22, dmsum);
		
			complex<double> S_tmp = Von_Neumann_Entropy(dmsum);
			outfile << fixed << setprecision(6) << c1 << "\t" << c2 << "\t" << S_tmp << endl;
		}
	}
	ts_EE = wall_clock() - ts_EE;

	outfile << "\n Total time for EE calc  " << (ts_EE)/60.0 << " minutes\n";
	outfile.close();

	return 0;
}
