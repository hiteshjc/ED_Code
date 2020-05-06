
#include "global.h"

void minE(	vector<double> &tmpE,
		vector<double> &tmpSz,
		double &tE,
		double &tSz)
{
	tE = tmpE[0]; tSz = tmpSz[0];

	for(int i = 1; i < tmpE.size(); i++)
	{
		if(tmpE[i] < tE)
		{
			tE = tmpE[i];
			tSz = tmpSz[i];
		}
	}

}


int main(int argc, char *argv[])
{
	////////////////////////// INPUT /////////////////////////////////////////
	string filein,fileout;

	if (argc <= 1)
	{
		cout << "Usage: " << argv[0] << " <Filename>" << endl;
		exit(1);
	}

	filein=argv[1];
	fileout=argv[2];

	ifstream infile(filein.c_str());
	
	string in; int n_spins = 0;
	string start_read("py)");
	string n_sites("sites/spins");

	vector<double> Energies, Spinz;

	while(infile >> in)
	{
		if(in.compare(n_sites) == 0)
		{
			char tmp; 
			infile >> tmp >> n_spins;
		}

		if(in.compare(start_read) == 0)
		{
//			cout << in << endl;
			
			int64_t id; double En, Sz; complex<double> k;
//			double nSz = 0.000;

			while(infile >> id >> En >> Sz >> k)
			{
//				cout << En << endl;
//				if(Sz == nSz || Sz == (nSz + 0.5) ) 
//				{	
					Energies.push_back(En);										
					Spinz.push_back(Sz);
//					nSz += 1.00;
//				}
			}		
		}
	}

	double tot_Sz = double(n_spins)/2.0;
	int tot_mag = 100;	
//	vector<double> minE, minSz;	
	ofstream outfile(fileout.c_str());
	
	for(int i = 0; i <= tot_mag; i++)
	{	
		int n = Energies.size();
		vector<double> tmpE(n);
		
		double h_mag = 1.25*double(i)/double(tot_mag);
	
		for(int j = 0; j < n; j++)
		{			
			tmpE[j] = Energies[j] - h_mag*Spinz[j];			
		}

		double tE, tSz;
		minE(tmpE, Spinz, tE, tSz);
		
//		cout << h_mag << " " << tE << "\t" << tSz << endl;
		outfile << fixed << setprecision(8) << tE << "\t" << tSz/tot_Sz << "\t" << h_mag << endl;
	}

	outfile.close();
	return 0;
}
