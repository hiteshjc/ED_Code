#include "read_vector.h"

void read_pEV(	int64_t &n_p,
		double &px,
		double &py,
		vector<int64_t> &pblock_states,
		vector<double> &norm,
		vector< complex <double> > &v_coeff, 
		string &filepEV)
{
	ifstream infile(filepEV.c_str());
	string line;
	string p1p2("np 		 px 		py");	

	while(getline(infile, line))
	{
		if(line.compare(p1p2) == 0)
		{
			infile >> n_p >> px >> py;
			break;
		}
	}	
	infile.close();
		
	string phead("pblock_states 		 Norm 			 v_coeff");
	pblock_states.resize(n_p);
	norm.resize(n_p);
	v_coeff.resize(n_p);	

	int64_t jb = 0;
	
	ifstream infile1(filepEV.c_str());
	while(getline(infile1, line))
	{
		complex<double> tmp_coeff;
 
		if(line.compare(phead) == 0)
		{
			while(infile1 >> pblock_states[jb] >> norm[jb] >> tmp_coeff)
			{	
				v_coeff[jb] = tmp_coeff;
				jb++;			
			}
		}
	}
	infile1.close();
}

void readEV(	vector<ids> &psi,
		string &fileEV)
{
	ifstream infile(fileEV.c_str());
	string line;
	string eigenvector("EIGENVECTOR");

	while(getline(infile, line))
	{
		if(line.compare(eigenvector) == 0)
		{
			ids tmp;
			while(infile >> tmp.stateID >> tmp.coeff)
			{	psi.push_back(tmp);}
		}

	}

	infile.close();
}


