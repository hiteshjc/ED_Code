#ifndef GUARD_EV_FROM_FILE_H
#define GUARD_EV_FROM_FILE_H

#include "global.h"
#include "idsClass.h"

void readCorr_List(	string &filein,
			bool &SS,
			vector<int> &list,
			bool &BB,
			vector<links> &bonds,
			vector<int> &T1,
			vector<int> &T2,
			string &filepEV1)
{	
	ifstream infile(filein.c_str());
	string line;
	string sSS("SS=");
	string slist("list=");
	string sBB("BB=");
	string sbonds("bonds=");
	string sEV1("pEV1_File_Path=");
	string sT1("T1=");
	string sT2("T2=");

	while(getline(infile, line))
	{
		string tmp;
		char r;
		istringstream linestream(line);

		while(linestream >> r)
		{
			tmp.push_back(r);

			if(tmp.compare(slist) == 0)
			{
				char s;
				while(linestream >> s)
				{
					if(s == '[' || s == ',') 
					{
						int tmp;
						linestream >> tmp;
						list.push_back(tmp);
					}
				}
				tmp.clear();
			}

			if(tmp.compare(sT1) == 0)
			{
				char s;
				while(linestream >> s)
				{
					if(s == '[' || s == ',')
					{
						int x; linestream >> x;
						T1.push_back(x);
					}
				}
				tmp.clear();
			}
					
			if(tmp.compare(sT2) == 0)
			{
				char s;
				while(linestream >> s)
				{
					if(s == '[' || s == ',')
					{
						int x; linestream >> x;
						T2.push_back(x);
					}
				}
				tmp.clear();
			}

			if(tmp.compare(sbonds) == 0)
			{
				char r_tmp;
				while(linestream >> r_tmp)
				{
					if(r_tmp == '[' || r_tmp == ',')
					{
						char op_b, comma, cl_b;
						int x, y;
						linestream >> op_b >> x >> comma >> y >> cl_b;

						links tmp_bond;  
						tmp_bond.x = x; tmp_bond.y = y; 
						bonds.push_back(tmp_bond);	    
					}
					else if (r_tmp == ']')
					break;
				}
				tmp.clear();
			}							

			if(tmp.compare(sSS) == 0)
			{
				linestream >> SS;
				tmp.clear();
			}

			if(tmp.compare(sBB) == 0)
			{
				linestream >> BB;
				tmp.clear();
			}

			if(tmp.compare(sEV1) == 0)
			{	
				char s;
				while(linestream >> s)
					filepEV1.push_back(s); 
				tmp.clear();
			}
			
		}
	}


	infile.close();
}

/*
void read_pEV(	vector<ids> &psi,
		int64_t &n_p,
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

*/
#endif
