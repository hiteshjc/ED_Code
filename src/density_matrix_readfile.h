#ifndef GUARD_EV_FROM_FILE_H
#define GUARD_EV_FROM_FILE_H

#include "global.h"
#include "idsClass.h"


void readCuts(	string &filein,
		vector<int> &Cut,
		string &fileEV1,
		string &fileEV2)
{	
	ifstream infile(filein.c_str());
	string line;
	string sCut("Cut=");
	string EV1("EV1_file=");
	string EV2("EV2_file=");

	while(getline(infile, line))
	{
		string tmp;
		char r;
		istringstream linestream(line);

		while(linestream >> r)
		{
			tmp.push_back(r);

			if(tmp.compare(sCut) == 0)
			{
				char s;
				while(linestream >> s)
				{
					if(s == '[' || s == ',')
					{
						int x; linestream >> x;
						Cut.push_back(x);
					}
				}
				tmp.clear();
			}


			if(tmp.compare(EV1) == 0)
			{	
				char s;
				while(linestream >> s)
					fileEV1.push_back(s); 
			}
			
			if(tmp.compare(EV2) == 0)
			{	
				char s;
				while(linestream >> s)
					fileEV2.push_back(s); 
			}
		}
	}



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
}


#endif
