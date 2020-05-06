#ifndef GUARD_EV_FROM_FILE_H
#define GUARD_EV_FROM_FILE_H

#include "global.h"
#include "idsClass.h"

void read_dmmes(	string &filein,
			vector<int> &Cut,
			string &file_dm11,
			string &file_dm12,
			string &file_dm21,
			string &file_dm22,
			double &c_min,
			double &interval,
			double &delta)
{	
	ifstream infile(filein.c_str());
	string line;
	string sCut("Cut=");
	string sdm11("dm11_file=");	
	string sdm12("dm12_file=");
	string sdm21("dm21_file=");
	string sdm22("dm22_file=");
	string sc_min("c_min=");	
	string sinterval("interval=");	
	string sdelta("delta=");

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

			if(tmp.compare(sdm11) == 0)
			{	
				char s;
				while(linestream >> s)
					file_dm11.push_back(s); 
			}
			

			if(tmp.compare(sdm12) == 0)
			{	
				char s;
				while(linestream >> s)
					file_dm12.push_back(s); 
			}

			if(tmp.compare(sdm21) == 0)
			{	
				char s;
				while(linestream >> s)
					file_dm21.push_back(s); 
			}

			if(tmp.compare(sdm22) == 0)
			{	
				char s;
				while(linestream >> s)
					file_dm22.push_back(s); 
			}

			if(tmp.compare(sc_min) == 0)
			{	
				linestream >> c_min; 
			}

			if(tmp.compare(sinterval) == 0)
			{	
				linestream >> interval; 
			}

			if(tmp.compare(sdelta) == 0)
			{	
				linestream >> delta; 
			}
		}
	}



}


#endif
