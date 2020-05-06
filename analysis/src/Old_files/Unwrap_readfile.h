#ifndef GUARD_EV_FROM_FILE_H
#define GUARD_EV_FROM_FILE_H

#include "global.h"
#include "idsClass.h"


void readunwrap(	string &filein,
			vector<int> &T1,
			vector<int> &T2,
			string &filepEV)
{	
	ifstream infile(filein.c_str());
	string line;
	string t1("T1=");
	string t2("T2=");
	string pEV("pEV_file_path=");

	while(getline(infile, line))
	{
		string tmp;
		char r;
		istringstream linestream(line);

		while(linestream >> r)
		{
			tmp.push_back(r);

			if(tmp.compare(t1) == 0)
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
					
			if(tmp.compare(t2) == 0)
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

			if(tmp.compare(pEV) == 0)
			{	
				char s;
				while(linestream >> s)
					filepEV.push_back(s); 
			}
			
		}
	}



}


#endif
