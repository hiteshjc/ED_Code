#ifndef GUARD_EV_FROM_FILE_H
#define GUARD_EV_FROM_FILE_H

#include "global.h"
#include "idsClass.h"


void read_mes(	string &filein,
		vector<int> &Cut,
		string &fileEV1,
		string &fileEV2,
		int &alpha,
		double &c_min,
		double &interval,
		double &delta)
{	
	ifstream infile(filein.c_str());
	string line;
	string sCut("Cut=");
	string EV1("EV1_file=");
	string EV2("EV2_file=");
	string sc_min("c_min=");	
	string sinterval("interval=");	
	string sdelta("delta=");
	string salpha("alpha=");

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

			if(tmp.compare(salpha) == 0)
			{	
				linestream >> alpha; 
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
void read_mes(	string &filein,
		vector<int> &Cut,
		string &filepEV1,
		string &filepEV2,
		vector<int> &T1,
		vector<int> &T2,
		int &alpha,
		double &c_min,
		double &interval,
		double &delta)
{	
	ifstream infile(filein.c_str());
	string line;
	string sCut("Cut=");
	string EV1("pEV1_file=");
	string EV2("pEV2_file=");
	string sc_min("c_min=");	
	string sinterval("interval=");	
	string sdelta("delta=");
	string salpha("alpha=");
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

			if(tmp.compare(salpha) == 0)
			{	
				linestream >> alpha; 
			}

			if(tmp.compare(EV1) == 0)
			{	
				char s;
				while(linestream >> s)
					filepEV1.push_back(s); 
			}
			
			if(tmp.compare(EV2) == 0)
			{	
				char s;
				while(linestream >> s)
					filepEV2.push_back(s); 
			}
		}
	}



}

#endif
