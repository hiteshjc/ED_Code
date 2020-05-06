#include "dm_readfile.h"


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


void readCuts(	string &filein,
		vector<int> &Cut,
		string &fileEV1,
		string &fileEV2,
		int &n_ones)
{	
	ifstream infile(filein.c_str());
	string line;
	string sCut("Cut=");
	string EV1("EV1_file=");
	string EV2("EV2_file=");
	string sn_ones("n_ones=");

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

			if(tmp.compare(sn_ones) == 0)
			{	
				linestream >> n_ones;
			}
		}
	}



}


void readpCuts(	string &filein,
		vector<int> &Cut,
		vector<int> &T1,
		vector<int> &T2,
		string &filepEV1,
		string &filepEV2)
{	
	ifstream infile(filein.c_str());
	string line;
	string sCut("Cut=");
	string EV1("pEV1_file=");
	string EV2("pEV2_file=");
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
