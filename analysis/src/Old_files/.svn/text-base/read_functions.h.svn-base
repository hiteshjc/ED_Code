#ifndef GUARD_READ_FUNCTIONS
#define GUARD_READ_FUNCTIONS

#include "global.h"

template <class T>
void read_Parameter(	string filein,
			string sParameter,
			T &Parameter)

{
	ifstream infile(filein.c_str() );
	string line;
	bool scan = true;

	cout << sParameter << endl;

	while( getline(infile, line) && scan )
	{
		string tmp;
		char r;
		istringstream linestream(line);

		while(linestream >> r)
		{
			tmp.push_back(r);
			if(tmp.compare(sParameter) == 0)
			{
				linestream >> Parameter; 
				scan = false;
			}
		}	
	}

	infile.close();
}

template <class T>
void read_list_of_sites(	string &filein,
				string slist,
				vector<T> &list)
{	
	ifstream infile(filein.c_str());
	string line;
	bool scan = true;

	while( getline(infile, line) && scan )
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
						T x; linestream >> x;
						list.push_back(x);
					}
				}
				tmp.clear();
				scan = false;
			}
		}
	}
	infile.close();
}

void read_file_string(	string &filein,
			string &key_string,
			string &file_string)
{	
	ifstream infile(filein.c_str());
	string line;
	bool scan = true;

	while( getline(infile, line) && scan )
	{
		string tmp;
		char r;
		istringstream linestream(line);

		while(linestream >> r)
		{
			tmp.push_back(r);
			if(tmp.compare(key_string) == 0)
			{	
				char s;
				while(linestream >> s)
					file_string.push_back(s); 
				scan = false;
			}
		}
	}
	infile.close();
}

#endif
