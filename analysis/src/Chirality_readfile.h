#ifndef GUARD_EV_FROM_FILE_H
#define GUARD_EV_FROM_FILE_H

#include "global.h"
#include "idsClass.h"

void readTriangles(	string &filein,
			vector<triangles> &tri,
			string &fileEV1,
			string &fileEV2)
{	
	ifstream infile(filein.c_str());
	string line;
	string sTriangles("Triangles=");
	string sEV1("EV1_File_Path=");
	string sEV2("EV2_File_Path=");

	while(getline(infile, line))
	{
		string tmp;
		char r;
		istringstream linestream(line);

		while(linestream >> r)
		{
			tmp.push_back(r);

			if(tmp.compare(sTriangles) == 0)
			{
				char s;
				while(linestream >> s)
				{
					if(s == '[' || s == ',') 
					{
						char ob, comma, cb;
						triangles tri_tmp; 
						linestream >> ob >> tri_tmp.i >> comma >> tri_tmp.j >> comma >> tri_tmp.k >> cb;
						tri.push_back(tri_tmp);
					}
				}
				tmp.clear();
			}
					
			if(tmp.compare(sEV1) == 0)
			{	
				char s;
				while(linestream >> s)
					fileEV1.push_back(s); 
			}
			
			if(tmp.compare(sEV2) == 0)
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
