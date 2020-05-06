#include "global.h"
#include "ED_2D_Lattice.h"

void ED_2D_Lattice_Geometry(	string input_filename, 
				int &n_sites,
				double &Jz,
				double &J1,
				double &J2,
				vector <coordinates> &adj_list,
				double &Ch,
				vector <triangles> &ch_list, 
				vector<int> &T1,
				vector<int> &T2,
				vector<int> &Cut_Pts,
				double &Sz,
				int &kx,
				int &ky,
				bool &All_Sectors,				
				double &Specific_Sz_Sector,
				bool &Specific_pSector,
				int &Num_Eigenvalues,
				bool &Find_Vectors,
				bool &Unwrap_pVector,
				string &Vector_File_Path,
				ofstream &outfile)
{
	ifstream infile(input_filename.c_str());
	stringstream infss;
	string line;
//	bool read = 0;
	int read_pairs = 0;
	int read_coor = 0;
//	double J1_tmp, J2_tmp;

// Input strings to be read from file
	string sJz("Jz=");
	string sJ1("J1=");
	string sPairs1("Pairs1=");
	string sJ2("J2=");
	string sPairs2("Pairs2=");

	string sCh("Ch=");
	string sTriangles("Triangles=");
	

	string sT1("T1=");
	string sT2("T2=");
	string sAll_Sectors("All_Sectors=");
	string sSpecific_Sz_Sector("Specific_Sz_Sector=");
	string sSpecific_pSector("Specific_pSector=");
	string sSz("Sz=");
	string skx("kx=");
	string sky("ky=");
	
	string sFind_Vectors("Find_Vectors=");	
	string sUnwrap_pVector("Unwrap_pVector=");
	string sNum_Eigenvalues("Num_Eigenvalues=");
	string sVector_File_Path("Vector_File_Path=");
	
	while(getline(infile, line))
	{
		istringstream linestream(line);
		char r;
		string read;
		
		while(linestream >> r)
		{
			if(r == '[')
				break;

			read.push_back(r);

			if(read.compare(sJz) == 0)
			{
				linestream >> Jz;
				read.clear();
			}
			
		// Reading in J1 and Pairs1 and storing them in the adjacency list	
			if(read.compare(sJ1) == 0)
			{	linestream >> J1;}
			
			if(read.compare(sPairs1) == 0)
			{
				char r_tmp;
				while(linestream >> r_tmp)
				{
					if(r_tmp == '[' || r_tmp == ',')
					{
						char op_b, comma, cl_b;
						int x, y;
						linestream >> op_b >> x >> comma >> y >> cl_b;

						coordinates alist; 
						alist.x = x; alist.y = y; alist.J = J1;
						adj_list.push_back(alist);	    
					}
					else if (r_tmp == ']')
					break;
				}
				read.clear();
			}		
	
		// Reading in J1 and Pairs1 and storing them in the adjacency list	
		
			if(read.compare(sJ2) == 0)
			{	linestream >> J2;}
			
			if(read.compare(sPairs2) == 0)
			{
				char r_tmp;
				while(linestream >> r_tmp)
				{
					if(r_tmp == '[' || r_tmp == ',')
					{
						char op_b, comma, cl_b;
						int x, y;
						linestream >> op_b >> x >> comma >> y >> cl_b;

						coordinates alist; 
						alist.x = x; alist.y = y; alist.J = J2;
						adj_list.push_back(alist);	    
					}
					else if (r_tmp == ']')
					break;
				}
				read.clear();
			}

		// Reading in ch and Triangles
			if(read.compare(sCh) == 0)
			{	linestream >> Ch;}
			
			if(read.compare(sTriangles) == 0)
			{
				char r_tmp;
				while(linestream >> r_tmp)
				{
					if(r_tmp == '[' || r_tmp == ',')
					{
						char op_b, comma, cl_b;
						int x, y, z;
						linestream >> op_b >> x >> comma >> y >> comma >> z >> cl_b;
						
				// Store all three possible combinations from each triangle
				// This way i and j represent the bond and k the Sz part
						triangles t1, t2, t3;
						t1.i = x; t1.j = y; t1.k = z;
						t2.i = y; t2.j = z; t2.k = x;
						t3.i = z; t3.j = x; t3.k = y;

						ch_list.push_back(t1);
						ch_list.push_back(t2);
						ch_list.push_back(t3);	    
					}
					else if (r_tmp == ']')
					break;
				}
				read.clear();
			}
			
			if(read.compare(sT1) == 0)
			{
				char s;
				
				while(linestream >> s)
				{
					if(s == '[' || s == ',')
					{
						int x; linestream >> x;	  
						T1.push_back(x);
					}	
					else if(s == ']')
					break;
				}      
				read.clear();
			}    


			if(read.compare(sT2) == 0)
			{
				char s;
			
				while(linestream >> s)
				{
					if(s == '[' || s == ',')
					{
						int x; linestream >> x;	  
						T2.push_back(x);
					}	
					else if(s == ']')
					break;
				}      
				read.clear();
			}    
			
			if(read.compare(sAll_Sectors) == 0)
			{
				linestream >> All_Sectors;
				read.clear();
			}
			
			if(read.compare(sSpecific_Sz_Sector) == 0)
			{
				linestream >> Specific_Sz_Sector;
				read.clear();
			}


			if(read.compare(sSpecific_pSector) == 0)
			{
				linestream >> Specific_pSector;
				read.clear();
			}

			if(read.compare(sSz) == 0)
			{
				linestream >> Sz;
				read.clear();
			}

			if(read.compare(skx) == 0)
			{
				linestream >> kx;
				read.clear();
			}
		
			if(read.compare(sky) == 0)
			{
				linestream >> ky;
				read.clear();
			}

			if(read.compare(sNum_Eigenvalues) == 0)
			{
				linestream >> Num_Eigenvalues;
				read.clear();
			}

			if(read.compare(sFind_Vectors) == 0)
			{
				linestream >> Find_Vectors;
				read.clear();
			}

			if(read.compare(sUnwrap_pVector) == 0)
			{
				linestream >> Unwrap_pVector;
				read.clear();
			}

			if(read.compare(sVector_File_Path) == 0)
			{
				char r_tmp; Vector_File_Path.clear();
				while(linestream >> r_tmp)
				{	Vector_File_Path.push_back(r_tmp);}
				read.clear();
			}
		}

	}


	if(T1.size() != T2.size())
	{
		outfile << "\nERROR: input T1 != T2\n";
		exit(1);
	}
	else
		n_sites = T1.size();
}	

