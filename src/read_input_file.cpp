#include "read_input_file.h"

complex<double> get_expphase(triangles t, vector <coordinates> &adj_list)
{
	for (int n=0;n<adj_list.size();n++)
	{
// RECALL: What enters chirality formula:  i/2 * ( exp (i theta_ij) S^i + S^j - exp (-theta_ij) S^i- Sj^+ ) Sk^z
// In Krishna's notation x-y (his adjacency list) = bond
// In Krishna's notation Chirality single term (his triangle) is i/2 * ( exp (i theta_ij) S^i + S^j - exp (-theta_ij) S^i- Sj^+ ) Sk^z
		if (t.i==adj_list[n].x and t.j==adj_list[n].y) return adj_list[n].expphase;
		if (t.i==adj_list[n].y and t.j==adj_list[n].x) return conj(adj_list[n].expphase);
	}
	return complex<double>(1.0,0.0);
}

bool neighbor_tri(triangles tri1, triangles tri2, vector <coordinates> &adj_list)
{
	// If triangles are nearest neighbors 
	for (int n=0;n<adj_list.size();n++)
        {
                int x=adj_list[n].x;
                int y=adj_list[n].y;

                if (x==tri1.i and y==tri2.i) return true;
                if (x==tri1.i and y==tri2.j) return true;
                if (x==tri1.i and y==tri2.k) return true;
                if (x==tri1.j and y==tri2.i) return true;
                if (x==tri1.j and y==tri2.j) return true;
                if (x==tri1.j and y==tri2.k) return true;
                if (x==tri1.k and y==tri2.i) return true;
                if (x==tri1.k and y==tri2.j) return true;
                if (x==tri1.k and y==tri2.k) return true;

                if (y==tri1.i and x==tri2.i) return true;
                if (y==tri1.i and x==tri2.j) return true;
                if (y==tri1.i and x==tri2.k) return true;
                if (y==tri1.j and x==tri2.i) return true;
                if (y==tri1.j and x==tri2.j) return true;
                if (y==tri1.j and x==tri2.k) return true;
                if (y==tri1.k and x==tri2.i) return true;
                if (y==tri1.k and x==tri2.j) return true;
                if (y==tri1.k and x==tri2.k) return true;
        }
        return false;
}

bool tri_common_element(triangles t1,triangles t2)
{
        if (t1.i==t2.i or t1.i==t2.j or t1.i==t2.k) return true;
        if (t1.j==t2.i or t1.j==t2.j or t1.j==t2.k) return true;
        if (t1.k==t2.i or t1.k==t2.j or t1.k==t2.k) return true;
        return false;
}



////////////////////////////////////////////////

void routine_parameters(	string input_filename,
				double &tol,
				bool &Use_ARPACK,
				int &Lanczos_iterations,
				int &Lanczos_ncycles)
{
	ifstream infile(input_filename.c_str());
	string line;

	string stol("tol=");
	string sUse_ARPACK("Use_ARPACK=");
	string sLanczos_iterations("Lanczos_iterations=");
	string sLanczos_ncycles("Lanczos_ncycles=");

	while( getline(infile,line) )
	{
		istringstream linestream(line);
		char r;
		string read;

		while(linestream >> r)
		{
			read.push_back(r);
	
			if(read.compare(stol) == 0)
			{	linestream >> tol;}
			

			if(read.compare(sUse_ARPACK) == 0)
			{
				linestream >> Use_ARPACK;
				read.clear();
			}

			if(read.compare(sLanczos_iterations) == 0)
			{
				linestream >> Lanczos_iterations;
				read.clear();
			}
			if(read.compare(sLanczos_ncycles) == 0)
			{
				linestream >> Lanczos_ncycles;
				read.clear();
			}
		
		}
		

	}
}

void ED_2D_Lattice_Geometry(	string input_filename, 
				int &n_sites,
				double &Jz,
				double &J1,
				double &J2,
				vector <coordinates> &adj_list,
				double &Ch,
				vector <triangles> &ch_list, 
				double &ChCh,
				vector <bowties> &chch_list, 
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
	string sChCh("ChCh=");
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

	std::vector<triangles> tri_list;	
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
						double ph;
						linestream >> op_b >> x >> comma >> y >> comma >> ph >> cl_b;
						outfile <<"   "<< x << "   " << y << "   " << ph << endl;
						
						coordinates alist; 
						alist.x = x; alist.y = y; alist.J = J1; alist.expphase = exp(complex<double>(0.0, ph));
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
			
			if(read.compare(sChCh) == 0)
			{	linestream >> ChCh;}
			
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
						t1.i = x; t1.j = y; t1.k = z; t1.expphase=get_expphase(t1,adj_list); 
						t2.i = y; t2.j = z; t2.k = x; t2.expphase=get_expphase(t2,adj_list);
						t3.i = z; t3.j = x; t3.k = y; t3.expphase=get_expphase(t3,adj_list);

						//outfile<<"t1.exphhase = "<<t1.expphase<<endl;
						//outfile<<"t2.exphhase = "<<t2.expphase<<endl;
						//outfile<<"t3.exphhase = "<<t3.expphase<<endl;
						tri_list.push_back(t1);
						ch_list.push_back(t1);
						ch_list.push_back(t2);
						ch_list.push_back(t3);	    
					}
					else if (r_tmp == ']')
					break;
				}
				read.clear();
			}
			
			outfile<<" Triangles "<<endl;
			for (int m=0;m<tri_list.size();m++)
			{
				outfile<< tri_list[m].i <<"  "<< tri_list[m].j <<"  "<< tri_list[m].k <<"  "<<endl;
			}	
			chch_list.clear();	
			for (int m=0;m<tri_list.size();m++)
			{
				for (int n=0;n<tri_list.size();n++)
				{
					triangles tri1=tri_list[m];	
					triangles tri2=tri_list[n];	
					if (m!=n and neighbor_tri(tri1,tri2,adj_list)==true)
					{
						bowties b1,b2,b3,b4,b5,b6,b7,b8,b9;
						b1.x = tri1.i; b1.y = tri1.j; b1.z = tri1.k; b1.p=tri2.i; b1.q=tri2.j ; b1.r = tri2.k ;
						b2.x = tri1.i; b2.y = tri1.j; b2.z = tri1.k; b2.p=tri2.k; b2.q=tri2.i ; b2.r = tri2.j ;
						b3.x = tri1.i; b3.y = tri1.j; b3.z = tri1.k; b3.p=tri2.j; b3.q=tri2.k ; b3.r = tri2.i ;
						
						b4.x = tri1.k; b4.y = tri1.i; b4.z = tri1.j; b4.p=tri2.i; b4.q=tri2.j ; b4.r = tri2.k ;
						b5.x = tri1.k; b5.y = tri1.i; b5.z = tri1.j; b5.p=tri2.k; b5.q=tri2.i ; b5.r = tri2.j ;
						b6.x = tri1.k; b6.y = tri1.i; b6.z = tri1.j; b6.p=tri2.j; b6.q=tri2.k ; b6.r = tri2.i ;
						
						b7.x = tri1.j; b7.y = tri1.k; b7.z = tri1.i; b7.p=tri2.i; b7.q=tri2.j ; b7.r = tri2.k ;
						b8.x = tri1.j; b8.y = tri1.k; b8.z = tri1.i; b8.p=tri2.k; b8.q=tri2.i ; b8.r = tri2.j ;
						b9.x = tri1.j; b9.y = tri1.k; b9.z = tri1.i; b9.p=tri2.j; b9.q=tri2.k ; b9.r = tri2.i ;
						
						chch_list.push_back(b1);	
						chch_list.push_back(b2);	
						chch_list.push_back(b3);	
						chch_list.push_back(b4);	
						chch_list.push_back(b5);	
						chch_list.push_back(b6);	
						chch_list.push_back(b7);	
						chch_list.push_back(b8);	
						chch_list.push_back(b9);	
					}
				}
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


void ED_3D_Lattice_Geometry(	string input_filename, 
				int &n_sites,
				double &Jz,
				double &J1,
				double &J2,
				vector <coordinates> &adj_list,
				double &Ch,
				vector <triangles> &ch_list, 
				double &ChCh,
				vector <bowties> &chch_list, 
				vector<int> &T1,
				vector<int> &T2,
				vector<int> &T3,
				vector<int> &Cut_Pts,
				double &Sz,
				int &kx,
				int &ky,
				int &kz,
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
	string sChCh("ChCh=");
	string sTriangles("Triangles=");
	

	string sT1("T1=");
	string sT2("T2=");
	string sT3("T3=");
	string sAll_Sectors("All_Sectors=");
	string sSpecific_Sz_Sector("Specific_Sz_Sector=");
	string sSpecific_pSector("Specific_pSector=");
	string sSz("Sz=");
	string skx("kx=");
	string sky("ky=");
	string skz("kz=");
	
	string sFind_Vectors("Find_Vectors=");	
	string sUnwrap_pVector("Unwrap_pVector=");
	string sNum_Eigenvalues("Num_Eigenvalues=");
	string sVector_File_Path("Vector_File_Path=");

	std::vector<triangles> tri_list;	
	tri_list.clear();
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
						double ph;
						linestream >> op_b >> x >> comma >> y >> comma >> ph >> cl_b;
						outfile <<"   "<< x << "   " << y << "   " << ph << endl;
						
						coordinates alist; 
						alist.x = x; alist.y = y; alist.J = J1; alist.expphase = exp(complex<double>(0.0, ph));
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
			
			if(read.compare(sChCh) == 0)
			{	linestream >> ChCh;}
			
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
						t1.i = x; t1.j = y; t1.k = z; t1.expphase=get_expphase(t1,adj_list); 
						t2.i = y; t2.j = z; t2.k = x; t2.expphase=get_expphase(t2,adj_list);
						t3.i = z; t3.j = x; t3.k = y; t3.expphase=get_expphase(t3,adj_list);

						//outfile<<"t1.exphhase = "<<t1.expphase<<endl;
						//outfile<<"t2.exphhase = "<<t2.expphase<<endl;
						//outfile<<"t3.exphhase = "<<t3.expphase<<endl;
						tri_list.push_back(t1);
						ch_list.push_back(t1);
						ch_list.push_back(t2);
						ch_list.push_back(t3);	    
					}
					else if (r_tmp == ']')
					break;
				}
				read.clear();
			}
			
			outfile<<" Triangles (3D)"<<endl;
			for (int m=0;m<tri_list.size();m++)
			{
				outfile<< tri_list[m].i <<"  "<< tri_list[m].j <<"  "<< tri_list[m].k <<"  "<<endl;
			}	
			chch_list.clear();	
			for (int m=0;m<tri_list.size();m++)
			{
				for (int n=0;n<tri_list.size();n++)
				{
					triangles tri1=tri_list[m];	
					triangles tri2=tri_list[n];	
					if (m!=n and neighbor_tri(tri1,tri2,adj_list)==true)
					{
						bowties b1,b2,b3,b4,b5,b6,b7,b8,b9;
						b1.x = tri1.i; b1.y = tri1.j; b1.z = tri1.k; b1.p=tri2.i; b1.q=tri2.j ; b1.r = tri2.k ;
						b2.x = tri1.i; b2.y = tri1.j; b2.z = tri1.k; b2.p=tri2.k; b2.q=tri2.i ; b2.r = tri2.j ;
						b3.x = tri1.i; b3.y = tri1.j; b3.z = tri1.k; b3.p=tri2.j; b3.q=tri2.k ; b3.r = tri2.i ;
						
						b4.x = tri1.k; b4.y = tri1.i; b4.z = tri1.j; b4.p=tri2.i; b4.q=tri2.j ; b4.r = tri2.k ;
						b5.x = tri1.k; b5.y = tri1.i; b5.z = tri1.j; b5.p=tri2.k; b5.q=tri2.i ; b5.r = tri2.j ;
						b6.x = tri1.k; b6.y = tri1.i; b6.z = tri1.j; b6.p=tri2.j; b6.q=tri2.k ; b6.r = tri2.i ;
						
						b7.x = tri1.j; b7.y = tri1.k; b7.z = tri1.i; b7.p=tri2.i; b7.q=tri2.j ; b7.r = tri2.k ;
						b8.x = tri1.j; b8.y = tri1.k; b8.z = tri1.i; b8.p=tri2.k; b8.q=tri2.i ; b8.r = tri2.j ;
						b9.x = tri1.j; b9.y = tri1.k; b9.z = tri1.i; b9.p=tri2.j; b9.q=tri2.k ; b9.r = tri2.i ;
						
						chch_list.push_back(b1);	
						chch_list.push_back(b2);	
						chch_list.push_back(b3);	
						chch_list.push_back(b4);	
						chch_list.push_back(b5);	
						chch_list.push_back(b6);	
						chch_list.push_back(b7);	
						chch_list.push_back(b8);	
						chch_list.push_back(b9);	
					}
				}
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
			
			if(read.compare(sT3) == 0)
			{
				char s;
			
				while(linestream >> s)
				{
					if(s == '[' || s == ',')
					{
						int x; linestream >> x;	  
						T3.push_back(x);
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

			if(read.compare(skz) == 0)
			{
				linestream >> kz;
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


	if(T1.size() != T2.size() or T1.size()!=T3.size() or T2.size()!=T3.size())
	{
		outfile << "\nERROR: input T1 != T2 or equivalent problem\n";
		exit(1);
	}
	else
		n_sites = T1.size();
}	

