#ifndef GUARD_read_input_file_h
#define GUARD_read_input_file_h

#include "global.h"


void routine_parameters(	string input_filename,
				double &tol,
				bool &Use_ARPACK,
				int &Lanczos_iterations,
				int &Lanczos_ncycles);

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
				ofstream &outfile);


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
				ofstream &outfile);


#endif
