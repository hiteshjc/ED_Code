
#include "global.h"			// Header file with Matrix Class
#include "read_input_file.h"		// Functions that reads in the input parameters from input file
#include "ARPACK_full_eigs.h"		// Functions for ARPACK computation for eigs for all sectors
#include "ARPACK_sector.h"		// Functions for computing eigs and eigvecs for specific sectors
#include "Lanczos_sector.h"
#include "idsClass.h"

int main(int argc, char *argv[])
{
	string filenameED,filenameP;

	if (argc <= 1)
	{
		cout << "Usage: " << argv[0] << " <Filename>" << endl;
		exit(1);
	}

	filenameED=argv[1];
	filenameP=argv[2];

////////////////////////////// ARPACK ROUTINE OPTIONS //////////////////////////////

	bool Use_ARPACK = false;
	int maxiter = 5000; 			// Max Arnoldi iterations (re-orthongonalization iterations) 
////////////////////////////// LANCZOS PARAMETERS //////////////////////////////

	int Lanczos_ncycles = 10;
	int Lanczos_iterations = 100;

////////////////////////////////////////////////////////////////////////////////

	int n_ev = 8; 				// Default # of eigenvalues to compute per sector
	double tol = 10e-12; 			// Tolerance/Accuracy of eigenvalues computed

	routine_parameters( filenameED, tol, Use_ARPACK, Lanczos_iterations, Lanczos_ncycles);	

////////////////////////////// DEFAULT/INPUT_VALUES ////////////////////////////

	double Jz = 1;	 			// Default lambda = Jz component in XXZ model
	double J1 = 1.0;			// Default value for J1 
	double J2 = 0.0;			// Default value for J2
	double Ch = 0.0;			// Default value for term ch
	double ChCh = 0.0;			// Default value for term chch
	int n_spins = 0; 			// #of spins/sites in model
	vector <coordinates> adj_list; 		// Contains the list of bonds between vertices i & j and their strengths
	vector <triangles> ch_list;		// Contains list of ordered triangles for chirality term
	vector <bowties> chch_list;		// Contains list of ordered triangles for chiral-chiral term
	vector <int> T1, T2, T3; 		// Vectors to store translations along a_1 and a_2 directions
	vector <int> Cut;			// Vector that list the Cut for the Entanglement entropy
	
	// Sector input values
	double sSz = -1;			// Sector Sz 
	int skx = -1, sky = -1, skz =-1;	// Sector (kx, ky)
	int Num_Eigenvalues = 2;		// # of eigs and eigvecs to compute in sector
	bool Find_Vectors =0;			// Don't compute evecs by default
	bool Unwrap_pVector=1;			// Unwrap momentum eigenvector by default
	string Vector_File_Path;		// Path were eigvecs in sector are stored

	// Which sectors to compute
	bool All_Sectors=0;			// Compute ARPACK eigs for all (Sz and p) sectors
	double Sz_Sector = -1.0;		// Specify Sz sector		
	bool pSector=0;				// Specify Sz and p sector
	
	
	ofstream outfile(filenameP.c_str());

	// Function to read input parameters from input file //
	int dim=2;
	if (dim==2)
	{
	ED_2D_Lattice_Geometry(	filenameED, n_spins, Jz, J1, J2, adj_list, Ch, ch_list, ChCh, chch_list, T1, T2,
				Cut, sSz, skx, sky, All_Sectors, Sz_Sector, pSector,
				Num_Eigenvalues, Find_Vectors, Unwrap_pVector, Vector_File_Path, outfile);
	}
	if (dim==3)
	{
	outfile<<"Dimension = "<<dim<<endl;
	ED_3D_Lattice_Geometry(	filenameED, n_spins, Jz, J1, J2, adj_list, Ch, ch_list, ChCh, chch_list, T1, T2, T3,
				Cut, sSz, skx, sky, skz, All_Sectors, Sz_Sector, pSector,
				Num_Eigenvalues, Find_Vectors, Unwrap_pVector, Vector_File_Path, outfile);
	}

/////////////////////////////////////////////////////////////////////////////////////////


	outfile << "\nDiagonalizing XXZ Hamiltonian\n";	
	outfile << " J1 = " << J1 << endl;
	outfile << " J2 = " << J2 << endl;	
	outfile << " Jz = lambda = " << Jz << endl;
	outfile << " Ch = " << Ch << endl;

	outfile << " # of sites/spins = " << n_spins << endl;

	outfile << " Tolerance = " << tol << endl << endl;

	outfile << " Use ARPACK = " << Use_ARPACK << endl;

	outfile << " Lanczos iterations = " << Lanczos_iterations << endl;
	outfile << " Lanczos ncycles = " << Lanczos_ncycles << endl;
/*
	// Computes the eigenvalues for all the sectors //
	if(All_Sectors || Sz_Sector > -1.0){

		double time_start = wall_clock();
		pVector Eij;
		ARPACK(Eij, Jz, n_spins, Sz_Sector, adj_list, Ch, ch_list, T1, T2, n_ev, tol, maxiter, outfile);
		double time_end = wall_clock();
		outfile << "\n\n Total time for all sectors " << double(time_end - time_start)/60.0 << " minutes\n\n";

		Eij.sort_by_En();
		Eij.outfile(outfile, n_spins, Jz, J2);
	}
*/
	// Computes the eigenvalues and eigenvectors of a specific sector //
	if(pSector)
	{
		if(Use_ARPACK)
		{
			sector_eigs_eigvecs(	skx, sky, sSz, n_spins, Jz, J2, adj_list, Ch, ch_list, ChCh, chch_list, T1, T2, Num_Eigenvalues,
						Find_Vectors, Unwrap_pVector, Vector_File_Path, tol, maxiter, outfile);
		}
		else
		{
			if (dim==2) 
			{
				sector_lanczos(	skx, sky, sSz, n_spins, Jz, J2, adj_list, Ch, ch_list, ChCh, chch_list, T1, T2, Num_Eigenvalues,
						Find_Vectors, Unwrap_pVector, Vector_File_Path, tol, Lanczos_iterations, Lanczos_ncycles, outfile);
			}
			else
			{
				outfile<<"Here "<<endl;
				sector_lanczos(	skx, sky, skz, sSz, n_spins, Jz, J2, adj_list, Ch, ch_list, ChCh, chch_list, T1, T2, T3, Num_Eigenvalues,
						Find_Vectors, Unwrap_pVector, Vector_File_Path, tol, Lanczos_iterations, Lanczos_ncycles, outfile);


			}
		}
	}

	outfile.close();

	return 0;

}
