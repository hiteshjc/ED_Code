#ifndef GUARD_Lanczos_sector_H
#define GUARD_Lanczos_sector_H

#include "global.h"
#include "LAPACK_matrix_functions.h"
#include "H_v.h"  
#include "pClasses.h"
#include "lanczos.h"
#include "Lin_Tables.h"

void sector_lanczos(	int &kx,
			int &ky,
			double &cSz,
			int &n_spins,
			double &lambda,
			double &J2,
			vector<coordinates> &adj_list,
			double &Ch,
			vector<triangles> &ch_list,
			double &ChCh,
			vector<bowties> &chch_list,
			vector<int> &T1,
			vector<int> &T2,
			int &n_sec_ev,
			bool &find_ev,
			bool &unwrap_peigvec,
			string &filepath_EV,
			double &tol,
			int &iterations,
			int &ncycles,
			ofstream &outfile);

void sector_lanczos(	int &kx,
			int &ky,
			int &kz,
			double &cSz,
			int &n_spins,
			double &lambda,
			double &J2,
			vector<coordinates> &adj_list,
			double &Ch,
			vector<triangles> &ch_list,
			double &ChCh,
			vector<bowties> &chch_list,
			vector<int> &T1,
			vector<int> &T2,
			vector<int> &T3,
			int &n_sec_ev,
			bool &find_ev,
			bool &unwrap_peigvec,
			string &filepath_EV,
			double &tol,
			int &iterations,
			int &ncycles,
			ofstream &outfile);

#endif

