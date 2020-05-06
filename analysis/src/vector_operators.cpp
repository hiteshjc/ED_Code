#include "vector_operators.h"

void rotate_vector(	vector<ids> &psi,
			vector<int> &Rot,
			vector<ids> &rpsi)
{

	int64_t nstates = psi.size();
	rpsi.resize(nstates);

	#pragma omp parallel for
	for(int64_t i = 0; i < nstates; i++)
	{
		int64_t tmpid = psi[i].stateID;

		int64_t newid = 0;

		for(int64_t j = 0; j < Rot.size(); j++)
		{
			bool bval = btest64(tmpid, Rot[j]);
			if(bval == 1)
				newid = ibset64(newid, j);

		}

		int64_t loc = find_ids_state(psi, newid);

		if (loc >= 0)
		{
			rpsi[loc].stateID = newid;
			rpsi[loc].coeff = psi[i].coeff;
		}
	}

}


