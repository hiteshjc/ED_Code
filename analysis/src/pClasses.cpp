#include "pClasses.h"
// Vector for storing Quantum numbers and energy eigenvalues

void pVector::swapVij(	vector<nQuantum>& E,
			int64_t i,
			int64_t j)
{
	nQuantum tmp = E[i];
	E[i] = E[j];
	E[j] = tmp;
}
 
void pVector::merge(	vector<nQuantum>& E,
			vector<nQuantum> E1,
			int64_t n1,
			vector<nQuantum> E2,
			int64_t n2)
{
  for (int i = 0, i1 = 0, i2 = 0; i < n1 + n2; i++)
  {
    if (i1 < n1 && i2 < n2)
    {
      if (E1[i1].En > E2[i2].En)
	E[i] = E2[i2++];
      else if(E1[i1].En == E2[i2].En)
      {
	if(E1[i1].Sz > E2[i2].Sz)
	  E[i] = E2[i2++];
	else
	  E[i] = E1[i1++];
      }
      else
	E[i] = E1[i1++];
    }
    else if (i1 < n1)
      E[i] = E1[i1++];
    else if (i2 < n2)
      E[i] = E2[i2++];	
    
  }
}

void pVector::merge_sort_En(	vector<nQuantum> &E,
				int64_t m)
{
  int64_t n1 = m/2, n2;
  if (m % 2 == 0)
    n2 = m/2;
  else
    n2 = m/2 + 1;
    
  vector<nQuantum> E1(n1), E2(n2);
  
  for (int64_t i = 0; i < n1; i++)
    E1[i] = E[i];
  for(int64_t i = 0; i < n2; i++)
    E2[i] = E[i + n1]; 

  if (n1 > 2)
    merge_sort_En(E1, n1);
  else if(n1 == 2)		// Base case with just two elements
    if (E1[0].En > E1[1].En)
      swapVij(E1,0,1);
    
  if (n2 > 2)
    merge_sort_En(E2, n2);
  else if (n2 == 2)		// Base case with just two elements
    if (E2[0].En > E2[1].En)
      swapVij(E2, 0, 1);
    
  merge(E, E1, n1, E2, n2);
    
  E1.clear();
  E2.clear();
}


void pVector::sort_by_Sz(	vector< pVector::nQuantum >& E)
{ 
  for(int64_t i = 0; i < E.size(); i++)
  {
    double tE = E[i].En;
    double tSz = E[i].Sz;
    int64_t j = i, k = i;
    while(abs(E[j].En-tE) < pres)
    {
      if(tSz > E[j].Sz)
      {	
	k = j;
	tSz = E[j].Sz;
      }
      j++;
    }    
    swapVij(E, k, i);    
  }
}


