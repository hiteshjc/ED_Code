#include "sort.h"

void swap_ids(vector<ids> &E, size_t i, size_t j)
{
  ids tmp = E[i];
  E[i] = E[j];
  E[j] = tmp;
}
 
void merge_ids(	vector<ids>& E, size_t m, 
		vector<ids> &E1, size_t n1,
		vector<ids> &E2, size_t n2)
{

	for (int64_t i = 0, i1 = 0, i2 = 0; i < n1 + n2; i++)
  	{
    		if (i1 < n1 && i2 < n2)
    		{
      			if (E1[i1].stateID > E2[i2].stateID)
				E[m + i] = E2[i2++];
			else
				E[m + i] = E1[i1++];
    		}
		else
    		{
     			if (i1 < n1 )
				E[m + i] = E1[i1++];
      
		      if (i2 < n2)
				E[m + i] = E2[i2++];	
    		}
  	}
}

void merge_sort_ids(vector< ids> &E, int64_t m, int64_t n)
{
	int64_t N = (n-m+1);
 
	if(N > 1)
	{
		int64_t n1 = N/2, n2;
		if (N % 2 == 0)
			n2 = N/2;
		else
			n2 = N/2 + 1;
  
		vector<ids> E1(n1), E2(n2);
  
		for (size_t i = 0; i < n1; i++)
    			E1[i] = E[i + m];
  		for(size_t i = 0; i < n2; i++)
    			E2[i] = E[i + m + n1]; 
		

		if (n1 > 2)
			merge_sort_ids(E1, 0, n1 - 1);
		else if(n1 == 2)		// Base case with just two elements
    			if (E1[0].stateID > E1[1].stateID)
      				swap_ids(E1,0,1);
    
		if (n2 > 2)
    			merge_sort_ids(E2, 0, n2 - 1);
  		else if (n2 == 2)		// Base case with just two elements
    		if (E2[0].stateID > E2[1].stateID)
      			swap_ids(E2, 0, 1);
    
		merge_ids(E, m, E1, n1, E2, n2);
    
  		E1.clear();
		E2.clear();
		
	}
}

