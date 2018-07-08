#include "Simple_Math_Func.h"

using namespace std;

long long int_power(long long x, long long y)
{
	if (y == 0)
		return 1;
	else if (y == 1)
		return x;
	else return x*int_power(x, y-1);
}

int64_t int64_t_power(int64_t x, int64_t y)
{
	if (y == 0)
		return 1;
	else if (y == 1)
		return x;
	else return x*int64_t_power(x, y-1);
}


size_t sz_t_power(size_t x, size_t y)
{
	if (y == 0)
		return 1;
	else if (y == 1)
		return x;
	else return x*sz_t_power(x, y-1);
}

size_t sz_t_factorial(size_t n)
{
	if (n == 1 || n == 0)
		return 1;
	else 
		return n*sz_t_factorial(n-1);  
}

int64_t factorial(int64_t n)
{
	if (n == 1 || n == 0)
    		return 1;
	else
		return n*factorial(n-1);  
}

int64_t n_choose_k(int n_sites, int n_ones)
{
	if (n_ones == 0)
		return 1;
    	else
		return (n_sites * n_choose_k(n_sites - 1, n_ones - 1)) / n_ones;
}
