#ifndef GLOBAL_H
#define GLOBAL_H

#include <iostream>
#include <vector>
#include <complex>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "time.h"
#include "stdio.h"
#include "stdlib.h"
#include <omp.h>
#include <sys/time.h>
#include <sys/types.h>

# define M_PI           3.14159265358979323846  /* pi */
# define pres           10e-16  /* pi */

using namespace std;


struct coordinates
{
	int x, y;  
	double J;
};

struct triangles{
	int i, j, k;
};

double wall_clock();

struct links{
	int x, y;
};



#endif
