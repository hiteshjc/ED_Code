#include "global.h"

// Function that returns the current clock time
double wall_clock()
{
	struct timeval tp;
	double sec, usec, time;
	gettimeofday(&tp, NULL);
		
	sec = static_cast<double> (tp.tv_sec);
	usec = static_cast<double> (tp.tv_usec)/1E6;
	time = sec + usec;
	return time;
}	

