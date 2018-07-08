#include "Basis_States.h"

using namespace std;


// btest64 tells you whether a bit is 0 or 1  at position pos

bool btest64(int64_t a,int pos)
{
	bitset<64> b;
	b=a;
	return b[pos];
}

// Takes an integer and sets a bit to 1 at location pos

int64_t ibset64(int64_t a,int pos)
{	return a| int64_t(1)<<pos;}

// Takes an integer and sets a bit to 0 at location pos

int64_t ibclr64(int64_t a,int pos)
{	return a & ~(int64_t(1)<<pos);}


int64_t flip(int64_t v, int i, int j,
	     bool ti, bool tj)
{
	if(ti == 1)
		v = ibclr64(v, i);
  	else
		v = ibset64(v, i);

	if(tj == 1)
		v = ibclr64(v, j);
	else
		v = ibset64(v, j);
  
	return v;  
}

void constrained_dets(int num_sites,int num_ones,std::vector<int64_t> &dets)
{
	// ----------------------------------------------------------------------------------------------
	//    ! Description   : Gives us a neat way of generating all configs with a fixed
	//    !                 particle number
	//    ! Author        : F. Petruzielo's code used by H.J. Changlani (ref. Bit Twiddling Hacks)
	//    ! ----------------------------------------------------------------------------------------------

	int64_t temp,temp1,temp2,temp3,temp4,temp5,temp6;
	int64_t i,num_configs;

	if (num_sites<num_ones)
	{
		cout<<"ERROR: Num sites < Num_ones"<<endl;
		return;
	}

	dets.clear();
	num_configs=n_choose_k(num_sites,num_ones);

	dets.push_back(pow(2,num_ones) - 1);
	for( i=1;i<num_configs;i++)
	{
		temp  = (dets[i-1] | dets[i-1] - 1) + 1;
		temp2 = (temp) & (-temp);
		temp3=  (dets[i-1]) & (-dets[i-1]);
		temp4=temp2/temp3;
		temp5=temp4>>1;
		temp6=temp5-1;
		dets.push_back(temp|temp6);
	}
}

