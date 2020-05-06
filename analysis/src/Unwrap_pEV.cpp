
#include "global.h"
#include "idsClass.h"
#include "read_inputfile.h"
#include "unwrap_pstate.h"

int main(int argc, char *argv[])
{
	string filein, filepEV, fileout;

	if (argc <= 1)
	{
		cout << "Usage: " << argv[0] << " <Filename>" << endl;
		exit(1);
	}

	filein=argv[1];
	fileout=argv[2];

	ofstream outfile(fileout.c_str());

	vector<int> T1, T2;
	read_T1T2(filein, T1, T2);
	read_pEV(filein, filepEV);

	/*outfile << " Input pEV file: " << filepEV << endl;

	outfile << "\nT1\t";	
	for(int i = 0; i < T1.size(); i++)
		outfile << T1[i] << " ";
	outfile << endl;

	outfile << "\nT2\t";	
	for(int i = 0; i < T2.size(); i++)
		outfile << T2[i] << " ";*/
	outfile << endl;
	vector<ids> psi;

	double ts_uw = wall_clock();
	read_n_unwrap(filein, filepEV, T1, T2, psi, outfile);
	ts_uw = wall_clock() - ts_uw;

	//outfile << "\n Size of unwrapped vector " << psi.size() << endl;
	//outfile << "\n Total time for unwrapping p-vector  " << (ts_uw)/60.0 << " minutes\n";
		
	//outfile << "\nEIGENVECTOR\n";
	Matrix mat(T1.size(),T1.size());
	zMatrix ninj(T1.size(),T1.size());
	for (int i=0;i<T1.size()*T1.size();i++) {mat[i]=0.0;ninj[i]=0.0;}
	for(int64_t jb = 0; jb < psi.size(); jb++)
	{
		int ctr=0;
		int el1,el2;
		for (int k=0;k<T1.size();k++)
		{
			if (btest64(psi[jb].stateID,k)==0 and ctr==1) {el2=k;ctr=ctr+1;}
			if (btest64(psi[jb].stateID,k)==0 and ctr==0) {el1=k;ctr=ctr+1;}
			//outfile << btest64(psi[jb].stateID,k); //<< "\t" << fixed << setprecision(15) << psi[jb].coeff << endl;

		}
		if (abs(psi[jb].coeff)>0.1 and real(psi[jb].coeff)>0) {mat(el1,el2)=+1;}
		if (abs(psi[jb].coeff)>0.1 and real(psi[jb].coeff)<0) {mat(el1,el2)=-1;}
		mat(el2,el1)=mat(el1,el2);
		ninj(el1,el2)+=1.0*psi[jb].coeff*conj(psi[jb].coeff);
		ninj(el2,el1)=ninj(el1,el2);	
		//outfile << "\t"<<fixed << setprecision(15) << real(psi[jb].coeff) << "  "<< imag(psi[jb].coeff)<<endl;
		outfile << fixed<< setprecision(15) << real(psi[jb].coeff) << "   "<< fixed << setprecision(15)<< imag(psi[jb].coeff)<<endl;
		//outfile << psi[jb].stateID << "\t" << fixed << setprecision(15) << psi[jb].coeff << endl;
	}
	outfile<<endl;
	for (int i=0; i<T1.size();i++)
	{
		for (int j=0; j< T1.size();j++)
		{
			if (mat(i,j)<0) outfile<<int(mat(i,j)*1.0001)<<" ";
			else outfile<<"+"<<int(mat(i,j)*1.0001)<<" ";
			
		}
		outfile<<endl;
	}
	outfile<<endl;	
	for (int i=0; i<T1.size();i++)
	{
		for (int j=0; j< T1.size();j++)
		{
			outfile<<real(ninj(i,j))<<" ";
			
		}
		outfile<<endl;
	}	
	outfile<<endl;	
	for (int i=0; i<T1.size();i++)
	{
		for (int j=0; j< T1.size();j++)
		{
			outfile<<i<<" "<<j<<" "<<real(ninj(i,j))<<endl;
			
		}
		outfile<<endl;
	}	
	outfile.close();


	return 0;
}
