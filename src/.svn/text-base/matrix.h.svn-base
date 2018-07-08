#ifndef MATRIX_H
#define MATRIX_H

// Evgenii B. Rudnyi, http:// MatrixProgramming.com

#include "global.h"
#include "idsClass.h"

using namespace std;

class Matrix : public std::vector<double>
{
	int64_t m;
	int64_t n;

public:
	Matrix() : m(0), n(0) {}
	Matrix(int64_t m_, int64_t n_) : std::vector<double>(m_*n_), m(m_), n(n_) {}
	Matrix(int64_t m_, int64_t n_, double val) : std::vector<double>(m_*n_, val), m(m_), n(n_) {}

	void resize(int64_t m_, int64_t n_)
	{	m = m_; n = n_; std::vector<double>::resize(m_*n_);}

	void reserve(int64_t m_, int64_t n_)
	{	std::vector<double>::reserve(m_*n_);}

	void clear()
	{	m = n = 0; std::vector<double>::clear();}

	int64_t NRows() const {return m;}
	int64_t NCols() const {return n;}

	double& operator()(int64_t i, int64_t j)
	{	return operator[](i + j*m);}

	const double& operator()(int64_t i, int64_t j) const
	{	return operator[](i + j*m);}

	void swap( Matrix &y)
	{
		std::vector<double>::swap(y);
		std::swap(n, y.n);
		std::swap(m, y.m);
	}

	void clearMemory()
	{
		Matrix empty;
		swap(empty);
	}
};

class zMatrix : public std::vector< complex< double > >
{
	int64_t m;
	int64_t n;

public:
	zMatrix() : m(0), n(0) {}
	zMatrix(int64_t m_, int64_t n_) : std::vector< complex<double> >(m_*n_), m(m_), n(n_) {}
	zMatrix(int64_t m_, int64_t n_, double val) : std::vector< complex<double> >(m_*n_, val), m(m_), n(n_) {}

	void resize(int64_t m_, int64_t n_)
	{m = m_; n = n_; std::vector< complex<double> >::resize(m_*n_);}
	void reserve(int64_t m_, int64_t n_)
	{std::vector< complex<double> >::reserve(m_*n_);}
	void clear()
	{m = n = 0; std::vector< complex<double> >::clear();}

	int64_t NRows() const {return m;}
	int64_t NCols() const {return n;}

	complex<double>& operator()(int64_t i, int64_t j)
	{	return operator[](i + j*m);}

	const complex<double>& operator()(int64_t i, int64_t j) const
	{	return operator[](i + j*m);}
	/*  
	void swap( Matrix &y)
	{
	vector< complex<double> >::swap(y);
	swap(n, y.n);
	std::swap(m, y.m);
	}
	*/  
	void clearMemory()
	{
		zMatrix empty;
	//    swap(empty);
	}	
};



// Sparse Matrix that stores entries as ids (stateID + coeff)
class sMatrix{
private:
	vector< vector<ids> > states;
	int64_t Nr, Nc;

public:
	sMatrix(): Nr(0), Nc(0) {}
	sMatrix(int64_t m, int64_t n){
		states.resize(m);
		Nr = m;	Nc = n;
	}

	void clear()
	{	states.clear(); }

	void resize(int64_t m, int64_t n){
		states.clear();	states.resize(m);
		Nr = m;	Nc = n;
	}

	int64_t NRows(){return Nr; }
	int64_t NCols(){return Nc; }
	int64_t RowSize(int64_t i){return states[i].size(); }

	ids me(int64_t i, int64_t j)	
	{
		return states[i][j];
	}

	void push_back(int64_t id, int64_t stateID, complex<double> coeff){
		if(id < Nr)
		{
			ids tmp; tmp.stateID = stateID; tmp.coeff = coeff;
			states[id].push_back(tmp);
		}
		else
			cout << "\nError in Sparse Matrix Push_Back \n";
	}

	vector<ids> return_sVector(int64_t id)
	{ return states[id];	}

};


#endif
