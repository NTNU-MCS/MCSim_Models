/*************************************************************************
 *
 * FILE		: $Archive: $ 
 * DATE		: 1999.02.04
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT	: A matrix class-library and appropriate functions.
 *			  Inspired by Matlab and Numerical Recipes.
 *
 *			  The most important properties of this library are:
 *				1) No differentiation between matrices and vectors.
 *				2) The matrices' sizes are kept internally, so the
 *				   user does not have to worry about getting them right.
 *				3) Operator overloading. Regular + - * etc. supported.
 *				4) Access to individual elements by using the Matlab-like
 *				   parenthesis notation. E.g:  a23=A(3,3)
 *				5) No error codes returned by the different functions.
 *				   Instead, an hierarchy of exceptions is used. The user
 *				   is them responsible for "catch"-ing these exceptions.
 *            
 * NOTES	: 
 *
 * AUTHOR(s): Karl-Petter Lindegaard
 *
 ************************************************************************/ 

#ifndef MTRX_H
#define MTRX_H


//#pragma warning( disable : 4290 ) 
#ifndef __cplusplus
#error "Must be C++ to use mtrx.h"
#endif

 
/*********************************************************************

  Constant definitions

*********************************************************************/

#define MAXITERATIONS 30		// Maximum number of iterations
#define EPSILON       1e-10		// Numeric resolution

#ifndef VARTYPE					// Define "double" as default type
#define VARTYPE double			// unless VARTYPE allready has been set
#endif



/*********************************************************************

  Include files

*********************************************************************/


#include "mmalloc.h"
#include "mtrxex.h"
//#include <iostream.h>
#include <iostream>;
using namespace std;

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif



/*********************************************************************

  Define a more readable matrix type

*********************************************************************/

#define matrix			CMatrix



/*********************************************************************

  Define the matrix class hierarchy

*********************************************************************/

/***********************************************
  CMatrixBase
***********************************************/

class CMatrixBase {
public:
	
	////////////////////////////////////////////////////
	//    Override New & Delete methods
	////////////////////////////////////////////////////

	void* operator new (size_t size) throw(EMatrixMemory) 
	{
		void *tmp = NULL;

		#ifdef COUTMESSAGES
		cout << "Operator new called." << endl;
		#endif

		if ( (tmp = mmalloc(size)) == NULL)
			throw EMatrixMemory();
		return mmalloc(size);
	}
	
	void  operator delete (void* p)
	{
		#ifdef COUTMESSAGES
		cout << "Operator delete called." << endl;
		#endif
		mfree (p);
	}
};


/***********************************************
  CMatrix
***********************************************/


class CMatrix : public CMatrixBase {

private:
	int    m, n;
	bool   empty;
	VARTYPE **X;
	
	////////////////////////////////////////////////////
	//    Allocate() and Free()
	////////////////////////////////////////////////////

	void Allocate(int rows, int cols);
	void Free();

public:
	////////////////////////////////////////////////////
	//    Constructors
	////////////////////////////////////////////////////
	
	CMatrix()
	{
		Allocate (0, 0);

		#ifdef COUTMESSAGES
		cout << "Default CMatrix constructor" << endl;
		#endif
	}
	
	CMatrix(int rows)
	{
		int i;

		Allocate (rows, 1);

		// Initialise vector to zero
		for (i=1;i<=m;i++)
			X[i][1] = 0;

		#ifdef COUTMESSAGES
		cout << "CMatrix constructor for vector, size " << rows << endl;
		#endif
	}
	
	CMatrix(int rows, int cols)
	{
		int i,j;

		Allocate(rows, cols);

		// Initialise matrix to zero
		for (i=1;i<=m;i++)
			for (j=1;j<=n;j++)
				X[i][j] = 0;

		#ifdef COUTMESSAGES
		cout << "CMatrix constructor called with m=" << m << " and n=" << n << endl;
		#endif
	}
	
	CMatrix(const CMatrix& rhs)
	{
		int i,j;
		
		#ifdef COUTMESSAGES
		cout << "CMatrix constructor for duplication called" << endl;
		#endif

		if ( !rhs.empty )
		{
			Allocate (rhs.m, rhs.n);
			
			for (i=1;i<=m;i++)
				for (j=1;j<=n;j++)
					X[i][j] = rhs.X[i][j];
		}
		else
		{
			m = n = 0;
			X = NULL;
			empty = true;
		}
	}
	
	
	~CMatrix()
	{
		#ifdef COUTMESSAGES
		cout << "Desctructor called for CMatrix with m=" << m << " and n=" << n << endl;
		#endif

		Free();
	}
	

	////////////////////////////////////////////////////
	//    Property describing methods
	////////////////////////////////////////////////////
	
	bool Isempty()	const	{ return empty; } 
	bool Iscolumn()	const	{ return ( (!empty) && (n==1) ); }
	bool Isrow()	const	{ return ( (!empty) && (m==1) ); }
	bool Isvector()	const	{ return ( (!empty) && ((m==1) || (n==1)) ); }
	bool Issquare()	const	{ return ( (!empty) && (m==n) ); }
	
	int Rows()		const	{ return m; }
	int Cols()		const	{ return n; }
	
	int Length()	const
	{
		if (Isvector())
			return m+n-1;
		else
			return 0;
	}
	
	
	////////////////////////////////////////////////////
	//    Access operators
	////////////////////////////////////////////////////
	
	VARTYPE& operator () (int row, int col) const
	{
		VARTYPE null = 0;
		if ( empty )
			throw EMatrixEmpty();
		else if ( ! ((row>0) && (row<=m) ) )
			throw EMatrixRow();
		else if ( ! ((col>0) && (col<=n) ) )
			throw EMatrixCol();
		else 
			return X[row][col];
		
		#ifdef COUTMESSAGES
		cout << "Reached end of operator()" << endl;		// Should not be reached!
		#endif

		//return null; 
	}
	
	VARTYPE& operator () (int elem) const
	{
		VARTYPE null = 0;
		if ( Isrow() )
		{
			if ( (elem>0) && (elem<=n) )
				return X[1][elem];
			else
				throw EMatrixCol();
		}
		else if ( Iscolumn() )
		{
			if ( (elem>0) && (elem<=m) )
				return X[elem][1];
			else
				throw EMatrixRow();
		}
		else
			throw EMatrixNovector();
		
		#ifdef COUTMESSAGES
		cout << "Reached end of operator()" << endl;		// Should not be reached!
		#endif

		//return null;
	}
	
	
	////////////////////////////////////////////////////
	//    Assignment operator
	////////////////////////////////////////////////////
	
	CMatrix& operator = (const CMatrix& rhs);
	
	
	////////////////////////////////////////////////////
	//    Basic operations
	////////////////////////////////////////////////////
	
	friend CMatrix operator + (const CMatrix& lhs, const CMatrix& rhs);
	friend CMatrix operator - (const CMatrix& lhs, const CMatrix& rhs);
	friend CMatrix operator * (const CMatrix& lhs, const CMatrix& rhs);
	friend CMatrix operator * (const CMatrix& lhs, const VARTYPE& rhs);
	friend CMatrix operator * (const VARTYPE& lhs, const CMatrix& rhs)
	{
		return rhs*lhs;
	}
	friend CMatrix operator | (const CMatrix& lhs, const CMatrix& rhs);	// Matlab's .*
	

	////////////////////////////////////////////////////
	//    Other operations
	////////////////////////////////////////////////////
	
	CMatrix Transp() const; 
	
	
	////////////////////////////////////////////////////
	//    Stream capabilities
	////////////////////////////////////////////////////
	
	friend ostream& operator << (ostream &ostrm, CMatrix& rhs);


	////////////////////////////////////////////////////
	//    Matlab support (only if Matlab is target)
	////////////////////////////////////////////////////

#ifdef MATLAB_MEX_FILE
	CMatrix(double *rhs, int rows, int cols);
	
	void Map2Matlab(double *rhs) const;
	// Assumes that enough space allready has been allocated with mxCreateDoubleMatrix()
#endif

	
};






/*********************************************************************

  Define different matrix functions ala Matlab.

*********************************************************************/

matrix Eye (int size);
matrix Ones (int rows, int cols);
matrix Zeros (int rows, int cols);
matrix Transpose(const matrix& A);

matrix Svd (const matrix& A, matrix& U, matrix &S, matrix &V);
matrix Inv (const matrix& A);

matrix Balance(const matrix& A);
matrix Hess(const matrix& A);
matrix Eig(const matrix& A);

VARTYPE Det(const matrix A);

matrix Interp1(const matrix& X, const matrix& Y, const matrix& xi);
matrix Interp2(const matrix& X, const matrix& Y, const matrix& Z, const matrix& XI, const matrix& YI);


#endif