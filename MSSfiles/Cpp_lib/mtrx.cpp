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

#include "mtrx.h"
#include <math.h>


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//	CMatrix		Member functions
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

void CMatrix::Allocate(int rows, int cols)
{
	int i;
	int elements;
	VARTYPE*  Buf;

	elements = rows*cols;

	if (elements > 0)
	{
		// Allocate mem once for all rows
		Buf  = (VARTYPE*)mmalloc(sizeof(VARTYPE)*elements);

		// Initialise matrix to zero
		// for (i=0;i<elements;i++)
		// 	Buf[i] = 0;

		// Allocate mem for row-pointers (one-indexed)
		X = (VARTYPE**)mmalloc(sizeof(VARTYPE*)*rows);
		X--;

		// Initialise row-pointers
		X[1] = Buf-1;
		for (i=2;i<=rows;i++)
			X[i] = X[i-1] + cols;

		empty = false;
		m = rows;
		n = cols;
	}
	else
	{
		empty = true;
		m = n = 0;
	}
}


void CMatrix::Free()
{
	if (empty == false)
	{
		mfree(X[1]+1);
		mfree(X+1);

		empty = true;
		m = n = 0;
	}
}



CMatrix& CMatrix::operator = (const CMatrix& rhs)
{
	int i, j;
	
#ifdef COUTMESSAGES
	cout << "The operator= method was called" << endl;
#endif
	
	Free();						// Delete current contents
	Allocate (rhs.m, rhs.n);
	
	if (!empty)
	{
		for (i=1;i<=m;i++)
			for (j=1;j<=n;j++)
				X[i][j] = rhs.X[i][j];
	}
	return *this;
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//	CMatrix		Operators
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

CMatrix operator + (const CMatrix& lhs, const CMatrix& rhs)
{
	int i;

	int elements;
	VARTYPE *dest, *src1, *src2;
	
#ifdef COUTMESSAGES
	cout << "ADD function" << endl;
#endif
	
	// Verify non-empty quantities
	if ( !lhs.empty && !rhs.empty )
	{
		
		// Verify identical sizes
		if ( (lhs.m==rhs.m) && (lhs.n==rhs.n) )
		{
			CMatrix A(rhs.m,rhs.n);

			elements = A.m*A.n;					// Number of elements in matrices
			dest = A.X[1];						// First element of A (result)
			src1 = lhs.X[1];					// First element of Left
			src2 = rhs.X[1];					// First element of Right
			for (i=1;i<=elements;i++)
				dest[i] = src1[i] + src2[i];	// For each element: A = Left + Right
				
			return A;
		}
		else
			throw EMatrixCorrespond();
	}
	else
		throw EMatrixEmpty();
	
	CMatrix A;	// Empty matrix
	return A;
}



CMatrix operator - (const CMatrix& lhs, const CMatrix& rhs)
{
	int i;

	int elements;
	VARTYPE *dest, *src1, *src2;

#ifdef COUTMESSAGES
	cout << "SUB function" << endl;
#endif	

	// Verify non-empty quantities
	if ( !lhs.empty && !rhs.empty )
	{
		
		// Verify identical sizes
		if ( (lhs.m==rhs.m) && (lhs.n==rhs.n) )
		{
			CMatrix A(rhs.m,rhs.n);

			elements = A.m*A.n;					// Number of elements in matrices
			dest = A.X[1];						// First element of A (result)
			src1 = lhs.X[1];					// First element of Left
			src2 = rhs.X[1];					// First element of Right
			for (i=1;i<=elements;i++)
				dest[i] = src1[i] - src2[i];	// For each element: A = Left - Right
				
			return A;
		}
		else
			throw EMatrixCorrespond();
		
	}
	else
		throw EMatrixEmpty();
	
	CMatrix A;	// Empty matrix
	return A;
}


CMatrix operator * (const CMatrix& lhs, const CMatrix& rhs)
{
	int row;
	int col;
	int i;
	
	VARTYPE *dest;
	VARTYPE *src1;
	VARTYPE *src2;

	int m,n,end;


#ifdef COUTMESSAGES
	cout << "MULT function" << endl;
#endif
	
	// Verify non-empty quantities
	if ( !lhs.empty && !rhs.empty )
	{
		
		// Verify legal size
		if ( (lhs.n==rhs.m)  )
		{
			CMatrix A(lhs.m,rhs.n);
			dest = &A.X[1][1];							// Element (1,1) of result matrix
			m = A.m;
			n = A.n;
			end = lhs.n;


			// Perform multiplication Dest = Src1 * Src2
			// Src1 is traversed row-wise and Src2 column-wise as allways.
			// Pointer arithmetics is used to optimize speed.
			
			for (row=1;row<=m;row++)
			{
				// src1   = &lhs.X[row][1]-1;				// First element of Src1's row
				src1 = lhs.X[row];

				for (col=1;col<=n;col++)
				{
					*dest  = 0;							// Reset Dest's element to zero
					src2   = &rhs.X[1][col];			// First element of Src2's column

					for (i=1;i<=end;i++)				
					{
						*dest += src1[i] * *src2;		// Perform multiplication
						src2  += n;						// Point to next element in Src2
					}
					dest++;								// Point to next element in Dest
				}
			}
			return A;
			
		}
		else
			throw EMatrixCorrespond();
		
	}
	else
		throw EMatrixEmpty();
	
	CMatrix A;	// Empty matrix
	return A;
}


CMatrix operator * (const CMatrix& lhs, const VARTYPE& rhs)
{
	int i, elements;

	VARTYPE *dest;

#ifdef COUTMESSAGES
	cout << "MULT with scalar" << endl;
#endif
	
	// Verify non-empty matrix
	if ( !lhs.empty )
	{
		CMatrix A = lhs;
		elements = lhs.m*lhs.n;			// Number of elements in matrix
		dest = &A.X[1][1];				// Point to first element
		for (i=1;i<=elements;i++)		
			*dest++ *= rhs;				// For each element, multiplicate with "rhs"
		return A;
	}
	else
		throw EMatrixEmpty();
	
	CMatrix A;	// Empty matrix
	return A;
}


CMatrix operator | (const CMatrix& lhs, const CMatrix& rhs)
{
	int i;

	int elements;
	VARTYPE *dest, *src1, *src2;
	
#ifdef COUTMESSAGES
	cout << "Element multiplication .*" << endl;
#endif
	
	// Verify non-empty quantities
	if ( !lhs.empty && !rhs.empty )
	{
		
		// Verify identical sizes
		if ( (lhs.m==rhs.m) && (lhs.n==rhs.n) )
		{
			CMatrix A(rhs.m,rhs.n);

			elements = A.m*A.n;					// Number of elements in matrices
			dest = A.X[1];						// First element of A (result)
			src1 = lhs.X[1];					// First element of Left
			src2 = rhs.X[1];					// First element of Right
			for (i=1;i<=elements;i++)
				dest[i] = src1[i] * src2[i];	// For each element: A = Left * Right
				
			return A;
		}
		else
			throw EMatrixCorrespond();
	}
	else
		throw EMatrixEmpty();
	
	CMatrix A;	// Empty matrix
	return A;
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//	CMatrix		Functions
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


CMatrix CMatrix::Transp() const
{
	int i;
	int j;

	if (!empty)
	{
		CMatrix A(n,m);

		for (i=1;i<=m;i++)
			for (j=1;j<=n;j++)
				A.X[j][i] = X[i][j];
		return A;
	}
	else
		throw EMatrixEmpty();

	CMatrix A;
	return A;
}



//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//	CMatrix		Stream capabilities
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


ostream& operator << (ostream &ostrm, CMatrix& rhs)
{
	int i, j;
	
	if ( rhs.empty == false )
		for (i=1; i<=rhs.m; i++)
		{
			for (j=1;j<=rhs.n; j++)
				ostrm << rhs.X[i][j] << ' ';
			ostrm << endl;
		}
		
		return ostrm;
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//	CMatrix		Matlab support (only if Matlab is target)
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

#ifdef MATLAB_MEX_FILE
CMatrix::CMatrix(double *rhs, int rows, int cols)
{
	int i,j,k;
	
	Allocate(rows, cols);
	
	// Copy matrix (Matlab stores'em column-wise)
	k = 0;
	for (j=1;j<=cols;j++)
	{
		for (i=1;i<=rows;i++)
		{
			X[i][j] = rhs[k];
			k++;
		}
	}
	
	#ifdef COUTMESSAGES
	cout << "CMatrix MATLAB constructor called with m=" << rows << " and n=" << cols << endl;
	#endif
}


void CMatrix::Map2Matlab(double *rhs) const
{
	int i,j;
	
	// Assumes that enough space allready has been allocated with mxCreateDoubleMatrix()
	
	// Copy contents
	for (j=1;j<=n;j++)
	{
		for (i=1;i<=m;i++)
		{
			*rhs = X[i][j];
			rhs++;
		}
	}
}
#endif
