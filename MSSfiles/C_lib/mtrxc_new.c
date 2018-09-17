#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "mmalloc.h"
#include "mtrxc.h"

//////////////////////////////////////////////////////////////////
//                      Macro definititions:                    //
//////////////////////////////////////////////////////////////////

#define MATRIXID 1
#define VECTORID 2


//////////////////////////////////////////////////////////////////
//                       Type definititions:                    //
//////////////////////////////////////////////////////////////////

/*
union vectorheader {
	int n;
	Vector x;
};
*/

union matrixheader {
	struct {
		unsigned short m; 
		unsigned short n;
		unsigned char type; 
	} s;
	Matrix x;
};

// typedef union vectorheader VectorHeader;
typedef union matrixheader MatrixHeader;


//////////////////////////////////////////////////////////////////
//                       Size information                      //
//////////////////////////////////////////////////////////////////

int MtrxGetM(void *x)
{
	MatrixHeader *mh = ((MatrixHeader*)x)-1;
	return mh->s.m;
}

int MtrxGetN(void *x)
{
	MatrixHeader *mh = ((MatrixHeader*)x)-1;
	return mh->s.n;
}

int MtrxIsMatrix(void *M)
{
	MatrixHeader *mh = ((MatrixHeader*)M)-1;
	return mh->s.type == MATRIXID;
}

int MtrxIsVector(void *M)
{
	MatrixHeader *mh = ((MatrixHeader*)M)-1;
	return mh->s.type == VECTORID;
}

/*int MtrxGetVectorLen(Vector x)
{
	VectorHeader *vh = ((VectorHeader*)x)-1;
	return vh->n;
}
*/

//////////////////////////////////////////////////////////////////
//                   Allocation/Deallocation                    //
//////////////////////////////////////////////////////////////////

Vector VectorNew(int m)
{
	int nbytes;
	Vector v;
	MatrixHeader *mh;

	nbytes = sizeof(MatrixHeader)+sizeof(MtrxType)*(m+1);
	mh = (MatrixHeader*)mmalloc(nbytes);
	mh->s.m = (unsigned short)m;
	mh->s.n = 1;
	mh->s.type = VECTORID;
	v = (Vector)(mh+1);
	return v;
}

/*
Vector VectorNew(int n)
{
	int nbytes;
	Vector v;
	VectorHeader *vh;

	nbytes = sizeof(VectorHeader)+sizeof(MtrxType)*(n+1);
	vh = (VectorHeader*)mmalloc(nbytes);
	vh->n = n;
	v = (Vector)(vh+1);
	return v;
}
*/

Matrix MtrxNew(int m, int n)
{
	int nbytes, k;
	MatrixHeader *mh;
	Matrix matr;

	nbytes = sizeof(MatrixHeader)+sizeof(Vector)*(m+1);
	mh = (MatrixHeader*)mmalloc(nbytes);

	matr = (Matrix)(mh+1);
	nbytes = (m+1)*(n+1)*sizeof(MtrxType);
	matr[0] = (Vector)mmalloc(nbytes);
	for (k=1;k<=m;k++)
		matr[k] = matr[k-1]+(n+1) /* *sizeof(MtrxType) */ ;

	//for (k=0;k<=m;k++1)
	//	matr[k] = (Vector*)mmalloc((n+1)*sizeof(MtrxType));

	mh->s.m = (unsigned short)m;
	mh->s.n = (unsigned short)n;
	mh->s.type = MATRIXID;

	return matr;
}

/*
void VectorFree(Vector v)
{
	VectorHeader *vh;
	vh = ((VectorHeader*)v)-1;
	mfree(vh);
}
*/

void VectorFree(Vector v)
{
	MatrixHeader *mh;
	mh = ((MatrixHeader*)v)-1;
	mfree(mh);
}


void MtrxFree(Matrix v)
{
	// int k, m;
	MatrixHeader *mh;

	//m = MtrxGetM(v);
	//for (k=m;k>=0;k--)
	//	mfree(m[k]);

	mfree(v[0]);

	mh = ((MatrixHeader*)v)-1;
	mfree(mh);
}

void MtrxClear(void *A)
{
	int i, j, m, n;
	Matrix M;
	Vector v;

	m = MtrxGetM(A);
	if (MtrxIsMatrix(A))
	{
		n = MtrxGetN(A);
		M = (Matrix)A;
		for (i=1;i<=m;i++)
			for (j=1;j<=n;j++)
				M[i][j] = 0.0;
	}
	else
	{
		v = (Vector)A;
		for (i=1;i<=m;i++)
			v[i] = 0.0;
	}

}

/************************************************************************
 *
 * FUNCTION		: MtrxEye
 * DATE			: 2001.03.06
 *
 * Copyright Karl-Petter Lindegaard 2001
 *
 * ABSTRACT		: Allocate a new matrix and fill it with ones along the
 *                diagonal: An identity matrix.
 *				  
 * PARAMETERS	:
 *	size		-I	int			Dimension n
 *
 * RETURNS		:
 *				-O	Matrix		Identity matrix I (nxn)
 *
 * NOTES		:
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
Matrix MtrxEye(int size)
{
	int i;
	Matrix A;

	if (size > 0)
	{
		A = MtrxNew(size,size);
		for (i=1;i<=size;i++)
			A[i][i] = 1;
		return A;
	}
	return NULL;
}


/************************************************************************
 *
 * FUNCTION		: MtrxOnes
 * DATE			: 2001.03.06
 *
 * Copyright Karl-Petter Lindegaard 2001
 *
 * ABSTRACT		: Allocate a new matrix and fill it with ones.
 *				  
 * PARAMETERS	:
 *	rows		-I	int			Number of rows
 *	cols		-I	int			Number of columns
 *
 * RETURNS		:
 *				-O	Matrix		Matrix (rows x cols)
 *
 * NOTES		:
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
Matrix MtrxOnes(int rows, int cols)
{
	int i,j;
	Matrix A;

	if ( (rows > 0) && (cols > 0) )
	{
		A = MtrxNew(rows,cols);
		for (i=1;i<=rows;i++)
			for (j=1;j<=cols;j++)
				A[i][j] = 1.0;
		return A;
	}
	return NULL;
}


/************************************************************************
 *
 * FUNCTION		: MtrxZeros
 * DATE			: 2001.03.06
 *
 * Copyright Karl-Petter Lindegaard 2001
 *
 * ABSTRACT		: Allocate and return a Matrix filled with zeros
 *				  
 * PARAMETERS	:
 *	rows		-I	int			Number of rows
 *	cols		-I	int			Number of columns
 *
 * RETURNS		:
 *				-O	CMatrix		Matrix (rows x cols)
 *
 * NOTES		:
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
Matrix MtrxZeros(int rows, int cols)
{
	int i,j;
	Matrix A;

	if ( (rows > 0) && (cols > 0) )
	{
		A = MtrxNew(rows,cols);
		for (i=1;i<=rows;i++)
			for (j=1;j<=cols;j++)
				A[i][j] = 0.0;
		return A;
	}
	return NULL;
}

Matrix MtrxCopyContents(Matrix Y, Matrix A)
{
	int m, n, i, j;

	m = MtrxGetM(A);
	n = MtrxGetN(A);

	for (i=1;i<=m;i++)
		for (j=1;j<=n;j++)
			Y[i][j] = A[i][j];

	return Y;
}


//////////////////////////////////////////////////////////////////
// Matrix manipulation functions                                //
//////////////////////////////////////////////////////////////////


/************************************************************************
 *
 * FUNCTION		: MtrxSwapCols
 * DATE			: 2001.06.17
 *
 * Copyright Karl-Petter Lindegaard 2001
 *
 * ABSTRACT		: Swap two column of a of a matrix. Column i will be
 *				  interchanged with column j.
 *				  
 * PARAMETERS	:
 *	Y			-O	Matrix		Output matrix
 *	A			-I	Matrix		Input matrix
 *	i			-I	int			column #i
 *  j			-I	int			column #j
 *
 * RETURNS		:
 *				-O	Matrix		Pointer to matrix Y
 *
 * NOTES		:
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
Matrix MtrxSwapCols(Matrix Y, Matrix A, int i, int j)
{
	MtrxType tmp;
	int k, l, m, n;

	if (Y!=NULL)
	{
		m = MtrxGetM(A);
		n = MtrxGetN(A);
		
		// If Y!=A then duplicate A to Y
		if ( Y!=A )
		{
			for (k=1;k<=m;k++)
				for (l=1;l<=n;l++)
					Y[k][l] = A[k][l];
		}

		// Start swapping columns
		for (k=1;k<=m;k++)
		{
			tmp = Y[k][i];
			Y[k][i] = Y[k][j];
			Y[k][j] = tmp;
		}
	}

	return Y;
}


/************************************************************************
 *
 * FUNCTION		: MtrxSwapRows
 * DATE			: 2001.06.17
 *
 * Copyright Karl-Petter Lindegaard 2001
 *
 * ABSTRACT		: Swap two rows of a of a matrix. Row i will be
 *				  interchanged with row j.
 *				  
 * PARAMETERS	:
 *	Y			-O	Matrix		Output matrix
 *	A			-I	Matrix		Input matrix
 *	i			-I	int			row #i
 *  j			-I	int			row #j
 *
 * RETURNS		:
 *				-O	Matrix		Pointer to matrix Y
 *
 * NOTES		:
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
Matrix MtrxSwapRows(Matrix Y, Matrix A, int i, int j)
{
	MtrxType tmp;
	int k, l, m, n;

	if (Y!=NULL)
	{
		m = MtrxGetM(A);
		n = MtrxGetN(A);
		
		// If Y!=A then duplicate A to Y
		if ( Y!=A )
		{
			for (k=1;k<=m;k++)
				for (l=1;l<=n;l++)
					Y[k][l] = A[k][l];
		}

		// Start swapping rows
		for (l=1;l<=n;l++)
		{
			tmp = Y[i][l];
			Y[i][l] = Y[j][l];
			Y[j][l] = tmp;
		}
	}

	return Y;
}

//////////////////////////////////////////////////////////////////
// Elementary matrix functions                                  //
//////////////////////////////////////////////////////////////////

int MtrxEqualSizes(void *A, void *B)
{
	if ( (A!=NULL) && (B!=NULL) )
	{
		if ( (MtrxGetM(A) == MtrxGetM(B)) &&
			(MtrxGetN(A) == MtrxGetN(B)) )
			return 1;
		else
			return 0;
	}
	else
		return 0;
}


/************************************************************************
 *
 * FUNCTION		: MtrxTranspose
 * DATE			: 12.30.98
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT		: Return the transpose of a matrix. 
 *				  
 * PARAMETERS	:
 *  Y           -O  Matrix      Matrix where to put result
 *  A			-I  Matrix		Matrix to make the transpose of
 *
 * RETURNS		:
 *
 * NOTES		:
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
void *MtrxTranspose(Matrix Y, Matrix A)
{
	int i, j, m, n;
	
	if (MtrxIsMatrix(A))
	{
		m = MtrxGetM(A);
		n = MtrxGetN(A);
		if ( (MtrxGetM(Y) == n) &&
			(MtrxGetN(Y) == m))
			
		{
			for (i=1;i<=m;i++)
				for (j=1;j<=n;j++)
					Y[j][i] = A[i][j];
		}
	}
	return Y;
}


// Y = A+B
void *MtrxAdd(void *Y, void *A, void *B)
{
	int elements, i;
	Vector y, a, b;
	
	if (MtrxEqualSizes(A,B))
	{
		if (MtrxIsMatrix(A))
		{
			// In case this is a matrix
			elements = MtrxGetM(A)*(MtrxGetN(A)+1)-1;
			y = &(((Matrix)Y)[1][1]);
			a = &(((Matrix)A)[1][1]);
			b = &(((Matrix)B)[1][1]);
			//a = &A[1][1];
			//b = &B[1][1];
		}
		else
		{
			// Or A and B are vectors
			elements = MtrxGetM(A);
			y = &(((Vector)Y)[1]);
			a = &(((Vector)A)[1]);
			b = &(((Vector)B)[1]);
		}

		// Perform the addition
		for (i=0;i<elements;i++)
			y[i] = a[i]+b[i];
	}

	return Y;
}


// Y = A-B
void *MtrxSub(Matrix Y, Matrix A, Matrix B)
{
	int elements, i;
	Vector y, a, b;
	
	if (MtrxEqualSizes(A,B))
	{
		if (MtrxIsMatrix(A))
		{
			// In case this is a matrix
			elements = MtrxGetM(A)*(MtrxGetN(A)+1)-1;
			y = &(((Matrix)Y)[1][1]);
			a = &(((Matrix)A)[1][1]);
			b = &(((Matrix)B)[1][1]);
			//a = &A[1][1];
			//b = &B[1][1];
		}
		else
		{
			// Or A and B are vectors
			elements = MtrxGetM(A);
			y = &(((Vector)Y)[1]);
			a = &(((Vector)A)[1]);
			b = &(((Vector)B)[1]);
		}

		// Perform the addition
		for (i=0;i<elements;i++)
			y[i] = a[i]-b[i];
	}

	return Y;
}

// Y = c*A
void *MtrxScl(void *Y, void *A, MtrxType c)
{
	int elements, i;
	Vector y, a;
	
	if (MtrxIsMatrix(A))
	{
		// In case A is a matrix
		elements = MtrxGetM(A)*(MtrxGetN(A)+1)-1;
		y = &(((Matrix)Y)[1][1]);
		a = &(((Matrix)A)[1][1]);
	}
	else
	{
		// Or A is a vector
		elements = MtrxGetM(A);
		y = &(((Vector)Y)[1]);
		a = &(((Vector)A)[1]);
	}
	
	// Perform the addition
	for (i=0;i<elements;i++)
		y[i] = c*a[i];
	
	return Y;
}

// Y = A*B or y=A*b
void *MtrxMult(void *Y, Matrix A, void *B)
{
	int i, j, k, m, n, end; // row, col;
	Vector dest, src1; //, src2;
	Matrix Ym, Bm;
	
	// What kind of input arguments do we have?
	if (MtrxIsMatrix(Y) && MtrxIsMatrix(B))
	{
		m = MtrxGetM(A);
		n = MtrxGetN(B);
		end = MtrxGetN(A);

		Ym = (Matrix)(Y);
		Bm = (Matrix)(B);

		for (i=1;i<=m;i++)
			for (j=1;j<=n;j++)
			{
				Ym[i][j] = 0.0;
				for (k=1;k<=end;k++)
					Ym[i][j] += A[i][k]*Bm[k][j];
			}


		/* Optimized multiplication (still bugsy)
		dest = (Vector)&(((Matrix)Y)[1][1]);
		
		for (row=1;row<=m;row++)
		{
			
			src1 = (MtrxType*)&A[row][1];
			
			for (col=1;col<=n;col++)
			{
				*dest = 0.0;
				src2  = (Vector)&(((Matrix)B)[1][col]);
				
				for (i=1;i<=end;i++)
				{
					*dest += src1[i] * *src2;
					src2  += n;
				}
				dest++;
			}
		}
		*/
	}
	else if (MtrxIsVector(Y) && MtrxIsVector(B))
	{
		m = MtrxGetM(A);
		n = MtrxGetN(A);

		dest = (Vector)Y;
		src1 = (Vector)B;
		for (i=1;i<=m;i++)
		{
			dest[i] = 0.0;
			for (j=1;j<=n;j++)
				dest[i] += A[i][j]*src1[j];
		}
	}

	return Y;
}

/*
// y = A*b
void MtrxVec(Vector y, Matrix A, Vector b)
{
	int i, j, m, n;

	m = MtrxGetM(A);
	n = MtrxGetN(A);

	for (i=1;i<=m;i++)
	{
		y[i] = 0.0;
		for (j=1;j<=n;j++)
			y[i] += A[i][j]*b[j];
	}
}

*/

#define min(a, b)  (((a) < (b)) ? (a) : (b))
#define max(a, b)  (((a) > (b)) ? (a) : (b))


MTRXTYPE sign(MTRXTYPE a, MTRXTYPE b) 
/* returns Sign(b)*Abs(a) */
{
	if (b >= (MTRXTYPE)0.0)
		return fabs(a);
	else 
		return (-fabs(a));
}

	
MTRXTYPE pythag(MTRXTYPE a, MTRXTYPE b)
/* Computes sqrt(a^2 + b^2) without destructive underflow or overflow */
{
	MTRXTYPE at = fabs(a);
	MTRXTYPE bt = fabs(b);

	if (at > bt)
		return (at*sqrt(1.0+(bt*bt)/(at*at)));
	else if (bt == 0.0)
		return 0.0;
	else 
		return (bt*sqrt(1.0+(at*at)/(bt*bt)));
}


/************************************************************************
 *
 * FUNCTION		: MtrxSvdComp
 * DATE			: 12.30.98
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT		: Compute the singular value decomposition
 *				  (SVD) of a matrix, A = U * W * V(transpose)
 *				  where U and V are both orthonormal matrixes,
 *				  and W is a diagonal matrix containing the
 *				  singular values.
 *
 * PARAMETERS	:
 *	U			-O	Matrix	First output matrix (m x m)
 *	S			-O	Matrix	Diagonal matrix with singular values (m x n)
 *	V			-O	Matrix	Last output matrix (n x n).
 *	A			-I	Matrix	Matrix (m x n) to decompose
 *
 * RETURNS		:
 *
 * NOTES		: REQUIRES THAT m>n for [m,n]=size(A) !!!!!!
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
void MtrxSvdComp (Matrix U, Matrix S, Matrix V, Matrix A)
{
	int m, n, mnmin;
	int nm, l, k, j, jj, its, i, flag;
	MTRXTYPE c, f, g, h, s, x, y, z, anorm, scale;
	Vector rv1, W; //, W2;

	l = 0;
	nm = 0;

	m = MtrxGetM(A);
	n = MtrxGetN(A);
	mnmin = min(m,n);

	rv1 = VectorNew(m); MtrxClear(rv1);
	W   = VectorNew(m); MtrxClear(W);
	//W2  = VectorNew(mnmin);
	
	
	// Initalise output matrices
	MtrxCopyContents(U, A); // both dims of U must be larger than A's
	MtrxClear (S);
	MtrxClear (V);


	// Start SVD decomposition
	g     = 0;
	scale = 0;
	anorm = 0;
	

	// Householder reduction to bidiagonal form
	for (i=1; i<=n; i++)
	{
		l = i+1;
		rv1[i] = scale*g;
		g     = 0;
		s     = 0;
		scale = 0;
		if (i<=m)
		{
			for (k=i;k<=m;k++)
				scale += fabs(U[k][i]);
			if (scale != 0.0)
			{
				for (k=i;k<=m;k++)
				{
					U[k][i] /= scale;
					s += U[k][i]*U[k][i];
				}
				f = U[i][i];
				g = -sign(sqrt(s),f);
				h = f*g - s;
				U[i][i] = f-g;
				for (j=l; j<=n; j++)
				{
					s = 0;
					for (k=i;k<=m;k++)
						s += U[k][i]*U[k][j];
					f = s/h;
					for (k=i;k<=m;k++)
						U[k][j] += f*U[k][i];
				}
				for (k=i;k<=m;k++)
					U[k][i] *= scale;
			}
		}
		W[i] = scale*g;
		g     = 0;
		s     = 0;
		scale = 0;
		
		if ((i<=m) && (i != n))
		{
			for (k=l;k<=n;k++)
				scale += fabs(U[i][k]);
			if (scale != 0.0)
			{
				for (k=l;k<=n;k++)
				{
					U[i][k] /= scale;
					s += U[i][k]*U[i][k];
				}
				f = U[i][l];
				g = -sign(sqrt(s),f);
				h = f*g - s;
				U[i][l] = f-g;
				for (k=l;k<=n;k++)
					rv1[k] = U[i][k] / h;
				for (j=l; j<=m; j++)
				{
					s = 0;
					for (k=l;k<=n;k++)
						s += U[j][k]*U[i][k];
					for (k=l;k<=n;k++)
						U[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++)
					U[i][k] *= scale;
			}
		}
		anorm = max (anorm, (fabs(W[i]) + fabs(rv1[i])) );
	}
	
	
	// Accumulation of right-hand transformations
	for (i=n;i>=1;i--)
	{
		if (i<n)
		{
			if (g!=0)
			{
				for (j=l;j<=n;j++)
					V[j][i] = (U[i][j]/U[i][l]) / g;
				for (j=l;j<=n;j++)
				{
					s = 0;
					for (k=l;k<=n;k++)
						s += U[i][k]*V[k][j];
					for (k=l;k<=n;k++)
						V[k][j] += s*V[k][i];
				}
			}
			for (j=l;j<=n;j++)
			{
				V[i][j] = 0;
				V[j][i] = 0;
			}
		}
		V[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}
	
	// mnmin = min(m,n);
	
	
	// Accumulation of left-hand transformations
	for (i=mnmin;i>=1;i--)
	{
		l = i + 1;
		g = W[i];
		for (j=l;j<=n;j++)
			U[i][j] = 0;
		if (g != 0.0)
		{
			g = 1.0/g;
			for (j=l;j<=n;j++)
			{
				s = 0;
				for (k=l;k<=m;k++)
					s += U[k][i]*U[k][j];
				f = (s/U[i][i])*g;
				for (k=i;k<=m;k++)
					U[k][j] += f*U[k][i];
			}
			for (j=i;j<=m;j++)
				U[j][i] *= g;
		}
		else
			for (j=i;j<=m;j++)
				U[j][i] = 0;
			U[i][i] += 1;
	}
	
	
	// Diagonalization of bidiagonal form: Loop over singular values
	// and over allowed iterations.
	
	for (k=n;k>=1;k--)
	{
		for (its=1;its<=MTRXMAXITERATIONS;its++)
		{
			flag = 1;
			for (l=k;l>=1;l--)			// Test for splitting
			{
				nm = l-1;				// Note that rv1(1) is allways zeros
				if ( (fabs(rv1[l])+anorm) == anorm)
				{
					flag = 0;
					break;
				}
				if ( (fabs(W[nm])+anorm) == anorm)
					break;				// Goto step 1
			}
			
			/*Step 1:*/
			if (flag)
			{
				c = 0;
				s = 0;
				for (i=1;i<=k;i++)
				{
					f = s*rv1[i];
					rv1[i] = c*rv1[i];
					if ( (fabs(f)+anorm) == anorm)
						break;			// Goto step 2
					g = W[i];
					h = pythag(f,g);
					W[i] = h;
					h = 1/h;
					c = g*h;
					s = -(f*h);
					for (j=1;j<=m;j++)
					{
						y = U[j][nm];
						z = U[j][i];
						U[j][nm] = (y*c) + (z*s);
						U[j][i]  = -(y*s) + (z*c);
					}
				}
			} // if (flag)
			z = W[k];
			if (l == k)				// Convergence
			{
				if (z < 0.0)		// Sing. values are made nonnegative
				{
					W[k] = -z;
					for (j=1;j<=n;j++)
						V[j][k] = -V[j][k];
				}
				break;
			}
			if (its == MTRXMAXITERATIONS)
			{
				// No convergence in x SVDCMP iterations
				VectorFree(W);
				VectorFree(rv1);
				//assert(0);
				return;
				// throw EMatrixIteration();
			}
			
			x  = W[l];				// Shift from bottom 2 by 2 minor
			nm = k - 1;
			y  = W[nm];
			g  = rv1[nm];
			h  = rv1[k];
			f  = ((y - z)*(y + z) + (g-h)*(g+h)) / (2.0*h*y);
			g  = pythag(f, 1.0);
			f  = ((x - z)*(x + z) + h*((y/(f + sign(g,f))) - h))/x;
			c  = 1.0;
			s  = 1.0;
			for (j=l;j<=nm;j++)
			{
				i  = j + 1;
				g  = rv1[i];
				y  = W[i];
				h  = s*g;
				g  = c*g;
				z  = pythag(f,h);
				rv1[j] = z;
				c  = f/z;
				s  = h/z;
				f  = (x*c) + (g*s);
				g  = -(x*s) + (g*c);
				h  = y*s;
				y  *= c;
				for (jj=l;jj<=n;jj++)
				{
					x = V[jj][j];
					z = V[jj][i];
					V[jj][j] = (x*c) + (z*s);
					V[jj][i] = -(x*s) + (z*c);
				}
				z = pythag(f,h);
				W[j] = z;
				if (z != 0.0)
				{
					z = 1.0/z;
					c = f*z;
					s = h*z;
				}
				f = (c*g) + (s*y);
				x = -(s*g) + (c*y);
				for (jj=l;jj<=m;jj++)
				{
					y = U[jj][j];
					z = U[jj][i];
					U[jj][j] = (y*c) + (z*s);
					U[jj][i] = -(y*s) + (z*c);
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			W[k]   = x;
		}
	}


/*
	// Householder reduction to bidiagonal form
	for (i=1; i<=n; i++)
	{
		l = i+1;
		rv1[i] = scale*g;
		g     = 0;
		s     = 0;
		scale = 0;
		if (i<=m)
		{
			for (k=i;k<=m;k++)
				scale += fabs(U[k][i]);
			if (scale != 0.0)
			{
				for (k=i;k<=m;k++)
				{
					U[k][i] /= scale;
					s += U[k][i]*U[k][i];
				}
				f = U[i][i];
				g = -sign(sqrt(s),f);
				h = f*g - s;
				U[i][i] = f-g;
				for (j=l; j<=n; j++)
				{
					s = 0;
					for (k=i;k<=m;k++)
						s += U[k][i]*U[k][j];
					f = s/h;
					for (k=i;k<=m;k++)
						U[k][j] += f*U[k][i];
				}
				for (k=i;k<=m;k++)
					U[k][i] *= scale;
			}
		}
		W[i] = scale*g;
		g     = 0;
		s     = 0;
		scale = 0;
		
		if ((i<=m) && (i != n))
		{
			for (k=l;k<=n;k++)
				scale += fabs(U[i][k]);
			if (scale != 0.0)
			{
				for (k=l;k<=n;k++)
				{
					U[i][k] /= scale;
					s += U[i][k]*U[i][k];
				}
				f = U[i][l];
				g = -sign(sqrt(s),f);
				h = f*g - s;
				U[i][l] = f-g;
				for (k=l;k<=n;k++)
					rv1[k] = U[i][k] / h;
				for (j=l; j<=m; j++)
				{
					s = 0;
					for (k=l;k<=n;k++)
						s += U[j][k]*U[i][k];
					for (k=l;k<=n;k++)
						U[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++)
					U[i][k] *= scale;
			}
		}
		anorm = max (anorm, (fabs(W[i]) + fabs(rv1[i])) );
	}
	
	
	// Accumulation of right-hand transformations
	for (i=n;i>=1;i--)
	{
		if (i<n)
		{
			if (g!=0)
			{
				for (j=l;j<=n;j++)
					V[j][i] = (U[i][j]/U[i][l]) / g;
				for (j=l;j<=n;j++)
				{
					s = 0;
					for (k=l;k<=n;k++)
						s += U[i][k]*V[k][j];
					for (k=l;k<=n;k++)
						V[k][j] += s*V[k][i];
				}
			}
			for (j=l;j<=n;j++)
			{
				V[i][j] = 0;
				V[j][i] = 0;
			}
		}
		V[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}
	
	// mnmin = min(m,n);
	
	
	// Accumulation of left-hand transformations
	for (i=mnmin;i>=1;i--)
	{
		l = i + 1;
		g = W[i];
		for (j=l;j<=n;j++)
			U[i][j] = 0;
		if (g != 0.0)
		{
			g = 1.0/g;
			for (j=l;j<=n;j++)
			{
				s = 0;
				for (k=l;k<=m;k++)
					s += U[k][i]*U[k][j];
				f = (s/U[i][i])*g;
				for (k=i;k<=m;k++)
					U[k][j] += f*U[k][i];
			}
			for (j=i;j<=m;j++)
				U[j][i] *= g;
		}
		else
			for (j=i;j<=m;j++)
				U[j][i] = 0;
			U[i][i] += 1;
	}
	
	
	// Diagonalization of bidiagonal form: Loop over singular values
	// and over allowed iterations.
	
	for (k=n;k>=1;k--)
	{
		for (its=1;its<=MTRXMAXITERATIONS;its++)
		{
			flag = 1;
			for (l=k;l>=1;l--)			// Test for splitting
			{
				nm = l-1;				// Note that rv1(1) is allways zeros
				if ( (fabs(rv1[l])+anorm) == anorm)
				{
					flag = 0;
					break;
				}
				if ( (fabs(W[nm])+anorm) == anorm)
					break;				// Goto step 1
			}
			
			// Step 1:
			if (flag)
			{
				c = 0;
				s = 0;
				for (i=1;i<=k;i++)
				{
					f = s*rv1[i];
					rv1[i] = c*rv1[i];
					if ( (fabs(f)+anorm) == anorm)
						break;			// Goto step 2
					g = W[i];
					h = pythag(f,g);
					W[i] = h;
					h = 1/h;
					c = g*h;
					s = -(f*h);
					for (j=1;j<=m;j++)
					{
						y = U[j][nm];
						z = U[j][i];
						U[j][nm] = (y*c) + (z*s);
						U[j][i]  = -(y*s) + (z*c);
					}
				}
			} // if (flag)
			z = W[k];
			if (l == k)				// Convergence
			{
				if (z < 0.0)		// Sing. values are made nonnegative
				{
					W[k] = -z;
					for (j=1;j<=n;j++)
						V[j][k] = -V[j][k];
				}
				break;
			}
			if (its == MTRXMAXITERATIONS)
			{
				// No convergence in x SVDCMP iterations
				VectorFree(W);
				VectorFree(rv1);
				return;
				// throw EMatrixIteration();
			}
			
			x  = W[l];				// Shift from bottom 2 by 2 minor
			nm = k - 1;
			y  = W[nm];
			g  = rv1[nm];
			h  = rv1[k];
			f  = ((y - z)*(y + z) + (g-h)*(g+h)) / (2.0*h*y);
			g  = pythag(f, 1.0);
			f  = ((x - z)*(x + z) + h*((y/(f + sign(g,f))) - h))/x;
			c  = 1.0;
			s  = 1.0;
			for (j=l;j<=nm;j++)
			{
				i  = j + 1;
				g  = rv1[i];
				y  = W[i];
				h  = s*g;
				g  = c*g;
				z  = pythag(f,h);
				rv1[j] = z;
				c  = f/z;
				s  = h/z;
				f  = (x*c) + (g*s);
				g  = -(x*s) + (g*c);
				h  = y*s;
				y  *= c;
				for (jj=l;jj<=n;jj++)
				{
					x = V[jj][j];
					z = V[jj][i];
					V[jj][j] = (x*c) + (z*s);
					V[jj][i] = -(x*s) + (z*c);
				}
				z = pythag(f,h);
				W[j] = z;
				if (z != 0.0)
				{
					z = 1.0/z;
					c = f*z;
					s = h*z;
				}
				f = (c*g) + (s*y);
				x = -(s*g) + (c*y);
				for (jj=l;jj<=m;jj++)
				{
					y = U[jj][j];
					z = U[jj][i];
					U[jj][j] = (y*c) + (z*s);
					U[jj][i] = -(y*s) + (z*c);
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			W[k]   = x;
		}
	}
*/
	// Generate return vector and matrix with sing. val along diagonal
	for (i=1;i<=mnmin;i++)
		S[i][i] = /*W2[i] =*/ W[i];


	VectorFree(W);
	VectorFree(rv1);
	//return W2;	// Return vector of singular values.
}

/************************************************************************
 *
 * FUNCTION		: MtrxSvd
 * DATE			: 12.30.98
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT		: Compute the singular value decomposition
 *				  (SVD) of a matrix, A = U * W * V(transpose)
 *				  where U and V are both orthonormal matrixes,
 *				  and W is a diagonal matrix containing the
 *				  singular values.
 *
 * PARAMETERS	:
 *	U			-O	Matrix	First output matrix (m x m)
 *	S			-O	Matrix	Diagonal matrix with singular values (m x n)
 *	V			-O	Matrix	Last output matrix (n x n).
 *	A			-I	Matrix	Matrix (m x n) to decompose
 *
 * RETURNS		:
 *
 * NOTES		: 
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
void MtrxSvd (Matrix U, Matrix S, Matrix V, Matrix A)
{
	int m, n;
	Matrix A2, U2, S2, V2;

	m = MtrxGetM(A);
	n = MtrxGetN(A);

	if (m>=n)
		MtrxSvdComp (U, S, V, A);
	else
	{
		// SvdComp requires m>=n. Thus, run SVD on A' instead.
		// A = U*S*V' leads to A' = V*S'*U'
		A2 = MtrxNew(n, m);
		U2 = MtrxNew(n, n);
		S2 = MtrxNew(n, m);
		V2 = MtrxNew(m, m);

		MtrxTranspose(A2, A);
		MtrxSvdComp (U2, S2, V2, A2);

		MtrxCopyContents(U, V2);
		MtrxTranspose(S, S2);
		MtrxCopyContents(V, U2);

		MtrxFree(V2);
		MtrxFree(S2);
		MtrxFree(U2);
		MtrxFree(A2);
	}
}


/************************************************************************
 *
 * FUNCTION		: MtrxSvdSort
 * DATE			: 16.06.01
 *
 * Copyright Karl-Petter Lindegaard 2001
 *
 * ABSTRACT		: Sort the singular values of a sv-decompositon A=U*S*V'
 *				  such that the sigmas in S appear in decreasing order
 *                along the diagonal (as in Matlab).
 *
 * PARAMETERS	:
 *	U			-IO	Matrix	First output matrix (m x m)
 *	S			-IO	Matrix	Diagonal matrix with singular values (m x n)
 *	V			-IO	Matrix	Last output matrix (n x n).
 *
 * RETURNS		:
 *
 * NOTES		:
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
void MtrxSvdSort(Matrix U, Matrix S, Matrix V)
{
	int m, n, nmax;
	int i, j, isigmamax;
	MtrxType sigmamax, tmp;

	// Get dimension of S
	m = MtrxGetM(S);
	n = MtrxGetN(S);

	// Max-number of singular values.
	nmax = min(m,n);

	// Bobble-sort this thing
	for (i=1;i<=nmax-1;i++)
	{
		isigmamax = i;
		sigmamax  = S[i][i];

		for (j=i+1;j<=nmax;j++)
		{
			if (S[j][j] > sigmamax)
			{
				sigmamax = S[j][j];
				isigmamax = j;
			}
		}

		if ( (isigmamax != i) && (sigmamax > MTRXEPSILON) )
		{
			// Swap columns in U and V
			MtrxSwapCols (U, U, i, isigmamax);
			MtrxSwapCols (V, V, i, isigmamax);

			// Manual swap in S (less operations)
			tmp = S[i][i];
			S[i][i] = S[isigmamax][isigmamax];
			S[isigmamax][isigmamax] = tmp;
		}
	}
}



/************************************************************************
 *
 * FUNCTION		: MtrxInv
 * DATE			: 12.30.98
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT		: Compute the inverse of a matrix.
 *
 * PARAMETERS	:
 *  Y           -O  Matrix      Inverse
 *	A			-I	matrix		Matrix to compute inverse of.
 *
 * RETURNS		:
 *
 * NOTES		:
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
Matrix MtrxInv(Matrix Y, Matrix A)
{
	int m,n,i;
	Matrix U, S, V, Ut, TEMP;

	if (MtrxEqualSizes(Y,A))
	{
		m = MtrxGetM(A);
		n = MtrxGetN(A);

		// Allocate 
		U    = MtrxNew(m,m);
		Ut   = MtrxNew(m,m);
		S    = MtrxNew(m,n);
		V    = MtrxNew(n,n);
		TEMP = MtrxNew(n,n);

		// Compute SVD
		MtrxSvd(U,S,V,A);

		// Check for singularity
		for (i=1;i<=n;i++)
		{
			if ( fabs(S[i][i]) < MTRXEPSILON)
				S[i][i] = MTRXEPSILON;
			S[i][i] = 1.0/S[i][i];		//  Make reciprocal
		}

		// Y = U*S*V' where S contains 1/sigma
		MtrxMult(TEMP, V, S); // TEMP = V*S
		MtrxTranspose(Ut, U); 
		MtrxMult(Y, TEMP, Ut);

		MtrxFree(TEMP);
		MtrxFree(V);
		MtrxFree(S);
		MtrxFree(Ut);
		MtrxFree(U);
	}

	return Y;
}


/************************************************************************
 *
 * FUNCTION		: MtrxCross
 * DATE			: 21.05.01
 *
 * Copyright Karl-Petter Lindegaard 2001
 *
 * ABSTRACT		: Compute cross product of two vectors.
 *
 * PARAMETERS	:
 *  y           -O  Vector      Cross-product of a and b
 *	a			-I	Vector		Input vector
 *  a			-I	Vector		Input vector
 *
 * RETURNS		:
 *              -O  Vector		Cross-product, same as y
 *
 * NOTES		:
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
Vector MtrxCross(Vector y, Vector a, Vector b)
{
	Matrix S;

	if ((MtrxGetM(a) == 3) && (MtrxGetM(b) == 3))
	{
		// Allocate memory for S
		S = MtrxNew(3,3);

		// Create skew-symmetric S matrix for cross product
		MtrxCrossS(S, a);

		// Calculate cross product
		MtrxMult (y, S, b);

		// De-allocate S matrix
		MtrxFree(S);
	}

	return y;
}

/************************************************************************
 *
 * FUNCTION		: MtrxCrossS
 * DATE			: 21.05.01
 *
 * Copyright Karl-Petter Lindegaard 2001
 *
 * ABSTRACT		: Create skew-symmetric S matrix used for cross product.
 *                       | 0  -a3  a2 |
 *                  S =  | a3  0  -a1 |
 *                       |-a2  a1  0  |
 *
 * PARAMETERS	:
 *  S           -O  Matrix      The skew-symmetric S
 *	a			-I	Vector		Input vector
 *
 * RETURNS		:
 *              -O  Matrix		S
 *
 * NOTES		:
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
Matrix MtrxCrossS(Matrix S, Vector a)
{
	if (MtrxGetM(a) == 3)
	{
		S[1][1] = 0.0;   S[1][2] = -a[3]; S[1][3] = a[2];
		S[2][1] = a[3];  S[2][2] = 0.0;   S[2][3] = -a[1];
		S[3][1] = -a[2]; S[3][2] = a[1];  S[3][3] = 0.0;
	}

	return S;
}


/************************************************************************
 *
 * FUNCTION		: MtrxNorm2
 * DATE			: 21.05.01
 *
 * Copyright Karl-Petter Lindegaard 2001
 *
 * ABSTRACT		: Compute the 2-norm of a matrix or vector
 *
 * PARAMETERS	:
 *	X			-I	void*		Input vector/matrix
 *
 * RETURNS		:
 *              -O  MTRXTYPE	The 2-norm.
 *
 * NOTES		: Matrix norm not yet supported!
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
MTRXTYPE MtrxNorm2(void *X)
{
	MTRXTYPE r;
	Vector x;
	int n, k;

	r = 0.0;
	if (MtrxIsVector(X))
	{
		x = (Vector)X;
		n = MtrxGetM(x);

		// Calculate norm for vector
		for (k=1;k<=n;k++)
			r = r + x[k]*x[k];
	}
	else
	{
		// Not implemented yet.
	}

	r = sqrt(r);
	return r;
}


/************************************************************************
 *
 * FUNCTION		: MtrxNormalize2
 * DATE			: 21.05.01
 *
 * Copyright Karl-Petter Lindegaard 2001
 *
 * ABSTRACT		: Normalize a matrix or vector using the 2-norm
 *
 * PARAMETERS	:
 *  Y			-O  void*		Normalized vector/matrix
 *	X			-I	void*		Input vector/matrix
 *
 * RETURNS		:
 *              -O  void*		The normalized vector/matrix. Same as Y
 *
 * NOTES		: Matrices not yet supported!
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
void *MtrxNormalize2(void *Y, void *X)
{
	MTRXTYPE r;

	r = MtrxNorm2(X);
	if (r>MTRXEPSILON)
		MtrxScl(Y, X, 1.0/r);
	else
		MtrxScl(Y, X, 1);
	
	return Y;
}


/************************************************************************
 *
 * FUNCTION		: MtrxTrace
 * DATE			: 21.05.01
 *
 * Copyright Karl-Petter Lindegaard 2001
 *
 * ABSTRACT		: Compute trace of a diagonal matrix
 *
 * PARAMETERS	:
 *	X			-I	void*		Input vector/matrix
 *
 * RETURNS		:
 *              -O  MTRXTYPE	The trace
 *
 * NOTES		: 
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
MTRXTYPE MtrxTrace(Matrix X)
{
	MTRXTYPE r;
	int m, n, k;

	r = 0.0;
	m = MtrxGetM(X);
	n = MtrxGetN(X);

	if ( m == n)
	{
		// Calculate trace
		for (k=1;k<=n;k++)
			r = r + X[k][k];
	}
	return r;
}


/************************************************************************
 *
 * FUNCTION		: MtrxRotationEulerAngles
 * DATE			: 21.05.01
 *
 * Copyright Karl-Petter Lindegaard 2001
 *
 * ABSTRACT		: Compute a rotation matrix from Euler angles.
 *
 * PARAMETERS	:
 *	R			-O	Matrix		Rotation matrix
 * Theta		-I	Vector		Vector of 3 Euler angles
 *
 * RETURNS		:
 *              -O  Matrix	The rotation R
 *
 * NOTES		: 
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *				
 * HISTORY		:
 *
 ***********************************************************************/
Matrix MtrxRotationEulerAngles(Matrix R, Vector Theta)
{
	MTRXTYPE spsi, cpsi, sphi, cphi, stheta, ctheta;
	
	if (MtrxGetM(Theta) == 3)
	{
		sphi = sin(Theta[1]);
		cphi = cos(Theta[1]);
		stheta = sin(Theta[2]);
		ctheta = cos(Theta[2]);
		spsi = sin(Theta[3]);
		cpsi = cos(Theta[3]);
		
		R[1][1]=cpsi*ctheta; R[1][2]=-spsi*cphi+cpsi*stheta*sphi; R[1][3]=spsi*sphi+cpsi*cphi*stheta;
		R[2][1]=spsi*ctheta; R[2][2]=cpsi*cphi+spsi*stheta*sphi;  R[2][3]=-cpsi*sphi+spsi*cphi*stheta;
		R[3][1]=-stheta;     R[3][2]=ctheta*sphi;                 R[3][3]=cphi*ctheta;
	}

	return R;
}


/************************************************************************
 *
 * FUNCTION		: MtrxRotationQuaternions
 * DATE			: 21.05.01
 *
 * Copyright Karl-Petter Lindegaard 2001
 *
 * ABSTRACT		: Compute a rotation matrix from Quaternions.
 *
 * PARAMETERS	:
 *	R			-O	Matrix		Rotation matrix
 * q			-I	Vector		Vector of quaternions
 *
 * RETURNS		:
 *              -O  Matrix	The rotation R
 *
 * NOTES		: 
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
Matrix MtrxRotationQuaternions(Matrix R, Vector q)
{
	MTRXTYPE e1, e2, e3, eta;

	if (MtrxGetM(q) == 4)
	{
		e1   = q[1];
		e2   = q[2];
		e3   = q[3];
		eta  = q[4];
		
		R[1][1] = e1*e1-e2*e2-e3*e3+eta*eta;
		R[1][2] = 2.0*(e1*e2-e3*eta);
		R[1][3] = 2.0*(e1*e3+e2*eta);

		R[2][1] = 2.0*(e1*e2+e3*eta);
		R[2][2] = -e1*e1+e2*e2-e3*e3+eta*eta;
		R[2][3] = 2.0*(e2*e3-e1*eta);

		R[3][1] = 2.0*(e1*e3-e2*eta);
		R[3][2] = 2.0*(e2*e3+e1*eta);
		R[3][3] = -e1*e1-e2*e2+e3*e3+eta*eta;
	}

	return R;
}


/************************************************************************
 *
 * FUNCTION		: MtrxQuaternions2EulerAngles
 * DATE			: 21.05.01
 *
 * Copyright Karl-Petter Lindegaard 2001
 *
 * ABSTRACT		: Computes Euler angles Theta from a vector q of containing
 * 				  unit quaternions. Assumes theta<>90 deg.
 *
 * PARAMETERS	:
 *	Theta		-O	Vector		Vector of 3 Euler angles
 * q			-I	Vector		Vector of quaternions.
 *
 * RETURNS		:
 *              -O  Vector	Euler angles (same as Theta)
 *
 * NOTES		: 
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
Vector MtrxQuaternions2EulerAngles (Vector Theta, Vector q)
{
	Matrix R;

	if ( (MtrxGetM(Theta)==3) && (MtrxGetM(q)==4) )
	{
		R = MtrxNew(3,3);

		// Create rotation matrix from q
		MtrxRotationQuaternions(R, q);

		// TIF says:
		// phi   =  atan2( cos(theta)*sin(phi), cos(theta)*cos(phi) );
		// theta = -asin( -sin(theta) );
		// psi   =  atan2( sin(psi)*cos(theta), cos(psi)*cos(theta) ); 

		// phi, theta, psi is calculated from R as follows
		Theta[1] = atan2(R[3][2], R[3][3]);   // phi
		Theta[2] = -asin(R[3][1]);            // theta
		Theta[3] = atan2(R[2][1], R[1][1]);   // psi
		
		MtrxFree(R);
	}

	return Theta;
}

/************************************************************************
 *
 * FUNCTION		: MtrxQuest
 * DATE			: 21.05.01
 *
 * Copyright Karl-Petter Lindegaard 2001
 *
 * ABSTRACT		: Computes the best fit of quaternions "q" that represents
 *					  the rotation from frame "v" to "w" such that 
 *                  r^w = R(q)^w_v * r^v
 *               Handles an arbitrary number of vectors r_i. Simply stack
 *               the vectors columnwise in W and V.
 *
 * PARAMETERS	:
 * q			-O	Vector		Vector of quaternions.
 * W			-I	Matrix		A 3xn matrix containing r^w_i
 * V			-I Matrix		A 3xn matrix containing r^v_i
 *
 * RETURNS		:
 *              -O  Vector	of quaternions (same as q)
 *
 * NOTES		: Uses svd instead of eig() because eig() not supported by mtrxc.
 *				  It is still an open question whether this approach actually
 *				  works. It will fail if abs(lambdamax) < abs(lambdamin) !!!!
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
Vector MtrxQuest(Vector q, Matrix W, Matrix V)
{
	int i, n, index;
	Matrix B, S, K;
	Matrix Us, Ss, Vs;
	Vector tmp, w, v, z;
	Matrix TMP, I4;
	MTRXTYPE rho, lambda;

	if ( (MtrxGetM(q)==4) && 
		(MtrxGetM(W)==3) && 
		MtrxEqualSizes(W, V))
	{
		n = MtrxGetN(W);

		// Allocate matrices
		B = MtrxZeros(3,3);
		S = MtrxZeros(3,3);
		K = MtrxNew(4,4);
		Us = MtrxNew(4,4);
		Ss = MtrxNew(4,4);
		Vs = MtrxNew(4,4);
		TMP = MtrxNew(3,3);
		z  = VectorNew(3);
		I4 = MtrxEye(4);

		
		// Construct B
		w   = VectorNew(3);
		v   = VectorNew(3);
		tmp = VectorNew(3);
		MtrxClear(z);
		for (i=1;i<=n;i++)
		{
			w[1] = W[1][i];
			w[2] = W[2][i];
			w[3] = W[3][i];
			v[1] = V[1][i];
			v[2] = V[2][i];
			v[3] = V[3][i];

			// B = B + W(:,i)*V(:,i)';
			TMP[1][1] = w[1]*v[1]; TMP[1][2] = w[1]*v[2]; TMP[1][3] = w[1]*v[3]; 
			TMP[2][1] = w[2]*v[1]; TMP[2][2] = w[2]*v[2]; TMP[2][3] = w[2]*v[3]; 
			TMP[3][1] = w[3]*v[1]; TMP[3][2] = w[3]*v[2]; TMP[3][3] = w[3]*v[3]; 
			// HERE WAS A MAJOR BUG. TMP[3][2] read TMP[2][2] previously!!!!
			MtrxAdd(B, B, TMP);

			// z = z+cross(W(:,i)*V(:,i));
			MtrxCross(tmp, w, v);
			MtrxAdd (z, z, tmp);
		}
		VectorFree(tmp);			
		VectorFree(v);
		VectorFree(w);

		// S   = B+B';
		// rho = trace(B);
		MtrxAdd(S, B, MtrxTranspose(TMP, B));
		rho = MtrxTrace(B);

		// Construct K
		K[1][1]=S[1][1]-rho; K[1][2]=S[1][2];     K[1][3]=S[1][3];     K[1][4]=z[1];
		K[2][1]=S[2][1];     K[2][2]=S[2][2]-rho; K[2][3]=S[2][3];     K[2][4]=z[2];
		K[3][1]=S[3][1];     K[3][2]=S[3][2];     K[3][3]=S[3][3]-rho; K[3][4]=z[3];
		K[4][1]=z[1];        K[4][2]=z[2];        K[4][3]=z[3];        K[4][4]=rho;

		// Find largest eigenvalue of K
		MtrxSvd (Us, Ss, Vs, K);

		// We known that lambda=3 allways, so the following
		// two lines are not required:
		// MtrxSvdSort (Us, Ss, Vs);
		// lambda = Ss[1][1];
		lambda = 3;
		
		// K = K-lambda*I;
		MtrxScl(I4, I4, lambda);
		MtrxSub(K, K, I4);

		// Perform svd of K-lambda*I to find eigenvector
		MtrxSvd(Us, Ss, Vs, K);

		MtrxSvdSort (Us, Ss, Vs);
		if (Ss[3][3] < MTRXEPSILON)
			index = 3;
		else
			index = 4;

		/******* OLDSTUFF *********
		// Find index of Ss where singular value 0 is
		lambda = 1000.0;		// Some large number
		index  = 0;
		for (i=1;i<=4;i++)	
		{
			if ( Ss[i][i] < lambda)
			{
				lambda = Ss[i][i];
				index = i;
			}
		}
		******** OLDSTUFF ********/

		// Extract column #index from Vs
		// (Could be that we should use opposite indexing)
		if (Vs[4][index] > 0)
		{
			q[1] = Vs[1][index];
			q[2] = Vs[2][index];
			q[3] = Vs[3][index];
			q[4] = Vs[4][index];
		}
		else
		{
			q[1] = -Vs[1][index];
			q[2] = -Vs[2][index];
			q[3] = -Vs[3][index];
			q[4] = -Vs[4][index];
		}

//		q[1] = Vs[index][1];
//		q[2] = Vs[index][2];
//		q[3] = Vs[index][3];
//		q[4] = Vs[index][4];

		// Normalize q
		MtrxNormalize2 (q, q);
	

		// Free space
		MtrxFree(I4);
		VectorFree(z);
		MtrxFree(TMP);
		MtrxFree(Vs);
		MtrxFree(Ss);
		MtrxFree(Us);
		MtrxFree(K);
		MtrxFree(S);
		MtrxFree(B);
	}
	return q;
}
