/*************************************************************************
 *
 * FILE		: $Archive: $ 
 * DATE		: 1999.02.04
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT	: This file contains Matlab-like functions for matrix
 *			  operations. Allthough, the input/output parameters are
 *			  somewhat different, the function names are kept identical
 *			  apart from that the libraries exported functions are 
 *			  capitalised.
 *            
 * NOTES	: 
 *
 * AUTHOR(s): Karl-Petter Lindegaard
 *
 ************************************************************************/ 

#include "mtrx.h"
#include <math.h>



/*********************************************************************

  Internal functions

*********************************************************************/


#define min(a, b)  (((a) < (b)) ? (a) : (b))
#define max(a, b)  (((a) > (b)) ? (a) : (b))


VARTYPE sign(VARTYPE a, VARTYPE b) 
/* returns Sign(b)*Abs(a) */
{
	if (b >= (VARTYPE)0.0)
		return fabs(a);
	else 
		return (-fabs(a));
}

	
VARTYPE pythag(VARTYPE a, VARTYPE b)
/* Computes sqrt(a^2 + b^2) without destructive underflow or overflow */
{
	VARTYPE at = fabs(a);
	VARTYPE bt = fabs(b);

	if (at > bt)
		return (at*sqrt(1.0+(bt*bt)/(at*at)));
	else if (bt == 0.0)
		return 0.0;
	else 
		return (bt*sqrt(1.0+(at*at)/(bt*bt)));
}




/*********************************************************************

  Exported Matlab-like functions

*********************************************************************/


/************************************************************************
 *
 * FUNCTION		: Eye
 * DATE			: 12.30.98
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT		: Return an identity matrix of a particular dimension
 *				  
 * PARAMETERS	:
 *	size		-I	int			Dimension n
 *
 * RETURNS		:
 *				-O	CMatrix		Identity matrix I (nxn)
 *
 * NOTES		:
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
matrix Eye(int size)
{
	int i;

	if (size > 0)
	{
		matrix A(size,size);
		for (i=1;i<=size;i++)
			A(i,i) = 1;
		return A;
	}
	matrix A;		// Empty matrix
	return A;
}


/************************************************************************
 *
 * FUNCTION		: Ones
 * DATE			: 12.30.98
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT		: Return a CMatrix filled with 1
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
matrix Ones(int rows, int cols)
{
	int i,j;

	if ( (rows > 0) && (cols > 0) )
	{
		matrix A(rows,cols);
		for (i=1;i<=rows;i++)
			for (j=1;j<=cols;j++)
				A(i,j) = 1;
		return A;
	}
	matrix A;		// Empty matrix
	return A;
}


/************************************************************************
 *
 * FUNCTION		: Zeros
 * DATE			: 12.30.98
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT		: Return a CMatrix filled with zeros
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
matrix Zeros(int rows, int cols)
{
	int i,j;

	if ( (rows > 0) && (cols > 0) )
	{
		matrix A(rows,cols);
		for (i=1;i<=rows;i++)
			for (j=1;j<=cols;j++)
				A(i,j) = 0;
		return A;
	}
	matrix A;		// Empty matrix
	return A;
}


/************************************************************************
 *
 * FUNCTION		: Tranpose
 * DATE			: 12.30.98
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT		: Return the transpose of a matrix. 
 *				  
 * PARAMETERS	:
 *  A			-I  CMatrix		Matrix to make the transpose of
 *
 * RETURNS		:
 *				-O	CMatrix		The transpose of A.
 *
 * NOTES		:
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
matrix Transpose(const matrix& A)
{
	return A.Transp();
}


/************************************************************************
 *
 * FUNCTION		: Det
 * DATE			: 12.30.98
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT		: Return the determinant of a real matrix. 
 *				  
 * PARAMETERS	:
 *  A			-I  CMatrix		Matrix to take the determinant of
 *
 * RETURNS		:
 *				-O	VARTYPE		The determinant
 *
 * NOTES		:
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
VARTYPE Det(matrix A)
{
	// Not implemented yet: Throw a general exception
	throw EMatrixException();

	return 0.0;
}


/************************************************************************
 *
 * FUNCTION		: Svd
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
 *	A			-I	CMatrix		Matrix to decompose
 *	U			-O	CMatrix&	First output matrix (m x n)
 *	S			-O	CMatrix&	Diagonal matrix with singular values
 *	V			-O	CMatrix&	Last output matrix.
 *
 * RETURNS		:
 *				-O	CMatrix		A vector (CMatrix of size nx1) containing
 *								singular values.
 *
 * NOTES		:
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
matrix Svd (const matrix& A, matrix& U, matrix &S, matrix &V)
{
	// Confirm non-empty input matrix.

	if (A.Isempty())
		throw EMatrixEmpty();


	// Local variables

	int m = A.Rows();
	int n = A.Cols();
	
	int nm, l, k, j, jj, its, i, flag;
	VARTYPE c, f, g, h, s, x, y, z, anorm, scale;
	int mnmin = min(m,n);
	
	matrix rv1 = Zeros(n,1);
	matrix W   = Zeros(n,1);
	matrix W2(mnmin,1);
	
	
	// Initalise output matrices
	
	U = A;
	S = Zeros(n,n);
	V = Zeros(n,n);


	// Start SVD decomposition
	
	g     = 0;
	scale = 0;
	anorm = 0;
	

	// Householder reduction to bidiagonal form
	for (i=1; i<=n; i++)
	{
		l = i+1;
		rv1(i) = scale*g;
		g     = 0;
		s     = 0;
		scale = 0;
		if (i<=m)
		{
			for (k=i;k<=m;k++)
				scale += fabs(U(k,i));
			if (scale != 0.0)
			{
				for (k=i;k<=m;k++)
				{
					U(k,i) /= scale;
					s += U(k,i)*U(k,i);
				}
				f = U(i,i);
				g = -sign(sqrt(s),f);
				h = f*g - s;
				U(i,i) = f-g;
				for (j=l; j<=n; j++)
				{
					s = 0;
					for (k=i;k<=m;k++)
						s += U(k,i)*U(k,j);
					f = s/h;
					for (k=i;k<=m;k++)
						U(k,j) += f*U(k,i);
				}
				for (k=i;k<=m;k++)
					U(k,i) *= scale;
			}
		}
		W(i) = scale*g;
		g     = 0;
		s     = 0;
		scale = 0;
		
		if ((i<=m) && (i != n))
		{
			for (k=l;k<=n;k++)
				scale += fabs(U(i,k));
			if (scale != 0.0)
			{
				for (k=l;k<=n;k++)
				{
					U(i,k) /= scale;
					s += U(i,k)*U(i,k);
				}
				f = U(i,l);
				g = -sign(sqrt(s),f);
				h = f*g - s;
				U(i,l) = f-g;
				for (k=l;k<=n;k++)
					rv1(k) = U(i,k) / h;
				for (j=l; j<=m; j++)
				{
					s = 0;
					for (k=l;k<=n;k++)
						s += U(j,k)*U(i,k);
					for (k=l;k<=n;k++)
						U(j,k) += s*rv1(k);
				}
				for (k=l;k<=n;k++)
					U(i,k) *= scale;
			}
		}
		anorm = max (anorm, (fabs(W(i)) + fabs(rv1(i))) );
	}
	
	
	// Accumulation of right-hand transformations
	for (i=n;i>=1;i--)
	{
		if (i<n)
		{
			if (g!=0)
			{
				for (j=l;j<=n;j++)
					V(j,i) = (U(i,j)/U(i,l)) / g;
				for (j=l;j<=n;j++)
				{
					s = 0;
					for (k=l;k<=n;k++)
						s += U(i,k)*V(k,j);
					for (k=l;k<=n;k++)
						V(k,j) += s*V(k,i);
				}
			}
			for (j=l;j<=n;j++)
			{
				V(i,j) = 0;
				V(j,i) = 0;
			}
		}
		V(i,i) = 1.0;
		g = rv1(i);
		l = i;
	}
	
	// mnmin = min(m,n);
	
	
	// Accumulation of left-hand transformations
	for (i=mnmin;i>=1;i--)
	{
		l = i + 1;
		g = W(i);
		for (j=l;j<=n;j++)
			U(i,j) = 0;
		if (g != 0.0)
		{
			g = 1.0/g;
			for (j=l;j<=n;j++)
			{
				s = 0;
				for (k=l;k<=m;k++)
					s += U(k,i)*U(k,j);
				f = (s/U(i,i))*g;
				for (k=i;k<=m;k++)
					U(k,j) += f*U(k,i);
			}
			for (j=i;j<=m;j++)
				U(j,i) *= g;
		}
		else
			for (j=i;j<=m;j++)
				U(j,i) = 0;
			U(i,i) += 1;
	}
	
	
	// Diagonalization of bidiagonal form: Loop over singular values
	// and over allowed iterations.
	
	for (k=n;k>=1;k--)
	{
		for (its=1;its<=MAXITERATIONS;its++)
		{
			flag = 1;
			for (l=k;l>=1;l--)			// Test for splitting
			{
				nm = l-1;				// Note that rv1(1) is allways zeros
				if ( (fabs(rv1(l))+anorm) == anorm)
				{
					flag = 0;
					break;
				}
				if ( (fabs(W(nm))+anorm) == anorm)
					break;				// Goto step 1
			}
			
			/*Step 1:*/
			if (flag)
			{
				c = 0;
				s = 0;
				for (i=1;i<=k;i++)
				{
					f = s*rv1(i);
					rv1(i) = c*rv1(i);
					if ( (fabs(f)+anorm) == anorm)
						break;			// Goto step 2
					g = W(i);
					h = pythag(f,g);
					W(i) = h;
					h = 1/h;
					c = g*h;
					s = -(f*h);
					for (j=1;j<=m;j++)
					{
						y = U(j,nm);
						z = U(j,i);
						U(j,nm) = (y*c) + (z*s);
						U(j,i)  = -(y*s) + (z*c);
					}
				}
			} // if (flag)
			z = W(k);
			if (l == k)				// Convergence
			{
				if (z < 0.0)		// Sing. values are made nonnegative
				{
					W(k) = -z;
					for (j=1;j<=n;j++)
						V(j,k) = -V(j,k);
				}
				break;
			}
			if (its == MAXITERATIONS)
				// No convergence in x SVDCMP iterations
				throw EMatrixIteration();
			
			x  = W(l);				// Shift from bottom 2 by 2 minor
			nm = k - 1;
			y  = W(nm);
			g  = rv1(nm);
			h  = rv1(k);
			f  = ((y - z)*(y + z) + (g-h)*(g+h)) / (2.0*h*y);
			g  = pythag(f, 1.0);
			f  = ((x - z)*(x + z) + h*((y/(f + sign(g,f))) - h))/x;
			c  = 1.0;
			s  = 1.0;
			for (j=l;j<=nm;j++)
			{
				i  = j + 1;
				g  = rv1(i);
				y  = W(i);
				h  = s*g;
				g  = c*g;
				z  = pythag(f,h);
				rv1(j) = z;
				c  = f/z;
				s  = h/z;
				f  = (x*c) + (g*s);
				g  = -(x*s) + (g*c);
				h  = y*s;
				y  *= c;
				for (jj=l;jj<=n;jj++)
				{
					x = V(jj,j);
					z = V(jj,i);
					V(jj,j) = (x*c) + (z*s);
					V(jj,i) = -(x*s) + (z*c);
				}
				z = pythag(f,h);
				W(j) = z;
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
					y = U(jj,j);
					z = U(jj,i);
					U(jj,j) = (y*c) + (z*s);
					U(jj,i) = -(y*s) + (z*c);
				}
			}
			rv1(l) = 0.0;
			rv1(k) = f;
			W(k)   = x;
		}
	}

	// Generate return vector and matrix with sing. val along diagonal
	for (i=1;i<=mnmin;i++)
		S(i,i) = W2(i) = W(i);
		
	return W2;	// Return vector of singular values.
}


/************************************************************************
 *
 * FUNCTION		: Inv
 * DATE			: 12.30.98
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT		: Compute the inverse of a matrix.
 *
 * PARAMETERS	:
 *	A			-I	matrix		Matrix to compute inverse of.
 *
 * RETURNS		:
 *				-O	matrix		The inverse.
 *
 * NOTES		:
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
matrix Inv(const matrix& A)
{
	int i,n;

	if (A.Issquare())				// Verify square matrix
	{
		matrix U, S, V;

		try
		{
			Svd(A, U, S, V);
		}
		catch (EMatrixException Ex)
		{
			throw;							// Reraise the exception......
		}

		// Check for singularity
		n = S.Rows();
		for (i=1;i<=n;i++)
		{
			if ( fabs(S(i,i)) < EPSILON)
				throw EMatrixSingular();
			else
				S(i,i) = 1.0/S(i,i);		//  Make reciprocal
		}

		matrix B = V * S * U.Transp();		// Calculate inverse matrix
		return B;
	}
	else 
		throw EMatrixNonsquare();

	matrix dummy(0);		// Never reached
	return dummy;
}


/************************************************************************
 *
 * FUNCTION		: Balance
 * DATE			: 03.11.99
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT		: Given a matrix, this routine produces a balanced
 *				  matrix with identical eigenvalues. A symmetric is 
 *				  allready balanced and is unaffected by this procedure.
 *				  The constant radix should be the machine's floating 
 *				  point radix.
 *
 * PARAMETERS	:
 *	A			-I	matrix		Matrix to produce balanced matrix of.
 *
 * RETURNS		:
 *				-O	matrix		The balanced matrix.
 *
 * NOTES		:
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
matrix Balance(const matrix& A)
{
	//local constants and variables
	int last, i, j;
	VARTYPE s, r, g, f, c, sqrdx;
	VARTYPE radix = 2.0;				// TEMPORARY

	if ( !A.Issquare() )
		throw EMatrixNonsquare();

	int n = A.Rows();

	/* Make a result matrix */
	matrix B = A;	

	/* start main */
	sqrdx = radix*radix;
	last  = 0;
	while (last==0) {
		last = 1;

		/*Calculate row and column norms*/
		for(i=1;i<=n;i++){					
			c = r = 0.0;
			for(j=1;j<=n;j++)
				if (j != i){
					c += fabs(B(j,i));
					r += fabs(B(i,j));
				}
			if((c != 0.0) && (r != 0.0)){
				g = r/radix;
				f = 1.0;
				s = c + r;
				while (c < g){				/* Find the integer power of the machine */
					f *= radix;				/* radix that comes closest to balancing the */
					c *= sqrdx;				/* matrix */
				}
				g = r*radix;
				while (c > g){
					f /= radix;
					c /= sqrdx;
				}
				if ((c+r)/f < 0.95*s){
					last = 0;
					g = 1.0/f;
					/* Apply similarity transformation */
					for(j=1;j<=n;j++)		
						B(i,j) *= g;
					for(j=1;j<=n;j++)
						B(j,i) *= f;
				}
			}
		}
	}
	return B;
}


/************************************************************************
 *
 * FUNCTION		: ElmHess
 * DATE			: 03.11.99
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT		: Reduction to Hessenberg form by the
 *				  elimination method. Recommended, but not required
 *				  is that this routine is preceeded by Balance.
 *
 * PARAMETERS	:
 *	A			-I	matrix		Matrix to produce an upper Hessenberg matrix
 *								with identical eigenvalues from. 
 *
 * RETURNS		:
 *				-O	matrix		The Hessenberg matrix. (On output
 *								from, the Hessenberg matrix is in elements
 *								B[i,j] with i < j+2, while elements with
 *								i > j+1 are returned with random values to
 *								be thought of as zero)
 *
 * NOTES		:
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
matrix ElmHess(const matrix& A)
{
	/* local variables */
	int m, j, i;
	VARTYPE x, y;

	if ( !A.Issquare() )
		throw EMatrixNonsquare();

	int n = A.Rows();

	/* Make a result matrix */
	matrix B = A;

	/* start main ***************************************************/
 	for(m=2;m<n;m++){
		x = 0.0;
		i = m; 
		/* find the pivot */
		for(j=m;j<=n;j++){						
			if (fabs(B(j,m-1)) > fabs(x)){
				x = B(j,m-1);
				i = j; 
			}
		}
		/*interchange the rows */
		if (i != m){							
			for(j=m-1;j<=n;j++){
				y = B(i,j);
				B(i,j) = B(m,j);
				B(m,j) = y;
			}
			for(j=1;j<=n;j++){
				y = B(j,i);
				B(j,i) = B(j,m);
				B(j,m) = y;
			}
		}
		/* carry out the elimination */
		if ( x!= 0.0) {							
			for(i=m+1;i<=n;i++){
				y = B(i,m-1);
				if (y != 0.0) {
					y /= x;
					B(i,m-1) = y;
					for(j=m;j<=n;j++)
						B(i,j) -= y*B(m,j);
					for(j=1;j<=n;j++)
						B(j,m) += y*B(j,i);
				}
			}
		}
	}

	return B;
}


/************************************************************************
 *
 * FUNCTION		: Hess
 * DATE			: 03.11.99
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT		: Reduction to Hessenberg form by the
 *				  elimination method. Recommended, but not required
 *				  is that this routine is preceeded by Balance.
 *
 * PARAMETERS	:
 *	A			-I	matrix		Matrix to produce an upper Hessenberg matrix
 *								with identical eigenvalues from. 
 *
 * RETURNS		:
 *				-O	matrix		The Hessenberg matrix.
 *
 * NOTES		:
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
matrix Hess(const matrix& A)
{
	/* local variables */
	int j, i;

	if ( !A.Issquare() )
		throw EMatrixNonsquare();

	int n = A.Rows();


	/* Call ElmHess() and zero out the entries below the subdiagonal */

	matrix B = ElmHess(A);	/* ElmHess() returns random values below sub-diag */
	for (i=3;i<=n;i++)
	{
		for (j=1;j<=i-2;j++)
			B(i,j) = 0;
	}

	return B;
}



/************************************************************************
 *
 * FUNCTION		: HQR
 * DATE			: 03.11.99
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT		: Find all eigenvalues of an upper Hessenberg matrix.
 *
 * PARAMETERS	:
 *	A			-I	matrix		Upper Hessenberg matrix to calcualte 
 *								eigenvalues from. 
 *
 * RETURNS		:
 *				-O	bool		True if eigenvalues are complex.
 *
 * NOTES		:
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
bool HQR(const matrix& A, matrix& Wr, matrix& Wi)
{
	/* Local variables  */
	int nn, m, l, k, j, its, i, mmin;
	VARTYPE z, y, x, w, v, u, t, s, r, q, p, anorm;
	bool complexeigs = false;

	if ( !A.Issquare() )
		throw EMatrixNonsquare();

	int n = A.Rows();
	Wr = Zeros(n,1);		/* Allocate result vectors. At the same time */
	Wi = Zeros(n,1);		/* their previous contents will be lost.     */


	/* start main *********************************************/
	
	/* Compute matrix norm for possible use in locating singl small subdiagonal
	elemnt */
	anorm = fabs(A(1,1));
	for(i=2;i<=n;i++)
		for(j=(i-1);j<=n;j++)
			anorm +=fabs(A(i,j));

	nn = n;
	t =0.0;
	
	/* Begin search for next eigenvalue */
	while (nn >= 1){
		its = 0;
		do{
			for(l=nn;l>=2;l--){							/*Begin iterations: look for */
				s = fabs(A(l-1,l-1)) + fabs(A(l,l));	/*single small subdiagonal element*/
				if (s == 0.0)
					s = anorm;
				if ((VARTYPE)(fabs(A(l,l-1)) + s) == s)
					break;
			}
			x = A(nn,nn);
			if (l == nn) {								/* One root found. */
				Wr(nn) = x+t;
				Wi(nn--)= 0.0;
			}
			else {
				y = A(nn-1,nn-1);
				w = A(nn,nn-1)*A(nn-1,nn);
				if (l == (nn-1)){						/* Two roots found..... */
					p = 0.5*(y-x);
					q = p*p +w;
					z = sqrt(fabs(q));
					x += t;
					if (q >= 0.0){						/* A real pair */ 
						z = p + sign(z,p);
						Wr(nn-1) = Wr(nn) = x + z;
						if (z)
							Wr(nn) = x - w/z;
						Wi(nn-1) = Wi(nn) = 0.0;
					}
					else{								/* A complex pair */
						Wr(nn-1) = Wr(nn) = x + p;
						Wi(nn-1) = -(Wi(nn) = z);
						complexeigs = true;
					}
					nn -= 2;
				}
				else {									/* No roots found. Continue */
					if (its == MAXITERATIONS) {			/* iterations */
						throw EMatrixIteration();
					}
					if ((its == 10)|| (its == 20)){		/* Form exceptional shift */
						t +=x;
						for(i=1;i<=nn;i++)
							A(i,i) -= x;
						s = fabs(A(nn,nn-1)) + fabs(A(nn-1,nn-2));
						y = x = 0.75*s;
						w = - 0.4375*s*s;
					}
					++its;
					for(m=(nn-2);m>=1;m--) {			/* Form shift and then look for */
						z = A(m,m);					/* 2 consecutive subdiagonal */
						r = x - z;						/* elements */
						s = y - z;
						p = (r*s - w)/A(m+1,m) + A(m,m+1); 
						q = A(m+1,m+1) - z - r - s;
						r = A(m+2,m+1);
						s = fabs(p) + fabs(q) + fabs(r);/* Scale to prevent underflow or */
						p /= s;							/* overflow */
						q /= s;
						r /= s;
						if (m == 1) 
							break;
						u = fabs(A(m,m-1))*(fabs(q) + fabs(r));
						v = fabs(p)*(fabs(A(m-1,m-1)) + fabs(z) + fabs(A(m+1,m+1)));
						if((double)(u+v) == v)
							break;
					}
					for(i=m+2;i<=nn;i++) {
						A(i,i-2) = 0.0;
						if (i != (m+2))
							A(i,i-3) = 0.0;
					}
					
					/* Double QR steps on rows 1 to nn and columns m to nn */
					for(k=m;k<=nn-1;k++) {
						if (k != m) {					/* Begin setup of householder */ 
							p = A(k,k-1);				/* vector */
							q =A(k+1,k-1);
							r =0.0;
							if (k != (nn-1)) 
								r = A(k+2,k-1);
							if ((x = fabs(p) + fabs(q) + fabs(r)) != 0.0) {
								p /= x;					/* Scale to prevent underflow */
								q /= x;					/* or overflow */
								r /= x;
							}
						}
						if ((s = sign(sqrt(p*p+q*q+r*r), p)) != 0.0) {
							if (k == m) {
								if (l != m)
									A(k,k-1) = - A(k,k-1);
							}
							else 
								A(k,k-1) =  -s*x;
							p += s;
							x = p/s;
							y = q/s;
							z = r/s;
							q /= p;
							r /= p;
							for(j=k;j<=nn;j++) {		/* Row modification */			
								p = A(k,j) + q*A(k+1,j);
								if (k != (nn-1)) {
									p += r*A(k+2,j);
									A(k+2,j) -= p*z;
								}
								A(k+1,j) -= p*y;
								A(k,j) -= p*x;
							}
							mmin = nn<k+3 ? nn: k+3;
							for(i=l;i<=mmin;i++) {		/* Column modification */
								p = x*A(i,k) + y*A(i,k+1);
								if (k != (nn-1)) {
									p += z*A(i,k+2);
									A(i,k+2) -= p*r;
								}
								A(i,k+1) -= p*q;
								A(i,k) -= p;
							}
						}
					}
				}
			}
		} while (l < nn-1);
	}
	return complexeigs;
}



/************************************************************************
 *
 * FUNCTION		: Eig
 * DATE			: 03.29.99
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT		: Find all eigenvalues of a matrix.
 *
 * PARAMETERS	:
 *	A			-I	matrix		Matrix to calcualte eigenvalues from. 
 *
 * RETURNS		:
 *				-O	matrix		If real eigenvalues only, the returned 
 *								matrix is in fact a nx1 column vector.
 *								If comples eigenvalues, the imaginary parts
 *								are stored in the second column, thus a nx2 
 *								matrix is returned.
 *
 * NOTES		:
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
matrix Eig(const matrix& A)
{
	int i;

	if ( !A.Issquare() )
		throw EMatrixNonsquare();

	bool comp;
	matrix Wr, Wi;
	int n = A.Rows();


	matrix B = Balance(A);		/* Balancing of matrix */
	matrix C = Hess(B);			/* Reduction to Hessenberg form */
	comp = HQR(C, Wr, Wi);		/* Calculate eigenvalues */

	matrix E;					/* Empty matrix */
	if (comp == true)
	{
		E = Zeros(Wr.Rows(), 2);
		for (i=1;i<=n;i++)		/* Copy imaginary parts */
			E(i,2) = Wi(i);
	}
	else
		E = Zeros(Wr.Rows(), 1);

	for (i=1;i<=n;i++)			/* Copy real parts */
		E(i,1) = Wr(i);

	return E;
}	





/************************************************************************
 *
 * FUNCTION		: Interp1
 * DATE			: 99.03.04
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT		: 1 dimentional linear table interpolation.
 *
 *				  YI = INTERP1(x,Y,xi) interpolates to find YI, the values of
 *				  the underlying function Y at the points in the vector xi.
 *				  The vector x specifies the points at which the data Y is
 *				  given. If Y is a matrix, then the interpolation is performed
 *				  for each column of Y and YI will be length(xi)-by-size(Y,2).
 *				  Out of range indexes in xi are returned as Y's boundary values.  
 *
 * PARAMETERS	:
 *	X			-I	matrix		Points at which Y is given (vector).
 *  Y			-I  matrix		Matrix of samples.
 *  XI			-I  matrix		Vector of points to interpolate for.
 *
 * RETURNS		:
 *  			-O	matrix		An mxi x ny matrix of interpolated values.
 *								where mxi = rows of XI, ny = cols of Y.
 *
 * NOTES		: Matlab interp1 - function
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
matrix Interp1(const matrix& X, const matrix& Y, const matrix& XI)
{
	double	xhigh, xlow, yhigh, ylow, xval, alpha;
    int  i, j, k;

	int mxi = XI.Length();
	int mx  = X.Length();
	int my  = Y.Rows();
	int ny  = Y.Cols();


	/* Check sizes of input arguments */

	if (mx != my)
		throw EMatrixCorrespond();


	/* Allocate space for result matrix */

	matrix YI(mxi,ny);


	/* for each element in XI */

	for (k=1;k<=mxi;k++)
    {  
	    xval = XI(k);

		/* find first element in tab which is larger than xval */

		if (xval <= X(1))
		{
			/* xval < smallest value to interpolate */
			for (i=1;i<=ny;i++)
				YI(k,i) = Y(1,i);
		}

		else if (xval >= X(my))
		{
			/* xval >= largest value to interpolate */	
			for (i=1;i<=ny;i++)
				YI(k,i) = Y(my,i);
		}
		else
		{

			/* xval in interval <X(1) X(my)> */

			for (j=1;(j<=my) && (X(j) <= xval);j++);
			xhigh = X(j);
			xlow  = X(j-1);

			for (i=1;i<=ny;i++)
			{
				yhigh = Y(j,i);
				ylow  = Y(j-1,i);
				alpha = (xhigh-xval)/(xhigh-xlow);
				YI(k,i) = alpha*ylow + (1.0-alpha)*yhigh;
			} /* for i...*/
	   		
		} /* else */

	} /* for k ....*/
	
	return YI;		
}



/************************************************************************
 *
 * FUNCTION		: Interp2
 * DATE			: 99.03.04
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT		: 2 dimentional linear table interpolation.
 *
 *				  ZI = INTERP2(X,Y,Z,XI,YI) interpolates to find ZI, the values of the
 *				  underlying 2-D function Z at the points in matrices XI and YI.
 *				  Matrices X and Y specify the points at which the data Z is given.
 *				  Out of range indexes in xi and yi are returned as Z's boundary values. 
 *
 * PARAMETERS	:
 *	X			-I	matrix		Row values at which Z is given (vector).
 *  Y			-I	matrix		Column values at which Z is given (vector)
 *  Z			-I	matrix		Matrix of samples
 *  XI			-I	matrix		Vector of row values to interpolate for
 *  YI			-I	matrix		Vector of column values to interpolate for
 *
 * RETURNS		:
 *  			-O	matrix		An myi x mxi matrix of interpolated values.
 *								where mxi = length of XI, myi = length of Y.
 *
 * NOTES		: Matlab interp2 - function
 *
 * AUTHOR(s)	: Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
matrix Interp2(const matrix& X, const matrix& Y, const matrix& Z, const matrix& XI, const matrix& YI)
{
	/* Pseudo-code -------------------------
	      1. Run Interp1 wrt yi
		  2. Transpose
		  3. Run Interp1 wrt xi
		  4. Transpose
		  5. Return result from 4.
	---------------------------------------*/

	/* Step 1: */
	matrix STEP1 = Interp1(Y, Z, YI);

	/* Step 2: */
	matrix STEP1T = STEP1.Transp();

	/* Step 3: */
	matrix STEP2 = Interp1(X, STEP1T, XI);

	/* Step 4 and 5: */
	return STEP2.Transp();
}
