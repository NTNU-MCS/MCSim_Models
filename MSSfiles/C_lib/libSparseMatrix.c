/* 
    libSparseMatrix.c

	Library file for sparse matrix functions.

	References:
	[1]	E. Horowitz, S. Sahni, S. Anderson-Freed, Fundamentals of Data Structures in C, 
		Computer Science Press, New York, 1993, ISBN 0-7167-8250-2.
	[2]	W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery, Numerical Recipes in C - The Art of Scientific Computing, 2nd ed. 
		Cambridge University Press, New York, USA, 1992.

	Copyright: 		Roger Skjetne, NTNU
    Date created: 	2011.03.31 Roger Skjetne
    Revised:      	

*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* #################################  SPARSE MATRIX FUNCTIONS	############################### */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

#include "libSparseMatrix.h"
#include "libGenericFunctions.h"

void transpose(SimStruct *S, sSparseMtrx A[], sSparseMtrx B[])
/* B is set to the transpose of A */
/* Note that zero-indexing is used in matrices.*/
{
	int n,i,j, currentB;
	n = (int)A[0].value;			// total number of elements
	B[0].row = A[0].col;	// rows in B = columns in A
	B[0].col = A[0].row;	// columns in B = rows in A
	B[0].value = (float)n;
	if(n>0) {	//non-zero matrix
		currentB = 1;
		for(i=0; i<A[0].col; i++)
		/* transpose by the columns in A */
			for(j=1; j<= n; j++)
			/* find elements from the current column */
				if(A[j].col == i) {
				/* elements is in current column, add it to B */
					B[currentB].row = A[j].col;
					B[currentB].col = A[j].row;
					B[currentB].value = A[j].value;
					currentB++;
				}
	}
}

void mmult(SimStruct *S, sSparseMtrx A[], sSparseMtrx B[], sSparseMtrx C[])
/* Multiply two sparse matrices C = A*B */
{
	int			i, j, column, flag = 0; 
	int			totalA = (int)A[0].value, totalB = (int)B[0].value; 
	int			totalC = 0;
	int			rows_A = A[0].row, cols_A = A[0].col;
	int			cols_B = B[0].col;
	int			row_begin = 1, row = A[1].row; 
	float		sum = 0.0;
	char		msg[MAX_STRINGCHAR+1];
	sSparseMtrx new_B[MAX_ELEMENTS+1];

	if(cols_A != B[0].row) {
		sprintf(msg,"Incompatible dimension of input matrices for Matrix multiplication.");
		SendMsg(S,ERROR_LOG,msg,"mmult");
	}

	transpose(S, B, new_B);
	
	/* Set boundary condition */
	if(totalA >= MAX_ELEMENTS-1) {
		sprintf(msg,"Out of elements in struct array sSparseMtrx A. Please increase MAX_ELEMENTS.");
		SendMsg(S,ERROR_LOG,msg,"mmult");
	}
	A[totalA+1].row = rows_A;
	new_B[totalB+1].row = cols_B;
	new_B[totalB+1].col = 0;
	new_B[totalB+1].value = 0.0;

	for(i=1; i<= totalA;){
		column = new_B[1].row;
		for(j=1; j<=totalB+1;) {
			/* multiply row of A by column of B */
			if(A[i].row != row) {
				storesum(S,C,&totalC,row,column,&sum,&flag);
				i = row_begin;
				for(; new_B[j].row == column; j++)
					;
				column = new_B[j].row;
			}
			else if(new_B[j].row != column) {
				storesum(S,C,&totalC,row,column,&sum,&flag);
				i = row_begin;
				column = new_B[j].row;
			}
			else switch(compare(S,A[i].col, new_B[j].col)) {
				case -1: /* go to next term in A */
					i++;
					break;
				case 0:  /* add terms, go to next term in  A and B */
					sum += (A[i++].value * new_B[j++].value);
					flag = 1;
					break;
				case 1:  /* advance to next term in B */
					j++;
					break;
			}
		} /* end of for j <= totalB+1 */
		for(; A[i].row == row; i++)
			;
		row_begin = i; row = A[i].row;
	} /* end of for i <= totalA */
	C[0].row	= rows_A;
	C[0].col	= cols_B;
	C[0].value	= (float)totalC;
}


void madd(SimStruct *S, sSparseMtrx A[], sSparseMtrx B[], sSparseMtrx C[])
/* Add two sparse matrices C = A + B */
{
	int		i = 1, j = 1, k, row_flag, col_flag;
	int		totalA = (int)A[0].value, totalB = (int)B[0].value; 
	int		totalC = totalA + totalB;
	char	msg[MAX_STRINGCHAR+1];


	if((A[0].row != B[0].row) || (A[0].col != B[0].col)) {
		sprintf(msg,"Incompatible dimension of input matrices for Matrix addition operation.");
		SendMsg(S,ERROR_LOG,msg,"madd");
	}
	else {
		/* Set boundary condition */
		if(totalA >= MAX_ELEMENTS-1) {
			sprintf(msg,"Out of elements in struct array sSparseMtrx A. Please increase MAX_ELEMENTS.");
			SendMsg(S,ERROR_LOG,msg,"madd");
		}
		if(totalB >= MAX_ELEMENTS-1) {
			sprintf(msg,"Out of elements in struct array sSparseMtrx B. Please increase MAX_ELEMENTS.");
			SendMsg(S,ERROR_LOG,msg,"madd");
		}
		A[totalA+1].row = A[0].row;
		B[totalB+1].row = B[0].row;

		for(k=1;k<=totalC;k++) {
			row_flag = compare(S, A[i].row, B[j].row);
			switch(row_flag) {
				case -1: /* A has elements in lowest row */
					C[k] = A[i];
					i++;
					break;
				case 0:  /* A and B has elements at the same starting row */
					col_flag = compare(S, A[i].col, B[j].col);
					switch(col_flag) {
						case -1: /* A has elements in lowest column */
							C[k] = A[i];
							i++;
							break;
						case 0:  /* A and B has elements at the same starting column */
							C[k] = A[i];
							C[k].value = A[i].value + B[j].value;
							i++; j++;
							totalC--;
							break;
						case 1:  /* B has elements in lowest column */
							C[k] = B[j];
							j++;
							break;
					}
					break;
				case 1:  /* B has elements in lowest row */
					C[k] = B[j];
					j++;
					break;
			}
		}
		C[0].row	= A[0].row;
		C[0].col	= A[0].col;
		C[0].value	= (float)totalC;
	}
}


void msub(SimStruct *S, sSparseMtrx A[], sSparseMtrx B[], sSparseMtrx C[])
/* Sparse matrix subtraction C = A - B */
{
	int			i;
	char		msg[MAX_STRINGCHAR+1];
	sSparseMtrx new_B[MAX_ELEMENTS+1];

	if((A[0].row != B[0].row) || (A[0].col != B[0].col)) {
		sprintf(msg,"Incompatible dimension of input matrices for Matrix subtraction operation.");
		SendMsg(S,ERROR_LOG,msg,"msub");
	}

	new_B[0].row	= B[0].row;
	new_B[0].col	= B[0].col;
	new_B[0].value	= B[0].value;

	for(i=1; i<=(int)B[0].value; i++) {
		new_B[i].row	= B[i].row;
		new_B[i].col	= B[i].col;
		new_B[i].value	= B[i].value*(-1);
	}
	madd(S, A, new_B, C);
}


double rowsum(SimStruct *S, sSparseMtrx A[], int Idx)
{
	int i, N = (int)A[0].value;
	double sum = 0.0;

	for(i=1;i<=N;i++) {
		if(A[i].row == Idx) 
			sum += (double)A[i].value;
		else if(A[i].row > Idx)
			break;
	}
	return sum;
}

void zerorow(SimStruct *S, sSparseMtrx A[], int Idx) 
{
	int i, N = (int)A[0].value, h1 = 0, h2 = -1;

	for(i=1;i<=N;i++) {
		if(A[i].row == Idx) {
			h2 = i;
			if(h1==0) h1 = i;
		}
		else if(A[i].row > Idx) {
			if(h2 != -1) A[i-(h2-h1+1)] = A[i];
			else break;
		}
	}
	A[0].value = A[0].value - (float)(h2-h1+1);
}


void zeroelement(SimStruct *S, sSparseMtrx A[], int row, int col) 
{
	int i, N = (int)A[0].value, h1 = N+1;

	for(i=1;i<=N;i++) {
		if((A[i].row == row) && (A[i].col == col)) 
			h1 = i;
		else if(i > h1) 
			A[i-1] = A[i];
	}
	if(h1 < N+1) A[0].value--;
}


void updateelement(SimStruct *S, sSparseMtrx A[], int row, int col, float value) 
{	// Function that updates an existing element or adds a new element to A
	int   i, N = (int)A[0].value, h1 = N+2, row0, row00, col0, col00;
	float val0, val00;
	char  msg[MAX_STRINGCHAR+1];

	if(row > A[0].row || col > A[0].col) {
		// New element is outside dimension of A.
		sprintf(msg,"Trying to add new element outside dimension of A.");
		SendMsg(S,ERROR_LOG,msg,"updateelement");
	}
	else if((row == A[N].row && col > A[N].col) || (row > A[N].row)) {
		// New element to be added after last existing element of A.
		A[N+1].row   = row;
		A[N+1].col   = col;
		A[N+1].value = value;
		A[0].value++;
	}
	else if(N == 0) {
		// First element to be added.
		A[1].row   = row;
		A[1].col   = col;
		A[1].value = value;
		A[0].value++;
	}
	else {
		// Existing element to be updated or new element to be inserted.
		for(i=1;i<=N+1;i++) {
			if((A[i].row == row) && (A[i].col == col) && (h1>i)) {
				A[i].value = value; // Existing element updated.
				break;
			}
			else if((A[i].row == row) && (A[i].col > col) && (h1>i)) {
				row0		= A[i].row;
				col0		= A[i].col;
				val0		= A[i].value;
				A[i].row	= row;
				A[i].col	= col;
				A[i].value	= value;
				h1 = i;
			}
			else if((A[i].row > row) && (h1>i)) {
				row0		= A[i].row;
				col0		= A[i].col;
				val0		= A[i].value;
				A[i].row	= row;
				A[i].col	= col;
				A[i].value	= value;
				h1 = i;
			}
			else if(i > h1) {
				row00		= A[i].row;
				col00		= A[i].col;
				val00		= A[i].value;
				A[i].row	= row0;		row0 = row00;
				A[i].col	= col0;		col0 = col00;
				A[i].value	= val0;		val0 = val00;
			}
		}
		if(h1 < N+2) A[0].value++;
	}
}


float readelementvalue(SimStruct *S, sSparseMtrx A[], int row, int col) 
{ // Function that reads out the value of an element in A
	int   i, N = (int)A[0].value;
	float val = 0.0;
	char  msg[MAX_STRINGCHAR+1];

	// Specified element is outside dimension of A.
	if(row > A[0].row || col > A[0].col) {
		sprintf(msg,"Trying to read element outside dimension of A.");
		SendMsg(S,ERROR_LOG,msg,"updateelement");
		return 0.0;
	}

	// Reading element.
	for(i=1;i<=N;i++) {
		if((A[i].row == row) && (A[i].col == col)) {
			val = A[i].value; // Existing element read.
			break;
		}
		else if(((A[i].row == row) && (A[i].col > col)) || (A[i].row > row)) 
			break;
	}

	return val;
}


int ifnonzero(SimStruct *S, sSparseMtrx A[], int row, int col) 
/* Function that returns 1 if the element (row,col) exist in the sparse array A[], 
	i.e. if the element is nonzero in the matrix A. If no match is found, it returns 0. */
{
	int i, k = 0, N = (int)A[0].value, h1 = N+1;

	for(i=1;i<=N;i++) {
		if((A[i].row == row) && (A[i].col == col)) { 
			k = 1; 
			break;
		}
		else if((A[i].row >= row) && (A[i].col > col))  
			break;
		else if(A[i].row > row)
			break;
	}
	return k;
}

int compare(SimStruct *S, int x, int y)
/* Compares x and y; returns -1 for less than, 0 for equal, and 1 for greater than*/
{
	if(x<y) return -1;
	else if (x==y) return 0;
	else return 1;
}

void storesum(SimStruct *S, sSparseMtrx C[], int *totalC, int row, int column, float *sum, int *flag)
/* if *flag != 0, then the sum along with its row and column position is stored as the *totalC+1 entry in C */
{
	char msg[MAX_STRINGCHAR+1];
	
	if(*flag)
		if(*totalC < MAX_ELEMENTS) {
			C[++*totalC].row = row;
			C[*totalC].col = column;
			C[*totalC].value = *sum;
			*sum = 0.0;
			*flag = 0;
		}
		else {
			sprintf(msg,"Numbers of terms in product exceeds %d.", MAX_ELEMENTS);
			SendMsg(S,ERROR_LOG,msg,"storesum");
		}
}






/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* ###############################  NUMERICAL RECIPES FUNCTIONS	############################### */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

#define NRANSI
#define TINY 1.0e-20
void	sLUdcmp(SimStruct *S, double **a, int n, int *indx, double *d)
/* Modified version of the LU decomposition function 'ludcmp' in [1]. */
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=dvector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) SendMsg(S,DEBUG_LOG,"Singular matrix in routine.","sLUdcmp");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_dvector(vv,1,n);
}

void	sLUbksb(double **a, int n, int *indx, double b[])
/* Modified version of the LU backsubstitution function 'lubksb' in [1]. */
{
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}

#undef TINY
#undef NRANSI


/*################################## END #############################################*/
