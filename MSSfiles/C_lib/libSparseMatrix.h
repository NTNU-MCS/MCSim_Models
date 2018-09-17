/* 
    libSparseMatrix.h

	Header file for the sparse matrix functions.

	References:
	[1]	E. Horowitz, S. Sahni, S. Anderson-Freed, Fundamentals of Data Structures in C, 
		Computer Science Press, New York, 1993, ISBN 0-7167-8250-2.
	[2]	W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery, Numerical Recipes in C - The Art of Scientific Computing, 2nd ed. 
		Cambridge University Press, New York, USA, 1992.

	Copyright: 		Roger Skjetne, NTNU
    Date created: 	2011.03.31 Roger Skjetne
    Revised:      	

*/


#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "simstruc.h"


/* Constants: */
#define MAX_ELEMENTS	300		/* Maxnumber of nonzero elements in matrices */
#define MAX_STRINGCHAR	256		/* Max number of characters string arrays	*/


typedef struct {	
/* 
	Defining a sparse matrix as a structure containing triples (row,col,value).
	The [0] element must contain number of rows, columns, and nonzero elements.
	Zero-indexing is used for the matrix to be represented. 
	Ex: A = [0 10; 20 30] gives 
		A[0] = (2,2,3), A[1] = (0,1,10), A[2] = (1,0,20), A[3] = (1,1,30) 
*/
	int		row;	// element row index (zero-indexing)
	int		col;	// element col index (zero-indexing)
	float	value;	// element value
} sSparseMtrx;


void  transpose(SimStruct *S, sSparseMtrx A[], sSparseMtrx B[]);
void  mmult(SimStruct *S, sSparseMtrx A[], sSparseMtrx B[], sSparseMtrx C[]);
void  madd(SimStruct *S, sSparseMtrx A[], sSparseMtrx B[], sSparseMtrx C[]);
void  msub(SimStruct *S, sSparseMtrx A[], sSparseMtrx B[], sSparseMtrx C[]);
double rowsum(SimStruct *S, sSparseMtrx A[], int Idx);
void  zerorow(SimStruct *S, sSparseMtrx A[], int Idx);
void  zeroelement(SimStruct *S, sSparseMtrx A[], int row, int col);
void  updateelement(SimStruct *S, sSparseMtrx A[], int row, int col, float value);
float readelementvalue(SimStruct *S, sSparseMtrx A[], int row, int col);
int   ifnonzero(SimStruct *S, sSparseMtrx A[], int row, int col);
int   compare(SimStruct *S, int x, int y);
void  storesum(SimStruct *S, sSparseMtrx C[], int *totalC, int row, int column, float *sum, int *flag);

void   sLUdcmp(SimStruct *S, double **a, int n, int *indx, double *d);
void   sLUbksb(double **a, int n, int *indx, double b[]);


#endif
