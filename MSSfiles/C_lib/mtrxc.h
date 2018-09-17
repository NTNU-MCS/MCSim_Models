#include "mtrxtype.h"
#include "mmalloc.h"

//////////////////////////////////////////////////////////////////
//                       Type definititions:                    //
//////////////////////////////////////////////////////////////////

typedef MTRXTYPE MtrxType;
typedef MTRXTYPE* Vector;
typedef MTRXTYPE** Matrix;


//////////////////////////////////////////////////////////////////
//                       Public Functions:                      //
//////////////////////////////////////////////////////////////////

/*
int MtrxGetM(Matrix x);
int MtrxGetN(Matrix x);
int MtrxGetVectorLen(Vector x);
*/
int MtrxGetM(void *x);
int MtrxGetN(void *x);


int MtrxIsMatrix(void *M); // non-zero if M is a matrix
int MtrxIsVector(void *M); // non-zero if M is a vector


//////////////////////////////////////////////////////////////////
//                   Allocation/Deallocation                    //
//////////////////////////////////////////////////////////////////

Vector VectorNew(int n);
Matrix MtrxNew(int m, int n);
void VectorFree(Vector v);
void MtrxFree(Matrix v);

Matrix MtrxEye(int size);
Matrix MtrxOnes(int rows, int cols);
Matrix MtrxZeros(int rows, int cols);


//////////////////////////////////////////////////////////////////
// Matrix manipulation functions                                //
//////////////////////////////////////////////////////////////////

Matrix MtrxSwapCols(Matrix Y, Matrix A, int i, int j);
Matrix MtrxSwapRows(Matrix Y, Matrix A, int i, int j);

//////////////////////////////////////////////////////////////////
// Elementary matrix functions                                  //
//////////////////////////////////////////////////////////////////

void *MtrxTranspose(Matrix Y, Matrix A);
void *MtrxAdd(void *Y, void *A, void *B);   // Y = A+B
void *MtrxSub(void *Y, void *A, void *B);   // Y = A-B
void *MtrxScl(void *Y, void *A, MtrxType c); // Y = c*A
void *MtrxMult(void *Y, Matrix A, void *B);  // Y = A*B
void MtrxSvd(Matrix U, Matrix S, Matrix V, Matrix A); // U*S*V' = A
void MtrxSvdSort(Matrix U, Matrix S, Matrix V);
Matrix MtrxInv(Matrix Y, Matrix A);

/*
void MtrxSub(Matrix Y, Matrix A, Matrix B);   // Y = A-B
void MtrxScl(Matrix Y, Matrix A, MtrxType c); // Y = c*A
void MtrxMult(Matrix Y, Matrix A, Matrix B);  // Y = A*B
void MtrxVec(Vector y, Matrix A, Vector b);   // y = A*b
*/

Vector MtrxCross(Vector y, Vector a, Vector b);
Matrix MtrxCrossS(Matrix S, Vector a);
MTRXTYPE MtrxNorm2(void *X);
void *MtrxNormalize2(void *Y, void *X);
MTRXTYPE MtrxTrace(Matrix X);

Matrix MtrxRotationEulerAngles(Matrix R, Vector Theta);
Matrix MtrxRotationQuaternions(Matrix R, Vector q);
Vector MtrxQuaternions2EulerAngles (Vector Theta, Vector q);
Vector MtrxQuest(Vector q, Matrix W, Matrix V);

int MtrxGaussJ(Matrix a, Matrix b);
