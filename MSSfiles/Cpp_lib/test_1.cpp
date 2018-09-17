/*************************************************************************
 *
 * FILE		: $Archive: $ 
 * DATE		: 1999.02.04
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT	: 
 *            
 * NOTES	: 
 *
 * AUTHOR(s): Karl-Petter Lindegaard
 *
 ************************************************************************/ 

/*************************************************************************
 * $History: $
 * 
 ************************************************************************/

#include "mtrx.h"
#include "mmalloc.h"


void Test1_Interp1()
{
	int i;

	matrix A(5,5);
	matrix X(6,1);

	for (i=1;i<=5;i++)
	{
		X(i) = (double)i;

		A(i,1) = 1000.0*i + 100;
		A(i,2) = 1000.0*i + 200;
		A(i,3) = 1000.0*i + 300;
		A(i,4) = 1000.0*i + 400;
		A(i,5) = 1000.0*i + 500;
	}
	cout << "A =" << endl << A;


	/* Select an input vector */

	matrix XI(3,1);
	XI(1) = 2.5;
	XI(2) = 1.1;
	XI(3) = 5.1;


	matrix XI2;

	/* Interpolate */

	matrix YI;

	try 
	{
		matrix YI = Interp1(X,A,XI2);
	}
	catch (EMatrixException Ex)
	{
		throw ;		// Reraise the exception
	}


	/* Write results to screen */

	cout << "XI = " << endl << XI;
	cout << "YI = " << endl << YI;
}


void Test1_Balance()
{
	matrix A(3,3);

	A(1,1) =  2.0;  A(1,2) = -4.0;  A(1,3) = 1.0;
	A(2,1) =  1.0;  A(2,2) = 10.0;  A(2,3) = 0.0;
	A(3,1) =  2.0;  A(3,2) =  1.0;  A(3,3) = 1.0;

	matrix B = Balance(A);
	matrix C = Hess(B);
	matrix D = Hess(A);

	cout << "Matrix A" << endl << A << endl;
	cout << "B = Balance(A)" << endl << B << endl;
	cout << "C = Hess(B)" << endl << C << endl;
	cout << "D = Hess(A)" << endl << D << endl;
}

void Test1_Eig()
{
	matrix A(3,3);

	A(1,1) =  2.0;  A(1,2) = -4.0;  A(1,3) = 1.0;
	A(2,1) =  1.0;  A(2,2) = 10.0;  A(2,3) = 0.0;
	A(3,1) =  2.0;  A(3,2) =  1.0;  A(3,3) = 1.0;

	matrix B = Eig(A);

	matrix C = A;
	C(3,1) = 0; C(3,2) = 0; C(3,3) = 0;
	matrix D = Eig(C);


	cout << "Matrix A" << endl << A << endl;
	cout << "B = Eig(A)" << endl << B << endl;
	cout << "Matrix C" << endl << C << endl;
	cout << "D = Eig(C)" << endl << D << endl;
}

void Test1_Performance_Mult()
{
	matrix A = Ones(100,100);
	matrix B;
	int i;
	int n;

	cout << "Enter number of 100x100 matrix multiplications (500 is a nice number): ";
	cin  >> n;
	cout << n << " multiplications in progress" << endl;

	for (i=1;i<=n;i++)
		B = A*A;

	cout << "Done! Type 0 and press ENTER";
	cin >> n;
}

void Test1_Performance_Svd()
{
	matrix A(6,6);
	matrix U;
	matrix S;
	matrix V;
	int i, n;

 	A(1,1)=35; A(1,2)=1;  A(1,3)=6;  A(1,4)=26; A(1,5)=19; A(1,6)=24;
    A(2,1)=3;  A(2,2)=32; A(2,3)=7;  A(2,4)=21; A(2,5)=23; A(2,6)=25;
    A(3,1)=31; A(3,2)=9;  A(3,3)=2;  A(3,4)=22; A(3,5)=27; A(3,6)=20;
    A(4,1)=8;  A(4,2)=28; A(4,3)=33; A(4,4)=17; A(4,5)=10; A(4,6)=15;
    A(5,1)=30; A(5,2)=5;  A(5,3)=34; A(5,4)=12; A(5,5)=14; A(5,6)=16;
    A(6,1)=4;  A(6,2)=36; A(6,3)=29; A(6,4)=13; A(6,5)=18; A(6,6)=11;

	//cout << "Enter number of 6x6 matrix SVDs (50000 is a nice number): ";
	//cin  >> n; 
	n = 20000;
	cout << n << " SVDs in progress" << endl;

	for (i=1;i<=n;i++)
		Svd(A, U, S, V);

	cout << "Done! Type 0 and press ENTER";
	//cin >> n;
}

int main(int,int)
{
	int before, after;

	before = mbytesfree();

	try
	{
		// Test1_Interp1();
		// Test1_Balance();
		// Test1_Eig();
		Test1_Performance_Mult();
		//Test1_Performance_Svd();
	}
	catch (EMatrixException Ex)
	{
		cout << Ex.what() << endl;
	}

	after = mbytesfree();

	cout << "Bytes lost during execution : " << before-after << endl;

	return 0;
}