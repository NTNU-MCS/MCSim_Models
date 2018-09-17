/***********************************************************************

  DEMO PROGRAM: Basic features for a C++ Matrix Library
  -----------------------------------------------------

	1.	It is vital that we can control the memory allocation to ensure
		that there will be no memory leak whatsoever. 
		The solution: Override the new-operator and use NEWALLOC instead.

	2.	Design a class hierarchy.

	3.	Override basic operators to facilitate code readability.



  Karl-Petter Lindegaard
  ABB Industri AS
  17.12.98
***********************************************************************/

#include <iostream>
#include <complex>
#include "mtrx.h"


void testfunc()
{
	double a;

	// Testing matrices
	matrix A(9,9);
	A(5,5) = 1607.1;
	matrix *B = new matrix(2,3);

	// Testing vectors
	matrix x(5);
	x(2) = 2.123;


	// Testing accessability methods
	cout << "x(2)  =" << x(2)   << endl;
	cout << "A(5,5)=" << A(5,5) << endl;

	// Adding of matrices
	matrix D = A+A;
	cout << "D(5,5) = " << D(5,5) << endl;


	// Testing the = operator on an allready declared matrix
	cout << endl << "Testing the = operator" << endl;
	D = A;

	// Testing Eye
	matrix EYE = Eye(3);
	cout << endl << "Eye(3) returned the following" << endl;
	cout << EYE << endl;

	D = EYE;
	D(1,2) = 2;
	D = 1607.0 * D;
	cout << "D =" << endl;
	cout << D;


	delete B;
}

/*
void testadd()
{
	vector a(3);
	vector b(3);
	vector c(4);

	a(1) = 100;
	b(3) = 200;

	vector d = a + b;
	vector e = a + c;

	cout << "Vector adding performed" << endl;
}
*/


void testexception()
{
	double a;
	matrix A(5,5);


	// Accessing an illegal element
	try
	{
		a = A(6,5);
		cout << "Back in business" << endl;		// Never reached
	}
	catch (EMatrixException Ex)
	{
		cout << "There was an error: " << Ex.What() << endl;
	}

	// Adding matrices with different sizes.
//	try
//	{ 
		matrix B = Eye(4);
		matrix C = A + B;
//	}
//	catch (EMatrixException Ex)
//	{
//		cout << "Try #2: " << Ex.What() << endl;
//	}
	
	// Multiplying a matrix with a vector
	cout << endl << "Multiply" << endl;

	A = Eye(3);
	A(2,3) = 4;
	matrix x(3);
	x(1) = 1;
	x(2) = 2;
	x(3) = 3.5;

	matrix x_dot = A*x;
	cout << endl << A << endl;
	cout << endl << x << endl;
	cout << endl << x_dot;
}

void testspeed()
{
	int i;

	matrix A(100,100);
	matrix B(100,200);

	cout << "Bytes free " << nbytesfree() << endl;

	for (i=1;i<=20;i++)
	{
		matrix C = A*B;
	}

	// cout << "Elapsed time was " << after.GetSecond() << " seconds" << endl;
}

void testsvd()
{
	char a;
	matrix A = Eye(4);

	matrix U, S, V, Sigma;
	matrix B, C;


	A(1,3) = 100;
	A(3,4) = -20;
	// A(4,4) = 0;
	cout << "A = " << A << endl;


	try {
		/*Sigma =*/ Svd(A,U,S,V);
		// cout << "U" << endl << U << endl << endl;
		// cout << "S" << endl << S << endl << endl;
		// cout << "V" << endl << V << endl << endl;
		// cin >> a;

	}
	catch (EMatrixException Ex) {
		cout << "Error while doing svd: " << Ex.What() << endl;
		return;
	}
	cout << "Singular values are " << endl << S/*igma*/ << endl;


	try {
		B = Inv(A);
	}
	catch (EMatrixException Ex) {
		cout << Ex.What() << endl;
		return;
	}
	cout << endl << "Inverse of A is" << endl << B << endl;

	C = A*B;
	cout << "A*inv(A) = " << C << endl;


	/***** BENCH-TEST ***************************************************
	cout << endl << "Starting bench...." << endl;
	for (int i=1;i<=10000;i++)
	{
		matrix D = Inv(A);
	}
	cout << i-1 << " inversions of a 4x4 matrix performed!" << endl;
	********************************************************************/
}


/*********************************************************************
  Main program
*********************************************************************/

void main() 
{
	int i;

	try {
		//testfunc();
		//testadd();
		//testexception();
		//testspeed();
		testsvd();
	}
	catch (EMatrixException Ex) {
		cout << "An error was catched by main()" << endl;
	}

	i = nmemleak();
	if (i)
		cout << "There was a memory leak!!!" << endl;

}
