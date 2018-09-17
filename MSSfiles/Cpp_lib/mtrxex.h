/*************************************************************************
 *
 * FILE		: $Archive: $ 
 * DATE		: 1999.02.04
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT	: Defines the matrix library's exception hierarchy. 
 *
 *			  The hierarchy is diveded into two sub-trees. The first
 *			  one consists of exceptions that occurs due to some 
 *			  programmatically caused bug, i.e. syntactic error. Typical
 *			  examples are addition of two matrices of different sizes.
 *			  So by tracking these exceptions individually, the developer
 *			  has a tool for tracking his/hers logic errors.
 *
 *			  The second sub-tree contains exceptions that may occur at
 *			  run-time. Such errors are e.g. trying to take the inverse 
 *			  of a singluar (but square) matrix and "out of memory" errors.
 *			  The common factor of all these exceptions are that the are
 *			  not necessarily caused by misspellings etc.
 *
 *			  This exception hierarchy is deduced from the common (ANSI?)
 *			  C++ exception class "exception".
 *            
 * NOTES	: 
 *
 * AUTHOR(s): Karl-Petter Lindegaard
 *
 ************************************************************************/ 

#ifndef MTRX_EX
#define MTRX_EX

#include <exception>
#include <string.h>



/***********************************************************

  Define our own super class in the exception hierarchy
  directly deduced from the (ANSI?) C++ exception class.

***********************************************************/

class EMatrixException // : public exception
{
protected:
	char reason[100];
	
public:
	EMatrixException()
	{
		strcpy (reason, "General matrix error");
	}
	virtual ~EMatrixException() {};

	virtual char *what()
	{
		return reason;
	}
};



/***********************************************************

  Define hierarchy of syntax caused matrix exceptions.
  
  If either one of these exceptions are raised, the source
  code is most likely incorrect.

***********************************************************/


class EMatrixLogic : public EMatrixException
{
public:
	EMatrixLogic()
	{
		strcpy (reason, "General matrix logic error. Check syntax.");
	}
};


class EMatrixEmpty : public EMatrixLogic
{
public:
	EMatrixEmpty()
	{
		strcpy (reason, "Matrix empty");
	}
};


class EMatrixIndex : public EMatrixLogic
{
public:
	EMatrixIndex()
	{
		strcpy (reason, "Indexing error");
	}
};


class EMatrixRow : public EMatrixIndex
{
public:
	EMatrixRow()
	{
		strcpy (reason, "Illegal row");
	}
};


class EMatrixCol : public EMatrixIndex
{
public:
	EMatrixCol()
	{
		strcpy (reason, "Illegal column");
	}
};


class EMatrixDimension : public EMatrixLogic
{
public:
	EMatrixDimension()
	{
		strcpy (reason, "Dimension error");
	}
};


class EMatrixCorrespond : public EMatrixDimension
{
public:
	EMatrixCorrespond()
	{
		strcpy (reason, "Dimensions does not correspond");
	}
};


class EMatrixNonsquare : public EMatrixDimension
{
public:
	EMatrixNonsquare()
	{
		strcpy (reason, "Nonsquare matrix");
	}
};


class EMatrixNovector : public EMatrixDimension
{
public:
	EMatrixNovector()
	{
		strcpy (reason, "Not a vector");
	}
};



/***********************************************************

  Define hierarchy of exceptions that track run-time
  errors that must be compensated for. 

***********************************************************/


class EMatrixRuntime : public EMatrixException
{
public:
	EMatrixRuntime()
	{
		strcpy (reason, "General matrix runtime exception");
	}
};

class EMatrixProperty : public EMatrixRuntime
{
public:
	EMatrixProperty()
	{
		strcpy (reason, "Matrix property");
	}
};

class EMatrixSingular : public EMatrixProperty
{
public:
	EMatrixSingular()
	{
		strcpy (reason, "Matrix singular");
	}
};

class EMatrixIllcondition : public EMatrixProperty
{
public:
	EMatrixIllcondition()
	{
		strcpy (reason, "Matrix is ill conditioned");
	}
};


class EMatrixMemory : public EMatrixRuntime
{
public:
	EMatrixMemory()
	{
		strcpy (reason, "Not enough memory");
	}
};

class EMatrixIteration : public EMatrixRuntime
{
public:
	EMatrixIteration()
	{
		strcpy (reason, "Max number of iterations reached");
	}
};



#endif