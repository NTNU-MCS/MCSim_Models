/*
 * sfun_forcealloc.cpp: S-Function for force allocation Cybership II.
 *--------------------------------------------------------------------
 * This S-Function implements the first version of nullspace based
 * force allocation for Cybership II.
 *
 * See \academic\lab\thralloc\nullspace for further details.
 *
 * To compile: mex -I\academic\lab\mtrx sfun_forcealloc.cpp \academic\lab\mtrx\mtrx.lib
 *
 * Copyright (c) 2000 Karl-Petter Lindegaard
 *
 */


/*
 * You must specify the S_FUNCTION_NAME as the name of your S-function
 * (i.e. replace sfuntmpl with the name of your S-function).
 */

#define S_FUNCTION_LEVEL 2
#define S_FUNCTION_NAME  sfun_forcealloc

/*
 * Need to include simstruc.h for the definition of the SimStruct and
 * its associated macro definitions.
 */
#include "simstruc.h"


/*
 * Include matrix library
 */
#include "mtrx.h"


/*
 * Include vessel data
 */
#include "cs2data.h"


#ifdef __cplusplus
extern "C" { // use the C fcn-call standard for all functions  
#endif       // defined within this scope                     


static double k1 = (T1LX-T2LX)/(T3LY-T2LY);
static double k2 = (T1LX-T3LX)/(T2LY-T3LY);

/*
 * Define macros for simpler access to persistent matrix objects
 */
#define getA1dagger(s) (ssGetPWork(s)[0])
#define getA2dagger(s) (ssGetPWork(s)[1])
#define getN1(s) (ssGetPWork(s)[2])
#define getN2(s) (ssGetPWork(s)[3])
#define getQ1(s) (ssGetPWork(s)[4])
#define getQ2(s) (ssGetPWork(s)[5])

// #define U(element) (*uPtrs[element])  /* Pointer to Input Port0 */


matrix nullsub1(const matrix&, const matrix&, const matrix&, const matrix&, double);
matrix nullsub2(const matrix&, const matrix&, const matrix&, const matrix&, double);




/* Error handling
 * --------------
 *
 * You should use the following technique to report errors encountered within
 * an S-function:
 *
 *       ssSetErrorStatus(S,"Error encountered due to ...");
 *       return;
 *
 * Note that the 2nd argument to ssSetErrorStatus must be persistent memory.
 * It cannot be a local variable. For example the following will cause
 * unpredictable errors:
 *
 *      mdlOutputs()
 *      {
 *         char msg[256];         {ILLEGAL: to fix use "static char msg[256];"}
 *         sprintf(msg,"Error due to %s", string);
 *         ssSetErrorStatus(S,msg);
 *         return;
 *      }
 *
 * See matlabroot/simulink/src/sfunctmpl.doc for more details.
 */

/*====================*
 * S-function methods *
 *====================*/

/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{
    /* See sfuntmpl.doc for more details on the macros below */

    ssSetNumSFcnParams(S, 2);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        /* Return if number of expected != number of actual parameters */
        return;
    }

    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    if (!ssSetNumInputPorts(S, 1)) return;
    ssSetInputPortWidth(S, 0, 3);
    ssSetInputPortDirectFeedThrough(S, 0, 3);

    if (!ssSetNumOutputPorts(S, 1)) return;
    ssSetOutputPortWidth(S, 0, 5);

    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 6);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    ssSetOptions(S, 0);
}



/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    This function is used to specify the sample time(s) for your
 *    S-function. You must register the same number of sample times as
 *    specified in ssSetNumSampleTimes.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);

}



#define MDL_INITIALIZE_CONDITIONS   /* Change to #undef to remove function */
#if defined(MDL_INITIALIZE_CONDITIONS)
  /* Function: mdlInitializeConditions ========================================
   * Abstract:
   *    In this function, you should initialize the continuous and discrete
   *    states for your S-function block.  The initial states are placed
   *    in the state vector, ssGetContStates(S) or ssGetRealDiscStates(S).
   *    You can also perform any other initialization activities that your
   *    S-function may require. Note, this routine will be called at the
   *    start of simulation and if it is present in an enabled subsystem
   *    configured to reset states, it will be call when the enabled subsystem
   *    restarts execution to reset the states.
   */
  static void mdlInitializeConditions(SimStruct *S)
  {
  }
#endif /* MDL_INITIALIZE_CONDITIONS */



#define MDL_START  /* Change to #undef to remove function */
#if defined(MDL_START) 
  /* Function: mdlStart =======================================================
   * Abstract:
   *    This function is called once at start of model execution. If you
   *    have states that should be initialized once, this is the place
   *    to do it.
   */
  static void mdlStart(SimStruct *S)
  {
	  matrix Q1 = Eye(4);
	  matrix Q2 = Eye(4);

	  matrix A1 = Zeros(3,4);
	  matrix n1(4,1);

	  matrix A2 = Zeros(3,4);
	  matrix n2(4,1);

	  // Configure A1, n1 and A1dagger etc.
	  //-----------------------------------
	  A1(1,2) = 1.0; A1(1,4) = 1.0;
	  A1(2,1) = 1.0; A1(2,3) = 1.0;
	  A1(3,1) = T1LX; A1(3,2)=-T2LY; A1(3,3)=T2LX; A1(3,4)=-T3LY;

	  n1(1)=-1.0; n1(2)=k1; n1(3)=1.0; n1(4)=-k1;

	  matrix A1T = A1.Transp();
	  matrix A1D = A1T*Inv(A1*A1T);


	  // Configure A2, n2 and A2dagger etc.
	  //-----------------------------------
	  A2(1,2) = 1.0; A2(1,3) = 1.0;
	  A2(2,1) = 1.0; A2(2,4) = 1.0;
	  A2(3,1) = T1LX; A2(3,2)=-T2LY; A2(3,3)=-T3LY; A2(3,4)=T3LX;

	  n2(1)=-1.0; n2(2)=-k2; n2(3)=k2; n2(4)=1.0;

	  matrix A2T = A2.Transp();
	  matrix A2D = A2T*Inv(A2*A2T);


	  // Store persistent variables
	  //---------------------------
	  getA1dagger(S) = (void *)new matrix(A1D);
	  getN1(S)       = (void *)new matrix(n1);
	  getQ1(S)       = (void *)new matrix(Q1);

	  getA2dagger(S) = (void *)new matrix(A2D);
	  getN2(S)       = (void *)new matrix(n2);
	  getQ2(S)       = (void *)new matrix(Q2);
  }
#endif /*  MDL_START */



/* Function: mdlOutputs =======================================================
 * Abstract:
 *    In this function, you compute the outputs of your S-function
 *    block. Generally outputs are placed in the output vector, ssGetY(S).
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
    real_T            *y    = ssGetOutputPortRealSignal(S,0);
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0);

	// Get parameters
	double k2alfa = *mxGetPr(ssGetSFcnParam(S, 0));
	double k3alfa = *mxGetPr(ssGetSFcnParam(S, 1));

	// Access persistent variables
	matrix *A1dagger = (matrix*)getA1dagger(S);
	matrix *N1       = (matrix*)getN1(S);
	matrix *Q1       = (matrix*)getQ1(S);

	matrix *A2dagger = (matrix*)getA2dagger(S);
	matrix *N2       = (matrix*)getN2(S);
	matrix *Q2       = (matrix*)getQ2(S);


	// Create matrices that's used for solving the sector problems
	matrix H1 = Zeros(3,3);
	matrix H2 = Zeros(3,3);

	H1(1,1)=k2alfa; H1(1,2)=-1.0;
	H1(2,1)=k2alfa;               H1(2,3)=-1;
	H1(3,1)=1.0;                  H1(3,3)=-k1;

	H2(1,1)=k3alfa; H2(1,2)=1.0;
	H2(2,1)=1.0;                  H2(2,3)=-k2;
	                H2(3,2)=1.0;  H2(3,3)=-1.0;

	matrix H1inv = Inv(H1);   // Rather slow....
	matrix H2inv = Inv(H2);   // Rather slow....

	// Map inputs (commanded forces Tauc)
	matrix Tauc(3,1);
	Tauc(1) = *uPtrs[0];
	Tauc(2) = *uPtrs[1];
	Tauc(3) = *uPtrs[2];

	// Call subroutines for each rudder
	matrix x1 = nullsub1(Tauc, *A1dagger, H1inv, *N1, k2alfa);
	matrix x2 = nullsub2(Tauc, *A2dagger, H2inv, *N2, k3alfa);

	// Compare results and pick the best solution
	matrix J1 = x1.Transp()*(*Q1)*x1;
	matrix J2 = x2.Transp()*(*Q2)*x2;
	if (J1(1) <= J2(1))
	{
		// Use: u = [x1(1) x1(2) x1(3) x1(4) 0];
		y[0]=x1(1); y[1] = x1(2); y[2]=x1(3); y[3]=x1(4); y[4]=0.0;
	}
	else
	{
		// u = [x2(1) x2(2) 0 x2(3) x2(4)];
		y[0]=x2(1); y[1] = x2(2); y[2]=0.0; y[3]=x2(3); y[4]=x2(4);
	}
}



#define MDL_UPDATE  /* Change to #undef to remove function */
#if defined(MDL_UPDATE)
  /* Function: mdlUpdate ======================================================
   * Abstract:
   *    This function is called once for every major integration time step.
   *    Discrete states are typically updated here, but this function is useful
   *    for performing any tasks that should only take place once per
   *    integration step.
   */
  static void mdlUpdate(SimStruct *S, int_T tid)
  {
  }
#endif /* MDL_UPDATE */



#define MDL_DERIVATIVES  /* Change to #undef to remove function */
#if defined(MDL_DERIVATIVES)
  /* Function: mdlDerivatives =================================================
   * Abstract:
   *    In this function, you compute the S-function block's derivatives.
   *    The derivatives are placed in the derivative vector, ssGetdX(S).
   */
  static void mdlDerivatives(SimStruct *S)
  {
  }
#endif /* MDL_DERIVATIVES */



/* Function: mdlTerminate =====================================================
 * Abstract:
 *    In this function, you should perform any actions that are necessary
 *    at the termination of a simulation.  For example, if memory was
 *    allocated in mdlStart, this is the place to free it.
 */
static void mdlTerminate(SimStruct *S)
{
	matrix *A1dagger = (matrix*)getA1dagger(S);
	matrix *N1       = (matrix*)getN1(S);
	matrix *Q1       = (matrix*)getQ1(S);
	matrix *A2dagger = (matrix*)getA2dagger(S);
	matrix *N2       = (matrix*)getN2(S);
	matrix *Q2       = (matrix*)getQ2(S);

	delete A1dagger;
	delete N1;
	delete Q1;
	delete A2dagger;
	delete N2;
	delete Q2;
}


/*======================================================*
 * See sfuntmpl.doc for the optional S-function methods *
 *======================================================*/


matrix nullsub1(const matrix& Tauc, 
				const matrix& A1dagger, 
				const matrix& H1inv, 
				const matrix& n1,
				double k2alfa)
{
	matrix x0 = A1dagger*Tauc;

	double x2 = x0(2);
	double y2 = x0(3);

	double l1 = y2 - k2alfa*x2;
	double l2 = -k1*y2 + x2;
	double lambda1;

	// Check sectors
	if ((y2<0) || (l2<0))
		lambda1 = -y2;
	else if ((l1>=0.0) && (l2>=0.0))
	{
		lambda1 = H1inv(3,2)*x0(3) + H1inv(3,3)*x0(2);
	}
	else
		lambda1 = 0;


	matrix x1;
	x1 = x0 + lambda1*n1;

	return x1;
}


matrix nullsub2(const matrix& Tauc, 
				const matrix& A2dagger, 
				const matrix& H2inv, 
				const matrix& n2,
				double k3alfa)
{
	matrix x0 = A2dagger*Tauc;

	double x3 = x0(3);
	double y3 = x0(4);

	double l1 = -y3 - k3alfa*x3;
	double l2 = -k2*y3 + x3;
	double lambda2;

	// Check sectors
	if ((y3>0) || (l2<0))
		lambda2 = -y3;
	else if ((l1>=0.0) && (l2>=0.0))
	{
		lambda2 = H2inv(3,2)*x0(3) + H2inv(3,3)*x0(4);
	}
	else
		lambda2 = 0;


	matrix x2;
	x2 = x0 + lambda2*n2;

	return x2;
}




/*=============================*
 * Required S-function trailer *
 *=============================*/


#ifdef __cplusplus
} // end of extern "C" scope
#endif

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif



