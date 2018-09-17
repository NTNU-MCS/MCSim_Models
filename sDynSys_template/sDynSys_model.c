/* 
    sDynSys_model.c

	Library file for functions defining the model of the dynamical system.

	References:		None

    Compile with:  > mex sDynSys_main.c sDynSys_model.c  ..\Libraries\C_lib\libGenericFunctions.c  

	Copyright: 		Roger Skjetne, NTNU
    Author: 		Roger Skjetne
    Date created: 	2011.04.13  RS Created for efficient modeling and simulation of control system simulations.
    Revised:      	

*/


#include "sDynSys_model.h"


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* #####################################  MODEL FUNCTIONS	 ############################################## */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void ModelOutputs(SimStruct *S, sSystemCFG *SystemStruct, real_T *xc, real_T *xd, real_T *u, real_T *y)
/*	Function that calculates the output y = h(xc,xd(k*Ts),u) of the system. */
{
	/*	Parameters:
			SystemStruct->Param[i].tag			i,j zero-indexed.
			SystemStruct->Param[i].size
			SystemStruct->Param[i].pvec[j],		The row array 'pvec[]' stores a matrix row-by-row. 
		Inputs: 
			u[0],  u[1],  u[2],  etc.
		Continuous states:
			xc[0], xc[1], xc[2], etc.
		Discrete states:
			xd[0], xd[1], xd[2], etc.

		Objective of function: To assign the outputs
			y[0],  y[1],  y[2],  etc.

	*/

	/*	Example model: 
				x_dot  = A*(x-xi) + xi_th*vs, 
				th_dot = vs - mu*W_th 
				y1 = x1, y2 = x2, y3 = xi1, y4 = xi2
	*/ 

	y[0] = xc[0];		// y1 = x1
	y[1] = xc[1];		// y2 = x2
	y[2] = cos(xc[2]);	// y3 = xi1[th]
	y[3] = -sin(xc[2]);	// y4 = xi2[th]

}

void ModelXcDerivatives(SimStruct *S, sSystemCFG *SystemStruct, real_T *xc_dot, real_T *xc, real_T *xd, real_T *u)
/*	Function that calculates the continuous state derivatives xc_dot = fc(xc,xd(k*Ts),u) of the system. */
{
	/*	Parameters:
			SystemStruct->Param[i].tag			i,j zero-indexed.
			SystemStruct->Param[i].size
			SystemStruct->Param[i].pvec[j],		The row array 'pvec[]' stores a matrix row-by-row. 
		Inputs: 
			u[0],  u[1],  u[2],  etc.
		Continuous states:
			xc[0], xc[1], xc[2], etc.
		Discrete states:
			xd[0], xd[1], xd[2], etc.

		Objective of function: To assign the continuous state derivatives
			xc_dot[0],	xc_dot[1],	xc_dot[2],  etc. 
	*/

	/*	Example model: 
				x_dot  = A*(x-xi) + xi_th*vs, 
				th_dot = vs - mu*W_th 
	*/ 

	real_T		mu, A[2][2], P[2][2], W_th, th, x1, x2, xi1, xi2, xi1_th, xi2_th, vs;

	A[0][0] = SystemStruct->Param[0].pvec[0];	// Element 1
	A[0][1] = SystemStruct->Param[0].pvec[1];	// Element 2
	A[1][0] = SystemStruct->Param[0].pvec[2];	// Element 3
	A[1][1] = SystemStruct->Param[0].pvec[3];	// Element 4
	mu		= SystemStruct->Param[1].pvec[0];
	P[0][0] = SystemStruct->Param[2].pvec[0];	// Element 1
	P[0][1] = SystemStruct->Param[2].pvec[1];	// Element 2
	P[1][0] = SystemStruct->Param[2].pvec[2];	// Element 3
	P[1][1] = SystemStruct->Param[2].pvec[3];	// Element 4

	x1 = xc[0];
	x2 = xc[1];
	th = xc[2];
	
	xi1 = cos(th);	xi1_th = -sin(th);
	xi2 = -sin(th);	xi2_th = -cos(th);

	vs   = u[0];
	W_th = -2*(x1-xi1)*(P[0][0]*xi1_th + P[0][1]*xi2_th) - 2*(x2-xi2)*(P[1][0]*xi1_th + P[1][1]*xi2_th);

	xc_dot[0] = A[0][0]*(x1-xi1) + A[0][1]*(x2-xi2) + xi1_th*vs;
	xc_dot[1] = A[1][0]*(x1-xi1) + A[1][1]*(x2-xi2) + xi2_th*vs;
	xc_dot[2] = vs - mu*W_th;

}

void ModelXdUpdate(SimStruct *S, sSystemCFG *SystemStruct, real_T *xc, real_T *xd, real_T *u)
/*	Function that calculates the discrete state updates xd(k+1) = fd(xd(k),xc,u) of the system. */
{
	/*	Parameters:
			SystemStruct->Param[i].tag			i,j zero-indexed.
			SystemStruct->Param[i].size
			SystemStruct->Param[i].pvec[j],		The row array 'pvec[]' stores a matrix row-by-row. 
		Inputs: 
			u[0],  u[1],  u[2],  etc.
		Continuous states:
			xc[0], xc[1], xc[2], etc.

		Objective of function: To assign the discrete state updates
			xd[0], xd[1], xd[2], etc.

	*/

	/*	Example model: 
				No discrete states 
	*/ 

	;	// No discrete states.
}



