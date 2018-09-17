
/* 
    sDynSys_model.h

	Header file for functions defining the model of the dynamical system.

	References:		None

	Copyright: 		Roger Skjetne, NTNU
    Author: 		Roger Skjetne
    Date created: 	2011.04.13  RS Created for efficient modeling and simulation of control system simulations.
    Revised:      	

*/


#ifndef S_DYNSYS_MODEL_H
#define S_DYNSYS_MODEL_H

#include "simstruc.h"
#include <time.h>
#include <math.h>
#include "..\Libraries\C_lib\libSparseMatrix.h"
#include "..\Libraries\C_lib\libDataStructures.h"
#include "..\Libraries\C_lib\libGenericFunctions.h"
//#include "..\Libraries\C_lib\mtrxc.h"


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* #####################################  Function declarations	################################### */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void ModelOutputs(SimStruct *S, sSystemCFG *SystemStruct, real_T *xc, real_T *xd, real_T *u, real_T *y);
void ModelXcDerivatives(SimStruct *S, sSystemCFG *SystemStruct, real_T *xc_dot, real_T *xc, real_T *xd, real_T *u);
void ModelXdUpdate(SimStruct *S, sSystemCFG *SystemStruct, real_T *xc, real_T *xd, real_T *u);


#endif