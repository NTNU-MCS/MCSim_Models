/*
    libGenericFunctions.h

	Header file for generic simulation functions related to the sDynSys s-function.

	References:
	[1]	W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery, Numerical Recipes in C - The Art of Scientific Computing, 2nd ed.
		Cambridge University Press, New York, USA, 1992.

	Copyright: 		Roger Skjetne, NTNU
    Date created: 	2011.04.13 Roger Skjetne
    Revised:		2011.07.12 RS Increased the MAX_xxx macros.

*/


#ifndef GENERIC_FUNCTIONS_H
#define GENERIC_FUNCTIONS_H

#include "simstruc.h"
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <Windows.h>

#define pi				3.141592653589793
#define epsilon			1.0e-10	// Small number

#define MAX_MSG_CHAR	64				/* Max number of characters in message */
#define MAX_MSG			200				/* Max number of messages to be sent with SendMsg */
#define MAX_STATUS		2				/* Max quarantine status given to message. Must be 2 or higher. */
#define ERROR_LOG		0
#define INFO_LOG		4
#define EVENT_LOG		5
#define DEBUG_LOG		9
#define FILE_LOG		10

#define MAX_STRING_SIZE	128					/* Max number of characters in strings */
#define MAX_LINE_SIZE	(8*MAX_STRING_SIZE)	/* Max number of characters to be read from a string line	*/
#define MAX_PARAMS		128					/* Max number of specified parameters */
#define MAX_PARAM_SIZE	256					/* Max number of elements in a parameter vector */


typedef struct {
/* Defining a structure to hold all the parameter configuration data */
	char		tag[MAX_STRING_SIZE+1];	// Tag of parameter
	int			size;					// Size of parameter vector
	double		pvec[MAX_PARAM_SIZE];	// Parameter vector
} sParams;
#define PARAM_FIELDS	3				/* Number of parameter fields in sParams struct */

typedef struct {
/* Defining super structure to hold config data with main info */
	double		SampleTime;				// Main sample time of s-function
	int			NumContStates;			// Number of continuous states
	int			NumDiscStates;			// Number of discrete states
	int			NumInputSignals;		// Number of control input signals
	int			NumOutputSignals;		// Number of output signals
	int			NumParameters;			// Number of parameters
	sParams		Param[MAX_PARAMS+1];	// Parameter data
} sSystemCFG;

typedef struct {
/* Defining a structure to message information */
	char	subfunc[MAX_MSG_CHAR+1];	// Name of subfunction sending message
	char	msg[8*MAX_MSG_CHAR+1];		// Message string
	int		msg_length;					// Length of message
	int		level;						// Priority level of message
	int		status;						// Activity status of message (0: Inactive, >0: Active)

} sMsgList;

int    Param2Struct(SimStruct *S, sSystemCFG *SystemStruct, char *filename);
void   SendMsg(SimStruct *S, int level, char *msg, char *subfunc);
double ticker(SimStruct *S);
int	   comp(const void *i, const void *j);
double rad2pipi(double x);

#endif