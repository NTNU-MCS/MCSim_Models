/* 
    sDynSys_main.c

	The sDynSys_main module implements an s-function for a dynamical system. The function calls the
	subfunctions in sDynSys_model.c where the models of the system is implemented.


    Input parameters:

	Parameter 0 -	Simulation parameters: [SampleTime, NumContStates, NumDiscStates, NumInputSignals, NumOutputSignals, NumParameters]
	Parameter 1 -	File name: Configuration data initialization file.
	Parameter 2 -	Xc0: Initial conditions continuous states.
	Parameter 3 -	Xd0: Initial conditions discrete states.


	Input signals (size):
	
	I0 - Mode control (Not implemented functionality yet):
						0: Debug information report flag; See function DebugReport for description of the 
							input values and generated reports.
						1: Reset; Resets module to initialization w.r.t. configuration. 
	I1 - Control input u



    Output signals (size):
	
	O0 - System outputs y
	O1 - Discrete states Xd
	O2 - Continuous states Xc 
	O3 - Continuous state derivatives Xc_dot 



    Compile with:  > mex sDynSys_main.c sDynSys_model.c  ..\Libraries\C_lib\libGenericFunctions.c  

	Copyright: 		Roger Skjetne, NTNU
    Author: 		Roger Skjetne
    Date created: 	2011.04.13  RS Created for efficient modeling and simulation of control system simulations.
    Revised:		2011.07.11  RS Updated with an extra output port for continuous state derivatives Xc_dot. 
			
*/


#define S_FUNCTION_NAME  sDynSys_main
#define S_FUNCTION_LEVEL 2


/*
 * Need to include simstruc.h for the definition of the SimStruct and
 * its associated macro definitions.
 */


#include "simstruc.h"
#include "sDynSys_model.h"
#include "..\Libraries\C_lib\libSparseMatrix.h"
#include "..\Libraries\C_lib\libDataStructures.h"
#include "..\Libraries\C_lib\libGenericFunctions.h"
#include "..\Libraries\C_lib\mtrxc.h"

#define FILE_NAME		"sDynSys_main"


/****************************************************************************
 * Input parameters       		  						                    *
 *		-PS: Specify the parameters as one column vector, row after row!    *
 ****************************************************************************/
#define SimParam_IDX 		0	// [SampleTime, NumContStates, NumDiscStates, NumInputSignals, NumOutputSignals, NumParameters]
#define SimParam_PARAM(S) ssGetSFcnParam(S,SimParam_IDX)

#define InitFname_IDX 		1	// Filename initialization (configuration) file.
#define InitFname_PARAM(S) ssGetSFcnParam(S,InitFname_IDX)

#define Xc0_IDX 			2	// Initial conditions continuous states
#define Xc0_PARAM(S) ssGetSFcnParam(S,Xc0_IDX)

#define Xd0_IDX 			3	// Initial conditions discrete states
#define Xd0_PARAM(S) ssGetSFcnParam(S,Xd0_IDX)


/*********************************
 * Vector and array sizes		 *
 *********************************/



#define NP 		4	/* Number of INPUT PARAMETERS							*/

#define NIP     2	/* Number of Input Ports								*/
#define NU0   	2	/* Input 0: Mode control input vector					*/

#define NOP     4	/* Number of Output Ports								*/



/* ========================================================== *
	Module special functions not included in library files
 * ========================================================== */

/*
void DebugReport(SimStruct *S, sSystemData *SystemData, sLocationData *LocationData, sPipeData *PipeData, sInputData *InputData, sPipestr_ptr *TouchDetect, int *ProximityDetect, int *InfoFlag)
*/







/*====================*
 * S-function methods *
 *====================*/


#define MDL_CHECK_PARAMETERS
#if defined(MDL_CHECK_PARAMETERS) && defined(MATLAB_MEX_FILE)
  /* Function: mdlCheckParameters =============================================== 
   * Abstract: Validate our parameters to verify they are okay. */
	static void mdlCheckParameters(SimStruct *S)
	{
		double	*simPtr = mxGetPr(SimParam_PARAM(S));
		int		NumContStates = (int)simPtr[1];
		int		NumDiscStates = (int)simPtr[2];

		/* Check parameters: */
		{
			if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
				/* Return if number of expected != number of actual parameters */
				ssSetErrorStatus(S, "The number of parameters to S-function does "
									"not match the expected number.");
				return;
			}
		}
	 
		/* Check 1st parameter: CtrlParam */
		{
			if (!mxIsNumeric(SimParam_PARAM(S)) ||
				 mxGetNumberOfElements(SimParam_PARAM(S)) != 6 ||	
				 simPtr[0] < -1.0 ||		// SampleTime
				 simPtr[1] <  0 ||			// NumContStates
				 simPtr[2] <  0 ||			// NumDiscStates
				 simPtr[3] <= 0 ||			// NumInputSignals
				 simPtr[4] <= 0 ||			// NumOutputSignals
				 simPtr[5] <= 0) {			// NumParameters
					 ssSetErrorStatus(S,"1st parameter to S-function must be a vector of size 6 with the simulation parameters: "
									"SampleTime	    Ts = -1 (variable) or Ts > 0.0 (fixed), "
									"NumContStates    >= 0, "
									"NumDiscStates    >= 0, "
									"NumInputSignals  >= 1, "
									"NumOutputSignals >= 1, and "
									"NumParameters    >= 1. ");
				return;
			}
		}
		  /* Check 2nd parameter: Filename to config file. */
		{
			if (!mxIsChar(InitFname_PARAM(S))) {
				ssSetErrorStatus(S, "2nd parameter to S-function must be a string "
									"specifying the filename with path to the file "
									"giving the initialization data.");
				return;
			}
		}
		/* Check 3rd parameter: Xc0 - Initial conditions continuous states */
		{
			if (!mxIsNumeric(Xc0_PARAM(S)) ||
				 mxGetNumberOfElements(Xc0_PARAM(S)) != NumContStates) {
				 ssSetErrorStatus(S,"3rd parameter to S-function, Xc0, must be a vector of initial "
									"conditions for the continuous states of the dynamic system.");
				return;
			}
		}
		/* Check 4th parameter: Xc0 - Initial conditions discrete states */
		{
			if (!mxIsNumeric(Xd0_PARAM(S)) ||
				 mxGetNumberOfElements(Xd0_PARAM(S)) != NumDiscStates) {
				 ssSetErrorStatus(S,"4th parameter to S-function, Xd0, must be a vector of initial "
									"conditions for the discrete states of the dynamic system.");
				return;
			}
		}
	}

#endif /* MDL_CHECK_PARAMETERS */
 


/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{

    ssSetNumSFcnParams(S, NP);  /* Number of expected parameters */
	#if defined(MATLAB_MEX_FILE)
		if (ssGetNumSFcnParams(S) == ssGetSFcnParamsCount(S)) {
			mdlCheckParameters(S);
			if (ssGetErrorStatus(S) != NULL) {
				return;
			}
		} 
		else {
			return; /* Parameter mismatch will be reported by Simulink */
		}
	#endif

    ssSetNumContStates(S, DYNAMICALLY_SIZED);
    ssSetNumDiscStates(S, DYNAMICALLY_SIZED);

    if (!ssSetNumInputPorts(S, NIP)) return;
	ssSetInputPortWidth(S, 0, NU0);
	ssSetInputPortWidth(S, 1, DYNAMICALLY_SIZED);

    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortDirectFeedThrough(S, 1, 1);

    if (!ssSetNumOutputPorts(S, NOP)) return;
	ssSetOutputPortWidth(S, 0, DYNAMICALLY_SIZED);
	ssSetOutputPortWidth(S, 1, DYNAMICALLY_SIZED);
	ssSetOutputPortWidth(S, 2, DYNAMICALLY_SIZED);
	ssSetOutputPortWidth(S, 3, DYNAMICALLY_SIZED);

	ssSetNumSampleTimes(S, 1);

    ssSetNumRWork(S, DYNAMICALLY_SIZED);
    ssSetNumIWork(S, DYNAMICALLY_SIZED);
    ssSetNumPWork(S, DYNAMICALLY_SIZED);

	ssSetNumModes(S, DYNAMICALLY_SIZED);
    ssSetNumNonsampledZCs(S, DYNAMICALLY_SIZED);

    /* Take care when specifying exception free code - see sfuntmpl.doc */
    //ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}

#if defined(MATLAB_MEX_FILE)


	#define MDL_SET_DEFAULT_PORT_DIMENSION_INFO
	void mdlSetDefaultPortDimensionInfo(SimStruct *S)
	{
		double	*simPtr = mxGetPr(SimParam_PARAM(S));
		int		NumContStates	 = (int)simPtr[1];
		int		NumDiscStates	 = (int)simPtr[2];
		int		NumInputSignals	 = (int)simPtr[3];
		int		NumOutputSignals = (int)simPtr[4];

		if (ssGetInputPortWidth(S, 0) == DYNAMICALLY_SIZED) {
			ssSetInputPortWidth(S, 0, NU0);
		}
		if (ssGetInputPortWidth(S, 1) == DYNAMICALLY_SIZED) {
			if (NumInputSignals < 1) 
				ssSetInputPortWidth(S, 0, 1);
			else
				ssSetInputPortWidth(S, 0, NumInputSignals);
		}
		if (ssGetOutputPortWidth(S, 0) == DYNAMICALLY_SIZED) {
			if (NumOutputSignals < 1) 
				ssSetOutputPortWidth(S, 0, 1);
			else
				ssSetOutputPortWidth(S, 0, NumOutputSignals);
		}
		if (ssGetOutputPortWidth(S, 1) == DYNAMICALLY_SIZED) {
			if (NumDiscStates < 1) 
				ssSetOutputPortWidth(S, 1, 1);
			else
				ssSetOutputPortWidth(S, 1, NumDiscStates);
		}
		if (ssGetOutputPortWidth(S, 2) == DYNAMICALLY_SIZED) {
			if (NumContStates < 1) 
				ssSetOutputPortWidth(S, 2, 1);
			else
				ssSetOutputPortWidth(S, 2, NumContStates);
		}
		if (ssGetOutputPortWidth(S, 3) == DYNAMICALLY_SIZED) {
			if (NumContStates < 1) 
				ssSetOutputPortWidth(S, 3, 1);
			else
				ssSetOutputPortWidth(S, 3, NumContStates);
		}
	}

	# define MDL_SET_INPUT_PORT_WIDTH
	static void mdlSetInputPortWidth(SimStruct *S, int_T port, int_T inputPortWidth)
	{
		double	*simPtr = mxGetPr(SimParam_PARAM(S));
		int		NumInputSignals	 = (int)simPtr[3];


		  if(port == 1) {
			  if ((inputPortWidth != NumInputSignals) && (inputPortWidth != 1)) {
				  ssSetErrorStatus(S,"sDynSys_main: Width for input port 1 (Control input u) should be equal to 1 or "
								   "the value given by parameter 0 (SimParam), element 3, i.e. NumInputSignals.");
				return;
			  }
			  else {
				if (NumInputSignals < 1) 
					ssSetInputPortWidth(S,port,1);				/* Control input signals	*/
				else
					ssSetInputPortWidth(S,port,inputPortWidth);	/* Control input signals	*/ 
			  }
		  }
		  else {
			ssSetErrorStatus(S,"sDynSys_main: mdlSetInputPortWidth called with erroneous port number.");
			return;
		  }
	}

	# define MDL_SET_OUTPUT_PORT_WIDTH
	static void mdlSetOutputPortWidth(SimStruct *S, int_T port, int_T outputPortWidth)
	{
		double	*simPtr = mxGetPr(SimParam_PARAM(S));
		int		NumContStates	 = (int)simPtr[1];
		int		NumDiscStates	 = (int)simPtr[2];
		int		NumOutputSignals = (int)simPtr[4];

		  if ((port == 0) && ((outputPortWidth != NumOutputSignals) && (outputPortWidth != 1))) {
			  ssSetErrorStatus(S,"sDynSys_main: Width of output port 0 (Output signals y) should be equal to 1 or "
								 "the value given by parameter 0 (SimParam), element 4, i.e. NumOutputSignals.");
			  return;
		  }
		  else if ((port == 1) && ((outputPortWidth != NumDiscStates) && (outputPortWidth != 1))) {
			  ssSetErrorStatus(S,"sDynSys_main: Width of output port 1 (Discrete states Xd) should be equal to 1 or "
								 "the value given by parameter 0 (SimParam), element 2, i.e. NumDiscStates.");
			  return;
		  }
		  else if ((port == 2) && ((outputPortWidth != NumContStates) && (outputPortWidth != 1))) {
			  ssSetErrorStatus(S,"sDynSys_main: Width of output port 2 (Continuous states Xc) should be equal to 1 or "
								 "the value given by parameter 1 (SimParam), element 1, i.e. NumContStates.");
			  return;
		  }
		  else if ((port == 3) && ((outputPortWidth != NumContStates) && (outputPortWidth != 1))) {
			  ssSetErrorStatus(S,"sDynSys_main: Width of output port 3 (Continuous state derivatives Xc_dot) should be equal to 1 or "
								 "the value given by parameter 1 (SimParam), element 1, i.e. NumContStates.");
			  return;
		  }

		  ssSetOutputPortWidth(S,port,outputPortWidth);

	}

	# define MDL_SET_WORK_WIDTHS
	static void mdlSetWorkWidths(SimStruct *S)
	{
		double	*simPtr = mxGetPr(SimParam_PARAM(S));
		int		NumContStates	 = (int)simPtr[1];
		int		NumDiscStates	 = (int)simPtr[2];

		if (NumContStates < 1) 
			ssSetNumContStates(S, 0);
		else
			ssSetNumContStates(S, NumContStates);

		if (NumDiscStates < 1) 
			ssSetNumDiscStates(S, 0);
		else
			ssSetNumDiscStates(S, NumDiscStates);


		ssSetNumRWork(S, 0);

		// I-work vectors
		//0 FAIL_FLAG
		ssSetNumIWork(S, 1);

		// R-work vectors
		//0 Xc_dot
		if (NumContStates < 1) 
			ssSetNumRWork(S, 0);
		else
			ssSetNumRWork(S, NumContStates);

		// P-work vectors
		//0 MsgList
		//1 tick
		//2 System
		ssSetNumPWork(S, 3);

 		ssSetNumModes(S, 0);
		ssSetNumNonsampledZCs(S, 0);
	}

#endif



/* Function: mdlInitializeSampleTimes =========================================
 * Abstract: S-function is comprised of only continuous sample time elements */
static void mdlInitializeSampleTimes(SimStruct *S)
{
	double	*simPtr = mxGetPr(SimParam_PARAM(S));
	double	Ts		= simPtr[0];

	if(Ts <= 0.0) {
		ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
		ssSetOffsetTime(S, 0, 0.0);
	}
	else {
		ssSetSampleTime(S, 0, Ts);
		ssSetOffsetTime(S, 0, 0.0);
	}
}



#define MDL_INITIALIZE_CONDITIONS
/* Function: mdlInitializeConditions ========================================
 * Abstract:
 *    If the initial condition parameter (X0) is not an empty matrix,
 *    then use it to set up the initial conditions, otherwise,
 *    set the intial conditions to all 0.0
 */
static void mdlInitializeConditions(SimStruct *S)
{
	double	*simPtr = mxGetPr(SimParam_PARAM(S));
	int		NumContStates = (int)simPtr[1];
	int		NumDiscStates = (int)simPtr[2];
	int		i;

    real_T  *xc0  = ssGetContStates(S);
	real_T	*xd0  = ssGetDiscStates(S);
    int_T   numXc = ssGetNumContStates(S);
	int_T	numXd = ssGetNumDiscStates(S);
	real_T *pr;

	if((NumContStates <= 0) || (mxGetM(Xc0_PARAM(S)) == 0)) {
		for (i=0; i<numXc; i++) {
			xc0[i] = 0.0;
		}
	}
	else {
		pr = mxGetPr(Xc0_PARAM(S));
		for (i=0; i<numXc; i++) {
			xc0[i] = pr[i];
		}
	}

	if((NumDiscStates <= 0) || (mxGetM(Xd0_PARAM(S)) == 0)) {
		for (i=0; i<numXd; i++) {
			xd0[i] = 0.0;
		}
	}
	else {
		pr = mxGetPr(Xd0_PARAM(S));
		for (i=0; i<numXd; i++) {
			xd0[i] = pr[i];
		}
	}
}


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
		int		i, FAIL_FLAG = 1, status;
		double	*tick;

		double	*simPtr = mxGetPr(SimParam_PARAM(S));
		int		NumContStates	 = (int)simPtr[1];
		int		NumDiscStates	 = (int)simPtr[2];
		int		NumInputSignals	 = (int)simPtr[3];
		int		NumOutputSignals = (int)simPtr[4];
		int		NumParameters	 = (int)simPtr[5];

		int		nu = (int)mxGetNumberOfElements(InitFname_PARAM(S));
		char	filename[MAX_STRING_SIZE], msg[MAX_LINE_SIZE];


		/* Message function initialization */
		const char_T	*sfunc		= ssGetModelName(S);
		const char_T	*sfunc_path	= ssGetPath(S);
		char			logfile[MAX_MSG_CHAR+1];
		sMsgList		*MsgList;
		FILE			*fp;

		sSystemCFG			*System;


		real_T				*Xd		= ssGetDiscStates(S);
		real_T				*Xc     = ssGetContStates(S);
		real_T				*Xc_dot;

		int_T				y0Width	= ssGetOutputPortWidth(S,0);
		int_T				y1Width	= ssGetOutputPortWidth(S,1);
		int_T				y2Width	= ssGetOutputPortWidth(S,2);

		int_T				u0Width	= ssGetInputPortWidth(S,0);
		int_T				u1Width	= ssGetInputPortWidth(S,1);
		InputRealPtrsType 	u1_U 	= ssGetInputPortRealSignalPtrs(S,1);


		// Initializing work pointers to NULL pointer:
		// P-work vectors
		//0 MsgList
		//1 tick
		//2 System
		MsgList				= NULL;
		tick				= NULL;
		System				= NULL;
		Xc_dot				= NULL;

		MsgList = (sMsgList *)malloc((MAX_MSG+1)*sizeof(sMsgList));
		if(!MsgList) {
			FAIL_FLAG = -3;		// Memory allocation problems
			ssSetErrorStatus(S, "sDynSys_main-mdlStart: Memory request failed for \'MsgList\' struct in mdlStart.\n");
		}
		else {
			MsgList[0].msg_length = -1;
			ssSetPWorkValue(S,0,MsgList);
			
			strcpy_s(logfile,MAX_MSG_CHAR+1,sfunc);
			strcat_s(logfile,MAX_MSG_CHAR+1,".log");
			status = fopen_s(&fp,logfile,"w+");
			if( status != 0 ) {
				sprintf_s(msg,MAX_LINE_SIZE,"The log file %s could not be opened.",logfile);
				ssSetErrorStatus(S, msg);
			}
			else {
				fprintf(fp,">>>>>>>>>>>>>>>>>>>>>> LOGFILE <<<<<<<<<<<<<<<<<<<<<<<<\n");
				fprintf(fp,">>>> %s/%s \n\n",sfunc_path,sfunc);
				fclose(fp);    
			}
			/* End of message function initialization */

			mxGetString(InitFname_PARAM(S),filename,nu+1);  // Getting input filename 
		}

		if(FAIL_FLAG > 0) {
			tick = (double *)malloc(3*sizeof(double));  
			if(!tick) {
				FAIL_FLAG = -3;		// Memory allocation problems
				sprintf_s(msg,MAX_LINE_SIZE,"Memory request failed for tick.");
				SendMsg(S,ERROR_LOG,msg,"mdlStart");
			}
			else {
				tick[0]	= ticker(S);
				tick[1]	= tick[0];
				tick[2]	= tick[0];
			}
		}

		if(FAIL_FLAG > 0) {
			/* Assigning system data to sSystemData struct *System */
			System = (sSystemCFG *)malloc(sizeof(sSystemCFG));
			if(!System) {
				FAIL_FLAG = -3;		// Memory allocation problems
				sprintf_s(msg,MAX_LINE_SIZE,"Memory request failed for \'System\' struct in mdlStart.");
				SendMsg(S,ERROR_LOG,msg,"mdlStart");
				System = NULL;
			}
			else {
				// Reading in configuration data:
				FAIL_FLAG = Param2Struct(S, System, filename);

				if(FAIL_FLAG > 0) {
					/* Assigning an Xc_dot vector */
					Xc_dot = (real_T *)malloc(NumContStates*sizeof(real_T));
					if(!Xc_dot) {
						FAIL_FLAG = -3;		// Memory allocation problems
						sprintf_s(msg,MAX_LINE_SIZE,"Memory request failed for \'Xc_dot\' struct in mdlStart.");
						SendMsg(S,ERROR_LOG,msg,"mdlStart");
						Xc_dot = NULL;
					}
					else {
						// Reading in configuration data:
						ModelXcDerivatives(S, System, Xc_dot, Xc, Xd, (real_T *)u1_U[0]);
					}
				}
			}
		}


		if(FAIL_FLAG > 0) {
			if(System->NumContStates != NumContStates) {
				FAIL_FLAG = -1;
				sprintf_s(msg,MAX_LINE_SIZE,"The number of continuous states (Xc) assigned in initialization file is %d and does not match the number as specified by mask parameter 0, element 1, i.e. NumContStates = %d.",System->NumContStates, NumContStates);
				SendMsg(S,ERROR_LOG,msg,"mdlStart");
			}
			else if(System->NumDiscStates != NumDiscStates) {
				FAIL_FLAG = -1;
				sprintf_s(msg,MAX_LINE_SIZE,"The number of discrete states (Xd) assigned in initialization file is %d and does not match the number as specified by mask parameter 0, element 2, i.e. NumDiscStates = %d.",System->NumDiscStates, NumDiscStates);
				SendMsg(S,ERROR_LOG,msg,"mdlStart");
			}
			else if(System->NumInputSignals != NumInputSignals) {
				FAIL_FLAG = -1;
				sprintf_s(msg,MAX_LINE_SIZE,"The number of control input signals (u) assigned in initialization file is %d and does not match the number as specified by mask parameter 0, element 3, i.e. NumInputSignals = %d.",System->NumInputSignals, NumInputSignals);
				SendMsg(S,ERROR_LOG,msg,"mdlStart");
			}
			else if(System->NumOutputSignals != NumOutputSignals) {
				FAIL_FLAG = -1;
				sprintf_s(msg,MAX_LINE_SIZE,"The number of output signals (y) assigned in initialization file is %d and does not match the number as specified by mask parameter 0, element 4, i.e. NumOutputSignals = %d.",System->NumOutputSignals, NumOutputSignals);
				SendMsg(S,ERROR_LOG,msg,"mdlStart");
			}
			else if(System->NumParameters != NumParameters) {
				FAIL_FLAG = -1;
				sprintf_s(msg,MAX_LINE_SIZE,"The number of parameters assigned in initialization file is %d and does not match the number as specified by mask parameter 0, element 5, i.e. NumParameters = %d.",System->NumParameters, NumParameters);
				SendMsg(S,ERROR_LOG,msg,"mdlStart");
			}
		}

		if(FAIL_FLAG > 0) {
			sprintf_s(msg,MAX_LINE_SIZE,"Successful initialization of S-function in mdlStart.");
			SendMsg(S,DEBUG_LOG,msg,"mdlStart");
		}
		else if(FAIL_FLAG == 0) {
			sprintf_s(msg,MAX_LINE_SIZE,"No fields were assigned from initialization data file.");
			SendMsg(S,ERROR_LOG,msg,"mdlStart");
		}
		else if(FAIL_FLAG == -1) {
			sprintf_s(msg,MAX_LINE_SIZE,"Error in parameters from initialization file.");
			SendMsg(S,ERROR_LOG,msg,"mdlStart");
		}
		else if(FAIL_FLAG == -2) {
			sprintf_s(msg,MAX_LINE_SIZE,"Cannot open initialization file.");
			SendMsg(S,ERROR_LOG,msg,"mdlStart");
		}
		else if(FAIL_FLAG == -3) {
			sprintf_s(msg,MAX_LINE_SIZE,"Memory allocation problems in S-function.");
			SendMsg(S,ERROR_LOG,msg,"mdlStart");
		}
		else if(FAIL_FLAG == -4) {
			sprintf_s(msg,MAX_LINE_SIZE,"Error in dynamic sizing of input/output ports for S-function.");
			SendMsg(S,ERROR_LOG,msg,"mdlStart");
		}
		else {
			sprintf_s(msg,MAX_LINE_SIZE,"Unknown configuration error of data structures. Check input parameters and initialization file.");
			SendMsg(S,ERROR_LOG,msg,"mdlStart");
		}
		// I-work vectors
		//0 FAIL_FLAG
		ssSetIWorkValue(S,0,FAIL_FLAG);

		// R-work vectors
		//0 Xc_dot
		if(!Xc_dot) {
			for(i=0;i<NumContStates;i++) 
				ssSetRWorkValue(S, (int_T)i, 0.0);
		}
		else {
			for(i=0;i<NumContStates;i++) 
				ssSetRWorkValue(S, (int_T)i, Xc_dot[i]);
		}
		free(Xc_dot);
		// P-work vectors
		//0 MsgList
		//1 tick
		//2 System
		ssSetPWorkValue(S,1,tick);
 		ssSetPWorkValue(S,2,System);

	}
#endif /*  MDL_START */



/* Function: mdlOutputs ======================================================= 
 */

static void mdlOutputs(SimStruct *S, int_T tid)
{
	int			i, InfoFlag = 0; 
	int			FAIL_FLAG	= ssGetIWorkValue(S,0);
	double		*tack, tick[20];
	char		msg[15*MAX_LINE_SIZE], str[MAX_LINE_SIZE];

	sSystemCFG	*System;

    real_T				*Xd			= ssGetDiscStates(S);
    real_T            	*Xc       	= ssGetContStates(S);
    real_T				Xc_dot;

    real_T            	*y0_Y		= ssGetOutputPortRealSignal(S,0);	
	int_T				y0Width		= ssGetOutputPortWidth(S,0);
    real_T            	*y1_Xd		= ssGetOutputPortRealSignal(S,1);	
	int_T				y1Width		= ssGetOutputPortWidth(S,1);
	real_T            	*y2_Xc		= ssGetOutputPortRealSignal(S,2);	
	int_T				y2Width		= ssGetOutputPortWidth(S,2);
	real_T            	*y3_XcDot	= ssGetOutputPortRealSignal(S,3);	
	int_T				y3Width		= ssGetOutputPortWidth(S,3);

    InputRealPtrsType 	u0_Control	= ssGetInputPortRealSignalPtrs(S,0);
	int_T				u0Width		= ssGetInputPortWidth(S,0);
    InputRealPtrsType 	u1_U 		= ssGetInputPortRealSignalPtrs(S,1);
	int_T				u1Width		= ssGetInputPortWidth(S,1);

	const char_T		*sfunc		= FILE_NAME;//ssGetModelName(S);
	const char_T		*sfunc_path	= ssGetPath(S);
	sMsgList			*MsgList;

	// Initializing work pointers:
	// P-work vectors
	//0 MsgList
	//1 tick
	//2 System
	MsgList				= ssGetPWorkValue(S,0);
	System				= ssGetPWorkValue(S,2);
	
	/* Calculating and setting outputs */
	if(FAIL_FLAG == 1) { // Continuing only if configuration of system structure has been successful 

		if(((int)(*u0_Control[0]) > 99)) tick[0] = ticker(S);

		if(System->NumOutputSignals > 0)
			ModelOutputs(S, System, Xc, Xd, (real_T *)u1_U[0], y0_Y);
		else {
			for(i=0;i<y0Width;i++)
				y0_Y[i]	 = 0.0;		
		}

		if(((int)(*u0_Control[0]) > 99)) tick[1] = ticker(S);

		/* Assigning outputs y0 - y1 */
		/* ========================= */
		if(System->NumDiscStates > 0)
			for(i=0;i<y1Width;i++)
				y1_Xd[i] = Xd[i];
		else {
			for(i=0;i<y1Width;i++)
				y1_Xd[i] = 0.0;
		}
		if(System->NumContStates > 0)
			for(i=0;i<y2Width;i++) {
				y2_Xc[i]	= Xc[i];
				y3_XcDot[i] = ssGetRWorkValue(S, (int_T)i);;
			}
		else {
			for(i=0;i<y2Width;i++) {
				y2_Xc[i]	= 0.0;
				y3_XcDot[i] = 0.0;
			}
		}

		if(((int)(*u0_Control[0]) > 99)) tick[2] = ticker(S);

		///* Debugging reporting. */
		//InfoFlag = (int)(*u0_Control[0]);
		//if((InfoFlag != (int)(x[0])) && ((InfoFlag > 0) && (InfoFlag <= 5))) {
		//	DebugReport(S, System, Location, Pipe, Input, TouchDetect, ProximityDetect, &InfoFlag);
		//}

		// Timing reporting
		if(((int)(*u0_Control[0]) > 99)) {
			// P-work vectors
			//1 tick
			tack	 = ssGetPWorkValue(S,1);
			tick[19] = ticker(S);
		
			if(((int)(*u0_Control[0]) == 100)) {
				if ssIsMajorTimeStep(S) {
					sprintf_s(str,MAX_LINE_SIZE,"mdlOutput - MAJOR timestep: Elapsed CPU time: %8.4f msec. since last. Execution time: %8.4f msec.\n",tick[19]-tack[0],tick[19]-tick[0]);
					printf("%s",str);
				}
				else {
					sprintf_s(str,MAX_LINE_SIZE,"mdlOutput - MINOR timestep: Elapsed CPU time: %8.4f msec. since last. Execution time: %8.4f msec.\n",tick[19]-tack[0],tick[19]-tick[0]);
					printf("%s",str);
				}
			}
			else if(((int)(*u0_Control[0]) == 102)) {
				sprintf_s(str,MAX_LINE_SIZE,"mdlOutput: Elapsed CPU time: %8.6f msec. since last. Execution time: %8.6f msec.\n",tick[19]-tack[0],tick[19]-tick[0]);
				printf("%s",str); 

				sprintf_s(str,MAX_LINE_SIZE,"   1- 0: %8.6f msec. - Calculating outputs y0\n",											tick[1] -tick[0]);	printf("%s",str); 
				sprintf_s(str,MAX_LINE_SIZE,"   2- 1: %8.6f msec. - Assigning outputs y1 - y3\n",										tick[2] -tick[1]);	printf("%s",str); 
				sprintf_s(str,MAX_LINE_SIZE,"  19- 2: %8.6f msec. - Debug reporting; DebugReport()\n",									tick[19]-tick[2]);	printf("%s",str); 

			}
			tack[0] = tick[19];
		}

		// Decreasing Message quarantine status of all messages in list.
		for(i = 0; i<=MAX_MSG; i++) {
			if(MsgList[i].msg_length < 0) {
				break;
			}
			if(MsgList[i].status > 0) {
				if(MsgList[i].status == 1) { 
					sprintf_s(msg,15*MAX_LINE_SIZE,"%d %s/%s", -i-1, sfunc_path,sfunc);	// Contra message for message going low. Last in mdlOutpouts.
					printf(msg);
				}
				MsgList[i].status = MsgList[i].status - 1;
			}
		}
	}
	else if(FAIL_FLAG != -3) {
		sprintf_s(str,MAX_LINE_SIZE,"Error in initialization of configuration data (Code %d). Sending fail-safe values to the outputs.", FAIL_FLAG);
		SendMsg(S,INFO_LOG,str,"mdlOutputs");
		/* Assigning outputs y0 - y2 */
		for(i=0;i<y0Width;i++)
			y0_Y[i]	 = 0.0;
		for(i=0;i<y1Width;i++)
			y1_Xd[i] = 0.0;
		for(i=0;i<y2Width;i++)
			y2_Xc[i] = 0.0;
		for(i=0;i<y3Width;i++)
			y3_XcDot[i] = 0.0;
	}
	else {
		sprintf_s(str,MAX_LINE_SIZE,"sDynSys-mdlOutputs: Error in execution of S-function (Code %d). Sending fail-safe values to the outputs.", FAIL_FLAG);
		ssSetErrorStatus(S, str);
		/* Assigning outputs y0 - y2 */
		for(i=0;i<y0Width;i++)
			y0_Y[i]	 = 0.0;
		for(i=0;i<y1Width;i++)
			y1_Xd[i] = 0.0;
		for(i=0;i<y2Width;i++)
			y2_Xc[i] = 0.0;
		for(i=0;i<y3Width;i++)
			y3_XcDot[i] = 0.0;
	}
}



#define MDL_DERIVATIVES
/* Function: mdlDerivatives =================================================
 */
static void mdlDerivatives(SimStruct *S)
{
	int			i, FAIL_FLAG	= ssGetIWorkValue(S,0);
	sSystemCFG	*System;
		
    real_T				*Xd		= ssGetDiscStates(S);
    real_T				*Xc     = ssGetContStates(S);
    real_T				*Xc_dot	= ssGetdX(S);
	int_T				numXc	= ssGetNumContStates(S);

    InputRealPtrsType 	u1_U 	= ssGetInputPortRealSignalPtrs(S,1);
	int_T				u1Width	= ssGetInputPortWidth(S,1);

	// Initializing work pointers:
	// P-work vectors
	//0 MsgList
	//1 tick
	//2 System
	System				= ssGetPWorkValue(S,2);
	
	/* Calculating and setting derivatives */
	if((FAIL_FLAG == 1) && (System->NumContStates > 0)) { 
		ModelXcDerivatives(S, System, Xc_dot, Xc, Xd, (real_T *)u1_U[0]);
	}
	else {
		for(i=0;i<numXc;i++)
			Xc_dot[i] = 0.0;
	}

	// R-work vectors
	//0 Xc_dot
	for(i=0;i<numXc;i++) 
		ssSetRWorkValue(S, (int_T)i, Xc_dot[i]);

}
 
#define MDL_UPDATE
/* Function: mdlUpdate ======================================================
 * Abstract:
 *      xdot = Ax + Bu
 */
static void mdlUpdate(SimStruct *S, int_T tid)
{
	int			FAIL_FLAG	= ssGetIWorkValue(S,0);
	sSystemCFG	*System;
		
    real_T				*Xd		= ssGetDiscStates(S);
    real_T				*Xc     = ssGetContStates(S);
	int_T				numXd	= ssGetNumDiscStates(S);

    InputRealPtrsType 	u1_U 	= ssGetInputPortRealSignalPtrs(S,1);
	int_T				u1Width	= ssGetInputPortWidth(S,1);

	// Initializing work pointers:
	// P-work vectors
	//0 MsgList
	//1 tick
	//2 System
	System				= ssGetPWorkValue(S,2);
	
	/* Calculating and setting derivatives */
	if((FAIL_FLAG == 1) && (System->NumDiscStates > 0)) { 
		ModelXdUpdate(S, System, Xc, Xd, (real_T *)u1_U[0]);
	}
	else {
		;	// No update
	}
}



 
/* Function: mdlTerminate =====================================================
 * Abstract:
 *    No termination needed, but we are required to have this routine.
 */
static void mdlTerminate(SimStruct *S)
{
	int				i, num; 
	void**			ptr;
	sSystemCFG		*System;

	System			= ssGetPWorkValue(S,2);

	num	= ssGetNumPWork(S);
	ptr = ssGetPWork(S);
	for(i=0;i<num;i++) {
		free(ptr[i]);
	}
}

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
