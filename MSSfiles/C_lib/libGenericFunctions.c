/* 
    libGenericFunctions.c

	Library file for generic simulation functions related to the sDynSys s-function.

	References:
	[1]	W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery, Numerical Recipes in C - The Art of Scientific Computing, 3rd ed. 
		Cambridge University Press, New York, USA, 2007.

	Copyright: 		Roger Skjetne, NTNU
    Date created: 	2011.04.13 Roger Skjetne
    Revised:      	

*/


#include "libGenericFunctions.h"
//#include "../../Lib/nr.h"

int		Param2Struct(SimStruct *S, sSystemCFG *SystemStruct, char *filename)
/* Function that assigns parameter configuration data to data structure SystemStruct. Returns:
	  1: successful
	  0: no fields were assigned
	 -1: error in configuration (see screen printout)
	 -2: filename was not found. 
*/
{
	char	Buff[MAX_LINE_SIZE+1];
	char	*s00, *s01, *ptr;
	char	s10[MAX_LINE_SIZE+1], s11[MAX_STRING_SIZE+1], msg[MAX_LINE_SIZE+1];
    FILE	*fp;

	int		i, j, k=-1, k_old=-1, p=-1, p_old=-1, d=0, m, MaxVectorSize = 0;
	size_t	h1, h2;
	int		N=0, M=0, num, count=0, status, n1=MAX_PARAMS, n2=0;
	int		NumContStates=0, NumDiscStates = 0, NumInputSignals=0, NumOutputSignals=0, NumParameters=0;
	double	fnum;
	

    /* Open file for input */
    fp = fopen(filename,"r");
	if( fp == NULL )
	{
        status = -2;
		sprintf(msg,"The parameter configuration file %s could not be opened.",filename);
		SendMsg(S,DEBUG_LOG,msg,"Param2Struct");
		SystemStruct = NULL;
	}
	else
	{
		while(fgets(Buff,sizeof(Buff),fp) != NULL) {
			M++;
			if(M > 2*MAX_PARAMS*(PARAM_FIELDS+2)+100) break; // Exiting while loop if iterations exceed threshold.
			strcpy(s10,Buff);
			ptr = strtok(s10, " ");
			if(strncmp(ptr, "%",1)) {
				for(i=0;i<sizeof(Buff);i++){
					s10[i] = tolower(Buff[i]);
				}

				s00 = strstr(s10, "system.param");
				if(s00 != NULL) {
					d++; // Number of hits on "System.Param" in text file.
					h1 = strcspn(Buff, "{");	// Finding where index begins in Buff string
					s01 = Buff+h1+1;
					h2 = strcspn(s01, "}");		// Finding where index ends in Buff string
					for(i=0; i<h2; i++) {
						s11[i] = Buff[i+h1+1];
					}
					s11[i] = '\0';	// String s11 = <index digits>
					k = atoi(s11);	// Assigning integer k to index in parameter structure
					if(k<=n1) n1 = k;	// Finding start index n1
					if(k>=n2) n2 = k;	// Finding end index n2

					if(strncmp(s01+h2+1, ".tag", 4)) {	// s01+h2+1 gives address to where ".tag" begins
						if(strncmp(s01+h2+1, ".spare1", 7)) {	// s01+h2+1 gives address to where ".spare1" begins
							if(strncmp(s01+h2+1, ".spare2", 7)) {	// s01+h2+1 gives address to where ".spare2" begins
								if(strncmp(s01+h2+1, ".spare3", 7)) {	// s01+h2+1 gives address to where ".spare3" begins
									if(strncmp(s01+h2+1, ".spare4", 7)) {	// s01+h2+1 gives address to where ".spare4" begins
										if(strncmp(s01+h2+1, ".size", 5)) {	// s01+h2+1 gives address to where ".size" begins
											if(strncmp(s01+h2+1, ".values", 7)) {	// s01+h2+1 gives address to where ".values" begins
												sprintf(msg,"Unknown field <%s> in data struct System.Param{%d} in file: %s.",s01+h2+1,k,filename);
												SendMsg(S,FILE_LOG,msg,"Param2Struct");
											}
											else {
												// Values:
												if(SystemStruct->Param[k].size > 0){
													j = 0;
													h1 = strcspn(Buff, "[");
													s01 = Buff+h1+1;
													h2 = strcspn(s01, "]")+1;
													for(m=0; m<SystemStruct->Param[k].size;m++) {
														for(i=0; i<h2; i++) {
															if(isdigit(Buff[i+h1+1]))	
																s11[j++] = Buff[i+h1+1];
															else if(Buff[i+h1+1] == '.')
																s11[j++] = Buff[i+h1+1];
															else if(Buff[i+h1+1] == '+')
																s11[j++] = Buff[i+h1+1];
															else if(Buff[i+h1+1] == '-')
																s11[j++] = Buff[i+h1+1];
															else if(Buff[i+h1+1] == 'e')
																s11[j++] = Buff[i+h1+1];
															else {
																s11[j] = '\0';	// String s11 = <n digits>
																fnum = (float)atof(s11); 
																SystemStruct->Param[k].pvec[m] = fnum;
																count++; // Structure field assignment made, count incremented.
																sprintf(msg,"SystemStruct->Param[%d].pvec[%d] = %.3f",k,m, SystemStruct->Param[k].pvec[m]);
																SendMsg(S,FILE_LOG,msg,"Param2Struct");
																h1 = i+h1+1;
																s01 = Buff+h1+1;
																h2 = strcspn(s01, "]")+1;
																j = 0;
																break;
															}
														}
													}
												}
											}
										}
										else {
											j = 0;
											h1 = strcspn(Buff, "=");
											s01 = Buff+h1+1;
											h2 = strcspn(s01, ";");
											for(i=0; i<h2; i++) {
												if(isdigit(Buff[i+h1+1])) 
													s11[j++] = Buff[i+h1+1];
												else if(Buff[i+h1+1] == '+')
													s11[j++] = Buff[i+h1+1];
												else if(Buff[i+h1+1] == '-')
													s11[j++] = Buff[i+h1+1];
												else
													s11[j++] = ' ';
											}
											s11[j] = '\0';	// String s11 = <size of parameter vector>
											num = atoi(s11);
											if(num > MaxVectorSize) MaxVectorSize = num;
											SystemStruct->Param[k].size = num;
											sprintf(msg,"SystemStruct->Param[%d].size = %d.",k, SystemStruct->Param[k].size);
											SendMsg(S,FILE_LOG,msg,"Param2Struct");
											count++; // Structure field assignment made, count incremented.
										}
									}
									else {
										; // spare4 is available
									}
								}
								else {
									; // spare3 is available
								}
							}
							else {
								; // spare2 is available
							}
						}
						else {
							; // spare1 is available
						}
					}
					else {
						h1 = strcspn(Buff, "\'");
						s01 = Buff+h1+1;
						h2 = strcspn(s01, "\'");
						for(i=0; i<h2; i++) {
							s11[i] = Buff[i+h1+1];
						}
						s11[i] = '\0';	// String s11 = <tag of equipment>
						strcpy(SystemStruct->Param[k].tag, s11);
						sprintf(msg,"SystemStruct->Param[%d].tag = \'%s\'",k, SystemStruct->Param[k].tag);
						SendMsg(S,FILE_LOG,msg,"Param2Struct");
						count++; // Structure field assignment made, count incremented.
					}
					if(k!=k_old) N++;
					k_old = k;
				}
				else {
					s00 = strstr(s10, "system.numcontstates");
					if(s00 != NULL) {
						s01 = strstr(s00, "=");
						if(s01 != NULL) {
							j = 0;
							h1 = strcspn(Buff, "=");
							s01 = Buff+h1+1;
							h2 = strcspn(s01, ";");
							for(i=0; i<h2; i++) {
								if(isdigit(Buff[i+h1+1])) 
									s11[j++] = Buff[i+h1+1];
								else if(Buff[i+h1+1] == '+')
									s11[j++] = Buff[i+h1+1];
								else if(Buff[i+h1+1] == '-')
									s11[j++] = Buff[i+h1+1];
								else
									s11[j++] = ' ';
								}
							s11[j] = '\0';	// String s11 = <number of continuous states>
							num = atoi(s11);
							SystemStruct->NumContStates = num;
							sprintf(msg,"SystemStruct->NumContStates = %d",SystemStruct->NumContStates);
							SendMsg(S,FILE_LOG,msg,"Param2Struct");
							count++; // Structure field assignment made, count incremented.
						}
					}
					else {
						s00 = strstr(s10, "system.numdiscstates");
						if(s00 != NULL) {
							s01 = strstr(s00, "=");
							if(s01 != NULL) {
								j = 0;
								h1 = strcspn(Buff, "=");
								s01 = Buff+h1+1;
								h2 = strcspn(s01, ";");
								for(i=0; i<h2; i++) {
									if(isdigit(Buff[i+h1+1])) 
										s11[j++] = Buff[i+h1+1];
									else if(Buff[i+h1+1] == '+')
										s11[j++] = Buff[i+h1+1];
									else if(Buff[i+h1+1] == '-')
										s11[j++] = Buff[i+h1+1];
									else
										s11[j++] = ' ';
									}
								s11[j] = '\0';	// String s11 = <number of discrete states>
								num = atoi(s11);
								SystemStruct->NumDiscStates = num;
								sprintf(msg,"SystemStruct->NumDiscStates = %d",SystemStruct->NumDiscStates);
								SendMsg(S,FILE_LOG,msg,"Param2Struct");
								count++; // Structure field assignment made, count incremented.
							}
						}
						else {
							s00 = strstr(s10, "system.numinputsignals");
							if(s00 != NULL) {
								s01 = strstr(s00, "=");
								if(s01 != NULL) {
									j = 0;
									h1 = strcspn(Buff, "=");
									s01 = Buff+h1+1;
									h2 = strcspn(s01, ";");
									for(i=0; i<h2; i++) {
										if(isdigit(Buff[i+h1+1])) 
											s11[j++] = Buff[i+h1+1];
										else if(Buff[i+h1+1] == '+')
											s11[j++] = Buff[i+h1+1];
										else if(Buff[i+h1+1] == '-')
											s11[j++] = Buff[i+h1+1];
										else
											s11[j++] = ' ';
										}
									s11[j] = '\0';	// String s11 = <number of control input signals>
									num = atoi(s11);
									SystemStruct->NumInputSignals = num;
									sprintf(msg,"SystemStruct->NumInputSignals = %d",SystemStruct->NumInputSignals);
									SendMsg(S,FILE_LOG,msg,"Param2Struct");
									count++; // Structure field assignment made, count incremented.
								}
							}
							else {
								s00 = strstr(s10, "system.numoutputsignals");
								if(s00 != NULL) {
									s01 = strstr(s00, "=");
									if(s01 != NULL) {
										j = 0;
										h1 = strcspn(Buff, "=");
										s01 = Buff+h1+1;
										h2 = strcspn(s01, ";");
										for(i=0; i<h2; i++) {
											if(isdigit(Buff[i+h1+1])) 
												s11[j++] = Buff[i+h1+1];
											else if(Buff[i+h1+1] == '+')
												s11[j++] = Buff[i+h1+1];
											else if(Buff[i+h1+1] == '-')
												s11[j++] = Buff[i+h1+1];
											else
												s11[j++] = ' ';
											}
										s11[j] = '\0';	// String s11 = <number of measured output signals>
										num = atoi(s11);
										SystemStruct->NumOutputSignals = num;
										sprintf(msg,"SystemStruct->NumOutputSignals = %d",SystemStruct->NumOutputSignals);
										SendMsg(S,FILE_LOG,msg,"Param2Struct");
										count++; // Structure field assignment made, count incremented.
									}
								}
								else {
									s00 = strstr(s10, "system.numparameters");
									if(s00 != NULL) {
										s01 = strstr(s00, "=");
										if(s01 != NULL) {
											j = 0;
											h1 = strcspn(Buff, "=");
											s01 = Buff+h1+1;
											h2 = strcspn(s01, ";");
											for(i=0; i<h2; i++) {
												if(isdigit(Buff[i+h1+1])) 
													s11[j++] = Buff[i+h1+1];
												else if(Buff[i+h1+1] == '+')
													s11[j++] = Buff[i+h1+1];
												else if(Buff[i+h1+1] == '-')
													s11[j++] = Buff[i+h1+1];
												else
													s11[j++] = ' ';
												}
											s11[j] = '\0';	// String s11 = <number of system parameters>
											num = atoi(s11);
											SystemStruct->NumParameters = num;
											sprintf(msg,"SystemStruct->NumParameters = %d",SystemStruct->NumParameters);
											SendMsg(S,FILE_LOG,msg,"Param2Struct");
											count++; // Structure field assignment made, count incremented.
										}
									}
									else {
										s00 = strstr(s10, "system.sampletime");
										if(s00 != NULL) {
											s01 = strstr(s00, "=");
											if(s01 != NULL) {
												j = 0;
												h1 = strcspn(Buff, "=");
												s01 = Buff+h1+1;
												h2 = strcspn(s01, ";");
												for(i=0; i<h2; i++) {
													if(isdigit(Buff[i+h1+1])) 
														s11[j++] = Buff[i+h1+1];
													else if(Buff[i+h1+1] == '.')
														s11[j++] = Buff[i+h1+1];
													else if(Buff[i+h1+1] == '+')
														s11[j++] = Buff[i+h1+1];
													else if(Buff[i+h1+1] == '-')
														s11[j++] = Buff[i+h1+1];
													else if(Buff[i+h1+1] == 'e')
														s11[j++] = Buff[i+h1+1];
												}
												s11[j] = '\0';	// String s11 = <sample time>
												fnum = atof(s11);
												SystemStruct->SampleTime = fnum;
												sprintf(msg,"SystemStruct->SampleTime = %.3f",SystemStruct->SampleTime);
												SendMsg(S,FILE_LOG,msg,"Param2Struct");
												count++; // Structure field assignment made, count incremented.
											}
										}
									}
								}
							}
						}
					}
				}
			}
			else {
				//sprintf(msg,"Commented textline in input file %s, line %d.",filename,M);
				//SendMsg(S,DEBUG_LOG,msg,"Param2Struct");
			}
		}
		fclose(fp);    
		sprintf(msg,"There were %d hits on \'System.Param\' in file: %s.",d,filename);
		SendMsg(S,FILE_LOG,msg,"Param2Struct");
		sprintf(msg,"There were %d assignments made to parameter data structure.",count);
		SendMsg(S,FILE_LOG,msg,"Param2Struct");  

		if(N == 0) {
			n1		= 0; 
			n2		= 0;
			status	= 0;
			sprintf(msg,"None parameters configured.");
			SendMsg(S,DEBUG_LOG,msg,"Param2Struct");
		}
		else if(n1 != 0) {
			status = -1;
			sprintf(msg,"Parameter data must be configured sequentially, starting with index i = 0 in input file: %s",filename);
			SendMsg(S,DEBUG_LOG,msg,"Param2Struct");
		}
		else if(N != n2+1) {
			status = -1;
			sprintf(msg,"Parameter data must be configured sequentially, in input file: %s",filename);
			SendMsg(S,DEBUG_LOG,msg,"Param2Struct");
		}
		else if(N != SystemStruct->NumParameters) {
			status = -1;
			sprintf(msg,"Configuration error: Number of parameters %d in input file: %s does not match the number of parameter assignments %d in the data structure.",SystemStruct->NumParameters,filename,N);
			SendMsg(S,DEBUG_LOG,msg,"Param2Struct");
		}
		else if(SystemStruct->NumParameters > MAX_PARAMS) {
			status = -1;
			sprintf(msg,"Configuration error: Number of parameters %d in input file: %s exceeds the macro MAX_PARAMS = %d in s-function.",SystemStruct->NumParameters,filename,MAX_PARAMS);
			SendMsg(S,DEBUG_LOG,msg,"Param2Struct");
		}
		else if(MaxVectorSize > MAX_PARAM_SIZE) {
			status = -1;
			sprintf(msg,"Configuration error: Maximum number of elements assigned to a parameter vecor is counted to %d in input file: %s, and this exceeds the macro MAX_PARAM_SIZE = %d in s-function.",MaxVectorSize,filename,MAX_PARAM_SIZE);
			SendMsg(S,DEBUG_LOG,msg,"Param2Struct");
		}
		else {
			status = 1;
			sprintf(msg,"Successful configuration of parameters from input file: %s.\n\n",filename);
			SendMsg(S,FILE_LOG,msg,"Param2Struct");
		}
	}

	return status;	
}


void	SendMsg(SimStruct *S, int level, char *msg, char *subfunc)
{
	int				i, SendIdx = -1;
	char			buf[4*MAX_STRING_SIZE+1], string[4*MAX_STRING_SIZE+1], sub_func[MAX_STRING_SIZE+1]; 
	char			timestr[22], logfile[MAX_MSG_CHAR+1];
	const char_T	*sfunc		= ssGetModelName(S);
	//const char_T	*sfunc		= FILE_NAME; 
	const char_T	*sfunc_path	= ssGetPath(S);
	//char			sfunc[20]	= "Sim_msg";
	//char			sfunc_path[20] = "path";
    FILE			*fp;
	struct tm		*ltime;
	time_t			t;
	sMsgList		*MsgList;


	strcpy(sub_func,subfunc);
	MsgList	= ssGetPWorkValue(S,0);

	strcpy(logfile,sfunc);
	strcat(logfile,".log");
	fp = fopen(logfile,"a+");
	if( fp == NULL ) {
		level = DEBUG_LOG;
		sprintf_s(msg,MAX_STRING_SIZE+1,"The log file %s could not be opened.",logfile);
		strcpy(sub_func,"SendMsg");
	}
	t = time(NULL);
	ltime = localtime(&t);
	strftime(timestr,21,"%Y-%m-%d %H:%M:%S",ltime);

	if(level == FILE_LOG && fp != NULL) {
		fprintf(fp,"%s %s - %s: %s\n",timestr,sfunc_path,sub_func,msg);
	}
	if((level != FILE_LOG)) {
		for(i = 0; i <= MAX_MSG; i++) {
			if(i >= MAX_MSG) {
				level = INFO_LOG;
				sprintf_s(msg,MAX_STRING_SIZE+1,"Message buffer of %s/%s has reached maximum entries %d.",sfunc_path,sfunc,MAX_MSG);
				strcpy(sub_func,"SendMsg");
				SendIdx = MAX_MSG;
			}
			if(MsgList[i].msg_length < 0) {
				strcpy(MsgList[i].subfunc, sub_func);
				strcpy(MsgList[i].msg, msg);
				MsgList[i].msg_length	= (int)strlen(msg);
				MsgList[i].level		= level;
				MsgList[i+1].msg_length	= -1;
				SendIdx = i;
				break;
			}
			else {
				if(strcmp(sub_func,MsgList[i].subfunc)) {
					// No hit.
				}
				else {	// Hit on subfunction
					if(strncmp(msg,MsgList[i].msg,MsgList[i].msg_length)) {
						// No hit.
					}
					else {	// Hit on message
						if(MsgList[i].status == 0) {
							SendIdx = i;
						}
						else if(MsgList[i].level > level) {
							SendIdx = i;
						}
						else {
							MsgList[i].status = MAX_STATUS; // Set active again
						}
						break;
					}
				}
			}
		}
		if(SendIdx >= 0) {
			sprintf_s(string,4*MAX_STRING_SIZE+1," %d %s/%s-%s: %s",SendIdx+1,sfunc_path,sfunc,sub_func,msg);

			/* Generate the message */
			switch (level) {
			   case ERROR_LOG:  /* Reports a critical error */
					sprintf_s(buf,4*MAX_STRING_SIZE+1,"ERROR (%s): %s \n",timestr,string);
					break;
			   case INFO_LOG:   /* Reports general information */
					sprintf_s(buf,4*MAX_STRING_SIZE+1,"INFO  (%s): %s \n",timestr,string);   
					break;
			   case EVENT_LOG:   /* Reports a simulation event */
					sprintf_s(buf,4*MAX_STRING_SIZE+1,"EVENT (%s): %s \n",timestr,string);   
					break;
			   case DEBUG_LOG:	/* Reports detailed debugging information */
					sprintf_s(buf,4*MAX_STRING_SIZE+1,"DEBUG (%s): %s \n",timestr,string);
					break;            
			   case FILE_LOG:	/* Reports message only to file */
					break;            
			}
		    
			/* Report to standard output (display) */
			printf(buf);
		    
			if( fp != NULL ) {
				fprintf(fp,buf);
			}
			MsgList[SendIdx].level = level;
			MsgList[SendIdx].status = MAX_STATUS;
			SendIdx = -1;

			if(level == ERROR_LOG)
				ssSetErrorStatus(S, buf);
		}
	}
	if( fp != NULL ) {
		fclose(fp);    
	}
}


double	ticker(SimStruct *S)
/* Function that returns a double giving point in time in milliseconds. */
{  
	char msg[MAX_LINE_SIZE+1];
	LARGE_INTEGER ticksPerSecond;
	LARGE_INTEGER tick;   // A point in time

	// Get the high resolution counter's accuracy
	if (!QueryPerformanceFrequency(&ticksPerSecond)) {
		sprintf(msg,"No go - QueryPerformance not present.");
		SendMsg(S,DEBUG_LOG,msg,"ticker");
		return -1.0; 
	}

	// what time is it?
	if (!QueryPerformanceCounter(&tick) ){ 
		sprintf(msg,"No go - Counter not installed.");
		SendMsg(S,DEBUG_LOG,msg,"ticker");
		return -1.0; 
	}

	return 1000*((double)tick.QuadPart/(double)ticksPerSecond.QuadPart);
}
int		comp(const void *i, const void *j)
/* Compares the integers (used with function qsort) */
{
	return *(int *)i - *(int *)j;
}


double	rad2pipi(double x)
{
	double y;

	if(x < 0.0) {
		y = fmod(x-pi,2*pi);
		y = y + pi;
	}
	else {
		y = fmod(x+pi,2*pi);
		y = y - pi;
	}

	return y;
}


