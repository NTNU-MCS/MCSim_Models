/* 
    libDataStructures.c

	Library file for a set of data structure fundamental functions related to stacks, queues, 
	lists, etc., in C.

	References:
	[1]	E. Horowitz, S. Sahni, S. Anderson-Freed, Fundamentals of Data Structures in C, 
		Computer Science Press, New York, 1993, ISBN 0-7167-8250-2.

	Copyright: 		Roger Skjetne, NTNU
    Date created: 	2011-04-03 Roger Skjetne
    Revised:      	 

*/


#include "libDataStructures.h"
#include "libGenericFunctions.h"


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* #####################################  DATA STRUCTURE FUNCTIONS	################################### */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void AddElementToIntStack(int stack[], int element) 
/* Adds new element at end of stack. If element already exists it does nothing. Care must be taken to avoid filling up the stack. */
{
	int i, found=1;
	for(i=0;i<stack[0];i++) {
		if(stack[i+1] == element) {
			found = -1;
			break;
		}
	}
	if(found>0) {
		stack[++i] = element;
		stack[0] = i;
	}
}


void ReplaceEntireIntStack(int target_stack[], int source_stack[]) 
/* Replaces target_stack[] with the source_stack[]. */
{
	int i;
	for(i=0;i<=source_stack[0];i++) {
		target_stack[i] = source_stack[i];
	}
}


void RemoveElementFromIntStack(int stack[], int element) 
/* Removes an element from the stack. If element does not exist, it does nothing. */
{
	int i, found=-1;
	for(i=1;i<=stack[0];i++) {
		if(found>0) {
			stack[found++]	= stack[i];
		}
		if((found<0) && (stack[i] == element)) {
			found = i;
		}
	}
	if(found>0)
		stack[0] = stack[0]-1;
}


int	AddToStack(SimStruct *S, s_ptr *top, sElement item)
/* Function that adds an element of type sElement to top of stack. It returns 1 if ok, and 0 if not. */
{
	char msg[128];

	s_ptr temp = (s_ptr)malloc(sizeof(sStack));
	if(!temp) {
		sprintf_s(msg,128,"Memory request failed for \'s_ptr\' struct.");
		SendMsg(S,ERROR_LOG,msg,"AddToStack");
		return -1; 
	}
	else {
		temp->item	= item;
		temp->next	= *top;
		*top		= temp;
		return 1;
	}
}


int	AddStackToStack(SimStruct *S, s_ptr *top1, s_ptr *top2, int MaxElements)
/* Function that adds a stack pointed to by top2 to top of stack pointed to by top1. It returns 1 if ok, and -1 if not. */
{
	int			status = 1, count=0;
	sElement	data;
	char		msg[128];

	while(*top2) {
		if((++count)>MaxElements) {
			sprintf_s(msg,128,"While adding stack items, the while-loop has reached maximum number of iterations. Breaking while-loop.");
			SendMsg(S,DEBUG_LOG,msg,"AddStackToStack");
			status = -1;
			break;
		}

		data = GetFromStack(S, top2);
		status = AddToStack(S, top1, data);
		if(status<1) break;
	}
	return status;
}


sElement GetFromStack(SimStruct *S, s_ptr *top)
/* Retrieves top element from stack. */
{
	char msg[128];
	s_ptr temp = *top;
	sElement item;

	if(!temp) {
		sprintf_s(msg,128,"Pointer top points to NULL. No element to retrieve.");
		SendMsg(S,DEBUG_LOG,msg,"GetFromStack");
		item.key = NULL;
	}
	else {
		item = temp->item;
		*top = temp->next;
		free(temp);
	}
	return item; 
}


int	DeleteMatchingElementInStack(SimStruct *S, s_ptr *top, sElement item, int MaxElements)
/* Finds the element in the stack that matches item and deletes it. The function returns 1: if the element is found and deleted,
   0: if the element is not found, and -1: if an error occurs.
*/
{
	char msg[256];
	int count=0, status = 0;
	s_ptr temp = *top, prev=NULL;

	if(CompareElements(temp->item, item)) { // Top element matches
		*top = temp->next;
		free(temp);
		status = 1;
	}
	else {
		prev = temp;
		temp = temp->next;
		while(temp) {
			if(CompareElements(temp->item, item)) {
				prev->next = temp->next;
				free(temp);
				status = 1;
				break;
			}
			else {
				prev = temp;
				temp = temp->next;
			}
			if((++count)>MaxElements) {
				sprintf_s(msg,256,"While looking for the element to delete, the while-loop has reached maximum number of iterations. Breaking while-loop.");
				SendMsg(S,DEBUG_LOG,msg,"DeleteMatchingElementInStack");
				status = -1;
				break;
			}
		}
	}

	return status;
}


void PurgeStackFromMemory(SimStruct *S, s_ptr *top, int MaxElements)
/* Purges all elements of a stack from memory. */
{
	char msg[256];
	int count=0;
	s_ptr temp = *top, next=NULL;

	while(temp) {
		next = temp->next;
		free(temp);
		temp = next;
		if((++count)>MaxElements) {
			sprintf_s(msg,256,"While trying to purge stack elements from memory, the while-loop has reached maximum number of iterations. Breaking while-loop.");
			SendMsg(S,DEBUG_LOG,msg,"PurgeStackFromMemory");
			break;
		}
	}
	*top = NULL;
}


int AddToListAtFront(SimStruct *S, l_ptr *first, sElement item)
/* Function that adds a new element at the front of the list. */
{
	char msg[128];

	l_ptr temp = (l_ptr)malloc(sizeof(sList));
	if(!temp) {
		sprintf_s(msg,256,"Memory request failed for \'sList\' struct.");
		SendMsg(S,ERROR_LOG,msg,"AddToListAtFront");
		return -1; 
	}

	temp->item		= item;
	temp->prev		= NULL;
	temp->next		= *first;
	if((*first)) {
		(*first)->prev	= temp;
	}
	*first			= temp;

	return 1;
}


int AddToListAtEnd(SimStruct *S, l_ptr *last, sElement item)
/* Function that adds a new element at the end of the list. */
{
	char msg[128];

	l_ptr temp = (l_ptr)malloc(sizeof(sList));
	if(!temp) {
		sprintf_s(msg,256,"Memory request failed for \'sList\' struct.");
		SendMsg(S,ERROR_LOG,msg,"AddToListAtEnd");
		return -1; 
	}

	temp->item		= item;
	temp->next		= NULL;
	temp->prev		= *last;
	if((*last)) {
		(*last)->next	= temp;
	}
	*last			= temp;

	return 1;
}


sElement RemoveElementFromList(SimStruct *S, l_ptr *first, l_ptr *last, l_ptr *ptr)
/* Retrieves the list element ptr from the list. */
{
	char msg[128];
	l_ptr temp = *ptr;
	sElement item;

	if(!temp) {
		sprintf_s(msg,128,"Pointer ptr points to NULL. No element to retrieve from list.");
		SendMsg(S,DEBUG_LOG,msg,"RemoveElementFromList");
		item.key = -1;
		return item;
	}

	item				= temp->item;

	if(temp->next)
		temp->next->prev	= temp->prev;
	else
		(*last) = temp->prev;

	if(temp->prev)
		temp->prev->next	= temp->next;
	else
		(*first) = temp->next;

	free(temp);

	return item; 
}


sElement GetFromListAtFront(SimStruct *S, l_ptr *first)
/* Retrieves first element in list. */
{
	sElement item;
	l_ptr dummy;

	item = RemoveElementFromList(S, first, &dummy, first);
	return item; 
}


sElement GetFromListAtEnd(SimStruct *S, l_ptr *last)
/* Retrieves last element in list. */
{
	sElement item;
	l_ptr dummy;

	item = RemoveElementFromList(S, &dummy, last, last);
	return item; 
}


int DeleteMatchingElementInList(SimStruct *S, l_ptr *first, l_ptr *last, sElement item, int MaxElements)
/* Finds the element in the list that matches item and deletes it. The function returns: 
		 1: if the element is found and deleted,
		 0: if the element is not found, and 
		-1: if an error occurs.
*/
{
	char msg[256];
	int count=0, status = 0;
	l_ptr temp = *first;
	sElement data;

	while(temp) {
		if((++count)>MaxElements) {
			sprintf_s(msg,256,"While looking for the element to delete, the while-loop has reached maximum number of iterations. Breaking while-loop.");
			SendMsg(S,DEBUG_LOG,msg,"DeleteMatchingElementInList");
			status = -1;
			break;
		}
		if(CompareElements(temp->item, item)) {
			data = RemoveElementFromList(S, first, last, &temp);
			status = 1;
			break;
		}
		else {
			temp = temp->next;
		}
	}
	return status; 
}


void PurgeListFromMemory(SimStruct *S, l_ptr *first, int MaxElements)
/* Purges all elements of a list from memory. */
{
	char msg[256];
	int count=0;
	l_ptr temp = *first, next=NULL;

	while(temp) {
		next = temp->next;
		free(temp);
		temp = next;
		if((++count)>MaxElements) {
			sprintf_s(msg,256,"While trying to purge list elements from memory, the while-loop has reached maximum number of iterations. Breaking while-loop.");
			SendMsg(S,DEBUG_LOG,msg,"PurgeListFromMemory");
			break;
		}
	}
	*first = NULL;
}


int	CompareElements(sElement a, sElement b)
/* Function that compares the fields of two elements of type sElement and returns 1 if these match. Otherwise it returns 0. */
{
	if(a.key == b.key)
		return 1;
	else
		return 0;
}




/*################################## END #############################################*/
