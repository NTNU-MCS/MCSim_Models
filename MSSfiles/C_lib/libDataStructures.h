/* 
    libDataStructures.h

	Header file for a set of data structure fundamental functions related to stacks, queues, 
	lists, etc., in C. 
		* Stacks are represented both as:
			- Simple integer stacks
			- Single-linked stack elements, dynamically created and freed in memory.
		* Lists are represented by doubly-linked list elements, dynamically created and freed in memory.

	References:
	[1]	E. Horowitz, S. Sahni, S. Anderson-Freed, Fundamentals of Data Structures in C, 
		Computer Science Press, New York, 1993, ISBN 0-7167-8250-2.

	Copyright: 		Roger Skjetne, NTNU
    Date created: 	2011-04-03 Roger Skjetne
    Revised:      	 

*/


#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

#include "simstruc.h"


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* #####################################  DATA STRUCTURE DEFINITIONS	################################### */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

typedef struct {	
/* Defining an element as part of a stack, queue, list, etc., in the csDataStructures library.
*/
	int		key;	// integer element
					// other fields???
} sElement;


struct sStack_element {
/* Structure that defines a stack element, containing items of type sElement. 
*/
	sElement				 item;	// Stack item of type sElement
	struct sStack_element	*next;	// Pointer to next element.
};
typedef struct sStack_element sStack;
typedef sStack *s_ptr;


struct sList_element {
/* Structure that defines a doubly-linked list element, containing items of type sElement. 
*/
	sElement				 item;	// Stack item of type sElement
	struct sList_element	*prev;	// Pointer to previous element.
	struct sList_element	*next;	// Pointer to next element.
};
typedef struct sList_element sList;
typedef sList *l_ptr;




/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* #####################################  DATA STRUCTURE DEFINITIONS	################################### */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* csDataStructureLib library functions */
void		AddElementToIntStack(int stack[], int element);
void		ReplaceEntireIntStack(int target_stack[], int source_stack[]);
void		RemoveElementFromIntStack(int stack[], int element);

int			AddToStack(SimStruct *S, s_ptr *top, sElement item);
int			AddStackToStack(SimStruct *S, s_ptr *top1, s_ptr *top2, int MaxElements);
sElement	GetFromStack(SimStruct *S, s_ptr *top);
int			DeleteMatchingElementInStack(SimStruct *S, s_ptr *top, sElement item, int MaxElements);
void		PurgeStackFromMemory(SimStruct *S, s_ptr *top, int MaxElements);

int			AddToListAtFront(SimStruct *S, l_ptr *first, sElement item);
int			AddToListAtEnd(SimStruct *S, l_ptr *last, sElement item);
sElement	RemoveElementFromList(SimStruct *S, l_ptr *first, l_ptr *last, l_ptr *ptr);
sElement	GetFromListAtFront(SimStruct *S, l_ptr *first);
sElement	GetFromListAtEnd(SimStruct *S, l_ptr *last);
int			DeleteMatchingElementInList(SimStruct *S, l_ptr *first, l_ptr *last, sElement item, int MaxElements);
void		PurgeListFromMemory(SimStruct *S, l_ptr *first, int MaxElements);
int			CompareElements(sElement a, sElement b);



#endif
