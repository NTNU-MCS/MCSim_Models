/*************************************************************************
 *
 * FILE		: $Archive: $ 
 * DATE		: 1999.02.19
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT	: MMALLOC is a dynamic memory allocation manager system.
 *            Instead of relying of the operation system to handle
 *            allocation and de-allocation, a new set of "alloc" and "free"
 *            functions had been implemented. 
 *            Main advantages: 
 *               1) The algoritm does not fragment memory
 *               2) Easier detection of memory leaks.
 *               3) Increased speed.
 *               4) New functionality can easily be included.
 *
 *			  Doing operations on a global storage as malloc and free does
 *			  is very critcal. To make sure that only one thread at a time
 *			  has access to the crical region, WIN32 API calls are used
 *			  for thread synchronization. Thus, this library can only be
 *			  used on WIN32-platforms.
 *
 * NOTES	: The algoritm is very much inspired by Kernighan & Ritchie
 *
 * AUTHOR(s): MAD/Karl-Petter Lindegaard
 *
 ************************************************************************/ 

#ifndef MMALLOC_H
#define MMALLOC_H


void *mmalloc(unsigned int nbytes);
void mfree(void *ap);
int mmemleak();
long mbytesfree();


#endif








