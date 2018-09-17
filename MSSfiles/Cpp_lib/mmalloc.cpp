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
 * AUTHOR(s): Karl-Petter Lindegaard
 *
 ************************************************************************/ 


#define WIN32THREADSAFE    // Undefine this line for non-Windows targets

#ifdef WIN32THREADSAFE
#include <windows.h>
#define ENTERCRITICALSECTION EnterCriticalSection (&MMVAR::CriticalSection);
#define LEAVECRITICALSECTION LeaveCriticalSection (&MMVAR::CriticalSection);

#else
#define ENTERCRITICALSECTION 
#define LEAVECRITICALSECTION 
#endif


/*************************************************************************

	Define a maximum memory buffer to allocate on

*************************************************************************/

#define MAXSIZE 64*1024	/* Default: 64 kb */


/*************************************************************************

	Type definitions

*************************************************************************/

typedef double Align;

union header {
	struct {
		union header *ptr;
		unsigned int size;
	} s;
	Align x;
};

typedef union header Header;


/*************************************************************************

	Variables

*************************************************************************/

namespace MMVAR
{
	static Header *base;			/* empty list to get started */
	static Header *freep = NULL;	/* start of free list */

	static char buffer[MAXSIZE];

#ifdef WIN32THREADSAFE
	CRITICAL_SECTION CriticalSection;
#endif
}



/*********************************************************************
 *
 * FUNCTION		: mmalloc
 * DATE			: 1999.02.19
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT		: General storage allocator. Similar to the 
 *                ANSI C "malloc". 
 *				  
 * PARAMETERS	:
 * nbytes		-I	int		Number of bytes to make room for.
 *
 * RETURNS		:
 *				-O	void*	Pointer to assigned memory block.
 *
 * NOTES		: If the NULL pointer is returned, no memory left.
 *
 * AUTHOR(s)	: MAD/Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
void *mmalloc(unsigned int nbytes)
{
	Header *p, *prevp;
	unsigned int nunits;

	nunits = (nbytes+sizeof(Header)-1)/sizeof(Header) + 1;

	if ((prevp = MMVAR::freep) == NULL) /* no free list yet */
	{
		/* THIS IS FOR WINDOWS NT ------------------------------------------------------
		base = (Header*)VirtualAlloc (NULL, MAXSIZE, 
			MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE); // Allocate once and for all 
		------------------------------------------------------------------------------*/
#ifdef WIN32THREADSAFE
		InitializeCriticalSection (&MMVAR::CriticalSection);
#endif
		ENTERCRITICALSECTION

		MMVAR::base = (Header*)MMVAR::buffer;  /* Instead of NT system call */

		p = MMVAR::base+1;
		MMVAR::base->s.ptr = p;
		MMVAR::base->s.size = 0;
		MMVAR::freep = prevp = MMVAR::base;

		p->s.ptr = MMVAR::freep;
		p->s.size = MAXSIZE/sizeof(Header) - 2;

		LEAVECRITICALSECTION	
	}

	ENTERCRITICALSECTION
	for (p = prevp->s.ptr ; ; prevp = p, p = p->s.ptr)
	{
		if (p->s.size >= nunits)		/* big enough */
		{
			if (p->s.size == nunits)	/* exactly */
				prevp->s.ptr = p->s.ptr;
			else						/* allocate tail end */
			{
				p->s.size -= nunits;
				p += p->s.size;
				p->s.size = nunits;
			}
			MMVAR::freep = prevp;

			LEAVECRITICALSECTION	
			return (void *)(p+1);
		}
		if (p == MMVAR::freep)		/* Wrapped around free list */
		{
			LEAVECRITICALSECTION	
			return NULL;
		}
	}
}





/*********************************************************************
 *
 * FUNCTION		: mfree
 * DATE			: 1999.02.19
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT		: General storage de-allocator. Similar to the 
 *                ANSI C "free". 
 *				  
 * PARAMETERS	:
 * ap			-I	void*	Pointer to block to free.
 *
 * RETURNS		:
 *
 * NOTES		:
 *
 * AUTHOR(s)	: MAD/Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
void mfree(void *ap)
{
	Header *bp, *p;

	ENTERCRITICALSECTION

	bp = (Header *)ap - 1;	/* Point to block header */
	for (p = MMVAR::freep; !(bp > p && bp < p->s.ptr); p = p->s.ptr)
		if (p >= p->s.ptr && (bp > p || bp < p->s.ptr))
			break;	/* Freed block at start or end of arena */

	if (bp + bp->s.size == p->s.ptr)	/* Join to upper nbr */
	{
		bp->s.size += p->s.ptr->s.size;
		bp->s.ptr = p->s.ptr->s.ptr;
	}
	else
		bp->s.ptr = p->s.ptr;

	if (p + p->s.size == bp)			/* Join to lower nbr */
	{
		p->s.size += bp->s.size;
		p->s.ptr = bp->s.ptr;
	}
	else
		p->s.ptr = bp;
	MMVAR::freep = p;

	LEAVECRITICALSECTION	
}



/*********************************************************************
 *
 * FUNCTION		: mmemleak
 * DATE			: 1999.02.19
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT		: mmemleak() is a test whether everything inside the 
 *                allocation buffer is freed or not. Thus, it can be 
 *                only used to track memory leaks are the outer most 
 *                level.
 *				  
 * PARAMETERS	:
 *
 * RETURNS		:
 *				-O	int		0 = All has been deleted. 
 *                          1 = There is still allocated memory present
 *                              in the buffer.
 *
 * NOTES		:
 *
 * AUTHOR(s)	: MAD/Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
int mmemleak()	/* Returns  */
{
	int bytes;

	if (MMVAR::freep)
	{
		bytes = sizeof(Header)*(MMVAR::base->s.ptr->s.size + 2);
		return bytes < MAXSIZE;
	}
	else
		return 0;
}



/*********************************************************************
 *
 * FUNCTION		: mbytesfree
 * DATE			: 1999.02.19
 *
 * Copyright Karl-Petter Lindegaard 1999
 *
 * ABSTRACT		: mbytesfree() returns the number of free memory in bytes 
 *                inside the buffer. Well suited tool for tracking memory
 *                leaks. 
 *				  
 * PARAMETERS	:
 *
 * RETURNS		:
 *				-O	long		Number of bytes free.
 *
 * NOTES		:
 *
 * AUTHOR(s)	: MAD/Karl-Petter Lindegaard
 *
 * HISTORY		:
 *
 ***********************************************************************/
long mbytesfree()
{
	long i;
	Header *p;
	
	if (MMVAR::freep)	// Has mmalloc been run, at least once?
	{
		
		ENTERCRITICALSECTION
		
		i = 0;
		if (p = MMVAR::base)
		{
			p = MMVAR::base->s.ptr;
			do
			{
				i = i + sizeof(Header)*p->s.size;
				p = p->s.ptr;
			} while (p != MMVAR::base);
		}
		else
			i = MAXSIZE - 2*sizeof(Header);
		
		LEAVECRITICALSECTION	
	}
	else
		i = MAXSIZE - 2*sizeof(Header);;

	return i;
}
