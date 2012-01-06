//    file LIST.CPP
 
#include "list.h"
#include "femcmpnn.h"
#include "string.h"
#include "debug.def"
#include <stdio.h>
#include <stdlib.h>


List ::  List (int s)
   // Constructor : creates a list of size s.
{
   register int   i ;
   FEMComponent** p ;

   size = s ;
   if (size) {
      values = new FEMComponent* [size] ;
      p      = values ;
      i      = size ;
      while (i--)
	 *p++ = NULL ;}                              // initialize 'values'
   else
      values = NULL ;
}


List :: ~List ()
// Destructor.
{
  int i = size;

  if (size) {
	 while (i--)
		delete (values[i]) ;
	 delete [] values ;
  }
}


void  List :: growTo (int newSize)
   // Expands the receiver from its current size to newSize, in order to acco-
   // modate new entries.
{
   register int i ;
   FEMComponent **newValues,**p1,**p2 ;

#ifdef DEBUG
   if (newSize <= size) {
      printf ("New list size (%d) not larger than current size (%d) ! \n",
	       newSize,size) ;
      exit(0) ; }
#endif

   newValues = new FEMComponent* [newSize] ;
   p1        = values ;
   p2        = newValues ;
   for (i=0 ; i<size ; i++)
      *p2++ = *p1++ ;
   for (i=size ; i<newSize ; i++)
      *p2++ = NULL ;

   if (values)
      delete [] values ;

   values = newValues ;
   size   = newSize ;
}


int  List :: includes (int i)
   // Returns True if the receiver has a non-null i-th entry, else returns
   // False.
{
   if (i > size)
      return FALSE ;
   else
      return (values[i-1]!=NULL) ;
}


void  List :: printYourself ()
   // Prints the receiver on screen.
{
   register int i ;

   printf ("list of components of size %d\n",size) ;
   for (i=1 ; i<=size ; i++) {
      if (values[i-1] == NULL)
	 printf ("%d  Nil \n",i) ;
      else
	 printf ("%d %lX\n",i,(long int)(values[i-1])) ;}
}


void  List :: put (int i, FEMComponent* anObject)
   // Stores anObject at position i. Enlarge the receiver if too small.
{
   if (size < i)
      this -> growTo(i) ;
   values[i-1] = anObject ;
}


      









