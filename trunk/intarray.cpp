//   file INTARRAY.CPP
 
#include "intarray.h"
#include "freestor.h"
#include <stdio.h>
#include "assert.h"


IntArray :: IntArray (int n)
   // Constructor : creates an array of size n (filled with garbage).
{
   size = n ;
   if (n)
      values = allocInt(size) ;
   else
      values = NULL ;
}


#ifdef DEBUG
int&  IntArray :: at (int i)
   // Returns the i-th coefficient of the receiver. Slow but safe.
{
   this -> checkBounds(i) ;
   return values[i-1] ;
}
#endif


#ifdef DEBUG
void  IntArray :: checkBounds (int i)
   // Checks that the receiver includes an index i.
{
   if (i<0) {
      printf ("array error on index : %d < 0 \n",i) ;
      assert(false) ; 
	}
   if (i>size) {
      printf ("array error on index : %d > %d \n",i,size) ;
      assert(false) ; 
	}
}
#endif


IntArray*  IntArray :: followedBy (IntArray* b)
   // Appends the array 'b' the receiver. Returns the receiver.
{
   register int i ;
   int          newSize ;
   int          *newValues,*p1,*p2,*p3 ;

   newSize = size + b->size ;
   if (newSize == size)
      return this ;

   newValues = allocInt(newSize) ;
   p1        = values ;
   p2        = b->values ;
   p3        = newValues ;

   i = size ;
   while (i--)
      *p3++ = *p1++ ;

   i = b->size ;
   while (i--)
      *p3++ = *p2++ ;

   if (values)
      freeInt (values) ;
   values = newValues ;
   size   = newSize ;
   return this ;
}


void  IntArray :: printYourself ()
   // Prints the receiver on screen.
{
   printf ("IntArray of size : %d\n",size) ;
   for (int i=1 ; i<=size ; ++i) {
      if (i>15) {
	 printf ("   (other components not printed)") ;
	 break ;}
      else
	 printf ("%d  ",this->at(i)) ;}
   printf ("\n") ;
}


