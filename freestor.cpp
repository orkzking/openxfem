//   file FREESTOR.C

#include "freestor.h"
#include <stdio.h>
#include <stdlib.h>


double*  allocDouble (int n)
   // Allocates and returns an array of n doubles. Issues an error message if
   // the free store is exhausted.
{
   double* answer ;

#  ifdef DEBUG
      if (! n) {
	 printf ("cannot allocate 0 coefficients in fucntion allocDouble\n");
	 exit(0) ;}
#  endif

   answer = (double*) calloc(n,sizeof(double)) ;
   if (answer)
      return answer ;
   else
      freeStoreError() ;
   return NULL ; // SC
}


int*  allocInt (int n)
   // Allocates and returns an array of n integers. Issues an error message
   // if the free store is exhausted.
{
   int* answer ;

#  ifdef DEBUG
      if (! n) {
	 printf ("cannot allocate 0 coefficients in function 'allocInt'\n") ;
	 exit(0) ;}
#  endif

   answer = (int*) calloc(n,sizeof(int)) ;
   if (answer)
      return answer ;
   else
      freeStoreError() ;
   return NULL ; //SC
}


void  freeStoreError ()
   // This function is called whenever operator "new" is unable to allocate 
   // memory.
{
   printf ("Sorry : free store exhausted \n") ;
   exit(0) ;
}


void  freeDouble (double* a)
   // Deallocates the list of decimals 'a'.
{
#ifdef SUN_STATION
   free ((char*)a) ;        /* the Sun compiler is really lousy */
#else
   free (a) ;
#endif
}


void  freeInt (int* a)
   // Deallocates the list of integers 'a'.
{
#ifdef SUN_STATION
   free ((char*)a) ;        /* see above */
#else
   free (a) ;
#endif
}
