//   ********************************************
//   *** DYNAMIC MEMORY ALLOCATION PROCEDURES ***
//   ********************************************


#ifndef _FREESTOR_H_
#define _FREESTOR_H_


/*
   This file does not define a class. Rather, it provides a few procedures
   related with dynamic memory allocation.
*/

   double*  allocDouble (int) ;
   int*     allocInt (int) ;
   void     freeStoreError () ;
   void     freeInt (int*) ;
   void     freeDouble (double*) ;

#endif // _FREESTOR_H_
