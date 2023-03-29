//   file INITIAL.CPP
 
#include "initial.h"
#include <stdlib.h>


double  InitialCondition :: give (char u)
   // Returns the value of the prescribed value of the kinematic unknown 'u'.
{
   if (this -> hasConditionOn(u))
      return  initialValueDictionary -> at(u) ;
   else
      return  0. ;
}


int  InitialCondition :: hasConditionOn (char u) 
   // Returns True if the receiver submits the unknown 'u' to an initial
   // condition, else returns False.
{
   if (initialValueDictionary == NULL)
      this -> readConditions() ;
   return  (initialValueDictionary->includes(u)) ;
}


void  InitialCondition :: printYourself ()
   // Prints the receiver on screen.
{
   printf ("Initial condition %d\ninitial values :\n",number) ;
   if (initialValueDictionary)
      initialValueDictionary -> printYourself() ;
   else
      printf ("   NULL\n") ;
}


void  InitialCondition :: readConditions ()
   // Sets up the dictionary where the receiver stores the conditions it
   // imposes.
{
   char   key [1] ;
   double value ;
   int    nCond,i,j ;

   initialValueDictionary = new Dictionary() ;
   i     = 0 ;
   nCond = this->readInteger("conditions",++i) ;

   for (j=1 ; j<=nCond ; j++) {
      this -> readString("conditions",++i,key) ;
      value = this->read("conditions",++i) ;
      initialValueDictionary -> add(key[0],value) ;}
}


