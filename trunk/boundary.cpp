//   file BOUNDARY.CPP
 
#include "boundary.h"

#include "timestep.h"
#include "loadtime.h"


double  BoundaryCondition :: give (char u, TimeStep* stepN)
   // Returns the value at stepN of the prescribed value of the kinematic 
   // unknown 'u'. Returns 0 if 'u' has no prescribed value.
{
   double value,factor ;
 
   if (! prescribedValueDictionary)
      this -> readConditions() ;

   if (prescribedValueDictionary -> includes(u)) {
      value  = prescribedValueDictionary -> at(u) ;
      factor = this -> giveLoadTimeFunction() -> at(stepN->giveTime()) ;
      return value*factor ;}
   else
      return 0. ;
}


void  BoundaryCondition :: readConditions ()
   // Sets up the dictionary where the receiver stores the conditions it
   // imposes.
{
   char   key [32] ; // 'key' is eventually of size 1, but some words that are
   double value ;    // read in the meantime in the data file can be larger !
   int    nCond,i,j ;

   prescribedValueDictionary = new Dictionary() ;
   i     = 0 ;
   nCond = this->readInteger("conditions",++i) ;

   for (j=1 ; j<=nCond ; j++) {
      this -> readString("conditions",++i,key) ;
      value = this->read("conditions",++i) ;
      prescribedValueDictionary -> add(key[0],value) ;}
}

