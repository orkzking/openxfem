//   file LOADTIME.CPP
 
#include "loadtime.h"
#include "constant.h"
#include "peak.h"
#include "piecewis.h"
#include "domain.h"
#include <stdlib.h>
#include <stdio.h>


LoadTimeFunction*  LoadTimeFunction :: ofType (char* aClass)
   // Returns a new function, which has the same number than the receiver,
   // but which also is typed (ConstantFunction, PeakFunction,...).
{
   LoadTimeFunction* newLTF ;

   if (! strncmp(aClass,"ConstantFunction",5))
      newLTF = new ConstantFunction(number,domain) ;
   else if (! strncmp(aClass,"PeakFunction",5))
      newLTF = new PeakFunction(number,domain) ;
   else if (! strncmp(aClass,"PiecewiseLinFunction",5))
      newLTF = new PiecewiseLinFunction(number,domain) ;
   else {
      printf("%s : unknown type of load-time function \n",aClass) ;
      exit(0) ;}

   return newLTF ;
}


LoadTimeFunction*  LoadTimeFunction :: typed ()
   // Returns a new load, which has the same number than the receiver,
   // but is typed.
{
   LoadTimeFunction*  newLTF ;
   char               type[32] ;

   this -> readString("class",type) ;
   newLTF = this -> ofType(type) ;

   return newLTF ;
}











