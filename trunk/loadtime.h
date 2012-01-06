//   ********************************
//   *** CLASS LOAD-TIME FUNCTION ***
//   ********************************
 

#ifndef loadtime_h

#include "femcmpnn.h"
#include "domain.h"


class LoadTimeFunction : public FEMComponent
/*!
   This abstract class is the superclass of the classes that implement
   functions  y = f(t) , where t is the time. These function are used for
   weighing the loads (see TJR Hughes, "The Finite Element Method", p 677).
   A load-time function is attribute of the domain. It is usually also at-
   tribute of one or more loads.
 DESCRIPTION
   The parameters which describe the function are defined by the subclasses
   of this class.
 TASK
   Returning the value 'y' at any abscissa 't'.
 */
{
   public :
      LoadTimeFunction (int i,Domain* d) : FEMComponent(i,d) {} 
      virtual ~LoadTimeFunction ()  {}

      // computations 
      virtual double     at (double)            { return 0. ;}

      // definition of a function
      LoadTimeFunction*  typed () ;
      LoadTimeFunction*  ofType (char*) ;
      char*              giveClassName (char* s) 
                                     { return strcpy(s,"LoadTimeFunction") ;}
      virtual void       instanciateYourself () {}

} ;

#define loadtime_h
#endif









