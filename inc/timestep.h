//   ***********************
//   *** CLASS TIME STEP ***
//   ***********************


#ifndef _TIMESTEP_H_
#define _TIMESTEP_H_

#include "femcmpnn.h"
#include <stdio.h>
class TimeIntegrationScheme ;


//! A time step 
/*!
   This class implements a step in a time history.
 DESCRIPTION :
   A time step n is characterized by its date 't', and by the time increment
   'dt' between step n-1 and step n.
   The time step knows the time history 'scheme' it belongs to. Its attribute
   'newLhs' is True if the LHS of the linear system must be (re)computed,
   else it is False ; for example, in quasi-static analysis, 'newLhs; is True
   if the step's number is 1 (the first time, the LHS must be assembled), and
   False otherwise (the same LHS is used - no need for recalculating it).
   The time obtains the value of 'newLhs' from its scheme.
 TASKS :
   - returning its time 't' and its time increment 'dt'.
   - returning the boolean 'newLhs' (method 'requiresNewLhs').
*/
class TimeStep : public FEMComponent{
   enum { FALSE , TRUE , UNKNOWN } ;

   protected :
      TimeIntegrationScheme*  scheme ;
      double                  t ;
      double                  deltaT ;
      int                     newLhs ;

   public :
      TimeStep (int,TimeIntegrationScheme*) ;       // constructors
      TimeStep (TimeIntegrationScheme*,TimeStep*) ;

      char*      giveClassName (char* s)    { return strcpy(s,"TimeStep") ;}
      TimeStep*  givePreviousStep() ;
      double     giveTime ()                { return t ;}
      double     giveTimeIncrement ()       { return deltaT ;}
      int        isNotTheLastStep () ;
      int        isTheCurrentTimeStep () ;
      int        requiresNewLhs () ;
} ;


#endif //_TIMESTEP_H_
