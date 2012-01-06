//   file TIMESTEP.CXX

#include "timestep.h"
#include "domain.h"
#include "timinteg.h"
#include <stdio.h>


TimeStep :: TimeStep (int n, TimeIntegrationScheme* s)
	  : FEMComponent(n,s->giveDomain())
   // Constructor. Creates a new time step, with number n, and belonging to
   // the time history of s. Used for the initial step (0 or 1).
{
   scheme = s ;
   newLhs = UNKNOWN ;
   if (n >= 1) {
      deltaT = read("dt") ;
      t      = deltaT ;}
   else {                            // case of a starting step (0)
      deltaT = 0. ;
      t      = 0. ;}
}


TimeStep :: TimeStep (TimeIntegrationScheme* s, TimeStep* stepN)
	  : FEMComponent(stepN->giveNumber()+1,s->giveDomain())
   // Standard constructor. Creates a new step n+1, following stepN, and
   // belonging to the time history of s.
{
   scheme = s ;
   newLhs = UNKNOWN ;
   deltaT = this -> read("dt") ;
   t      = stepN->giveTime() + deltaT ;
}


TimeStep*  TimeStep :: givePreviousStep ()
   // Not accepted in-line.
{
   return scheme->givePreviousStep() ;
}


int  TimeStep :: isNotTheLastStep ()
   // Returns True if the time history contains steps after the receiver,
   // else returns False.
{
   return  (number != scheme->giveNumberOfSteps()) ;
}

			  
int  TimeStep :: isTheCurrentTimeStep ()
   // Not accepted in-line.
{
   return this==scheme->giveCurrentStep() ;
}


int  TimeStep :: requiresNewLhs ()
   // Returns True if the linear system's LHS must be recomputed, else
   // returns False.
{
   if (newLhs == UNKNOWN)
      newLhs = domain->giveTimeIntegrationScheme()->requiresNewLhsAt(this) ;
   return newLhs ;
}



