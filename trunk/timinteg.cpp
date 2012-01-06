//   file TIMINTEG.CXX

#include "timinteg.h"
#include "newmark.h"
#include "static.h"
#include <stdio.h>
#include <stdlib.h>


TimeIntegrationScheme :: TimeIntegrationScheme (int n, Domain* d)
		       : FEMComponent (n,d)
   // Constructor. Creates a scheme with number n, belonging to domain d.
{
   numberOfSteps = 0 ;
   previousStep  = NULL ;
   currentStep   = NULL ;
}


TimeStep*  TimeIntegrationScheme :: giveNextStep ()
   // Returns the next time step in the time history.
{
   if (previousStep == NULL)
      currentStep = this -> GiveInitialStep() ;
   else {
      if (previousStep->giveNumber() == this->giveNumberOfSteps())
	 currentStep = NULL ;
      else
	 currentStep = new TimeStep(this,previousStep) ;}

   return currentStep ;
}


TimeStep*  TimeIntegrationScheme :: GiveInitialStep ()
   // Should be implemented by subclass.
{
   printf ("error : TIS::GiveInitialStep() is implemented by subclasses\n") ;
   exit(0) ;

   return NULL ;
}


int  TimeIntegrationScheme :: giveNumberOfSteps ()
   // Returns the number of time steps in the time history.
{
   if (! numberOfSteps)
      numberOfSteps = this -> readNumberOf("TimeStep") ;
   return numberOfSteps ;
}


int  TimeIntegrationScheme :: requiresNewLhsAt (TimeStep* stepN)
   // Returns True if the time increment deltaT of stepN differs from that
   // of the previous step, else returns False. This is the default imple-
   // mentation of the method.
{
   TimeStep* pastStep ;
   int       n,answer ;

   n = stepN->giveNumber() ;
   if (n <= 1)
      answer = TRUE ;
   else {
      pastStep = new TimeStep(n-1,this) ;
      answer = (pastStep->giveTimeIncrement() != stepN->giveTimeIncrement());
      delete pastStep ;}

   return answer ;
}


TimeIntegrationScheme*  TimeIntegrationScheme :: typed ()
   // Returns a new element, which has the same number than the receiver,
   // but is typed (Newmark or Static).
{
   TimeIntegrationScheme* scheme ;
   char                   type[32] ;

   this -> readString("class",type) ;

   if (! strcmp(type,"Static"))
      scheme = (TimeIntegrationScheme*) new Static(number,domain) ;
   else if (! strcmp(type,"Newmark"))
      scheme = (TimeIntegrationScheme*) new Newmark(number,domain) ;
   else {
      printf ("%s : unknown time integration scheme type \n",type) ;
      exit(0) ;}

   return  scheme ;
}


void  TimeIntegrationScheme :: updateYourself ()
   // Updates the receiver at end of step.
{
   delete previousStep ;
   previousStep = currentStep ;
   currentStep  = NULL ;
}
