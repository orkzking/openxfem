
//   *************************************
//   *** CLASS TIME INTEGRATION SCHEME ***
//   *************************************


#ifndef timinteg_h

#include "femcmpnn.h"
#include "timestep.h"
#include "string.h"


class TimeIntegrationScheme : public FEMComponent
/*!
   This abstract class is the superclass of the classes that implement time
   integration algorithms (Newmark, quasi-static,...). A time integration is
   an attribute of the domain.
 DESCRIPTION :
   'numberOfSteps' is the number of time steps the problem consists in.
   'currentStep' is the time step currently solved. 'previousStep' is the
   former time step.
 TASKS :
   - managing the time history of the problem :
       .returning (including updating) the current time step and the previous
	one ;
       .storing and returning the coefficients of the finite difference
	formula ;
   - saying if elemental left-hand sides
       .are diagonal ;
       .need being recalculated and reassembled at any time step.
*/
{
   protected :
      int            numberOfSteps ;
      TimeStep*      currentStep ;
      TimeStep*      previousStep ;

      enum { FALSE , TRUE } ;

   public :
      TimeIntegrationScheme (int,Domain*) ;
      virtual ~TimeIntegrationScheme ()
			       { delete currentStep ; delete previousStep ;}

      // definition of the type (Newmark, Static, etc)
      virtual int             isNewmark ()         { return FALSE ;}
      virtual int             isStatic ()          { return FALSE ;}
      TimeIntegrationScheme*  typed () ;

      // coefficients of the finite difference formula
      virtual double          giveAlpha ()         { return 0. ;}
      virtual double          giveBeta ()          { return 0. ;}
      virtual double          giveGamma ()         { return 0. ;}

      // management of the time steps
      TimeStep*               giveCurrentStep ()   { return currentStep ;}
      virtual TimeStep*       GiveInitialStep () ;
      TimeStep*               giveNextStep () ;
      int                     giveNumberOfSteps () ;
      TimeStep*               givePreviousStep ()  { return previousStep ;}

      // effect on elemental left-hand side (diagonal, invariant, etc)
      virtual int             hasDiagonalLhs ()    { return FALSE ;}
      virtual int             requiresNewLhsAt (TimeStep*) ;

      // miscellaneous
      char*                   giveClassName (char* s)
			       { return strcpy(s,"TimeIntegrationScheme") ;}
      void                    updateYourself () ;
} ;

#define timinteg_h
#endif








