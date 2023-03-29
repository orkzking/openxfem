//   ***************************
//   *** CLASS LINEAR SYSTEM ***
//   ***************************
 

#ifndef _LINSYST_H_
#define _LINSYST_H_

#include "flotarry.h"
#include "skyline.h"


class LinearSystem
/*!
 DESCRIPTION :
   This class defines the attribute 'leftHandSide', 'rightHandSide' and
   'solutionArray'. Attribute 'currentSize' is used during the equation
   numbering process (see below).
 TASKS :
   - solving itself (implemented by the subclasses of LinearSystem) ;
   - assigning an equation number to any degree of freedom which asks for
     one (method 'giveUpdatedCurrentSize'). The system does so by returning
     its attribute 'currentSize' incremented every times by 1.
 */
{
   protected :
      Skyline*     leftHandSide ;
      FloatArray*  rightHandSide ;
      FloatArray*  solutionArray ;
      int          currentSize ;

   public :
      LinearSystem () ;
      ~LinearSystem () ;

      void          carveYourselfFor (Domain*) ;
      Skyline*      giveLhs ()                    { return leftHandSide ;}
      FloatArray*   giveRhs ()                    { return rightHandSide ;}
      FloatArray*   giveSolutionArray ()          { return solutionArray ;}
      int           giveUpdatedSize ()            { return ++currentSize ;}
      void          solveYourself () ;
      void          solveYourselfAndCheckSolution () ;
      void          updateYourself () ;
      void          updateYourselfExceptLhs () ;

	   void          setLHSTo (Skyline* aLHS);
	   void          setRHSTo (FloatArray* aFloatArray);
		void          resetCurrentSize(){currentSize = 0 ;} // NVP 2005-09-05

		void			  initializeLinearSystem();
} ;


#endif // _LINSYST_H_
