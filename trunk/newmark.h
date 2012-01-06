//   *********************
//   *** CLASS NEWMARK ***
//   *********************
 

#include "timinteg.h"


class Newmark : public TimeIntegrationScheme
/*!
   This class implements a Newmark beta-gamma predictor-corrector time inte-
   gration scheme, with accelerations as primary unknowns.
   .A Newmark scheme assigns to the dofs the following 5 unknowns :
       - the displacement 'd' and its predictor 'D';
       - the velocity 'v' and its predictor 'V' ;
       - the acceletion 'a'.
    Their expressions are calculated by the dofs in class Dof.
   .A Newmark scheme assigns to the element an equivalent mass matrix as
    left-hand side, and an equivalent right-hand side as right-hand side.
    Their expressions are claculated by the element in class Element.
   .A Newmark scheme enforces the linear system's left-hand side to be recom-
    puted at every time step where the time increment deltaT is modified.
 DESCRIPTION :
   A Newmark scheme is characterized by two coefficients 'beta' and 'gamma'.
   These are stored as pointers, rather than numbers, so that their state
   (initialized vs NULL) can be checked.
 */
{
   enum { FALSE , TRUE } ;

   private :
      double*  beta ;
      double*  gamma ;

   public :
      Newmark (int,Domain*) ;
      ~Newmark ()                     {delete beta ; delete gamma ;}

      double     giveBeta () ;
      double     giveGamma () ;
      TimeStep*  GiveInitialStep ()   {return new TimeStep(0,this);}
      int        hasDiagonalLhs ()    {return (this->giveBeta() == 0.) ;}
      int        isNewmark ()         {return TRUE ;}
} ;








