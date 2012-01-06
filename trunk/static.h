//   ********************
//   *** CLASS STATIC ***
//   ********************


#include "timinteg.h"
#include "timestep.h"

class Static : public TimeIntegrationScheme
/*!
   This class implements a quasi-static time integration scheme. Multiple
   time steps in quasi-statics are useful for computing a structure subjected
   to various load cases.
   The static scheme assigns to the dofs a unique unknown :
      - the displacement 'd' .
   The static scheme assumes that the linear system's left-hand side remains
   unmodified all along the time history.
*/
{
   public :
      Static (int i,Domain* d) : TimeIntegrationScheme(i,d)  { }

      TimeStep*  GiveInitialStep ()  { return new TimeStep(1,this);}
      int        isStatic ()         { return TRUE ;}
      int        requiresNewLhsAt (TimeStep* stepN)
				     { return (stepN->giveNumber() == 1) ;}
} ;











