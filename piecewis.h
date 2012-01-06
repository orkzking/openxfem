//   ***************************************
//   *** CLASS PIECEWISE LINEAR FUNCTION ***
//   ***************************************
 

#include "loadtime.h"


class PiecewiseLinFunction : public LoadTimeFunction
/*!
   This class implements a piecewise linear function.
 DESCRIPTION
   The function is defined by 'numberOfPoints' points. 'dates' and 'values'
   store respectively the abscissas (t) and the values (f(t)) of the points
*/
{
   private :
      int      numberOfPoints ;
      double*  dates ;
      double*  values ;

   public :
      PiecewiseLinFunction (int i,Domain* d) : LoadTimeFunction(i,d)
			    { numberOfPoints=0 ; dates=NULL ; values=NULL ;}
      ~PiecewiseLinFunction ()             { delete dates ; delete values ;}

      double  at (double) ;
      void    getPoints () ;
      void    instanciateYourself ()       { this->getPoints() ;}
} ;

