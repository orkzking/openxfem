//   ***************************
//   *** CLASS PEAK FUNCTION ***
//   ***************************
 

#include "loadtime.h"


class PeakFunction : public LoadTimeFunction
/*!
   This class implements a function that is 0 everywhere, except in a single
   point.
 DESCRIPTION
   't' is the only abscissa (time) where the function is not 0. 'value' is
   the value of the function 't'. Both 't' and 'value' are pointers, rather
   than numbers, so that their state (initialized or not) can be checked.
*/
{
   private :
      double*  t ;
      double*  value ;

   public :
      PeakFunction (int i,Domain* d) : LoadTimeFunction(i,d)
					     { t=NULL ; value=NULL ;}
      ~PeakFunction ()                       { delete t ; delete value ;}

      double  at (double) ;
      void    getCoefficients () ;
      void    instanciateYourself ()         { this->getCoefficients() ;}
} ;








