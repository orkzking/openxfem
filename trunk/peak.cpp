//   file PEAK.CPP
 
#include "peak.h"
#include <math.h>


double  PeakFunction :: at (double time)
   // Returns the value of the receiver at time 'time'.
{
   const double precision = 0.000001 ;

   if (! value)
      this -> getCoefficients() ;

   if (fabs(*t - time) < precision)
      return *value ;
   else
      return 0. ;
}


void  PeakFunction :: getCoefficients ()
   // Reads the date anf the time increment of the receiver in the data file.
{
   t     = new double[1] ;
   value = new double[1] ;

   *t     = this -> read("t") ;
   *value = this -> read("f(t)") ;
}


