//   file NEWMARK.CPP
 
#include "newmark.h"


Newmark :: Newmark (int n, Domain* aDomain)
	 : TimeIntegrationScheme (n,aDomain)
   // Constructor.
{
   beta  = NULL ;
   gamma = NULL ;
}


double  Newmark :: giveBeta ()
   // Returns the beta coefficient of the receiver. Reads it if it does not
   // exist yet.
{
   if (beta)
      return *beta ;
   else {
      beta    = new double[1] ;
      (*beta) = this -> read("beta") ;
      return *beta ;}
}


double  Newmark :: giveGamma ()
   // Returns the gamma coefficient of the receiver. Reads it if it does not
   // exist yet.
{
   if (gamma)
      return *gamma ;
   else {
      gamma    = new double[1] ;
      (*gamma) = this -> read("gamma") ;
      return *gamma ;}
}


