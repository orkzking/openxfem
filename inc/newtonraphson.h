//   ****************************
//   *** CLASS NEWTON-RAPHSON ***
//   ****************************


#ifndef NewtonRaphson_h
#define NewtonRaphson_h
#include "nlsolver.h"

class NewtonRaphson : public NLSolver
/*!
   This class implements a Newton-Raphson solver
 */
{	public :
		NewtonRaphson (int n,Domain* d) ;	// constructor
       ~NewtonRaphson () ;			// destructor
	   char* giveClassName (char*) ;
	   void instanciateYourself () ;
} ;

#endif
