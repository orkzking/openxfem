//   **********************
//   *** CLASS NLSOLVER ***
//   **********************


#ifndef _NLSOLVER_H_
#define _NLSOLVER_H_

#include "femcmpnn.h"
#include "flotarry.h"

class Dictionary; class Domain; class LinearSystem;
//! A class of a Nonlinear Solver 
/*!
 This class implements a NonLinear solver (subclasses are NEWTON-RAPHSON,
 BFGS, ...)
 */
class NLSolver : public FEMComponent
{  
   protected :
      Dictionary* propertyDictionary ;
		Domain* domain;
		LinearSystem* linearSystem;
		int maxIterations;
		double tolerance;
		int numberOfIterations;
		int convergenceStatus;
		int consistentDep;
		int currentIteration;

   public :
	   LinearSystem* giveLinearSystem() { return linearSystem; } ;
	   int giveNumberOfIterations() { return	numberOfIterations; } ;
	   int giveCurrentIteration() { return	currentIteration; } ;
	   int giveConvergenceStatus() { return	convergenceStatus; } ;
	   int giveConsistentDep() { return	consistentDep; } ;
	   void updateYourself();
       NLSolver (int n,Domain* d) ;                  // constructor
       ~NLSolver () ;								        // destructor
	   double  give (char) ;
	   NLSolver* ofType (char*) ;
	   NLSolver* typed () ;
	   void  giveKeyword (char*) ;
	   char* giveClassName (char*) ;
	   virtual void instanciateYourself () { } ;
	   FloatArray* Solve();

} ;

#endif // _NLSOLVER_H_
