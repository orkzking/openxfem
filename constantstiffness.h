//   ********************************
//   *** CLASS CONSTANT STIFFNESS ***
//   ********************************


#ifndef _CONSTANTSTIFFNESS_H_
#define _CONSTANTSTIFFNESS_H_

#include "nlsolver.h"

//! Constant stiffness solver
class ConstantStiffness : public NLSolver  
/*!
   This class implements a Constant stiffness solver
 */
{	public :
		ConstantStiffness (int n,Domain* d) ;	// constructor
       ~ConstantStiffness () ;			// destructor
	   char* giveClassName (char*) ;
	   void instanciateYourself () ;
} ;

#endif // _CONSTANTSTIFFNESS_H_
