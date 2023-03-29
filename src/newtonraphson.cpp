//   file NEWTONRAPHSON.CPP

#include "newtonraphson.h"
#include "dictionr.h"

NewtonRaphson :: NewtonRaphson (int n,Domain* aDomain) : NLSolver (n, aDomain) 
   // Constructor. Creates a NewtonRaphson solver.
{	;
}


NewtonRaphson :: ~NewtonRaphson ()
   // Destructor.
{  ;
}

char* NewtonRaphson :: giveClassName (char* s)
{	return strcpy(s,"NewtonRaphson") ;
}

void  NewtonRaphson :: instanciateYourself ()
{
   double value ;

#  ifdef VERBOSE
      printf ("\ninstanciating Newton Raphson solver %d\n",number) ;
#  endif

	//maximum number of iterations
   if (this->readWhetherHas("n")) {
      value = this -> read("n") ;
      propertyDictionary -> add('n',value) ;}

   //tolerance for RHS norm
   if (this->readWhetherHas("t")) {
      value = this -> read("t") ;
      propertyDictionary -> add('t',value) ;}
}


