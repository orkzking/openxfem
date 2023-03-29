//   file CONSTANTSTIFFNESS.CPP

#include "constantstiffness.h"
#include "dictionr.h"

ConstantStiffness :: ConstantStiffness (int n,Domain* aDomain) : NLSolver (n, aDomain) 
   // Constructor. Creates a ConstantStiffness solver.
{	;
}


ConstantStiffness :: ~ConstantStiffness ()
   // Destructor.
{  ;
}

char* ConstantStiffness :: giveClassName (char* s)
{	return strcpy(s,"ConstantStiffness") ;
}

void  ConstantStiffness :: instanciateYourself ()
{
   double value ;

#  ifdef VERBOSE
      printf ("\ninstanciating Constant Stiffness solver %d\n",number) ;
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


