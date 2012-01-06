//   file NLSOLVER.CPP

#include "nlsolver.h"
#include "skyline.h"
#include "flotarry.h"
#include "dictionr.h"
#include "domain.h"
#include "linsyst.h"
#include "timinteg.h"
#include "constantstiffness.h"
#include "newtonraphson.h"
#include <stdio.h>
#include <iostream>

using namespace std;

NLSolver :: NLSolver (int n,Domain* d) : FEMComponent(n,d)
// Constructor. Creates a NLSolver.
{  
  //maybe this should be read from the dat file?
  propertyDictionary = new Dictionary() ;
  this->numberOfIterations = 0;
  this->currentIteration = 0;
  this->convergenceStatus = 0;
  this->domain = d;
  this->linearSystem = new LinearSystem();
}

NLSolver :: ~NLSolver ()
// Destructor.
{  
  delete propertyDictionary; // SC-Purify-10.09.97
}

char* NLSolver :: giveClassName (char* s)
{	
  return strcpy(s,"NLSolver") ;
}

void  NLSolver :: giveKeyword (char* keyword)
{
  strcpy(keyword,"NLSolver") ;
}

NLSolver*  NLSolver :: ofType (char* aClass)
// Returns a new nonlinear solver, which has the same number than the receiver,
// but belongs to aClass (NewtonRaphson, or ConstantStiffness,...).
{
  NLSolver* newNLSolver ;	    

  if (! strcmp(aClass,"NewtonRaphson"))
	 newNLSolver = new NewtonRaphson(number,domain) ;
  else if (! strcmp(aClass,"ConstantStiffness"))
	 newNLSolver = new ConstantStiffness(number,domain) ;
  else {
	 printf ("%s : unknown nonlinear solver type \n",aClass) ;
	 exit(0) ;}

  return newNLSolver ;
}

NLSolver*  NLSolver :: typed ()
// Returns a new nonlinear solver, which has the same number than the receiver,
// but belongs to aClass (NewtonRaphson, or ModifiedNewtonRaphson,...).

{
  NLSolver* newNLSolver ;
  char     type[32] ;

  this -> readString("class",type) ;
  newNLSolver = this -> ofType(type) ;

  return newNLSolver ;
}

double  NLSolver :: give (char aProperty)
// Returns the value of the property aProperty (e.g. maxIterations, ...)
{
  double  value ;

  if (propertyDictionary -> includes(aProperty))
	 value = propertyDictionary -> at(aProperty) ;
  else {                                         // read and store the value
	 value = this -> read(aProperty) ;
	 propertyDictionary -> add(aProperty,value) ;}
  return value ;
}

FloatArray* NLSolver :: Solve ()
// Solves the equation y=f(x)=0, where f is given by the domain.
// Returns the solution x.
{
  Skyline      *jacobian;
  FloatArray   *x,*y,*dx, *dxacc;
  double       initialNorm, norm, tolerance, maxIterations;
  const double PRECISION = 1.e-9;
  int          hasConverged = 0 ;
  FILE         *logFile;

  // MP-Purify-2.10.05
   this->linearSystem->initializeLinearSystem();

  //reading tolerance and maxIterations
  tolerance = this->give('t');
  maxIterations = this->give('n');

  // starting point
  x = this->domain->GiveInitialGuess();  
  
  logFile = fopen(this->domain->giveLogFileName(),"a");

  //write info on the screen and into the log file
  printf ("Solving step: %d\n",this->domain->giveTimeIntegrationScheme()->giveCurrentStep()->giveNumber()) ;
  fprintf (logFile,"Solving step: %d\n",this->domain->giveTimeIntegrationScheme()->giveCurrentStep()->giveNumber()) ;
  cout << endl ;

  // deltax accumulated
  dxacc = new FloatArray(this->domain->giveNumberOfFreeDofs());
  int i;
  // corrections
  for(i = 1 ; i <= maxIterations ; i++) 
  {
	 this->currentIteration = i;

	 // compute right-hand side
	 
	 y = this->domain->ComputeFunctionalAt(dxacc); 
	 jacobian = this->domain->ComputeJacobian();
	 this->linearSystem->setRHSTo(y);
	 this->linearSystem->setLHSTo(jacobian);

    /*// ----------------------  DEBUG --------------------------------------
	 // Export the LHS to Matlab 
	 FloatMatrix* copy = jacobian -> AsFloatMatrix() ;
	 FILE       *mlbFile;
	 char        mlbFileName[40] = "LHS" ; // name of the file
	 strcat(mlbFileName,".m");             // extension .m ( Matlab M file)
	 mlbFile = fopen(mlbFileName, "w");

	 fprintf (mlbFile,"lhs = [ ");
	 double val ;

	 for(size_t u = 1 ; u <= copy->giveNumberOfRows() ; u++ )
	 {
		for(size_t v = 1 ; v <= copy->giveNumberOfColumns() ; v++ )
		{
		  val = copy->at(u,v) ;
		  fprintf (mlbFile,"%8.5e ",val);
		}
		fprintf(mlbFile,"\n");
	 }
	 fprintf (mlbFile," ];\n ");
	 fclose(mlbFile);
	 // -----------------------------------------------------------------*/

	 // check the convergence
	 norm = y->giveNorm();

	 if (i == 1)
		initialNorm = norm;

	 //to avoid divisions by zero
	 if (initialNorm == 0.) initialNorm = 1.;

	 //write info on the screen and into the log file
	 printf ("Iteration: %d - RHS Norm = %12.4e\n",i,(float)(norm/initialNorm));
	 fprintf (logFile,"Iteration: %d - RHS Norm = %12.4e\n",i,(float)(norm/initialNorm));

	 //checking the convergence
	 if (norm/initialNorm < tolerance || norm < PRECISION) 
	 {
		hasConverged = 1;
		linearSystem->updateYourself();
		break; 
	 }

	 //checking the divergence
	 if (norm/initialNorm > 10) 
	 {
		//write info on the screen and into the log file
		printf ("ATTENTION: STEP HAS DIVERGED!\n");
		fprintf (logFile,"ATTENTION: STEP HAS DIVERGED!\n");
		hasConverged = 0;
		linearSystem->updateYourself();
		break; 
	 }

	 // set up and solve the linear system
	 this->linearSystem->solveYourself();
	 dx = this->linearSystem->giveSolutionArray();
	 // update the solution
	 x->add(dx);

	 // update deltax accumulated
	 dxacc->add(dx);

	 // clear the memory
	 linearSystem->updateYourself();
	 delete y;	// MP-Purify-2.10.05
  }

  this->numberOfIterations = i;
  this->convergenceStatus = hasConverged;

  fclose(logFile);

  delete dxacc; // SC-Purify-10.09.97
  delete jacobian;
  delete y;
  
  return x;
}

void NLSolver :: updateYourself()
// ******************************
// Mofified by NVP for XFEM implementation
// 2005-09-05
{ 
  currentIteration = 0 ; 
  numberOfIterations = 0 ; 
  convergenceStatus = 0 ; 

  if(this->domain->isXFEMorFEM())
	 this->giveLinearSystem()->resetCurrentSize();
} 