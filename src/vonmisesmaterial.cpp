//   file VONMISESMATERIAL.CPP

#include "vonmisesmaterial.h"
#include "dictionr.h"
#include "domain.h"
#include "gausspnt.h"
#include "element.h"
#include "mathutil.h"
#include "timinteg.h"
#include <stdlib.h>
#include <math.h>
#include <iostream> 

class MathUtil; class TimeIntegrationScheme ;

char* VonMisesMaterial :: giveClassName (char* s)
{	return strcpy(s,"VonMisesMaterial") ;
}

void  VonMisesMaterial :: instanciateYourself ()
{
  double value ;

#  ifdef VERBOSE
  printf ("\ninstanciating Von Mises material %d\n",number) ;
#  endif

  if (this->readWhetherHas("E")) {
	 value = this -> read("E") ;
	 propertyDictionary -> add('E',value) ;}

  if (this->readWhetherHas("n")) {
	 value = this -> read("n") ;
	 propertyDictionary -> add('n',value) ;}

  if (this->readWhetherHas("d")) {
	 value = this -> read("d") ;
	 propertyDictionary -> add('d',value) ;}

  if (this->readWhetherHas("A")) {
	 value = this -> read("A") ;
	 propertyDictionary -> add('A',value) ;}

  if (this->readWhetherHas("t")) {
	 value = this -> read("t") ;
	 propertyDictionary -> add('t',value) ;}

  if (this->readWhetherHas("k")) {
	 value = this -> read("k") ;
	 propertyDictionary -> add('k',value) ;}

}

void  VonMisesMaterial :: printYourself ()
// Prints the receiver on screen.
{
  printf ("Von Mises Material with properties : \n") ;
  propertyDictionary -> printYourself() ;
}

void  VonMisesMaterial :: ComputeStress (FloatArray* dElem, Element* elem, GaussPoint* gp)

{
  FloatArray *sigmaTrial, *deltaEpsilon;
  FloatMatrix *Del;

  Del = elem->giveConstitutiveMatrix()->GiveCopy();

  // strain increment corresponding to dxacc = delta x accumulated
  // (since the last converged step)
  deltaEpsilon = elem->computeStrainIncrement(gp,dElem);

  // delta sigma
  sigmaTrial = Del->Times(deltaEpsilon); // SC-Purify-10.09.97
  // sigma(prev)+delta sigma
  sigmaTrial->add(gp->givePreviousStressVector());
  // computation of f(sigmaTrial)
  double fSigTr = this->computeYieldFunctionFor(sigmaTrial);

  // check if the trial stress is outside the yield surface
  if (fSigTr > 0) {
	 // if true: plastic correction
	 double deltaGamma;
	 FloatArray *dFDSigma = this->computeDFDSigma(sigmaTrial);
	 deltaGamma = fSigTr / (dFDSigma->transposedTimes(Del->Times(dFDSigma)));
	 sigmaTrial->minus((Del->Times(dFDSigma))->times(deltaGamma));
	 //tell the GaussPoint that it has reached the yield surface
	 gp->isPlastic();
	 gp->setDeltaGamma(deltaGamma);
	 delete dFDSigma;
  }

  gp->letStressVectorBe(sigmaTrial); //gp stores sigmaTrial (corrected), not a copy!

  // deletion of temporary arrays
  delete Del;
  delete deltaEpsilon;

}

FloatMatrix*  VonMisesMaterial :: ComputeConstitutiveMatrix
(GaussPoint* gp, Element* elem)
{
  FloatMatrix *Del, *Dep;

  Del = elem->giveConstitutiveMatrix()->GiveCopy(); // SC-Purify-10.09.97

  //for debugging 

  /*
  if (elem->giveDomain()->giveTimeIntegrationScheme()->giveCurrentStep()->giveTime() == 5) {
  if (elem->giveDomain()->giveNLSolver()->giveCurrentIteration() == 2) {
  if (elem->giveNumber() == 90) {
  if (gp->giveNumber() == 1) {
  int dum = 0;
  }
  }
  }
  }
  */

  //for first iteration, always return Del
  if (elem->giveDomain()->giveNLSolver()->giveCurrentIteration() == 1) {
	 return Del;
  }

  //for Full Newton-Raphson, return Dep if plastic, Del else
  if (gp->givePlasticCode() == 1) {
	 //GaussPoint is plastic
	 FloatArray *sigma, *DFDSigma;
	 sigma = gp->giveStressVector();
	 DFDSigma = this->computeDFDSigma(sigma);
	 int aSize = sigma->giveSize();
	 //test
	 NLSolver* nlSolver;
	 nlSolver    = this->domain->giveNLSolver();
	 //if consistent, compute DelStar, else Del if enough
	 if (nlSolver->give('c') == 1) {
		FloatArray *kronDelta, *a1;
		double deltaGamma, j2, A1, A2, A3, e, nu, A4;
		//deltaGamma computation
		deltaGamma =  gp->giveDeltaGamma();
		//factors
		j2 = sigma->computeInvariantJ2();
		double sqrtJ2 = sqrt(j2);
		A1 = deltaGamma*(-.25/(j2*sqrtJ2));
		A2 = deltaGamma*(0.5/sqrtJ2);
		e = elem -> giveMaterial() -> give('E') ;
		nu = elem -> giveMaterial() -> give('n') ;
		A3 = - nu / e;
		A4 = A3 - (A2 / 3.0e0);
		//a1 computation (dj2/dsigma)
		a1 = sigma->computeDJ2DSigma();
		//kronDelta definition
		kronDelta = new FloatArray(aSize);
		kronDelta->at(1) = 1;
		kronDelta->at(2) = 1;
		kronDelta->at(3) = 0;
		kronDelta->at(4) = 1;
		//T
		FloatMatrix* T;
		T = new FloatMatrix(aSize,aSize);
		T->at(1,1) = 1;
		T->at(2,2) = 1;
		T->at(3,3) = 2;
		T->at(4,4) = 1;
		//Ginv
		FloatMatrix* Ginv;
		double G;
		G = e / (2.0e0*(1.0e0+nu));
		Ginv = T->Times(1.0e0/(2.0e0*G));
		//Ninv
		FloatMatrix *N, *Ninv;
		N = Ginv->Plus(T->Times(A2));
		Ninv = new FloatMatrix(aSize,aSize);
		for (int i = 1; i<=aSize; i++) {
		  Ninv->at(i,i) = 1 / N->at(i,i);
		}
		//MathUtil
		MathUtil* aMathUtil = new MathUtil();
		//Kinv
		FloatMatrix* Kinv;
		Kinv = aMathUtil->giveShermanMorrisonInverse(Ninv, kronDelta, A4);
		//DelStar
		Del = aMathUtil->giveShermanMorrisonInverse(Kinv, a1, A1);
		//delete
		delete kronDelta;
		delete a1;
		delete T;
		delete Ginv;
		delete N;
		delete Ninv;
		delete aMathUtil;
	 }
	 //Finally, Dep
	 FloatArray *DelTimesDFDSigma;
	 FloatMatrix *aContribution;
	 double aFactor;
	 DelTimesDFDSigma = Del->Times(DFDSigma);
	 aFactor = 1. / DFDSigma->transposedTimes(DelTimesDFDSigma);
	 aContribution = DelTimesDFDSigma->timesTransposed(DelTimesDFDSigma);
	 Dep = Del->Minus(aContribution->times(aFactor));
	 //delete
	 delete Del;
	 delete DFDSigma;
	 //delete sigma;
	 delete DelTimesDFDSigma;
	 delete aContribution;
	 //return consistent (or not) tangent matrix Dep
	 return Dep;
  } else {
	 //GaussPoint is Elastic
	 return Del;
  } 
}

FloatArray* VonMisesMaterial :: computeDFDSigma (FloatArray* r) 
{
  FloatArray *answer = new FloatArray(r->giveSize());

  double a11 = (*r)[0];
  double a22 = (*r)[1];
  double a12 = (*r)[2];
  double a33 = (*r)[3];

  double j2 = r->computeInvariantJ2();
  if (j2 == 0) return answer;

  double factor = 1.0e0 / (6.0e0*sqrt(j2));

  (*answer)[0] = factor * ((2.0e0*a11) - a22 - a33);
  (*answer)[1] = factor * ((2.0e0*a22) - a11 - a33);
  (*answer)[2] = factor * (6.0e0 * a12);
  (*answer)[3] = factor * ((2.0e0*a33) - a22 - a11);

  return answer;
}

double VonMisesMaterial :: computeYieldFunctionFor (FloatArray* sigma) 
{
  double answer;
  double k = this->give('k');
  double j2 = sigma->computeInvariantJ2();
  answer = sqrt(j2) - k;

  return answer;
}

double VonMisesMaterial :: computeStressLevelFor (GaussPoint* gp)
{
  double answer;
  double k = this->give('k');
  FloatArray* sigma = gp->giveStressVector();
  double j2 = sigma->computeInvariantJ2();
  answer = sqrt(j2) / k;
  return answer;
}


