//   file VONMISESMATERIAL_H.CXX

#include "vonmisesmaterial_h.h"
#include "dictionr.h"
#include "domain.h"
#include "gausspnt.h"
#include "element.h"
#include "mathutil.h"
#include "timinteg.h"
#include <stdlib.h>
#include <math.h>


class MathUtil; class TimeIntegrationScheme ;

char* VonMisesMaterial_H :: giveClassName (char* s)
{	return strcpy(s,"VonMisesMaterial_H") ;
}

void  VonMisesMaterial_H :: instanciateYourself ()
{
  double value ;

#  ifdef VERBOSE
  printf ("\ninstanciating Von Mises material with Hardening %d\n",number) ;
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

  if (this->readWhetherHas("H")) {
	 value = this -> read("H") ;
	 propertyDictionary -> add('H',value) ;}

  if (this->readWhetherHas("b")) {
	 value = this -> read("b") ;
	 propertyDictionary -> add('b',value) ;}

}

void  VonMisesMaterial_H :: printYourself ()
// Prints the receiver on screen.
{
  printf ("Von Mises Material with Hardening with properties : \n") ;
  propertyDictionary -> printYourself() ;
}

void  VonMisesMaterial_H :: ComputeStress (FloatArray* dElem, Element* elem, GaussPoint* gp)

{
  FloatArray *sigmaTrial, *deltaEpsilon, *alphaPrev ; 
  FloatMatrix *Del;
  double radius;

  double H = this->give('H'); //plastic modulus
  double b = this->give('b'); //hardening parameter

  Del = elem->giveConstitutiveMatrix()->GiveCopy();


  // for debugging
  if (elem->giveDomain()->giveTimeIntegrationScheme()->giveCurrentStep()->giveTime() == 5) {
	 if (elem->giveDomain()->giveNLSolver()->giveCurrentIteration() == 3) {
		if (elem->giveNumber() == 90) {
		  if (gp->giveNumber() == 1) {
			 int dum = 0;
		  }
		}
	 }
  }



  // strain increment corresponding to dxacc = delta x accumulated
  // (since the last converged step)
  deltaEpsilon = elem->computeStrainIncrement(gp,dElem);

  // delta sigma
  sigmaTrial = Del->Times(deltaEpsilon); // SC-Purify-10.09.97
  // sigma(prev)+delta sigma
  sigmaTrial->add(gp->givePreviousStressVector());
  //hardening terms
  alphaPrev = gp->givePreviousBackStressVector()->GiveCopy();
  if (elem->giveDomain()->giveTimeIntegrationScheme()->giveCurrentStep()->giveNumber() == 1) {
	 radius = (this->give('k'))*sqrt(2.);
  } else {
	 radius = gp->givePreviousYieldRadius();
  }
  // computation of f(sigmaTrial)
  double fSigTr = this->computeYieldFunctionFor(sigmaTrial, alphaPrev, radius);

  // check if the trial stress is outside the yield surface
  if (fSigTr > 0) { 
	 // if true: plastic correction
	 double deltaGamma;
	 FloatArray *dFDSigma = this->computeDFDSigma(sigmaTrial, alphaPrev);
	 //consistency...
	 deltaGamma = fSigTr / ((dFDSigma->transposedTimes(Del->Times(dFDSigma))) + (2./3.)*H);
	 //correction on the stress state
	 sigmaTrial->minus((Del->Times(dFDSigma))->times(deltaGamma));
	 //correction on the hardening parameters
	 FloatArray *ksi = this->computeKsi(sigmaTrial, alphaPrev);
	 ksi->times(1. / ksi->computeTensorialNorm());
	 alphaPrev->minus((ksi->Times(deltaGamma*(2./3.)*H*(1.-b)))->negated());
	 radius += (deltaGamma*(2./3.)*H*b);
	 //tell the GaussPoint that it has reached the yield surface
	 gp->isPlastic();
	 gp->setDeltaGamma(deltaGamma);
	 delete dFDSigma;
	 delete ksi;
  }

  gp->letStressVectorBe(sigmaTrial); //gp stores sigmaTrial (corrected), not a copy!
  gp->letBackStressVectorBe(alphaPrev); //gp stores alphaPrev (corrected), not a copy!
  gp->setYieldRadius(radius);

  // deletion of temporary arrays
  delete Del;
  delete deltaEpsilon;
}

FloatMatrix*  VonMisesMaterial_H :: ComputeConstitutiveMatrix
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
	 FloatArray *sigma, *alpha, *DFDSigma;
	 sigma = gp->giveStressVector();
	 alpha = gp->givePreviousBackStressVector();
	 DFDSigma = this->computeDFDSigma(sigma, alpha);
	 FloatArray *DelTimesDFDSigma;
	 FloatMatrix *aContribution;
	 double aFactor;
	 double H = this->give('H'); //plastic modulus
	 DelTimesDFDSigma = Del->Times(DFDSigma);
	 aFactor = 1. / (DFDSigma->transposedTimes(DelTimesDFDSigma) + (2./3.)*H);
	 aContribution = DelTimesDFDSigma->timesTransposed(DelTimesDFDSigma);
	 Dep = Del->Minus(aContribution->times(aFactor));
	 //delete
	 delete Del;
	 delete DFDSigma;
	 //delete sigma;
	 //delete alpha;
	 delete DelTimesDFDSigma;
	 delete aContribution;
	 //return tangent (not yet consistent) matrix Dep
	 return Dep;
  } else {
	 //GaussPoint is Elastic
	 return Del;
  } 
}

FloatArray* VonMisesMaterial_H :: computeKsi (FloatArray* sigma, FloatArray* alpha) 
{
  FloatArray* ksi;
  ksi = sigma->computeDeviatoricPart();
  ksi->minus(alpha);
  return ksi;
}

FloatArray* VonMisesMaterial_H :: computeDFDSigma (FloatArray* sigma, FloatArray* alpha) 
{	
  FloatArray* ksi;
  ksi = this->computeKsi(sigma, alpha);
  ksi->times(1. / ksi->computeTensorialNorm());
  double ksi12 = ksi->at(3);
  ksi->at(3) = 2.0*ksi12; //twice the shear term in the derivative
  return ksi;
}

double VonMisesMaterial_H :: computeYieldFunctionFor (FloatArray* sigma, FloatArray *alpha, double r) 
{
  FloatArray* ksi;
  ksi = this->computeKsi(sigma, alpha);
  double answer = ksi->computeTensorialNorm() - r;
  delete ksi;
  return answer;
}

double VonMisesMaterial_H :: computeStressLevelFor (GaussPoint* gp) 
{
  FloatArray* sigma = gp->giveStressVector();
  FloatArray* alpha = gp->giveBackStressVector();
  FloatArray* ksi = this->computeKsi(sigma, alpha);
  double radius = gp->giveYieldRadius();
  double answer = ksi->computeTensorialNorm() / radius;
  delete ksi;
  return answer;
}


