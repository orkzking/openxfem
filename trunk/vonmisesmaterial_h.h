//   ********************************
//   *** CLASS VONMISESMATERIAL_H ***
//   ********************************
 
#ifndef _VONMISESMATERIAL_H_H
#define _VONMISESMATERIAL_H_H

#include "vonmisesmaterial.h"


//!Von Mises material with linear kinematic and isotropic hardening.
/*!
 This class implements a Von Mises material with linear
 kinematic and isotropic hardening.

 The hardening parameters are H', the plastic modulus,
 and beta (= 0 => kinematic, = 1 => isotropic)
*/
class VonMisesMaterial_H : public VonMisesMaterial{
public :
	   VonMisesMaterial_H (int n,Domain* d) : VonMisesMaterial(n, d) { ; } ;
      ~VonMisesMaterial_H () { ; } ;    

      char*    giveClassName (char* s) ;
      void     instanciateYourself () ;
      void     printYourself () ;
	  void     ComputeStress (FloatArray* ,Element*, GaussPoint*) ;
	  FloatMatrix* ComputeConstitutiveMatrix (GaussPoint* , Element*);

	  FloatArray* computeKsi (FloatArray*, FloatArray*) ;
	  FloatArray* computeDFDSigma (FloatArray*, FloatArray*) ;
	  double computeYieldFunctionFor (FloatArray*, FloatArray*, double) ;
	  double computeStressLevelFor(GaussPoint*) ;

	  bool  isElastic(){return false;}

} ;
#endif // _VONMISESMATERIAL_H_H
