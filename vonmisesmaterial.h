//   ******************************
//   *** CLASS VONMISESMATERIAL ***
//   ******************************
 
#ifndef _VONMISESMATERIAL_H_
#define _VONMISESMATERIAL_H_

#include "material.h"

//! A Von Mises material 
/*!
   This class implements a material in a finite element problem. A material
   is an attribute of a domain. It is usually also attribute of many elements.
 DESCRIPTION
   The attribute 'propertyDictionary' contains all the properties of a mate-
   rial, like its Young modulus, its mass density, its area or its thickness.
 TASK
   Returning a material property (method 'give') ;
*/
class VonMisesMaterial : public Material{
public :
	  VonMisesMaterial (int n,Domain* d) : Material(n, d) { ; } ;
     ~VonMisesMaterial () { ; } ;    

     char*    giveClassName (char* s) ;
     void     instanciateYourself () ;
     void     printYourself () ;
	  void     ComputeStress (FloatArray* ,Element*, GaussPoint*) ;
	  FloatMatrix* ComputeConstitutiveMatrix (GaussPoint* , Element*);

	  FloatArray* computeDFDSigma (FloatArray*) ;
	  double computeYieldFunctionFor (FloatArray*) ;
	  double computeStressLevelFor(GaussPoint*) ;

	  bool  isElastic(){return false;}


} ;
#endif //_VONMISESMATERIAL_H_
