//   *****************************
//   *** CLASS ELASTICMATERIAL ***
//   *****************************
 
#ifndef _ELASTICMATERIAL_H_
#define _ELASTICMATERIAL_H_

#include "material.h"

//! Elastic material
class ElasticMaterial : public Material
/*!
   This class implements a material in a finite element problem. A material
   is an attribute of a domain. It is usually also attribute of many elements.
 DESCRIPTION
   The attribute 'propertyDictionary' contains all the properties of a mate-
   rial, like its Young modulus, its mass density, its area or its thickness.
 TASK
   Returning a material property (method 'give') ;
 */
{
   public :
	   ElasticMaterial (int n,Domain* d) : Material(n, d){}
      ~ElasticMaterial () {}   

      char*    giveClassName (char* s) ;
      void     instanciateYourself () ;
      void     printYourself () ;

		bool     isElastic(){return true;}

	   void     ComputeStress (FloatArray* ,Element*, GaussPoint*) ;
	   FloatMatrix* ComputeConstitutiveMatrix (GaussPoint* , Element*);
	   double   computeStressLevelFor(GaussPoint*) { return 0; } ;

} ;
#endif // _ELASTICMATERIAL_H_
