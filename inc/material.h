//   **********************
//   *** CLASS MATERIAL ***
//   **********************


#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include "femcmpnn.h"
#include "flotarry.h"

class Dictionary ; class GaussPoint; class Element;
class FloatMatrix;


//! An abstract class for Material
/*!
This class implements a material in a finite element problem. A material
is an attribute of a domain. It is usually also attribute of many elements.
DESCRIPTION
The attribute 'propertyDictionary' contains all the properties of a mate-
rial, like its Young modulus, its mass density, its area or its thickness.
TASK
Returning a material property (method 'give') ;
*/
class Material : public FEMComponent
{
protected :
  Dictionary*    propertyDictionary ;

public :
  Material (int n,Domain* d) ;
  ~Material () ;    

  virtual void giveKeyword (char*) ;
  virtual void instanciateYourself () { ; } ;
  double   give (char) ;
  char*    giveClassName (char* s) ;
  void     printYourself () ;
  Material* ofType (char* aClass) ;
  Material* typed () ;

  // check if material is elastic or not
  // virtual method that requires derived classes implementing this method!
  virtual bool  isElastic(){return NULL;}

  //new-SC
  virtual void  ComputeStress(FloatArray* ,Element*, GaussPoint*){;} ;
  virtual FloatMatrix* ComputeConstitutiveMatrix (GaussPoint* , Element*){return NULL ;} ;
  virtual double computeStressLevelFor(GaussPoint*){ return 0; } ;
} ;
#endif // _MATERIAL_H_

