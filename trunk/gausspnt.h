//   *************************
//   *** CLASS GAUSS POINT ***
//   *************************
 

#ifndef _GAUSSPNT_H_
#define _GAUSSPNT_H_

#include "stressarray.h"
#include "flotarry.h"
#include "geometry_base.h"
#include <stdio.h>
class Element ;


//! a Gauss point
/*!
   This class implements a point for a gaussian quadrature. A Gauss point is
   usually attribute of an element.
 DESCRIPTION
   A Gauss point is identified by its 'number' and by 'element' - the element
   it belongs to. The array 'coordinates' defines its position in an axis sys-
   tem defined by its element. 'weight' is the pondaration of the point in the
   numerical integration ; 'weight' is naturally positive and smaller than 2.
   'strainVector' and 'stressVector' store the current values of strains and
   stresses ; the size of these arrays and the nature of their values is de-
   termined by the element.
 TASKS
   - returning its coordinates, its weight and the strains and stresses ;
   - reading/writing its attribute in a file.
 REMARK
   The Gauss point is a rather passive object : it does not compute its strains
   and stresses - it just stores them. They are computed by the element. Ifever
   materially nonlinearity is implemented, the role of the Gauss point is 
   likely to grow.
 */
class GaussPoint
{
private :
  int          number ;
  int		      sizeOfVector ;
  Element*     element ;
  Mu::Point*   coordinates ;
  double       weight ;
  FloatArray*  strainVector ;
  FloatArray*  stressVector ;
  FloatArray*  previousStressVector ;
  int          plasticCode ;
  double       stressLevel;
  double       deltaGamma;
  // for hardening purposes
  FloatArray*  backStressVector ;
  FloatArray*  previousBackStressVector ;
  double	      yieldRadius;
  double	      previousYieldRadius;

public :
  GaussPoint  (Element*,int,Mu::Point*,double,int) ;    //<! constructor
  ~GaussPoint  () ;                                      //<! destructor

  double       giveCoordinate (int i)const    { return (*coordinates)[i] ;}
  Mu::Point*   giveCoordinates ()const        { return coordinates ;}
  FloatArray*  giveStrainVector ()            { return strainVector ;}
  FloatArray*  giveStressVector ()				 { return stressVector ;}
  FloatArray*  givePreviousStressVector ()	 { return previousStressVector ;}
  FloatArray*  giveBackStressVector ()			 { return backStressVector ;}
  FloatArray*  givePreviousBackStressVector (){ return previousBackStressVector ;}
  double       giveWeight ()               { return weight ;}
  int          giveSize ()                 { return sizeOfVector ;}
  void         letStrainVectorBe (FloatArray* v);
  void         letStressVectorBe (FloatArray* v);
  void         letPreviousStressVectorBe (FloatArray* v) { previousStressVector = v ;}
  void         letBackStressVectorBe (FloatArray* v);
  void         letPreviousBackStressVectorBe (FloatArray* v) { previousBackStressVector = v ;}
  void         printOutput (FILE*) ;
  void         printBinaryResults(FILE*) ;
  void         getBinaryRecord(float*);
  void         updateYourself () ;
  void		   isPlastic() {plasticCode = 1;} ;
  int          givePlasticCode()const {return plasticCode;} ;
  double       giveStressLevel() {return stressLevel;} ;
  void         computeStressLevel() ;
  double	      giveDeltaGamma()const {return deltaGamma;} ;
  int		      giveNumber()const {return number;} ;
  void		   setDeltaGamma(double aDouble) { deltaGamma = aDouble ;};
  double	      giveYieldRadius() {return yieldRadius;} ;
  void		   setYieldRadius(double aDouble) { yieldRadius = aDouble ;};
  double	      givePreviousYieldRadius() {return previousYieldRadius;} ;
  void		   setPreviousYieldRadius(double aDouble) { previousYieldRadius = aDouble ;};
  /*/
   compute the equivalent von Mises stress. 
	NVP 2005, XFEM implementation
	sigma_e = \sqrt{3}(II_s)^{1/2}
   */
  double       computeVonMisesStress();
  /*!
   Update components of GPs at the end of each step.
	In present implementation, step 2 is made alike to step 1.
	2005-09-13. This could be not efficient.
   */
  void         updateYourselfForXFEM();
  /*!
   Returns the element to which the GaussPoint is associated
   */
  Element*     giveMyElement(){return element;}
} ;

#endif //_GAUSSPNT_H_
