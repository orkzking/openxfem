

#ifndef _MICT4_H_
#define _MICT4_H_

#include "ansplate.h"
#include "fei2dquadlin.h"
#include "planelast.h"
#include "enrichmentitem.h"
#include <vector>


//! A popular element of the family of ANS elements. The MITC-4 element
/*!
 * Inputs: 
 *  - number of nodes : 4 nodes
 *  - number and location of sampling points : 4 sampling points
 *  - shape functions used for interpolating the natural shear strain 
 *    -# \f$N_1 = 0.5*(1-\eta) N_2 = 0.5*(1+\eta)\f$
 *	   -# \f$N_3 = 0.5*(1-\ksi) N_4 = 0.5*(1+\ksi)\f$
 *
 * Remarks: for enriched elements (cracked plate), attention with the integration scheme!!!
 * Start implemented: 1-12-2005 by Nguyen Vinh Phu, some portions are taken from work
 * of Mourad Belgasmia.
*/

class MITC4 : public ANSPlate
{
public :

  MITC4 (int,Domain*) ;              //!< constructor
  ~MITC4 ()  {;}                     //!< destructor

  /**
  Define the FE interpolation used for standard approximation
  */
  FEI2dQuadLin*      giveFEInterpolation();
  /**
  Define the FE interpolation used for for enriched approximation
  */
  FEI2dQuadLin*      giveXFEInterpolation();

  FloatMatrix*  ComputeAssumedStrainShapeFunction(GaussPoint*);
  void          setttingSamplingPoints();
  
  /*!
  Computes the area of the element
  */
  double   area();
  void               computeGaussPoints () ;
  GaussPoint**  setGaussQuadForJ_Integral();
} ;

#endif // _MITC4_H_
