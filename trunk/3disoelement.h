//   **************************
//   *** CLASS PlaneProblem ***
//   **************************

#ifndef _3DISOELEMENT_H_
#define _3DISOELEMENT_H_

#include "element.h"
#include "enrichmentitem.h"


//! A 3D isoparametric elements for three dimensional problems
/*!
This class implements a geneal 3D isoparametric finite element. 
Each node has 3 degrees of freedom: u,v,w.
This is the super class of all classes implementing the four-noded tetrahedron element,
ten-noded tetrahedron elements...
*/

class Iso3dElement : public Element
{
public :

  Iso3dElement(int,Domain*) ;              //!< constructor
 ~Iso3dElement ()  {;}                     //!< destructor

  /** @name  Numerical integration of stiffness matrix
  *  Functions to perform the numerical integration of stiffness matrix
  */
  //@{
  FloatMatrix*       ComputeBuMatrixAt (GaussPoint*) ;
  FloatMatrix*       ComputeNmatrixAt (GaussPoint* );
  double             computeVolumeAround (GaussPoint*) ;
  //@}
  /**
  Compute the constituitive matrix 
  */
  FloatMatrix*       computeConstitutiveMatrix(); 
  /*!
	Set the enrichment for the receiver.
   */
  void               setEnrichmentForMyNodes(EnrichmentItem*);
 
} ;

#endif // _3DISOELEMENT_H_
