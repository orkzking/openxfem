//   **************************
//   *** CLASS QUAD_U       ***
//   **************************

#ifndef _QUAD_U_H_
#define _QUAD_U_H_

#include "planeproblem.h"
#include "fei2dquadlin.h"
#include "enrichmentitem.h"
#include <vector>



//! 4 noded isoparametric quadrilateral elements for plane elasticity
/*!
This class implements an isoparametric four-node quad plane-
strain elasticity finite element. Each node has 2 degrees of freedom.
DESCRIPTION :
One single additional attribute is needed for Gauss integration purpose :
'jacobianMatrix'. This 2x2 matrix contains polynomials.
TASKS :
- calculating its Gauss points ;
- calculating its B,D,N matrices and dV.
*/

class Quad4_U : public PlaneProblem
{
public :

  Quad4_U (int,Domain*) ;             //!< constructor
 ~Quad4_U ()  {;}                     //!< destructor

  void               computeGaussPoints () ;
  /**
  Define the FE interpolation used for standard approximation
  */
  FEI2dQuadLin*      giveFEInterpolation();
  /**
  Define the FE interpolation used for for enriched approximation
  */
  FEI2dQuadLin*      giveXFEInterpolation();
  
  /*!
  Partition the element into subtriangles for numerical integration
  @returns the vector containing the subtriangles which defined in 
  local coordinate of the parent element.
  \todo : will be implemented in the near future
  */
  std::vector<DelaunayTriangle *>*  PartitionMySelf();
  /*!
  Computes the area of the element
  */
  double   area();
  /*!
  Computes the area of portion above the enrichment item
  */
  double   computeAreaAboveEnrItem(EnrichmentItem *enrItem);
  GaussPoint**  setGaussQuadForJ_Integral();

} ;

#endif // _QUAD_U_H_
