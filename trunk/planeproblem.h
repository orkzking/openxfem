//   **************************
//   *** CLASS PlaneProblem ***
//   **************************

#ifndef _PLANEPROBLEM_H_
#define _PLANEPROBLEM_H_

#include "element.h"
#include "planelast.h"
#include "enrichmentitem.h"


//! An 2D isoparametric elements for plane elasticity
/*!
 * This class implements an isoparametric finite element for plane elasticity problem. 
 * This is the super class of the 3 node and 6 node triangle elements, 4 node quadrilateral elements
 * Each node has 2 degrees of freedom which are the displacement along X and Y directions.
 *
 * Tasks:
 *  - Compute the elastic constitutive matrix D
 *  - Compute the discretized strain matrix B, only the classical part, not the enriched part
 *  - Set enrichment for its nodes
 *  - Compute the shape function matrix N
 *  - Compute the detJ*wt*t for the integration transformation
*/

class PlaneProblem : public Element
{
public :

  PlaneProblem (int,Domain*) ;              //!< constructor
  ~PlaneProblem ()  {;}                     //!< destructor

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
 
private:

  PlaneElasticity * planeElast; 

} ;

#endif // _PLANEPROBLEM_H_
