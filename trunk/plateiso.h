//   **************************
//   *** CLASS PlaneProblem       ***
//   **************************

#ifndef _PLATEISO_H_
#define _PLATEISO_H_

#include "element.h"
#include "fei2dquadlin.h"
#include "planelast.h"
#include "enrichmentitem.h"
#include <vector>



//! An isoparametric quadrilateral elements for bending plate
/*!
This class implements a 2D isoparametric elements for bending plate according 
to Mindlin-Reissener theory. 
Each node has three DOFs: tranverse displacement w and two rotations theta_x and theta_y
Selective reduced integration is used to avoid shear locking:
 + Bending stiffness matrix  : 
 + Shearing stiffness matrix :
Remarks: for enriched elements (cracked plate), attention with the integration scheme!!!
Start implemented: 1-12-2005 by Nguyen Vinh Phu, some portions are taken from work
of Mourad Belgasmia.
*/

class PlateIso : public Element
{
public :

  PlateIso (int,Domain*) ;              //!< constructor
  ~PlateIso ()  {;}                     //!< destructor

  /** @name  Numerical integration of stiffness matrix
  *  Functions to perform the numerical integration of stiffness matrix
  */
  //@{
  FloatMatrix*         ComputeBuMatrixAt (GaussPoint*) ;
  FloatMatrix*         ComputeBBmatrixAt (GaussPoint*) ; // bending stiffness matrix
  virtual FloatMatrix* ComputeBSmatrixAt (GaussPoint*) ; // shearing stiffness matrix
  FloatMatrix*         ComputeNmatrixAt (GaussPoint* );
  void                 computeGaussPoints () ;
  virtual void         computeGaussPointsS () ;          // Gauss quad used for shear terms
  double               computeVolumeAround (GaussPoint*) ;
  //@}

  /** @name  Constitutive matrices
  *  Functions to compute the total constitutive matrix from its contributions 
  */
  //@{
  FloatMatrix*       computeConstitutiveMatrix() ;
  FloatMatrix*       computeConstitutiveMatrixS() ;
  FloatMatrix*       computeConstitutiveMatrixB() ;
  //@}

  /** @name  Stiffness matrix
  *  Functions to compute the total stiffness matrix from its contributions 
  */
  //@{
  FloatMatrix*       computeStiffnessMatrix();            // [Ke]
  FloatMatrix*       computeBendingStiffnessMatrix();     // [Kb]
  FloatMatrix*       computeShearingStiffnessMatrix();    // [Ks]
  //@}
  
  /*!
	Set the enrichment for the receiver.
   */
  void               setEnrichmentForMyNodes(EnrichmentItem*);
  /*!
  Partition the element into subtriangles for numerical integration
  @returns the vector containing the subtriangles which defined in 
  local coordinate of the parent element.
  \todo : will be implemented in the near future
  */
  std::vector<DelaunayTriangle *>*  PartitionMySelf();
  
  /*!
  Computes the area of portion above the enrichment item
  */
  double   computeAreaAboveEnrItem(EnrichmentItem *enrItem);
  
private:
  GaussPoint**  gaussPointArrayS ;
} ;

#endif // _PLATEISO_H_
