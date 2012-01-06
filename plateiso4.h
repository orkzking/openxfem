//   **************************
//   *** CLASS PlaneProblem       ***
//   **************************

#ifndef _PLATEISO4_H_
#define _PLATEISO4_H_

#include "plateiso.h"
#include "fei2dquadlin.h"
#include "planelast.h"
#include "enrichmentitem.h"
#include <vector>



//! Four-noded isoparametric quadrilateral elements for bending plate
/*!
This class implements a four-noded isoparametric quadrilateral elements 
for bending plate according to Mindlin-Ressener theory. 
Each node has three DOFs: tranverse displacement w and two rotations theta_x and theta_y
Selective reduced integration is used to avoid shear locking:
 + Bending stiffness matrix  : 2x2 Gauss quadrature rule
 + Shearing stiffness matrix : 1 Gauss quadrature rule
Remarks: for enriched elements(cracked plate), attention with the integration scheme!!!
Start implemented: 1-12-2005 by Nguyen Vinh Phu, some portions are taken from work
of Mourad Belgasima.
*/

class PlateIsoQ4 : public PlateIso
{
public :

  PlateIsoQ4 (int,Domain*) ;              //!< constructor
  ~PlateIsoQ4 ()  {;}                     //!< destructor

  /**
  Define the FE interpolation used for standard approximation
  */
  FEI2dQuadLin*      giveFEInterpolation();
  /**
  Define the FE interpolation used for for enriched approximation
  */
  FEI2dQuadLin*      giveXFEInterpolation();

  void   computeGaussPoints () ;
  void   computeGaussPointsS () ; // Gauss quad used for shear terms
 
  /*!
  Computes the area of the element
  */
  double   area();
  GaussPoint**  setGaussQuadForJ_Integral();

private:
  GaussPoint**  gaussPointArrayS ;
} ;

#endif // _PLATEISO4_H_
