/************************************************************************************ 

Copyright (C) 2005
Stephane BORDAS, Cyrille DUNANT, Vinh Phu NGUYEN, Quang Tri TRUONG, Ravindra DUDDU

This file is part of the XFEM C++ Library (OpenXFEM++) written 
and maintained by above authors.

This program is free software; you can redistribute it and/or modify it.

This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
license text for more details.

Any feedback is welcome. Emails : nvphu80@yahoo.com, ...

*************************************************************************************/

#ifndef _TETRA4_H_
#define _TETRA4_H_

#include "3disoelement.h"
#include "fei3dlineartetra.h"
#include "geometry_2D.h"
#include <vector>

//! Four-noded isoparametric tetrahedra elements for 3D elasticity problems.
/*!
This class implements an isoparametric four-noded tetrahedra finite element. 
Each node has 3 degrees of freedom:u,v,w.
DESCRIPTION :
One single additional attribute is needed for Gauss integration purpose :
'jacobianMatrix'. This 2x2 matrix contains polynomials.
TASKS :
- calculating its Gauss points ;
- calculating its B,D,N matrices and dV.
*/

class Tetra4 : public Iso3dElement
{
public :

  Tetra4 (int,Domain*) ;    //!< constructor
  ~Tetra4 ()  {}            //!< destructor

  void               computeGaussPoints () ;

  /**
  Define the FE interpolation used for standard approximation
  */
  FEI3dLinearTetrahedra*      giveFEInterpolation();
  /**
  Define the FE interpolation used for for enriched approximation
  */
  FEI3dLinearTetrahedra*     giveXFEInterpolation();
  /*!
  Partition the element into subtriangles for numerical integration
  @returns the vector containing the subtriangles which defined in 
  local coordinate of the parent element.
  */
  std::vector<DelaunayTriangle *>*  PartitionMySelf();
  /*!
  Computes the area of portion above the enrichment item
  */
  double            computeAreaAboveEnrItem(EnrichmentItem*);
  /*!
  Computes the volume of the element
  */
  double   volume();

  GaussPoint**  setGaussQuadForJ_Integral();

} ;

#endif //_TETRA4_H_
