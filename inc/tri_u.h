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

#ifndef _TRI_U_H_
#define _TRI_U_H_

#include "planeproblem.h"
#include "fei2dtrilin.h"
#include "geometry_2D.h"
#include <vector>

//! 3 noded isoparametric triangle elements for plane elasticity
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

class Tri3_U : public PlaneProblem
{
public :

  Tri3_U (int,Domain*) ;    //!< constructor
  ~Tri3_U ()  {}            //!< destructor

  void               computeGaussPoints () ;

  /**
  Define the FE interpolation used for standard approximation
  */
  FEI2dTriLin*       giveFEInterpolation();
  /**
  Define the FE interpolation used for for enriched approximation
  */
  FEI2dTriLin*       giveXFEInterpolation();
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
  Computes the area of the element
  */
  double   area();
  /*!
  Make a triangle from three nodes of the receiver so that can use Cyrille's facilities
  */
  Mu::Triangle*     makeTriangle(); 
  GaussPoint**  setGaussQuadForJ_Integral();
} ;

#endif //_TRI_U_H_
