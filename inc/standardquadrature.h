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

#ifndef _STANDARDQUADRATURE_H_
#define _STANDARDQUADRATURE_H_

#include "integrationrule.h"
#include "geometry_2D.h"



//! Standard Gauss-Legendre Quadrature
/*!
 Class representing Gaussian-quadrature integration rule. 
 The number of integration points and their coordinates and integration weights depends on
 integration rule type (rule for integration in 1d, 2d, 3d) and required  acurracy.\n
DESCRIPTION:
   Implements integration rule class.
  Stores integration points used for integration
  of necesary terms (for example computation of  stiffness matrix 
  or computation of element nodal force vector ) 
  and it  corresponds to some local strains 
  on finite element level. Finite element can have many 
  integration rules corresponding to  different strains.\n

TASKS:
  returning number of integration points used
  returning requested integration point - method getIntegrationPoint
  */
class StandardGaussLegendreQuadrature : public IntegrationRule
{

public:

  StandardGaussLegendreQuadrature (){}
  ~StandardGaussLegendreQuadrature(){}
 /*!
  Returns requred number of integration points to exactly integrate
  polynomial of order approxOrder on given domain.
  When approxOrder is too large and is not supported by implementation
  method returns -1.
  */
 int getRequiredNumberOfIntegrationPoints (IntDomain dType, int approxOrder) ;

protected:
 /*!
  Sets up receiver's  integration points on unit line integration domain.
  @returns the array of Gauss points 
  */
  void  SetUpPointsOnLine    (size_t) ;
 /*!
  Sets up receiver's  integration points on triangular (area coords) integration domain.
  @returns the array of Gauss points
  */
	void  SetUpPointsOnTriangle (size_t) ;
 /*!
  Sets up receiver's  integration points on unit square integration domain.
  @returns the array of Gauss points
  */
	void  SetUpPointsOnSquare  (size_t) ;
 /*!
  Sets up receiver's  integration points on tetrahedra domain.
  @returns the array of Gauss points
  */
	void  SetUpPointsOnTetrahedra (size_t);
  
};

#endif //_STANDARDQUADRATURE_H_
