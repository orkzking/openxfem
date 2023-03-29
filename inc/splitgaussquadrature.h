/************************************************************************************ 

   Copyright (C) 2005
   Stephane BORDAS, Cyrille DUNANT, Vinh Phu NGUYEN, Quang Tri TRUONG, Ravindra DUDDU

   This file is part of the XFEM C++ Library (XFEMLIB) written 
   and maintained by above authors.

   This program is free software; you can redistribute it and/or modify it.

   This program is distributed in the hope that it will be useful, 
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.
   
	Any feedback is welcome. Emails : nvinhphu@gmail.com, ...

*************************************************************************************/

#ifndef _SPLITGAUSSQUADRATURE_H_
#define _SPLITGAUSSQUADRATURE_H_

#include "integrationrule.h"
#include "standardquadrature.h"
#include "feinterpol.h"
#include "element.h"


//! Gauss quadrature for discontinuous functions
/*!
 This class implementes a discontinuous Gauss Quadrature
 used to numerically integrate the weak form in elements cut
 by the discontinuities.
DESCRIPTION:
 First, we compute the intersections betwwen the element edges
 and the discontinuities. Then, the Delaunay triangulation is applied
 to get the subtriangles. The Gauss points of each subtriangle are converted to
 those defined in parent T3 element.
TASKS:
  returning number of integration points used
  */
class SplitGaussLegendreQuadrature : public IntegrationRule
{
public:

  SplitGaussLegendreQuadrature  (){};//!< Constructor
  ~SplitGaussLegendreQuadrature (){};//!< Destructor

  /**
	Initializes the receiver. Receiver integration points are created acording to given parameters.
	@param elem describes the finite element
	@param nPoints is the necessary number of Gauss Points for each sub-triangle.
	*/
  void setUpIntegrationPoints (Element* elem,size_t nPoints);

  /*! 
	Sets up receiver's  integration points on 2D integration domain.
	@returns the array of Gauss points
	*/
  void  SetUpPointsOn2dDomain(Element*,size_t nPoints);
 
  /*! 
	Sets up receiver's  integration points on 3D integration domain.
	@returns the array of Gauss points
	\todo 
	*/
  void  SetUpPointsOn3dDomain(Element*,size_t nPoints){} 
};

#endif //_SPLITGAUSSQUADRATURE_H_
