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
   
	Any feedback is welcome. Emails : nvphu80@yahoo.com, ...

*************************************************************************************/
#ifndef _FEINTERPOL_H_
#define _FEINTERPOL_H_

#include "flotmtrx.h"
#include "flotarry.h"
#include "node.h"
#include "domain.h"
#include "geometry_base.h"


//! Abstract class, interface for FE interpolation
/*!
 This class is an abstract class, the superclass of all classes that imple-
 ment the finite element interpolations such as FE2dquadlin, FE2dtrilin...
 */
class FEInterpolation
{
public:
   FEInterpolation(int o):order(o){}    //!< constructor
	virtual ~FEInterpolation(){}         //!< destructor

	int giveInterpolationOrder()const {return order;}
	/*!
	 Computes the shape functions at point with coordinates
	 given by lcoords
	 */
	virtual FloatArray* evalN(Mu::Point* lcoords)=0;
	/*!
	 Computes the derivatives of shape functions at point with coordinates
	 given by lcoords
	 */
	virtual FloatMatrix* evaldNdx(Domain* d,IntArray* nodes,Mu::Point *lcoords)=0;
	/*!
	 Compute the Jacobian matrix 
	 */
	virtual FloatMatrix* giveJacobianMatrixAt (Domain*,IntArray*,Mu::Point*)=0 ;
	/*!
    Returns the global coordinates of a point having local coordinates
	 Used to compute the global coord. of a Gauss point.
	 */
	virtual Mu::Point*   local2Global(Domain*,IntArray*,Mu::Point* )=0;
	/*!
    Returns the local coordinates of a point having global coordinates
	 */
	virtual Mu::Point*   global2Local(Domain* d,IntArray* nodes,Mu::Point* p)=0;

protected:
	int order;
};


#endif //_FEINTERPOL_H_
