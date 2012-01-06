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

#ifndef _FEI2DQUADLIN_H_
#define _FEI2DQUADLIN_H_

#include "feinterpol2d.h"
#include "flotarry.h"
#include "domain.h"
#include "intarray.h"
#include "geometry_base.h"
#include <stdio.h>
class Node;


//! FE interpolation for 4 node quadrilateral elements
/*!
 This class implementes the finite element interpolation for 
 one 4 nodes quadrilateral element.
 */
class FEI2dQuadLin:public FEInterpolation2d
{
public:
   /*! Constructor
	 */
	FEI2dQuadLin (int ind1, int ind2) : FEInterpolation2d (1) {xind = ind1; yind=ind2;}
	/*!
	 Compute the shape function matrix at Gauss points
	 */
	FloatArray*   evalN(Mu::Point*);
	/*!
	 Compute the derivatives of shape functions w.r.t physical coord. 
	 at Gauss points
	 */
	FloatMatrix*  evaldNdx(Domain*,IntArray*,Mu::Point*);
	/** @name Derivations of shape functions w.r.t local coord. system
	 *  Compute the shape function w.r.t parent coordinate system (xi,eta)
	 *  at Gauss points
	 */
   //@{
	FloatArray*   giveDerivativeKsi (double=0,double=0);
	FloatArray*   giveDerivativeEta (double=0,double=0); //@}
	/*!
	 Compute the Jacobian matrix 
	 */
	FloatMatrix*  giveJacobianMatrixAt (Domain*,IntArray*,Mu::Point*) ;
	/*!
	 Compute the global coordinates of a point with local coord. 
	 */
	Mu::Point*    local2Global(Domain* d,IntArray*nodes,Mu::Point* p);
   /*!
	 Compute the local coordinates of a point. \TODO 
	 */
	Mu::Point*   global2Local(Domain* d,IntArray* nodes,Mu::Point* p);
	

protected:
	int xind,yind;
};


#endif  // _FEI2DQUADLIN_H_
