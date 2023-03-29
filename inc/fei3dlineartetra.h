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

#ifndef _FEI3DLINEARTETRA_H_
#define _FEI3DLINEARTETRA_H_

#include "feinterpol3d.h"
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
class FEI3dLinearTetrahedra:public FEInterpolation3d
{
public:
   /*! Constructor
	 */
  FEI3dLinearTetrahedra():FEInterpolation3d (1){}
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
	 *  Compute the shape function w.r.t parent coordinate system (xi,eta,zeta)
	 *  at Gauss points
	 */
   //@{
	FloatArray*   giveDerivativeKsi ();
	FloatArray*   giveDerivativeEta (); 
	FloatArray*   giveDerivativeZeta();
	//@}
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
	
};


#endif  // _FEI3DLINEARTETRA_H_
