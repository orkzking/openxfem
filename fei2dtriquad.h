// **************************
// *** CLASS FEI2DTRIQUA  ***
// **************************

#ifndef _FEI2DTRIQUA_H_
#define _FEI2DTRIQUA_H_

#include "feinterpol2d.h"
#include "flotarry.h"
#include "domain.h"
#include "intarray.h"
#include "geometry_base.h"
#include <stdio.h>
class Node;


//! FE interpolation for six noded triangle elements
/*!
 This class implementes the quadratic finite interpolation for
 six noded triangle elements.
 */
class FEI2dTriQuad : public FEInterpolation2d
{
public:
   /*! Constructor
	 */
	FEI2dTriQuad (int ind1, int ind2) : FEInterpolation2d (2) {xind = ind1; yind=ind2;}
   /*!
	 Compute the shape function matrix at Gauss points
	 */
	FloatArray*  evalN(Mu::Point*);
	/*!
	 Compute the derivatives of shape functions w.r.t physical coord. 
	 at Gauss points
	 */
	FloatMatrix* evaldNdx(Domain*,IntArray*,Mu::Point*);
	/** @name Derivations of shape functions w.r.t local coord. system
	 *  Compute the shape function w.r.t parent coordinate system (xi,eta)
	 *  at Gauss points
	 */
   //@{
	FloatArray*  giveDerivativeKsi (double,double);
	FloatArray*  giveDerivativeEta (double,double);
	//@}
	/*!
	 Compute the Jacobian matrix 
	 */
	FloatMatrix* giveJacobianMatrixAt (Domain*,IntArray*,Mu::Point*) ;
	/*!
	 Compute the global coordinates of a point. 
	 */
	Mu::Point*   local2Global(Domain*,IntArray*,Mu::Point*);
	/*!
	 Compute the local coordinates of a point. 
	 */
	Mu::Point*   global2Local(Domain* d,IntArray* nodes,Mu::Point* p);
protected:
	int xind, yind;
};


#endif // _FEI2DTRIQUA_H_
