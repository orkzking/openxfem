//   *********************************************
//   ***        CLASS LEVELSETDESCRIPTION       ***
//   *********************************************

 
#ifndef _VECTORLSDESCRIPTION_H_
#define _VECTORLSDESCRIPTION_H_

#include "geometrydescription.h"
#include "vertex.h"
#include "flotarry.h"
#include "domain.h"


//! Vector level set description of geometry entities  
/*!
 The geometry of the discontinuities such as PiecewiseLinear, Circle, Vertex
 is represented by the vector level sets
 see Ventura, Xu and Belytschko 2002
 */

class VectorLevelSetDescription:public GeometryDescription
{

public:
  VectorLevelSetDescription(){}             //!< Constructor
	~VectorLevelSetDescription(){}           //!< Destructor
} ;


#endif // _VECTORLSDESCRIPTION_H_