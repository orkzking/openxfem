//  *********************************************
//   ***       CLASS STANDARDDESCRIPTION       ***
//   *********************************************


#ifndef _STANDARDDESCRIPTION_H_
#define _STANDARDDESCRIPTION_H_

#include "geometrydescription.h"
#include "flotarry.h"


//! Standard description of geometry entities  
/*!
 The geometry of the discontinuities such as PiecewiseLinear, Circle, Vertex
 is traditionally represented, i.e., to find out enriched nodes, we use the
 geometry predicate, see Sukumar and Prevost 2003
 */

class StandardDescription:public GeometryDescription
{
public:
  StandardDescription(){}             //!< Constructor
  ~StandardDescription(){}           //!< Destructor

  /*!
  check if element interacts with GeometryEntity ot not
  using the pure geometry predication
  */
  bool interactsWith(Element*,GeometryEntity*);
} ;


#endif // _STANDARDDESCRIPTION_H_
