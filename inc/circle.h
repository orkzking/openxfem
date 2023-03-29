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
   
	Any feedback is welcome. Emails : nvphu80@yahoo.com,stephane.bordas@epfl.ch

*************************************************************************************/
 
#ifndef  _CIRCLE_H_
#define  _CIRCLE_H_

#include "geometryentity.h"
#include "vertex.h"
#include "flotarry.h"
#include "domain.h"
#include "geometry_2D.h"
#include <iostream>


/*! A Circle is geometry of holes
 This class implements a circle with center P and radius R
 */


class Cercle:public GeometryEntity
{
public:

  /*!
   Constructor
   */
   Cercle(int n,Domain* d):GeometryEntity(n,d),radius(0),center(0){};
	~Cercle(){;}           //!< Destructor

	/*!
	 Return the center of the receiver 
	 */
	Vertex*    giveCenter();
	/*!
	 Return the radius of the receiver 
	 */
	double     giveRadius();
	/*!
	 Check if a point defined by coord is inside (negative) or
	 outside the circle (positive value)
	 */
	double     computeSignedDistanceOfPoint(Mu::Point*);
	/*!
	 Returns the position of the point compared to the receiver
	 @return 1 if point is inside and @return 0 otherwise
	 */
	int        givePositionComparedTo(Mu::Point*);
	/*!
	 Check the intersection between circle and triangle
	 @return true if intersects 
	 */
	bool       intersects(const Mu::Triangle* t);
	bool       intersects(const Mu::Segment* s );
	/*!
	 Computes the intersection points of circle and triangle
	 @param t the triangle 
	 @return intersection points stored in vector<Mu::Point>
	 */
	std::vector<Mu::Point> intersection(const Mu::Triangle* t);
   std::vector<Mu::Point> intersection(const Mu::Segment* s) ;
	/*!
    Export geometry to Matlab file for plotting.
	 Not yet implemented at this point.\todo
    */
	void  exportToMatlab(std::string&);

	
private:
	int      center ;
	double   radius ;
} ;

#endif // _CIRCLE_H_