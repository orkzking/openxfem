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
   
	Any feedback is welcome. Emails : nvinhphu@gmail.com, ...

*************************************************************************************/

 
#ifndef _PIECEWISELINEAR_H_
#define _PIECEWISELINEAR_H_

#include "geometryentity.h"
#include "domain.h"
#include "geometry_2D.h"
#include <list>
#include <vector>
#include <iostream>
class Vertex ; class GeometryDescription; class CrackInterior;

//! A Polyline is geometry of a straight crack
/*!
 This class implements a polyline composed by several vertices   
 A PiecewiseLinear is geometry of a straight crack 
 */



class PiecewiseLinear:public GeometryEntity
{

public:
   
	PiecewiseLinear(int,Domain*) ; //!< Constructor
	~PiecewiseLinear(); //!< Destructor

	/** @name List of vertices 
    *  Methods to get and give vertices of the polyline
	 *  Replace Vertex by Mu::Point makes us easier to insert points into the PiecewiseLinear.
    */
   //@{
	
	Vertex*                  giveVertex (int) ;
	std::list<Mu::Point*>*   getListOfVertices();
	std::list<Mu::Point*>*   giveMyListOfVertices();
   //@}

	/*!
    Find the segment of polyline closest to a given point	  
	 returns two vertices of this segment
	 */
	Mu::Segment FindSegmentClosestTo(Mu::Point*);
	/*!
    Compute the signed distance of a given point to
	 the receiver.
	 */
	double       computeSignedDistanceOfPoint(Mu::Point*); 
	 /*!
	  Check the position of a point w.r.t the receiver
	  @return  1 if p is above the receiver
	  @return -1 if p is below the receiver
	  @return  0 if p is on the receiver
	  */
	int          givePositionComparedTo(Mu::Point* p);

	/*!
	 Make a vector of segments
	 */
	std::vector<Mu::Segment*>*   makeSegments();
	/*!
	 Return vector of segments of the receiver
	 */
	std::vector<Mu::Segment*>*   giveSegments();
	/*!
	 Check the intersection of the receiver with a triangle
	 */
	bool   intersects(const Mu::Triangle* tri);
	bool   intersects(const Mu::Segment* s ) ;
	/*!
	 Compute the intersection of the receiver with a triangle
	 @param tri the triangle 
	 @return 
	 */
	std::vector<Mu::Point> intersection(const Mu::Triangle* tri);
   std::vector<Mu::Point> intersection(const Mu::Segment* s) ;
	/*!
	 Returns the extreme segments whose end is the vertex v( often v is the TIP)
	 @param v the Vertex
	 @return segment containing the \c Vertex v
	 */
	 Mu::Segment*          giveSegmentContain(Vertex* v) ;

	 /*!
	  Defines the bounding box of the reveiver
	  Used to "quickly" find out elements interacting with the receiver
	  */
	 //void                  makeBoundingBox();
	 //Mu::Rectangle*        giveBoundingBox();
	 //std::vector<Element*> giveElementsAroundMe();

	 /*!
     Export geometry to Matlab file for plotting.
     */
	 void  exportToMatlab(std::string&);

	 /** @name Methods for updating the geometry of PiecewiseLinear
	  *  1. Update the vertices
	  *  2. Update the segments. 
	  *  2005-09-04
     */
    //@{
	 void  insertNewVertexInFront(Vertex* v);
    void  insertNewVertexAtBack(Vertex* v);
	 void  insertNewSegmentInFront(Mu::Segment *s){segmentList->insert(segmentList->begin(),s);}
	 void  insertNewSegmentAtBack(Mu::Segment *s){segmentList->push_back(s);}
	 //@}

	 /*!
	  Check if a point belongs to the receiver
	  */
	 bool  isBelongToMe(Mu::Point* p);

private:
    std::list<Mu::Point*>*     vertexList ;   //!< list of vertices of the PiecewiseLinear
	 std::vector<Mu::Segment*>* segmentList;   //!< segments of the Piecewiselinear
	 //Mu::Rectangle*             myBoundingBox; //!< bounding box of the PL
    
} ;

#endif // _PIECEWISELINEAR_H_
