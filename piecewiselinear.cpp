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

#include "piecewiselinear.h"

#include "vertex.h"
#include "flotarry.h"
#include "geometry_2D.h"
#include "node.h"
#include <list>
#include <map>
#include <vector>
#include <math.h>
#include <algorithm>
#include <typeinfo>
#include <limits>
#include <iostream>


PiecewiseLinear :: PiecewiseLinear(int n, Domain* aDomain)
:GeometryEntity(n,aDomain)
{
  vertexList       = NULL ;
  segmentList      = NULL ;
  //myBoundingBox    = NULL ;
}

PiecewiseLinear::~PiecewiseLinear()
{
  delete vertexList    ;
  delete segmentList   ;
  //delete myBoundingBox ;
}

std::list<Mu::Point*>* PiecewiseLinear :: getListOfVertices()
// **********************************************************
// read from the input file, the vertices of the receiver.
{
  size_t numberOfVertices = this->readInteger("numOfVertices") ;

  vertexList = new std::list<Mu::Point*>;

  for (size_t i = 1 ; i <= numberOfVertices ; i++)
  {
	 size_t n = this->readInteger("vertices",i) ;
	 Vertex *aVertex = static_cast<Vertex*>(domain->giveGeoEntity(n)) ;
	 aVertex->myParent = this; // parent of ith Vertex is this (PiecewiseLinear)
	 Mu::Point *p = new Mu::Point(aVertex->giveCoordinate(1),aVertex->giveCoordinate(2));
	 vertexList->push_back(p)  ;
  }

  return vertexList ;
}

std::list<Mu::Point*>* PiecewiseLinear :: giveMyListOfVertices()
// *************************************************************
// return the vertexList of the receiver. Creates it if it doesn't exist yet
{
  if(vertexList == NULL)
	 vertexList = this->getListOfVertices();

  return vertexList ;
}

Mu::Segment PiecewiseLinear :: FindSegmentClosestTo(Mu::Point* p)
// *************************************************************************
// For a given point, find the segment of the receiver which is closest to p
{
  std::map<double, Mu::Point*>  myMap;

  vertexList = this->giveMyListOfVertices() ;

  // if this polyline composed of 2 vertices, then return itself, otherwise, continue...
  if (vertexList->size() == 2)
  {
	 return Mu::Segment(**vertexList->begin(),**++vertexList->begin());
  }

  for(std::list<Mu::Point*> ::iterator i = vertexList->begin() ; i != vertexList->end(); i++)
  {
	 double x1 =(*i)->x ;       // the x coord. of ith vertex
	 double y1 =(*i)->y ;       // the y coord. of ith vertex
	 double d  =(x1-p->x)*(x1-p->x) + (y1-p->y)*(y1-p->y) ; // square distance between two points
	 myMap[d] = *i ;
  }

  // corrected 2005-08-10
  std::map<double, Mu::Point*>::iterator bIter = myMap.begin();
  double x0 = ((*bIter).second)->x  ;
  double y0 = ((*bIter).second)->y  ;

  ++bIter ;

  double x1 = ((*bIter).second)->x  ;
  double y1 = ((*bIter).second)->y  ;

  return Mu::Segment(Mu::Point(x0,y0),Mu::Point(x1,y1));

}
/*   NO NEED THIS METHOD , NEED FOR VECTOR LEVEL SET METHOD
FloatArray* PiecewiseLinear::FindClosestProjectPointOf(FloatArray* nodeCoord)
{
pair<Vertex*, Vertex*> segment ;
double x1,y1,x2,y2,x3,y3,d,u;
FloatArray* answer;

x3 = nodeCoord->at(1) ;
y3 = nodeCoord->at(2) ;

segment = this ->FindSegmentClosestTo(nodeCoord) ;
x1 = segment.first->giveCoordinate(1) ;
y1 = segment.first->giveCoordinate(2) ;
x2 = segment.second->giveCoordinate(1) ;
y2 = segment.second->giveCoordinate(2) ;

d = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
u = ((x3-x1)*(x2-x1) + (y3-y1)*(y2-y1))/d;

answer->at(1) = x1 + u*(x2-x1);
answer->at(2) = y1 + u*(y2-y1);

return answer ;

} */

double PiecewiseLinear::computeSignedDistanceOfPoint(Mu::Point* coord)
// *********************************************************************
// compute the signed distance from a point whose coord. given by coord
// to the receiver.
{
  int sign = 0 ;

  Mu::Segment segment = this->FindSegmentClosestTo(coord);

  double x1 = segment.first()->x ;
  double y1 = segment.first()->y ;
  double x2 = segment.second()->x ;
  double y2 = segment.second()->y ;

  double delta = (x1-coord->x)*(y2-coord->y) - (x2-coord->x)*(y1-coord->y) ;

  double tol = std::numeric_limits<double>::epsilon();
  if (delta > tol)
	 sign =  1;
  else if (delta < -tol)
	 sign = -1;

  // base P(b) of the perpendicular from P to the extended line L
  // may be outside the range of the segment. In this case, the actual
  // shortest distance is from the point P to one of the endpoints of segment.

  Mu::Point v(x2 - x1, y2 - y1);
  Mu::Point w(coord->x - x1, coord->y - y1);

  double c1 = w*v;

  if ( c1 <= 0 )
	 return (w.norm())*sign;

  double c2 = v*v;
  Mu::Point w1(coord->x - x2, coord->y - y2) ;

  if ( c2 <= c1 )
	 return (w1.norm())*sign;

  // base P(b) of the perpendicular from P to the extended line L
  // is inside the range of the segment

  double b = c1 / c2;

  Mu::Point Pb(x1 + b * v.x, y1 + b * v.y) ;
  double d = dist(Pb, (*coord)) ;

  return d*sign;

}

int PiecewiseLinear:: givePositionComparedTo(Mu::Point* coord)
// ***********************************************************
// Using the "Orientation test" to detect the position of point
// coord w.r.t the segment.
// Result : 1(above), -1(below) and 0(on).
{
  double const EPSILON = 0.000001;

  Mu::Segment segment = this->FindSegmentClosestTo(coord);

  double x1 = segment.first()->x ;
  double y1 = segment.first()->y ;
  double x2 = segment.second()->x ;
  double y2 = segment.second()->y ;

  double delta ;
  if(x2 > x1)
    delta = (x1-coord->x)*(y2-coord->y) - (x2-coord->x)*(y1-coord->y) ;
  else
    delta = (x2-coord->x)*(y1-coord->y) - (x1-coord->x)*(y2-coord->y) ;

  if (delta > EPSILON)        // coord is above the PLine
	 return 1;
  else if (delta < -EPSILON)  // coord is below the PLine
	 return -1;
  else
	 return 0;                 // coord is on the PLine
}

std::vector<Mu::Segment*>* PiecewiseLinear::makeSegments()
// *******************************************************
// from the vertices define the connected segments so that
// can use routines of Mu::Segment
// 2005-09-08 :
// Corrected version for multiple segments ( than 2 vertices, for curved cracks )
{
  segmentList = new std::vector<Mu::Segment*>;

  this->giveMyListOfVertices();

  if (vertexList->size() == 2)
  {
	 Mu::Segment* seg = new Mu::Segment(**vertexList->begin(),**++vertexList->begin());
	 segmentList->push_back(seg);
	 return segmentList;
  }

  for(std::list<Mu::Point*> ::iterator i = vertexList->begin(); i != --vertexList->end() ; i++)
  {
	 std::list<Mu::Point*> ::iterator j = i ;
	 Mu::Segment* seg = new Mu::Segment(**i,**++j);
	 segmentList->push_back(seg);
  }
  /*
  // Debug only
  // -------------------------------------------------------------
  for(size_t k = 0 ; k < segmentList->size() ; k++)
  {
	 (*segmentList)[k]->print();
	 std::cout << std::endl ;
  }
  // -------------------------------------------------------------
  */
  return segmentList;
}

std::vector<Mu::Segment*>* PiecewiseLinear::giveSegments()
// *******************************************************
// return the segments of the receiver. Creates them if they do not exist yet.
{
  if (!segmentList)
	 this->makeSegments();

  return segmentList;
}

bool PiecewiseLinear::intersects(const Mu::Triangle* t)
// ****************************************************
// Check the intersection between the receiver and the triangle t
// Corrected 2005-08-10 for multiple segments crack !!!
{
  this->giveSegments();
  std::vector<Mu::Segment*>::iterator i;
  for(i = segmentList->begin() ; i != segmentList->end(); i++)
  {
	 if((*i)->intersects(t))
		return true ;
  }
}

bool PiecewiseLinear::intersects(const Mu::Segment* s)
// ***************************************************
// check if the receiver intersects with Segment s or not
// using Mu::Segment::intersects(s)
{
  this->giveSegments();
  std::vector<Mu::Segment*>::iterator i;
  for(i = segmentList->begin() ; i != segmentList->end(); i++)
  {
	 if((*i)->intersects(s))
		return true ;
  }
}

std::vector<Mu::Point> PiecewiseLinear::intersection(const Mu::Triangle* t)
// ************************************************************************
{
  std::vector<Point> ret ;
  std::vector<Point> result ;
  std::vector<Mu::Segment*> ::iterator i;

  std::vector<Mu::Segment*>* segmentList = this->giveSegments();

  for(i = segmentList->begin(); i != segmentList->end(); i++)
  {
	 if ( (*i)->intersects(t) )      // this segment cuts the triangle
	 {
		ret = (*i)->intersection(t) ; // compute the intersections
		result.insert(result.end(),ret.begin(),ret.end());
	 }
  }

  // check for case the segment goes throught the corners of the triangle t
  // corrected 2005-08-09 when got trouble with structured triangle mesh.

  if(result.size() == 1)
  {
	 for(i = segmentList->begin(); i != segmentList->end(); i++)
	 {
		if((*i)->on(t->getBoundingPoint(0)))
		  result.push_back(*(t->getBoundingPoint(0)));
		if((*i)->on(t->getBoundingPoint(1)))
		  result.push_back(*(t->getBoundingPoint(1)));
		if((*i)->on(t->getBoundingPoint(2)))
		  result.push_back(*(t->getBoundingPoint(2)));
	 }
  }

  return result ;
}

std::vector<Mu::Point> PiecewiseLinear::intersection(const Mu::Segment* s)
// ***********************************************************************
{
  std::vector<Point> result ;
  std::vector<Mu::Segment*> ::iterator i;

  std::vector<Mu::Segment*>* segmentList = this->giveSegments();

  for(i = segmentList->begin(); i != segmentList->end(); i++)
  {
	 if ( (*i)->intersects(s) )      // this segment cuts the triangle
	 {
		result.push_back((*i)->intersection(s));
	 }
  }

  return result ;
}

Mu::Segment* PiecewiseLinear::giveSegmentContain(Vertex* tip)
// ***********************************************************
// give segment contains the crack tip. If tip is the first tip, then
// need to change the direction of the segment.
// The tip segment is VERY IMPORTANT : used to compute the inclination angle of the
// crack, used for the tip to do the element partitioning.
{

  std::vector<Mu::Segment*> :: iterator i = segmentList->begin();
  std::vector<Mu::Segment*> :: iterator j = segmentList->end();

  double x = tip->giveCoordinate(1);
  double y = tip->giveCoordinate(2);
  Mu::Point* point = new Mu::Point(x,y);

  // check first segment
  if ( *point == *((*i)->first())) // point is the first TIP
  {
	 //std::cout << " This is first TIP : " << std::endl ; // debug only
	 //(*i)->first()->print();  std::cout << std::endl ;
	 Mu::Segment * seg = new Mu::Segment(*((*i)->second()),*((*i)->first()));
	 delete point ;
	 return seg ;
  }

  if ( *point == *((*i)->second()))  // point is the second TIP
  {
	 Mu::Segment * seg = new Mu::Segment(*((*i)->first()),*((*i)->second())); // return a copy !!!
	 delete point ;
	 return seg ;
  }

  // check last segment
  std::vector<Mu::Segment*> :: iterator k = j ;
  Mu::Segment *s = *--k ; // last segment
  if ( *point == *((*(--j))->second()))
  {
	 Mu::Segment *seg = new Mu::Segment( *( s->first() ),*( s->second() ) );
	 delete point ;
    return seg ;
  }
}

void PiecewiseLinear :: exportToMatlab(string& theString)
// *****************************************************
// format of file :
// crack = [x1 y1 x2 y2 ... xn yn]
{
  char value[30];
  string space(" ");
  string newline("\n");
  double x,y ;
  for (std::list<Mu::Point*> ::iterator i = vertexList->begin() ; i != vertexList->end(); i++)
  {
	 x =(*i)->x ;  // the X coord. of ith vertex
	 y =(*i)->y ;  // the Y coord. of ith vertex

	 _gcvt(x,6,value);     // transforms the double param1 in char param3 with 3 digits
	 theString += value;
	 theString += space;
	 _gcvt(y,6,value);
	 theString += value;
	 theString += space;
	 theString += newline;
  }

  theString += newline;
}

void PiecewiseLinear :: insertNewVertexInFront(Vertex *v)
// *****************************************************
{

  Mu::Point *p = new Mu::Point(v->giveCoordinate(1),v->giveCoordinate(2));
  vertexList->push_front(p);
}

void PiecewiseLinear :: insertNewVertexAtBack(Vertex *v)
// *****************************************************
{
  Mu::Point *p = new Mu::Point(v->giveCoordinate(1),v->giveCoordinate(2));
  vertexList->push_back(p);
}

bool PiecewiseLinear :: isBelongToMe(Mu::Point* p)
{
  this->giveSegments();
  bool found = false;
  for(size_t i = 0 ; i < segmentList->size() ; i++)
  {
	 if ((*segmentList)[i]->on(p))
	 {
		found = true ;
		i = segmentList->size() ;
	 }
  }
  return found ;
}

/*
void  PiecewiseLinear :: makeBoundingBox()
//****************************************
// just consider the case of two vertices !!!
{
  vertexList = this->giveMyListOfVertices();

  double x0 = (*vertexList->begin())->giveCoordinate(1);
  double y0 = (*vertexList->begin())->giveCoordinate(2);

  double x1 = (*++vertexList->begin())->giveCoordinate(1);
  double y1 = (*++vertexList->begin())->giveCoordinate(2);

  double originX = 0.5 * (x0 + x1);
  double originY = 0.5 * (y0 + y1);
  double length = x1 - x0 ;

  double height;
  // check if the receiver is horizontal or oblique
  if(y0 == y1)
	 height = 1.5 ;     // NEED BE IMPROVED LATER !!!
  else
	 height = y1 - y0 ;

  myBoundingBox = new Mu::Rectangle(length,height,originX,originY);

}

Mu::Rectangle*  PiecewiseLinear::giveBoundingBox()
//************************************************
{
  if(myBoundingBox == NULL)
	 this->makeBoundingBox() ;

  return myBoundingBox ;
}

std::vector<Element*> PiecewiseLinear :: giveElementsAroundMe()
//*************************************************************
{
  Mu::Rectangle* myBoundingBox = this->giveBoundingBox();

  std::map<Node*,std::vector<Element*> > nodeElemMap = this->domain->giveNodalSupports();

  std::vector<Element*> ret ;

  for(size_t i = 0 ; i < domain->giveNumberOfNodes() ; i++)
  {
	 Node *node = domain->giveNode(i+1) ;
	 if (myBoundingBox->in(*node->makePoint()))
	 {
		ret.insert(ret.end(),nodeElemMap[node].begin(),nodeElemMap[node].end()) ;
	 }
  }

  std::sort(ret.begin(),ret.end());
  std::unique(ret.begin(),ret.end());

  return ret ;
}*/

