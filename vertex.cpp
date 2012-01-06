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

#include "vertex.h"

#include "domain.h"
#include "flotarry.h"
#include "node.h"
#include "element.h"



Vertex :: Vertex(int n, Domain* aDomain)
     :GeometryEntity(n,aDomain)
{
  coordinates = NULL ;
  myParent    = NULL ;
}

Vertex::~Vertex()
{
  delete coordinates ;
  //delete myParent ;
}

void Vertex :: getCoordinates()
// ****************************
// Get from the data file all of the coordinates of the receiver.
{
   size_t numberOfCoordinates = this->readInteger("coord") ;
   coordinates = new FloatArray(numberOfCoordinates) ;
   for (size_t i = 1 ; i <= numberOfCoordinates ; i++)
	{
	   coordinates->at(i) = this->read("coord",i+1) ;
	}
}

FloatArray*  Vertex :: giveCoordinates() 
{
   if (!coordinates)
      this -> getCoordinates() ;
	return coordinates ;
}

double  Vertex :: giveCoordinate(int i) 
// ************************************
// Returns the i-th coordinate of the receiver.
{
   if (!coordinates)
      this -> getCoordinates() ;
	return coordinates->at(i) ;
}

void Vertex :: setCoordinates(double x, double y)
// **********************************************
{
  coordinates->at(1) = x ; 
  coordinates->at(2) = y ; 
}

bool Vertex :: interactsWith(Element* elem)
// ****************************************
// check if element elem contains the Vertex using the geometry 
// predicates : orientation test
{
	 size_t count = 0;
	 double const EPSILON = 1e-8 ; // this tolerence has strong effect on the tip !!!   

	 double x0 = this->giveCoordinate(1); // coord. of this vertex
	 double y0 = this->giveCoordinate(2);

	 double delta,x1,y1,x2,y2  ;
	
	 size_t numNodes = elem->giveNumberOfNodes();

	 for (size_t i = 1 ; i <= numNodes ; i++)
	 {
		  // coordinates of first node of edge
		  x1 = elem->giveNode(i)->giveCoordinate(1);
		  y1 = elem->giveNode(i)->giveCoordinate(2);

        // coordinates of second node of edge
		  if(i != numNodes)
		  {
			 x2 = elem->giveNode(i+1)->giveCoordinate(1);
			 y2 = elem->giveNode(i+1)->giveCoordinate(2);
		  }
		  else
		  {
			 x2 = elem->giveNode(1)->giveCoordinate(1);
			 y2 = elem->giveNode(1)->giveCoordinate(2);
		  }

		  delta = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0) ;

		  if (delta > EPSILON)
				count += 1 ;
	 }
	 
	 return (count == numNodes); 		
}

bool Vertex::intersects(const Mu::Triangle * t)
// ********************************************
// check the intersection of segment containing the receiver and triangle t.
// Attention : if the crack goes through one of element's nodes => the method 
// Segment::intersects(t) returns FALSE !!!
{
  Mu::Segment * seg = myParent->giveSegmentContain(this);

  bool temp = seg->intersects(t);
  if(temp)
	 return true ;

  // else, need to check more ...
  if(seg->on(t->getBoundingPoint(0)))
	 return true ;
  if(seg->on(t->getBoundingPoint(1)))
	 return true ;
  if(seg->on(t->getBoundingPoint(2)))
	 return true ;

  return false ;
}

bool Vertex::intersects(const Mu::Segment * s)
// ********************************************
// check the intersection of segment containing the receiver and triangle t.
// Attention : if the crack goes through one of element's nodes => the method 
// Segment::intersects(t) returns FALSE !!!
{
  Mu::Segment * seg = myParent->giveSegmentContain(this);
  return seg->intersects(s);
}

std::vector<Mu::Point> Vertex :: intersection(const Mu::Triangle* t)
// *****************************************************************
// two cases may occur : 
//    1. t cuts the crack segment(normal one),then return the intersection point 
//       and the Vertex itself. 
//    2. crack segment is entirely inside the triangle => !!! Corrected : 2005-09-23
{
  std::vector<Mu::Point> ret ;
  std::vector<Mu::Point> temp ;

  Mu::Segment * seg = myParent->giveSegmentContain(this); 
  
  if(this->intersects(t))       // normal case 1
  {
	 temp = seg->intersection(t) ;
	 ret.insert(ret.end(), temp.begin(),temp.end());
	 // check for case crack goes throught node. Corrected 2005-08-09 when deal with 
	 // structured triangle mesh .
	 if(ret.size() == 0)
	 {
		for(size_t i = 0 ; i < 3 ; i++)
		{
		  if(seg->on(t->getBoundingPoint(i)))
			 ret.push_back(*(t->getBoundingPoint(i)));
		}
	 }
  }
  else   // crack segment within element. Occurred when crack grows with small increment.
  {
	 std::vector<Mu::Point> intersect = myParent->intersection(t);
	 ret.insert(ret.end(),intersect.begin(),intersect.end());
  }

  // insert the CrackTip into the result
  Mu::Point point((*coordinates)[0],(*coordinates)[1]);
  ret.push_back(point);

  if(ret.size() != 2) 
  {
	 std::cout <<" Total number of intersection points :" << ret.size() << endl;
	 assert(false);
  }

  return ret ;
}

std::vector<Mu::Point> Vertex :: intersection(const Mu::Segment* s)
// ****************************************************************
{
  std::vector<Mu::Point> ret ;
  
  Mu::Segment * seg = myParent->giveSegmentContain(this);
  // compute the intersection point, if any
  if (this->intersects(s))
  {
	 ret.push_back(seg->intersection(s));

	 // insert the CrackTip into the result
	 double xTip = (*coordinates)[0];
	 double yTip = (*coordinates)[1];
	 Mu::Point point(xTip,yTip);

	 ret.push_back(point);
  }

  return ret;
}

void Vertex :: print()
// *******************
{
  std::cout << "(" << (*coordinates)[0] << "," << (*coordinates)[1] << ")" ;
  std::cout << std::endl ;
}