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

#include "circle.h"

#include "vertex.h"
#include "geometry_2D.h"
#include <math.h>
#include <iostream>
#include <string>


double Cercle :: giveRadius()
{
  if(radius == 0)
	 radius = this->read("radius") ;
  return radius ;
}

Vertex* Cercle :: giveCenter()
{
  if (center == 0)
	 center = this ->readInteger("center") ;

  Vertex *myCenter = static_cast<Vertex*>(domain->giveGeoEntity(center)) ;

  return myCenter;
}

double Cercle ::  computeSignedDistanceOfPoint(Mu::Point * coord)
// compute the difference between the distance from center xc to 
// point xi and the radius r.
// return ||xi-xc|| - r
{
  double xc,yc,distance;

  xc = this->giveCenter()->giveCoordinate(1) ;
  yc = this->giveCenter()->giveCoordinate(2) ;

  distance = sqrt((coord->x-xc)*(coord->x-xc)+(coord->y-yc)*(coord->y-yc));

  return (distance - (this->giveRadius()));

}

int Cercle:: givePositionComparedTo(Mu::Point* coord)
{
  double d = this->computeSignedDistanceOfPoint(coord);
  return ( d>0 ? 1 : 0 );

}

bool Cercle::intersects(const Mu::Triangle* t)
{
  Segment s0((*(*t->getBoundingPoints())[0]), (*(*t->getBoundingPoints())[1])) ;
  Segment s1((*(*t->getBoundingPoints())[1]), (*(*t->getBoundingPoints())[2])) ;
  Segment s2((*(*t->getBoundingPoints())[2]), (*(*t->getBoundingPoints())[0])) ;

  double xc = this->giveCenter()->giveCoordinate(1) ;
  double yc = this->giveCenter()->giveCoordinate(2) ;
  Mu::Circle* c = new Mu::Circle(this->giveRadius(),new Mu::Point(xc,yc));

  return s0.intersects(c) || s1.intersects(c) || s2.intersects(c) ;
}

bool Cercle::intersects(const Mu::Segment* s)
{
  double xc = this->giveCenter()->giveCoordinate(1) ;
  double yc = this->giveCenter()->giveCoordinate(2) ;
  Mu::Circle* c = new Mu::Circle(this->giveRadius(),xc,yc);
  bool ret = s->intersects(c);
  delete c ;
  return ret ;
}

std::vector<Point> Cercle::intersection(const Triangle *t) 
{
  double xc = this->giveCenter()->giveCoordinate(1) ;
  double yc = this->giveCenter()->giveCoordinate(2) ;
  Mu::Circle* c = new Mu::Circle(this->giveRadius(),xc,yc);

  std::vector<Point> ret ;
  Segment s0((*(*t->getBoundingPoints())[0]), (*(*t->getBoundingPoints())[1])) ;

  if(s0.intersects(c))
	 ret.insert(ret.end(),s0.intersection(c).begin(),s0.intersection(c).end()) ;

  Segment s1((*(*t->getBoundingPoints())[1]), (*(*t->getBoundingPoints())[2])) ;

  if(s1.intersects(c))
	 ret.insert(ret.end(),s0.intersection(c).begin(),s0.intersection(c).end()) ;

  Segment s2((*(*t->getBoundingPoints())[2]), (*(*t->getBoundingPoints())[0])) ;

  if(s2.intersects(c))
	 ret.insert(ret.end(),s0.intersection(c).begin(),s0.intersection(c).end()) ;

  delete c ;

  return ret ;
}

std::vector<Point> Cercle::intersection(const Segment *s) 
{
  double xc = this->giveCenter()->giveCoordinate(1) ;
  double yc = this->giveCenter()->giveCoordinate(2) ;
  Mu::Circle* c = new Mu::Circle(this->giveRadius(),xc,yc);

  std::vector<Point> ret ;

  if(s->intersects(c))
	 ret = s->intersection(c) ;
  delete c ;
  return ret ;
}

void Cercle::exportToMatlab(std::string& theString)
{
  char value[30];	
  std::string space(" ");
  std::string newline("\n");

  double x = this->giveCenter()->giveCoordinate(1) ;  
  double y = this->giveCenter()->giveCoordinate(2) ; 
  double r = this->giveRadius() ;

#ifdef WIN32
  _gcvt(x,6,value);     // transforms the double param1 in char param3 with 3 digits
  theString += value;
  theString += space;
  _gcvt(y,6,value);
  theString += value;
  theString += space;
  _gcvt(r,6,value);
  theString += value;
  theString += space;
#else
  qgcvt(x,6,value);     // transforms the double param1 in char param3 with 3 digits
  theString += value;
  theString += space;
  qgcvt(y,6,value);
  theString += value;
  theString += space;
  qgcvt(r,6,value);
  theString += value;
  theString += space;
#endif
  theString += newline;
}

