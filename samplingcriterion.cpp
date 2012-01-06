//
// C++ Implementation: samplingcriterion
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "samplingcriterion.h"

#include <algorithm>
#include <functional>
#include <cmath>  

#ifndef M_PI 
#define M_PI  3.14159265358979323846
#endif


SamplingCriterion::SamplingCriterion()
{
	toAdd = true ;
}


bool SamplingCriterion::add() const
{
	return toAdd ;
}

SamplingCriterion::~SamplingCriterion()
{
}

bool SamplingCriterion::meetsCriterion(const DelaunayTriangle * t)  
{
	return true ;
}

MinimumAngle::MinimumAngle(double a) 
{ 
	angle = a ;
}
	
MinimumAngle::~MinimumAngle() 
{
	
}
	
bool MinimumAngle::meetsCriterion(const DelaunayTriangle * t)  
{
	Point v0 = Point(t->second->x - t->first->x, t->second->y - t->first->y) ;
	Point v1 = Point(t->third->x - t->second->x, t->third->y - t->second->y) ;
	Point v2 = Point(t->first->x - t->third->x, t->first->y - t->third->y) ;
	std::vector<double> ang;
	ang.push_back(
	               asin(
	                     ( v0.x*v1.y - v0.y*v1.x )/( v0.norm()*v1.norm() )
	                   )
	             ) ;
	ang.push_back(
	               asin(
	                     ( v1.x*v2.y - v1.y*v2.x )/( v1.norm()*v2.norm() )
	                   )
	             ) ;
	ang.push_back(
	               asin(
	                     ( v2.x*v0.y - v2.y*v0.x )/( v2.norm()*v0.norm() )
	                   )
	             ) ;
	
	std::sort(ang.begin(), ang.end()) ;

	if(ang[0]+ang[1]+ang[2] > 1.000001 * M_PI)
	{
		std::cout << "what a wierd triangle..." ; t->print() ;
		std::cout << "angle sum is " << ang[0]+ang[1]+ang[2] << std::endl ;
		assert(false) ;
	}

	return ang[0] > angle ;
}

void MinimumAngle::reset() 
{
}

Point * MinimumAngle::suggest(const DelaunayTriangle * t) const
{
	Point v0 = Point(t->second->x - t->first->x, t->second->y - t->first->y) ;
	Point v1 = Point(t->third->x - t->second->x, t->third->y - t->second->y) ;
	Point v2 = Point(t->first->x - t->third->x, t->first->y - t->third->y) ;
	std::map<double, Point> sides;
	sides[v0.norm()] = v0 ;
	sides[v1.norm()] = v1 ;
	sides[v2.norm()] = v2 ;
	
	if(t->getRadius() < sides.rbegin()->first)
		return new Point(t->getCircumCenter()) ;
	
	if(sides.rbegin()->second == v0)
		return new Point(t->second->x*0.5+ t->first->x*0.5, t->second->y*0.5+ t->first->y*0.5 ) ;
	if(sides.rbegin()->second == v1)
		return new Point(t->second->x*0.5+ t->third->x*0.5, t->second->y*0.5+ t->third->y*0.5 ) ;
	if(sides.rbegin()->second == v2)
		return new Point(t->first->x*0.5+ t->third->x*0.5, t->third->y*0.5+ t->first->y*0.5 ) ;
	
	assert(false) ;
	return NULL ;
}

MaximumLength::MaximumLength(double l)
{
	length = l ;
}
	
MaximumLength::~MaximumLength()
{
}
	
bool MaximumLength::meetsCriterion(const DelaunayTriangle * t) 
{
	Point v0 = Point(t->second->x - t->first->x, t->second->y - t->first->y) ;
	Point v1 = Point(t->third->x - t->second->x, t->third->y - t->second->y) ;
	Point v2 = Point(t->first->x - t->third->x, t->first->y - t->third->y) ;
	
	std::vector<double> len;
	len.push_back(v0.norm()) ;
	len.push_back(v1.norm()) ;
	len.push_back(v2.norm()) ;
	
	std::sort(len.begin(), len.end()) ;
	
	return len[0] < length ;
}

void MaximumLength::reset() 
{
}

Counter::Counter(size_t nit)
{
	this->c = 0 ;
	this->s = nit ;
	toAdd = false ;
}
	
Counter::~Counter()
{
}
	
bool Counter::meetsCriterion(const DelaunayTriangle * t) 
{
	c++ ;
	
	return c < s ;
}

void Counter::reset()
{
	c = 0 ;
}
