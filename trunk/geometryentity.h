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
#ifndef	_GEOMETRYENTITY_H_
#define	_GEOMETRYENTITY_H_

#include	"femcmpnn.h"
#include	"geometrydescription.h"
#include	"geometry_2D.h"
#include <stdio.h>
#include <iostream>

class Element;

//! A	geometry	entity is geometry of an enrichment	item.
/*!
This	class	is	an	abstract	class, the superclass of all classes that	imple-
ment	different types of geometry entities: Vertex, PiecewiseLinear,	PicewiseParabolic,
Circle and	Ellipse.\n
DESCRIPTION \n
Enrichment items will use this geometry to find	enriched	nodes.
*/


class	GeometryEntity:public FEMComponent
{

public:

  //!< Constructor
  GeometryEntity(int n ,Domain* aDomain):FEMComponent(n,aDomain){
	 geoDescription = 0 ;
  }
  virtual ~GeometryEntity(){}//!< Destructor

  /** @name	Definition of an element
   *  Functions	to	define an element
   */
  //@{
  GeometryEntity*				  typed () ;		  
  GeometryEntity*				  ofType	(char*) ;
  char*						 giveClassName	(char* s)
  { return	strcpy(s,"GeometryEntity")	;}
  int						 giveNumber	()
  { return	FEMComponent::giveNumber()	;}
  //@}
  /*!
   Do	the interaction between	the geometry entity and	the finite elements
   This function depends on the Geometry Description of this geometry
   Actual code	is	implemented	in	classes of GeometryDescription
  */
  virtual bool interactsWith(Element*);
  /*!
   Compute	the signed distance of a given point to the receiver.
   Virtual	function	 
  */
  virtual double	  computeSignedDistanceOfPoint(Mu::Point*){ return	NULL ;}
  /*!
   Read from the input file the type	of	description	of	the receiver
   Possible	geometry	descriptions are standard,	level	set and vector	
   level	set.
	@return : StandardDescription (geoDescription=1)
	          LevelSetDescription (geoDescription=2)
				 VectorLevelSetDescription (geoDescription=3)
  */
   GeometryDescription*	giveMyGeoDescription();
  /*!
   Returns	the position of the point compared to the	receiver
   1	if	point	is	inside (circle) above (line) and	0 otherwise
   Virtual	function, so derived	class	have to implemente this	method  
   */
  virtual int	 givePositionComparedTo(Mu::Point*)	{ return	NULL ;}

  /*!
   Check the intersection	between the	geometry	entity and	the triangle
   */
  virtual bool intersects(const Mu::Triangle* t){return	NULL;}
  virtual bool intersects(const Mu::Segment* s ){return	NULL;}
  /*!
   Computes the intersection	points between	the receiver with	the triangle
   */
  virtual std::vector<Mu::Point> intersection(const Mu::Triangle* t)
  {
	 return std::vector<Mu::Point>();
  }
  virtual std::vector<Mu::Point> intersection(const Mu::Segment* s)
  {
	 return std::vector<Mu::Point>();
  }
  /*!
  Computes the intersection points between the receiver with a given element
  Used to partition split elements into sub-cells
   */
  virtual std::vector<Mu::Point> intersection(Element* e);
  /*!
   Export geometry to Matlab file for plotting.
   */
  virtual void  exportToMatlab(std::string&){;}

protected:

  size_t	   geoDescription ;//!< Description of Geometry entity
} ;


#endif // _GEOMETRYENTITY_H_
