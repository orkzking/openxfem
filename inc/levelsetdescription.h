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
   
	Any feedback is welcome. Emails : nvphu80@yahoo.com, ...

*************************************************************************************/
 
#ifndef _LEVELSETDESCRIPTION_H_
#define _LEVELSETDESCRIPTION_H_

#include "geometrydescription.h"
#include "domain.h"
class Element; class GeometryEntity;


//! Level set description of geometry entities  
/*!
 The geometry of the discontinuities such as PiecewiseLinear, Circle, Vertex
 is represented by the level sets defined at nodes.
 see Storlaska 2001
 */
class LevelSetDescription:public GeometryDescription
{
public:
  LevelSetDescription(){}          //!< Constructor
	~LevelSetDescription(){}        //!< Destructor
  /*!
	Check if element interacts with GeometryEntity ot not
	using the level set at nodes of element
	 */
	bool interactsWith(Element*,GeometryEntity*);
} ;


#endif // _LEVELSETDESCRIPTION_H_
