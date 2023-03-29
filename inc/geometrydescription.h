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
 
#ifndef _GEOMETRYDESCRIPTION_H_
#define _GEOMETRYDESCRIPTION_H_


class Element; class GeometryEntity;

// The description type of geometry entities
/*!
 This class is an abstract class, the superclass of all classes that imple-
 ment different types of geometry descriptions such as LevelSetDescription, 
 VectorLevelSetDescription, and StandardDescription.
 */
class GeometryDescription
{
public:
	 GeometryDescription(){}   //!< Constructor
    virtual ~GeometryDescription(){}  //!< Destructor

	 /*!
	  Do the interaction between element and geometry entity
	  Virtual function that requires the derived classes must
	  implement this function.
	  */
	 virtual bool interactsWith(Element*,GeometryEntity*)=0;

} ;


#endif //_GEOMETRYDESCRIPTION_H_
