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
 
#ifndef _HOLE_H_
#define _HOLE_H_

#include "enrichmentitem.h"

//! Holes
/*!
 This class implements a hole.
 */
class Hole:public EnrichmentItem
{

public:
	  Hole(int,Domain*);   //!< constructor
	 ~Hole(){}             //!< destructor
	  void  instanciateYourself();

	 /*!
	  Time step termination
	  */
	 void     printOutputAt (TimeStep*, FILE*, FILE*){}

	 void     updateEnrichment(){;}

	 void     printYourSelf();
     
} ;


#endif //_HOLE_H_








