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

 
#ifndef _JUNCTIONFUNCTION_H_
#define _JUNCTIONFUNCTION_H_

#include "enrichmentfunction.h"
#include "flotarry.h"
#include "geometry_base.h"
#include "element.h"
#include "node.h"


class JunctionFunction : public EnrichmentFunction
{
public:
   /*!Constructor
	 */
   JunctionFunction(Domain* aDomain,int n) : EnrichmentFunction(aDomain,n){}
	virtual ~JunctionFunction(){}  //!< Destructor
	/*!
	 Compute the step function at point with coordinates
	 given by FloatArray* coord
	 Most of the times, coord is a Gauss Point
	 */
	double         EvaluateYourSelfAt(GaussPoint* gp);
	double         EvaluateYourSelfAt(Element*,Node*);
	/*! 
	 Compute the gradients of discontinuous function
	 */
	FloatArray*    EvaluateYourGradAt(GaussPoint* gp);
} ;


#endif //_JUNCTIONFUNCTION_H_
