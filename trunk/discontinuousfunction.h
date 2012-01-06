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

 
#ifndef _DISCONTINUOUSFUNCTION_H_
#define _DISCONTINUOUSFUNCTION_H_

#include "enrichmentfunction.h"
#include "flotarry.h"
#include "element.h"
#include "node.h"
#include "geometry_base.h"
#include "gausspnt.h"

//! Discontinuous function : Heaviside function H(x)
/*!
 This class implements a discontinuous function H(x)
 used to enrich the nodes whose support is cut by the discontinuities
 */
class DiscontinuousFunction:public EnrichmentFunction
{
public:
   /*!Constructor
	 */
   DiscontinuousFunction(Domain* aDomain,int n) : EnrichmentFunction(aDomain,n){}
	virtual ~DiscontinuousFunction(){}  //!< Destructor
	/*!
	 Compute the step function at \c GaussPoint gp 
	 */
	double         EvaluateYourSelfAt(GaussPoint* gp);
	double         EvaluateYourSelfAt(Element*,Node*);
	/*! 
	 Compute the gradients of discontinuous function at \c GaussPoint gp 
	 */
	FloatArray*    EvaluateYourGradAt(GaussPoint* gp);
} ;


#endif //_DISCONTINUOUSFUNCTION_H_
