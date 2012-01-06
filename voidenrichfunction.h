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

#ifndef _VOIDFUNCTION_H_
#define _VOIDFUNCTION_H_

#include "enrichmentfunction.h"
#include "flotarry.h"
#include "geometry_base.h"


//! Enrichment function for voids : function V(x)
/*!
 This class implements a discontinuous function V(x), see Daux et al. 2000
 used to enrich the nodes whose support is cut by holes
 */
class VoidFunction:public EnrichmentFunction{
public:
   /*!Constructor
	 */
   VoidFunction(Domain* aDomain,int n) : EnrichmentFunction(aDomain,n){}
	virtual ~VoidFunction(){}  //!< Destructor
	/*!
	 Compute the void function at point with coordinates
	 given by @param FloatArray* p
	 @return 1 if p is outside the hole and 0 otherwise
	 */
	double         EvaluateYourSelfAt(Mu::Point* p);
	/*! 
	 Compute the gradients of void function, this value is always zero !
	 */
	FloatArray*    EvaluateYourGradAt(Mu::Point*);
} ;


#endif //_DISCONTINUOUSFUNCTION_H_
