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

#ifndef _BIMATCRACKASYMPTOTIC_H_
#define _BIMATCRACKASYMPTOTIC_H_

#include "crackasymptotic.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "geometry_base.h"
#include "element.h"
#include "node.h"
#include "gausspnt.h"

//! Asymptotic functions of interfacial cracks
/*!
   This class implements the asymptotic displacement field of homogeneous
   media.\n
 TASKS
   - Compute the value of asymptotic functions.
   - Compute the derivatives of asymptotic functions w.r.t global coord.
 */
class BiMatCrackAsymp : public CrackAsymptotic
{
public:

	BiMatCrackAsymp(Domain*,int,int);   //!< Constructor
	~BiMatCrackAsymp(){;}               //!< Destructor
	/*! 
	 Compute the asymptotic functions
	 */
	double      EvaluateYourSelfAt(GaussPoint *gp);
   double      EvaluateYourSelfAt(Element*,Node*);
	/*! 
	 Compute the gradients of asymptotic functions
	 */
	FloatArray* EvaluateYourGradAt(GaussPoint *gp);
	
private:
	 int order; //!< the order of asymp. function
	 double  computeBiMatConst();//!< Compute the constant Epsilon
} ;

#endif // _BIMATCRACKASYMPTOTIC_H_
