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
#ifndef _FEINTERPOL2D_H_
#define _FEINTERPOL2D_H_

#include "feinterpol.h"

//! 2D FE interpolation
/*!
 This class is an abstract class, the superclass of all classes that imple-
 ment the 2D finite element interpolations such as FE2dquadlin, FE2dtrilin...
 */
class FEInterpolation2d:public FEInterpolation
{

public:

  FEInterpolation2d(int o):FEInterpolation(o){}

  virtual FloatArray* giveDerivativeKsi (double=0,double=0) = 0 ;
  virtual FloatArray* giveDerivativeEta (double=0,double=0) = 0 ;
	
};


#endif //_FEINTERPOL2D_H_
