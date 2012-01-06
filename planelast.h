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

#ifndef _PLANELAST_H_
#define _PLANELAST_H_

#include "flotmtrx.h"
#include "material.h"

//! Compuatation of constituive matrix of plane elasticity
/*!
 This class implementes the computation ofconstituive matrix for plane stress problem.
 It's used in class Quad_U or Tri_U to build the constituive matrix 
 */

class PlaneElasticity
{
public:
  virtual FloatMatrix* computeConstitutiveMatrix(Material* mat) = 0;
};

class PlanStrain : public PlaneElasticity
{
public:
  FloatMatrix* computeConstitutiveMatrix(Material* mat);
};

class PlanStress : public PlaneElasticity
{
public:
  FloatMatrix* computeConstitutiveMatrix(Material* mat);
};


#endif // _PLANELAST_H_
