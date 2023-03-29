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


#include "planelast.h"

FloatMatrix* PlanStrain :: computeConstitutiveMatrix(Material* mat)
// ****************************************************************
// compute the constituive matrix of plane strain  
{
   double e     = mat -> give('E') ;
   double nu    = mat -> give('n') ;
   double ee    = e / ((1.+nu) * (1.-nu-nu)) ;
   double shear = e / (2.+nu+nu) ;

   FloatMatrix* constitutiveMatrix = new FloatMatrix(4,4) ;

   constitutiveMatrix->at(1,1) = (1.-nu) * ee ;
   constitutiveMatrix->at(1,2) =     nu  * ee ;
   constitutiveMatrix->at(2,1) =     nu  * ee ;
   constitutiveMatrix->at(2,2) = (1.-nu) * ee ;
   constitutiveMatrix->at(3,3) =  shear ;
   constitutiveMatrix->at(1,4) =     nu  * ee ;
   constitutiveMatrix->at(2,4) =     nu  * ee ;
   constitutiveMatrix->at(4,1) =     nu  * ee ;
   constitutiveMatrix->at(4,2) =     nu  * ee ;
   constitutiveMatrix->at(4,4) = (1.-nu) * ee ;

   return constitutiveMatrix ;
}

FloatMatrix* PlanStress :: computeConstitutiveMatrix(Material* mat)
// ****************************************************************
// compute the constituive matrix of plane stress  
{
   double e     = mat -> give('E') ;
   double nu    = mat -> give('n') ;
   double shear = e / (1. - nu*nu) ;

   FloatMatrix*  constitutiveMatrix = new FloatMatrix(4,4) ;

   constitutiveMatrix->at(1,1) = shear ;
   constitutiveMatrix->at(1,2) = nu  * shear ;
   constitutiveMatrix->at(2,1) = nu  * shear ;
   constitutiveMatrix->at(2,2) = shear ;
   constitutiveMatrix->at(3,3) = shear * (1-nu)/2.  ;
   
   return constitutiveMatrix ;
}
