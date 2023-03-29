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

#include "voidenrichfunction.h"
#include "geometry_base.h"
#include "circle.h"


double VoidFunction::EvaluateYourSelfAt(Mu::Point* p)
// ***************************************************
// Compute the value of void function V(x)
{
  Cercle * myGeo = static_cast<Cercle *>(activeEnrItem->giveMyGeo());

  double xc = myGeo->giveCenter()->giveCoordinate(1) ;
  double yc = myGeo->giveCenter()->giveCoordinate(2) ;

  double radi = myGeo->giveRadius();

  Mu::Circle* c = new Mu::Circle(radi,new Mu::Point(xc,yc));

  if (c->in(p))
	 return 0.0 ;
  else
	 return 1.0;
}

FloatArray* VoidFunction ::EvaluateYourGradAt(Mu::Point* coord)
// ************************************************************
// Compute the derivatives of void function V(x)
// This derivative is always zero 
{
	return new FloatArray(2);
}

