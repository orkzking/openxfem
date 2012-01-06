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

#include "discontinuousfunction.h"
#include "element.h"
#include "node.h"
#include "geometryentity.h"
#include "feinterpol.h"

double DiscontinuousFunction::EvaluateYourSelfAt(GaussPoint* gp)
// *************************************************************
// Compute the value of discontinuous function H(x) at coord.
{
  Element *e = gp->giveMyElement();
  Mu::Point *globalCoord = 
	 e->giveFEInterpolation()->local2Global(domain,e->giveNodeArray(),gp->giveCoordinates());
  int check = this->activeEnrItem->giveMyGeo()->givePositionComparedTo(globalCoord);
  delete globalCoord ;
  if (check == 1)
	 return 1.0;
  else if (check == -1)
	 return -1.0;
  else
  {
	 return 0.0 ;    // corrected 2005-08-11. Need to discuss more on this !!!
  }
}
/*
double DiscontinuousFunction :: EvaluateYourSelfAt(Element* e ,Node* n)
// ********************************************************************
// for nodes lie on the crack, the H(x) at them are multi-valued(+1,-1)
// then need to check these nodes belong to element whose center above or below 
// the crack to choose the correct H(x) value.
{
Point* p = n->makePoint() ; 	
int check = this->activeEnrItem->giveMyGeo()->givePositionComparedTo(p);
delete p ; 
if (check == 1)         // above the crack
return 1.0;
if (check == -1)        // below the crack 
return -1.0;

// node lies on crack

Point *center = e->giveMyCenter(); 
int temp = this->activeEnrItem->giveMyGeo()->givePositionComparedTo(center);
delete center ;

if (temp == 1)              // center is above the crack
return 1.0;
else if (temp == -1)        // center is below the crack
return -1.0;
else
return 0.0 ;              // center is also on the crack, 2005-08-09

} */


double DiscontinuousFunction :: EvaluateYourSelfAt(Element* e ,Node* n)
// ********************************************************************
// Nodes on the crack then H = zero !!!
// Corrected 2005-08-11
{
  Point* p = n->makePoint() ; 

  int check = this->activeEnrItem->giveMyGeo()->givePositionComparedTo(p);

  if (check == 1)         // above the crack
	 return 1.0;
  if (check == -1)        // below the crack 
	 return -1.0;
  else
  {
	 //std::cout << " Split-enriched node belong to the crack !!! " << std::endl ;
	 return 0.0 ;
  }
  delete p ; 
} 

FloatArray* DiscontinuousFunction ::EvaluateYourGradAt(GaussPoint* gp)
// ********************************************************************
// Compute the derivatives of discontinuous function H(x)
// This derivative is zero 
{
  return new FloatArray(2);
}
