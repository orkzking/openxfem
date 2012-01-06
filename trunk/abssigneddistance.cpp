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

#include "abssigneddistance.h"
#include "element.h"
#include "node.h"
#include "gausspnt.h"
#include "feinterpol.h"
#include "mathfem.h"
#include "geometryentity.h"


double AbsSignedDistance::EvaluateYourSelfAt(GaussPoint* gp)
// *********************************************************
// Compute the ramp function used to model the material interface
// Definition: abs(signed distance to the interface)
{
  Element *e = gp->giveMyElement();
  Mu::Point* p = gp->giveCoordinates();
  FloatArray* N = e->giveFEInterpolation()->evalN(p);
  double ls; double ret = 0.0 ; 
  for(size_t i = 0 ; i < e->giveNumberOfNodes() ; i++)
  {
	 ls = e->giveNode(i+1)->giveLevelSet(activeEnrItem) ;
	 if(ls == 0) // signed distance not yet computed at this node
	 {
		ls = activeEnrItem->giveMyGeo()->computeSignedDistanceOfPoint(e->giveNode(i+1)->makePoint());
		//cout<<ls<<endl;
	 }
	 ret += (*N)[i]*ls ;
  }

  return abs(ret) ;
}


double AbsSignedDistance :: EvaluateYourSelfAt(Element* e ,Node* n)
// ****************************************************************
{
  double ret = n->giveLevelSet(activeEnrItem);
  if(ret == 0) // signed distance not yet computed at this node
  {
	 ret = activeEnrItem->giveMyGeo()->computeSignedDistanceOfPoint(n->makePoint());
  }
  return abs(ret) ;
} 

FloatArray* AbsSignedDistance ::EvaluateYourGradAt(GaussPoint* gp)
// ***************************************************************
// Compute the derivatives of discontinuous function H(x)
// This derivative is zero 
{
  Element *e = gp->giveMyElement(); 
  Mu::Point* p = gp->giveCoordinates();
  FloatArray*     N = e->giveFEInterpolation()->evalN(p);
  FloatMatrix* dNdx = e->giveFEInterpolation()->evaldNdx(domain,e->giveNodeArray(),p);
  double ls;
  double phiPt = 0.0 ; double dFdx = 0.0 ; double dFdy = 0.0 ;
  



  for(size_t i = 0 ; i < e->giveNumberOfNodes() ; i++)
  {
	 ls = e->giveNode(i+1)->giveLevelSet(activeEnrItem) ;
	 if(ls == 0) // signed distance not yet computed at this node
		ls = activeEnrItem->giveMyGeo()->computeSignedDistanceOfPoint(e->giveNode(i+1)->makePoint());

	 phiPt += (*N)[i]*ls ;
	 dFdx  += dNdx->at(i+1,1)*ls ; // dFdx
	 dFdy  += dNdx->at(i+1,2)*ls ; // dFdy
  }

  FloatArray* answer = new FloatArray(2);
  (*answer)[0] = dFdx * sgn(phiPt) ; 
  (*answer)[1] = dFdy * sgn(phiPt) ;

  return answer;
}
