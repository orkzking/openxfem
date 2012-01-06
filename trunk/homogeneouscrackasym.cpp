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

#include "homogeneouscrackasymp.h"

#include "crackasymptotic.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "cracktip.h"
#include "assert.h"
#include <stdlib.h>
#include <math.h>
#include "feinterpol.h"

#ifndef M_PI 
#define M_PI  3.14159265358979323846
#endif


HomoElastCrackAsymp :: HomoElastCrackAsymp(Domain* aDomain,int n,int m)
: CrackAsymptotic (aDomain,n),order(m)
{
}
/*
// computation of asymptotic functions and their derivatives.  
// sqrt(r)*sin(theta/2),sqrt(r)*cos(theta/2),sqrt(r)*sin(theta/2)*sin(theta),sqrt(r)*cos(theta/2)*sin(theta)

double HomoElastCrackAsymp::EvaluateYourSelfAt(Mu::Point* coord)
// **************************************************************
// compute the i-th asymptotic function 
{
FloatArray * localCoord = 
dynamic_cast<CrackTip*>(this->activeEnrItem)->computeLocalCoordOf(coord) ;

double xloc = (*localCoord)[0];
double yloc = (*localCoord)[1];

double r = sqrt( xloc * xloc + yloc * yloc ); 
double theta = atan2(yloc,xloc) ;  

// for sure, -Pi <= theta <= Pi
assert(theta <= M_PI && theta >= -M_PI);

double answer ;
switch(this->order)
{
case 1:
answer = sqrt(r) * sin(theta/2.) ; 
break ;
case 2:
answer = sqrt(r) * cos(theta/2.);
break ;
case 3:
answer = sqrt(r) * sin(theta/2.) * sin(theta);
break;
case 4:
answer = sqrt(r) * cos(theta/2.) * sin(theta);
break ;
}

return answer ;
}

double HomoElastCrackAsymp::EvaluateYourSelfAt(Element* e ,Node* n)
// ****************************************************************
// compute the i-th asymptotic function at Nodes.
// the first asymptotic function is multi-valued for nodes lie on the crack.
{
Mu::Point* p = n->makePoint() ; 	
FloatArray * localCoord = 
dynamic_cast<CrackTip*>(this->activeEnrItem)->computeLocalCoordOf(p) ;

double xloc = (*localCoord)[0];
double yloc = (*localCoord)[1];

double r = sqrt( xloc * xloc + yloc * yloc ); 
double theta = atan2(yloc,xloc) ;  

// for sure, -Pi <= theta <= Pi
assert(theta <= M_PI && theta >= -M_PI);
delete p ;

double answer ;
switch(this->order)
{
case 1:
if(theta != M_PI)                     // node does not lie on the crack
answer = sqrt(r) * sin(theta/2.) ; 
else                                  // node lies on the crack !!!
{
Mu::Point *center = e->giveMyCenter(); 
Mu::Segment *tipSeg = dynamic_cast<CrackTip*>(this->activeEnrItem)->giveTipSegment();
// use orientation test to check the center's location w.r.t the tip segment 
double x1 = tipSeg->first()->x ;
double y1 = tipSeg->first()->y ;
double x2 = tipSeg->second()->x ;
double y2 = tipSeg->second()->y ;

double delta = (x1-center->x)*(y2-center->y) - (x2-center->x)*(y1-center->y) ;
delete center ;

if (delta > 0.000001)         // center is above crack
answer = sqrt(r) ; 
else                          // center is below crack
answer = -sqrt(r) ; 
}
break ;
case 2:
answer = sqrt(r) * cos(theta/2.);
break ;
case 3:
answer = sqrt(r) * sin(theta/2.) * sin(theta);
break;
case 4:
answer = sqrt(r) * cos(theta/2.) * sin(theta);
break ;
}

return answer ;
}

FloatArray* HomoElastCrackAsymp::EvaluateYourGradAt(Mu::Point* p)
// **************************************************************
// compute the grad of i-th asymp. function
// First, derivatives w.r.t local crack tip are computed. Then, transformed
// to the global coordinate system.
// result = [dPhidx dPhidy]
{
// local crack tip coordinates of Point p
FloatArray *localCoord = 
static_cast<CrackTip*>(activeEnrItem)->computeLocalCoordOf(p) ;

// compute the inclination angle of the local crack tip coord. sys w.r.t the x axe
Mu::Segment *seg = static_cast<CrackTip*>(activeEnrItem)->giveTipSegment();
double x1 = seg->first()->x ;
double y1 = seg->first()->y ;
double x2 = seg->second()->x ;
double y2 = seg->second()->y ;

double alpha = atan2(y2-y1,x2-x1) ; // inclination angle

double xloc = (*localCoord)[0];
double yloc = (*localCoord)[1];

double r = sqrt( xloc * xloc + yloc * yloc );
double theta = atan2(yloc,xloc) ;

// for sure, -Pi <= theta <= Pi
assert(theta <= M_PI && theta >= -M_PI);

// !!! How about the case r = 0 ???
double fac  = 1./(2.* sqrt(r)) ;

double dPhidx1,dPhidx2 ;
FloatArray* answer = new FloatArray(2); 
switch(this->order)
{
case 1:
dPhidx1 = -fac * sin(theta/2.);
dPhidx2 =  fac * cos(theta/2.);
(*answer)[0] = dPhidx1 * cos(alpha) - dPhidx2 * sin(alpha)  ;
(*answer)[1] = dPhidx1 * sin(alpha) + dPhidx2 * cos(alpha)  ;
break;
case 2:
dPhidx1 = fac * cos(theta/2.); 
dPhidx2 = fac * sin(theta/2.); 
(*answer)[0] = dPhidx1 * cos(alpha) - dPhidx2 * sin(alpha)  ;
(*answer)[1] = dPhidx1 * sin(alpha) + dPhidx2 * cos(alpha)  ;
break;
case 3:
dPhidx1 = -fac * sin(3.* theta/2.) * sin(theta);
dPhidx2 =  fac * (sin(theta/2.) + sin(3.* theta/2.)*cos(theta));
(*answer)[0] = dPhidx1 * cos(alpha) - dPhidx2 * sin(alpha)  ;
(*answer)[1] = dPhidx1 * sin(alpha) + dPhidx2 * cos(alpha)  ;
break;
case 4:
dPhidx1 = -fac * cos(3.* theta/2.) * sin(theta);
dPhidx2 =  fac * (cos(theta/2.) + cos(3.* theta/2.)*cos(theta));
(*answer)[0] = dPhidx1 * cos(alpha) - dPhidx2 * sin(alpha)  ;
(*answer)[1] = dPhidx1 * sin(alpha) + dPhidx2 * cos(alpha)  ;
break;
}

return answer;
}
*/

// Alternative asymptotic functions implementation !!!
// sqrt(r)*sin(theta/2),sqrt(r)*cos(theta/2),
// sqrt(r)*sin(theta/2)*cos(theta),sqrt(r)*cos(theta/2)*cos(theta)


double HomoElastCrackAsymp::EvaluateYourSelfAt(GaussPoint* gp)
// ***********************************************************
// compute the i-th asymptotic function 
{
  Element *e = gp->giveMyElement();
  Mu::Point *globalCoord = 
	 e->giveFEInterpolation()->local2Global(domain,e->giveNodeArray(),gp->giveCoordinates());
  std::valarray<double> *localCoord = 
	 dynamic_cast<CrackTip*>(this->activeEnrItem)->computePolarCoordOf(globalCoord);
  
  double r     = (*localCoord)[0];
  double theta = (*localCoord)[1];

  // for sure, -Pi <= theta <= Pi
  assert(theta <= M_PI && theta >= -M_PI);
  delete localCoord ; // Purify 27-10-2005

  double answer ;
  switch(this->order)
  {
  case 1:
	 answer = sqrt(r) * sin(theta/2.) ; 
	 break ;
  case 2:
	 answer = sqrt(r) * cos(theta/2.);
	 break ;
  case 3:
	 answer = sqrt(r) * sin(theta/2.) * cos(theta); 
	 break;
  case 4:
	 answer = sqrt(r) * cos(theta/2.) * cos(theta); 
	 break ;
  }
  delete globalCoord ;
  return answer ;
}
/*
double HomoElastCrackAsymp::EvaluateYourSelfAt(Element* e ,Node* n)
// ****************************************************************
// compute the i-th asymptotic function at Nodes.
// the first asymptotic function is multi-valued for nodes lie on the crack.
{
Mu::Point* p = n->makePoint() ; 	
FloatArray * localCoord = 
dynamic_cast<CrackTip*>(this->activeEnrItem)->computeLocalCoordOf(p) ;

double xloc = (*localCoord)[0];
double yloc = (*localCoord)[1];

double r = sqrt( xloc * xloc + yloc * yloc ); 
double theta = atan2(yloc,xloc) ;  

// for sure, -Pi <= theta <= Pi
assert(theta <= M_PI && theta >= -M_PI);
delete p ;

double answer ;
switch(this->order)
{
case 1:
if(theta != M_PI)                     // node does not lie on the crack
answer = sqrt(r) * sin(theta/2.) ; 
else                                  // node lies on the crack !!!
{
Mu::Point *center = e->giveMyCenter(); 
Mu::Segment *tipSeg = dynamic_cast<CrackTip*>(this->activeEnrItem)->giveTipSegment();
// use orientation test to check the center's location w.r.t the tip segment 
double x1 = tipSeg->first()->x ;
double y1 = tipSeg->first()->y ;
double x2 = tipSeg->second()->x ;
double y2 = tipSeg->second()->y ;

double delta = (x1-center->x)*(y2-center->y) - (x2-center->x)*(y1-center->y) ;
delete center ;

if (delta > 0.000001)         // center is above crack
answer = sqrt(r) ; 
else                          // center is below crack
answer = -sqrt(r) ; 
}
break ;
case 2:
answer = sqrt(r) * cos(theta/2.);
break ;
case 3:
answer = sqrt(r) * sin(theta/2.) * cos(theta);
break;
case 4:
answer = sqrt(r) * cos(theta/2.) * cos(theta);
break ;
}

return answer ;
}*/

double HomoElastCrackAsymp::EvaluateYourSelfAt(Element* e ,Node* n)
// ****************************************************************
// compute the i-th asymptotic function at Nodes.
// corrected 2005-08-11. No distinguish for case theta = PI .
{
  Mu::Point* p = n->makePoint() ; 	
  std::valarray<double> *localCoord = 
	 dynamic_cast<CrackTip*>(this->activeEnrItem)->computePolarCoordOf(p) ;

  double r     = (*localCoord)[0];
  double theta = (*localCoord)[1];
  // for sure, -Pi <= theta <= Pi
  assert(theta <= M_PI && theta >= -M_PI);
  delete p ; delete localCoord ; // Purify 27-10-2005

  double answer ;
  switch(this->order)
  {
  case 1:
	 answer = sqrt(r) * sin(theta/2.) ; 
	 // debug only
	 //if(theta == M_PI)
	 // std::cout << " Tip-enriched node belong to the crack !!! " << std::endl ;
	 break ;
  case 2:
	 answer = sqrt(r) * cos(theta/2.);
	 break ;
  case 3:
	 answer = sqrt(r) * sin(theta/2.) * cos(theta);
	 break;
  case 4:
	 answer = sqrt(r) * cos(theta/2.) * cos(theta);
	 break ;
  }

  return answer ;
}

FloatArray* HomoElastCrackAsymp::EvaluateYourGradAt(GaussPoint* gp)
// ****************************************************************
// compute the grad of i-th asymp. function
// 1. derivatives w.r.t local crack tip coordinate system are computed. 
// 2. transformed to the global coordinate system.
// result = [dPhidx dPhidy]
{
  // inclination angle of the local crack tip coord. sys w.r.t the x axe
  double alpha = static_cast<CrackTip*>(activeEnrItem)->giveInclinationAngle();

  Element *e = gp->giveMyElement();
  Mu::Point *globalCoord = 
	 e->giveFEInterpolation()->local2Global(domain,e->giveNodeArray(),gp->giveCoordinates());
  std::valarray<double> *localCoord = 
	 dynamic_cast<CrackTip*>(this->activeEnrItem)->computePolarCoordOf(globalCoord) ;

  double r     = (*localCoord)[0];
  double theta = (*localCoord)[1];

  // for sure, -Pi <= theta <= Pi
  assert(theta <= M_PI && theta >= -M_PI);

  double r2 = sqrt(r) ; // r = 0 ???
  double fac  = 0.5/r2 ;

  double dPhidx1,dPhidx2 ;
  FloatArray* answer = new FloatArray(2); 
  switch(this->order)
  {
  case 1:
	 dPhidx1 = -fac * sin(theta/2.);
	 dPhidx2 =  fac * cos(theta/2.);
	 (*answer)[0] = dPhidx1 * cos(alpha) - dPhidx2 * sin(alpha)  ;
	 (*answer)[1] = dPhidx1 * sin(alpha) + dPhidx2 * cos(alpha)  ;
	 break;
  case 2:
	 dPhidx1 = fac * cos(theta/2.); 
	 dPhidx2 = fac * sin(theta/2.); 
	 (*answer)[0] = dPhidx1 * cos(alpha) - dPhidx2 * sin(alpha)  ;
	 (*answer)[1] = dPhidx1 * sin(alpha) + dPhidx2 * cos(alpha)  ;
	 break;
  case 3:
	 dPhidx1 = fac * sin(theta/2.) * (2*sin(theta)*sin(theta)-cos(theta)) ;        
	 dPhidx2 = fac * cos(theta)* cos(3*theta/2.);
	 (*answer)[0] = dPhidx1 * cos(alpha) - dPhidx2 * sin(alpha)  ;
	 (*answer)[1] = dPhidx1 * sin(alpha) + dPhidx2 * cos(alpha)  ;
	 break;
  case 4:
	 dPhidx1 =  fac * cos(theta/2.) * (cos(theta) + 2*sin(theta)*sin(theta)) ;
	 dPhidx2 = -fac * cos(theta) * sin(3*theta/2.);
	 (*answer)[0] = dPhidx1 * cos(alpha) - dPhidx2 * sin(alpha)  ;
	 (*answer)[1] = dPhidx1 * sin(alpha) + dPhidx2 * cos(alpha)  ;
	 break;
  }
  delete localCoord ; // Purify 27-10-2005
  delete globalCoord ;
  return answer;
}






