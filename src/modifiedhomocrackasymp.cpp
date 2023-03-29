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

#include "modifiedhomocrackasymp.h"

#include "crackasymptotic.h"
#include "flotarry.h"
#include "piecewiselinear.h"
#include "cracktip.h"
#include "crackinterior.h"
#include "assert.h"
#include <stdlib.h>
#include <math.h>


#ifndef M_PI 
#define M_PI  3.14159265358979323846
#endif


ModifiedHomoCrackAsymp :: ModifiedHomoCrackAsymp(Domain* aDomain,int n)
       : CrackAsymptotic (aDomain,n)
{
}

std::valarray<double> ModifiedHomoCrackAsymp::computeRampFunctionAt(Mu::Point* p)
{
   CrackInterior* myCrack = 
	   dynamic_cast<CrackTip*>(this->activeEnrItem)->giveMyCrackInterior();
	
	double a = myCrack->giveCrackLength() ; // half crack length

	FloatArray * localCoord = 
	  dynamic_cast<CrackTip*>(this->activeEnrItem)->computeLocalCoordOf(p) ;

	double xloc = (*localCoord)[0];
	double yloc = (*localCoord)[1];

	// compute ramp function and its local derivatives
	double ramp,dRdx1 ;
	if(xloc <= -2.0*a)
	{
		ramp  = 0.0 ;
		dRdx1 = 0.0 ;
	}
	else if(xloc <= 0.0 && xloc >= -2.0*a)
	{
	   ramp = 0.25*(xloc+2.0*a)*(xloc+2.0*a)*(a-xloc)/(a*a*a) ;
      dRdx1 = -0.75*(xloc+2.0*a)*xloc/(a*a*a) ;
	}
	else
	{
	  ramp  = 1.0 ;
	  dRdx1 = 0.0 ;
	}

	std::valarray<double> result(3);
	result[0] = ramp ;
	result[1] = dRdx1 ;
	result[2] = 0.0 ;

	return  result ;
}

double ModifiedHomoCrackAsymp::EvaluateYourSelfAt(Mu::Point* coord)
// **************************************************************
// compute the modified first branch function:
// \sqrt{r}sin(theta/2)*ramp function
// used to model crack nucleation, see Dolbow and Bellec, 2003 for details.
{
   std::valarray<double> *localCoord = 
	  dynamic_cast<CrackTip*>(this->activeEnrItem)->computePolarCoordOf(coord) ;

	double r     = (*localCoord)[0];
	double theta = (*localCoord)[1];

	// for sure, -Pi <= theta <= Pi
	assert(theta <= M_PI && theta >= -M_PI);

	std::valarray<double> ramp = this->computeRampFunctionAt(coord);

	double ret = sqrt(r) * sin(theta/2.) * ramp[0] ; 
	return ret ;
}

FloatArray* ModifiedHomoCrackAsymp::EvaluateYourGradAt(Mu::Point* p)
// **************************************************************
// compute the grad of the modified asymp. function: phi1(x)*r(x)
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

	double r2 = sqrt(r) ; // r = 0 ???
	double fac  = 0.5/r2 ;

   // local derivatives of first branch function
	double dPhi1dx1 = -fac * sin(theta/2.);
	double dPhi1dx2 =  fac * cos(theta/2.);

	// local derivatives of ramp function
	
	double Phi1 = sqrt(r) * sin(theta/2.) ;
   std::valarray<double> ramp = this->computeRampFunctionAt(p);

	double dPhidx1 = dPhi1dx1 * ramp[0] + Phi1 * ramp[1] ;
	double dPhidx2 = dPhi1dx2 * ramp[0] + Phi1 * ramp[2] ;

	FloatArray* answer = new FloatArray(2); 
	(*answer)[0] = dPhidx1 * cos(alpha) - dPhidx2 * sin(alpha)  ;
	(*answer)[1] = dPhidx1 * sin(alpha) + dPhidx2 * cos(alpha)  ;

	return answer;
}

double ModifiedHomoCrackAsymp::EvaluateYourSelfAt(Element* e ,Node* n)
// *******************************************************************
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
  
	std::valarray<double> ramp = this->computeRampFunctionAt(p);
	double answer  = sqrt(r) * sin(theta/2.) * ramp[0] ; 
   delete p ; delete localCoord ;
	return answer ;
}
	
	
	


