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


#include "bimatcrackasymp.h"

#include "crackasymptotic.h"
#include "auxiliaryfield.cpp"
#include "auxiliaryfield.h"
#include "material.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "cracktip.h"
#include "feinterpol.h"
#include <math.h>

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

class CrackTip;

BiMatCrackAsymp :: BiMatCrackAsymp(Domain* aDomain,int n,int m)
: CrackAsymptotic (aDomain,n),order(m){}


double BiMatCrackAsymp::computeBiMatConst()
// *****************************************
// Compute the bi material constant
{
  std::vector<Material*>*  matArray = 
	 dynamic_cast<CrackTip*>(activeEnrItem)->giveMatArray();

  // Material constants of material 1 and 2
  double E1  = (*matArray)[0]->give('E');
  double nu1 = (*matArray)[0]->give('n');
  double E2  = (*matArray)[1]->give('E');
  double nu2 = (*matArray)[1]->give('n');

  double shear1 = E1/(2. + nu1 + nu1);
  double shear2 = E2/(2. + nu2 + nu2);

  double K1 = Field<PlaneStrain>::K(nu1);     // the Kolosov coeficient for mat1
  double K2 = Field<PlaneStrain>::K(nu2);     // the Kolosov coeficient for mat2
  double beta    = ( shear1*(K2 - 1.0) - shear2*(K1 - 1.0) )/( shear1*(K2 + 1.0) + shear2*(K1 + 1.0) ) ;
  double epsilon = 1.0/(2.0*M_PI)*log( (1.0 - beta)/(1.0 + beta) );

  return epsilon;
}

double BiMatCrackAsymp::EvaluateYourSelfAt(GaussPoint* gp)
// ******************************************************
// Compute the i-th crack asymptotic function
{
  double epsilon = this->computeBiMatConst();

  Element *e = gp->giveMyElement();
  Mu::Point *globalCoord = 
	 e->giveFEInterpolation()->local2Global(domain,e->giveNodeArray(),gp->giveCoordinates());
  std::valarray<double> *localCoord = 
	 dynamic_cast<CrackTip*>(this->activeEnrItem)->computePolarCoordOf(globalCoord);
  double r     = (*localCoord)[0];
  double theta = (*localCoord)[1];

  double epsLogr = epsilon*log(r) ;
  double sqrtR = sqrt(r);
  double A = sqrtR * cos(epsLogr);
  double B = exp(epsilon*theta);
  double C = sqrtR * sin(epsLogr);

  double answer ;
  switch(this->order)
  { 
  case 1:
	 answer = A/B*sin(theta/2.);		break ;
  case 2:
	 answer = A/B*cos(theta/2.);		break ;
  case 3:
	 answer = A*B*sin(theta/2.);		break ;
  case 4:
	 answer = A*B*cos(theta/2.);		break ;
  case 5:
	 answer = A*B*sin(theta)*sin(theta/2.);	break ;
  case 6:
	 answer = A*B*sin(theta)*cos(theta/2.);	break ;
  case 7:
	 answer = C/B*sin(theta/2.);		break ;
  case 8:
	 answer = C/B*cos(theta/2.);		break ;
  case 9:
	 answer = C*B*sin(theta/2.);		break ;
  case 10:
	 answer = C*B*cos(theta/2.);		break ;
  case 11:
	 answer = C*B*sin(theta)*sin(theta/2.);	break ;
  case 12:
	 answer = C*B*sin(theta)*cos(theta/2.);	break ;
  }
  return answer ;
}

double BiMatCrackAsymp::EvaluateYourSelfAt(Element* e,Node* n)
// ***********************************************************
// Compute the i-th crack asymptotic function
// The first version, inorge the discontinous of asymtotic function throung the crack
// ham naobat lientuc (sin(theta/2) thi sua lai!!! 
{
  double epsilon = this->computeBiMatConst();
  Mu::Point* p = n->makePoint() ; 	

  std::valarray<double> *localCoord = 
	 dynamic_cast<CrackTip*>(this->activeEnrItem)->computePolarCoordOf(p) ;

  double r     = (*localCoord)[0];
  double theta = (*localCoord)[1];

  double epsLogr = epsilon*log(r) ;
  double sqrtR = sqrt(r);
  double A = sqrtR * cos(epsLogr);
  double B = exp(epsilon*theta);
  double C = sqrtR * sin(epsLogr);

  double answer ;
  switch(this->order)
  { 
  case 1:
	 answer = A/B*sin(theta/2.);		break ;
  case 2:
	 answer = A/B*cos(theta/2.);		break ;
  case 3:
	 answer = A*B*sin(theta/2.);		break ;
  case 4:
	 answer = A*B*cos(theta/2.);		break ;
  case 5:
	 answer = A*B*sin(theta)*sin(theta/2.);	break ;		
  case 6:
	 answer = A*B*sin(theta)*cos(theta/2.);	break ;
  case 7:
	 answer = C/B*sin(theta/2.);		break ;
  case 8:
	 answer = C/B*cos(theta/2.);		break ;
  case 9:
	 answer = C*B*sin(theta/2.);		break ;
  case 10:
	 answer = C*B*cos(theta/2.);		break ;
  case 11:
	 answer = C*B*sin(theta)*sin(theta/2.);	break ;
  case 12:
	 answer = C*B*sin(theta)*cos(theta/2.);	break ;
  }
  delete p;
  return answer ;
}
/*
double BiMatCrackAsymp::EvaluateYourSelfAt(Element* e,Node* n)
// Compute the i-th crack asymptotic function
// The second version, modify asymtotic function wich discontinous throung the crack
// ham naobat lientuc (sin(theta/2) thi sua lai!!! 
{
double epsilon = this->computeBiMatConst();
//std::cout<<"BiMatCrackAsymp::EvaluateYourSelfAt: epsilon = "<<epsilon<<std::endl; Debug only TQT
Mu::Point* p = n->makePoint() ; 	
FloatArray * localCoord = 
dynamic_cast<CrackTip*>(activeEnrItem)->computeLocalCoordOf(p) ;

double xloc = localCoord->at(1);
double yloc = localCoord->at(2);

double r = sqrt( xloc * xloc + yloc * yloc );
double theta = atan2(yloc,xloc) ;  

double delta;

//double discon = floor(1.0 + (theta - M_PI)*1000000000.0);
if (fabs(theta) == M_PI)
{		
//Above crack:	delta > 0.000001			EPSILON
Mu::Point *center = e->giveMyCenter(); 
Mu::Segment *tipSeg = dynamic_cast<CrackTip*>(this->activeEnrItem)->giveTipSegment();
// use orientation test to check the center's location w.r.t the tip segment 
double x1 = tipSeg->first()->x ;
double y1 = tipSeg->first()->y ;
double x2 = tipSeg->second()->x ;
double y2 = tipSeg->second()->y ;

//delta = (x1-center->x)*(y2-center->y) - (x2-center->x)*(y1-center->y) ;
if(x2 > x1)
delta = (x1-center->x)*(y2-center->y) - (x2-center->x)*(y1-center->y) ;
else
delta = (x2-center->x)*(y1-center->y) - (x1-center->x)*(y2-center->y) ;   

delete center ;


if (delta < -0.000001)		//This Node belong to the lower element!!! Does it synonym to the lower-half? TQT
{
theta = -M_PI;
//if ( ((e->giveNumber()) == 53 )&( (n->giveNumber()) == 59 ) )
//	std::cout<<"Element : "<<e->giveNumber()<<" Node "<<n->giveNumber()<<" is Low"<<std::endl;
//if ( (e->giveNumber()) >50 )
//	std::cout<<" E"<<e->giveNumber()<<" N"<<n->giveNumber();
}
else
{
theta = M_PI;
//if ( (e->giveNumber()) == 53 )
//	std::cout<<"Element : "<<e->giveNumber()<<" Node "<<n->giveNumber()<<" is up"<<std::endl;
//if ( (e->giveNumber()) >50 )
//	std::cout<<" E"<<e->giveNumber()<<" N"<<n->giveNumber();
}

}
double A = sqrt(r)*cos(epsilon*log(r));
double B = exp(epsilon*theta);
double C = sqrt(r)*sin(epsilon*log(r));

//Debug only TQT
if (e->giveNumber() == 212)
{
std::cout<<"( "<<n->giveCoordinate(1)<<" , "<<n->giveCoordinate(2)<<" ) - ";
std::cout<<"(xloc, yloc) = ("<<xloc<<"; "<<yloc<<" ) - (r, theta) = ("<<r<<" ; "<<theta*180.0/M_PI<<" )"<<std::endl;
}
//End debug!!!

double answer ;

switch(this->order)
{ 
case 1:
answer = A/B*sin(theta/2.);				break ;		
case 2:
answer = A/B*cos(theta/2.);				break ;
case 3:
answer = A*B*sin(theta/2.);				break ;
case 4:
answer = A*B*cos(theta/2.);				break ;
case 5:
answer = A*B*sin(theta)*sin(theta/2.);	break ;
case 6:
answer = A*B*sin(theta)*cos(theta/2.);	break ;
case 7:
answer = C/B*sin(theta/2.);				break ;
case 8:
answer = C/B*cos(theta/2.);				break ;
case 9:
answer = C*B*sin(theta/2.);				break ;
case 10:
answer = C*B*cos(theta/2.);				break ;
case 11:
answer = C*B*sin(theta)*sin(theta/2.);	break ;
case 12:
answer = C*B*sin(theta)*cos(theta/2.);	break ;
}
delete p;
return answer ;
}
*/
FloatArray* BiMatCrackAsymp::EvaluateYourGradAt(GaussPoint* gp)
// **************************************************************
// For The crack align to the Bi-Material interface
// Derivatives of 12 enrichment functions 
{
  // inclination angle of the local crack tip coord. sys w.r.t the x axe
  double alpha = static_cast<CrackTip*>(activeEnrItem)->giveInclinationAngle();
  // polar coordinates r and theta  
  Element *e = gp->giveMyElement();
  Mu::Point *globalCoord = 
	 e->giveFEInterpolation()->local2Global(domain,e->giveNodeArray(),gp->giveCoordinates());
  std::valarray<double> *localCoord = 
	 dynamic_cast<CrackTip*>(this->activeEnrItem)->computePolarCoordOf(globalCoord);

  double r     = (*localCoord)[0];
  double theta = (*localCoord)[1];

  double epsilon = this->computeBiMatConst();

  double epsLogr = epsilon*log(r) ;
  double sqrtR = sqrt(r);

  double A = cos(epsLogr);
  double B = sin(epsLogr);
  double C = 2.*epsilon*cos(epsLogr - theta);
  double D = 2.*epsilon*sin(epsLogr - theta);
  double E = 2.*epsilon*cos(epsLogr + theta);
  double F = 2.*epsilon*sin(epsLogr + theta);
  double G = exp(epsilon*theta);

  double ST   = sin(theta)    ;
  double CT   = cos(theta)    ;
  double ST2  = sin(theta/2.) ;
  double CT2  = cos(theta/2.) ;
  double S3T2 = sin(3. * theta/2.) ;				  
  double C3T2 = cos(3. * theta/2.) ;
  double fac  = 1./( 2. * sqrtR ) ;

  double dPhidx1,dPhidx2 ;
  FloatArray* answer = new FloatArray(2); 

  switch(this->order)
  {
  case 1:
	 dPhidx1 = - fac/G * ST2 * (A + D); 
	 dPhidx2 =   fac/G * (A * CT2 - C * ST2);
	 (*answer)[0] = dPhidx1 * cos(alpha) - dPhidx2 * sin(alpha)  ;
	 (*answer)[1] = dPhidx1 * sin(alpha) + dPhidx2 * cos(alpha)  ;
	 break;
  case 2:
	 dPhidx1 =   fac/G * CT2 * (A - D);
	 dPhidx2 =   fac/G * (A * ST2 - C * CT2);
	 (*answer)[0] = dPhidx1 * cos(alpha) - dPhidx2 * sin(alpha)  ;
	 (*answer)[1] = dPhidx1 * sin(alpha) + dPhidx2 * cos(alpha)  ;
	 break;
  case 3:
	 dPhidx1 = - fac * G * ST2 * (A + F);
	 dPhidx2 =   fac * G * (A * CT2 + E * ST2);
	 (*answer)[0] = dPhidx1 * cos(alpha) - dPhidx2 * sin(alpha)  ;
	 (*answer)[1] = dPhidx1 * sin(alpha) + dPhidx2 * cos(alpha)  ;
	 break;
  case 4:
	 dPhidx1 =   fac * G * CT2 * (A - F);
	 dPhidx2 =   fac * G * (A * ST2 + E * CT2);
	 (*answer)[0] = dPhidx1 * cos(alpha) - dPhidx2 * sin(alpha)  ;
	 (*answer)[1] = dPhidx1 * sin(alpha) + dPhidx2 * cos(alpha)  ;
	 break;
  case 5:
	 dPhidx1 = - fac * G * ST * (A * S3T2 + F * ST2);
	 dPhidx2 =   fac * G * (A * (ST2 + S3T2 * CT) + E * ST2 * ST);
	 (*answer)[0] = dPhidx1 * cos(alpha) - dPhidx2 * sin(alpha)  ;
	 (*answer)[1] = dPhidx1 * sin(alpha) + dPhidx2 * cos(alpha)  ;
	 break;
  case 6:
	 dPhidx1 = - fac * G * ST * (A * C3T2 + F * CT2);
	 dPhidx2 =   fac * G * (A * (CT2 + C3T2 * CT) + E * CT2 * ST);
	 (*answer)[0] = dPhidx1 * cos(alpha) - dPhidx2 * sin(alpha)  ;
	 (*answer)[1] = dPhidx1 * sin(alpha) + dPhidx2 * cos(alpha)  ;
	 break;
  case 7:
	 dPhidx1 = - fac/G * ST2 * (B - C);
	 dPhidx2 =   fac/G * (B * CT2 - D * ST2);
	 (*answer)[0] = dPhidx1 * cos(alpha) - dPhidx2 * sin(alpha)  ;
	 (*answer)[1] = dPhidx1 * sin(alpha) + dPhidx2 * cos(alpha)  ;
	 break;
  case 8:
	 dPhidx1 =  fac/G * CT2 * (B + C);
	 dPhidx2 =  fac/G * (B * ST2 - D * CT2);
	 (*answer)[0] = dPhidx1 * cos(alpha) - dPhidx2 * sin(alpha)  ;
	 (*answer)[1] = dPhidx1 * sin(alpha) + dPhidx2 * cos(alpha)  ;
	 break;
  case 9:
	 dPhidx1 = -fac * G * ST2 * (B - E);
	 dPhidx2 =  fac * G * (B * CT2 + F * ST2);
	 (*answer)[0] = dPhidx1 * cos(alpha) - dPhidx2 * sin(alpha)  ;
	 (*answer)[1] = dPhidx1 * sin(alpha) + dPhidx2 * cos(alpha)  ;
	 break;
  case 10:
	 dPhidx1 =  fac * G * CT2 * (B + E);
	 dPhidx2 =  fac * G * (B * ST2 + F * CT2);
	 (*answer)[0] = dPhidx1 * cos(alpha) - dPhidx2 * sin(alpha)  ;
	 (*answer)[1] = dPhidx1 * sin(alpha) + dPhidx2 * cos(alpha)  ;
	 break;
  case 11:
	 dPhidx1 = -fac * G * ST * (B * S3T2 - E *ST2);
	 dPhidx2 =  fac * G * (B * (ST2 + S3T2 * CT) + F * ST2 * ST);
	 (*answer)[0] = dPhidx1 * cos(alpha) - dPhidx2 * sin(alpha)  ;
	 (*answer)[1] = dPhidx1 * sin(alpha) + dPhidx2 * cos(alpha)  ;
	 break;
  case 12:
	 dPhidx1 = -fac * G * ST * (B * C3T2 - E * CT2);
	 dPhidx2 =  fac * G * (B * (CT2 + C3T2 * CT) + F * CT2 * ST);
	 (*answer)[0] = dPhidx1 * cos(alpha) - dPhidx2 * sin(alpha)  ;
	 (*answer)[1] = dPhidx1 * sin(alpha) + dPhidx2 * cos(alpha)  ;
	 break;
  }
  return answer;
}