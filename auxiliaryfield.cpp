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

#include "auxiliaryfield.h"

#include "flotmtrx.h"
#include "cracktip.h"
#include "material.h"
#include "assert.h"
#include <iostream>

#ifndef M_PI 
#define M_PI  3.14159265358979323846
#endif


/*
template<class M1, class M2,const FieldType field>
AuxiliaryField<M1,M2,field>::AuxiliaryField()
{
  if(mode == 0)
	 std::cout << " Auxiliary field for mode I is being created " << std::endl ;
  else if(mode == 1)
	 std::cout << " Auxiliary field for mode II is being created " << std::endl ;
}
*/

template<class M1, class M2, const FieldType field>
void AuxiliaryField<M1,M2,field>::ComputeComponentsOfAuxField(CrackTip* tip, Mu::Point* p,ModeType mode ,FloatMatrix& AuxGradDisp,FloatArray& AuxEps,FloatArray& AuxStress)
// *************************************************************************************************************************************************************
// compute the auxiliary field 
{

  switch (this->giveNumberOfMaterials())
  {
  case 1:
	 this->ComputeComponentsOfOneMat(tip,p,mode,AuxGradDisp,AuxEps,AuxStress);
	 break;
  case 2:
	 this->ComputeComponentsOfBiMat(tip,p,mode,AuxGradDisp,AuxEps,AuxStress);
	 break;
  }
}


template<class M1, class M2,const FieldType field>
void AuxiliaryField<M1,M2,field>::
ComputeComponentsOfOneMat(CrackTip* tip, Mu::Point* p, ModeType mode ,FloatMatrix& AuxGradDisp,FloatArray& AuxEps,FloatArray& AuxStress)
// ***********************************************************************************************************************
// computes the derivatives of displacement field of the auxiliary field 
// for homogeneous elastic cracks.
// AuxGradDisp = [du1dx1 du1dx2
//                du2dx1 du2dx2]
// AuxStress   = [sigma_xx sigma_yy sigma_xy]
{
  
  std::vector<Material*>*  matArray = tip->giveMatArray();

  // Material constants of material 
  double E  = (*matArray)[0]->give('E');
  double nu = (*matArray)[0]->give('n');

  // get the coords in the crack tip local coordinate system

  std::valarray<double> *localCoord = tip->computePolarCoordOf(p) ;

  double r     = (*localCoord)[0];
  double theta = (*localCoord)[1];

  // for sure, -Pi <= theta <= Pi
  assert(theta <= M_PI && theta >= -M_PI);

  double SQR = sqrt(r);  
  double CT  = cos(theta);
  double ST  = sin(theta);
  double CT2 = cos(theta/2.);
  double ST2 = sin(theta/2.);
  double C3T2= cos(3*theta/2.);  
  double S3T2= sin(3*theta/2.);

  // derivatives of r, theta w.r.t local coordinate system
  double drdx = CT;
  double drdy = ST;
  double dtdx = -ST/r;
  double dtdy = CT/r;

  double K  = Field<field>::K(nu);     // the Kolosov coefficient 
  double mu = E/(2.+ nu + nu);
  
  double  FACDisp   = sqrt(1/(2*M_PI))/(2*mu);
  double  FACStress = 1./( sqrt(2.*(M_PI)*r) );   

  
  double du1dr,du1dt,du2dr,du2dt;

  if(mode == 0)      // mode I
  {
	 du1dr = FACDisp*0.5/SQR*CT2*(K - CT);
	 du1dt = FACDisp*SQR*(-0.5*ST2*(K - CT) + CT2*ST);

	 du2dr = FACDisp*0.5/SQR*ST2*(K - CT);
	 du2dt = FACDisp*SQR*(0.5*CT2*(K - CT) + ST2*ST);

	 // compute derivative of displacement field of auxiliary field

	 AuxGradDisp.at(1,1) = du1dr*drdx + du1dt*dtdx ;
	 AuxGradDisp.at(1,2) = du1dr*drdy + du1dt*dtdy ;
	 AuxGradDisp.at(2,1) = du2dr*drdx + du2dt*dtdx ;
	 AuxGradDisp.at(2,2) = du2dr*drdy + du2dt*dtdy ;

	 // compute strain vector of auxiliary field

	 AuxEps.at(1) = AuxGradDisp.at(1,1) ; // epsilon_xx
	 AuxEps.at(2) = AuxGradDisp.at(2,2) ; // epsilon_yy
	 AuxEps.at(3) = (0.5)*(AuxGradDisp.at(1,2) + AuxGradDisp.at(2,1)); // epsilon_xy

	 // compute stress vector

	 AuxStress.at(1) = FACStress * CT2 * (1. - ST2 * S3T2) ; // sigma_xx
	 AuxStress.at(2) = FACStress * CT2 * (1. + ST2 * S3T2) ; // sigma_yy
	 AuxStress.at(3) = FACStress * ST2 * CT2 * C3T2 ;        // sigma_xy

  }
  else if(mode == 1) // mode II
  {
	 du1dr = FACDisp * 0.5/SQR * ST2 * (K + 2 + CT);
	 du1dt = FACDisp * SQR * (0.5 * CT2 * (K + 2 + CT) - ST2*ST);

	 du2dr = -FACDisp * 0.5*(1/SQR)*CT2*(K - 2 + CT);
	 du2dt = -FACDisp * SQR*(-0.5*ST2*(K - 2 + CT) - CT2*ST);

	 // compute derivative of displacement field of auxiliary field

	 AuxGradDisp.at(1,1) = du1dr*drdx + du1dt*dtdx ;
	 AuxGradDisp.at(1,2) = du1dr*drdy + du1dt*dtdy ;
	 AuxGradDisp.at(2,1) = du2dr*drdx + du2dt*dtdx ;
	 AuxGradDisp.at(2,2) = du2dr*drdy + du2dt*dtdy ;

	 // compute strain vector of auxiliary field

	 AuxEps.at(1) = AuxGradDisp.at(1,1) ; // epsilon_xx
	 AuxEps.at(2) = AuxGradDisp.at(2,2) ; // epsilon_yy
	 AuxEps.at(3) = (0.5)*(AuxGradDisp.at(1,2) + AuxGradDisp.at(2,1)); // epsilon_xy

	 AuxStress.at(1) = -FACStress * ST2 * (2. + CT2 * C3T2) ;  // sigma_xx
	 AuxStress.at(2) =  FACStress * ST2 * CT2 * C3T2 ;         // sigma_yy
	 AuxStress.at(3) =  FACStress * CT2 * (1. - ST2 * S3T2) ;  // sigma_xy
  }
  else
	 assert(false);

  delete localCoord ;

}

template<class M1, class M2, const FieldType field>
void AuxiliaryField<M1,M2,field>::
ComputeComponentsOfBiMat(CrackTip* tip, Mu::Point* p, ModeType mode ,FloatMatrix& AuxGradDisp,FloatArray& AuxEps,FloatArray& AuxStress)
// **********************************************************************************************************************
// computes the auxiliary field for interfacial cracks
{
  
  std::vector<Material*>*  matArray = tip->giveMatArray();

  // Material constants of material 1 and 2

  double E1  = (*matArray)[0]->give('E');
  double nu1 = (*matArray)[0]->give('n');
  double E2  = (*matArray)[1]->give('E');
  double nu2 = (*matArray)[1]->give('n');

  // get the coords in the crack tip local coordinate system  

  FloatArray * localCoord = tip->computeLocalCoordOf(p) ;


  double xloc = localCoord->at(1);
  double yloc = localCoord->at(2);

  double r = sqrt( xloc * xloc + yloc * yloc );
  double theta = atan2(yloc,xloc) ;
  
  // for sure, -Pi <= theta <= Pi
  assert(theta <= M_PI && theta >= -M_PI);

  double shear1 = E1 / (2.+nu1+nu1);
  double shear2 = E2 / (2.+nu2+nu2);	

  double K1 = Field<field>::K(nu1);     // the Kolosov coeficient 
  double K2 = Field<field>::K(nu2);     // the Kolosov coeficient
  
  double epsilon = 0.5/M_PI*log((shear1 + shear2 * K1)/(shear2 + shear1 * K2));
  double epLogr = epsilon*log(r);
  double cosEpLogr = cos(epLogr);
  double sinEpLogr = sin(epLogr);
  double temp = 0.25 + epsilon * epsilon ;
  double beta  = (0.5 * cosEpLogr + epsilon * sinEpLogr)/temp ;
  double beta_ = (0.5 * sinEpLogr - epsilon * cosEpLogr)/temp ;
  double phi = epLogr + theta/2.;

  double delta, gamma,gamma_ , A;

  double ST2 = sin(theta/2.);
  double CT2 = cos(theta/2.);
  double CT  = cos(theta);
  double ST  = sin(theta);
  double SP  = sin(phi) ;
  double CP  = cos(phi) ;

  //Need to differentiate the upper and lower half plane

  Mu::Segment *tipSeg = tip->giveTipSegment();
  // use orientation test to check the center's location w.r.t the tip segment 
  double x1 = tipSeg->first()->x ;
  double y1 = tipSeg->first()->y ;
  double x2 = tipSeg->second()->x ;
  double y2 = tipSeg->second()->y ;

  //delta = (x1-center->x)*(y2-center->y) - (x2-center->x)*(y1-center->y) ;
  double deltaSig;
  if(x2 > x1)
	  deltaSig = (x1 - p->x)*(y2 - p->y) - (x2 - p->x)*(y1 - p->y) ;
  else
	  deltaSig = (x2 - p->x)*(y1 - p->y) - (x1 - p->x)*(y2 - p->y) ;   
  // End Here

  //if((theta>=0)&&(theta<=M_PI))	// Belong to the upper-half: Material M1
  if (deltaSig > 0.0000001)
  { 
	 delta = exp((theta - M_PI)*epsilon);
	 gamma = K1 * delta - 1./delta ;
	 gamma_= K1 * delta + 1./delta ;
	 A = 1./(4.*shear1*cosh(M_PI*epsilon));
  }
  else if (deltaSig < -0.0000001)
  {                          // Belong to the lower-half: Material M2
	 delta = exp((theta + M_PI)*epsilon);
	 gamma = K2 * delta - 1./delta;
	 gamma_= K2 * delta + 1./delta;
	 A = 1./(4.*shear2*cosh(M_PI*epsilon));
  }
  else
	 assert(false);

  double T1 = 2. * delta * ST * SP ;
  double T2 = 2. * delta * ST * CP ;
  double T3 = 2. * delta * CT * SP ;
  double T4 = 2. * delta * CT * CP ;

  double C = beta_* gamma * CT2 - beta * gamma_* ST2;
  double D = beta * gamma * CT2 + beta_* gamma_* ST2;
  double E = beta_* gamma_* CT2 - beta * gamma * ST2;		
  double F = beta * gamma_* CT2 + beta_* gamma * ST2;

  double dT1dr = epsilon*T2/r;
  double dT1dt = epsilon*T1 + T2/2.0 + T3;
  double dT2dr = - epsilon*T1/r;
  double dT2dt = epsilon*T2 - T1/2.0 + T4;

  double dCdr = epsilon*D/r;
  double dCdt = - F/2.0 + epsilon*E;		
  double dDdr = -epsilon*C/r;
  double dDdt = E/2.0 + epsilon*F;		

  double f1,f2,df1dr,df1dt,df2dr,df2dt;

  switch (mode){ // mode I or mode II	
	case 0:                           // Mode I
	  f1 = D + T1;
	  f2 = - C - T2;

	  df1dr = dDdr + dT1dr ;
	  df1dt = dDdt + dT1dt ;
	  df2dr = - dCdr - dT2dr ;
	  df2dt = - dCdt - dT2dt ;			 
	  break;
	case 1:                           // Mode II	
	  f1 = - C + T2;
	  f2 = - D + T1;

	  df1dr = - dCdr + dT2dr;
	  df1dt = - dCdt + dT2dt;			 
	  df2dr = - dDdr + dT1dr;
	  df2dt = - dDdt + dT1dt;

	  break;
  }

  double df1dx = df1dr*CT - df1dt*ST/r;
  double df1dy = df1dr*ST + df1dt*CT/r;
  double df2dx = df2dr*CT - df2dt*ST/r;
  double df2dy = df2dr*ST + df2dt*CT/r;

  double B     = sqrt(r/(M_PI+M_PI));

  AuxGradDisp.at(1,1) = A*(B*df1dx + (f1*CT)/(4.0*M_PI*B));
  AuxGradDisp.at(1,2) = A*(B*df1dy + (f1*ST)/(4.0*M_PI*B));
  AuxGradDisp.at(2,1) = A*(B*df2dx + (f2*CT)/(4.0*M_PI*B));
  AuxGradDisp.at(2,2) = A*(B*df2dy + (f2*ST)/(4.0*M_PI*B));

  // compute strain vector of auxiliary field
  
  AuxEps.at(1) = AuxGradDisp.at(1,1) ; // epsilon_xx
  AuxEps.at(2) = AuxGradDisp.at(2,2) ; // epsilon_yy
  AuxEps.at(3) = (0.5)*(AuxGradDisp.at(1,2) + AuxGradDisp.at(2,1)); // epsilon_xy

  // compute stress vector of auxiliary field 

  // Plane strain
  // Should use Template!! TQT
  double ee, shear, nu ;
  if (deltaSig > 0.0000001)		//Material 1
  {
	 ee    = E1 / ((1.+nu1) * (1.-nu1-nu1)) ;
	 shear = E1 / (2.+nu1+nu1) ;
	 nu    = nu1;
  }
  else if (deltaSig < -0.0000001)
  {
	 ee    = E2 / ((1.+nu2) * (1.-nu2-nu2)) ;
	 shear = E2 / (2.+nu2+nu2) ;
	 nu    = nu2;
  }
  else
	 assert(false);

   FloatMatrix* constitutiveMatrix = new FloatMatrix(3,3) ;

   constitutiveMatrix->at(1,1) = (1.-nu) * ee ;
   constitutiveMatrix->at(1,2) =     nu  * ee ;
   constitutiveMatrix->at(2,1) =     nu  * ee ;
   constitutiveMatrix->at(2,2) = (1.-nu) * ee ;
   constitutiveMatrix->at(3,3) =  shear ;

   //AuxStress = constitutiveMatrix->Times(&AuxEps) ;
   AuxStress.at(1) = ee*( (1.0 - nu)*AuxEps.at(1) + nu*AuxEps.at(2) ) ;
   AuxStress.at(2) = ee*( nu*AuxEps.at(1) + (1.0 - nu)*AuxEps.at(2) ) ;
   AuxStress.at(3) = shear * AuxEps.at(3) ;

   //End Plane Strain

/*  // Plane stress
  double ee, shear, nu ;  
  if (deltaSig > 0.0000001)		//Material 1
  {
	  nu    = nu1;
	  ee    = E1 / (1. - nu*nu) ;
	  shear = E1 / (2. + nu+nu) ;	 
  }
  else
  {
	  nu    = nu2;
	  ee    = E2 / (1. - nu*nu);
	  shear = E2 / (2. + nu+nu) ;	 
  }
   FloatMatrix* constitutiveMatrix = new FloatMatrix(3,3) ;

   constitutiveMatrix->at(1,1) = ee ;
   constitutiveMatrix->at(1,2) =     nu * ee ;
   constitutiveMatrix->at(2,1) =     nu * ee ;
   constitutiveMatrix->at(2,2) = ee ;
   constitutiveMatrix->at(3,3) = shear ;	
   //AuxStress = constitutiveMatrix->Times(&AuxEps) ;
   AuxStress.at(1) = ee*( AuxEps.at(1)    + nu*AuxEps.at(2) ) ;
   AuxStress.at(2) = ee*( nu*AuxEps.at(1) + AuxEps.at(2)    ) ;
   AuxStress.at(3) = ee*0.5*(1. - nu)*AuxEps.at(3) ;

   // End Plane Stress 
*/
}
