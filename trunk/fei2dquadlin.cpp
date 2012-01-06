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
   
	Any feedback is welcome. Emails : nvphu80@yahoo.com, ...

*************************************************************************************/

#include "fei2dquadlin.h"

#include "flotmtrx.h"
#include "flotarry.h"
#include "geometry_base.h"

FloatArray*  FEI2dQuadLin::evalN(Mu::Point* p)
// *******************************************
// computes the 4 shape functions at point p.
{
	FloatArray* answer = new FloatArray(4) ;

	(*answer)[0] = (1. - p->x) * (1. - p->y) * 0.25 ;
	(*answer)[1] = (1. + p->x) * (1. - p->y) * 0.25 ;
	(*answer)[2] = (1. + p->x) * (1. + p->y) * 0.25 ;
	(*answer)[3] = (1. - p->x) * (1. + p->y) * 0.25 ;

	return answer ;
}

FloatMatrix* FEI2dQuadLin :: evaldNdx(Domain* d,IntArray* Nodes,Mu::Point* lcoords)
// ********************************************************************************
// computes the derivatives of shape functions w.r.t global coordinates
// dNdx = [N1,x N1,y
//         N2,x N2,y
//         N3,x N3,y
//         N4,x N4,y]
{

	FloatMatrix * jacobianMatrix = this->giveJacobianMatrixAt(d,Nodes,lcoords);
	FloatMatrix * inv = jacobianMatrix->GiveInverse();

	FloatArray *dNdKsi = this->giveDerivativeKsi(lcoords->y); // derivative of N w.r.t Ksi
	FloatArray *dNdEta = this->giveDerivativeEta(lcoords->x); // derivative of N w.r.t Eta

	FloatMatrix *answer = new FloatMatrix(4,2);

	for (size_t i = 0 ; i < 4 ; i++) 
	{
		answer->at(i+1,1) = (*dNdKsi)[i]*inv->at(1,1) + (*dNdEta)[i]*inv->at(1,2);
		answer->at(i+1,2) = (*dNdKsi)[i]*inv->at(2,1) + (*dNdEta)[i]*inv->at(2,2);
	}

	delete jacobianMatrix ;
	delete inv ;
	delete dNdKsi ; delete dNdEta ; // Purify, 11-10-05
	 
	return answer;
}

FloatArray*  FEI2dQuadLin :: giveDerivativeKsi(double eta,double ksi)
// ******************************************************************
// computes the derivatives of shape functions w.r.t Ksi 
{

	FloatArray* answer = new FloatArray(4);

	(*answer)[0] = (-1. + eta) * 0.25;
	(*answer)[1] = ( 1. - eta) * 0.25;
	(*answer)[2] = ( 1. + eta) * 0.25;
	(*answer)[3] = (-1. - eta) * 0.25;

	return answer ;
}

FloatArray*  FEI2dQuadLin :: giveDerivativeEta(double ksi,double eta)
// ******************************************************************
// computes the derivatives of shape functions w.r.t Eta 
{
	FloatArray* answer = new FloatArray(4);

	(*answer)[0] = (-1. + ksi) * 0.25;
	(*answer)[1] = (-1. - ksi) * 0.25;
	(*answer)[2] = ( 1. + ksi) * 0.25;
	(*answer)[3] = ( 1. - ksi) * 0.25;

	return answer ;
}

FloatMatrix* FEI2dQuadLin :: giveJacobianMatrixAt(Domain* d,IntArray* nodes,Mu::Point* lcoords)
// ********************************************************************************************
// computes the Jacobian matrix J
{
	double ksi = lcoords->x;
	double eta = lcoords->y;

	FloatArray *dNdKsi = this->giveDerivativeKsi(eta);
	FloatArray *dNdEta = this->giveDerivativeEta(ksi);

	FloatMatrix *jacobianMatrix = new FloatMatrix(2,2);

	for (size_t i = 0 ; i < 4 ; i++)
	{
		double x = d -> giveNode((*nodes)[i]) -> giveCoordinate(xind);
		double y = d -> giveNode((*nodes)[i]) -> giveCoordinate(yind);

		jacobianMatrix->at(1,1) += (*dNdKsi)[i]*x;
		jacobianMatrix->at(1,2) += (*dNdKsi)[i]*y;
		jacobianMatrix->at(2,1) += (*dNdEta)[i]*x;
		jacobianMatrix->at(2,2) += (*dNdEta)[i]*y;
	}

	delete dNdKsi ; delete dNdEta ; 

	return jacobianMatrix ;
}

Mu::Point* FEI2dQuadLin::local2Global(Domain* d,IntArray* nodes,Mu::Point* lcoords)
// ********************************************************************************
// Compute the global coordinates of the local point lcoords using the usual finite
// element interpolation: \f$ x=\sumN_{i}x_{i} \f$
{
	Mu::Point* answer = new Mu::Point();

	for (int i = 0 ; i < 4 ; i++)
	{
		double x = d -> giveNode((*nodes)[i])->giveCoordinate(1);
		double y = d -> giveNode((*nodes)[i])->giveCoordinate(2);
		
		answer->x += (*this->evalN(lcoords))[i]*x;
		answer->y += (*this->evalN(lcoords))[i]*y;
	}
	return answer;
}

Mu::Point* FEI2dQuadLin::global2Local(Domain* d,IntArray* nodes,Mu::Point* p)
// **************************************************************************
// Could be used with Newton Raphson method or solve by hand as follow
// The equations system need to be solved :
//             a1-a2*r-a3*s+a4*r*s = 4*xo
//             b1-b2*r-b3*s+b4*r*s = 4*yo
// (xo,yo) are coordinates of point p
// (r,s)  are unknows to be determined, which are the local coord. of point p
// Written in 2005-08-22 by NVP and Truong Quang Tri.
{
   Node* node1 = d -> giveNode((*nodes)[0]) ;
   Node* node2 = d -> giveNode((*nodes)[1]) ;
   Node* node3 = d -> giveNode((*nodes)[2]) ;
   Node* node4 = d -> giveNode((*nodes)[3]) ;

   double x1 = node1 -> giveCoordinate(1) ;
   double x2 = node2 -> giveCoordinate(1) ;
   double x3 = node3 -> giveCoordinate(1) ;
   double x4 = node4 -> giveCoordinate(1) ;

   double y1 = node1 -> giveCoordinate(2) ;
   double y2 = node2 -> giveCoordinate(2) ;
   double y3 = node3 -> giveCoordinate(2) ;
   double y4 = node4 -> giveCoordinate(2) ;

	double a1 = x1 + x2 + x3 + x4 ;  // coefficients of the first equation.
	double a2 = x1 - x2 - x3 + x4 ;
	double a3 = x1 + x2 - x3 - x4 ;
	double a4 = x1 - x2 + x3 - x4 ;

	double b1 = y1 + y2 + y3 + y4 ;  // coefficients of the second equation.
	double b2 = y1 - y2 - y3 + y4 ;
	double b3 = y1 + y2 - y3 - y4 ;
	double b4 = y1 - y2 + y3 - y4 ;

	// from these two equations, deduce the quadratic equation : a*r^2 - b*r + c = 0
	double a = a2 * b4 - b2 * a4 ;
	double b = a1 * b4 + a2 * b3 - a3 * b2 - b1 * a4 - b4 * 4.0 * p->x + a4 * 4.0 * p->y ;
	double c = a1 * b3 - a3 * b1 - 4.0 * p->x * b3 + 4.0 * p->y * a3 ;

	// solve this quadratic equation for \ksi ...
	double ksi1, ksi2 ;   // two solutions
	double ksi,eta ;      // correct solution
	double delta = b * b - 4.0 * a * c ;

   if(fabs(a) < 10e-6)      // a considered to be zero, linear equation
	{
	  ksi1 = c/b ;
	  ksi2 = c/b ;
	  ksi  = ksi1 ;
	}                  // quadratic equation  
	else if(delta > 0)    // two solutions, i.e., two points
	{
	  ksi1 = ( b + sqrt(delta) )/(2.0 * a);
     ksi2 = ( b - sqrt(delta) )/(2.0 * a);
	  // choose the internal points,i.e.,  -1 <= ksi <=1
	  if( (ksi1 >= -1) && (ksi1 <=1))
		 ksi = ksi1 ;
	  else if( (ksi2 >= -1) && (ksi2 <=1))
		 ksi = ksi2 ;
	  else
	  {
		 std::cout << " Hmm. Two points are both outside !!! " << std::endl ;
		 assert(false);
	  }
	}
	else if(delta == 0)   // two real coincident solutions
	{
	  ksi1 = b/(2.0 * a) ;
	  ksi2 = b/(2.0 * a) ;
	  ksi  = ksi1 ;
	}
	else
	{
	  std::cout << " Two complex conjugate solutions :" << std::endl ;
	}

	// then, for \eta
	
	double denom = b3 + ksi * b4;
	
	if (fabs(denom) <= 10e-10)
	{
		eta = (4.0 * p->x - a1 - a2 * ksi)/(a4 * ksi + a3);
	}
	else
	{
        eta = (4.0 * p->y - b1 - b2 * ksi)/denom;
	};
	
   Mu::Point* answer = new Mu::Point(ksi,eta);

	return answer;
}


