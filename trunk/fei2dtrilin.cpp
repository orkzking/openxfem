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

#include "fei2dtrilin.h"

#include "flotmtrx.h"
#include "flotarry.h"
#include "geometry_base.h"


FloatArray*  FEI2dTriLin::evalN(Mu::Point* lcoords)
// ************************************************
// computes the value of shape functions at a given point with coord. lcoords
{
	double n1  =  1. - lcoords->x - lcoords->y;
	double n2  =  lcoords->x;
	double n3  =  lcoords->y;

	FloatArray* answer = new FloatArray(3) ;

	(*answer)[0] = n1 ; 
	(*answer)[1] = n2 ;
	(*answer)[2] = n3 ;

   return answer ;
}

FloatMatrix* FEI2dTriLin :: evaldNdx(Domain* d,IntArray* Nodes,Mu::Point* lcoords)
// ********************************************************************************
// computes the derivatives of shape functions w.r.t global coordinates
// dNdx = [N1,x N1,y
//         N2,x N2,y
//         N3,x N3,y]
{
	FloatMatrix *jacobianMatrix = giveJacobianMatrixAt(d,Nodes,lcoords);
	
	FloatMatrix *inv = jacobianMatrix->GiveInverse();

	FloatArray *dNdKsi = giveDerivativeKsi();
	FloatArray *dNdEta = giveDerivativeEta();

	FloatMatrix * answer = new FloatMatrix(3,2);

	for (size_t i = 0 ; i < 3 ; i++)
	{
		answer->at(i+1,1) = (*dNdKsi)[i]*inv->at(1,1) + (*dNdEta)[i]*inv->at(1,2);
		answer->at(i+1,2) = (*dNdKsi)[i]*inv->at(2,1) + (*dNdEta)[i]*inv->at(2,2);
	}

	delete jacobianMatrix ;
	delete inv ;
	delete dNdKsi ; // Purify, 11-10-05
	delete dNdEta ; // Purify, 11-10-05

	return answer;
}

FloatArray*  FEI2dTriLin :: giveDerivativeKsi(double ksi,double eta)
// *****************************************************************
// computes the derivative of shape functions w.r.t \ksi 
{
   FloatArray* answer = new FloatArray(3);

	(*answer)[0] = -1. ;
	(*answer)[1] =  1. ; 
	(*answer)[2] =  0. ; 

   return answer ;
}

FloatArray*  FEI2dTriLin :: giveDerivativeEta(double ksi,double eta)
// *****************************************************************
// computes the derivative of shape functions w.r.t \eta
{
   FloatArray* answer = new FloatArray(3);

	(*answer)[0] = -1. ;
	(*answer)[1] =  0. ; 
	(*answer)[2] =  1. ; 

   return answer ;
}

FloatMatrix* FEI2dTriLin :: giveJacobianMatrixAt(Domain* d,IntArray* nodes,Mu::Point* lcoords)
// *******************************************************************************************
// computes the Jacobian matrix.
{
	FloatMatrix *jacobianMatrix = new FloatMatrix(2,2);

	FloatArray *dNdKsi = giveDerivativeKsi();
	FloatArray *dNdEta = giveDerivativeEta();

	for (size_t i = 0 ; i < 3 ; i++)
	{
		double x = d -> giveNode((*nodes)[i])->giveCoordinate(1);
		double y = d -> giveNode((*nodes)[i])->giveCoordinate(2);

		jacobianMatrix->at(1,1) += (*dNdKsi)[i]*x;
		jacobianMatrix->at(1,2) += (*dNdKsi)[i]*y;
		jacobianMatrix->at(2,1) += (*dNdEta)[i]*x;
		jacobianMatrix->at(2,2) += (*dNdEta)[i]*y;
	}

	delete dNdKsi ; delete dNdEta ;

	return jacobianMatrix ;
}

Mu::Point* FEI2dTriLin::local2Global(Domain* d,IntArray* nodes,Mu::Point* lcoords)
// *******************************************************************************
// returns the global coord. of a given point with local coord. (compute the global
// coord. of a Gauss point, for example).
{
	Mu::Point* answer = new Mu::Point();
	FloatArray *N = this->evalN(lcoords);

	for (size_t i = 0 ; i < 3 ; i++)
	{
		double x = d -> giveNode((*nodes)[i])->giveCoordinate(1);
		double y = d -> giveNode((*nodes)[i])->giveCoordinate(2);
		
		answer->x += (*N)[i]*x;
		answer->y += (*N)[i]*y;
	}

	delete N ;

	return answer;
}

Mu::Point* FEI2dTriLin::global2Local(Domain* d,IntArray* nodes,Mu::Point* coords)
// ******************************************************************************
// P : point with global(physical) coordinates x.
// computes its local coordinates (xi,eta).
// From the FE interpolation x = N_{i}(xi,eta) * x_{i}, solve this equa. with 
// variables are xi, eta using Crame method.
{
	
	double x1 = d -> giveNode((*nodes)[0])->giveCoordinate(1);
	double y1 = d -> giveNode((*nodes)[0])->giveCoordinate(2);

	double x2 = d -> giveNode((*nodes)[1])->giveCoordinate(1);
	double y2 = d -> giveNode((*nodes)[1])->giveCoordinate(2);

	double x3 = d -> giveNode((*nodes)[2])->giveCoordinate(1);
	double y3 = d -> giveNode((*nodes)[2])->giveCoordinate(2);

	double x = coords->x ;
	double y = coords->y ;
	
   double A  = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1) ;
	double A1 = (x-x1) *(y3-y1)-(x3-x1)*(y-y1)  ;
	double A2 = (x2-x1)*(y-y1) -(x-x1) *(y2-y1) ;

	Mu::Point* answer = new Mu::Point(A1/A,A2/A);

   // for sure, check if the transformed point are really in the local coord. system
	//assert(answer->x >= 0 && answer->x <= 1 && answer->y >= 0 && answer->y <= 1);

	return answer;
}
