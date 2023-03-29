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

#include "fei3dlineartetra.h"

#include "flotmtrx.h"
#include "flotarry.h"
#include "geometry_base.h"


FloatArray*  FEI3dLinearTetrahedra::evalN(Mu::Point* p)
// ****************************************************
// computes the value of shape functions at a given point p
// Result:  N=[1-xi-eta-zeta ; xi ; eta ; zeta]
{    
	double n1  =  1.0 - p->x - p->y - p->z;
	double n2  =  p->x ;
	double n3  =  p->y ;
	double n4  =  p->z ;

	FloatArray* answer = new FloatArray(3) ;

	(*answer)[0] = n1 ; 
	(*answer)[1] = n2 ;
	(*answer)[2] = n3 ;
	(*answer)[3] = n4 ;

   return answer ;
}

FloatMatrix* FEI3dLinearTetrahedra :: evaldNdx(Domain* d,IntArray* Nodes,Mu::Point* lcoords)
// ********************************************************************************
// computes the derivatives of shape functions w.r.t global coordinates
// dNdx = [N1,x N1,y N1,z
//         N2,x N2,y N2,z
//         N3,x N3,y N3,z
//         N4,x N4,y N4,z]
{
	FloatMatrix *jacobianMatrix = giveJacobianMatrixAt(d,Nodes,lcoords);
	
	FloatMatrix *inv = jacobianMatrix->GiveInverse();

	FloatArray *dNdKsi  = giveDerivativeKsi();
	FloatArray *dNdEta  = giveDerivativeEta();
	FloatArray *dNdZeta = giveDerivativeZeta();

	FloatMatrix * answer = new FloatMatrix(4,3);

	for (size_t i = 0 ; i < 4 ; i++)
	{
		answer->at(i+1,1) = (*dNdKsi)[i]*inv->at(1,1) + (*dNdEta)[i]*inv->at(1,2) + (*dNdZeta)[i]*inv->at(1,3);
		answer->at(i+1,2) = (*dNdKsi)[i]*inv->at(2,1) + (*dNdEta)[i]*inv->at(2,2) + (*dNdZeta)[i]*inv->at(2,3);
		answer->at(i+1,3) = (*dNdKsi)[i]*inv->at(3,1) + (*dNdEta)[i]*inv->at(3,2) + (*dNdZeta)[i]*inv->at(3,3);
		
	}

	delete jacobianMatrix ;
	delete inv ;
	delete dNdKsi ; // Purify, 11-10-05
	delete dNdEta ; // Purify, 11-10-05
	delete dNdZeta ; // Purify, 11-10-05

	return answer;
}

FloatArray*  FEI3dLinearTetrahedra :: giveDerivativeKsi()
// ******************************************************
// computes the derivative of shape functions w.r.t \ksi 
{
   FloatArray* answer = new FloatArray(4);

	(*answer)[0] = -1. ;
	(*answer)[1] =  1. ; 
	(*answer)[2] =  0. ; 
	(*answer)[3] =  0. ; 

   return answer ;
}

FloatArray*  FEI3dLinearTetrahedra :: giveDerivativeEta()
// ******************************************************
// computes the derivative of shape functions w.r.t \eta
{
   FloatArray* answer = new FloatArray(4);

	(*answer)[0] = -1. ;
	(*answer)[1] =  0. ; 
	(*answer)[2] =  1. ; 
	(*answer)[3] =  0. ; 

   return answer ;
}

FloatArray*  FEI3dLinearTetrahedra :: giveDerivativeZeta()
// ******************************************************
// computes the derivative of shape functions w.r.t \zeta
{
   FloatArray* answer = new FloatArray(4);

	(*answer)[0] = -1. ;
	(*answer)[1] =  0. ; 
	(*answer)[2] =  0. ; 
	(*answer)[3] =  1. ; 

   return answer ;
}

FloatMatrix* FEI3dLinearTetrahedra :: giveJacobianMatrixAt(Domain* d,IntArray* nodes,Mu::Point* lcoords)
// *******************************************************************************************
// computes the Jacobian matrix [3x3].
{
	FloatMatrix *jacobianMatrix = new FloatMatrix(3,3);

	FloatArray *dNdKsi  = giveDerivativeKsi();
	FloatArray *dNdEta  = giveDerivativeEta();
	FloatArray *dNdZeta = giveDerivativeZeta();

	for (size_t i = 0 ; i < 4 ; i++)
	{
		double x = d -> giveNode((*nodes)[i])->giveCoordinate(1);
		double y = d -> giveNode((*nodes)[i])->giveCoordinate(2);
		double z = d -> giveNode((*nodes)[i])->giveCoordinate(3);

		jacobianMatrix->at(1,1) += (*dNdKsi)[i]*x;
		jacobianMatrix->at(1,2) += (*dNdKsi)[i]*y;
		jacobianMatrix->at(1,3) += (*dNdKsi)[i]*z;

		jacobianMatrix->at(2,1) += (*dNdEta)[i]*x;
		jacobianMatrix->at(2,2) += (*dNdEta)[i]*y;
		jacobianMatrix->at(2,3) += (*dNdEta)[i]*z;

		jacobianMatrix->at(3,1) += (*dNdZeta)[i]*x;
		jacobianMatrix->at(3,2) += (*dNdZeta)[i]*y;
		jacobianMatrix->at(3,3) += (*dNdZeta)[i]*z;
	}

	delete dNdKsi ; delete dNdEta ; delete dNdZeta;

	return jacobianMatrix ;
}

Mu::Point* FEI3dLinearTetrahedra::local2Global(Domain* d,IntArray* nodes,Mu::Point* lcoords)
// *******************************************************************************
// returns the global coord. of a given point with local coord. (compute the global
// coord. of a Gauss point, for example).
{
	Mu::Point* answer = new Mu::Point();
	FloatArray *N = this->evalN(lcoords);

	for (size_t i = 0 ; i < 4 ; i++)
	{
		double x = d -> giveNode((*nodes)[i])->giveCoordinate(1);
		double y = d -> giveNode((*nodes)[i])->giveCoordinate(2);
		double z = d -> giveNode((*nodes)[i])->giveCoordinate(3);
		
		answer->x += (*N)[i]*x;
		answer->y += (*N)[i]*y;
		answer->z += (*N)[i]*z;
	}

	delete N ;

	return answer;
}

Mu::Point* FEI3dLinearTetrahedra::global2Local(Domain* d,IntArray* nodes,Mu::Point* coords)
// ******************************************************************************
// P : point with global(physical) coordinates x.
// computes its local coordinates (xi,eta).
// From the FE interpolation x = N_{i}(xi,eta) * x_{i}, solve this equa. with 
// variables are xi, eta using Crame method.
{
	
	double x1 = d -> giveNode((*nodes)[0])->giveCoordinate(1);
	double y1 = d -> giveNode((*nodes)[0])->giveCoordinate(2);
	double z1 = d -> giveNode((*nodes)[0])->giveCoordinate(3);

	double x2 = d -> giveNode((*nodes)[1])->giveCoordinate(1);
	double y2 = d -> giveNode((*nodes)[1])->giveCoordinate(2);
	double z2 = d -> giveNode((*nodes)[1])->giveCoordinate(3);

	double x3 = d -> giveNode((*nodes)[2])->giveCoordinate(1);
	double y3 = d -> giveNode((*nodes)[2])->giveCoordinate(2);
	double z3 = d -> giveNode((*nodes)[2])->giveCoordinate(3);

   double x4 = d -> giveNode((*nodes)[3])->giveCoordinate(1);
	double y4 = d -> giveNode((*nodes)[3])->giveCoordinate(2);
	double z4 = d -> giveNode((*nodes)[3])->giveCoordinate(3);

	double x = coords->x ;
	double y = coords->y ;
	double z = coords->z ;
	
   FloatMatrix *a = new FloatMatrix(3,3);
	a->at(1,1) = x2 - x1 ; a->at(1,2) = x3 - x1 ; a->at(1,3) = x4 - x1 ;
	a->at(2,1) = y2 - y1 ; a->at(2,2) = y3 - y1 ; a->at(2,3) = y4 - y1 ;
	a->at(3,1) = z2 - z1 ; a->at(3,2) = z3 - z1 ; a->at(3,3) = z4 - z1 ;

	FloatMatrix *a1 = new FloatMatrix(3,3);
	a1->at(1,1) = x - x1 ; a1->at(1,2) = x3 - x1 ; a1->at(1,3) = x4 - x1 ;
	a1->at(2,1) = y - y1 ; a1->at(2,2) = y3 - y1 ; a1->at(2,3) = y4 - y1 ;
	a1->at(3,1) = z - z1 ; a1->at(3,2) = z3 - z1 ; a1->at(3,3) = z4 - z1 ;

	FloatMatrix *a2 = new FloatMatrix(3,3);
	a2->at(1,1) = x2 - x1 ; a2->at(1,2) = x - x1 ; a2->at(1,3) = x4 - x1 ;
	a2->at(2,1) = y2 - y1 ; a2->at(2,2) = y - y1 ; a2->at(2,3) = y4 - y1 ;
	a2->at(3,1) = z2 - z1 ; a2->at(3,2) = z - z1 ; a2->at(3,3) = z4 - z1 ;

	FloatMatrix *a3 = new FloatMatrix(3,3);
	a3->at(1,1) = x2 - x1 ; a3->at(1,2) = x3 - x1 ; a3->at(1,3) = x - x1 ;
	a3->at(2,1) = y2 - y1 ; a3->at(2,2) = y3 - y1 ; a3->at(2,3) = y - y1 ;
	a3->at(3,1) = z2 - z1 ; a3->at(3,2) = z3 - z1 ; a3->at(3,3) = z - z1 ;

	double A  = a->giveDeterminant();
   
	double A1 = a1->giveDeterminant();
	double A2 = a2->giveDeterminant();
	double A3 = a3->giveDeterminant();

	Mu::Point* answer = new Mu::Point(A1/A,A2/A,A3/A);

   // for sure, check if the transformed point are really in the local coord. system
	//assert(answer->x >= 0 && answer->x <= 1 && answer->y >= 0 && answer->y <= 1);

	return answer;
}


