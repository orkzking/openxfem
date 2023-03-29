// FILE FEI2DTRIQUA.CPP

#include "fei2dtriquad.h"

#include "flotmtrx.h"
#include "flotarry.h"
#include "polynoxy.h"

FloatArray*  FEI2dTriQuad :: evalN(Mu::Point* p)
// *********************************************
{
	FloatArray* answer = new FloatArray(6) ;

	double ksi = p->x ;
	double eta = p->y ;

	(*answer)[0] = 1 - 3*(ksi + eta) + 4*ksi*eta + 2*(ksi*ksi + eta*eta);
	(*answer)[1] = ksi * (2*ksi-1);
	(*answer)[2] = eta * (2*eta-1);
	(*answer)[3] = 4* ksi * (1-ksi-eta);
	(*answer)[4] = 4* ksi * eta;
	(*answer)[5] = 4* eta * (1-ksi-eta);
   
   return answer ;
}

FloatMatrix* FEI2dTriQuad :: evaldNdx(Domain* d,IntArray* Nodes,Mu::Point* p)
// *************************************************************************
// computes the derivatives of shape functions w.r.t global coordinates
// dNdx = [N1,x N1,y
//         N2,x N2,y
//         ...  ...
//         N6,x N6,y]
{

	FloatMatrix *jacobianMatrix= giveJacobianMatrixAt(d,Nodes,p);
	FloatMatrix *inv = jacobianMatrix->GiveInverse();

	FloatArray *dNdKsi = giveDerivativeKsi(p->x,p->y);
	FloatArray *dNdEta = giveDerivativeEta(p->x,p->y);

	FloatMatrix *answer = new FloatMatrix(6,2);

	for (size_t i = 0 ; i < 6 ; i++)
	{
		answer->at(i+1,1) = (*dNdKsi)[i]*inv->at(1,1) + (*dNdEta)[i]*inv->at(1,2);
		answer->at(i+1,2) = (*dNdKsi)[i]*inv->at(2,1) + (*dNdEta)[i]*inv->at(2,2);
	}

	delete jacobianMatrix ;
	delete inv ;

	return answer;
}

FloatArray*  FEI2dTriQuad :: giveDerivativeKsi (double eta,double ksi)
// ******************************************************************
// computes the derivative of shape functions w.r.t \ksi 
{
	FloatArray* answer = new FloatArray(6);

	(*answer)[0] =  4* (ksi+eta)-3;
	(*answer)[1] =  4* ksi - 1 ; 
	(*answer)[2] =  0 ; 
	(*answer)[3] =  4* (1-eta-2*ksi) ; 
	(*answer)[4] =  4* eta ; 
	(*answer)[5] = -4* eta ; 

   return answer ;
}

FloatArray*  FEI2dTriQuad :: giveDerivativeEta (double eta, double ksi)
// ******************************************************************
// computes the derivative of shape functions w.r.t \eta
{
   FloatArray* answer = new FloatArray(6);

	(*answer)[0] =  4* (ksi+eta)-3;
	(*answer)[1] =  0 ; 
	(*answer)[2] =  4* eta - 1 ; 
	(*answer)[3] = -4* ksi ; 
	(*answer)[4] =  4* ksi ; 
	(*answer)[5] =  4*(1-ksi-2*eta) ; 

   return answer ;
}

FloatMatrix* FEI2dTriQuad :: giveJacobianMatrixAt (Domain* d,IntArray* nodes, Mu::Point* lcoords)
// *********************************************************************************************
// computes the Jacobian matrix.
{
 
	FloatArray * dNdKsi = giveDerivativeKsi(lcoords->x,lcoords->y);
	FloatArray * dNdEta = giveDerivativeEta(lcoords->x,lcoords->y);

	FloatMatrix *jacobianMatrix = new FloatMatrix(2,2);

	for (size_t i = 0 ; i < 6 ; i++)
	{
		double x = d -> giveNode((*nodes)[i]) -> giveCoordinate(1);
		double y = d -> giveNode((*nodes)[i]) -> giveCoordinate(2);

		jacobianMatrix->at(1,1) += (*dNdKsi)[i]*x;
		jacobianMatrix->at(1,2) += (*dNdKsi)[i]*y;
		jacobianMatrix->at(2,1) += (*dNdEta)[i]*x;
		jacobianMatrix->at(2,2) += (*dNdEta)[i]*y;
	}

	delete dNdKsi ; delete dNdEta ;
	
	return jacobianMatrix ;
}

Mu::Point* FEI2dTriQuad :: local2Global(Domain* d,IntArray* nodes,Mu::Point* lcoords)
// *******************************************************************************
// returns the global coord. of a given point with local coord. (compute the global
// coord. of a Gauss point, for example).
{
	Mu::Point* answer = new Mu::Point();

	for (size_t i = 0 ; i < 6 ; i++)
	{
		double x = d -> giveNode((*nodes)[i])->giveCoordinate(1);
		double y = d -> giveNode((*nodes)[i])->giveCoordinate(2);
		
		answer->x += (*this->evalN(lcoords))[i]*x;
		answer->x += (*this->evalN(lcoords))[i]*y;
	}
	return answer;
}

Mu::Point* FEI2dTriQuad :: global2Local(Domain* d,IntArray* nodes,Mu::Point* lcoords)
// **********************************************************************************
// Implemented later ...... need to use Newton Raphson method
{
	Mu::Point* answer = new Mu::Point(0,0);

	return answer;
}