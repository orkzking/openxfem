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


#ifndef _INTEGRATIONRULE_H_
#define _INTEGRATIONRULE_H_


#include "geometry_base.h"
#include "element.h"
#include <vector>


//! Finite element integration rules
/*!
This class implements integration rule class.\n
DESCRIPTION: \n
Stores integration points used for integration
of necesary terms (for example computation of  stiffness matrix 
or computation of element nodal force vector ) 
and it  corresponds to some local strains 
on finite element level. Finite element can have many 
integration rules corresponding to  different strains.\n
*/

typedef enum{LINE,TRIANGLE,SQUARE,TETRAHEDRA
} IntDomain;

class IntegrationRule
{

public:

	IntegrationRule ();         //<! Constructor   
	virtual ~IntegrationRule(); //<! Destructor
	
	//!  Returns number of integration points of receiver.
	size_t giveNumberOfIntegrationPoints () {return numberOfIntegrationPoints;}
	/*!
	Initializes the receiver. Receiver integration points are created acording to given parameters.
	@param dType describes integration domain
	@param nPoints required number of integration points of receiver
	*/
	virtual void setUpIntegrationPoints (IntDomain dType, size_t nPoints);
   /*!
	 Set up the integration points for elements cut by the discontinuities
	 Will be defined in only one derived class : SplitGaussLegendreQuadrature
	*/
	virtual void setUpIntegrationPoints (Element* elem,size_t nPoints){}
	/*!
	Abstract service.
	Returns requred number of integration points to exactly integrate
	polynomial of order approxOrder on given domain.
	When approxOrder is too large and is not supported by implementation
	method returns -1. Must be overloaded by derived classes.
	*/
	virtual int   getRequiredNumberOfIntegrationPoints (IntDomain dType, int approxOrder){return 0;}
	std::vector<double>*      giveWeightArray();
	std::vector<Mu::Point*>*  giveIntegrationPointVector();

protected:

	size_t  numberOfIntegrationPoints;
	std::vector<Mu::Point*>*   integrationPointVector; //!< vector of integration points  
	std::vector<double>*       weightArray ; //!< vector of weights
   /*!
    Sets up receiver's  integration points on unit line integration domain.
    Default implementaion does not sets up any integration points and returns 0.
    Must be overloaded by deived classes.
    @returns the array of Gauss points
	*/
	virtual void  SetUpPointsOnLine       (size_t) {}
	/*! 
	Sets up receiver's  integration points on triangular (area coords) integration domain.
	Default implementaion does not sets up any integration points and returns 0.
	Must be overloaded by deived classes.
	@returns the array of Gauss points
	*/
	virtual void  SetUpPointsOnTriangle    (size_t) {}
	/*! 
	Sets up receiver's  integration points on unit square integration domain.
	Default implementaion does not sets up any integration points and returns 0.
	Must be overloaded by deived classes.
	@returns the array of Gauss points
	*/
	virtual void  SetUpPointsOnSquare     (size_t) {}
	/*! 
	Sets up receiver's  integration points on tetrahedra domain.
	Default implementaion does not sets up any integration points and returns 0.
	Must be overloaded by deived classes.
	@returns the array of Gauss points
	*/
	virtual void  SetUpPointsOnTetrahedra (size_t) {}
	
};


#endif //_INTEGRATIONRULE_H_

