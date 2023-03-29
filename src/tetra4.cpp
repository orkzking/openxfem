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

#include "tetra4.h"


#include "gausspnt.h"
#include "domain.h"     
#include "delaunay.h"     
#include "standardquadrature.h" 
#include "splitgaussquadrature.h" 
#include "integrationrule.h" 
#include "fei3dlineartetra.h"
#include "enrichmentitem.h" 
#include "crackinterior.h"
#include "cracktip.h"
#include "geometryentity.h" 
#include "geometry_2D.h" 
#include "functors.h"
#include <vector>
#include <algorithm>


Tetra4 :: Tetra4 (int n, Domain* aDomain) : Iso3dElement(n,aDomain)
{
  numberOfNodes  = 4 ;
}

void  Tetra4 :: computeGaussPoints ()
// *********************************
// Sets up the array containing the Gauss points of the receiver.
// Three Gauss quadrature rules are used
// 1. Classical elements : 1 gauss points are used
// 2. Split and tip elements : 13 Gp for each sub-triangles
// 3. Near tip enriched elements : 13 GP are used (due to asymptotic functions)
// REMARK: \TODO for 3D X-FEM implementation !!!
{
  numberOfIntegrationRules = 1 ; 
  quadratureRuleArray = new IntegrationRule*;

  // choose the correct integration rule to be used ...

  if (this->isEnriched() == false )  // classical elements
  { 
	 quadratureRuleArray[0] = new StandardGaussLegendreQuadrature;
	 quadratureRuleArray[0]->setUpIntegrationPoints(::TETRAHEDRA,1);
  }
  else if (enrichmentItemListOfElem != NULL) // split, tip elements
  {
	 quadratureRuleArray[0] = new SplitGaussLegendreQuadrature;
	 quadratureRuleArray[0]->setUpIntegrationPoints(this,13);
  }
  else  // partially enriched elements 
  {
	 size_t count = 0 ;
	 for(size_t i = 0 ; i < numberOfNodes ; i++)
	 {
		if(this->giveNode(i+1)->isTipEnriched())
		  count += 1 ;
	 }
	 if(count != 0) // near tip enriched elements
	 {
		quadratureRuleArray[0] = new StandardGaussLegendreQuadrature;
		quadratureRuleArray[0]->setUpIntegrationPoints(::TETRAHEDRA,13);
	 }
	 else
	 {
		quadratureRuleArray[0] = new StandardGaussLegendreQuadrature;
		quadratureRuleArray[0]->setUpIntegrationPoints(::TETRAHEDRA,3);
	 }
  }

  // tranform the obtained integration points into the GaussPoint array of element

  vector<double>* weight = quadratureRuleArray[0]->giveWeightArray();
  vector<Mu::Point*>* coordArray = quadratureRuleArray[0]->giveIntegrationPointVector();	

  if(gaussPointArray)
  {
	 for (size_t i = 0 ; i < numberOfGaussPoints ; i++)
		delete gaussPointArray[i] ;
    delete gaussPointArray ;
  }

  numberOfGaussPoints = coordArray->size();

  gaussPointArray  = new GaussPoint* [numberOfGaussPoints] ;

  for(size_t i = 0 ; i < numberOfGaussPoints ; i++)
  {
	 gaussPointArray[i] = new GaussPoint(this,i+1,coordArray->at(i),weight->at(i),4) ;
  }
  delete weight ; delete coordArray ;   // Purify, 08-11-05
}


FEI3dLinearTetrahedra* Tetra4::giveFEInterpolation()
// *************************************************
// set up the FE interpolation used for standard FE approximation 
{
  if (!standardFEInterpolation)
	 standardFEInterpolation = new FEI3dLinearTetrahedra();

  return static_cast<FEI3dLinearTetrahedra*>(standardFEInterpolation);
}


FEI3dLinearTetrahedra* Tetra4::giveXFEInterpolation()
// **************************************************
// set up the FE interpolation used for enriched approximation 
{
  if (!enrichmentFEInterpolation)
	 enrichmentFEInterpolation = new FEI3dLinearTetrahedra;

  return static_cast<FEI3dLinearTetrahedra*>(enrichmentFEInterpolation);
}


double Tetra4 :: volume() 
//*********************
// NOT YET IMPPLEMENTED
{
  return 0. ;
}

GaussPoint**  Tetra4 :: setGaussQuadForJ_Integral()
// ***********************************************
// Two Gauss quadrature rules are used
// 1. Classical elements : 7 Gauss points are used
// 2. Split and tip elements : 13 Gps for each sub-triangles

{
  IntegrationRule *GaussQuad = new IntegrationRule ;

  if (!enrichmentItemListOfElem) // no need discontinuous Gauss quad 
  { 
	 GaussQuad = new StandardGaussLegendreQuadrature;
	 GaussQuad->setUpIntegrationPoints(::TRIANGLE,7);
  }
  else
  {
	 GaussQuad = new SplitGaussLegendreQuadrature;
	 GaussQuad->setUpIntegrationPoints(this,7);
  } 

  // tranform the obtained integration points into the GaussPoint array of element

  vector<double>     *weight     = GaussQuad->giveWeightArray();
  vector<Mu::Point*> *coordArray = GaussQuad->giveIntegrationPointVector();	

  numOfGPsForJ_Integral = coordArray->size() ;

  GaussPoint** ret  = new GaussPoint* [numOfGPsForJ_Integral] ;

  for(size_t i = 0 ; i < numOfGPsForJ_Integral ; i++)
  {
	 ret[i] = new GaussPoint(this,i+1,coordArray->at(i),weight->at(i),4) ;
  }

  delete weight ; delete coordArray ;  // Purify, 14-10-05

  return ret ;
}

std::vector<DelaunayTriangle *>* Tetra4 :: PartitionMySelf()
// *********************************************************
// Partitions split elements and tip elements into sub-triangles
// for numerical integration of the weak form.
// For computational cost, these sub-triangles are already defined in local coord. system
// of the parent three-noded triangle element.
{
  Mu::Point *P0 = new Mu::Point(0.,0.,0.); 
  Mu::Point *P1 = new Mu::Point(0.,0.,0.); 
  Mu::Point *P2 = new Mu::Point(0.,0.,0.); 
  DelaunayTree *dt = new DelaunayTree(P0,P1,P2) ;
  vector<DelaunayTriangle *> *tri = dt->getTriangles(); 
  return tri ;
}

double Tetra4 :: computeAreaAboveEnrItem(EnrichmentItem *enrItem)
// *************************************************************
// The enrItem splits Tri_U into two portions : compute the area of the
// portion above the enrItem.
// Rewritten version : 2005-09-20.
{
  return 0.0 ;
}