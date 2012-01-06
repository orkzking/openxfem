

//   file PLATEISO4.CPP
 
#include "plateiso4.h"

#include "node.h"
#include "material.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "diagmtrx.h"
#include "fei2dquadlin.h"
#include "intarray.h"
#include "domain.h"     
#include "standardquadrature.h" 
#include "splitgaussquadrature.h"
#include "integrationrule.h" 
#include "enrichmentitem.h"
#include "piecewiselinear.h"
#include "vertex.h"
#include "geometry_2D.h"
#include "geometryentity.h" 
#include "functors.h"
#include "crackinterior.h"
#include "cracktip.h"
#include <math.h>
#include <vector>



 
PlateIsoQ4 :: PlateIsoQ4 (int n, Domain* aDomain)
      : PlateIso(n,aDomain)
{
   numberOfNodes  = 4 ;	
	this->computeGaussPointsS();
}

void  PlateIsoQ4 :: computeGaussPoints()
// ************************************
// Sets up the array containing the Gauss points of the receiver.
// Modified to use class Integration Rule and discontinuous Gauss quadrature.
{
  numberOfIntegrationRules = 1 ;
  quadratureRuleArray = new IntegrationRule*;
  int numOfGPs = 16 ;
  if (this->isEnriched() == false )  // classical elements
  { 
	 quadratureRuleArray[0] = new StandardGaussLegendreQuadrature;
	 quadratureRuleArray[0]->setUpIntegrationPoints(SQUARE,4);
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
		quadratureRuleArray[0]->setUpIntegrationPoints(SQUARE,numOfGPs);
	 }
	 else           // step enriched elements
	 {
		quadratureRuleArray[0] = new StandardGaussLegendreQuadrature;
		quadratureRuleArray[0]->setUpIntegrationPoints(SQUARE,4);
	 }
  }

  // tranform the obtained integration points into the GaussPoint array of element

  vector<double>* weight = quadratureRuleArray[0]->giveWeightArray();
  vector<Mu::Point*>* coordArray = quadratureRuleArray[0]->giveIntegrationPointVector();	

  numberOfGaussPoints = coordArray->size();

  gaussPointArray  = new GaussPoint* [numberOfGaussPoints] ;

  for(size_t i = 0 ; i < numberOfGaussPoints ; i++)
  {
	 gaussPointArray[i] = new GaussPoint(this,i+1,coordArray->at(i),weight->at(i),4) ;
  }
  delete weight ; delete coordArray ;   // Purify, 08-11-05
}

void  PlateIsoQ4 :: computeGaussPointsS ()
// **************************************
{  
   gaussPointArrayS    = new GaussPoint* [1] ;
	Mu::Point *p = new Mu::Point(0,0);
   gaussPointArrayS[0] = new GaussPoint(this,1,p,4.,4) ;  
 }

GaussPoint**  PlateIsoQ4 :: setGaussQuadForJ_Integral()
// ************************************************
// Sets up the array containing the Gauss points used for the J integral
// computation.
{
  
  IntegrationRule *GaussQuad = new IntegrationRule ;
  size_t numOfGPs1 = 9 ;
  size_t numOfGPs2 = 7 ; // number of GPs for each sub-triangle
  
  // choose the correct integration rule to be used ...

  if (!enrichmentItemListOfElem) // no enrichment 
  { 
	 GaussQuad = new StandardGaussLegendreQuadrature;
	 GaussQuad->setUpIntegrationPoints(SQUARE,numOfGPs1); // 9,25,36,49,56,
  }
  else
  {
	 GaussQuad = new SplitGaussLegendreQuadrature;
	 GaussQuad->setUpIntegrationPoints(this,numOfGPs2); // 7 or 13
  } 

  vector<double>* weight;
  vector<Mu::Point*>* coordArray;	

  weight = GaussQuad->giveWeightArray();
  coordArray = GaussQuad->giveIntegrationPointVector();

  numOfGPsForJ_Integral = coordArray->size() ;
  GaussPoint** ret = new GaussPoint* [numOfGPsForJ_Integral] ;
  
  for(size_t i = 0 ; i < weight->size() ; i++)
	 ret[i] = new GaussPoint(this,i+1,coordArray->at(i),weight->at(i),4) ;

  delete weight ; delete coordArray ;   // Purify, 14-10-05

  return ret ;
}

FEI2dQuadLin* PlateIsoQ4::giveFEInterpolation()
{

  if (!standardFEInterpolation)
	 standardFEInterpolation = new FEI2dQuadLin(1,2);

  return static_cast<FEI2dQuadLin*>(standardFEInterpolation);

}

 
FEI2dQuadLin* PlateIsoQ4::giveXFEInterpolation()
{
  if (!enrichmentFEInterpolation)
	 enrichmentFEInterpolation = new FEI2dQuadLin(1,2);

  return static_cast<FEI2dQuadLin*>(enrichmentFEInterpolation);
}

double PlateIsoQ4 :: area()
{
  double x0 = this->giveNode(1)->giveCoordinate(1);
  double y0 = this->giveNode(1)->giveCoordinate(2);

  double x1 = this->giveNode(2)->giveCoordinate(1);
  double y1 = this->giveNode(2)->giveCoordinate(2);

  double x2 = this->giveNode(3)->giveCoordinate(1);
  double y2 = this->giveNode(3)->giveCoordinate(2);

  double x3 = this->giveNode(4)->giveCoordinate(1);
  double y3 = this->giveNode(4)->giveCoordinate(2);

  double A1 = 0.5 * ((x0-x2)*(y1-y2) - (x1-x2)*(y0-y2)) ; // signed area 
  double A2 = 0.5 * ((x0-x3)*(y2-y3) - (x2-x3)*(y0-y3)) ; // signed area 

  return A1 + A2 ;
}

