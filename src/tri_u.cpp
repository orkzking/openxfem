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

#include "tri_u.h"


#include "gausspnt.h"
#include "domain.h"     
#include "delaunay.h"     
#include "standardquadrature.h" 
#include "splitgaussquadrature.h" 
#include "integrationrule.h" 
#include "fei2dtrilin.h"
#include "enrichmentitem.h" 
#include "crackinterior.h"
#include "cracktip.h"
#include "geometryentity.h" 
#include "geometry_2D.h" 
#include "piecewiselinear.h"
#include "vertex.h"
#include "functors.h"
#include <vector>
#include <algorithm>


Tri3_U :: Tri3_U (int n, Domain* aDomain) : PlaneProblem(n,aDomain)
{
  numberOfNodes  = 3 ;
}

void  Tri3_U :: computeGaussPoints ()
// *********************************
// Sets up the array containing the Gauss points of the receiver.
// Three Gauss quadrature rules are used
// 1. Classical elements : 3 gauss points are used
// 2. Split and tip elements : 13 Gp for each sub-triangles
// 3. Near tip enriched elements : 13 GP are used (due to asymptotic functions)
{
  numberOfIntegrationRules = 1 ; 
  quadratureRuleArray = new IntegrationRule*;

  // choose the correct integration rule to be used ...

  if (this->isEnriched() == false )  // classical elements
  { 
	 quadratureRuleArray[0] = new StandardGaussLegendreQuadrature;
	 quadratureRuleArray[0]->setUpIntegrationPoints(::TRIANGLE,1);
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
		quadratureRuleArray[0]->setUpIntegrationPoints(::TRIANGLE,13);
	 }
	 else
	 {
		quadratureRuleArray[0] = new StandardGaussLegendreQuadrature;
		quadratureRuleArray[0]->setUpIntegrationPoints(::TRIANGLE,3);
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


FEI2dTriLin* Tri3_U::giveFEInterpolation()
// **************************************
// set up the FE interpolation used for standard FE approximation 
{
  if (!standardFEInterpolation)
	 standardFEInterpolation = new FEI2dTriLin(1,2);

  return static_cast<FEI2dTriLin*>(standardFEInterpolation);
}


FEI2dTriLin* Tri3_U::giveXFEInterpolation()
// ***************************************
// set up the FE interpolation used for enriched approximation 
{
  if (!enrichmentFEInterpolation)
	 enrichmentFEInterpolation = new FEI2dTriLin(1,2);

  return static_cast<FEI2dTriLin*>(enrichmentFEInterpolation);
}

std::vector<DelaunayTriangle *>* Tri3_U::PartitionMySelf()
// ******************************************************
// Partitions split elements and tip elements into sub-triangles
// for numerical integration of the weak form.
// For computational cost, these sub-triangles are already defined in local coord. system
// of the parent three-noded triangle element.
{
  std::vector<Mu::Point>  temp;
  std::vector<Mu::Point>  intersectPoints;
  GeometryEntity* myGeo;

  Mu::Triangle *myTriangle = this->makeTriangle();
  
  // tip-elements interacted with both the tip and  the crack interior,
  // therefore, need to remove crack interior from the list of enrichment items. 
  list<EnrichmentItem*> ::iterator iter1,iter2;
  iter1=find_if(enrichmentItemListOfElem->begin(),enrichmentItemListOfElem->end(),IsType<CrackTip,EnrichmentItem>());
  iter2=find_if(enrichmentItemListOfElem->begin(),enrichmentItemListOfElem->end(),IsType<CrackInterior,EnrichmentItem>());

  if (iter1 != enrichmentItemListOfElem->end() && iter2 != enrichmentItemListOfElem->end()) 
	 enrichmentItemListOfElem->remove(*iter2);
  
  // loop on enrichment items and get the intersection points 
  for(std::list<EnrichmentItem*>::iterator i = enrichmentItemListOfElem->begin() ; i != enrichmentItemListOfElem->end() ; i++)
  {
	 if( (typeid(**i) == typeid(CrackTip)) && (dynamic_cast<CrackTip*>(*i)->amIOnElementEdge()) )
	 {
		intersectPoints.push_back(*this->giveMyCenter());
	 }
	 else
	 {
		myGeo = (*i)->giveMyGeo();
		temp = myGeo->intersection(this);
		if(temp.size() == 1)   // for the case crack touched the element edge
		{
		  Mu::Point *p1 = this->giveNode(1)->makePoint();
		  Mu::Point *p2 = this->giveNode(2)->makePoint();
		  Mu::Point *p3 = this->giveNode(3)->makePoint();
		  Mu::Segment *s1 = new Mu::Segment(*p1,*p2) ;
		  Mu::Segment *s2 = new Mu::Segment(*p2,*p3) ;
		  Mu::Segment *s3 = new Mu::Segment(*p3,*p1) ; 

		  PiecewiseLinear *geo = dynamic_cast<PiecewiseLinear*>(myGeo);
		  std::list<Mu::Point*> *vertices = geo->giveMyListOfVertices();

		  for (std::list<Mu::Point*> ::iterator i = vertices->begin() ; i != vertices->end(); i++)
		  {
			 if( s1->on(*i) || s2->on(*i) || s3->on(*i) )
			 {
				(*i)->print();
				temp.push_back(**i) ;
			 }
		  }
		  delete s1 ; delete s2 ; delete s3 ; 
		  delete p1 ; delete p2 ; delete p3 ; 
		}
		intersectPoints.insert(intersectPoints.end(),temp.begin(),temp.end());
      
      // -----------------------------------------------------------------------------------------------
		// look for kink point. 2005-09-06
		// just segments have kink points, i.e., geometry of CrackInterior only
		// So, use typeid() to test ...
		// Definition of kink point : 
		// 1. A vertex locates within element
		// 2. Is not collinear with intersection points !!! Current code does not check this condition !!!
      
		if( typeid(**i) == typeid(CrackInterior) )
		{
		  std::list<Mu::Point*> *vertices = dynamic_cast<PiecewiseLinear*>(myGeo)->giveMyListOfVertices();
        // attention, just greater than 2 vertices => kink point :) 2005-09-13
		  if(vertices->size() != 2)
		  {
			 for (std::list<Mu::Point*> ::iterator i = vertices->begin() ; i != vertices->end(); i++)
			 {
				if(this->isWithinMe(*i))
				{
				  std::cout << " Oh, found one kink point !!! " << endl ;
				  intersectPoints.push_back(**i) ;
				}
			 }
		  }
		}
		// -----------------------------------------------------------------------------------------------
		
	 }   // end of check CrackTip
  }     // end of loop on enrichment items

  // make the three nodes Mu::Point
  Mu::Point *p0 = new Mu::Point(this->giveNode(1)->giveCoordinate(1),this->giveNode(1)->giveCoordinate(2)) ;
  Mu::Point *p1 = new Mu::Point(this->giveNode(2)->giveCoordinate(1),this->giveNode(2)->giveCoordinate(2)) ;
  Mu::Point *p2 = new Mu::Point(this->giveNode(3)->giveCoordinate(1),this->giveNode(3)->giveCoordinate(2)) ;

  // transform to local coordinate of the parent element
  Mu::Point *localP0 = this->giveFEInterpolation()->global2Local(domain,this->giveNodeArray(),p0);
  Mu::Point *localP1 = this->giveFEInterpolation()->global2Local(domain,this->giveNodeArray(),p1);
  Mu::Point *localP2 = this->giveFEInterpolation()->global2Local(domain,this->giveNodeArray(),p2);

  DelaunayTree *dt = new DelaunayTree(localP0,localP1,localP2) ;
  
  // insert the intersection points + kink points into the DT
  for(size_t j = 0 ; j < intersectPoints.size() ; j++)
  {
	 dt->insert(standardFEInterpolation->global2Local(domain,this->giveNodeArray(),&intersectPoints[j]));
  }

  // do the Delaunay triangulation ...
  vector<DelaunayTriangle *> *tri = dt->getTriangles(); 
  
  delete myTriangle ; //delete dt ;
  delete p0 ; delete p1 ; delete p2 ;

  return tri;
}

double Tri3_U :: computeAreaAboveEnrItem(EnrichmentItem *enrItem)
// *************************************************************
// The enrItem splits Tri_U into two portions : compute the area of the
// portion above the enrItem.
// Rewritten version : 2005-09-20.
{
  GeometryEntity* myGeo = enrItem->giveMyGeo() ;
  Mu::Triangle *tri = this->makeTriangle();
  std::vector<Mu::Point>  intersect = myGeo->intersection(tri); 

  // Attention for the case the crack goes throught the node
  // these case because of numerical stability, Segment::intersection() returns FALSE !!!
  // Therefore, need to use Segment::on(Point* p) to check.
  PiecewiseLinear *geo = static_cast<PiecewiseLinear*>(enrItem->giveMyGeo());
  Mu::Point *p1 = this->giveNode(1)->makePoint();
  Mu::Point *p2 = this->giveNode(2)->makePoint();
  Mu::Point *p3 = this->giveNode(3)->makePoint();
  if (intersect.size() == 1)
  {
	 std::vector<Mu::Segment*> *seg =  geo->giveSegments();
	 for(size_t i = 0 ; i < seg->size() ; i++)
	 {
		if( (*seg)[i]->on(p1))
		  intersect.push_back(*p1) ;
		if( (*seg)[i]->on(p2))
		  intersect.push_back(*p2) ;
		if( (*seg)[i]->on(p3))
		  intersect.push_back(*p3) ;
	 }
  }
  // just do for the case enrItem completely splits the element => always have 2 intersection points
  if (intersect.size() != 2)
  {
	 std::cout << "enrItem must completely splits the element !!! " << std::endl;
	 assert(false);
  }
  
  DelaunayTree *dt = new DelaunayTree(p1,p2,p3);
  // insert intersection points into DelaunayTree
  dt->insert(&intersect[0]) ;
  dt->insert(&intersect[1]) ;

  // looking for kink points !!! 2005-09-18
  std::list<Mu::Point*> *vertices = geo->giveMyListOfVertices();
  if(vertices->size() != 2)
  {
	 for(std::list<Mu::Point*>::iterator i = vertices->begin() ; i != vertices->end(); i++)
	 {
		if(this->isWithinMe(*i) && (isAligned(*i,&intersect[0],&intersect[1]) == false))
		  dt->insert(*i) ;
	 }
  }

  double x1 = intersect[0].x ; // first intersection point
  double y1 = intersect[0].y ; // first intersection point

  double x2 = intersect[1].x ; // second intersection point
  double y2 = intersect[1].y ; // second intersection point

  vector<DelaunayTriangle *> *v  = dt->getTriangles(); 
  double A = 0.0  ; double delta ;
  // loop on delaunay triangles and compute the accumulated area
  for(size_t i = 0 ; i < v->size() ; i++)
  {
	 if((*v)[i]->isTriangle == true)
	 {
		Mu::Point *center = (*v)[i]->getCenter() ;
		if(x2 > x1)
		  delta = (x1-center->x)*(y2-center->y)-(x2-center->x)*(y1-center->y);
		else
		  delta = (x2-center->x)*(y1-center->y)-(x1-center->x)*(y2-center->y);

		if(delta > 0.000001) 
		  A  += (*v)[i]->area(); 
	 }
  }
  delete p1 ; delete p2 ; delete p3 ; 
  delete tri ; delete dt ; delete v ;  // Purify,14-10-2005

  return A ;
}

double Tri3_U :: area() 
//********************
{
  double x0 = this->giveNode(1)->giveCoordinate(1);
  double y0 = this->giveNode(1)->giveCoordinate(2);

  double x1 = this->giveNode(2)->giveCoordinate(1);
  double y1 = this->giveNode(2)->giveCoordinate(2);

  double x2 = this->giveNode(3)->giveCoordinate(1);
  double y2 = this->giveNode(3)->giveCoordinate(2);

  double A = 0.5 * ((x0-x2)*(y1-y2) - (x1-x2)*(y0-y2)) ; // signed area 

  if (A < 0) A = -A ;

  return A ;
}

Mu::Triangle*  Tri3_U::makeTriangle()
//**********************************
{
  double x1 = this->giveNode(1)->giveCoordinate(1);
  double y1 = this->giveNode(1)->giveCoordinate(2);

  double x2 = this->giveNode(2)->giveCoordinate(1);
  double y2 = this->giveNode(2)->giveCoordinate(2);

  double x3 = this->giveNode(3)->giveCoordinate(1);
  double y3 = this->giveNode(3)->giveCoordinate(2);

  Mu::Point *p1 = new Mu::Point(x1,y1) ;
  Mu::Point *p2 = new Mu::Point(x2,y2) ;
  Mu::Point *p3 = new Mu::Point(x3,y3) ;

  Mu::Triangle *tri = new Mu::Triangle(p1,p2,p3);

  //delete p1,p2,p3;

  return tri ;
}

GaussPoint**  Tri3_U :: setGaussQuadForJ_Integral()
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


