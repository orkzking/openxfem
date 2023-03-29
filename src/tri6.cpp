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

#include "tri6.h"

#include "material.h"
#include "gausspnt.h"
#include "flotarry.h"
#include "planelast.h"
#include "domain.h"     
#include "delaunay.h"     
#include "standardquadrature.h" 
#include "splitgaussquadrature.h" 
#include "integrationrule.h" 
#include "fei2dtrilin.h"
#include "fei2dtriquad.h"
#include "enrichmentitem.h" 
#include "crackinterior.h"
#include "cracktip.h"
#include "geometryentity.h" 
#include "geometry_2D.h" 
#include "piecewiselinear.h"
#include "vertex.h"
#include "functors.h"
#include <math.h>
#include <vector>
#include <algorithm>


Tri6_U :: Tri6_U (int n, Domain* aDomain) : PlaneProblem(n,aDomain)
{
  numberOfNodes  = 6 ;  
}

void  Tri6_U :: computeGaussPoints ()
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
	 quadratureRuleArray[0]->setUpIntegrationPoints(::TRIANGLE,7);
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
		quadratureRuleArray[0]->setUpIntegrationPoints(::TRIANGLE,7);
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
}


FEI2dTriQuad* Tri6_U::giveFEInterpolation()
// ****************************************
// set up the FE interpolation used for standard FE approximation 
{
  if (!standardFEInterpolation)
	 standardFEInterpolation = new FEI2dTriQuad(1,2);

  return static_cast<FEI2dTriQuad*>(standardFEInterpolation);
}


FEI2dTriLin* Tri6_U::giveXFEInterpolation()
// ****************************************
// set up the FE interpolation used for enriched approximation 
{
  if (!enrichmentFEInterpolation)
	 enrichmentFEInterpolation = new FEI2dTriLin(1,2);

  return static_cast<FEI2dTriLin*>(enrichmentFEInterpolation);
}

std::vector<DelaunayTriangle *>* Tri6_U::PartitionMySelf()
// *******************************************************
// Partitions split elements and tip elements into sub-triangles
// for numerical integration of the weak form.
// For computational cost, these sub-triangles are already defined in local coord. system
// of the parent 3 noded triangle element.
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
	 myGeo = (*i)->giveMyGeo();                // get the geometry of the enrichment item
	 temp  = myGeo->intersection(myTriangle);  // compute the intersection between element and myGeo 
	 if(temp.size() == 1)                      // for the case crack touched the element edge
	 {
		std::cout << "This is element numbered " << this->giveNumber()<< std::endl;
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
		  //delete p ;
		}
		delete s1 ; delete s2 ; delete s3 ; 
	 }
	 intersectPoints.insert(intersectPoints.end(),temp.begin(),temp.end());
  }
  // make the three nodes Mu::Point
  Mu::Point *p0 = new Mu::Point(this->giveNode(1)->giveCoordinate(1),this->giveNode(1)->giveCoordinate(2)) ;
  Mu::Point *p1 = new Mu::Point(this->giveNode(2)->giveCoordinate(1),this->giveNode(2)->giveCoordinate(2)) ;
  Mu::Point *p2 = new Mu::Point(this->giveNode(3)->giveCoordinate(1),this->giveNode(3)->giveCoordinate(2)) ;

  // transform to local coordinate of the parent element
  Mu::Point *localP0 = this->giveXFEInterpolation()->global2Local(domain,this->giveNodeArray(),p0);
  Mu::Point *localP1 = this->giveXFEInterpolation()->global2Local(domain,this->giveNodeArray(),p1);
  Mu::Point *localP2 = this->giveXFEInterpolation()->global2Local(domain,this->giveNodeArray(),p2);

  DelaunayTree *dt = new DelaunayTree(localP0,localP1,localP2) ;
  
  // insert the intersection points into the DT
  for(size_t j = 0 ; j < intersectPoints.size() ; j++)
  {
	 //intersectPoints[j].print(); // DEBUG ONLY. REMEMBER TO COMMENT THIS LINE !!!
	 dt->insert(standardFEInterpolation->global2Local(domain,this->giveNodeArray(),&intersectPoints[j]));
  }

  // look for the kink point. Is this still necessary if we use DT->insert(Segment*) ???
  
  // do the Delaunay triangulation ...
  vector<DelaunayTriangle *> *tri = dt->getTriangles(); 
  
  delete myTriangle ; 
  delete p0 ; delete p1 ; delete p2 ;
  
  return tri;
}

double Tri6_U :: computeAreaAboveEnrItem(EnrichmentItem *enrItem)
// **************************************************************
// The enrItem splits Tri_U into two portions : compute the area of the
// portion above the enrItem.
{
  GeometryEntity* myGeo = enrItem->giveMyGeo() ;
  Mu::Triangle *tri = this->makeTriangle();
  std::vector<Mu::Point>  intersect = myGeo->intersection(tri); 

  //std::cout << "Dealing with element " << this->giveNumber() << std::endl ; 

  // Attention for the case the crack touches edge of the triangle or goes throught the node
  // these case because of numerical stable, return false when use Segment :: intersection()!!!
  // Therefore, need to use Segment::on(Point* p) to check.
  if (intersect.size() == 1)
  {
	 Mu::Point *p1 = this->giveNode(1)->makePoint();
	 Mu::Point *p2 = this->giveNode(2)->makePoint();
	 Mu::Point *p3 = this->giveNode(3)->makePoint();

	 PiecewiseLinear *geo = static_cast<PiecewiseLinear*>(enrItem->giveMyGeo());
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
	 delete p1 ; delete p2 ; delete p3 ; 
  }

  // just do for the case enrItem completely splits the element => always have 2 intersection points
  if (intersect.size() != 2)
  {
	 std::cout << "enrItem must completely splits the element !!! " << std::endl;
	 assert(false);
  }

  double x1 = intersect[0].x ; // first intersection point
  double y1 = intersect[0].y ; // first intersection point

  double x2 = intersect[1].x ; // second intersection point
  double y2 = intersect[1].y ; // second intersection point

  double x,y ;   // contains coords of each node of element
  double delta ; // signed area of the triangle, ( orientation test)

  Node *nodeI ;

  DelaunayTree *dt;
  dt = NULL ;

  for(size_t i = 0 ; i < numberOfNodes ; i++)
  {         
	 nodeI = this -> giveNode(i+1) ;
	 x = nodeI -> giveCoordinate(1);
	 y = nodeI -> giveCoordinate(2);
	 delta = (x1-x)*(y2-y)-(x2-x)*(y1-y);
	 if ( delta > 32*std::numeric_limits<double>::epsilon())// nodeI is above the intersection segment 
	 {
		if(dt == NULL)
		  dt = new DelaunayTree(&intersect[0],&intersect[1],nodeI->makePoint());
		dt->insert(nodeI->makePoint());
	 }
  }	

  vector<DelaunayTriangle *> * v  = dt->getTriangles(); 

  double A = 0.0  ;
  for(size_t i = 0 ; i < v->size() ; i++)
	 if((*v)[i]->isTriangle == true)
	 {
		A  += (*v)[i]->area();
	 }

  return A ;
}

double Tri6_U :: area() 
//*********************
{
  double x0 = this->giveNode(1)->giveCoordinate(1);
  double y0 = this->giveNode(1)->giveCoordinate(2);

  double x1 = this->giveNode(2)->giveCoordinate(1);
  double y1 = this->giveNode(2)->giveCoordinate(2);

  double x2 = this->giveNode(3)->giveCoordinate(1);
  double y2 = this->giveNode(3)->giveCoordinate(2);

  double A = 0.5 * ((x0-x2)*(y1-y2) - (x1-x2)*(y0-y2)) ; // signed area 

  return A ;
}

Mu::Triangle*  Tri6_U::makeTriangle()
//***********************************
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

GaussPoint**  Tri6_U :: setGaussQuadForJ_Integral()
// ************************************************
// Two Gauss quadrature rules are used
// 1. Classical elements : 13 Gauss points are used
// 2. Split and tip elements : 13 Gps for each sub-triangles

{

  IntegrationRule *GaussQuad = new IntegrationRule ;

  if (!enrichmentItemListOfElem) // no need discontinuous Gauss quad 
  { 
	 GaussQuad = new StandardGaussLegendreQuadrature;
	 GaussQuad->setUpIntegrationPoints(::TRIANGLE,13);
  }
  else
  {
	 GaussQuad = new SplitGaussLegendreQuadrature;
	 GaussQuad->setUpIntegrationPoints(this,7);
  } 


  // tranform the obtained integration points into the GaussPoint array of element

  vector<double>* weight = GaussQuad->giveWeightArray();
  vector<Mu::Point*>* coordArray = GaussQuad->giveIntegrationPointVector();	

  GaussPoint** ret  = new GaussPoint* [coordArray->size()] ;

  for(size_t i = 0 ; i < coordArray->size() ; i++)
  {
	 ret[i] = new GaussPoint(this,i+1,coordArray->at(i),weight->at(i),4) ;
  }

  return ret ;

}