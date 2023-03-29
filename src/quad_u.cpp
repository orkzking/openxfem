/************************************************************************************

Copyright (C) 2005
Stephane BORDAS, Cyrille DUNANT, Vinh Phu NGUYEN, Quang Tri TRUONG, Ravindra DUDDU

This file is part of the XFEM C++ Library (XFEMLIB) written
and maintained by above authors.

This program is free software; you can redistribute it and/or modify it.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
license text for more details.

Any feedback is welcome. Emails : nvinhphu@gmail.com, ...

*************************************************************************************/
#include "quad_u.h"

#include "gausspnt.h"
#include "diagmtrx.h"
#include "fei2dquadlin.h"
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


Quad4_U :: Quad4_U (int n, Domain* aDomain)
: PlaneProblem(n,aDomain)
{
  numberOfNodes  = 4 ;
}


void  Quad4_U :: computeGaussPoints ()
// **********************************
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

GaussPoint**  Quad4_U :: setGaussQuadForJ_Integral()
// ************************************************
// Sets up the array containing the Gauss points used for the J integral
// computation.
{
  IntegrationRule *GaussQuad = new IntegrationRule ;

  size_t numOfGPs1 = 9 ;
  size_t	numOfGPs2 = 7 ; // number of GPs for each sub-triangle

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

FEI2dQuadLin* Quad4_U::giveFEInterpolation()
// *****************************************
// Returns the shape functions used to interpolate the "true" displacement
{

  if (!standardFEInterpolation)
	 standardFEInterpolation = new FEI2dQuadLin(1,2);

  return static_cast<FEI2dQuadLin*>(standardFEInterpolation);

}


FEI2dQuadLin* Quad4_U::giveXFEInterpolation()
// ******************************************
// Returns the shape functions multiplied with enrichment functions to form
// the PUM shape functions
{
  if (!enrichmentFEInterpolation)
	 enrichmentFEInterpolation = new FEI2dQuadLin(1,2);

  return static_cast<FEI2dQuadLin*>(enrichmentFEInterpolation);
}

std::vector<DelaunayTriangle *>* Quad4_U:: PartitionMySelf()
// ********************************************************
// 2005-08-23.
{
  std::vector<Mu::Point>  temp,pts,intersectPoints ;
  GeometryEntity* myGeo;
  // TIP-ELEMENTS INTERACTED WITH BOTH THE TIP AND  THE CRACK INTERIOR,
  // THEREFORE, NEED TO REMOVE CRACK INTERIOR FROM THE LIST OF ENRICHMENT ITEMS.

  list<EnrichmentItem*> ::iterator iter1,iter2;
  iter1=find_if(enrichmentItemListOfElem->begin(),enrichmentItemListOfElem->end(),IsType<CrackTip,EnrichmentItem>());
  iter2=find_if(enrichmentItemListOfElem->begin(),enrichmentItemListOfElem->end(),IsType<CrackInterior,EnrichmentItem>());

  if (iter1 != enrichmentItemListOfElem->end() && iter2 != enrichmentItemListOfElem->end())
	 enrichmentItemListOfElem->remove(*iter2);

  // LOOP ON ENRICHMENT ITEMS AND GET THE INTERSECTION POINTS
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
		if(temp.size() == 1)     // for the case crack touched the element edge
		{
		  Mu::Point *p1 = this->giveNode(1)->makePoint();
		  Mu::Point *p2 = this->giveNode(2)->makePoint();
		  Mu::Point *p3 = this->giveNode(3)->makePoint();
		  Mu::Point *p4 = this->giveNode(4)->makePoint();
		  Mu::Segment *s1 = new Mu::Segment(*p1,*p2) ;
		  Mu::Segment *s2 = new Mu::Segment(*p2,*p3) ;
		  Mu::Segment *s3 = new Mu::Segment(*p3,*p4) ;
		  Mu::Segment *s4 = new Mu::Segment(*p4,*p1) ;

		  PiecewiseLinear *geo = dynamic_cast<PiecewiseLinear*>(myGeo);
		  std::list<Mu::Point*> *vertices = geo->giveMyListOfVertices();

		  for (std::list<Mu::Point*> ::iterator k = vertices->begin() ; k != vertices->end(); k++)
		  {
			 if( s1->on(*k) || s2->on(*k) || s3->on(*k)|| s4->on(*k) )
			 {
				(*k)->print();
				temp.push_back(**k) ;
			 }
		  }
		  delete s1 ; delete s2 ; delete s3 ; delete s4 ;
		  delete p1 ; delete p2 ; delete p3 ; delete p4 ; // Purify, 15-10-2005
		}

		intersectPoints.insert(intersectPoints.end(),temp.begin(),temp.end());

		// ************************************************************************************
		// ********                    LOOK FOR KINK POINT                               ******
		// ************************************************************************************
		// just segments have kink points, i.e., geometry of CrackInterior only
		// So, use typeid() to test ...
		// Attention: just greater than 2 vertices => kink point :) 2005-09-13

		if( typeid(**i) == typeid(CrackInterior) )
		{
		  std::list<Mu::Point*> *vertices = dynamic_cast<PiecewiseLinear*>(myGeo)->giveMyListOfVertices();
		  if(vertices->size() != 2)
		  {
			 for(std::list<Mu::Point*> ::iterator k = vertices->begin() ; k != vertices->end(); k++)
			 {
				if( this->isWithinMe(*k) )
				  intersectPoints.push_back(**k) ;
			 }
		  }
		}          // end of looking for kink points  */
	 }            // end of check CrackTip
  }              // end of loop on enrichment items

  //for(size_t i = 0 ; i < intersectPoints.size() ; i++)
  //{
	// intersectPoints[i].print();
	// std::cout<<std::endl;
  //}

  // make the four nodes Mu::Point
  Mu::Point *p0 = new Mu::Point(this->giveNode(1)->giveCoordinate(1),this->giveNode(1)->giveCoordinate(2)) ;
  Mu::Point *p1 = new Mu::Point(this->giveNode(2)->giveCoordinate(1),this->giveNode(2)->giveCoordinate(2)) ;
  Mu::Point *p2 = new Mu::Point(this->giveNode(3)->giveCoordinate(1),this->giveNode(3)->giveCoordinate(2)) ;
  Mu::Point *p3 = new Mu::Point(this->giveNode(4)->giveCoordinate(1),this->giveNode(4)->giveCoordinate(2)) ;

  // transform to local coordinate of the parent element
  Mu::Point *localP0 = this->giveFEInterpolation()->global2Local(domain,this->giveNodeArray(),p0);
  Mu::Point *localP1 = this->giveFEInterpolation()->global2Local(domain,this->giveNodeArray(),p1);
  Mu::Point *localP2 = this->giveFEInterpolation()->global2Local(domain,this->giveNodeArray(),p2);
  Mu::Point *localP3 = this->giveFEInterpolation()->global2Local(domain,this->giveNodeArray(),p3);

  DelaunayTree *dt = new DelaunayTree(localP0,localP1,localP2) ;
  dt->insert(localP3);
  // insert the intersection points into the DT
  for(size_t j = 0 ; j < intersectPoints.size() ; j++)
  {
	 dt->insert(standardFEInterpolation->global2Local(domain,this->giveNodeArray(),&intersectPoints[j]));
  }

  // do the Delaunay triangulation ...
  vector<DelaunayTriangle *> *tri = dt->getTriangles();

  delete p0 ; delete p1 ; delete p2 ; delete p3 ;

  return tri;
}


double Quad4_U :: area()
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

double Quad4_U :: computeAreaAboveEnrItem(EnrichmentItem *enrItem)
// **************************************************************
// The enrItem splits Quad_U into two portions : compute the area of the
// portion above the enrItem.
// First, define the intersection points, insert into the DelaunayTree
// perform the tesselation => triangles => area
// Latest version : 2005-09-07, when solve with structured Q4 mesh.
{
  GeometryEntity* myGeo = enrItem->giveMyGeo() ;
  std::vector<Mu::Point> pts,intersect,temp ;

  // compute the intersection between element edges and enrItem
  // -------------------------------------------------------------------------
  for(size_t i = 0 ; i < this->giveNumberOfNodes() ; i++)
  {
	 Mu::Point *p = this->giveNode(i+1)->makePoint() ;
	 pts.push_back(*p);
	 delete p ;
  }

  for(size_t i = 0 ; i < pts.size() ;  i++)
  {
	 Mu::Segment s(pts[i],pts[(i+1)%pts.size()]) ;

	 if(myGeo->intersects(&s))
	 {
		temp = myGeo->intersection(&s) ;
		intersect.insert(intersect.end(),temp.begin(),temp.end()) ;
	 }
  }
  // -------------------------------------------------------------------------

  // If crack touched element edge , the above code just gives one intersection !!!
  // -------------------------------------------------------------------------
  PiecewiseLinear *geo = static_cast<PiecewiseLinear*>(enrItem->giveMyGeo());

  if (intersect.size() == 1)
  {
	 std::cout << " Crack tips touched element edge !!! " << std::endl ;
	 std::list<Mu::Point*> *vertices = geo->giveMyListOfVertices() ;
	 std::list<Mu::Point*>::iterator j = vertices->end();
	 Mu::Point *lastVertex = *(--j) ;
	 lastVertex->print() ;  std::cout << std::endl ;   // debug only
	 for(size_t k = 0 ; k < pts.size() ;  k++)
	 {
		Mu::Segment s(pts[k],pts[(k+1)%pts.size()]) ;
		s.print() ;
		if(s.on(lastVertex))
		{
		  std::cout << " found edge containing tip  " << endl ;
		  intersect.push_back(*lastVertex);
		  k = pts.size() ;
		}
		//for( std::list<Mu::Point*>::iterator j = vertices->begin() ; j != vertices->end() ; j++)
		//  if(s.on(*j))
		//  {
		//	 intersect.push_back(**j);
		//	 j = vertices->end() ;     // stop, no test anymore.
		//  }
	 }
  }               // end of if (intersect.size() == 1)
  // -------------------------------------------------------------------------

  // just do for the case enrItem completely splits the element => always have
  // 2 intersection points
  if (intersect.size() != 2)
  {
	 std::cout << "enrItem must completely splits the element !!! " << std::endl;
	 assert(false);
  }

  Mu::Point *p1 = this->giveNode(1)->makePoint();
  Mu::Point *p2 = this->giveNode(2)->makePoint();
  Mu::Point *p3 = this->giveNode(3)->makePoint();
  Mu::Point *p4 = this->giveNode(4)->makePoint();
  DelaunayTree *dt = new DelaunayTree(p1,p2,p3);
  // insert points into DelaunayTree
  dt->insert(p4) ;
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
		Mu::Point *center = (*v)[i]->getCenter();
		if(x2 > x1)
		  delta = (x1-center->x)*(y2-center->y)-(x2-center->x)*(y1-center->y);
		else
		  delta = (x2-center->x)*(y1-center->y)-(x1-center->x)*(y2-center->y);

		if(delta > 0.000001)
		  A  += (*v)[i]->area();
	 }
  }
  delete p1 ; delete p2 ; delete p3 ; delete p4 ;
  delete dt ; delete v ;  // Purify,14-10-2005

  return A ;

}

