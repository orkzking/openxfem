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

#include "crackinterior.h"

#include "cracktip.h"
#include "crackjunction.h"
#include "vertex.h"
#include "element.h"
#include "node.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <list>
#include <iostream>
#include <algorithm>

CrackInterior :: CrackInterior(int n,Domain* aDomain)
: EnrichmentItem(n,aDomain)
{
  a = 0.0 ;
  size_t check = this->readIfHas("isBranchOf");
  if(check)
  {
	 std::cout << "Hi, I am a branched crack" << endl;
	 isBranchedCrack = true ;
	 myMasterCrack = dynamic_cast<CrackInterior*>(domain->giveEnrichmentItem(check));
  }
  else
  {
	 std::cout << "Hi, I am a master crack" << endl;
	 isBranchedCrack = false ;
	 myMasterCrack = NULL ;
  }
}

void CrackInterior :: getMyTips()
// ******************************
{
  size_t numOfTips = this->readInteger("myTips");

  for(size_t i = 0 ; i < numOfTips ; i++)
  {
	 size_t tip = this->readInteger("myTips",i+2) ;
	 myTips.push_back(static_cast<CrackTip*>(domain->giveEnrichmentItem(tip)));
	 static_cast<CrackTip*>(domain->giveEnrichmentItem(tip))->setMyAssociatedCrackInterior(this);
  }
}

std::vector<CrackTip*> CrackInterior :: giveMyTips()
// **************************************************
{
  if(myTips.size() == 0)
	 this->getMyTips();
  return myTips ;
}

void CrackInterior :: resolveLinearDependency()
// *********************************************
{
  std::vector<Element*>*  elemList = this->giveElementsInteractWithMe(); 

  /*// remove tip elements from listOfInteractedElements
  std::vector<CrackTip*>  tips = this->giveMyTips(); 
  std::vector<Element*>   tipElements ;
  std::vector<Element*>   *temp ;
  for(size_t i = 0 ; i < tips.size() ; i++)
  {
  temp = tips[i]->giveElementsInteractWithMe();
  tipElements.insert(tipElements.end(),temp->begin(),temp->end());
  }

  for(size_t i = 0 ; i < tipElements.size() ; i ++)
  {
  std::vector<Element*>::iterator it = std::find(elemList->begin(),elemList->end(),tipElements[i]);
  if(it != elemList->end())
  elemList->erase(it);
  } */

#ifdef NVP5JULY05
  std::cout<< "Elements split by this CrackInterior : " << std::endl ;
  for(size_t i = 0 ; i < elemList->size() ; i++)
	 std::cout << (*elemList)[i]->giveNumber() << " " ;
  std::cout << std::endl ;
#endif

  Element *e  ;

  // loop on elements that interact with the receiver
  // to avoid doing the same operation on repeated nodes, at first get the
  // vector of nodes of interacted elements, then loop on it and perform ...
  std::list<Node*> checkedNodes;
  for(size_t i = 0 ; i < elemList->size() ; i++)
  {
	 e = (*elemList)[i];
	 for(size_t j = 0 ; j < e->giveNumberOfNodes(); j++)       		
		checkedNodes.push_back(e->giveNode(j+1));		 
  }

  // get rid of repeated nodes ...
  checkedNodes.sort();
  checkedNodes.unique();

  GeometryEntity *geoOfEnrItem = this->giveMyGeo() ;
  if(isBranchedCrack == false)
  {
	 for(std::list<Node*>::iterator it = checkedNodes.begin() ; it != checkedNodes.end() ; it++)
	 {
		if((*it)->isStepEnriched())
		{
		  // if node lies on the crack, always enriched by H(x),i.e., no need to check :)
		  Mu::Point *p = (*it)->makePoint() ;
		  int check = geoOfEnrItem->givePositionComparedTo(p);
		  if(check == 1 || check == -1)
			 (*it) -> resolveLinearDependency(this);
		  delete p ;
		}
	 }
  }
  else
  {
	 for(std::list<Node*>::iterator it = checkedNodes.begin() ; it != checkedNodes.end() ; it++)
	 {
		if((*it)->isStepEnriched() && (*it)->isJunctionEnriched() == false)
		{
		  // if node lies on the crack, always enriched by H(x),i.e., no need to check :)
		  Mu::Point *p = (*it)->makePoint() ;
		  int check = geoOfEnrItem->givePositionComparedTo(p);
		  if(check == 1 || check == -1)
			 (*it) -> resolveLinearDependency(this);
		  delete p ;
		}
	 }
  }
}

void CrackInterior :: resolveConflictsInNodalEnrichment()
// If CrackInterior and MaterialInterface coincide, not enriched
// by both of them, just enriched by H(x) => discontinuous strain 
// Remark: just called for interfacial crack cuts element !!! Otherwise,
// this functions is no need to be performed!!!
// 1-1-2006
{

  std::vector<Element*>*  elemList = this->giveElementsInteractWithMe();
  Element *elem  ;
  std::list<Node*> checkedNodes;

  for(size_t i = 0 ; i < elemList->size() ; i++)
  {
	 elem = (*elemList)[i];
	 for(size_t j = 0 ; j < elem->giveNumberOfNodes(); j++) 
		if(elem->giveNode(j+1)->isTipEnriched() == false)
		  checkedNodes.push_back(elem->giveNode(j+1));		 
  }

  for(std::list<Node*>::iterator it = checkedNodes.begin() ; it != checkedNodes.end() ; it++)
  {
	 (*it)->resolveConflictsInEnrichment2();
  }
}

void CrackInterior :: UpdateMyGeometry()
// **************************************
// When cracks grow, need to insert new tips into list of vertices correctly !!!
// Also insert new tip segment(s) into list of segments
// Code's convention : 
//        + Crack with two tips : [left tip right tip]
//        + Crack with one tip  : [right tip] 
// 2005-08-30
{
  Vertex* tip ;
  Mu::Segment *s,*seg ;
  PiecewiseLinear* geo = static_cast<PiecewiseLinear*>(myGeometry);

  if(myTips.size() == 1)                // CrackInterior with ONE tip
  {
	 tip = myTips[0]->giveTip();
	 s   = myTips[0]->giveTipSegment();
	 seg = new Mu::Segment(*s->first(),*s->second()); // insert a copy of tip segment not itself !!! 2005-09-26

	 //std::list<Mu::Point*> *vertices = geo->giveMyListOfVertices() ;
	 //std::list<Mu::Point*>::iterator it1, it2;
	 //it1 = vertices->begin() ; it2 = vertices->end() ;
	 //if( (*--it2)->x > (*it1)->x ) // right tip
	 //{
	 geo->insertNewSegmentAtBack(seg);
	 geo->insertNewVertexAtBack(tip) ; 
	 //}
	 //else                          // left tip
	 //{
	 //	geo->insertNewSegmentInFront(seg);
	 //	geo->insertNewVertexInFront(tip) ; 
	 //}
  }
  else if (myTips.size() == 2)          // CrackInterior with TWO tips
  {
	 for(size_t i = 0 ; i < myTips.size() ; i++)
	 {
		tip = myTips[i]->giveTip();
		s = myTips[i]->giveTipSegment();
		// DEBUG ONLY
		//std::cout << " Working with tip " << i << std::endl ;
		//tip->print() ;
		//s->print() ;
		if(i==0)      // first tip, i.e., left tip !!!
		{
		  geo->insertNewVertexInFront(tip) ; // insert at front of the list !!!
		  // insert a copy of tip segment not itself !!! 2005-09-26
		  seg = new Mu::Segment(*s->first(),*s->second());
		  geo->insertNewSegmentInFront(seg);
		}
		else if(i==1) // second tip, i.e., right tip !!!
		{
		  geo->insertNewVertexAtBack(tip) ; // insert at front of the list !!!
		  // insert a copy of tip segment not itself !!! 2005-09-26
		  seg = new Mu::Segment(*s->first(),*s->second());
		  geo->insertNewSegmentAtBack(seg) ;
		}
	 }
  }
}

void CrackInterior :: printYourSelf()
{
  std::cout << " CrackInterior numbered " << this->giveNumber() << endl ;
}

void CrackInterior::treatMeshGeoInteractionForMyTips()
// ****************************************************
// 2005-09-07
// Allow tips touch element edge.
// Not sure, it is efficient.
{
  std::vector<Element*> *interactedElems = this->giveElementsInteractWithMe(); 
  this->giveMyTips() ; // tips of the crack, maybe one or two.
  Element *tipElem = NULL ;
  CrackTip *tip ;

  for(size_t i = 0 ; i < myTips.size() ; i++ ) // loop on tips of crack
  {
	 tip = myTips[i] ;
	 for(size_t j = 0 ; j < interactedElems->size() ; j++ )
	 {
		if( tip->interactsWith((*interactedElems)[j]) )
		{
		  std::cout << " The tip is within element :) " << std::endl ;
		  tipElem = (*interactedElems)[j] ;
		  j = interactedElems->size() ;         // not to check more, Tri's help .
		}
	 }

	 // not yet find out tip element, do more tests
	 if(tipElem == NULL)
	 {
		for(size_t k = 0 ; k < this->domain->giveNumberOfElements() ; k++)
		{
		  Element *e = this->domain->giveElement(k+1) ;
		  if(tip->interactsWith(e))
		  {
			 tipElem = e ;
			 k = this->domain->giveNumberOfElements() ;
		  }
		}
	 }

	 // If the tips touch element edges then need check more ...
	 if(tipElem == NULL) // Ah, I do not like this at all :)
	 {
		std::cout << " The tip touched element edge !!! " << std::endl ;
		double x = tip->giveTip()->giveCoordinate(1) ;
		double y = tip->giveTip()->giveCoordinate(2) ;
		Mu::Point *p = new Mu::Point(x,y) ;
		for(size_t k = 0 ; k < interactedElems->size() ; k++ )
		{
		  if((*interactedElems)[k]->IsOnEdge(p))
		  {
			 std::cout << " Ah, I found you !!! " << std::endl ;
			 tipElem = (*interactedElems)[k] ;
			 std::cout << " Number of tip element  " << tipElem->giveNumber() << std::endl ;
			 k = interactedElems->size() ;      // not to check more, Tri's help .
		  }
		}

		tip->setIsOnElementEdge();

		bool found ;

		std::list<Element*> *supportOfTipElem = tipElem->giveNeighboringElements() ;
		std::list<Element*> :: iterator j = supportOfTipElem->begin() ;
		found = false ;
		while( (!found) && (j != supportOfTipElem->end()) )
		{
		  if((*j)->IsOnEdge(p))
		  {
			 tipElem = *j ;
			 std::cout << " Number of tip element  " << tipElem->giveNumber() << std::endl ; // DEBUG 
			 found = true ;
		  }
		  j++ ;
		}
		delete p ;
	 }                     // end of check ( tipElem = NULL )
	 tipElem->isEnrichedWith(tip) ;
	 tip->setListOfInteractedElements(tipElem);
	 tipElem = NULL ;
  }                       // end on loop on tips 
}

double CrackInterior::giveCrackLength()
{
  if(a == 0.0)
	 this->computeCrackLength();
  return a ;
}

void CrackInterior::computeCrackLength()
{
  PiecewiseLinear *crack_geo = dynamic_cast<PiecewiseLinear*>(this->giveMyGeo()); 
  std::list<Mu::Point*> *vertices  = crack_geo->giveMyListOfVertices();

  std::list<Mu::Point*>::iterator it;

  it = vertices->begin();
  double x1 = (*it)->x ; 
  double y1 = (*it)->y ;

  it = vertices->end();
  double x2 = (*--it)->x ; 
  double y2 = (*--it)->y ;

  a = 0.5*sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)); // half crack length
}

void CrackInterior::treatJunctionEnrichment()
{
  // geometry of master crack
  PiecewiseLinear *geo1 = dynamic_cast<PiecewiseLinear*>(myMasterCrack->giveMyGeo());
  // geometry of the receiver
  PiecewiseLinear *geo2 = dynamic_cast<PiecewiseLinear*>(this->giveMyGeo());
  std::list<Mu::Point*>* poi = geo2->giveMyListOfVertices();

  // compute the intersection (junction point)
  Mu::Point *intersect_point;
  for(std::list<Mu::Point*>::iterator i = poi->begin(); i != poi->end(); i++)
  {
	 //(*i)->print();
	 if( geo1->isBelongToMe(*i) )
	 {
		intersect_point = *i ;
	 }
  }
  // detect element containing this interection point
  std::vector<Element*>* elems = this->giveElementsInteractWithMe();
  Element *junction_elem ;
  for(size_t i = 0 ; i < elems->size() ; i++ )
  {
	 if( (*elems)[i]->isWithinMe(intersect_point) )
	 {
		junction_elem = (*elems)[i] ;
		i = elems->size() ;         // not to check more, Tri's help .
	 }
  }
  // define a CrackJunction object
  size_t num_of_enritems = this->domain->giveNumberOfEnrichmentItems();
  CrackJunction *crack_junction = new CrackJunction(num_of_enritems+1,this->domain);

  //this->domain->increaseNumberOfEnrichmentItems();
  // set components for CrackJunction
  crack_junction->components = new std::vector<CrackInterior*>;
  crack_junction->components->push_back(myMasterCrack);
  crack_junction->components->push_back(this);
  // set enrichment for nodes of "junction" element
  for(size_t i = 0 ; i < junction_elem->giveNumberOfNodes() ; i++ )
  {
	 Node *aNode = junction_elem->giveNode(i+1);
	 aNode->deleteEnrichmentItem(this); // nghi ngo thang ne lam ne!!!
	 aNode->addNewEnrichmentItem(crack_junction);
  }
}


void CrackInterior::treatMeshGeoInteraction(Element *e)
// ****************************************************
// Check if element interacts with enrichment item or not. If so, insert this element
// into the list of each enrichment item
// Modified at 2005-09-07 to make it more efficient than before.
// 28-12-2005: MATERIAL INTERFACE IMPLEMENTATION, ASSUMING 2 MATERIALS 
{
  if( this->interactsWith(e) )
  {
	 e->isEnrichedWith(this);
	 this->setListOfInteractedElements(e);
  }
}