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

Any feedback is welcome. Emails : nvphu80@yahoo.com,stephane.bordas@epfl.ch

*************************************************************************************/


#include "enrichmentdetector.h"

#include "enrichmentitem.h"
#include "domain.h"
#include "element.h"
#include "delaunay.h"
#include "geometry_2D.h"
#include "cracktip.h"
#include "crackinterior.h"
#include "vertex.h"
#include "node.h"
#include "assert.h"
#include <typeinfo>
#include <iostream>
#include <vector>
#include <list>


void ElementAroundPointDetect :: setEnrichedNodes(EnrichmentItem *enrItem)
// ***********************************************************************
// enrItem must be CrackTip !!!
// enrich all nodes of tip-element.
// NVP 2005 06 03 
{
  // checking for compatibility ...
  if (!(typeid(static_cast<CrackTip*>(enrItem)) == typeid(CrackTip*)))
  {
	 std :: cout << " Enrichment item and Enrichment detector is not compatible together "<< std :: endl ;
	 assert(false);
  }

  std::vector<Element*>* elemList = enrItem->giveElementsInteractWithMe(); 

  // a Tip just interacts with ONE element !!!
  assert(elemList->size() == 1);

  (*elemList)[0]->setEnrichmentForMyNodes(enrItem) ;
}

void BallAroundPointDetect :: setEnrichedNodes(EnrichmentItem *enrItem)
// ********************************************************************
// enrItem must be CrackTip !!!
// nodes belong to the circle are enriched (Fixed area enrichment).
// NVP 2005-06-03 
{
  // checking for compatibility ...
  if (!(typeid(static_cast<CrackTip*>(enrItem)) == typeid(CrackTip*)))
  {
	 std :: cout << " Enrichment item and Enrichment detector is not compatible together "<< std :: endl ;
	 assert(false);
  }

  double enrichRadius = static_cast<CrackTip*>(enrItem)->giveEnrichRadius();
  double xTip = static_cast<CrackTip*>(enrItem)->giveTip()->giveCoordinate(1);
  double yTip = static_cast<CrackTip*>(enrItem)->giveTip()->giveCoordinate(2);
  Mu::Point  *center = new Mu::Point(xTip,yTip);
  Mu::Circle *ball   = new Mu::Circle(enrichRadius,center);

  // find nodes belong to the ball ...
  vector<Element*> *interactedElements = enrItem->giveElementsInteractWithMe();
  Element *tipElem = (*interactedElements)[0] ; 

  std::list<Element*> ret = tipElem->conflicts(ball);

  Node *node; 
  Mu::Point *poi ;
  for(std::list<Element*> :: iterator i = ret.begin() ; i != ret.end() ; i++)
  {
	 for(size_t j = 0 ; j < (*i)->giveNumberOfNodes() ; j++)
	 {
		node = (*i)->giveNode(j+1);
		poi = node->makePoint();
		if(ball->in(poi))                    // node belongs to the ball, then enriched...
		  node->isEnrichedWith(enrItem);
		delete poi ;                         //Purify,14-10-2005
	 }
  }
  delete center ; delete ball ;            //Purify,14-10-2005
}

void SplitElementDetect :: setEnrichedNodes(EnrichmentItem *enrItem)
// *****************************************************************
// enrItem must be CrackInterior !!!
// NVP 2005-06-03 
{
  // checking for compatibility ...
  //if (!(typeid(static_cast<CrackInterior*>(enrItem)) == typeid(CrackInterior*)) )
  //{
	// std :: cout << " Enrichment item and Enrichment detector is not compatible together "<< std :: endl ;
	 //assert(false);
  //}

  std::vector<Element*>* elemList = enrItem->giveElementsInteractWithMe(); 

  for(size_t i = 0 ; i < elemList->size() ; i++)
  {
    (*elemList)[i]->setEnrichmentForMyNodes(enrItem) ;
  }
}

void BandOfWidthDetect :: setEnrichedNodes(EnrichmentItem *enrItem)
// ****************************************************************
// 
{
}
