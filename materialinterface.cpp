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
   
	Any feedback is welcome. Emails : nvphu80@yahoo.com, ...

*************************************************************************************/

#include "materialinterface.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "element.h"
#include "node.h"
#include "geometryentity.h"

MaterialInterface :: MaterialInterface(int n,Domain* aDomain)
              : EnrichmentItem(n,aDomain)
{
}

std::vector<size_t>  MaterialInterface::giveMaterialIDs()
// ******************************************************
// Read from the input file the material info
{
  if (materialIDs.size()==0)
  {
	 size_t nMats  = this->readInteger("Mat") ; // number of materials
	 for (size_t i = 0 ; i < nMats ; i++)
	 {
		size_t matNumber = this->readInteger("Mat",i+2) ;
		materialIDs.push_back(matNumber);
	 }
  }

  return materialIDs ;
}

void MaterialInterface::printYourSelf()
{
  std::cout << "I am a Material Interface " << std::endl;
}

void MaterialInterface::treatMeshGeoInteraction(Element *e)
// ********************************************************
// Check if element interacts with enrichment item or not. If so, insert this element
// into the list of each enrichment item
// Modified at 2005-09-07 to make it more efficient than before.
// 28-12-2005: MATERIAL INTERFACE IMPLEMENTATION, ASSUMING 2 MATERIALS 
{
  if( this->interactsWith(e) )
  {
	 e->isEnrichedWith(this);
	 this->setListOfInteractedElements(e);
	 e->computeNodalLevelSets(this);
	 e->setMultiMaterial();
	 e->setMaterialIDs(this->giveMaterialIDs());
  }
  else
  {
	 Mu::Point *p = e->giveNode(1)->makePoint();
	 double ls = this->giveMyGeo()->computeSignedDistanceOfPoint(p);
	 if(ls > 0)
		e->setMaterial(this->giveMaterialIDs()[0]);
	 else
		e->setMaterial(this->giveMaterialIDs()[1]);
  }
  
}
