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

#include "crackjunction.h"

#include "junctionfunction.h"
#include "element.h"
#include "node.h"
#include <stdio.h>
#include <vector>
#include <list>
#include <iostream>


CrackJunction :: CrackJunction(int n,Domain* aDomain)
              : EnrichmentItem(n,aDomain)
{
  size_t num_of_enrFns = this->domain->giveNumberOfEnrichmentFunctions();
  EnrichmentFunction *enrFn = new JunctionFunction(this->domain,num_of_enrFns+1);
  //this->domain->increaseNumberOfEnrichmentFunctions();
  myEnrichFns = new vector<EnrichmentFunction*> ;
  myEnrichFns->push_back(enrFn) ;
  enrFn->setMyEnrichmentItem(this) ;
  components = NULL ;
  numOfInstances = 0 ;
  ++numOfInstances;
}

void CrackJunction :: printYourSelf()
{
  std::cout<< " a CrackJunction " << " " ;
}

