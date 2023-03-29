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



#include "enrichmentfunction.h"
#include "enrichmentitem.h"
#include "discontinuousfunction.h"
#include "homogeneouscrackasymp.h"
#include "bimatcrackasymp.h"
#include "voidenrichfunction.h"
#include "modifiedhomocrackasymp.h"
#include "abssigneddistance.h"
#include <stdlib.h>
#include <vector>
#include <algorithm>

EnrichmentFunction :: EnrichmentFunction(Domain* aDomain,int n)
              : FEMComponent(n,aDomain)
{
  myEnrItems = NULL ;
  activeEnrItem = NULL ;
}

EnrichmentFunction :: ~EnrichmentFunction()
{
 ;
}

EnrichmentFunction*  EnrichmentFunction :: typed ()
// ************************************************
// Returns a new enrichment function,which has the same number than the receiver,
// but is typed (Discontinuous, or Asymptotic,..).
{
   EnrichmentFunction* newEnrichmentFunction ;
   char     type[32] ;

   this -> readString("class",type) ;
   newEnrichmentFunction = this -> ofType(type) ;

   return newEnrichmentFunction ;
}
EnrichmentFunction*  EnrichmentFunction :: ofType (char* aClass)
// **************************************************************
// Returns a new enrichment function,which has the same number than the receiver,
// but belongs to a class (Discontinuous, or Asymptotic,..).
{
   EnrichmentFunction* newEnrichmentFunction ;

   if (! strcmp(aClass,"DiscontinuousField"))
      newEnrichmentFunction = new DiscontinuousFunction(domain,number) ;
   else if (! strcmp(aClass,"HomoElastCrackAsymp1"))
      newEnrichmentFunction = new HomoElastCrackAsymp(domain,number,1) ;
   else if (! strcmp(aClass,"HomoElastCrackAsymp2"))
      newEnrichmentFunction = new HomoElastCrackAsymp(domain,number,2) ;
   else if (! strcmp(aClass,"HomoElastCrackAsymp3"))
      newEnrichmentFunction = new HomoElastCrackAsymp(domain,number,3) ;
   else if (! strcmp(aClass,"HomoElastCrackAsymp4"))
      newEnrichmentFunction = new HomoElastCrackAsymp(domain,number,4) ;
	else if (! strcmp(aClass,"ModifiedHomoCrackAsymp"))
      newEnrichmentFunction = new ModifiedHomoCrackAsymp(domain,number) ;
	else if (! strcmp(aClass,"BiMatCrackAsymp1"))
      newEnrichmentFunction = new BiMatCrackAsymp(domain,number,1) ;
	else if (! strcmp(aClass,"BiMatCrackAsymp2"))
      newEnrichmentFunction = new BiMatCrackAsymp(domain,number,2) ;
	else if (! strcmp(aClass,"BiMatCrackAsymp3"))
      newEnrichmentFunction = new BiMatCrackAsymp(domain,number,3) ;
	else if (! strcmp(aClass,"BiMatCrackAsymp4"))
      newEnrichmentFunction = new BiMatCrackAsymp(domain,number,4) ;
	else if (! strcmp(aClass,"BiMatCrackAsymp5"))
      newEnrichmentFunction = new BiMatCrackAsymp(domain,number,5) ;
	else if (! strcmp(aClass,"BiMatCrackAsymp6"))
      newEnrichmentFunction = new BiMatCrackAsymp(domain,number,6) ;
	else if (! strcmp(aClass,"BiMatCrackAsymp7"))
      newEnrichmentFunction = new BiMatCrackAsymp(domain,number,7) ;
	else if (! strcmp(aClass,"BiMatCrackAsymp8"))
      newEnrichmentFunction = new BiMatCrackAsymp(domain,number,8) ;
	else if (! strcmp(aClass,"BiMatCrackAsymp9"))
      newEnrichmentFunction = new BiMatCrackAsymp(domain,number,9) ;
	else if (! strcmp(aClass,"BiMatCrackAsymp10"))
      newEnrichmentFunction = new BiMatCrackAsymp(domain,number,10) ;
	else if (! strcmp(aClass,"BiMatCrackAsymp11"))
      newEnrichmentFunction = new BiMatCrackAsymp(domain,number,11) ;
	else if (! strcmp(aClass,"BiMatCrackAsymp12"))
      newEnrichmentFunction = new BiMatCrackAsymp(domain,number,12) ;
	else if (! strcmp(aClass,"VoidFunction"))
      newEnrichmentFunction = new VoidFunction(domain,number) ;
	else if (! strcmp(aClass,"AbsSignedDistance"))
      newEnrichmentFunction = new AbsSignedDistance(domain,number) ;
   else {
      printf ("%s : unknown enrichment function type \n",aClass) ;
      exit(0) ;}

   return newEnrichmentFunction ;
}

void EnrichmentFunction :: setMyEnrichmentItem(EnrichmentItem *enrItem)
// ********************************************************************
{
  
  if(myEnrItems == NULL) 
	 myEnrItems = new std::vector<EnrichmentItem*>;

  if( std::find(myEnrItems->begin(),myEnrItems->end(),enrItem)
	 == myEnrItems->end())
	 myEnrItems->push_back(enrItem);
}

void EnrichmentFunction :: findActiveEnrichmentItem(EnrichmentItem *enrItem)
{
  std::vector<EnrichmentItem*>::iterator it ;
  it = std::find(myEnrItems->begin(),myEnrItems->end(),enrItem);

  activeEnrItem = *it ;
}

void EnrichmentFunction :: printMyEnrichmentItem()
// **********************************************
{
  for(size_t i = 0 ; i < myEnrItems->size() ; i++)
   (*myEnrItems)[i]->printYourSelf() ;
}





