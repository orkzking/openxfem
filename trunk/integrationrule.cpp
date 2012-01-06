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

#include "integrationrule.h"
#include "gausspnt.h"
#include "assert.h"
#include <vector>


IntegrationRule::IntegrationRule ()
{
  numberOfIntegrationPoints = 0;
  weightArray = NULL ;    
  integrationPointVector = NULL ;  
}


IntegrationRule::~IntegrationRule ()
{
  delete weightArray ;
  delete integrationPointVector ;
}

std::vector<double>*  IntegrationRule ::  giveWeightArray()
// return a copy of weight so that we can delete 
{
  std::vector<double>* weight = new std::vector<double>(0) ;
  weight->insert(weight->end(),weightArray->begin(),weightArray->end());
  return weight ;
}

std::vector<Mu::Point*>* IntegrationRule :: giveIntegrationPointVector()
// return a copy version of integrationPointVector => can delete temporary variables
{
  std::vector<Mu::Point*>* coordArray = new std::vector<Mu::Point*>(0);
  coordArray->insert(coordArray->end(),integrationPointVector->begin(),integrationPointVector->end());
  return coordArray ;   
}

void IntegrationRule :: setUpIntegrationPoints (IntDomain mode, size_t nPoints)
{

  if (mode == LINE)
	 this->SetUpPointsOnLine (nPoints) ;
  else if (mode == ::TRIANGLE)
	 this->SetUpPointsOnTriangle (nPoints) ;
  else if (mode == SQUARE)
	 this->SetUpPointsOnSquare  (nPoints) ;
  else if (mode == TETRAHEDRA)
	 this->SetUpPointsOnTetrahedra(nPoints) ;
  else
	 assert(false);
}