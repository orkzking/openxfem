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

#include "junctionfunction.h"
#include "crackinterior.h"
#include "crackjunction.h"
#include "element.h"
#include "node.h"

double JunctionFunction :: EvaluateYourSelfAt(GaussPoint* gp)
// **********************************************************
// Compute the value of junction function J(x) at coord.
// J(x) is computed based on the value of H1(x) and H2(x).
{
  std::vector<CrackInterior*> *compo = 
  dynamic_cast<CrackJunction*>(this->activeEnrItem)->giveComponents();

  vector<EnrichmentFunction*>* enrFns1 = (*compo)[0]->giveEnrFuncVector();
  vector<EnrichmentFunction*>* enrFns2 = (*compo)[1]->giveEnrFuncVector();

  double H1 = (*enrFns1)[0]->EvaluateYourSelfAt(gp);
  double H2 = (*enrFns2)[0]->EvaluateYourSelfAt(gp);
 
  return ( H1 < 0 ) ? H2 : H1 ;
}

double JunctionFunction :: EvaluateYourSelfAt(Element* e,Node* n)
{
  Point* p = n->makePoint() ; 
  std::vector<CrackInterior*> *compo = 
  dynamic_cast<CrackJunction*>(this->activeEnrItem)->giveComponents();

  vector<EnrichmentFunction*>* enrFns1 = (*compo)[0]->giveEnrFuncVector();
  vector<EnrichmentFunction*>* enrFns2 = (*compo)[1]->giveEnrFuncVector();

  double H1 = (*enrFns1)[0]->EvaluateYourSelfAt(e,n);
  double H2 = (*enrFns2)[0]->EvaluateYourSelfAt(e,n);
  
  delete p ; 
  return ( H1 < 0 ) ? H2 : H1 ;
}

FloatArray* JunctionFunction :: EvaluateYourGradAt(GaussPoint* gp)
// ***************************************************************
// Compute the derivatives of discontinuous function H(x)
// This derivative is zero 
{
	return new FloatArray(2);
}
