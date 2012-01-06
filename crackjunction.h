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

Any feedback is welcome. Emails : nvinhphu@gmail.com,stephane.bordas@epfl.ch

*************************************************************************************/

#ifndef _CRACKJUNCTION_H_
#define _CRACKJUNCTION_H_

#include "enrichmentitem.h"
#include "element.h"
#include <vector>

class GeometryEntity; class CrackInterior; 

//! The Crack junction 
/*!
This class implements a crack interior.
*/
class CrackJunction : public EnrichmentItem
{
friend class CrackInterior;
public:

  CrackJunction(int,Domain*);   //!< constructor
  ~CrackJunction(){;}           //!< destructor
  /**
  Update the geometry of the crack   
  */
  void     UpdateMyGeometry(){}
  /*!
  Time step termination. 
  */
  void     printOutputAt (TimeStep*, FILE*, FILE*){}
  /*!
  Resolve the linear dependency for nodes enriched by H(x) function.
  */
  void     resolveLinearDependency(){}

  void     resolveConflictsInNodalEnrichment(){}

  /*!
  Update enrichment after cracks grow. Do no thing here.
  See class CrackTip. 
  */
  void     updateEnrichment(){;}

  void     printYourSelf();
  std::vector<CrackInterior*> * giveComponents(){return this->components;} 
  size_t   giveNumberOfInstances(){return numOfInstances;}
protected:
  std::vector<CrackInterior*> *components;
  size_t                numOfInstances; // number of CrackJunction object
} ;


#endif //_CRACKJUNCTION_H_
