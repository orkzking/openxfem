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

#ifndef _CRACKINTERIOR_H_
#define _CRACKINTERIOR_H_

#include "enrichmentitem.h"
#include "element.h"
#include <vector>

class Vertex; class GeometryEntity;class CrackTip; 


//! The Crack interior 
/*!
 This class implements a crack interior.
 */
class CrackInterior : public EnrichmentItem
{
public:

  CrackInterior(int,Domain*);   //!< constructor
  ~CrackInterior(){;}           //!< destructor
  /**
   Update the geometry of the crack   
   */
  void     UpdateMyGeometry(); //  start on 2005/08/31
  /*!
   Time step termination. 
   */
  void     printOutputAt (TimeStep*, FILE*, FILE*){}
  /*!
   Read from the input file the tips of this crack
   */
  void                   getMyTips();
  std::vector<CrackTip*> giveMyTips();
  /*!
   Resolve the linear dependency for nodes enriched by H(x) function.
   */
  void     resolveLinearDependency();

  void     resolveConflictsInNodalEnrichment();

  /*!
   Update enrichment after cracks grow. Do no thing here.
	See class CrackTip. 
   */
  void     updateEnrichment(){;}

  void     printYourSelf();
  /*!
   Find out elements interacting with tips of the receiver by 
	looping on set of elements interacting with the receiver which has been
	already found.
	2005-09-07.
   */
  void     treatMeshGeoInteractionForMyTips() ;
  double   giveCrackLength();
  void     computeCrackLength();
  bool     IsaBranchedCrack(){return isBranchedCrack;}
  void     treatJunctionEnrichment(); // branched crack, 2005-11-12
  void     treatMeshGeoInteraction(Element* e);
private:
  std::vector<CrackTip*>  myTips;//!< vector of associated crack tips
  double                  a; //!< half crack length
  bool                    isBranchedCrack; // branched crack, 2005-11-12
  CrackInterior*          myMasterCrack ;  // branched crack, 2005-11-12     
} ;


#endif //_CRACKINTERIOR_H_
