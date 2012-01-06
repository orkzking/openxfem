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
 


#ifndef _ENRICHMENTDETECTOR_H_
#define _ENRICHMENTDETECTOR_H_

#include "enrichmentitem.h"


//! Enrichment detector, abstract base class 
/*!
 An enrichment item based on the enrichment dectetor to enrich nodes.
 The kind of EnrichmentDetector will be chosen by read from the data file. It's necessary 
 to make sure that the EnrichmentDetector is compatible with the EnrichmentItem.
 */

class EnrichmentDetector
{
public:
  /**
   Virtual method, derived classes will have to implement this method
	Set enriched nodes caused by an EnrichmentItem
	@param enrItem test \c EnrichmentItem
   */
  virtual void setEnrichedNodes(EnrichmentItem *enrItem) = 0;
};


//! PointDetector is EnrichmentDetector for CrackTip
/*!
 An abstract class
 This class has two derived classes : \c ElementAroundPointDetect and BallAroundPointDetect
 */
class PointDetector : public EnrichmentDetector{};


//! ElementAroundPointDetect
/*!
 Set enriched nodes of element around the crack tip
 */
class ElementAroundPointDetect : public PointDetector
{
public:
  void setEnrichedNodes(EnrichmentItem *enrItem);
};


//! BallAroundPointDetect
/*!
 Set enriched nodes belong to the ball of radius r centered at the crack tip
 See Patrick Laborde et al. (2004) for fixed area enrichment 
 */
class BallAroundPointDetect : public PointDetector
{
public:
  void setEnrichedNodes(EnrichmentItem *enrItem);
};


//! PolyLineDetector is EnrichmentDetector for CrackInterior
/*!
 An abstract class
 This class has two derived classes : \c SplitElementDetect and BandOfWidthDetect
 */
class PolyLineDetector : public EnrichmentDetector{};

//! SplitElementDetect - CrackInterior
/*!
 Set enriched nodes of elements cut by the CrackInterior
 */
class SplitElementDetect : public PolyLineDetector
{
 public:
  void setEnrichedNodes(EnrichmentItem *enrItem);
};

//! \todo for Ravindra :)
class BandOfWidthDetect : public PolyLineDetector
{
 public:
  void setEnrichedNodes(EnrichmentItem *enrItem);
};



#endif // _ENRICHMENTDETECTOR_H_