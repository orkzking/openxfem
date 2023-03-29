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

#ifndef _CRACKTIP_H_
#define _CRACKTIP_H_

#include "enrichmentitem.h"
#include "auxiliaryfield.h"
#include "intarray.h"
#include "delaunay.h"
#include "material.h"
#include "timestep.h"
#include <vector>
#include <list>

class CrackGrowthIncrementLaw; class CrackGrowthDirectionLaw; 
class Vertex; class GeometryEntity; class CrackInterior ;

typedef enum{ // There are two crack types : homogeneous and bimaterial crack
  HomoElast,
  BiMatElast,
} CrackType;


//! This class implementes a crack tip 
/*! 
 This class together with the associated class CrackInterior, model 2D cracks.
 */
class CrackTip : public EnrichmentItem
{
public:

  CrackTip(int,Domain*);   //!< constructor
  ~CrackTip();             //!< destructor

  /*!
   Returns the tip of the receiver
   */
  Vertex*     giveTip();
  /*!
   Returns the field of the crack
   @return : PLANESTRAIN or PLANESTRESS 
   */
   FieldType  giveMyFieldForSIF() ;
  /*!
   Define the  integration domain for interaction integral.
   Be a set of elements that belong to a ball of radius r = h*sqrt(tip_element)
   centered at the crack tip. Often, h = 2.0.
   */
  std::list<Element*> buildIntegrationDomain1();
  /*!
   Set of elements intersecting with the boundary of the circle radius r = h*sqrt(tip_element)
   centered at the crack tip. Often, h = 2.0.
   */
  std::list<Element*> buildIntegrationDomain2();

  /** @name Computation of SIFs
   *  Computation of the stress intensity factors 
   */
  //@{
  void   computeSIFs(TimeStep* stepN);
  void   computeK_eq(TimeStep* stepN); 
  //@}

  /** @name Return the SIFs
   *  Give the the stress intensity factors 
   */
  //@{
  double giveK_i(){return K_i ;}     // used in MaxHoopStress to compute growth direction
  double giveK_ii(){return K_ii ;}
  //@}

  /*!
   Compute the interaction integrations I(1,2).
	@return I[0] is I(1,2) for mode I and I[1] is I(1,2) for mode II.
   Then used to compute the SIFs.
  */
  std::valarray<double> computeInteractionIntegral(TimeStep* stepN);
  /*!
   Return the state of the receiver.
   */
  bool giveState()const {return isActive;}
  /**
   Kill the crack tip. Whenever the tip reaches a free boundary, we kill it.
   If the tip reaches another crack, then we also kill this tip and introduce
   a junction. See Daux et all 2000 and Budyn 2003.
  */
  void kill(){isActive = false;}
  /**
   Read from the input file the type of the crack
   At this moment, just consider two cracks : homogeneous elastic cracks
   and interfacial cracks.
   */
  void  crackTypeInitialization();
  /**
   Check the current position of the crack tip, then deduce the type of 
   the crack, for example, if the tip is in the MaterialInterface => Bimaterial
   crack
   */
  void  crackTypeUpdate();
  /**
   Reads from the input file the materials associated with the reveiver
	@return one material if crack is HomoElast and two materials for 
	cracks of type BiMatElast
   */
  std::vector<Material*>*   giveMatArray();
  /*!
   Return the segment containing the TIP
   */
  Mu::Segment*   giveTipSegment();
  /*!
   Compute the coord in local crack tip coord system
	@param p : point needs to be computed the local coord
	@return FloatArray* ret , (*ret)[0] = xloc and (*ret)[1] = yloc
   */
  /*!
   @Return true if the updated tip still inside the tip-element.
   */
  bool           InOrOutTipElement(){return stillInTipElement ;}  
  /*!
   Compute the XY coordinates of Point p in the local crack tip coordinate system
   */
  FloatArray*    computeLocalCoordOf(Mu::Point* p);
  /*!
   Compute the polar coordinates of Point p in the local crack tip coordinate system
   */
  std::valarray<double>*    computePolarCoordOf(Mu::Point* p);
  /*!
	Time step termination
	*/
  void           printOutputAt (TimeStep*, FILE*, FILE*);
  /*!
   @return the radius of domain integration for SIFs computation
   */
  double         giveRadiusOfDomainIntegration();
  /*!
   @return the enrichment radius used in fixed enrichment area detector.
	Used with EnrichmentDetector : BallAroundMe only.
   */
  double         giveEnrichRadius();
  /*!
   Do no thing. Just for polymorphism. Only CrackInterior implement this method.
   */
  void           resolveLinearDependency(){}
  /*!
   Do not enrich nodes of tip elements with H(x).
	2005-09-11
   */
  void           resolveConflictsInNodalEnrichment();

  void           printYourSelf();

  /** @name Update geometry when crack grows
   *  Start from 2005-08-29.
   */
  //@{
  /**
   Update the geometry of the crack   
   */
   void          UpdateMyGeometry(); 
  /**
   Define the circle centered at the old tip with radius equal to distance between
	the old and new tips.
	This circle will be used to efficiently find out new enriched nodes. 
   */
	Mu::Circle*   DefineDomainForUpdatedEnrichment();
  /**
   Find elements which are risk to interact with the new crack segment.
	These elements are defined based on the circle built from the previous method.
   */
	std::list<Element*>  DefineUpdatedElements();
  /**
   Update the enrichment
	Just do here, no need to do with CrackInterior !!!
   */
	void   updateEnrichment(); 
  //@}

	void    setMyAssociatedCrackInterior(CrackInterior* crInt){ myAssociatedCrackInterior = crInt ;}
	CrackInterior* giveMyCrackInterior(){return myAssociatedCrackInterior;} 
	bool    amIOnElementEdge(){return isOnElementEdge ;}
	void    setIsOnElementEdge(){isOnElementEdge = true ;}
	/**
	 Inclination angle of the tip segment w.r.t the horizontal axe 
	 */
	double  giveInclinationAngle();
	void    computeInclinationAngle();
	void    treatMeshGeoInteraction(Element* e){}
private:
  CrackType  tipID;                  //!< Homo or interfacial cracks  
  FieldType  field;                  //!< PLANESTRAIN or PLANESTRESS 
  std::vector<Material*>*  matArray; //!< Array of materials
  double K_i,K_ii;                   //!< SIFs of mode I and II
  double K_eq  ;                     //!< equivalent SIF 
  Mu::Segment*       tipSegment ;    //!< crack segment containing the tip
  bool    stillInTipElement ;        //!< marker, true if the updated tip still locates in tip element
  bool                      isActive;    //!< switch isActive state on or off
  Mu::Circle*   makeBall();              //!< make ball centered at the tip, radius = 2*tip_element.area()
  double computeSIFMultiplier();
  double computeSIFMultiplierForOneMat();
  double computeSIFMultiplierForBiMat();
  double enrichRadius; //!< radius of enrichment area,often equals to 1/10 dimension of domain
  CrackInterior*   myAssociatedCrackInterior ; //!< Crack interior of this tip.
  bool             isOnElementEdge ; //!< marker, true if tip is on element edge.
  double           alpha ; //!< inclination angle of tip segment
} ;

bool   nodeIsInCircle(Node* ,Mu::Circle*);

#endif  // _CRACKTIP_H_
