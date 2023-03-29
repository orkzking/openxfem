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

#ifndef _ENRICHMENTITEM_H_
#define _ENRICHMENTITEM_H_

#include "femcmpnn.h"
#include "domain.h"
#include <vector>

class GeometryEntity; class EnrichmentFunction;
class Element; class Vertex;class TimeStep ; class EnrichmentDetector;

using namespace std;

//! Enrichment items : cracks, holes, material interfaces...
/*!
This class is an abstract class, the superclass of all classes that imple-
ment objects that cause enrichment : cracks, holes, material interfaces,
and biofilms...\n
*/
class EnrichmentItem:public FEMComponent
{

public:
  EnrichmentItem(int,Domain*); //!< Constructor
  virtual ~EnrichmentItem();   //!< Destructor

  /*! Get the initial geometry of enrichment item from the input file
   */
  void                 getGeometry(); 
  /*!
   Reads from the input file the enrichment functions used to
   model this enrichment item
   */
  vector< EnrichmentFunction* >* giveEnrFuncVector(); 
  /*!
   Check if enrichment item interacts with element or not 
   Use its geometry to check this interaction
   */
  bool interactsWith(Element*);
  /*!
   Enriched nodes of element with this enrichment item
   Based on the \c EnrichmentDetector
   */
  virtual void treatEnrichment();
  /** @name Concerning the enrichment detector of enrichment item
   *  Functions to define and give the enrichment detector.
   */
  //@{
  EnrichmentDetector*   defineMyEnrDetector();
  EnrichmentDetector*   giveMyEnrDetector();
  //@}

  /*!
   Returns the geometry of the receiver.
   */
  GeometryEntity*        giveMyGeo(); 
  /*!
   Set the list of elements interact with the receiver
   */
  void                   setListOfInteractedElements(Element *e); 
  /*!
   Returns a list of elements that interact with the receiver.
   This list is created and inserted new elements when doing the mesh geo interaction.
   */
  vector<Element*>*      giveElementsInteractWithMe();
  /** @name Definition of an enrichment item
   *  Functions to define an enrichment item
   */
  //@{
  EnrichmentItem*      typed () ;        
  EnrichmentItem*      ofType (char*) ;
  char*                giveClassName (char* s)
          { return strcpy(s,"EnrichmentItem") ;}
  int                  giveNumber ()
          { return FEMComponent::giveNumber() ;} 
  //@}
  /*!
  Time step termination
  */
  virtual void     printOutputAt (TimeStep*, FILE*, FILE*){}
  /*!
   If a node is enriched by the Heaviside function H(x) and it does not satisfy
   the area inclusion criteria ( see Dolbow et al. 1999) then this node is not
   enriched by H(x) any more => matrix is not singular.
   Virtual method, only derived class CrackInterior overrides this method.
  */
  virtual void     resolveLinearDependency(){}
  /*!
   A node should not be enriched by both H(x) and branch functions of the same
   crack. Remove H(x), just enriched by branch functions.
	Pure virtual method, only derived class CrackInterior overrides this method.
   New method coded in 2005-09-11.
   */
  virtual void     resolveConflictsInNodalEnrichment(){}
  /*!
   Update the geometry, virtual function that requires derived class (Crack,
   Hole,...) must implement this method.
	2005-09-02 : National day of Vietnam !!!
  */
  virtual void     UpdateMyGeometry(){;}
  /*!
   Update the enrichment. Just CrackTip has actual implementation since the update is
	just local around the crack tip.
   */
  virtual void     updateEnrichment(){;}
  /*!
   
   */
  virtual void     treatMeshGeoInteraction(Element* e){;}

  //=================== Methods used for debugging ==============================
  virtual void     printYourSelf(){;}

protected:
  GeometryEntity*      myGeometry;      //!< the original geometry
  int                  geometryID;      //!< 
  vector<EnrichmentFunction*>*  myEnrichFns ;  //!< the enrichment functions
  EnrichmentDetector*  myEnrDetector; //!< enrichment detector 
  vector<Element*>*    interactedElements ; //!< a list of elements interacting with me
} ;


#endif //_ENRICHMENTITEM_H_
