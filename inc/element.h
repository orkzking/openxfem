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
//   *********************
//   *** CLASS ELEMENT *** 
//   *********************
 

#ifndef _ELEMENT_H_
#define _ELEMENT_H_

#include "femcmpnn.h"
#include "domain.h"
#include "flotmtrx.h"
#include "geometry_base.h"
#include "delaunay.h"
#include <stdio.h>
#include <vector>
#include <list>
class TimeStep ; class Node ; class Material ; class GaussPoint ;
class FloatMatrix ; class FEInterpolation; 
class IntegrationRule ; class EnrichmentFunction; class EnrichmentItem;

//! A general finite element, base class of all element classes
/*!
 *  This abstract class is the most important class of the program. It is the
 *  superclass of all classes implementing finite elements (bar, shell, etc).
 *  An element is an attribute of a domain.
 * DESCRIPTION :
 *  The basic data of an element are the numbers of its 'numberOfNodes' nodes,
 *  stored in 'nodeArray', of its 'material', of its body loads (eg, the dead
 *  weight) stored in 'bodyLoadArray'. These data are obtained from the domain.
 *  The element also possesses 'numberOfGaussPoints' Gauss points, stored in
 *  'gaussPointArray'.
 *  The calculated data of an element are its 'massMatrix', its 'stiffnessMatrix',
 *  its 'locationArray'. Since the load vector is recalculated at every
 *  time step, it is not given the status of attribute.
 * TASKS :
 *  -defining itself :
 *    .typing itself (methods 'typed' and 'ofType'). When the domain creates
 *     an element, it actually creates a temporary instance of Element, then
 *     asks this element to transform itself into an element of the right type
 *     (PlaneStrain, Truss2D, etc) ;
 *    .obtaining its basic data : the element reads in the data file the num-
 *     ber of these objects, then obtains the data from the domain (methods
 *     'giveNode', 'giveMaterial',etc) ;
 *  -calculating its contribution to the problem :
 *    .calculating its stiffness matrix K, its mass Matrix M, its load vector
 *     f, its location array ;
 *    .calculating its contribution to the LHS and RHS of the linear system,
 *     using Static,Newmark,etc, formula. These contributions are usually
 *     combinations of M,K,f.
 *    .requesting from the LHS and the RHS of the linear system to assemble
 *     these contributions, by the means of the location array ;
 *  -performing end-of-step operations :
 *    .calculating the strains and stresses at its Gauss points ;
 *    .printing its output in the data file and updating itself ;
 */

class Element : public FEMComponent
{
   
   protected :
      size_t            numberOfNodes ;//!< number of nodes of one element
      IntArray*         nodeArray ;    //!< array contains nodes of element
      int               material ;     //!< material ID
      IntArray*         bodyLoadArray ;
      size_t            numberOfGaussPoints ; //!< number of Gauss points
		size_t            numOfGPsForJ_Integral ; //!< number of Gauss points 
      GaussPoint**      gaussPointArray ;          
      IntArray*         locationArray ;            //!< scatter matrix
	   FloatMatrix*      massMatrix ;               //!< element mass matrix
      FloatMatrix*      stiffnessMatrix ;          //!< stiffness matrix
		FloatMatrix*      constitutiveMatrix ; //!< constituive matrix of material
	   FEInterpolation*  standardFEInterpolation;   //!< Interpolation of standard part
	   FEInterpolation*  enrichmentFEInterpolation; //!< Interpolation of enrichement part
		size_t            numberOfIntegrationRules ; //!< number of used integration rules
	   IntegrationRule** quadratureRuleArray ;	   //!< Quadrature rule 
		std::list<EnrichmentItem*>* enrichmentItemListOfElem;//!< List of enrichment items acting on element  
		std::list<Element*>* neighbors;//!< neighboring elements of this element 
	   bool              checked; //!< Marker, used at method conflicts(Circle*) 
		bool              isUpdated ;       //!< marker for crack growth. 2005-09-03 
		bool              isMultiMaterial;  //!< one or multi-material element 
		std::vector<FloatMatrix*>*      constitutiveMatrices ; //!< constituive matrix of material
		std::vector<size_t> materialIDs ; //!< multi-materials element
   public :

      Element (int,Domain*) ;              //<! Constructor
		Element (){}                         //<! Default constructor
      virtual ~Element () ;                //<! Destructor

     	/** @name Assembly function 
	    *  Functions to perform the assembly procedure
		 */
      //@{
      void                  assembleYourselfAt (TimeStep*) ;
      void                  assembleLhsAt (TimeStep*) ;
      void                  assembleRhsAt (TimeStep*) ;
      IntArray*             giveLocationArray () ; 
		IntArray*             computeLocationArray () ;//@}
      
		/** @name Left-hand side
	    *  Functions to compute the left-hand side
		 */
      //@{
      FloatMatrix*          computeLhsAt (TimeStep*) ;
      FloatMatrix*          ComputeNewmarkLhsAt (TimeStep*) ;
      FloatMatrix*          ComputeStaticLhsAt (TimeStep*) ; //@}

      /** @name Right-hand side
	   *  Functions to compute the right-hand side
		*/
      //@{
      FloatArray*           computeRhsAt (TimeStep*, FloatArray*) ;
      FloatArray*           ComputeNewmarkRhsAt (TimeStep*, FloatArray*) ;
      FloatArray*           ComputeStaticRhsAt (TimeStep*, FloatArray*) ;//@}

       
		/** @name Mass and Stiffness matrix
	    *  Functions to compute the elementary matrices
		 */
      //@{
      FloatMatrix*          giveMassMatrix () ;
      FloatMatrix*          GiveStiffnessMatrix () ; // Purify, Matthias 2005
      FloatMatrix*          giveConstitutiveMatrix () ;
	   virtual FloatMatrix*  computeMassMatrix () ;
      virtual FloatMatrix*  ComputeConsistentMassMatrix () ;
		// NOT USED ANYMORE !!! @see ComputeTangentStiffness()
      virtual FloatMatrix*  computeStiffnessMatrix () ;
		virtual FloatMatrix*  computeConstitutiveMatrix (){ return NULL ;}   
		//virtual FloatMatrix*  computeConstitutiveMatrix (Material*){ return NULL ;}   
		//@}

      /** @name Load vector
	    *  Functions to compute the Load vector
		 */
      //@{
      FloatArray*           ComputeLoadVectorAt (TimeStep*) ;
      FloatArray*           ComputeBcLoadVectorAt (TimeStep*) ;
      virtual FloatArray*   ComputeBodyLoadVectorAt (TimeStep*) ;
	   virtual FloatArray*   ComputeResultingBodyForceAt (TimeStep*) ; //@}
  
		/** @name Strains and stresses
	    *  Functions to compute the strain and stresses
		 */
      //@{
      virtual FloatArray*   computeStrainVector (GaussPoint*,TimeStep*) ;
		void                  computeStrain (TimeStep* stepN);
	   FloatArray*			    computeStrainIncrement (GaussPoint* , FloatArray* );

      /** @name  Vector of nodal unknowns
	    *  Functions to compute the unknowns
		 */
      //@{
      virtual FloatArray*   ComputeVectorOf (char,TimeStep*) ;
      FloatArray*           ComputeVectorOfPrescribed (char,TimeStep*) ;
		virtual FloatArray*   ComputeVectorOfDisplacement (char,TimeStep*) ; // XFEM. NVP 07/05
      int                   computeNumberOfDofs () ; 
		size_t                computeNumberOfDisplacementDofs () ; // XFEM. NVP 2005
		//@}

		/** @name Interpolation, numerical integration
	    *  Functions to compute B matrix, N matrix, dV ...
		 */
      //@{
      virtual void          computeGaussPoints ()             {}
	   size_t                giveNumberOfGaussPoints(); // NVP 2005
		size_t                giveNumOfGPtsForJ_Integral(); 
		GaussPoint**          giveGaussPointArray();
      virtual FloatMatrix*  ComputeBmatrixAt (GaussPoint*);  
	   virtual FloatMatrix*  ComputeBuMatrixAt (GaussPoint*){ return NULL ;}
      virtual FloatMatrix*  ComputeNmatrixAt (GaussPoint*) { return NULL ;}
      virtual double        computeVolumeAround (GaussPoint*){ return NULL ;}
		//@}
		/*! 
	   Returns the finite interpolation which is used to
		approximate the continuous field such as displacement
	   */
	   virtual FEInterpolation*      giveFEInterpolation(){ return NULL ;}
		/*! 
	   Returns the finite interpolation which is used to
		approximate the discontinuous functions. 
		It's necessary to use linear interpolation for enriched
		approximation fields because of blending elements. For more
		details, see Chessa et all 2003
	   */
	   virtual FEInterpolation*      giveXFEInterpolation(){ return NULL ;}

      /** @name  Data management
	   *  Functions to manage the data of an element
		*/
      //@{
      Node*                 giveNode (int) ;
		IntArray*             giveNodeArray() ;
      Material*             giveMaterial () ;
      IntArray*             giveBodyLoadArray () ;
      void                  instanciateYourself () ;
		//@}

      // time step termination
      void                  printOutputAt (TimeStep*, FILE*, FILE*) ;
      void                  updateYourself () ;

      /** @name  Definition of an element
	   *  Functions to define an element
		*/
      //@{
      Element*              typed () ;
      Element*              ofType (char*) ;
      char*                 giveClassName (char* s) 
                                      { return strcpy(s,"Element") ;}
      int                   giveNumber ()  
				      { return FEMComponent::giveNumber() ;}
		//@}

	  //new-SC
	   FloatArray*			    ComputeInternalForces (FloatArray*);
	   FloatMatrix*          ComputeTangentStiffnessMatrix ();

	  /*! 
	   Returns true if this element is enriched and false otherwise.
	   */
	   bool                 isEnriched();
	  /*! 
	   If enriched with an enrichment item, then insert this
		item into the enrichmentItemListOfElem
	   */
	   void                 isEnrichedWith(EnrichmentItem*);
	  /*! 
	   Do the interaction between element and the enrichment items
	   */
     void                  treatGeoMeshInteraction();
	  /*!
	   Set the enrichment for the receiver. Virtual method since actual action depends 
		on the type of elements. For example, with linear elements, enrich all nodes. For 
		T6 elements but use enriched linear shape functions, just enrich three corner nodes.
	   */
	  virtual void          setEnrichmentForMyNodes(EnrichmentItem*){}
	  size_t                giveNumberOfNodes () const {return numberOfNodes;} ;
		/*!
		 Set up the integration rule will be used for numerical integration
		 for standard finite element, just use the continuous integration rule
		 for elements interacts with discontinuities, SplitGaussQuad should be used
		 */
	  void setUpIntegrationRule();
	  /*!
      Partition the element into subtriangles for numerical integration
      @returns the vector containing the subtriangles which defined in 
		local coordinate of the parent element.
      */
     virtual std::vector<DelaunayTriangle *>*  PartitionMySelf(){ return NULL ;}
	  /*!
	   compute the neighboring elements for a given element
	   */
	  std::list<Element*>*  ComputeNeighboringElements();
	  /*!
	   give the neighboring elements for a given element
	   */
	  std::list<Element*>*  giveNeighboringElements();
	  /*!
	   Print the neighboring elements of the receiver
		Useful for debugging purpose
	   */
	  void                 printMyNeighbors();
	  /*!
	   check if the receiver is belong to the Circle c or not.
	   */
	  bool                 in(Mu::Circle *c);
	  /*!
	   Operator < used to sort elements in List
	   */
	  bool  operator < (const Element& elem2)
	  {
		 return number < elem2.number ;
	  }
	  /*!
	   Computes the area of the element
	   */
	  virtual double   area() {return NULL ;}
	  /*!
	   Computes the area of portion above the enrichment item
	   */
	  virtual double  computeAreaAboveEnrItem(EnrichmentItem*){return NULL ;}
	  
	  /*!
	   Find elements conflicting with a Circle c with center belong to the receiver
		An element is called conflict with c if at least one node is inside c
	   */
	  std::list<Element*> conflicts(Mu::Circle *c);
	  /*!
	   Check the intersection between the receiver and the circle.
		2005-08-11, written so that will be used to compute the I integral.
	   */
	  bool                intersects(Mu::Circle *c);     
     /*!
	   Check the intersection between the receiver and the Segment.
		2005-08-29, code for crack growth simulation.
	   */
	  bool                intersects(Mu::Segment *seg); 

	  /*!
	   Check whether a \c Point p locates on edge of the receiver.
		2005-09-06
	   */
	  bool                IsOnEdge(Mu::Point* p);

	  /*!
	   Check if the \c Point p is within the receiver or not
		2005-09-06
	   */
	  bool                isWithinMe(Mu::Point* p);
     /*!
	   Check if the \c Point p is coincident to one of receiver's nodes or not
		2005-09-06
	   */
	  bool                isCoincidentToMyNode(Mu::Point* p);
	  /*!
	   Reset data member "checked" to false for the next checking ( used in Element::conflicts(Circle*)
	   */
	  void                clearChecked();
	  /*!
	   Just for debug purpose. Print all enrichment items acting on the receiver.
	   */
	  void                printMyEnrItems();
	  /*!
	   Return the center of the receiver. Used to compute the Heaviside function
		H(x) for nodes belong to the crack which is multi valued function(+1,-1),so
		compute the center of element at hand, if H(this center)=1 => H(node)=1.
		See Sukumar et Prevost Part I for more details.
	   */
	   Mu::Point*           giveMyCenter();
     /*!
	   Computes the stress points of element for plotting.
	   */
	  virtual void          computeStressPoints (){}
	  /*!
	   Defines the Gauss quadrature for elements used to compute the J integral.
		For split and tip-elements, discontinuous Gauss quadrature need to be used.
		For non-enriched elements, high order Gauss quadrature is used since we're integrating
		non-polynomial functions : sin(x), cos(x), sqrt(r) ...
		Implemented 2005-08-13.
	   */
	  virtual GaussPoint**  setGaussQuadForJ_Integral(){return NULL;}
	  /*!
	   Delete the list of EnrichmentItem of the receiver so that it can be reset
		when crack grows. 2005-08-30.
	   */
	  void    eraseListOfEnrichmentItems(){enrichmentItemListOfElem->clear() ;}
     /*!
	   Set isUpdated of the receiver to TRUE
		It means that this element need to recompute its stiffness matrix.
		Element changes their status due to crack growth !!! 
		2005-09-03
	   */
	  void    setStateOfElement(){isUpdated = true ; }
	  void    reinitializeStateOfElement(){isUpdated = false ; }
	  bool    isUpdatedElement(){return isUpdated ;}
     /*! 
	   If one element is enriched by both CrackInterior and MaterialInterface,
		then remove the MaterialInterface.
		Assuming that the CrackInterior is on the MaterialInterface !!!
		If it is not the case, no need this conflicts resolve.
	   */
	  void    resolveConflictsInEnrItems();
	  bool    isSplitElement(){return enrichmentItemListOfElem == NULL ? false : true ;}
	
	  /*! RB-SB-2004-10-29 begin*/
	  virtual void exportStressResultsToMatlab(std::string& theStringStress);
	  virtual void exportStrainResultsToMatlab(std::string& theStringStrain);
	  virtual void exportGaussPointsToMatlab(std::string& theString);
	  //SB2004-10-13
	  GaussPoint*  giveGaussPointNumber(int i)	;
	  FloatArray*  ComputeGlobalNodalCoordinates();

	  void		   reinitializeStiffnessMatrix();	// MP-purify 3.10.05
	  void         computeNodalLevelSets(EnrichmentItem*);
	  void         setMultiMaterial(){isMultiMaterial = true;}
	  bool         containsMultiMats(){return isMultiMaterial;} 
	  void         setMaterialIDs(std::vector<size_t> matIDs){materialIDs = matIDs;}
	  void         setMaterial(size_t matID){material = matID;}
	  std::vector<size_t> giveMatIDs(){return materialIDs;}
	  void         updateMaterialID();
	  std::list<EnrichmentItem*>* giveListOfEnrichmentItems(){return enrichmentItemListOfElem;}
} ;


#endif
