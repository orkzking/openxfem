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
 
#ifndef _ENRICHMENTFUNCTION_H_
#define _ENRICHMENTFUNCTION_H_

#include "enrichmentitem.h"
#include "flotarry.h"
#include "geometry_base.h"
#include "element.h"
#include "gausspnt.h"
#include "node.h"
#include <vector>

//! Abstract base class for definition of enrichment functions
/*!
 *  This class is an abstract class, the superclass of all classes that implement 
 *  the enrichment functions: discontinuos function, asymptotic functions...\n 
 *
 *  Main tasks is to evaluate its values and gradients at a given point. Since one
 *  enrichment function can be used to model numerous enrichment items, it must have
 *  a list containing all enrichment items being enriched by itself.
 */
class EnrichmentFunction : public FEMComponent
{
public:
	EnrichmentFunction(Domain*,int);    //!< Constructor
	~EnrichmentFunction();              //!< Destructor

	/*!
	 Computes the value of enrichment function at a point
	 For instance, computes this function at \c GaussPoint gp
	 to calculate the enriched stiffness matrix
	 */
	virtual double      EvaluateYourSelfAt(GaussPoint* gp){ return NULL ;}
	
	/*!
	 Computes the value of enrichment function at a \c Node.
	 Used for the nodes belong to the crack => multi-valued enrichment function.
	 */
	virtual double      EvaluateYourSelfAt(Element*,Node*){ return NULL ;}
	/*!
	 Computes the derivatives of enrichment function at a point
	 For instance, computes this function at \c Gausspointgp
	 to calculate the enriched stiffness matrix
	 */
	virtual FloatArray* EvaluateYourGradAt(GaussPoint* gp){ return NULL ;}
   /*!
	 Set the EnrichmentItem associated with the receiver.
	 */
	void                setMyEnrichmentItem(EnrichmentItem*);
	void                findActiveEnrichmentItem(EnrichmentItem*);
	void                printMyEnrichmentItem();

	/** @name  Definition 
	 *  Functions to define an enrichment function
	 */
   //@{
	EnrichmentFunction* typed();
	EnrichmentFunction* ofType(char*);
	char*               giveClassName (char* s)
                         { return strcpy(s,"EnrichmentFunction") ;}//@}

protected:
	int number; 
	std::vector<EnrichmentItem*> *myEnrItems; //!< enrichment items modeled by the receiver
	EnrichmentItem               *activeEnrItem;//!< current enrichment item modeled by this function
} ;


#endif //_ENRICHMENTFUNCTION_H_
