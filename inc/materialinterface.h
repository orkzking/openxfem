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

Any feedback is welcome. Emails : nvphu80@yahoo.com, ...

*************************************************************************************/

#ifndef _MATERIALINTERFACE_H_
#define _MATERIALINTERFACE_H_

#include "enrichmentitem.h"

//! Material Interface
/*!
* This class implements a material interface.
* Each \c Material Interface is associated to two \c Materials, the part below this
* interface is assumed belonging to the ith material, and the upper part is of the
* (i+1)th material.
* The corresponding enrichment function is the absolute of signed distance, which
* is implemented in class \c AbsoluteSignedDistance.
* (@see paper "Arbitrary discontinuities in finite element" of Prof. TB)
* The input file is :
*   1 MaterialInterface Mat 2 1 2 geometry 1 enrFuncs 1 1 enrichScheme 3
*/
class MaterialInterface:public EnrichmentItem
{
public:
  MaterialInterface(int,Domain*);   //!< constructor
  ~MaterialInterface(){}             //!< destructor

  /*!
   Gives the list of two materials stored in materials
	Read from the input file for the first time.
   */
  std::vector<size_t>  giveMaterialIDs();

  /*!
  Time step termination
  */
  void     printOutputAt (TimeStep*, FILE*, FILE*){}
  /*!
  Since the \c Material Interface is fixed, this function is empty
  */
  void     updateEnrichment(){;}
  void     printYourSelf();
  void     treatMeshGeoInteraction(Element* e);

private:
  std::vector<size_t> materialIDs;//!< associated materials
} ;


#endif //_MATERIALINTERFACE_H_
