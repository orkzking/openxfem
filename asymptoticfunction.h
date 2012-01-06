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
 
#ifndef _ASYMPTOTICFUNCTION_H_
#define _ASYMPTOTICFUNCTION_H_

#include "enrichmentfunction.h"
#include "domain.h"


//! Abstract class for asymptotic functions
/*!
 This class is an abstract class, the superclass of all classes that imple-
 ment the asymptotic functions, for instance, CrackAsymptotic.
 This is a virtual base class, it means that it just provide the interface
 No object of this class can be created !
 */
class AsymptoticFunction:public EnrichmentFunction
{
public:
   //!< Constructor
  AsymptoticFunction(Domain* aDomain,int n): EnrichmentFunction(aDomain,n){}
  virtual ~AsymptoticFunction(){}      //!< Destructor
} ;


#endif //_ASYMPTOTICFUNCTION_H_
