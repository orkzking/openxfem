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



#include "crackgrowthincrementlaw.h"
#include "fixedincrement.h"
#include "assert.h"
#include <stdio.h>
#include <iostream>

CrackGrowthIncrementLaw*  CrackGrowthIncrementLaw :: typed ()
// **********************************************************
// Returns a new CrackGrowthDirectionLaw,which has the same number than the receiver,
// but is typed (MaxHoopStress,...)
{
   CrackGrowthIncrementLaw *newIncrementLaw ;
   char     type[32] ;

   this -> readString("class",type) ;
   newIncrementLaw = this -> ofType(type) ;

   return newIncrementLaw ;
}

CrackGrowthIncrementLaw*  CrackGrowthIncrementLaw :: ofType (char* aClass)
// ***********************************************************************
// Returns a new CrackGrowthDirectionLaw, which has the same number than the receiver,
// but belongs to a class (MaxHoopStress,...)
{
   CrackGrowthIncrementLaw  *newIncrementLaw ;

   if (! strcmp(aClass,"FixedIncrement"))
      newIncrementLaw = new FixedIncrement(number,domain) ;
	//else if (! strcmp(aClass,"MultipleOfTipElementSize"))
   //   newIncrementLaw = new MultipleOfTipElementSize(number,domain) ;
   else 
	{
      printf ("%s : unknown Crack increment law \n",aClass) ;
      assert(false) ;
	}

   return newIncrementLaw ;
}

