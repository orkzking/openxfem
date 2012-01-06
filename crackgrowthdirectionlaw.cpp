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



#include "crackgrowthdirectionlaw.h"
#include "maxhoopstress.h"
#include "assert.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

CrackGrowthDirectionLaw*  CrackGrowthDirectionLaw :: typed ()
// **********************************************************
// Returns a new CrackGrowthDirectionLaw,which has the same number than the receiver,
// but is typed (MaxHoopStress,...)
{
   CrackGrowthDirectionLaw *newDirectionLaw ;
   char     type[32] ;

   this -> readString("class",type) ;
   newDirectionLaw = this -> ofType(type) ;

   return newDirectionLaw ;
}

CrackGrowthDirectionLaw*  CrackGrowthDirectionLaw :: ofType (char* aClass)
// ***********************************************************************
// Returns a new CrackGrowthDirectionLaw, which has the same number than the receiver,
// but belongs to a class (MaxHoopStress,...)
{
   CrackGrowthDirectionLaw  *newDirectionLaw ;

   if (! strcmp(aClass,"MaxHoopStress"))
      newDirectionLaw = new MaxHoopStress(number,domain) ;
   else 
	{
      printf ("%s : unknown Crack growth direction law \n",aClass) ;
      assert(false) ;
	}

   return newDirectionLaw ;
}

