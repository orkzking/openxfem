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

#include "hole.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

using namespace std;

Hole :: Hole(int n,Domain* aDomain)
              : EnrichmentItem(n,aDomain)
{
	 ;
}

void Hole::printYourSelf()
{
  cout << "I am a Hole " << endl;
}
