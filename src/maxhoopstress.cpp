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

#include "maxhoopstress.h"
#include "cracktip.h"
#include <math.h>



MaxHoopStress :: MaxHoopStress(int n, Domain *d)
:CrackGrowthDirectionLaw(n,d)
{
}

MaxHoopStress :: ~MaxHoopStress()
{
}

double MaxHoopStress :: computeGrowthAngle(CrackTip *aTip)
// *******************************************************
// compute the direction growth of crack 
// using maximum hoop stress criteria
{
	 
	 double KI  = aTip->giveK_i();
	 double KII = aTip->giveK_ii();

	 double ratio = KII/KI  ;
	 double A = -2.0*(ratio); 
	 double B = 1.0 + sqrt(1.0 + 8.0 * ratio * ratio);

	 double answer = 2.0 * atan2(A,B);

	 return answer;
}


