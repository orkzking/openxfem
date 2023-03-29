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

#ifndef _MAXHOOPSTRESS_H_
#define _MAXHOOPSTRESS_H_

#include "crackgrowthdirectionlaw.h"
#include "cracktip.h"

/*!
A class for determining the growth direction using the 
maximum hoop stress criteria of Erdogan and Sih.
*/
class MaxHoopStress : public CrackGrowthDirectionLaw
{
public:
  MaxHoopStress(int, Domain*); //!< Constructor
  ~MaxHoopStress();            //!< Destructor
  /*!
  Compute the direction of crack growth
  using the maximum hoop stress criteria of Erdogan and Sih.
  */
  double computeGrowthAngle(CrackTip*);
} ;


#endif // _MAXHOOPSTRESS_H_


