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

#ifndef _CRACKGROWTHDIRECTIONLAW_H_
#define _CRACKGROWTHDIRECTIONLAW_H_

#include "femcmpnn.h"
#include "cracktip.h"
#include <stdio.h>
#include <stdlib.h>


//! Law used to determine the direction of crack growth
/*!
 This class is an abstract class, the superclass of all classes that imple-
 ment the crack growth direction law : MaxHoopStress... 
*/
class CrackGrowthDirectionLaw : public FEMComponent
{
public:
  //<! Constructor
  CrackGrowthDirectionLaw(int n,Domain* aDomain):FEMComponent(n,aDomain){}   
  ~CrackGrowthDirectionLaw(){}             //!< destructor

  /*! 
   Computes the growth direction of the crack. Actual implementation is
   done in derived classes.
  */
  virtual double computeGrowthAngle(CrackTip*){return NULL ;}

  /** @name  Definition 
   *  Functions to define an CrackGrowthDirectionLaw
   */
  //@{
  CrackGrowthDirectionLaw* typed();
  CrackGrowthDirectionLaw* ofType(char*);
  char*               giveClassName (char* s)
  { return strcpy(s,"CrackGrowthDirectionLaw") ;}
  //@}

} ;


#endif // _CRACKGROWTHDIRECTIONLAW_H_