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

#ifndef _CRACKGROWTHINCREMENTLAW_H_
#define _CRACKGROWTHINCREMENTLAW_H_

#include "femcmpnn.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>


//! Law used to determine the crack increment length.
/*!
 This class is an abstract class, the superclass of all classes that imple-
 ment the crack growth increment law : Paris law, fixed increment and adaptive 
 increment( equal to tip-element size).
*/
class CrackGrowthIncrementLaw : public FEMComponent
{
public:
 
  CrackGrowthIncrementLaw(int n,Domain* aDomain):FEMComponent(n,aDomain){} //<! Constructor   
  ~CrackGrowthIncrementLaw(){}             //!< destructor

  /*! 
  Computes the increment length of the crack. Actual implement is
  done in derived classes.
  */
  virtual double computeIncrementLength(){return NULL ;}

  /** @name  Definition 
   *  Functions to define an CrackGrowthIncrementLaw
   */
  //@{
  CrackGrowthIncrementLaw* typed();
  CrackGrowthIncrementLaw* ofType(char*);
  char*               giveClassName (char* s)
  { return strcpy(s,"CrackGrowthIncrementLaw") ;}
  //@}

} ;

#endif //_CRACKGROWTHINCREMENTLAW_H_
