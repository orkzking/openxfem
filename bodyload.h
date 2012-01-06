//   ***********************
//   *** CLASS BODY LOAD ***
//   ***********************

 
#ifndef _BODYLOAD_H_
#define _BODYLOAD_H_

#include "load.h"
class Element ; class TimeStep ;


//! Body Load
/*!
 This abstract class is the superclass of all loads (e.g., the dead weight)
 that act on the whole volume of finite elements. A body load is usually
 attribute of many elements.
 */
class BodyLoad : public Load
{
   public :
      BodyLoad (int i,Domain* d) : Load (i,d) {}        // constructor

      virtual FloatArray*  ComputeForceOn (Element*,TimeStep*) = 0 ;
} ;


#endif








