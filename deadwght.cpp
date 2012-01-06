//   file DEADWGHT.CPP
 
#include "deadwght.h"
#include "element.h"
#include "timestep.h"
#include "loadtime.h"
#include "material.h"
#include "flotarry.h"



FloatArray*  DeadWeight :: ComputeForceOn (Element* elem, TimeStep* stepN)
// Returns the force (and not just the acceleration) induced at stepN by
// the receiver on the element 'elem'.
{
   double density,factor ;

   factor = this -> giveLoadTimeFunction() -> at(stepN->giveTime()) ;
   if (factor == 0.)
      return NULL ;
   else {
      density = elem -> giveMaterial() -> give('d') ;
      return  this -> giveComponentArray() -> Times(density*factor) ;}
}

