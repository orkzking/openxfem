//   file NODLOAD.CPP
 
#include "nodload.h"
#include "loadtime.h"
#include "flotarry.h"
#include "timestep.h"

FloatArray*  NodalLoad :: ComputeValueAt (TimeStep* stepN)
// Returns an array, the force induced at stepN by the receiver.
{
   FloatArray* force ;
   
   double factor = this -> giveLoadTimeFunction() -> at(stepN->giveTime()) ;
   if (factor == 0.)
      return NULL ;
   else 
	{
      force  = this -> giveComponentArray() ;
      return  force -> Times(factor) ;
	}
}


