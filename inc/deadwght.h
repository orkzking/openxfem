//   *************************
//   *** CLASS DEAD WEIGHT ***
//   *************************
 

#include "bodyload.h"


class DeadWeight : public BodyLoad
/*!
   This class implements a gravity-like load.
 DESCRIPTION
   The attribute 'componentArray' contains the components of an acceleration
   'a', expected (but not required) to be downwards vertical.
 TASK
   returning the body force  rho*a  acting on a given element.
 */
{
   public :
      DeadWeight (int i,Domain* d) : BodyLoad(i,d) {}         // constructor

      FloatArray*  ComputeForceOn (Element*,TimeStep*) ;
      void         saveAttributesIn (FILE*) ;
} ;









