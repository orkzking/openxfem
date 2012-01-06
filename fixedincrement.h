//   ***********************************
//   ***   CLASS FIXEDINCREMENT      ***
//   ***********************************

 
#ifndef _FIXEDINCREMENT_H_
#define _FIXEDINCREMENT_H_

#include "crackgrowthincrementlaw.h"


/*!
 Read from the input file the crack increment length
 used for crack growth simulation.
 */
class FixedIncrement : public CrackGrowthIncrementLaw
{

public:

  FixedIncrement(int,Domain*); //!< Constructor
  ~FixedIncrement();           //!< Destructor

  /*! 
   Computes the increment length of the crack.
   */
  double computeIncrementLength();

private:

  double   delta ; //!< the crack increment lenght
} ;

#endif //_FIXEDINCREMENT_H_
