//   **********************
//   *** CLASS MATHUTIL ***
//   **********************
 
#ifndef MathUtil_H
#define MathUtil_H

#include <vector>
#include "flotarry.h"

using namespace std;

class FloatMatrix;

class MathUtil
/*!
   This class implements a mathematical utilities (Sherman-Morrison, ...)
 */
{
   public :
	   MathUtil () { ; } ;
      ~MathUtil () { ; } ;

	  FloatMatrix* giveShermanMorrisonInverse (FloatMatrix* , FloatArray*, double);
	  

} ;
#endif









