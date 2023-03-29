//   **************************
//   *** CLASS STRESS ARRAY ***
//   **************************


#ifndef _STRESSARRAY_H_
#define _STRESSARRAY_H_

#include "stressarray.h"
#include "flotarry.h"
#include <stdlib.h>

class StressArray : public FloatArray  
/*!
   This class implements a stress array
*/
{	public :
		StressArray () ;	// constructor
		StressArray (int n):FloatArray(n) {;}; // constructor
		~StressArray () ;	// destructor
		
} ;

#endif //_STRESSARRAY_H_
