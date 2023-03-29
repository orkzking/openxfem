//   file MATHUTIL.CPP

#include "mathutil.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include <stdlib.h>
#include <math.h>

class FloatArray; class FloatMatrix;

FloatMatrix*  MathUtil :: giveShermanMorrisonInverse
                                       (FloatMatrix* Kinv, FloatArray* a, double x)
{	//this implements the Sherman-Morrison method to compute
	//the inverse of a matrix which can be put in the form:
	//
	//			    T
	//M = K + x * aa
	//
	//K must be diagonal
	//
	FloatMatrix *M, *aNumerator,*temp1,*temp3;
	FloatArray	*temp2;
	
	double aFactor;
	
	temp3		= a->timesTransposed(a);
	temp1		= Kinv->Times(temp3);
	aNumerator	= (temp1->Times(Kinv))->times(x);
	temp2		= Kinv->Times(a);
	aFactor = 1 / (1 + (x * (a->transposedTimes(temp2))));

	M = Kinv->Minus(aNumerator->times(aFactor));
	delete aNumerator;
	delete temp1;
	delete temp2;
	delete temp3;
	return M;
}



