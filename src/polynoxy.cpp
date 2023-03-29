//   file POLYNOXY.CXX
 
#include "polynoxy.h"
#include "geometry_base.h"
#include <stdlib.h>


PolynomialXY :: PolynomialXY (int n)
   // Constructor. Creates an (X,Y) polynomial of degree n.
{
   degree       = n ;
   coefficients = new FloatArray(++n * ++n / 2) ;
}


double  PolynomialXY :: evaluatedAt (Mu::Point * aPoint)
{
	if (degree == 0)
		return (*coefficients)[0] ;
	else if (degree == 1)
		return (*coefficients)[0] + (*coefficients)[1]*aPoint->x + (*coefficients)[2]*aPoint->x ;
	else
		assert(false) ;
}


void  PolynomialXY :: print()
   // Prints the receiver on screen.
{
	std::cout << "Polynomial(X,Y) of degree " << degree << std::endl ;
	
	for (size_t i = 0 ; i < coefficients->giveSize() ; i++)
		std::cout << (*coefficients)[i] << std::endl ;
}








