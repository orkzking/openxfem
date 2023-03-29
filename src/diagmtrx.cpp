//   file DIAGMTRX.CPP
 
#include "diagmtrx.h"
#include "flotarry.h"
#include <iostream>
#include <assert.h>


void  DiagonalMatrix :: checkBounds (size_t i,size_t j)
   // Checks that the receiver possesses a (i,j) coefficient.
{
	assert(i <= nRows) ;
	assert(j == i) ;
}


FloatMatrix*  DiagonalMatrix :: GiveCopy ()
   // Returns a copy of the receiver. For consistency with method giveCopy()
   // of class FloatMatrix, returns the answer as a FloatMatrix.
{
   DiagonalMatrix *answer ;
   double         *P1,*P2 ;
   int            i ;

   answer = new DiagonalMatrix(nRows) ;
   P1     = values ;
   P2     = answer->values ;
   i      = nRows ;
   while (i--)
      *P2++ = *P1++ ;
   return (FloatMatrix*)answer ;
}

FloatArray*  DiagonalMatrix :: Times (FloatArray* anArray)
   // Returns the product of the receiver and anArray.
{

   FloatArray* answer = new FloatArray(nRows) ;
	
   for (size_t i=0 ; i< nRows ; i++) {
	   (*answer)[i] = values[i] * (*anArray)[i];
   }
   return answer ;
}


void  DiagonalMatrix :: printYourself ()
   // Prints the receiver on screen.
{
	std::cout <<  "Diagonal matrix of size " << nRows << std::endl ;
   for (size_t i=0 ; i<nRows ; i++)
	   std::cout << values[i] << "  " << std::flush ;
  
	std::cout << std::endl ;
}


