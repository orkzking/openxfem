//   *****************************
//   *** CLASS DIAGONAL MATRIX ***
//   *****************************
 

#ifndef _DIAGMTRX_H_
#define _DIAGMTRX_H_


#include "flotmtrx.h"
#include "flotarry.h"

class DiagonalMatrix : public FloatMatrix
/*!
   This class implements a square matrix containing non zero coefficients
   only on the diagonal. Diagonal matrices are typically used as mass
   matrices of finite elements.
 DESCRIPTION :
   A diagonal matrix of size n has nRows=n and nColumns=1. 'values' stores
   only the n diagonal coefficients.
 TASKS :
   Similar to a standard FloatMatrix.
 */
{

   public :
      DiagonalMatrix (int n) : FloatMatrix(n,1)    { }

#     ifdef DEBUG
       double&      at (size_t i,size_t j)
			    {this->checkBounds(i,j) ; return values[i-1];}
#     else
      double&      at (int i,int j) const          { return values[i-1] ;}
#     endif
      void          checkBounds (size_t,size_t) ;
      FloatMatrix*  GiveCopy () ;
      int           isDiagonal () const            { return true ;}
      void          printYourself () ;
	   FloatArray*   Times (FloatArray*) ;
} ;

#endif // _DIAGMTRX_H_

