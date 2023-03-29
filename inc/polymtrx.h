//   *******************************
//   *** CLASS POLYNOMIAL MATRIX ***
//   *******************************
 
#ifndef _POLYMTRX_H_
#define _POLYMTRX_H

#include "matrix.h"
#include "flotarry.h"
#include "geometry_base.h"

class Polynomial ; class FloatMatrix ;

class PolynomialMatrix:public Matrix
/*!
   This class implements a matrix which contains polynomials. These matrices
   are typically used as jacobian matrices of finite elements.
 DESCRIPTION :
   The matrix stores its nRows*nColumns polynomials column by column.
 TASKS :
   - storing and returing coefficients (i.e., polynomials), in method 'at' ;
   - evaluating itself at a given point (method 'EvaluatedAt'). This prod-
     uces a matrix of the same size, containing numbers.
*/
{
   protected :
      Polynomial**  values ;

   public :
      PolynomialMatrix (int,int) ;           // constructor
      ~PolynomialMatrix () ;                 // destructor

#     ifdef DEBUG
	Polynomial*&     at (int,int) ;
#     else
	Polynomial*&     at (int i,int j) {return values[(j-1)*nRows+i-1] ;}
#     endif
	FloatMatrix*       EvaluatedAt (Mu::Point*) ;
      PolynomialMatrix*  initialized () ;
      void               printYourself () ;
} ;
#endif  // _POLYMTRX_H