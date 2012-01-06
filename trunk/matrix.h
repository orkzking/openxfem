//   ********************
//   *** CLASS MATRIX ***
//   ********************

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <cstdlib>

class FloatMatrix ; class IntArray ;

//! A general class for matrices
/*!
This abstract class is the superclass of the class that implement
matrices (FloatMatrix, PolynomialMatrix,...).
DESCRIPTION :
A matrix is characterized by its number of rows and columns.
TASKS :
Its tasks are defined in the subclasses.
*/
class Matrix
{
protected:
  size_t  nRows ;
  size_t  nColumns ;

public:
  Matrix ()             { }                         // default constructors
  Matrix (int n,int m)  { nRows=n ; nColumns=m ;}   // constructors
  ~Matrix ()            { }                         // destructor

  void          checkBounds (int,int) ;
  size_t        giveNumberOfRows ()    const  { return nRows ;}
  size_t        giveNumberOfColumns () const  { return nColumns ;}
  bool          isSquare () const             { return (nRows==nColumns) ;}
};

#endif // _MATRIX_H
