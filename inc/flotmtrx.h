//   **************************
//   *** CLASS FLOAT MATRIX ***
//   **************************


#ifndef _FLOTMTRX_H_
#define _FLOTMTRX_H_

#include "matrix.h"
#include "debug.def"
#include "freestor.h"
#include <stdio.h>
#include <stdlib.h>

class FloatArray ; class DiagonalMatrix ;



class FloatMatrix : public Matrix
/*!
   This class implements a matrix containing floating-point numbers.
 DESCRIPTION :
   The matrix stores its nRows*nColumns coefficients in 'values'. The coef-
   ficients are stored column by column, like in Fortran.
 TASKS :
   - storing and retrieveing a coefficient (method 'at') ;
   - performing standard operations : addition, multiplication, transposition,
     inversion, lumping, rotation, etc ;
 */
{
   enum { FALSE , TRUE } ;

   protected :
      double*  values ;

   public:
      FloatMatrix (int n,int m):Matrix(n,m) {values = allocDouble(n*m);}
		FloatMatrix () : Matrix(0,0)          {values = NULL ;}
      virtual ~FloatMatrix ()               {if(values) freeDouble(values);}

#ifdef DEBUG
      virtual double&   at (int,int) ;
#else
      virtual double&   at (int i,int j)    {return values[(j-1)*nRows+i-1];}
#endif
      virtual FloatMatrix*  GiveCopy () ;
      double        giveDeterminant () ;
      FloatMatrix*  GiveInverse () ;
      FloatMatrix*  GiveTransposition () ;
      virtual int   isDiagonal ()const      {return FALSE ;}
      FloatMatrix*  Lumped () ;
      FloatMatrix*  Over (double f)         {return this->Times(1./f) ;}
      FloatMatrix*  plus (FloatMatrix*) ;
      FloatMatrix*  plusDiagonalMatrix (DiagonalMatrix*) ;
      void          plusProduct (FloatMatrix*,FloatMatrix*,double) ;
      virtual void  printYourself () ;
      FloatMatrix*  rotatedWith (FloatMatrix*) ;
      FloatMatrix*  symmetrized () ;
      FloatMatrix*  times (double f) ;
      FloatMatrix*  Times (double f)     {return this->GiveCopy()->times(f);}
      FloatArray*   Times (FloatArray*) ;
      FloatMatrix*  Times (FloatMatrix*) ;
	   double*       givePointer ()     {return values;}

	  //SC
	  FloatMatrix*  Minus (FloatMatrix*) ;
	  FloatMatrix*  Plus (FloatMatrix*) ;
     // Implemented during coding X-FEM	
	  FloatMatrix* FollowedBy(FloatMatrix*);
	  FloatMatrix* SwapTwoRows(size_t,size_t);
	  FloatArray*  giveRow(size_t i);
} ;


#endif // _FLOTMTRX_H_

