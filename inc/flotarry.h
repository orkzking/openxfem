//   *************************
//   *** CLASS FLOAT ARRAY ***
//   *************************

/*
   This class implements an array of floating-point numbers.
 DESCRIPTION :
   A FloatArray stores its coefficients in an array 'values' of size 'size'.
 TASKS :
   - storing and returning a coefficient (method 'at') ;
   - expanding its size in order to store additional coefficients (method
     'growTo') ;
   - performing basic operations : summation, product, rotation, etc ;
   - assembling to itself another array, typically an elemental or nodal
     load vector (method 'assemble') ;
   - reading/writing its description on a given file.
 REMARKS :
   - for the sake of efficiency, the array 'values' is allocated using the
     C 'calloc' function rather than the 'new' operator.
   - method 'givePointer' is an encapsulation crime. It is used only for
     speeding up method 'dot' of class RowColumn.
*/ 

#ifndef _FLOTARRY_H_
#define _FLOTARRY_H_

#include "freestor.h"
#include "debug.def"
#include <stdlib.h>

class IntArray ; class FloatMatrix ;

class FloatArray
{
   enum { FALSE , TRUE } ;

   protected :
		int      size ;
		double*  values ;

   public :
		FloatArray ()      {size=0; values=NULL;}               // constructor
		FloatArray (int) ;                                       // constructor
		virtual ~FloatArray () {if(values) freeDouble(values);}  // destructor

#ifdef DEBUG
		double&      at (int i) ;
#else
		double&      at (int i)                { return values[i-1] ;}
#endif
		double&      operator[](const size_t i){return values[i];}
		FloatArray*  add (FloatArray*) ;
		void         assemble (FloatArray*,IntArray*) ;
		void         checkBounds (int) ;
		void         checkSizeTowards (IntArray*) ;
		int          containsOnlyZeroes () ;
		FloatArray*  GiveCopy ()               { return this->Times(1.) ;}
		double*      givePointer ()            { return values ;} // see above
		int          giveSize ()               { return size ;}
		void         growTo (int) ;
		int          isNotEmpty ()             { return (size != 0) ;}
		FloatArray*  negated () ;
		void         printYourself () ;
		FloatArray*  reinitialized() ;
		FloatArray*  rotatedWith (FloatMatrix*,char) ;
		FloatArray*  times (double) ;
		FloatArray*  Times (double) ;

		friend double  dotProduct (double*,double*,int) ;

		double   giveNorm ()                {return norm(values,size);}
		double   norm (double *p, int n) ;

		FloatArray* minus();
		FloatArray* minus (FloatArray*);
		FloatArray* Minus (FloatArray*);
		FloatArray* Plus (FloatArray*);

		FloatArray* Extract (IntArray*);

		FloatArray* plusProduct (FloatMatrix* , FloatArray* , double);

		// for Dep
		double transposedTimes (FloatArray*);
		FloatMatrix* timesTransposed (FloatArray*);

		// new from Milan J.
		double computeInvariantI1() ;
		double computeInvariantJ2() ;
		double computeTensorialNorm() ;
		FloatArray* computeDeviatoricPart() ;
		FloatArray* computeHydrostaticPart() ;
		FloatArray* computeDJ2DSigma () ;

} ;


#endif
