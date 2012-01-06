//   *********************
//   *** CLASS SKYLINE ***
//   *********************


#ifndef _SKYLINE_H_
#define _SKYLINE_H_

#include "lhs.h"
#include "column.h"
#include "flotarry.h"
#include <stdio.h>

class Domain ; class FloatMatrix ;
//! A Skyline matrix 
/*!
   This class implements a symmetric matrix stored in a compact (skyline)
   form. A skyline is usually an attribute of the linear system.
 DESCRIPTION :
   Attribute 'size' is number of columns of the skyline. Attribute 'columns'
   is a list of 'size' columns, each of them of any size. Attribute 'isFac-
   torized' is True if the skyline is already in U(T).D.U factorized form,
   else it is False.
 TASKS :
   - storing or returning a coefficient (method 'at').
   - assembling to itself a matrix, typically and elemental stiffness matrix
     (method 'assemble') ; this may include enlarging the profile (method
     'growTo') ;
   - when associated with an array b, solving the system Ax=b (methods 'fac-
     torized','forwardReductionWith',etc) ;
   - when associated with a domain, shaping its own profile (method 'carve-
     Yourself') ;
   - resetting to zero all of its coefficients (method 'reinitialized').
*/
class Skyline : public LHS
{
   enum { FALSE , TRUE } ;
   protected :
      size_t         size ;
      Column**       columns ;
      bool           isFactorized ;

   public :
      Skyline () ;                         // constructor
      ~Skyline () ;                        // destructor

      FloatMatrix*  AsFloatMatrix () ;
	   double&       at (int i,int j)     {return this->giveColumn(j)->at(i);}
      void          assemble (FloatMatrix*,IntArray*) ;
      FloatArray*   backSubstitutionWith (FloatArray*) ;
      void          carveYourselfFor (Domain*) ;
      void          checkSizeTowards (IntArray*) ;
      void          createColumn (int,int) ;
      Skyline*      diagonalScalingWith (FloatArray*) ;
      Skyline*      factorized () ;
      Skyline*      forwardReductionWith (FloatArray*) ;
      Column*       giveColumn (int j)               { return columns[j-1] ;}
      int           giveNumberOfColumns ()           { return size ;}
		void          growTo (int) ;
      void          growTo (size_t) ;
      void          printYourself () ;
      Skyline*      reinitialized () ;
} ;

#endif // _SKYLINE_H_
