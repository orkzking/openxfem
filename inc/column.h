//   ********************
//   *** CLASS COLUMN ***
//   ********************


#ifndef _COLUMN_H_
#define _COLUMN_H_

#include "flotarry.h"
class FloatMatrix ; class Skyline ;


class Column : public FloatArray
/*!
   This class implements a column in a matrix stored in segmented form
   (symmetric skyline). A column is a particular kind of FloatArray.
 DESCRIPTION :
   A column n stores in 'values' its 'size' coefficients, upwards :
   .values[0]      = diagonal coefficient A(n,n)
   .values[1]      = off-diagonal coefficient A(n-1,n)
   .values[size-1] = highest non-0 coefficient of the n-th column.
 TASKS :
   Those inherited from FloatArray.
 */
{
private :
	int       number ;
	Skyline*  matrix ;

public :
	Column (int n,int size,Skyline* m) : FloatArray(size)
					       { number=n ; matrix=m ; }
	~Column () {}

	double&  at (int i)       { return values[number-i] ;} 
	void     checkSizeTowards (IntArray*,int) ;
	double   dot (Column*,int,int) ;
	Column*  GiveCopy () ;
	int      giveHighestRow ()   { return number - size + 1 ;}
} ;

#endif // _COLUMN_H_
