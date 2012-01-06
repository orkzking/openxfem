//   file COLUMN.CXX

#include "column.h"
#include "skyline.h"
#include "intarray.h"
#include "flotmtrx.h"
#include "mathfem.h"
#include "freestor.h"
#include "debug.def"
#include <numeric>
#include <math.h>


void  Column :: checkSizeTowards (IntArray* loc, int matrixIsDiagonal)
   // Enlarges the receiver if 'loc' points to coefficients outside of the
   // receiver. This test if trivial if the matched matrix is diagonal, since
   // a column should always contin at least the diagonal coefficient.
{
   int i,coeff,highRow ;

   if (matrixIsDiagonal)
      return ;

   highRow = 32000 ;
   for (i=1 ; i<=loc->giveSize() ; i++) {   // find the most off-diagonal
      coeff = loc->at(i) ;                  //       coefficient in 'loc'
      if (coeff)
		  highRow = min(highRow,coeff) ;}

   if (highRow < this->giveHighestRow())
      this -> growTo(number-highRow+1) ;
}


double  Column :: dot (Column* b, int start, int stop)
   // Returns the scalar product of the receiver by 'b', from index 'start'
   // (included) to index 'stop' (included).
{
   double *P1,*P2 ;
   int    i ;

#  ifdef DEBUG
     if (start<=number-size || start<=b->number-b->size || stop>number ||
							   stop>b->number) {
       printf ("error : cannot make dot product from %d to %d\n",start,stop);
       printf ("        of Column %d (size %d) and Column %d (size %d)\n",
		number,size,b->number,b->size) ;
       exit(0) ;}
#  endif

   P1 = values + number - stop ;             // P1 points to a(stop)
   P2 = b->values + b->number - stop  ;      // P2 points to b(stop)
   i  = stop - start + 1 ;

   return  dotProduct(P1,P2,i) ;
}

Column*  Column :: GiveCopy ()
   // Returns a copy of the receiver.
{
   Column        *answer ;
   double        *p1,*p2 ;
   register int  i ;

   answer = new Column(number,size,matrix) ;
   p1 = answer->values ;
   p2 = values ;
   i  = size ;
   while (i--)
      *p1++ = *p2++ ;
   return answer ;
}