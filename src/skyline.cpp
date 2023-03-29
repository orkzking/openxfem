//   file SKYLINE.CPP

#include "skyline.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "domain.h"
#include "element.h"
#include "timinteg.h"
#include "mathfem.h"
#include "debug.def"
#include "verbose.def"
#include <math.h>
#include <stdlib.h>

enum  { NOT_DIAGONAL , DIAGONAL } ;


Skyline :: Skyline ()
// Constructor. Creates a skyline of size 0.
{
  size         = 0 ;
  columns      = NULL ;
  isFactorized = FALSE ;
}


Skyline :: ~Skyline ()
// Destructor.
{
  int j ;

  if (size)
  {
	 for (j=0 ; j<size ; j++)
		delete (columns[j]) ;
	 delete [] columns ;
  }
}


FloatMatrix*  Skyline :: AsFloatMatrix ()
// Returns a matrix, the receiver in a full storage form. This is useful
// for debugging and printings.
{
  FloatMatrix* answer ;
  Column*      columnJ ;
  int          i,j,start ;
  double       coeff ;

  answer = new FloatMatrix(size,size) ;
  for ( j = 1 ; j <= size ; j++)
  {
	 columnJ = this->giveColumn(j) ;
	 start   = columnJ->giveHighestRow() ;
	 for (i = start ; i <= j ; i++)
	 {
		coeff = columnJ->at(i) ;
		answer->at(i,j) = coeff ;
		answer->at(j,i) = coeff ;
	 }
  }
  return answer ;
}


void  Skyline :: assemble (FloatMatrix* mat, IntArray* loc)
// Assembles the elemental matrix 'mat' to the receiver, using 'loc' as a
// location array. The values in ke corresponding to a zero coefficient
// in loc are not assembled.
{
  int     i,j,ii,jj,dim ;
  Column* columnJJ ;

#  ifdef DEBUG
  dim = mat -> giveNumberOfRows() ;
  if (dim != loc->giveSize()) {
	 printf ("error : dimension of 'mat' and 'loc' mismatch \n") ;
	 exit(0) ;
  }
  this -> checkSizeTowards(loc) ;
#  endif

  dim = mat -> giveNumberOfRows() ;

  if (mat -> isDiagonal()) {               // mat is a diagonal matrix
	 for (j=1 ; j<=dim ; j++) {
		jj = loc->at(j) ;
		if (jj) {
		  columnJJ = this->giveColumn(jj) ;
#           ifdef DEBUG
		  columnJJ -> checkSizeTowards(loc,DIAGONAL) ;
#           endif
		  columnJJ->at(jj) += mat->at(j,j) ;
		}
	 }
  }

  else {                                   // mat is a full matrix
	 for (j = 1 ; j <= dim ; j++) {
		jj = loc->at(j) ;
		if (jj) {
		  columnJJ = this->giveColumn(jj) ;
#           ifdef DEBUG
		  columnJJ -> checkSizeTowards(loc,NOT_DIAGONAL) ;
#           endif
		  for (i = 1 ; i <= dim ; i++) {
			 ii = loc->at(i) ;
			 if (ii && ii<=jj)
				columnJJ->at(ii) += mat->at(i,j) ;
		  }
		}
	 }
  }
}


FloatArray*  Skyline :: backSubstitutionWith (FloatArray* y)
// Returns the solution x of the system U.x = y , where U is the receiver.
// nota : x overwrites y
{
  int     i,j,start ;
  double  yJ ;
  Column* columnJ ;

  for (j=size ; j>0 ; j--) {
	 columnJ = this->giveColumn(j) ;
	 yJ      = y->at(j) ;
	 start   = columnJ->giveHighestRow() ;
	 for (i=start ; i<j ; i++)
		y->at(i) -= columnJ->at(i) * yJ ;}

  return y ;
}


void  Skyline :: carveYourselfFor (Domain* aDomain)
// Instanciates the profile of the receiver and initilalizes all coeffi-
// cients to zero.
{
  IntArray *LocArray ;
  int      n,nElem,i,ii,test ;

  nElem = aDomain -> giveNumberOfElements() ;
  test  = aDomain -> giveTimeIntegrationScheme() -> hasDiagonalLhs() ;
  for (n=1 ; n<=nElem ; n++) {
	 LocArray = aDomain -> giveElement(n) -> giveLocationArray() ;
	 this -> checkSizeTowards(LocArray) ;        // gain new columns
	 for (i=1 ; i <= LocArray->giveSize() ; i++) {
		ii = LocArray->at(i) ;
		if (ii)
		  this -> giveColumn(ii) -> checkSizeTowards(LocArray,test) ;}}
}


void  Skyline :: checkSizeTowards (IntArray*loc)
// Increases the number of columns of the receiver if 'loc' points to
// out-of-range columns.
{
  int i,dim,maxCol ;

  maxCol = 0 ;
  dim    = loc->giveSize() ;
  for (i=1 ; i<=dim ; i++)
	 maxCol = std::max (maxCol,loc->at(i)) ;

  if (maxCol > size)                       // enlarge the matrix
	 this -> growTo(maxCol) ;
}


void  Skyline :: createColumn (int j,int s)
// Construct a new j-th column, of size 's'.
{
  columns[j-1] = new Column(j,s,this) ;
}


Skyline*  Skyline :: diagonalScalingWith (FloatArray* y)
// Scales y by the diagonal of the receiver. Returns the receiver.
{
  int j ;

  for (j=1 ; j<=size ; j++) {
#     ifdef DEBUG
    std :: cout << j << "   " << this->at(j,j) << "\n";
	 if (fabs(this->at(j,j)) < 0.000000000001) {
		printf ("error in Skyline::diagScaling : pivot %d is zero \n",j) ;
		exit(0) ;}
#     endif
	 y->at(j) /= this->at(j,j) ;}
  return this ;
}


Skyline*  Skyline :: factorized ()
// Returns the receiver in  U(transp).D.U  Crout factorization form.
{
  int    i,j,highRow,start ;
  double sum,temp ;
  Column *columnI,*columnJ ;

  if (isFactorized)
	 return this ;

  // for every column, process :
  for (j = 2 ; j <= size ; j++)
  {
	 // first, the off-diagonal terms
	 columnJ = columns[j-1] ;
	 highRow = columnJ->giveHighestRow() ;
	 for (i = highRow+1 ; i < j ; i++)
	 {
		columnI = columns[i-1] ;
		start   = std::max (highRow,columnI->giveHighestRow()) ;
		columnJ->at(i) -= columnI->dot(columnJ,start,i-1) ;
	 }
	 // next, the diagonal term
	 sum = 0. ;
	 for (i = highRow ; i < j ; i++)
	 {
		columnI = columns[i-1] ;
		temp    = columnJ->at(i) ;
		columnJ->at(i) = temp / columnI->at(i) ;
		sum           += temp * columnJ->at(i) ;
	 }
	 columnJ->at(j) -= sum ;
  }

  isFactorized = TRUE ;

  return this ;
}


Skyline*  Skyline :: forwardReductionWith (FloatArray* b)
// Computes y, the solution of the system  U(transp).y = b , where U is
// the receiver (assumed to be in a factorized form). Returns the receiver.
// nota : y overwrites b
{
  int     i,j,start ;
  double  sum ;
  Column* columnJ ;

#  ifdef DEBUG
  if (! isFactorized) {
	 printf ("error : a skyline must be factorized before reductions\n");
	 exit(0) ;}
#  endif

  if (! size)
	 return 0 ;                                        // null size system
  if (b->giveSize() != size) {
#     ifdef VERBOSE
	 printf ("too small LHS : enlarge it \n") ;
#     endif
	 b -> growTo(size) ;}

  for (j=1 ; j<=size ; j++)
  {
	 columnJ = this -> giveColumn(j) ;
	 sum     = 0. ;
	 start   = columnJ -> giveHighestRow() ;
	 for( i = start ; i < j ; i++)
		sum += columnJ->at(i) * b->at(i) ;
	 b->at(j) -= sum ;
  }

  return this ;
}


void  Skyline :: growTo (int m)
// Increases to 'm' the number of columns of the receiver.
{
  int      j ;
  Column** newColumns ;

  if (m <= size) {
	 printf ("Error : skyline [%d] cannot grow to %d\n",size,m) ;
	 exit(0) ;}

  newColumns = new Column* [m] ;                // copy old columns
  for (j=0 ; j<size ; j++)
	 newColumns[j] = columns[j] ;
  if (columns)
	 delete [] columns ;
  columns = newColumns ;

  for (j=size ; j<m ; j++)                      // create new columns
	 this -> createColumn(j+1,1) ;

  size = m ;
}


void  Skyline :: printYourself ()
// Prints the receiver on screen.
{
  FloatMatrix* copy ;

  copy = this -> AsFloatMatrix() ;
  copy -> printYourself() ;
  delete copy ;
}


Skyline*  Skyline :: reinitialized ()
// Returns the receiver with all coefficients set to zero.
{
  int j ;

  for (j=1 ; j<=size ; j++)
	 this -> giveColumn(j) -> reinitialized() ;
  isFactorized = FALSE ;

  return this ;
}

