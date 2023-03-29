//   file FLOTMTRX.CPP

#include "flotmtrx.h"
#include "diagmtrx.h"
#include "flotarry.h"
#include "mathfem.h"
#include <numeric>
#include <iostream>
#include <assert.h>


#ifdef DEBUG
double&  FloatMatrix :: at (int i,int j)
// Returns the coefficient (i,j) of the receiver. Safer but slower than
// the inline version of method 'at'.
{
  this->checkBounds(i,j) ;
  return values[(j-1)*nRows + i - 1] ;
}
#endif


FloatMatrix*  FloatMatrix :: GiveCopy ()
// Creates and returns a copy of the receiver.
{
  FloatMatrix *answer ;
  double      *P1,*P2 ;
  int         i ;

  answer = new FloatMatrix(nRows,nColumns) ;
  P1 = answer -> values ;
  P2 = values ;
  i  = nRows * nColumns ;
  while (i--)
	 *P1++ = *P2++ ;
  return answer ;
}


double  FloatMatrix :: giveDeterminant ()
// Returns the determinant of the receiver.
{
#  ifdef DEBUG
  if (! this->isSquare()) {
	 printf ("error : cannot compute determinant of a %d by %d matrix \n"
		,nRows,nColumns) ;
	 exit(0) ;}
#  endif

  if (nRows == 1)
	 return  values[0] ;
  else if (nRows == 2)
	 return  (values[0]*values[3] - values[1]*values[2]) ;
  else if (nRows == 3)   
	 return  (values[0]*values[4]*values[8]+values[1]*values[5]*values[6]
	 +values[2]*values[3]*values[7]-values[1]*values[3]*values[8]
	 -values[2]*values[4]*values[6]-values[0]*values[5]*values[7]) ;
  else {
	 printf ("sorry : cannot inverse %d by %d matrices \n",nRows,nColumns) ;
	 exit(0) ;
	 return 0.0 ; // SC
  }
}


FloatMatrix*  FloatMatrix :: GiveInverse ()
// Returns a new matrix, the inverse of the receiver. (implemented only
// for 1x1 and 2x2 matrices)
{
  FloatMatrix* answer ;
  double       det ;
  double*      p ;

#  ifdef DEBUG
  if (! this->isSquare()) {
	 printf ("error : cannot inverse a %d by %d matrix ! \n",
		nRows,nColumns);
	 exit(0) ;}
#  endif

  answer = new FloatMatrix(nRows,nRows) ;
  p      = answer->values ;

  if (nRows == 1)
	 p[0] = 1. / values[0] ;
  else if (nRows == 2) {
	 det  = values[0]*values[3] - values[1]*values[2] ;
	 p[0] =  values[3] / det ;
	 p[1] = -values[1] / det ;
	 p[2] = -values[2] / det ;
	 p[3] =  values[0] / det ;}
  else {
	 printf ("error : cannot inverse a %d by %d matrix ! \n",nRows,nColumns);
	 exit(0) ;}

  return answer ;
}


FloatMatrix*  FloatMatrix :: GiveTransposition ()
// Returns a new matrix, the transposition of the receiver.
{
  FloatMatrix  * answer = new FloatMatrix(nColumns,nRows) ;
  for (size_t i=0 ; i<nRows ; i++)
	 for (size_t j=0 ; j<nColumns ; j++)
		answer->at(j+1,i+1) = this->at(i+1,j+1) ;
  return answer ;
}



FloatMatrix*  FloatMatrix :: Lumped ()
// Returns a new diagonal matrix, which is the lumped receiver : all
// coefficients on a column are concentrated on the diagonal.
{

  DiagonalMatrix *answer = new DiagonalMatrix(nRows) ;
  double *p      = values ;

  for (size_t j=0 ; j< nColumns ; j++) 
  {
	 double sum = 0. ;
	 register int i   = nRows ;
	 while (i--)
		sum += *p++ ;

	 answer -> at(j+1,j+1) = sum ;
  }

  return answer ;
}


FloatMatrix*  FloatMatrix :: plus (FloatMatrix* aMatrix)
// Adds aMatrix to the receiver. If the receiver has a null size,
// adjusts its size to that of aMatrix. Returns the modified receiver.
{
  register int i ;
  int      n,m ;
  double   *P1,*P2 ;

  if (aMatrix -> isDiagonal())
	 return  this->plusDiagonalMatrix((DiagonalMatrix*)aMatrix) ;

  n = aMatrix -> nRows ;
  m = aMatrix -> nColumns ;
  if (nRows*nColumns == 0) {
	 if (values)
		freeDouble (values) ;
	 nRows    = n ;
	 nColumns = m ;
	 i        = n * m ;
	 values   = allocDouble(i) ;
	 P1       = values ;
	 P2       = aMatrix->values ;
	 while (i--)
		*P1++ = *P2++ ;}
  else {
#     ifdef DEBUG
	 if (n-nRows || m-nColumns) {
		printf ("dimensions mismatch : r1,c1,r2,c2 : %d %d %d %d\n",
		  nRows,n,nColumns,m) ;
		exit(0) ;}
#     endif

	 P1 = values ;
	 P2 = aMatrix->values ;
	 i  = n * m ;
	 while (i--)
		*P1++ += *P2++ ;}

  return this ;
}


FloatMatrix*  FloatMatrix :: plusDiagonalMatrix (DiagonalMatrix* aMatrix)
{
  size_t n = aMatrix -> giveNumberOfRows() ;

  if (nRows*nColumns == 0) 
  {
	 if (values)
		freeDouble (values) ;

	 nRows    = n ;
	 nColumns = n ;
	 values   = allocDouble(n*n) ;
  }

  assert(n == nRows) ;

  for (size_t i=0 ; i<nRows ; i++)
	 this->at(i+1,i+1) += aMatrix->at(i+1,i+1) ;

  return this ;
}


void  FloatMatrix :: plusProduct (FloatMatrix* a, FloatMatrix* b, double dV)
// Adds to the receiver the product  a(transposed).b dV .
// If the receiver has a null size, it is expanded.
// This method assumes that both the receiver and the product above are
// symmetric matrices, and therefore computes only the upper half of the
// receiver ; the lower half is not modified. Other advantage : it does
// not compute the transposition of matrix a.
{
  
  double coeff ;
  double *P1,*P2 ;

  if (nRows*nColumns == 0) {
	 nRows  = nColumns = a->nColumns ;
	 values = allocDouble(nRows*nColumns) ;
  }

  int p = a->nRows ;

  for (size_t i = 1 ; i <= nRows ; i++)
	 for (size_t j = i ; j <= nColumns ; j++) 
	 {
		P1 = &(a->at(1,i)) ;
		P2 = &(b->at(1,j)) ;
		coeff = dotProduct(P1,P2,p) ;
		this->at(i,j) += coeff * dV ;
	 }
}


void  FloatMatrix :: printYourself ()
// Prints the receiver on screen.
{
  std::cout << "FloatMatrix with dimensions :" << nRows << ", " << nColumns << std::endl ;

  if (nRows<=10 && nColumns<=8)
	 for (size_t i = 0 ; i < nRows ; i++) 
	 {
		for (size_t j = 0 ; j < nColumns && j < 10 ; j++)
		  std::cout << this->at(i+1,j+1) << "  " << std::flush ;

		std::cout << std::endl ;
	 }
  else 
	 std::cout << "   large matrix : coefficients not printed" << std::endl ;
}


FloatMatrix*  FloatMatrix :: rotatedWith (FloatMatrix* r)
// Returns the receiver 'a' rotated according the change-of-base matrix r.
// The method performs the operation  a = r(transp) * a * r .
{
  register int  i ;
  double        *p1,*p2 ;
  FloatMatrix   *rt,*rta,*rtar ;

  rt   = r -> GiveTransposition() ;          //  r(transp)
  rta  = rt -> Times(this) ;                 //  r(transp) . a
  rtar = rta -> Times(r) ;                   //  r(transp) . a . r

  p1 = values ;
  p2 = rtar -> values ;
  i  = nRows*nColumns ; 
  while (i--)
	 *p1++ = *p2++ ;                         // copy rtar into the receiver

  delete rt ;
  delete rta ;
  delete rtar ;
  return this ;
}


FloatMatrix*  FloatMatrix :: symmetrized ()
// Initializes the lower half of the receiver to the upper half.
{
  assert(nRows == nColumns) ;

  for (size_t i=1 ; i < nRows ; i++)
	 for (size_t j=0 ; j<i ; j++)
		this->at(i+1,j+1) = this->at(j+1,i+1) ;

  return this ;
}


FloatMatrix*  FloatMatrix :: times (double factor)
// Multiplies every coefficient of the receiver by factor. Answers the
// modified receiver.
{
  register int i ;
  double*      p ;

  p = values ;
  i = nRows * nColumns ;
  while (i--)
	 *p++ *= factor ;
  return this ;
}


FloatArray*  FloatMatrix :: Times (FloatArray* anArray)
// Returns the product of the receiver and anArray.
{
  assert(nColumns == anArray->giveSize()) ;

  FloatArray* answer = new FloatArray(nRows) ;
  for (size_t i=0 ; i<nRows ; i++) 
  {
	 double sum = 0. ;
	 for (size_t j=0 ; j<nColumns ; j++)
		sum += this->at(i+1,j+1) * (*anArray)[j] ;
	 (*answer)[i] = sum ;
  }

  return answer ;
}

FloatMatrix*  FloatMatrix :: Minus (FloatMatrix* aMatrix)
{
  FloatMatrix* answer = new FloatMatrix(nRows,nColumns) ;
  for (size_t i=0 ; i<nRows ; i++) 
  {
	 for (size_t j=0 ; j< nColumns ; j++)
	 {
		answer->at(i+1,j+1) = this->at(i+1,j+1) - (aMatrix->at(i+1,j+1));
	 }
  }
  return answer ;
}

FloatMatrix*  FloatMatrix :: Plus (FloatMatrix* aMatrix)
{
  FloatMatrix* answer = new FloatMatrix(nRows,nColumns) ;
  for (size_t i=0 ; i<nRows ; i++) 
  {
	 for (size_t j=0 ; j<nColumns ; j++)	
	 {
		answer->at(i+1,j+1) = this->at(i+1,j+1) + (aMatrix->at(i+1,j+1));
	 }
  }
  return answer ;
}


FloatMatrix*  FloatMatrix :: Times (FloatMatrix* aMatrix)
// Returns the product of the receiver and aMatrix. Easier to use than
// operator * . 
{
  assert(nColumns == aMatrix->nRows) ;

  size_t p = aMatrix -> nColumns ;
  FloatMatrix* answer = new FloatMatrix(nRows,p) ;
  for (size_t i=0 ; i<nRows ; i++)
  {
	 for (size_t j=0 ; j<p ; j++) 
	 {
		double coeff = 0. ;

		for (size_t k=0 ; k<nColumns ; k++)
		  coeff += this->at(i+1,k+1) * aMatrix->at(k+1,j+1) ;

		answer->at(i+1,j+1) = coeff ;
	 }
  }
  return answer ;
}

FloatMatrix* FloatMatrix::FollowedBy(FloatMatrix* matrixFollowingReceiver)
// ***********************************************************************
// NVP - 2006
// Used to compute the B matrix for enriched finite elements
// For enriched elements : B = [B_fem B_enriched]
{
  FloatMatrix* answer ;

  if(values == NULL)  // if the receiver is NULL
  {
	 if(matrixFollowingReceiver->values != NULL)  
		answer = matrixFollowingReceiver->GiveCopy() ; // Purify, 14-10-2005
	 else
		answer = new FloatMatrix();
  }
  else 
  {
	 if(matrixFollowingReceiver->values == NULL)
		answer = new FloatMatrix((*this));
	 else 
	 {
		assert(matrixFollowingReceiver->nRows == nRows) ;

		size_t Rows = matrixFollowingReceiver->nRows;
		size_t Cols = matrixFollowingReceiver->nColumns;

		answer = new FloatMatrix(nRows,nColumns + Cols);

		for (size_t i = 0 ; i < nRows ; i++)
		  for (size_t j = 0 ; j < nColumns ; j++)
			 answer->at(i+1,j+1) = this->at(i+1,j+1); 

		for (size_t i = 0 ; i < Rows ; i++)
		  for (size_t j = 0 ; j < Cols ; j++)
		  {
			 double temp = matrixFollowingReceiver->at(i+1,j+1);
			 answer->at(i+1,nColumns+j+1) = temp ;
		  }
	 }
  }

  return answer;

}

FloatMatrix* FloatMatrix :: SwapTwoRows(size_t row1,size_t row2)
// Swap two rows of the matrix
{
  assert(row2 < this->giveNumberOfRows()) ;

  for(size_t j=0 ; j< this->giveNumberOfColumns() ; j++)
  {
	 double temp = this->at(row1,j+1);
	 this->at(row1,j+1) = this->at(row2,j+1);
	 this->at(row2,j+1) = temp;
  }

  return this;
}

FloatArray* FloatMatrix :: giveRow(size_t n)
// returns the ith row of the receiver
{
  FloatArray *answer = new FloatArray(this->giveNumberOfColumns());

  for(size_t j = 0 ; j < this->giveNumberOfColumns() ; j++)
  {
	 (*answer)[j] = this->at(n,j+1);
  }

  return answer ;
}




