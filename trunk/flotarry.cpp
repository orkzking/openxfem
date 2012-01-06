//   file FLOTARRY.CPP

#include "flotarry.h"

#include "intarray.h"
#include "flotmtrx.h"
#include "mathfem.h"
#include <numeric>
#include <functional>
#include <math.h>
#include <string.h>
#include <assert.h>


FloatArray :: FloatArray (int n)
   // Constructor : creates an array of size n (filled with garbage).
{
   size = n ;
   if (size)
      values = allocDouble(size) ;
   else
      values = NULL ;
}

double FloatArray :: norm (double *p, int n) 
   // Returns the Euclidian norm of the n-size array p.
   // Remark: if norm**2 is needed, use rather function 'squaredNorm', which
   //         saves a square-root operation.
{
   return sqrt(dotProduct(p,p,n));
}

FloatArray*  FloatArray :: minus ()
   // Switches the sign of every coefficient of the receiver. Returns the
   // receiver.
{
   double* p;
   double  x;
   int     i;

   i = size;
   p = values;
   while (i--) {
      x    = - *p;
      *p++ = x; }
   return this;
}

FloatArray*  FloatArray :: Minus (FloatArray* aFloatArray)
{	//returns a new float array
  FloatArray *answer;
  if(this == NULL){
	 answer = new FloatArray(aFloatArray->giveSize());
	 for (int i = 1 ; i <= aFloatArray->giveSize() ; i++) 
		answer->at(i) = -aFloatArray->at(i);
  }
  else{
	 answer = new FloatArray(this->size);
	 for (int i=1; i<=this->size; i++) 
		answer->at(i) = this->at(i) - aFloatArray->at(i);
  }
  return answer;
}

FloatArray*  FloatArray :: Plus (FloatArray* aFloatArray)
{	//returns a new float array
	FloatArray *answer;
	answer = new FloatArray(this->size);
	for (int i=1; i<=this->size; i++) {
		answer->at(i) = this->at(i) + aFloatArray->at(i);
	}
	return answer;
}


FloatArray*  FloatArray :: minus (FloatArray* b)
   // Performs the operation a=a-b, where a stands for the receiver. 
   // - if size(a)=0, adjusts size(a) to size(b);
   // - if size(a)>size(b), OK. does as if the trailing coeffs of b were 0;
   // - if size(a)<size(b), OK if the trailing coeffs of b are 0.
   // Returns the receiver.
{
   double *p1,*p2;
   int    i;

   if (!b)
      return this;
   else if (! size) {                            
      size   = b->size;
      values = allocDouble(size); }
#  ifdef DEBUG
   else if (size < b->size) 
      for (i=size+1; i<=b->size; i++)
         if (b->at(i) != 0.) {
            printf ("\nError in DArray->minus(DArray)\n");
            exit(0); }
#  endif
      
   p1 = values;
   p2 = b->values;
	i  = min(size,b->size);
   while (i--)
      *p1++ -= *p2++;
   return this;
}

FloatArray*  FloatArray :: Extract (IntArray* loc)
    // Returns an array which contains only the coefficients of the receiver
    // that are pointed to by 'loc'.
{
    FloatArray *answer;
    int    pos,n,i;

#   ifdef DEBUG
   /*n = loc->giveSize();
   for (i=1; i<=n; i++) {
      int k = loc->at(i);
      if (k<1 || k>size) {
         printf ("\nError in DArray(%d)::Extract(loc) :\n",size);
         printf ("     loc(%d)=%d points out of array\n\n",i,k);
         exit(0); }} */
#   endif

    n      = loc->giveSize();
    answer = new FloatArray(n);
    for (i=1; i<=n; i++) {
		pos = loc->at(i);
		if (pos == 0) {//means bc!
			answer->values[i-1] = 0;
		} else {
		    answer->values[i-1] = this->at(pos);
		}
    }
    return answer;
}

FloatArray*  FloatArray :: plusProduct (FloatMatrix* A, FloatArray* y, double alpha)
   // Performs the operation  x = x + alpha A(transp).y , where x stands for
   // the receiver.
{
   double *pA,*pX,*pY,sum;
   int    n,m,i,j;

   n = A->giveNumberOfRows();
   m = A->giveNumberOfColumns();
   if (! size) {
      size   = m;
      values = allocDouble(m); }

#  ifdef DEBUG
   if (size != m) { 
      printf ("\nError in DArray::plusProduct :\n");
      printf ("    size(array)=%d not equal to nColumns(matrix)=%d\n\n",size,m);
      exit(0); }
   int p = y->giveSize();
   if (n != p) { 
      printf ("\nError in DArray::plusProduct :\n");
      printf ("      nRows(matrix)=%d not equal to size(array2)=%d\n\n",n,p);
      exit(0); }
   if (A->isDiagonal()) {
      printf ("\nError in DArray::plusProduct :\n");
      printf ("      matrix is diagonal: unexpected case\n\n");
      exit(0); }
#  endif

   pX = values;
   pA = A->givePointer();        // horrible but efficient
   i  = size;
   while (i--) {
      pY  = y->values;
      sum = 0.;
      j   = n;
      while (j--)
         sum += (*pA++) * (*pY++);
      *pX++ += alpha*sum; }
   
  return this;
}


FloatArray*  FloatArray :: add (FloatArray* b)
   // Performs the operation a=a+b, where a stands for the receiver. If the
   // receiver's size is 0, adjusts its size to that of b. Returns the
   // receiver.
{
   register int i ;
   double       *p1,*p2 ;

   if (!b || b->giveSize()==0)
      return this ;

   if (! size) {                              // null-sized array
      size   = b -> size ;
      values = allocDouble(size) ;}

#  ifdef DEBUG
      if (size != b->size) {                  // unmatching sizes
	 printf ("FloatArray dimension mismatch in a[%d]->add(b[%d])\n",
		  size,b->size) ;
	 exit(0) ;}
#  endif

   p1 = values ;
   p2 = b -> values ;
   i  = size ;
   while (i--)
      *p1++ += *p2++ ;
   return this ;
}


void  FloatArray :: assemble (FloatArray* fe, IntArray* loc)
// Assembles the array fe (typically, the load vector of a finite
// element) to the receiver, using loc as location array.
{
  int i,ii,n ;

#  ifdef DEBUG
  if ((n=fe->giveSize()) != loc->giveSize()) 
  {
	 printf ("dimensions of 'fe' and 'loc' mismatch \n") ;
	 assert(false) ; 
  }
  this -> checkSizeTowards(loc) ;
#  endif

  n = fe->giveSize() ;
  for (i = 1 ; i <= n ; i++) {
	 ii = loc->at(i) ;
	 if (ii)                            // if non 0 coefficient,
		this->at(ii) += fe->at(i) ;      // then assemble
  }         
}


#ifdef DEBUG
double&  FloatArray :: at (int i)
   // Returns the i-th coefficient of the receiver. Slow but safe.
{
   this -> checkBounds(i) ;
   return values[i-1] ;
}
#endif


#ifdef DEBUG
void  FloatArray :: checkBounds(int i)
   // Checks that the receiver's size is not smaller than 'i'.
{
   if (i<=0) {
      printf ("array error on index : %d <= 0 \n",i) ;
      assert(false) ; 
	}

   if (i>size) {
      printf ("array error on index : %d > %d \n",i,size) ;
      assert(false) ; 
	}
}
#endif


void  FloatArray :: checkSizeTowards (IntArray* loc)
   // Expands the receiver if loc points to coefficients beyond the size of
   // the receiver.
{
   int i,n,high ;

   high = 0 ;
   n    = loc -> giveSize() ;
   for (i=1 ; i<=n ; i++)
	  high = max(high,(loc->at(i))) ;
   if (high > size)                             // receiver must be expanded
      this -> growTo(high) ;
}


int  FloatArray :: containsOnlyZeroes ()
   // Returns True if all coefficients of the receiver are 0, else returns
   // False.
{
   register int i ;
   double       *p ;

   p = values ;
   i = size ;
   while (i--)
      if (*p++ != 0.)
	 return FALSE ;

   return TRUE ;
}


void  FloatArray :: growTo (int n)
   // Expands the receiver up to size n (n is assumed larger than 'size').
   // Initializes all new coefficients to zero.
{
   register int i ;
   double       *newValues,*p1,*p2 ;

#  ifdef DEBUG
    if (!n || n<=size) {
	   printf ("error in FloatArray(%d)::growTo(%d) \n",size,n) ;
	   assert(false) ;
	 }
#  endif

   newValues = allocDouble(n) ;

   p1 = values ;
   p2 = newValues ;
   i  = size ;
   while (i--)
      *p2++ = *p1++ ;

   if (values)
      freeDouble (values) ;
   values = newValues ;
   size   = n ;
}


FloatArray*  FloatArray :: negated ()
   // Switches the sign of every coefficient of the receiver. Returns the
   // receiver.
{
   register int i ;
   double       x ;
   double*      p ;

   i = size ;
   p = values ;
   while (i--) {
      x    = - *p ;
      *p++ = x ;}
   return this ;
}


void  FloatArray :: printYourself ()
   // Prints the receiver on screen.
{
   printf ("FloatArray of size : %d \n",size) ;
   for (int i=1 ; i<=size ; ++i)
      printf ("%.12f  ",this->at(i)) ;
   printf ("\n") ;
}


FloatArray*  FloatArray :: reinitialized ()
   // Returns the receiver with all coefficients set to 0.
{
   if (values)
      freeDouble (values) ;
   values = allocDouble(size) ;
   return this ;
}


FloatArray*  FloatArray :: rotatedWith (FloatMatrix* r, char mode)
   // Returns the receiver 'a' rotated according the change-of-base matrix r.
   // If mode = 't', the method performs the operation  a = r(transp) * a .
   // If mode = 'n', the method performs the operation  a = r * a .
{
   register int i ;
   double       *p1,*p2 ;
   FloatMatrix  *rot ;
   FloatArray   *rta ;

   if (mode == 't')
      rot = r -> GiveTransposition() ;
   else
      rot = r ;

   rta = rot -> Times(this) ;

   p1 = values ;
   p2 = rta -> values ;
   i  = size ;
   while (i--)
      *p1++ = *p2++ ;

   if (mode == 't')
      delete rot ;
   delete rta ;
   return this ;
}


FloatArray*  FloatArray :: times (double factor)
   // Multiplies every coefficient of the receiver by factor. Answers the
   // modified receiver.
{
   register int i ;
   double*      p ;

   p = values ;
   i = size ;
   while (i--)
      *(p++) *= factor ;
   return this ;
}


FloatArray*  FloatArray :: Times (double factor)
   // Returns a new array, whose components are those of the receicer, times
   // factor.
{
   register int i ;
   double       *p1,*p2 ;
   FloatArray*  answer ;

   answer = new FloatArray(size) ;
   p1     = values ;
   p2     = answer -> values ;
   i      = size ;
   while (i--)
      *p2++ = factor * (*p1++) ;
   return answer ;
}

double FloatArray :: transposedTimes (FloatArray* aFloatArray) {
	//returns a double (1xn * nx1)
	double answer = 0;
	for (int i=1; i<=this->size; i++) {
		answer+= this->at(i) * aFloatArray->at(i);
	}
	return answer;
}

FloatMatrix* FloatArray :: timesTransposed (FloatArray* aFloatArray) {
	//returns a float matrix (nx1 * 1xn)
	FloatMatrix* answer;
	answer = new FloatMatrix(this->size,this->size);
	for (int i=1; i<=this->size; i++) {
		for (int j=1; j<=this->size; j++) {
			answer->at(i,j) = this->at(i) * aFloatArray->at(j);
		}
	}
	return answer;
}

double FloatArray :: computeTensorialNorm()
//***************************************
// Without forgetting to double the s12 term!!!

{	double a11, a22, a33, a12, norm;

	a11 = this -> at(1);
	a22 = this -> at(2);
	a12 = this -> at(3);
	a33 = this -> at(4);

	norm = sqrt(a11*a11+a22*a22+a33*a33+2*a12*a12);
	return norm;
}

double FloatArray :: computeInvariantI1()
//***************************************
// Compute the invariant I1

{  double       I1;

   I1 = this->at(1)+this->at(2)+this->at(4);
   //this->at(3) = sigma12 !!

   return I1;
}

double FloatArray :: computeInvariantJ2()
//***************************************
// Compute the invariant J2 (0.5 * s  s  )
//								    ij ij

{  double   J2;
   double   p, a11, a22, a33, a12;

// Computation of the deviatoric part
   a11 = this -> at(1);
   a12 = this -> at(3);
   a22 = this -> at(2);
   a33 = this -> at(4);
   p = (a11+a22+a33)/3.;
   a11 -= p; // s11
   a22 -= p; // s22
   a33 -= p; // s33

// Computation of J2
   J2 = 0.5*(a11*a11+a22*a22+a33*a33+2.*a12*a12);

   return J2;
}

FloatArray* FloatArray :: computeDeviatoricPart()
//***********************************************
// Compute the deviatoric part of the second order tensor

{  double      p;
   FloatArray* answer;

   p = (this -> computeInvariantI1())/3.;
   answer = new FloatArray(this->size);
   answer -> at(1) = this->at(1) - p;
   answer -> at(2) = this->at(2) - p;
   answer -> at(3) = this->at(3);
   answer -> at(4) = this->at(4) - p;

   return answer;
}

FloatArray* FloatArray :: computeHydrostaticPart()
//***********************************************
// Compute the deviatoric part of the second order tensor

{  double      p;
   FloatArray* answer;

   answer = new FloatArray(this->size);

   p = (this -> computeInvariantI1())/3.;
   answer -> at(1) = p;
   answer -> at(2) = p;
   answer -> at(4) = p;

   return answer;
}

FloatArray* FloatArray :: computeDJ2DSigma () 
{
	FloatArray* answer;
	double a11, a12, a22, a33;
	double const aThird = 1.0f/3.0f;

    answer = new FloatArray(this->size);

	a11 = this->at(1);
	a12 = this->at(3);
	a22 = this->at(2);
	a33 = this->at(4);

	answer->at(1) = aThird * ((2*a11) - a22 - a33);
	answer->at(2) = aThird * ((2*a22) - a11 - a33);
	answer->at(3) = (2 * a12);
	answer->at(4) = aThird * ((2*a33) - a22 - a11);
	return answer;
}




