// This file contains some functions used in the finite element
// application.
// ref : Lippman p 104
// ! friend functions may not be declared inline (C++\Views p 17)
 

inline int min (int i,int j)
    { return (i<=j ? i : j) ; } 
// 
inline int max (int i,int j)
    { return (i>=j ? i : j) ; }

inline double  dotProduct (double* P1, double* P2, register int i)
   // A non-member function. Returns the dot product of the first 'i' coef-
   // ficienst of the two arrays P1 and P2. This method applies to many
   // situations, eg row*column products with matrices.
{
   double answer ;

   answer = 0. ;
   while (i--)
      answer += *P1++ * *P2++ ;
   return answer ;
}

/// Returns the signum of given value (if value is < 0 returns -1, otherwise returns 1)
inline double sgn (double i){ return (i< 0. ? -1.: 1.); }

/**
 Solves cubic equation for real roots.
 @param a,b,c,d - coefficients of equation in form: \f$ax^3 + bx^2 + cx + d = 0\f$
 @param r1,r2,r3 - roots (only first num roots is valid)
 @param num      - number of roots resolved
 Copyright: this functions is copied from OOFEM written by Borek Patzak
*/
void cubic (double a, double b, double c, double d, double*r1, double *r2, double *r3, int *num);





