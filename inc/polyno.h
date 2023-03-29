//   ************************
//   *** CLASS POLYNOMIAL ***
//   ************************


#ifndef polyno_h

#include "flotarry.h"
#include "geometry_base.h"
#include <stdio.h>


class Polynomial
  /*!
  This abstract class is the superclass of the classes that implement
  polynomials : polynomials with one variable (X), with two variables
  (X,Y), etc.
  DESCRIPTION :
  The degree of a polynomial is stored in attribute 'degree ; for example,
  P(X,Y) = 2 + 4 X + 5 XY has degree 2.
  The coefficients of the polynomial are stored in 'coefficients' ; their
  sequence depends on the the polynomial's type (PolynomialXY, etc).
  TASKS :
  The task of a polynomial is to return its value at a given point (method
  'evaluatedAt').
  */
{
protected :
  int          degree ;
  FloatArray*  coefficients ;

public :
  Polynomial ()                       { }
  virtual ~Polynomial ()              { delete coefficients;}

  double&         at (int i)          { return (*coefficients)[i-1] ;}
  virtual double  evaluatedAt (Mu::Point*)  = 0 ; // pure virtual method
  virtual void    printYourself ()           = 0 ;
} ;

#define polyno_h
#endif








