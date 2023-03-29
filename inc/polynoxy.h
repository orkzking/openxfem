//   ***************************
//   *** CLASS POLYNOMIAL XY ***
//   ***************************


#include "polyno.h"
#include "geometry_base.h"
#include "flotarry.h"


class PolynomialXY : public Polynomial
  /*!
  This class implements a polynomial with two variables X and Y.
  DESCRIPTION :
  The coefficients a_i of a polynomial are stored in the following sequence:
  2                2       3
  p(X,Y) = a1 + a2 X + a3 Y + a4 X  + a5 X Y + a6 Y  + a7 X  + ...
  REMARK :
  The degree of such polynomial is temporarily limited to 1 ; it can be
  easily increased whenever required.
  */
{
public :
  PolynomialXY (int) ;                 // constructor

  double  evaluatedAt (Mu::Point*) ;
  void    print() ;
} ;







