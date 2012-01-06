//   file CONSTANT.CPP

#include "constant.h"


double  ConstantFunction :: giveValue ()
   // Returns the constant value of the receiver. Reads 'value' in the data
   // file if it hasn't been done yet.
{
   if (! value) {
      value    = new double[1] ;
      value[0] = this -> read("f(t)") ;}

   return *value ;
}