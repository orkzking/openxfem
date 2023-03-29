//   *******************************
//   *** CLASS CONSTANT FUNCTION ***
//   *******************************

 
#include "loadtime.h"


class ConstantFunction : public LoadTimeFunction
/*!
   This class implement a function  y = f(t)  that is constant in time.
 DESCRIPTION
   'value' is the constant value of the function. It is a pointer, rather
   than a number, so that its state (initialized or not) can be checked.
 */
{
   private :
      double*  value ; 

   public :
      ConstantFunction (int i,Domain* d) : LoadTimeFunction(i,d) {value=0 ;}
      ~ConstantFunction ()                       { delete value ;}

      double  at (double)             { return this->giveValue();}
      double  giveValue () ;
      void    instanciateYourself ()  { this->giveValue() ;}
} ;








