//   ***********************
//   *** CLASS INT ARRAY ***
//   ***********************


#ifndef _INTARRAY_H_
#define _INTARRAY_H_

#include "freestor.h"
#include "debug.def"
#include <stdlib.h>


class IntArray
  /*
  This class implements an array of integers.
  DESCRIPTION :
  An IntArray stores its coefficients in an array 'value' of size 'size'.
  TASKS :
  - storing and reterning coefficients (method 'at') ;
  - appending another IntArray to itself ;
  REMARK :
  see Remark 2 in file "floatarry.hxx".
  */
{
private :
  int   size ;
  int*  values ;

public :
  IntArray (int) ;                                   // constructor
  ~IntArray ()  { if(values) freeInt(values) ;}      // destructor

#     ifdef DEBUG
  int&    at (int) ;
#     else
  int&    at (int i)                     { return values[i-1] ;}
#     endif
  int&       operator[](const size_t i){return values[i];}
  void       checkBounds (int) ;
  IntArray*  followedBy (IntArray*) ;
  int        giveSize ()                 { return size ;}
  int        isEmpty ()                  { return size==0 ;}
  void       printYourself () ;
  int*       givePointer ()              {return values;} 
} ;


#endif


