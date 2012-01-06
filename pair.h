//   ******************
//   *** CLASS PAIR ***
//   ******************


#include <stdio.h>

class Pair
/*! 
   This class implements key/value associations, e.g., the key 'E' and its
   value 210000. A pair is used as an entry in a dictionary.
 DESCRIPTION :
   A pair has three components : its name (a character), its value (a number),
   a pointer to the next pair in the dictionary.
 TASKS : 
   - Returning its key, or its value, or the next pair ;
   - Appending another pair to itself.
*/
{
   private :
      char    key ;
      double  value ;
      Pair*   next ;

   public :
      Pair (char k,double v)    { key=k ; value=v ; next=NULL ;}
      ~Pair ()                  { }

      void     append (Pair* p) { next = p ;}
      char     giveKey ()       { return key ;}
      Pair*    giveNext ()      { return next ;}
      double&  giveValue()      { return value ;}
      void     printYourself () { printf("   Pair (%c,%f)\n",key,value);}
} ;








