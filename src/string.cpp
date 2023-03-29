//   file STRING.CXX

#include "string.h"
#include "mathfem.h"                /*  for min()  */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>                  /*  for isInteger(s)  */


char*  concatenate (char* s,int i,char* target)
   // A non member function. Returns target, the concatenation of s and i.
{
   char s2[8] ;

   target[0] = '\0' ;
   strcat (target,s) ;
   strcat (target,itos(i,s2)) ;

   return target ;
}


char*  ctos (char c, char* s)
   // A non-member function. Converts the character c into a string s.
   // Returns s.
{
   s[0] = c ;
   s[1] = '\0' ;
   return s ;
}


bool  isInteger (char* s)
   // A non-member function. Returns True if s contains an integer, else
   // returns False. Leading and trailing blanks are not allowed.
   // This function always returns True if the DEBUG option is not used.
{
#  ifdef DEBUG
      int  len,i ;
      char c ;

      len = strlen(s) ;
      if (! len)
	 return false ;

      c = s[0] ;
      if (c=='+' || c=='-')
	 i=1 ;
      else if (isdigit(c))
	 i=0 ;
      else
	 return false ;

      for ( ; i<len ; i++)
	 if (! isdigit(s[i]))
	    return false ;
      return true ;
#  else
      return true ;
#  endif
}



bool  isNumber (char* s)
   // A non member function. Returns True is s contains a number, else
   // returns False. Leading and trailing blanks are not allowed.
   // s should have a sign (optional), an integer part, a dot and a decimal
   // part (both optional), and an exponent part (optional). The exponent
   // part has 'e' or 'E', a sign (optional) and digits.
   // This finction always returns True if the DEBUG option is not used.
{
#  ifdef DEBUG
      int  len,i ;
      char c ;

      len = strlen(s) ;
      if (! len)
	 return false ;

      c = s[0] ;                             // look for a sign
      if (c=='+' || c=='-')
	 i=1 ;
      else if (isdigit(c))
	 i=0 ;
      else
	 return false ;


      for (i=i+1 ; i<len ; i++)              // process integer part
	 if (! isdigit(c=s[i])) {
	    if (c=='.')
	       goto decimal ;
	    else if (c=='e' || c=='E')
	       goto exponent ;
	    else
	       return false ;}
      return true ;

      decimal :
      for (i=i+1 ; i<len ; i++)
	 if (! isdigit(c=s[i])) {
	    if (c=='e' || c=='E')
	       goto exponent ;
	    else
	       return false ;}
      return true ;

      exponent :
      c = s[i+1] ;
      if (c=='+'|| c=='-')
	 i++ ;
      if (i == len-1)
	 return false ;
      for (i=i+1 ; i<len ; i++)
	 if (! isdigit(s[i]))
	    return false ;
      return true ;
#  else
      return true ;
#  endif
}


char*  itos (int n, char* s)
   // A non member function, copied from Kernighan & Ritchie, p64. Implemen-
   // ted for portability (the PC function itoa(int,char*,int) is not ANSI
   // C). Returns s, the string conversion of n.
{
   int  i,j,sign ;
   char c ;

   if (n == 1)
      return strcpy(s,"1") ;
   if (n == 2)
      return strcpy(s,"2") ;

   if ((sign=n) < 0)
      n = -n ;
   i = 0 ;
   do
      s[i++] = n % 10 + '0' ;
   while ((n/=10) > 0) ;
   if (sign < 0)
      s[i++] = '-' ;
   s[i] = '\0' ;

   // reverse 's'
   j = strlen(s) - 1 ;
   for (i=0 ; i<j ; i++,j--) {
      c    = s[i] ;
      s[i] = s[j] ;
      s[j] = c ;}

   return s ;
}

