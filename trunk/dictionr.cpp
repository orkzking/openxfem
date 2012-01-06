//   file DICTIONR.CPP
 
#include "dictionr.h"
#include <stdlib.h>


Dictionary :: ~Dictionary ()
   // Destructor.
{
   Pair* Next ;

   while (first) {
      Next = first -> giveNext() ;
      delete first ;
      first = Next ;}
}


Pair*  Dictionary :: add (char k, double v)
   // Adds the pair (k,v) to the receiver. Returns this new pair.
{
   Pair* newPair ;

#  ifdef DEBUG
      if (this->includes(k))  {
	 printf ("Error in Dictionary::add(%c)",k," - key already exists\n");
	 exit(0) ;}
#  endif

   newPair = new Pair(k,v) ;
   if (last)
      last -> append(newPair) ;
   else                                   // empty dictionary
      first = newPair ;
   last = newPair ;

   return newPair ;
}


double&  Dictionary :: at (char aKey)
   // Returns the value of the pair which key is aKey. If such pair does
   // not exist, creates it and assign value 0.
{
   Pair *next,*newPair ;

   next = first ;
   while (next) 
   {
      if (next->giveKey() == aKey)
		return next->giveValue() ;
      next = next->giveNext() ;
   }

   newPair = this->add(aKey,0) ;          // pair does not exist yet
   return newPair->giveValue() ;
}


int  Dictionary :: includes (char aKey)
   // Returns True if the receiver contains a pair which key is aKey, else
   // returns False.
{
   Pair* next ;

   next = first ;
   while (next) 
	{
	  if (next->giveKey() == aKey)
		 return TRUE ;
	  next = next->giveNext() ;
	}

   return FALSE ;
}


void  Dictionary :: printYourself ()
   // Prints the receiver on screen.
{
   Pair* next ;

   printf ("Dictionary : \n") ;

   next = first ;
   while (next) {
      next -> printYourself() ;
      next = next -> giveNext() ;}
}


