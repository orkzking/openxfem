//   ******************
//   *** CLASS LIST ***
//   ******************
 

#ifndef _LIST_H_
#define _LIST_H_

class FEMComponent ;


class List
/*!
   This class implements an array which contains elements, nodes, materials,
   loads or load-time functions.
 DESCRIPTION :
   The objects are stored in 'values', an array of FEMComponents of size
   'size'.
 TASKS :
   - storing (method 'put') and returning (method 'at') FEMComponents ;
   - expanding itself, in order to accomodate more FEMComponents.
 */
{
   enum { FALSE } ;

   protected :
      int             size ;
      FEMComponent**  values ;

   public :
      List (int) ;                // constructor
      ~List () ;                  // destructor

      FEMComponent*   at (int i)           {return values[i-1] ;}
      int             giveSize ()          {return size ;}
      void            growTo (int) ;
      int             includes (int) ;
      int             isEmpty ()           {return (size==0) ;}
      int             isNotEmpty ()        {return (size!=0) ;}
      void            printYourself () ;
      void            put (int,FEMComponent*) ;
} ;


#endif  // _LIST_H_
