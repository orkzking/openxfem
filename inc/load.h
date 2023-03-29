//   ******************
//   *** CLASS LOAD ***
//   ******************
 

#ifndef _LOAD_H_
#define _LOAD_H_

#include "femcmpnn.h"
#include "domain.h"
#include "flotarry.h"
#include "dictionr.h"


//! Loads : Boundary conditions, Nodal loads, Dead loads
/*!
   This abstract class is the superclass of the classes that implement loads
   (body load, nodal load, boundary conditions, etc). A load is an attribute
   of the domain. It is usually also attribute of several elements, nodes or
   dofs.
 DESCRIPTION
   The load stores its values in 'componentArray'. The components of a load
   at a given time step is the product of 'componentArray' by the value of
   the function 'loadTimeFunction' at that time step.
 TASK
   Returning its components and its load-time function ;
 */
class Load : public FEMComponent
{
   protected :
      FloatArray*    componentArray ;
      int            loadTimeFunction ;

   public :
      Load (int,Domain*) ;                              //!< constructor
      virtual ~Load ()  { delete componentArray ;}      //!< destructor

      // computations
      LoadTimeFunction*  giveLoadTimeFunction () ;
      FloatArray*        giveComponentArray () ;

      // definition of a load
      Load*              typed () ;
      Load*              ofType (char*) ;
      char*              giveClassName (char* s)   {return strcpy(s,"Load");}
      virtual void       instanciateYourself () ;
} ;


#endif // _LOAD_H_

