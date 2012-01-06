//   *****************
//   *** CLASS DOF ***
//   *****************
 

#include "dictionr.h"
#include <stdio.h>
class Domain ; class Node ; class TimeStep ; class BoundaryCondition ; 
class InitialCondition ;


class Dof 
/*!
   This class implements a nodal degree of freedom. A dof is usually attri-
   bute of one node.
 DESCRIPTION
   'number' and 'node' are used for reading/writing data in the data file.
   'equationNumber' is the number of the associated equation in the linear
   system ; it is 0 if the dof is subjected to a boundary condition (b.c.).
   'bc' is the number of the b.c. the dof is subjected to, if any.
   'unknowns' and 'pastUnknowns' are the dictionaries where the dof stores
   its unknowns (e.g., the displacement 'd', the velocity 'v' and the acce-
   leration 'a'), at the current time step and at the previous one.
 TASKS
   - equation numbering, in method 'giveEquationNumber' ;
   - managing its b.c. and its i.c., if any (methods 'hasBc', 'giveBc', etc);
   - managing its unknowns. This includes retrieving the associated solution
     from the linear system, calculating the unknowns (through Newmark, sta-
     tic, etc, formulae), printing the unknowns in the data file ;
 REMARKS
   - class Dof is not a subclass of FEMComponent : a dof belongs to a single
     node, not to the domain ;
   - class Dof is not restricted to structural analysis problems. Unknowns
     may also be pressures, temperatures, etc.
 */
{
   enum { FALSE , TRUE } ;
   public :
      Dof (int,Node*) ;                                       //<! constructor
      ~Dof ()   { delete unknowns ; delete pastUnknowns ;}    //<! destructor.

      double              computeNewmarkUnknown (char,TimeStep*) ;
      double              computeStaticUnknown (char,TimeStep*) ;
      double              computeUnknown (char,TimeStep*) ;
      double              getSolution () ;
      BoundaryCondition*  giveBc () ;
      int                 giveEquationNumber () ;
      InitialCondition*   giveIc () ;
      double              givePastUnknown (char,TimeStep*) ;
      double              giveUnknown (char,TimeStep*) ;
	   double              computeInitialAcceleration (TimeStep* stepN) ; //SC 1.99
      int                 hasBc () ;
      int                 hasIc () ;
      int                 hasIcOn (char) ;
      void                print (char,TimeStep*, FILE*) ;
      void                print (char,char,char,TimeStep*,FILE*) ;
      void                printNewmarkOutputAt (TimeStep* stepN, FILE* disFile)
				      { this -> print('d','v','a',stepN,disFile) ;}
      void                printOutputAt (TimeStep*, FILE*) ;
      void                printStaticOutputAt (TimeStep* stepN, FILE* disFile)
				      { this -> print('d',stepN,disFile) ;}
      void                printYourself () ;
      void                updateYourself () ;

   private :

      int          number ;
      Node*        node ;
      int          equationNumber ;
      int          bc ;
      int          ic ;
      Dictionary*  unknowns ;
      Dictionary*  pastUnknowns ;
} ;

