//   ************************
//   *** CLASS NODAL LOAD ***
//   ************************
 

#include "load.h"
class TimeStep ;


class NodalLoad : public Load
/*!
   This class implements a concentrated load (force, moment,...) that acts
   directly on a node (not on an element). A nodal load is usually attribute
   of one or more nodes.
 DESCRIPTION
   The attribute 'componentArray' contains, for example for the case of a
   plane beam structure, 2 forces and 1 moment.
 REMARK
   class Node is not restricted to structural problems. For example, in ther-
   mal analysis, a nodal load would be a concentrated heat source.
 */
{
   public :
      NodalLoad (int i,Domain* d) : Load(i,d) {}           // constructor

      FloatArray*  ComputeValueAt (TimeStep*) ;
} ;









