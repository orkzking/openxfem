//   ******************
//   *** CLASS NODE ***
//   ******************

#ifndef _NODE_H_
#define _NODE_H_

#include "femcmpnn.h"
#include "geometry_base.h"
#include "domain.h"
#include <stdio.h>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
class Dof ; class NodalLoad ; class TimeStep ;
class EnrichmentItem;

using namespace std;

//! Node 
/*!
This class implements a node in a finite element mesh. A node is an attri-
bute of a domain. It is usually also attribute of a few elements.
DESCRIPTION
The node possesses 'numberOfDofs' degrees of freedom, stored in 'dofArray'.
In 'loadArray' it stores the number of every nodal load it is subjected to
(typically, concentrated forces and moments).
In 'locationArray' the node stores the equation number of each of its dofs.
This location array is used by the node for assembling its load vector to
the right-hand side of the linear system ; it is also used by elements for
calculating their own location arrays.
TASKS
- managing its position in space. In geometrically linear analysis, this
only amounts to managing its coordinates ;
- managing its degrees of freedom (method 'giveDof') ;
- calculating its nodal load vector and having it assembled to the linear
system ;
- printing and updating at end of step ;
- managing its swapping to and from disk.
*/
class Node : public FEMComponent
{
private:
  FloatArray*  coordinates ;      //!< coordinates of node
  int          numberOfDofs ;     //!< number of Degree Of Freedom 
  Dof**        dofArray ;         //!< array of Dofs
  IntArray*    loadArray ;        //!< loads acting on the node
  IntArray*    locationArray ;    //!< location array of node in the global matrix 
  int          isEnriched ;       //!< switch of enriched or not    
  list<EnrichmentItem*>*     enrichmentItemListOfNode; //!< list of enrichment items of node  
  bool         isUpdated ; //!< updated or not when crack growth, 2005-09-03.  
  std::map<EnrichmentItem*,double>  levelSet ;//!<
public:

  Node (int,Domain*) ;                        // constructor
  ~Node () ;                                  // destructor

  /** @name Node Coordinates
  *  Methods to get and give coordinates of node
  */
  //@{
  double       giveCoordinate (int) ;
  FloatArray*  giveCoordinates();
  void         getCoordinates () ; //@}


  /** @name Degree of freedom 
  *  Give the Dofs
  */
  //@{
  Dof*         giveDof (int) ;
  int          computeNumberOfDofs () ; // NVP 2005, for XFEM
  size_t       giveNumberOfTrueDofs () ; // XFEM. NVP 2005
  int          giveNumberOfDofs () ; //@}

  /** @name Nodal load vector
  *  Methods for computing location array
  *  compute the load vector and assemble this into the right hand side
  */
  //@{
  void         assembleYourLoadsAt (TimeStep*) ;
  FloatArray*  ComputeLoadVectorAt (TimeStep*) ;
  IntArray*    giveLoadArray () ;
  IntArray*    giveLocationArray () ;
  IntArray*    computeLocationArray () ;
  IntArray*    giveStandardLocationArray () ;
  IntArray*    giveEnrichedLocationArray () ; //@}

  // time step termination
  void         printOutputAt (TimeStep*, FILE*, FILE*) ;
  void         updateYourself () ;

  /** @name Miscellaneous
  */
  //@{
  char*        giveClassName (char* s)      { return strcpy(s,"Node") ;}
  void         instanciateYourself () ;
  void         printYourself () ;
  void         printBinaryResults (TimeStep*, FILE*);
  void         getBinaryRecord (TimeStep*, float*);   //@}

  /*!
  @return 1 if this is enriched node and 0 otherwise.
  */
  int                    getIsEnriched(){return isEnriched;}
  /*!
  Gives the list of EnrichmentItems of the receiver. Used in computation of B matrix.
  */
  list<EnrichmentItem*>* giveEnrItemListOfNode();  
  /*!
  Delete the list of EnrichmentItem of the receiver so that it can be reset
  when crack grows. 2005-08-30.
  */
  void    eraseListOfEnrichmentItems(){enrichmentItemListOfNode->clear() ;}
  /*!
  Insert a new EnrichmentItem objetc into list of enrichment items of node
  Coding for branched cracks. 2005-11-13
  */
  void    addNewEnrichmentItem(EnrichmentItem*);
  /*!
  Erase an EnrichmentItem from the list of enrichment items of node
  Coding for branched cracks. 2005-11-13
  */
  void    deleteEnrichmentItem(EnrichmentItem*);
  /*!
  Insert enrItem to the list of EnrichmentItems of the receiver.
  */
  void    isEnrichedWith(EnrichmentItem* enrItem); 
  /*!
  A node just can be enriched by either asym. funcs or discontinous func. Therefore
  if the enrichmentItemListOfNode contains both CrackInterior and CrackTip
  then we need remove the CrackInterior from this list, to get the final
  enrichmentItemListOfNode
  */
  void          resolveConflictsInEnrichment();
  /*!
  A node is enriched by both CrackInterior and MaterialInterface (coincident geometry)
  Or, just enrich this node by CrackInterior (modify the file CrackInterior.cpp to
  switch between two options)
  */
  void          resolveConflictsInEnrichment2();
  /*!
  Used to remove nodes from the set of Heaviside enriched nodes
  @See Dolbow et al (1999) for details.
  */
  void          resolveLinearDependency(EnrichmentItem*) ;
  /*!
  Just for checking whether the enrichment detection is correct or not.
  @return true if this is tip-enriched node.
  */
  bool          isTipEnriched();
  /*!
  Just for checking whether the enrichment detection is correct or not.
  @return true if this is step-enriched node.
  */
  bool          isStepEnriched();
  bool          isJunctionEnriched(); // branched crack implementation.
  bool          isInterfaceEnriched();// material interface implementation
  /*!
  From the coordinates of the receiver, make a Mu::Point in order to use
  facilities on computational geometry of Cyrille Dunant.
  */
  Mu::Point*    makePoint();
  /*!
  Operator < used to sort elements in List
  */
  bool  operator < (const Node& node2){
	 return number < node2.number ;
  }
  /*!
  Whenever the status of node is changed ( become enriched when crack grows), then
  set variable isUpdated = true .
  */
  void    setIsUpdated(){isUpdated = true ;} 

  void    setLevelSets(EnrichmentItem*,double);
  /*!
  Returns the level set at node which is the signed distance to the 
  \c EnrichmentItem enrItem
  */
  double  giveLevelSet(EnrichmentItem* enrItem);

  // ===============================================================
  //                            DEBUGGING ...
  /*!
  print enriched nodes in file to check
  */
  void         printEnrichedNode();  
  /*!
  Print the support of the receiver, usefull for debugging 
  */
  void         printNodalSupport();
  int          giveMyNumber(){return FEMComponent::giveNumber();}
} ;



#endif //_NODE_H_
