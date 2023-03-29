//   file NODE.CPP

#include "node.h"

#include "dof.h"
#include "nodload.h"
#include "timestep.h"
#include "linsyst.h"
#include "flotarry.h"
#include "intarray.h"
#include "enrichmentitem.h"
#include "element.h"
#include "tri_u.h"
#include "geometryentity.h"
#include "geometry_2D.h"
#include "functors.h"
#include "crackinterior.h"
#include "cracktip.h"
#include "crackjunction.h"
#include "materialinterface.h"
#include <stdlib.h>
#include <list>
#include <map>
#include <vector>
#include <algorithm>

using namespace std;

Node :: Node (int n, Domain* aDomain)
: FEMComponent (n,aDomain)
// Constructor. Creates a node with number n, belonging to aDomain.
{
  numberOfDofs  = 0 ;
  coordinates   = NULL ;
  dofArray      = NULL ;
  loadArray     = NULL ;
  locationArray = NULL ;
  isEnriched    = 0    ;
  enrichmentItemListOfNode = NULL ; 
  isUpdated  = false ; 
}

Node :: ~Node()
// Destructor.
{
  delete coordinates ;
  size_t i = numberOfDofs ;
  if (dofArray) 
  {
	 while (i--)
		delete dofArray[i] ;
	 delete  dofArray ;
  }
  delete loadArray ;
  delete locationArray ;
  delete enrichmentItemListOfNode ;
}


void  Node :: assembleYourLoadsAt (TimeStep* stepN)
// Forms at stepN the vector of the concentrated loads of the receiver ; 
// then, if it exists, assembles it to the system's right-hand side.
{
  LinearSystem *system ;
  FloatArray   *loadVector,*rhs ;
  IntArray     *loc ;

#  ifdef VERBOSE
  printf ("assembling loads of node %d\n",number) ;
#  endif

  loadVector = this -> ComputeLoadVectorAt(stepN) ;

  if (loadVector) 
  {
	 system = (LinearSystem*) (domain->giveNLSolver()->giveLinearSystem()) ;
	 rhs    = system -> giveRhs() ;
	 loc    = this -> giveLocationArray() ;
	 rhs -> assemble(loadVector,loc) ;
	 delete loadVector ;
  }
}


FloatArray*  Node :: ComputeLoadVectorAt (TimeStep* stepN)
// Computes the vector of the nodal loads of the receiver.
{
  int        i,n,nLoads ;
  NodalLoad  *loadN ;
  FloatArray *answer,*contribution ;

  if (this -> giveLoadArray() -> isEmpty())
	 return NULL ;
  else 
  {
	 answer = new FloatArray(0) ;
	 nLoads = loadArray->giveSize() ;          // the node may be subjected
	 for (i = 1 ; i <= nLoads ; i++)           // to more than one load
	 {            
		n            = loadArray -> at(i) ;
		loadN        = (NodalLoad*) domain->giveLoad(n) ;
		contribution = loadN -> ComputeValueAt(stepN) ;      // can be NULL
		answer -> add(contribution) ;
		delete contribution ;
	 }
	 if (answer->giveSize())
		return answer ;
	 else 
	 {
		delete answer ;
		return NULL ;
	 }
  }
}


void  Node :: getCoordinates ()
// Get from the data file all of the coordinates of the receiver.
{
  int numberOfCoordinates = this->readInteger("coord") ;
  coordinates         = new FloatArray(numberOfCoordinates) ;
  for (size_t i = 1 ; i <= numberOfCoordinates ; i++)
	 coordinates->at(i) = this->read("coord",i+1) ;
}


double  Node :: giveCoordinate (int i)
// Returns the i-th coordinate of the receiver.
{
  if (! coordinates)
	 this -> getCoordinates() ;
  return coordinates->at(i) ;
}

FloatArray*  Node :: giveCoordinates()
{
  if (! coordinates)
	 this->getCoordinates();

  return coordinates ;
}

Dof*  Node :: giveDof (int i)
// Returns the i-th degree of freedom of the receiver. Creates the array
// containing the dofs of the receiver, if such array does not exist yet.
{

  if (! dofArray) 
  {
	 dofArray = new Dof* [this->giveNumberOfDofs()] ;
	 for (size_t j = 0 ; j < numberOfDofs ; j++)
		dofArray[j] = new Dof(j+1,this) ;
  }

  return dofArray[i-1] ;
}


IntArray*  Node :: giveLoadArray ()
// Returns the list containing the number of every nodal loads that act on
// the receiver. If this list does not exist yet, constructs it. This list
// is not to be confused with the load vector.
{
  int nLoads ;

  if (! loadArray)    // the list does not exist yet
  {                        
	 nLoads    = this->readIfHas("loads") ;
	 loadArray = new IntArray(nLoads) ;
	 for (size_t i = 1 ; i <= nLoads ; i++)
		loadArray->at(i) = this->readInteger("loads",i+1) ;
  }

  return loadArray ;
}

int  Node :: computeNumberOfDofs ()
// ********************************
// Computes the number of degrees of freedom of the receiver.
// Need to include the additional Dofs of enriched nodes (XFEM implementation)
// NVP 2005
{
   numberOfDofs = this->readInteger("nDofs") ; // real Dofs ( standard ones )

  // if the receiver is a non-enriched node, stop. Otherwise it has 
  // ( NSD * number_of_enrichmentfunctions ) additional Dofs for each enrichment item.

  if ( isEnriched == 0 )  
	 return numberOfDofs ;

  size_t enrichedDofs = 0 ;
  list<EnrichmentItem*>* enrItemList = this->giveEnrItemListOfNode();

  list<EnrichmentItem*>::iterator iter ;
  for(iter = enrItemList->begin(); iter != enrItemList->end(); iter++)
  {
	 vector<EnrichmentFunction*>*	enrFnVector = (*iter)->giveEnrFuncVector();
	 enrichedDofs += 2 * enrFnVector->size();
  }

  numberOfDofs += enrichedDofs ;

  return numberOfDofs ;
}

size_t  Node :: giveNumberOfTrueDofs ()
// *************************************
// True Dofs are displacement Dofs, just read from input file
{
  size_t nDofs = this->readInteger("nDofs") ;
  return nDofs ;
}

int  Node :: giveNumberOfDofs ()
// *****************************
// Returns the number of degrees of freedom of the receiver.
{
  if (numberOfDofs == 0)
	 this->computeNumberOfDofs() ; 
  else if(isUpdated)
  {
	 this->computeNumberOfDofs() ; // recompute number of Dofs !!!
  }
  return numberOfDofs ;
}


void  Node :: instanciateYourself ()
// Gets from the data file all the data of the receiver.
{
  int i ;

#  ifdef VERBOSE
  printf ("instanciating node %d\n",number) ;
#  endif

  this -> getCoordinates() ;
  this -> giveLocationArray () ;
  this -> giveLoadArray() ;

  numberOfDofs = this -> readInteger("nDofs") ;
  for (i=1 ; i<=numberOfDofs ; i++)
	 this -> giveDof(i) -> hasBc() ;
}


void  Node :: printOutputAt (TimeStep* stepN, FILE* disFile, FILE* s00File)
{
 
#  ifdef VERBOSE
  printf ("node %d printing output\n",number) ;
#  endif

  for (size_t i= 1 ; i <= numberOfDofs ; i++) 
	 this -> giveDof(i) -> printOutputAt(stepN, disFile) ;

  this->printBinaryResults(stepN, s00File);

}

void  Node :: printBinaryResults (TimeStep* stepN, FILE* s00File)
// Prints the strains and stresses on the data file.
// Modified by NVP ( for XFEM ), additional Dofs !!!
{	
  float a[3];         // Why 3 ??? Is this a limitation of the code?

  if (s00File != NULL) {
	 //to put it at the end of the file
	 int pos = fseek(s00File, 0L, SEEK_END);
	 //writing
	 this->getBinaryRecord(stepN, a);
	 fwrite(a, sizeof(a), 1, s00File);
  }
}

void  Node :: getBinaryRecord (TimeStep* stepN, float* a)
// ******************************************************
// Modified by NVP ( for XFEM ), additional Dofs !!!
{
  size_t numberOfDofs = this->readInteger("nDofs") ; // just need "actual" displacement field

  for (size_t i = 1 ; i <= numberOfDofs ; i++) { 
	 a[i-1] = (float)(this -> giveDof(i) -> giveUnknown('d',stepN));
  }
  //water-pressure = 0
  a[2] = (float)0;
}

void  Node :: printYourself ()
// Prints the receiver on screen.
{
  double x = this->giveCoordinate(1) ;
  double y = this->giveCoordinate(2) ;

  printf ("Node %d    coord : x %f  y %f\n",number,x,y) ;
  for ( size_t i = 0 ; i < numberOfDofs ; i++) {
	 if (dofArray[i])
		dofArray[i] -> printYourself() ;
	 else
		printf ("dof %d is nil \n",i+1) ;
  }

  if (locationArray)
	 locationArray->printYourself() ;
  else
	 printf ("locationArray = nil \n") ;

  if (loadArray)
	 loadArray->printYourself() ;
  else
	 printf ("loadArray = nil \n") ;

  printf ("\n") ;
}


void  Node :: updateYourself()
// ***************************
// Updates the receiver at end of step.
// Modification made by NVP for XFEM. 2005-09-05
{
#  ifdef VERBOSE
  printf ("updating node %d\n",number) ;
#  endif

  if(this->domain->isXFEMorFEM() == false)      // FEM problems
  {
	 for (size_t i = 0 ; i < numberOfDofs ; i++)
		this -> giveDof(i+1) -> updateYourself() ;
  }
  else                                          // XFEM problems
  {
	 for (size_t i = 0 ; i < numberOfDofs ; i++)
		delete dofArray[i] ;
	 delete  dofArray ;
	 dofArray = NULL ;
    //numberOfDofs = 0 ;

	 delete locationArray ;
	 locationArray = NULL ;

	 delete loadArray ;
	 loadArray = NULL ;
  }
}


IntArray*  Node :: giveLocationArray ()
// ************************************
// Returns the location array of the receiver. 
// The location array contains the equation number of
// every  degree of freedom of the receiver.
{
  if (!locationArray)
	 this -> computeLocationArray() ;
  return locationArray ;
}


IntArray*  Node :: computeLocationArray()
// **************************************
// compute the location array of the receiver.
// Having form : u = [u1 u2 a1 a2 ] where u1,u2 are displacement DOFs
// and a1,a2 are enriched DOFs( assuming that this node is enriched by 1 function)
{
  locationArray = new IntArray(this->giveNumberOfDofs()) ;

  for (size_t i = 1 ; i <= numberOfDofs ; i++)
	 locationArray->at(i) = this->giveDof(i)->giveEquationNumber() ;
 
  /*// DEBUG ...
  if(isEnriched)
  {
	 std::cout << " location array of node " <<  this->giveNumber()<< endl ;
	 for(size_t i = 0 ; i < locationArray->giveSize() ; i++) 
		std::cout << (*locationArray)[i] << " " ;
	 std::cout << std::endl ; 
  }
  // ----------------------------------------------------------------------*/

  return locationArray ;
}

IntArray*  Node :: giveStandardLocationArray()
// *******************************************
// return the location of standard DOFs, i.e., displacement DOFs
{
  if(isEnriched == 0)                          // standard nodes ( non-enriched)
	 return this->giveLocationArray();

  IntArray *temp = this->giveLocationArray();

  IntArray *ret = new IntArray(2) ;

  ret->at(1) = temp->at(1);
  ret->at(2) = temp->at(2);

  return ret;
}

IntArray*  Node :: giveEnrichedLocationArray ()
// ********************************************
{
  // standard nodes ( non-enriched)
  if(isEnriched == 0)
	 return NULL ;

  // enriched nodes
  IntArray *temp = this->giveLocationArray();
  IntArray *ret = new IntArray(temp->giveSize()-2) ;

  for(size_t i = 1 ; i <= ret->giveSize() ; i++)
  {
	 ret->at(i) = temp->at(i+2);
  }

  return ret;
}


list<EnrichmentItem*>* Node::giveEnrItemListOfNode()
// *************************************************
// before give the list of enrichment items, make sure there
// is no conflicts in this !
{
  //this->resolveConflictsInEnrichment();
  return enrichmentItemListOfNode;
}


void Node :: isEnrichedWith(EnrichmentItem* enrItem)
// *************************************************
// If node is enriched with enrichment item enrItem, then insert 
// enrItem into the list of enrichment items
// also sets member isEnriched equals to 1.
{
  if(enrichmentItemListOfNode == NULL)
	 enrichmentItemListOfNode = new std::list<EnrichmentItem*>;

  if( find(enrichmentItemListOfNode->begin(),enrichmentItemListOfNode->end(),enrItem)
	 == enrichmentItemListOfNode->end())
	 enrichmentItemListOfNode->push_back(enrItem);

  isEnriched = 1 ; // this node is now enriched 

  this->domain->setEnrichedNodesList(this); 
}

void Node :: addNewEnrichmentItem(EnrichmentItem *enrItem)
{
  enrichmentItemListOfNode->push_back(enrItem);
}

void Node :: deleteEnrichmentItem(EnrichmentItem *enrItem)
{
  enrichmentItemListOfNode->remove(enrItem);
}

void Node :: resolveConflictsInEnrichment()
// ****************************************
// check if this  list contains both a CrackTip and a CrackInterior
// then remove the CrackInterior from the list.
// functor IsType<CrackTip,EnrichmentItem>() defined generically in file "functors.h"
{
  list<EnrichmentItem*> ::iterator iter1,iter2;

  iter1=find_if(enrichmentItemListOfNode->begin(),enrichmentItemListOfNode->end(),IsType<CrackTip,EnrichmentItem>());
  iter2=find_if(enrichmentItemListOfNode->begin(),enrichmentItemListOfNode->end(),IsType<CrackInterior,EnrichmentItem>());
  
  if (iter1 != enrichmentItemListOfNode->end() && iter2 != enrichmentItemListOfNode->end()) 
	 enrichmentItemListOfNode->remove(*iter2);
  /*
  // debug: check the list after being removed
  // ----------------------------------------------------------------------------------
  list<EnrichmentItem*> ::iterator it;
  for(it = enrichmentItemListOfNode->begin();it!=enrichmentItemListOfNode->end();it++)
  {
	 (*it)->printYourSelf();
  }
  // ----------------------------------------------------------------------------------
  */
}

void Node :: resolveConflictsInEnrichment2()
// *****************************************
{
  list<EnrichmentItem*> ::iterator it1,it2;

  it1 = find_if(enrichmentItemListOfNode->begin(),enrichmentItemListOfNode->end(),IsType<CrackInterior,EnrichmentItem>());
  it2 = find_if(enrichmentItemListOfNode->begin(),enrichmentItemListOfNode->end(),IsType<MaterialInterface,EnrichmentItem>());

  if (it1 != enrichmentItemListOfNode->end() && it2 != enrichmentItemListOfNode->end()) 
	 enrichmentItemListOfNode->remove(*it2);
}

void Node :: resolveLinearDependency(EnrichmentItem *enrItem)
// **********************************************************
// Using the area inclusion criterion to remove Heaviside enriched nodes
// NVP 2005-06-04
{
  // This tolerance is not constant, but depends on the mesh size h !!!
  // how to implement this?
  double eps = 0.0001 ;
  
  // get the support of this Node from the Domain
  map<Node*,vector<Element*> > nodeElemMap = this->domain->giveNodalSupports();
  vector<Element*>  mySupport = nodeElemMap[this] ; 

  double A   = 0.0 ;   // A is the area of the support
  double A1  = 0.0 ;   // the area above the CrackInterior
  double temp = 0.0 ;

  vector<Element*>* elemList = enrItem->giveElementsInteractWithMe(); 
  vector<Element*> :: iterator iter ;
  // Loop on elements in nodal support, say e. If e intersects with enrItem, then
  // compute the area of e above enrItem by e->giveAreaAboveEnrItem() ;
  for(size_t i = 0 ; i < mySupport.size() ; i++)
  {
	 A  += mySupport[i]->area();
	 iter = find(elemList->begin(),elemList->end(),mySupport[i]); 
	 if (iter != elemList->end()) 
	 {
		A1 += mySupport[i]->computeAreaAboveEnrItem(enrItem); // area of part of element above the enrItem
		temp += mySupport[i]->area() ; // area of elements in nodal support which are split by the crack
	 }
  }
  GeometryEntity *geoOfEnrItem = enrItem->giveMyGeo() ;
  Mu::Point *p = this->makePoint() ;
  int check = geoOfEnrItem->givePositionComparedTo(p);

  if(check == 1)       // node is above the crack 
	 A1 += A - temp ;
  double A2 = A - A1 ; // the area below the Crack

  // finally the area inclusion criteria is applied
  if ((A1/A < eps) || (A2/A < eps))
  {
	 enrichmentItemListOfNode->remove(enrItem); //do not enrich this node with H(x) any more 
	 // after being removed enrichment, if enrichmentItemListOfNode does not contain any
	 // enrichment item, then this node is a classical node !!! NVP 2005-07-12 
	 if(enrichmentItemListOfNode->size() == 0 )
	 {
		isEnriched = 0 ;
		// update the enriched-nodes-list of DOMAIN !!!
		this->domain->removeNodeFromEnrichedNodesList(this); 
	 }
  }

  delete p ;
  /*
  // ------------------------- DEBUG --------------------------------------------
  std::cout<< " checking node " << number << " " << std::endl ;
  std::cout<< " The first  ratio  :     " << A1/A    << std::endl;
  std::cout<< " The second ratio  :     " << A2/A   << std::endl;
  // ----------------------------------------------------------------------------
  */
}

bool Node :: isTipEnriched()
//**************************
// return true if the receiver enriched by the asymptotic functions
// just for plot :)
{
  if( isEnriched == 0 )
    return false ;

  list<EnrichmentItem*> ::iterator it;

  it = 
	 find_if(enrichmentItemListOfNode->begin(),enrichmentItemListOfNode->end(),IsType<CrackTip,EnrichmentItem>());
  
  return (it != enrichmentItemListOfNode->end()) ? true : false ;
}

bool Node :: isStepEnriched()
//***************************
// return true if the receiver enriched by the step function
// just for plot :)
{
  if( isEnriched == 0 )
    return false ;

  list<EnrichmentItem*> ::iterator it;

  it = 
	 find_if(enrichmentItemListOfNode->begin(),enrichmentItemListOfNode->end(),IsType<CrackInterior,EnrichmentItem>());
  
  return (it != enrichmentItemListOfNode->end()) ? true : false ;
}

bool Node :: isJunctionEnriched()
//*******************************
// return true if the receiver enriched by the junction function
// just for plot :)
{
  if( isEnriched == 0 )
    return false ;

  list<EnrichmentItem*> ::iterator it;

  it = 
	 find_if(enrichmentItemListOfNode->begin(),enrichmentItemListOfNode->end(),IsType<CrackJunction,EnrichmentItem>());
  
  return (it != enrichmentItemListOfNode->end()) ? true : false ;
}

bool Node :: isInterfaceEnriched()
//********************************
// return true if the receiver enriched by the material interface
// just for plot :)
{
  if( isEnriched == 0 )
    return false ;

  list<EnrichmentItem*> ::iterator it;

  it = 
	 find_if(enrichmentItemListOfNode->begin(),enrichmentItemListOfNode->end(),IsType<MaterialInterface,EnrichmentItem>());
  
  return (it != enrichmentItemListOfNode->end()) ? true : false ;
}

Mu::Point* Node::makePoint()
//**************************
{
  double x1 = this->giveCoordinate(1);
  double y1 = this->giveCoordinate(2);

  Mu::Point *result = new Mu::Point(x1,y1) ;

  return result ;
}

void Node :: printNodalSupport()
//******************************
// Only for debugging
{
  std::map<Node*,std::vector<Element*> > nodeElemMap = this->domain->giveNodalSupports();
  std::vector<Element*> mySupport = nodeElemMap[this] ;
  
  std::cout << "  Support of node " << number << " : " ;
  for(size_t i = 0 ; i < mySupport.size() ; i++)
	 std::cout << mySupport[i]->giveNumber() << " " ;
  std::cout << std::endl;
}

void Node :: printEnrichedNode()
//******************************
// only for debugging
{
  if(isEnriched)
  {
	 std::cout<< " Node " << number << " is enriched by: " << endl ;
	 list<EnrichmentItem*> ::iterator it;
	 for(it = enrichmentItemListOfNode->begin() ; it != enrichmentItemListOfNode->end() ; it++)
		(*it)->printYourSelf();
	 std::cout << std::endl ;
	 std::cout<< " Number of Dofs " << this->giveNumberOfDofs() << std::endl ;
  }
}

void Node::setLevelSets(EnrichmentItem* enrItem, double ls)
{
  pair<EnrichmentItem*,double> a_pair;
  a_pair = make_pair(enrItem,ls);
  levelSet.insert(a_pair);
}

double Node::giveLevelSet(EnrichmentItem *enrItem)
{
  return levelSet[enrItem];
}