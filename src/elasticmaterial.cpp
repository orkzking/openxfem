//   file ELASTICMATERIAL.CPP

#include "elasticmaterial.h"
#include "dictionr.h"
#include "domain.h"
#include "gausspnt.h"
#include "element.h"
#include "materialinterface.h"
#include "enrichmentitem.h"
#include "feinterpol.h"
#include "geometryentity.h"
#include "functors.h"
#include <stdlib.h>
#include <typeinfo>

#include "planelast.h" 

char* ElasticMaterial :: giveClassName (char* s)
{	return strcpy(s,"ElasticMaterial") ;
}

void  ElasticMaterial :: instanciateYourself ()
{
  double value ;

#  ifdef VERBOSE
  printf ("\ninstanciating material %d\n",number) ;
#  endif

  if (this->readWhetherHas("E")) {
	 value = this -> read("E") ;
	 propertyDictionary -> add('E',value) ;}

  if (this->readWhetherHas("n")) {
	 value = this -> read("n") ;
	 propertyDictionary -> add('n',value) ;}

  if (this->readWhetherHas("d")) {
	 value = this -> read("d") ;
	 propertyDictionary -> add('d',value) ;}

  if (this->readWhetherHas("A")) {
	 value = this -> read("A") ;
	 propertyDictionary -> add('A',value) ;}

  if (this->readWhetherHas("t")) {
	 value = this -> read("t") ;
	 propertyDictionary -> add('t',value) ;}
}


void  ElasticMaterial :: printYourself ()
// Prints the receiver on screen.
{
  printf ("ElasticMaterial with properties : \n") ;
  propertyDictionary -> printYourself() ;
}

void  ElasticMaterial :: ComputeStress(FloatArray* dxacc, Element* elem, GaussPoint* gp)
// MODIFIFED FOR INTERFACE IMPLEMENTATION
{  
  // strains
  FloatArray * deltaEpsilon = elem->computeStrainIncrement(gp,dxacc);

  // delta sigma
  FloatArray *sigma; 
  if(elem->containsMultiMats() == false)
	 sigma = elem->giveConstitutiveMatrix()->Times(deltaEpsilon); // SC-Purify-10.09.97
  else
  {
	 double sd ; FloatMatrix *d ;
	 PlaneElasticity  *planeElast = new PlanStrain; // LIMITATION !!!

	 std::list<EnrichmentItem*> *enrichmentItemList = elem->giveListOfEnrichmentItems();
	 Mu::Point *globalCoord = 
		elem->giveFEInterpolation()->local2Global(domain,elem->giveNodeArray(),gp->giveCoordinates());
	 std::vector<size_t> matIDs = elem->giveMatIDs();

	 std::list<EnrichmentItem*>::iterator it ;
	 for(it = enrichmentItemList->begin() ; it != enrichmentItemList->end() ; it++)
	 {
		if(dynamic_cast<MaterialInterface*>(*it))
		{
		  sd = (*it)->giveMyGeo()->computeSignedDistanceOfPoint(globalCoord);
		  if(sd>0)
		  {			 
			 Material *mat1 = domain->giveMaterial(matIDs[0]); 
			 d = planeElast->computeConstitutiveMatrix(mat1);	
			 sigma = d->Times(deltaEpsilon); 
		  }
		  else
		  {			 
			 Material *mat2 = domain->giveMaterial(matIDs[1]); 
			 d = planeElast->computeConstitutiveMatrix(mat2);		
			 sigma = d->Times(deltaEpsilon); 
		  }
		}
	 }
	 delete globalCoord ;
  }

  // sigma(prev)+delta sigma
  sigma->add(gp->givePreviousStressVector());

  delete deltaEpsilon;

  gp->letStressVectorBe(sigma); //gp stores sigma, not a copy!
}
/*
FloatMatrix*  ElasticMaterial :: ComputeConstitutiveMatrix(GaussPoint* gp, Element* e)
// ************************************************************************************
// Returns a copy of the elasticity matrix of the receiver.
{
return e->giveConstitutiveMatrix()->GiveCopy();
}*/


FloatMatrix*  ElasticMaterial :: ComputeConstitutiveMatrix(GaussPoint* gp, Element* e)
// ************************************************************************************
// Returns a copy of the elasticity matrix of the receiver.
// Modified for multi-materials element.
// X-FEM implementation for material interface.
{
  if(e->containsMultiMats() == false)
	 return e->giveConstitutiveMatrix()->GiveCopy();

  double sd ; FloatMatrix *d ;
  PlaneElasticity  *planeElast = new PlanStrain; // LIMITATION !!!

  std::list<EnrichmentItem*> *enrichmentItemList = e->giveListOfEnrichmentItems();
  Mu::Point *globalCoord = 
	 e->giveFEInterpolation()->local2Global(domain,e->giveNodeArray(),gp->giveCoordinates());
  std::vector<size_t> matIDs = e->giveMatIDs();

  std::list<EnrichmentItem*>::iterator it ;
  for(it = enrichmentItemList->begin() ; it != enrichmentItemList->end() ; it++)
  {
	 if(dynamic_cast<MaterialInterface*>(*it))
	 {
		sd = (*it)->giveMyGeo()->computeSignedDistanceOfPoint(globalCoord);
		//std::cout << sd << std::endl;
		if(sd>0)
		{		  
		  Material *mat1 = domain->giveMaterial(matIDs[0]); 
		  d = planeElast->computeConstitutiveMatrix(mat1);
		  //std::cout<<e->giveNumber()<< " E = " << mat1->give('E') << " ";
		  //globalCoord->print();std::cout<<endl;
		}
		else
		{
		  Material *mat2 = domain->giveMaterial(matIDs[1]); 
		  d = planeElast->computeConstitutiveMatrix(mat2);
		  //std::cout<<e->giveNumber()<< " E = " << mat2->give('E') << " " ;
		  //globalCoord->print();std::cout<<endl;
		}
	 }
  }
  delete globalCoord ;
  return d ;
}
