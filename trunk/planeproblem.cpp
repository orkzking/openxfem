//   file PlaneProblem.CPP
 
#include "planeproblem.h"

#include "node.h"
#include "material.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "domain.h"     
#include "standardquadrature.h" 
#include "splitgaussquadrature.h"
#include "integrationrule.h" 
#include "enrichmentitem.h"
#include "geometry_2D.h"
#include "geometryentity.h" 
#include <math.h>



 
PlaneProblem :: PlaneProblem (int n, Domain* aDomain)
      : Element (n,aDomain)
{  
	planeElast = new PlanStrain;  // better to have a if statement !!!
}

 
FloatMatrix*  PlaneProblem :: ComputeBuMatrixAt (GaussPoint *aGaussPoint)
//***********************************************************************
// Returns the [4 x 2*n] strain-displacement B matrix of the receiver,
// evaluated at aGaussPoint.
// Format of B = [N1,x 0    N2,x 0 ... Nn,x 0
//                 0   N1,y 0  N2,y... 0   Nn,y
//                N1,y N1,x N2,y N2,x ... 
//                0    0     0   0     ] 
{
	Mu::Point *Coord = aGaussPoint -> giveCoordinates() ;
   FloatMatrix *dNdx = this->giveFEInterpolation()->evaldNdx(domain,this->giveNodeArray(),Coord);

	FloatMatrix  *answer = new FloatMatrix(4,2*numberOfNodes);

	for (size_t i = 0 ; i < numberOfNodes ; i++)
	{
		answer->at(1,2*i+1) = dNdx->at(i+1,1);
		answer->at(2,2*i+2) = dNdx->at(i+1,2);
		answer->at(3,2*i+1) = dNdx->at(i+1,2);
		answer->at(3,2*i+2) = dNdx->at(i+1,1);
	}

   delete dNdx ;

   return answer ;
}

FloatMatrix*  PlaneProblem ::computeConstitutiveMatrix()
// *****************************************************
// computes the constituive matrix 
{
   constitutiveMatrix = planeElast->computeConstitutiveMatrix(this->giveMaterial()); 

   return constitutiveMatrix ;
}  

 
FloatMatrix*  PlaneProblem :: ComputeNmatrixAt (GaussPoint* aGaussPoint)
// *********************************************************************
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
// N = [N1 0 N2 0 N3 0 N4 0
//      0 N1 0 N2 0 N3 0  N4]
{
	Mu::Point *Coord = aGaussPoint -> giveCoordinates() ;
   FloatArray *n = this->giveFEInterpolation()->evalN(Coord);

   FloatMatrix *answer = new FloatMatrix(2,2*numberOfNodes);

	for (size_t i = 0 ; i < numberOfNodes ; i++)
	{
		answer->at(1,2*i+1) = (*n)[i];
		answer->at(2,2*i+2) = (*n)[i];
	}

   delete n;

   return answer ;
} 

 
double  PlaneProblem :: computeVolumeAround (GaussPoint* aGaussPoint)
// ******************************************************************
// Returns the portion of the receiver which is attached to aGaussPoint.
{
   FloatMatrix* jacob ;
   double       determinant,weight,thickness,volume ;

	Mu::Point *coord  = aGaussPoint -> giveCoordinates() ;
	jacob       = this -> giveFEInterpolation()->giveJacobianMatrixAt(domain,this->giveNodeArray(),coord);
   determinant = fabs (jacob->giveDeterminant()) ;
   weight      = aGaussPoint -> giveWeight() ;
   thickness   = this -> giveMaterial() -> give('t') ;

   volume      = determinant * weight * thickness ;

   delete jacob ;
   return volume ;
}  

void PlaneProblem::setEnrichmentForMyNodes(EnrichmentItem *enrItem)
// ****************************************************************
// High order XFEM : 
// If use linear shape functions for enriched part of the displacement approximation
// then just enrich three corner nodes. Otherwise, if also use quadratic shape funcs.
// enrich all of six  nodes.
{
  size_t intOrder = this->giveXFEInterpolation()->giveInterpolationOrder();

  if(intOrder == 1)     // linear PUM shape functions enrich all nodes
  {
	 for(size_t i = 0 ; i < numberOfNodes ; i++)
	 {         
		this->giveNode(i+1)->isEnrichedWith(enrItem)  ;
	 }
  } 
  else                  // just enrich three corner nodes
  {
	 for(size_t i = 0 ; i < 0.5*numberOfNodes ; i++)
	 {         
		this->giveNode(i+1)->isEnrichedWith(enrItem)  ;
	 }
  }
}
 