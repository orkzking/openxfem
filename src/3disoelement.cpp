//   file.CPP
 
#include "3disoelement.h"

#include "node.h"
#include "material.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "domain.h"     
#include "enrichmentitem.h"
#include "geometry_2D.h"
#include "standardquadrature.h" 
#include "splitgaussquadrature.h"
#include "integrationrule.h" 
//#include "geometryentity.h" 
//#include <math.h>

 
Iso3dElement :: Iso3dElement (int n, Domain* aDomain)
      : Element (n,aDomain)
{  
	this->computeConstitutiveMatrix();
}

 
FloatMatrix*  Iso3dElement :: ComputeBuMatrixAt (GaussPoint *aGaussPoint)
//***********************************************************************
// Returns the [6x8] strain-displacement B matrix of the receiver,
// evaluated at aGaussPoint.
// Format of B = [N1,x 0   0  ...
//                 0   N1,y 0 ...
//                 0   0  N1,z... 
//                N1,y N1,x 0 ... 
//                0 N1,z N1,y ... 
//                N1,z 0 N1,x ...] 
{
	Mu::Point *Coord = aGaussPoint -> giveCoordinates() ;
   FloatMatrix *dNdx = this->giveFEInterpolation()->evaldNdx(domain,this->giveNodeArray(),Coord);

	FloatMatrix  *answer = new FloatMatrix(6,3*numberOfNodes);

	for (size_t i = 0 ; i < numberOfNodes ; i++)
	{
		answer->at(1,3*i+1) = dNdx->at(i+1,1); // epsilon_xx
		answer->at(2,3*i+2) = dNdx->at(i+1,2); // epsilon_yy
		answer->at(3,3*i+3) = dNdx->at(i+1,3); // epsilon_zz
		answer->at(4,3*i+1) = dNdx->at(i+1,2); // epsilon_xy
		answer->at(4,3*i+2) = dNdx->at(i+1,1);
		answer->at(5,3*i+2) = dNdx->at(i+1,2); // epsilon_yz
		answer->at(5,3*i+3) = dNdx->at(i+1,3);
		answer->at(6,3*i+1) = dNdx->at(i+1,3); // epsilon_zx
		answer->at(6,3*i+3) = dNdx->at(i+1,1);
	}

   delete dNdx ;

   return answer ;
}

FloatMatrix*  Iso3dElement ::computeConstitutiveMatrix()
// *****************************************************
// computes the constituive matrix of isotropic, linear elastic material
{
   constitutiveMatrix = new FloatMatrix(6,6); 

	Material* mat   = this -> giveMaterial() ;

   double e     = mat -> give('E') ;
   double nu    = mat -> give('n') ;
   double ee    = e / ((1.+nu) * (1.-nu-nu)) ;
   
	constitutiveMatrix->at(1,1) = ee * (1.0 - nu) ;
	constitutiveMatrix->at(1,2) = ee * nu ;
	constitutiveMatrix->at(1,2) = ee * nu ;

	constitutiveMatrix->at(2,1) = ee * nu ;
	constitutiveMatrix->at(2,2) = ee * (1.0 - nu) ;
	constitutiveMatrix->at(2,3) = ee * nu ;

	constitutiveMatrix->at(3,1) = ee * nu ;
	constitutiveMatrix->at(3,2) = ee * nu ;
	constitutiveMatrix->at(3,3) = ee * (1.0 - nu) ;

	constitutiveMatrix->at(4,4) = ee * (1.0 - nu - nu)*0.5 ;
	constitutiveMatrix->at(5,5) = ee * (1.0 - nu - nu)*0.5 ;
	constitutiveMatrix->at(6,6) = ee * (1.0 - nu - nu)*0.5 ;

   return constitutiveMatrix ;
}  

 
FloatMatrix*  Iso3dElement :: ComputeNmatrixAt (GaussPoint* aGaussPoint)
// *********************************************************************
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
// N = [N1 0 0  ...
//      0 N1 0  ...
//      0 0  N1 ...]
{
	Mu::Point *Coord = aGaussPoint -> giveCoordinates() ;
   FloatArray *n = this->giveFEInterpolation()->evalN(Coord);

   FloatMatrix *answer = new FloatMatrix(3,3*numberOfNodes);

	for (size_t i = 0 ; i < numberOfNodes ; i++)
	{
		answer->at(1,3*i+1) = (*n)[i];
		answer->at(2,3*i+2) = (*n)[i];
		answer->at(3,3*i+3) = (*n)[i];
	}

   delete n;

   return answer ;
} 

 
double  Iso3dElement :: computeVolumeAround (GaussPoint* aGaussPoint)
// ******************************************************************
// Returns the portion of the receiver which is attached to aGaussPoint.
{
   FloatMatrix* jacob ;
   double       determinant,weight,volume ;

	Mu::Point *coord  = aGaussPoint -> giveCoordinates() ;
	jacob       = this -> giveFEInterpolation()->giveJacobianMatrixAt(domain,this->giveNodeArray(),coord);
   determinant = fabs (jacob->giveDeterminant()) ;
   weight      = aGaussPoint -> giveWeight() ;
   
   volume      = determinant * weight ;

   delete jacob ;
   return volume ;
}  

void Iso3dElement::setEnrichmentForMyNodes(EnrichmentItem *enrItem)
// ****************************************************************
// High order XFEM : 
// If use linear shape functions for enriched part of the displacement approximation
// then just enrich three corner nodes. Otherwise, if also use quadratic shape funcs.
// enrich all of six  nodes.
{
  size_t intOrder = this->giveXFEInterpolation()->giveInterpolationOrder();

  if(intOrder == 1)          // just enrich three corner nodes
  {
	 for(size_t i = 0 ; i < 3 ; i++)
	 {         
		this->giveNode(i+1)->isEnrichedWith(enrItem)  ;
	 }
  }
  else if(intOrder == 2)     // enrich all nodes
  {
	 for(size_t i = 0 ; i < numberOfNodes ; i++)
	 {         
		this->giveNode(i+1)->isEnrichedWith(enrItem)  ;
	 }
  }
  else
	 assert(false) ;     // can not have higher interpolation for quadratic elements !!!
}

 
