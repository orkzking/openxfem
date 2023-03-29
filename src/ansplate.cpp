

//   file ANSPLATE.CPP 

#include "ansplate.h"

#include "node.h"
#include "material.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "diagmtrx.h"
#include "fei2dquadlin.h"
#include "intarray.h"
#include "domain.h"     
#include "standardquadrature.h" 
#include "splitgaussquadrature.h"
#include "integrationrule.h" 
#include "enrichmentitem.h"
#include "geometry_2D.h"
#include <math.h>
#include <vector>




ANSPlate :: ANSPlate (int n, Domain* aDomain)
: PlateIso(n,aDomain)
{
  //samplingPoints = NULL ;
}

FloatMatrix*  ANSPlate :: ComputeBnsMatrixAt(Mu::Point* p)
// ******************************************************
// B = [N1,xi  -N1dydxi  N1dxdxi ...
//      N1,eta -N1dydeta N1dxdeta ...]
{
  FloatArray *dNdKsi = 
	 (dynamic_cast<FEInterpolation2d*>(this->giveFEInterpolation()))->giveDerivativeKsi(p->x,p->y);
  FloatArray *dNdEta = 
	 (dynamic_cast<FEInterpolation2d*>(this->giveFEInterpolation()))->giveDerivativeEta(p->x,p->y);

  double dXdKsi = 0.0 ; double dYdKsi = 0.0 ;
  double dXdEta = 0.0 ; double dYdEta = 0.0 ; 
  for(size_t i = 0 ; i < numberOfNodes ; i++)
  {
	 Node *aNode = this->giveNode(i+1);
	 double x = aNode->giveCoordinate(1);
	 double y = aNode->giveCoordinate(2);

	 dXdKsi +=  (*dNdKsi)[i] * x ; 
	 dYdKsi +=  (*dNdKsi)[i] * y ; 

	 dXdEta +=  (*dNdEta)[i] * x ; 
	 dYdEta +=  (*dNdEta)[i] * y ; 
  }
  FloatArray* N = this->giveFEInterpolation()->evalN(p) ;
  FloatMatrix *answer = new FloatMatrix(2,3*numberOfNodes);
  for(size_t i = 0 ; i < numberOfNodes ; i++)
  {
	 answer->at(1,3*i+1) = (*dNdKsi)[i] ; 
	 answer->at(1,3*i+2) = -(*N)[i] * dYdKsi ; 
	 answer->at(1,3*i+3) =  (*N)[i] * dXdKsi ; 

	 answer->at(2,3*i+1) = (*dNdEta)[i] ; 
	 answer->at(2,3*i+2) = -(*N)[i] * dYdEta ; 
	 answer->at(2,3*i+3) =  (*N)[i] * dXdEta ; 
  }

  delete dNdKsi ;
  delete dNdEta ;
  delete N ;

  return answer ;
}

FloatMatrix*  ANSPlate :: ComputeBSmatrixAt (GaussPoint *aGaussPoint)
//********************************************************************
{
  FloatMatrix *Bns;
  FloatArray  *row1 = new FloatArray(3*numberOfNodes);
  FloatArray  *row2 = new FloatArray(3*numberOfNodes);
  Mu::Point *Coord  = aGaussPoint -> giveCoordinates() ;   

  // Assumed strain shape functions evaluated at Gauss points !!!
  FloatMatrix *Ns  = ComputeAssumedStrainShapeFunction(aGaussPoint);

  FloatMatrix  *answer = new FloatMatrix(2,3*numberOfNodes);

  for(size_t s = 0 ; s < 0.5*samplingPoints.size() ; s++)
  {
	 // Discretized natural shear strain evaluated at sampling points S !!!
	 Bns = this->ComputeBnsMatrixAt(samplingPoints[s]);
	 row1->add(Bns->giveRow(1)->times(Ns->at(1,s+1) )) ;
	 row2->add(Bns->giveRow(2)->times(Ns->at(2,s+1))) ;
  }

  for(size_t i = 0 ; i < 3*numberOfNodes ; i++)
  {
	 answer->at(1,i+1) = (*row1)[i] ;
	 answer->at(2,i+1) = (*row2)[i] ;
  }

  delete Bns ;
  delete row1 ;
  delete row2 ;
  delete Ns ;

  return answer ;
}
