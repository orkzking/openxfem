

//   file PLATEISO.CPP
 
#include "plateiso.h"

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
#include "piecewiselinear.h"
#include "vertex.h"
#include "geometry_2D.h"
#include "geometryentity.h" 
#include "functors.h"
#include "crackinterior.h"
#include "cracktip.h"
#include <math.h>
#include <vector>



 
PlateIso :: PlateIso(int n, Domain* aDomain)
      : Element(n,aDomain)
{   
	this->computeGaussPointsS();
}

 
FloatMatrix*  PlateIso :: ComputeBuMatrixAt (GaussPoint *aGaussPoint)
//********************************************************************
// Returns the [4x8] strain-displacement B matrix of the receiver,
// evaluated at aGaussPoint.
// Format of B = [BB;BS]
{	
	FloatMatrix  *answer = new FloatMatrix(5,3*numberOfNodes);
   FloatMatrix  *BB = this->ComputeBBmatrixAt(aGaussPoint);
	FloatMatrix  *BS = this->ComputeBSmatrixAt(aGaussPoint);

	for (size_t i = 0 ; i < 5 ; i++)
	{
		answer->at(1,i+1) = BB->at(1,i+1);		
		answer->at(2,i+1) = BB->at(2,i+1);		
		answer->at(3,i+1) = BB->at(3,i+1);		
		answer->at(4,i+1) = BS->at(1,i+1);		
		answer->at(5,i+1) = BS->at(2,i+1);		
	}
	delete BB ; delete BS ; 
   return answer ;
}

FloatMatrix*  PlateIso :: ComputeBBmatrixAt (GaussPoint *aGaussPoint)
//********************************************************************
// Returns the [3x8] bending strain-displacement BB matrix of the receiver,
// evaluated at aGaussPoint.
// Format of BB = [0 0 N1,x     ...
//                 0 -N1,y 0    ...
//                 0 -N1,x N1,y ] 
{
	Mu::Point *Coord = aGaussPoint -> giveCoordinates() ;
   FloatMatrix *dNdx = this->giveFEInterpolation()->evaldNdx(domain,this->giveNodeArray(),Coord);

	FloatMatrix  *answer = new FloatMatrix(3,3*numberOfNodes);

	for (size_t i = 0 ; i < numberOfNodes ; i++)
	{
		answer->at(1,3*i+3) = dNdx->at(i+1,1);
		answer->at(2,2*i+2) = -dNdx->at(i+1,2);
		answer->at(3,2*i+2) = -dNdx->at(i+1,1);
		answer->at(3,3*i+3) = dNdx->at(i+1,2);
	}

   delete dNdx ;

   return answer ;
}

FloatMatrix*  PlateIso :: ComputeBSmatrixAt (GaussPoint *aGaussPoint)
//********************************************************************
// Returns the [2x8] shearing strain-displacement BB matrix of the receiver,
// evaluated at aGaussPoint.
// Format of BS = [N1,x 0  N1 ...
//                 N1,y -N1 0 ...] 
{
	Mu::Point *Coord  = aGaussPoint -> giveCoordinates() ;
   FloatMatrix *dNdx = this->giveFEInterpolation()->evaldNdx(domain,this->giveNodeArray(),Coord);
   FloatArray *n     = this->giveFEInterpolation()->evalN(Coord);

	FloatMatrix  *answer = new FloatMatrix(2,3*numberOfNodes);

	for (size_t i = 0 ; i < numberOfNodes ; i++)
	{
		answer->at(1,2*i+1) = dNdx->at(i+1,1);
		answer->at(1,3*i+3) = n->at(i+1);
		answer->at(2,2*i+1) = dNdx->at(i+1,2);		
		answer->at(2,2*i+2) = -n->at(i+1);
	}

   delete dNdx ;

   return answer ;
}

FloatMatrix*  PlateIso ::computeConstitutiveMatrix()
// *************************************************
// computes the constituive matrix D
// D = [Eb [0]
//      [0] Es]
{
   constitutiveMatrix = new FloatMatrix(5,5); 
	FloatMatrix* Eb = this->computeConstitutiveMatrixB();
	FloatMatrix* Es = this->computeConstitutiveMatrixS();
	for(size_t i = 1 ; i <= 3 ; i++)
	{
		constitutiveMatrix->at(1,i) = Eb->at(1,i);
		constitutiveMatrix->at(2,i) = Eb->at(2,i);
		constitutiveMatrix->at(3,i) = Eb->at(3,i);		
	}

	constitutiveMatrix->at(4,4) = Eb->at(1,1);
	constitutiveMatrix->at(5,5) = Eb->at(2,2);

   delete Eb ; delete Es ;
   return constitutiveMatrix ;
}  

FloatMatrix*  PlateIso ::computeConstitutiveMatrixB()
// ***************************************************
// computes the constituive matrix coresponding to bending part
{
  Material *mat = this -> giveMaterial() ;

  double h  = mat -> give('t') ;
  double e  = mat -> give('E') ;
  double nu = mat -> give('n') ;

  double D = e*h*h*h/(12.*(1.-nu*nu)); // flexural rigidity

  FloatMatrix *ret = new FloatMatrix(3,3);
  ret->at(1,1) = D ; ret->at(1,2) = D*nu ;
  ret->at(2,1) = D*nu  ; ret->at(2,2) = D ;
  ret->at(3,3) = 0.5*D*(1. - nu) ;

  return ret;
}  

FloatMatrix*  PlateIso ::computeConstitutiveMatrixS()
// ***************************************************
// computes the constituive matrix coresponding to shearing part
{
  Material *mat = this -> giveMaterial() ;

  double h  = mat -> give('t') ;
  double e  = mat -> give('E') ;
  double nu = mat -> give('n') ;
  
  double k = 5./6. ;                       // shear factor
  double fac = k*e*h/(2.*(1.0 + nu)); 

  FloatMatrix *ret = new FloatMatrix(2,2);
  ret->at(1,1) = fac ; 
  ret->at(2,2) = fac ;

  return ret;
}  

FloatMatrix*  PlateIso :: computeStiffnessMatrix()
// ************************************************
// Computes numerically the stiffness matrix of the receiver.
// [Ke] = [Kb] + [Ks]
{
  this->computeBendingStiffnessMatrix();
  FloatMatrix *stiffnessMatrixS = this->computeShearingStiffnessMatrix();
  stiffnessMatrix = stiffnessMatrix ->GiveCopy() ;
  stiffnessMatrix->plus(stiffnessMatrixS);

  delete stiffnessMatrixS ;
  return stiffnessMatrix ;
}

FloatMatrix* PlateIso :: computeBendingStiffnessMatrix()
//*******************************************************
// Computes numerically the bending stiffness matrix of the receiver, [Kb]
{
   FloatMatrix *b,*db,*d ;
   GaussPoint  *gp ;
	double dV ;
   stiffnessMatrix = new FloatMatrix() ;
   for(int i = 0 ; i < numberOfGaussPoints ; i++)
	{
      gp = gaussPointArray[i] ;
      b  = this -> ComputeBBmatrixAt(gp) ;
      d  = this -> computeConstitutiveMatrixB() ;
      dV = this -> computeVolumeAround(gp) ;
      db = d -> Times(b) ;
      stiffnessMatrix -> plusProduct(b,db,dV) ;
	  
      delete b ;
      delete db ;
	}

   return stiffnessMatrix -> symmetrized() ;
}

FloatMatrix* PlateIso :: computeShearingStiffnessMatrix()
//********************************************************
// Computes numerically the shearing stiffness matrix of the receiver, [Ks]
{
   FloatMatrix *b,*db,*d ;
   GaussPoint  *gp ;
	double dV ;
   FloatMatrix *ret = new FloatMatrix() ; 
   for(int i = 0 ; i < 1 ; i++)
	{
      gp = gaussPointArrayS[i] ;
      b  = this -> ComputeBSmatrixAt(gp) ;
      d  = this -> computeConstitutiveMatrixS() ;
      dV = this -> computeVolumeAround(gp) ;
      db = d -> Times(b) ;
      ret -> plusProduct(b,db,dV) ;
	  
      delete b ;
      delete db ;
	}

   return ret -> symmetrized() ;
}

FloatMatrix*  PlateIso :: ComputeNmatrixAt (GaussPoint* aGaussPoint)
// ***************************************************************
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
// N = [N1 0 N2 0 N3 0 N4 0
//      0 N1 0 N2 0 N3 0  N4]
{
	Mu::Point *Coord = aGaussPoint -> giveCoordinates() ;
   FloatArray *n = this->giveFEInterpolation()->evalN(Coord);

   FloatMatrix *answer = new FloatMatrix(2,8);

	for (size_t i = 0 ; i < 4 ; i++)
	{
		answer->at(1,2*i+1) = (*n)[i];
		answer->at(2,2*i+2) = (*n)[i];
	}

   delete n;

   return answer ;
} 

 
double  PlateIso :: computeVolumeAround (GaussPoint* aGaussPoint)
// ***************************************************************
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

 
void  PlateIso :: computeGaussPoints()
// ************************************
// Sets up the array containing the Gauss points of the receiver.
// Modified to use class Integration Rule and discontinuous Gauss quadrature.
{
  numberOfIntegrationRules = 1 ;
  quadratureRuleArray = new IntegrationRule*;

  int numOfGPs = 16 ;
  if (this->isEnriched() == false )  // classical elements
  { 
	 quadratureRuleArray[0] = new StandardGaussLegendreQuadrature;
	 quadratureRuleArray[0]->setUpIntegrationPoints(SQUARE,4);
  }
  else if (enrichmentItemListOfElem != NULL) // split, tip elements
  {
	 quadratureRuleArray[0] = new SplitGaussLegendreQuadrature;
	 quadratureRuleArray[0]->setUpIntegrationPoints(this,13);
  }
  else  // partially enriched elements 
  {
	 size_t count = 0 ;
	 for(size_t i = 0 ; i < numberOfNodes ; i++)
	 {
		if(this->giveNode(i+1)->isTipEnriched())
		  count += 1 ;
	 }
	 if(count != 0) // near tip enriched elements
	 {
		quadratureRuleArray[0] = new StandardGaussLegendreQuadrature;
		quadratureRuleArray[0]->setUpIntegrationPoints(SQUARE,numOfGPs);
	 }
	 else           // step enriched elements
	 {
		quadratureRuleArray[0] = new StandardGaussLegendreQuadrature;
		quadratureRuleArray[0]->setUpIntegrationPoints(SQUARE,4);
	 }
  }

  // tranform the obtained integration points into the GaussPoint array of element

  vector<double>* weight = quadratureRuleArray[0]->giveWeightArray();
  vector<Mu::Point*>* coordArray = quadratureRuleArray[0]->giveIntegrationPointVector();	

  numberOfGaussPoints = coordArray->size();

  gaussPointArray  = new GaussPoint* [numberOfGaussPoints] ;

  for(size_t i = 0 ; i < numberOfGaussPoints ; i++)
  {
	 gaussPointArray[i] = new GaussPoint(this,i+1,coordArray->at(i),weight->at(i),4) ;
  }
  delete weight ; delete coordArray ;   // Purify, 08-11-05
}

void  PlateIso :: computeGaussPointsS ()
// **************************************
{  
   gaussPointArrayS    = new GaussPoint* [1] ;
	Mu::Point *p = new Mu::Point(0,0);
   gaussPointArrayS[0] = new GaussPoint(this,1,p,4.,4) ;  
 }

void PlateIso::setEnrichmentForMyNodes(EnrichmentItem *enrItem)
// **********************************************************
// loop on its nodes and enrich them with enrItem
{
  for(size_t j = 0 ; j < numberOfNodes ; j++)
  {         
	 this->giveNode(j+1)->isEnrichedWith(enrItem)  ;
  }			 
}
 
std::vector<DelaunayTriangle *>* PlateIso:: PartitionMySelf()
// ********************************************************
// 2005-08-23.
{
  std::vector<Mu::Point>  temp,pts,intersectPoints, tam ;
  GeometryEntity* myGeo;
  
  // TIP-ELEMENTS INTERACTED WITH BOTH THE TIP AND  THE CRACK INTERIOR,
  // THEREFORE, NEED TO REMOVE CRACK INTERIOR FROM THE LIST OF ENRICHMENT ITEMS. 

  list<EnrichmentItem*> ::iterator iter1,iter2;
  iter1=find_if(enrichmentItemListOfElem->begin(),enrichmentItemListOfElem->end(),IsType<CrackTip,EnrichmentItem>());
  iter2=find_if(enrichmentItemListOfElem->begin(),enrichmentItemListOfElem->end(),IsType<CrackInterior,EnrichmentItem>());

  if (iter1 != enrichmentItemListOfElem->end() && iter2 != enrichmentItemListOfElem->end()) 
	 enrichmentItemListOfElem->remove(*iter2);

  // LOOP ON ENRICHMENT ITEMS AND GET THE INTERSECTION POINTS 
  for(std::list<EnrichmentItem*>::iterator i = enrichmentItemListOfElem->begin() ; i != enrichmentItemListOfElem->end() ; i++)
  {
	 if( (typeid(**i) == typeid(CrackTip)) && (dynamic_cast<CrackTip*>(*i)->amIOnElementEdge()) )
	 {
		intersectPoints.push_back(*this->giveMyCenter());
	 }
	 else
	 {
		myGeo = (*i)->giveMyGeo();                // get the geometry of the enrichment item

		for(size_t j = 0 ; j < this->giveNumberOfNodes() ; j++)
		{
		  Mu::Point *p = this->giveNode(j+1)->makePoint() ;
		  pts.push_back(*p);
		  delete p ;
		}

		for(size_t j = 0 ; j < pts.size() ;  j++)
		{
		  Mu::Segment s(pts[j],pts[(j+1)%pts.size()]) ;
		  if(myGeo->intersects(&s))
		  {
			 tam = myGeo->intersection(&s) ;
			 temp.insert(temp.end(),tam.begin(),tam.end()) ;
		  }
		}

		if(temp.size() == 1)     // for the case crack touched the element edge
		{
		  Mu::Point *p1 = this->giveNode(1)->makePoint();
		  Mu::Point *p2 = this->giveNode(2)->makePoint();
		  Mu::Point *p3 = this->giveNode(3)->makePoint();
		  Mu::Point *p4 = this->giveNode(4)->makePoint();
		  Mu::Segment *s1 = new Mu::Segment(*p1,*p2) ;
		  Mu::Segment *s2 = new Mu::Segment(*p2,*p3) ;
		  Mu::Segment *s3 = new Mu::Segment(*p3,*p4) ; 
		  Mu::Segment *s4 = new Mu::Segment(*p4,*p1) ; 

		  PiecewiseLinear *geo = dynamic_cast<PiecewiseLinear*>(myGeo);
		  std::list<Mu::Point*> *vertices = geo->giveMyListOfVertices();

		  for (std::list<Mu::Point*> ::iterator k = vertices->begin() ; k != vertices->end(); k++)
		  {
			 if( s1->on(*k) || s2->on(*k) || s3->on(*k)|| s4->on(*k) )
			 {
				(*k)->print();
				temp.push_back(**k) ;
			 }
		  }
		  delete s1 ; delete s2 ; delete s3 ; delete s4 ;
		  delete p1 ; delete p2 ; delete p3 ; delete p4 ; // Purify, 15-10-2005
		}

		intersectPoints.insert(intersectPoints.end(),temp.begin(),temp.end());

		// Nucleation crack : no intersection between crack & element!!!
		if(intersectPoints.size()==0)
		{
		  CrackInterior* myCrack = 
			 dynamic_cast<CrackTip*>(*i)->giveMyCrackInterior();
		  PiecewiseLinear *crack_geo = dynamic_cast<PiecewiseLinear*>(myCrack->giveMyGeo()); 
		  std::list<Mu::Point*> *vertices  = crack_geo->giveMyListOfVertices();

		  std::list<Mu::Point*>::iterator it;

		  it = vertices->begin();
		  double x1 = (*it)->x ; 
		  double y1 = (*it)->y ;

		  it = vertices->end();
		  double x2 = (*--it)->x ; 
		  double y2 = (*--it)->y ;
		  
		  Mu::Point *p1 = new Mu::Point(x1,y1);
		  Mu::Point *p2 = new Mu::Point(x2,y2);
		  intersectPoints.push_back(*p1);
		  intersectPoints.push_back(*p2);
		}

		// ************************************************************************************
		// ********                    LOOK FOR KINK POINT                               ******
		// ************************************************************************************
		// just segments have kink points, i.e., geometry of CrackInterior only
		// So, use typeid() to test ...
      // Attention: just greater than 2 vertices => kink point :) 2005-09-13

		if( typeid(**i) == typeid(CrackInterior) )
		{
		  std::list<Mu::Point*> *vertices = dynamic_cast<PiecewiseLinear*>(myGeo)->giveMyListOfVertices();
		  if(vertices->size() != 2)
		  {
			 for(std::list<Mu::Point*> ::iterator k = vertices->begin() ; k != vertices->end(); k++)
			 {
				if( this->isWithinMe(*k) )
				  intersectPoints.push_back(**k) ;
			 }
		  }
		}          // end of looking for kink points  */
	 }            // end of check CrackTip
  }              // end of loop on enrichment items

  // make the four nodes Mu::Point
  Mu::Point *p0 = new Mu::Point(this->giveNode(1)->giveCoordinate(1),this->giveNode(1)->giveCoordinate(2)) ;
  Mu::Point *p1 = new Mu::Point(this->giveNode(2)->giveCoordinate(1),this->giveNode(2)->giveCoordinate(2)) ;
  Mu::Point *p2 = new Mu::Point(this->giveNode(3)->giveCoordinate(1),this->giveNode(3)->giveCoordinate(2)) ;
  Mu::Point *p3 = new Mu::Point(this->giveNode(4)->giveCoordinate(1),this->giveNode(4)->giveCoordinate(2)) ;

  // transform to local coordinate of the parent element
  Mu::Point *localP0 = this->giveFEInterpolation()->global2Local(domain,this->giveNodeArray(),p0);
  Mu::Point *localP1 = this->giveFEInterpolation()->global2Local(domain,this->giveNodeArray(),p1);
  Mu::Point *localP2 = this->giveFEInterpolation()->global2Local(domain,this->giveNodeArray(),p2);
  Mu::Point *localP3 = this->giveFEInterpolation()->global2Local(domain,this->giveNodeArray(),p3);

  DelaunayTree *dt = new DelaunayTree(localP0,localP1,localP2) ;
  dt->insert(localP3);
  // insert the intersection points into the DT
  for(size_t j = 0 ; j < intersectPoints.size() ; j++)
  {
	 dt->insert(standardFEInterpolation->global2Local(domain,this->giveNodeArray(),&intersectPoints[j]));
  }

  // do the Delaunay triangulation ...
  vector<DelaunayTriangle *> *tri = dt->getTriangles(); 

  delete p0 ; delete p1 ; delete p2 ; delete p3 ; 
  //delete dt ;

  return tri;
}

double PlateIso :: computeAreaAboveEnrItem(EnrichmentItem *enrItem)
// **************************************************************
// The enrItem splits PlateIso4 into two portions : compute the area of the
// portion above the enrItem.
// First, define the intersection points, insert into the DelaunayTree
// perform the tesselation => triangles => area
// Latest version : 2005-09-07, when solve with structured Q4 mesh. 
{
  GeometryEntity* myGeo = enrItem->giveMyGeo() ;
  std::vector<Mu::Point> pts,intersect,temp ;
  
  // compute the intersection between element edges and enrItem
  // -------------------------------------------------------------------------
  for(size_t i = 0 ; i < this->giveNumberOfNodes() ; i++)
  {
	 Mu::Point *p = this->giveNode(i+1)->makePoint() ;
	 pts.push_back(*p);
	 delete p ;
  }
  for(size_t i = 0 ; i < pts.size() ;  i++)
  {
	 Mu::Segment s(pts[i],pts[(i+1)%pts.size()]) ;
	 if(myGeo->intersects(&s))
	 {
		temp = myGeo->intersection(&s) ;
		intersect.insert(intersect.end(),temp.begin(),temp.end()) ;
	 }
  }
  // -------------------------------------------------------------------------

  // If crack touched element edge , the above code just gives one intersection !!!
  // -------------------------------------------------------------------------
  PiecewiseLinear *geo = static_cast<PiecewiseLinear*>(enrItem->giveMyGeo());
  if (intersect.size() == 1) 
  {
	 std::cout << " Crack tips touched element edge !!! " << std::endl ;
	 std::list<Mu::Point*> *vertices = geo->giveMyListOfVertices() ;
	 std::list<Mu::Point*>::iterator j = vertices->end();
	 Mu::Point *lastVertex = *(--j) ;
	 lastVertex->print() ;  std::cout << std::endl ;   // debug only
	 for(size_t k = 0 ; k < pts.size() ;  k++)
	 {
		Mu::Segment s(pts[k],pts[(k+1)%pts.size()]) ;
		s.print() ;
		if(s.on(lastVertex))
		{
		  std::cout << " found edge containing tip  " << endl ;
		  intersect.push_back(*lastVertex);
		  k = pts.size() ;
		}
		//for( std::list<Mu::Point*>::iterator j = vertices->begin() ; j != vertices->end() ; j++)
		//  if(s.on(*j))
		//  {
		//	 intersect.push_back(**j);
		//	 j = vertices->end() ;     // stop, no test anymore.
		//  }
	 }
  }               // end of if (intersect.size() == 1) 
  // -------------------------------------------------------------------------

  // just do for the case enrItem completely splits the element => always have 
  // 2 intersection points
  if (intersect.size() != 2)
  {
	 std::cout << "enrItem must completely splits the element !!! " << std::endl;
	 assert(false);
  }

  Mu::Point *p1 = this->giveNode(1)->makePoint();
  Mu::Point *p2 = this->giveNode(2)->makePoint();
  Mu::Point *p3 = this->giveNode(3)->makePoint();
  Mu::Point *p4 = this->giveNode(4)->makePoint();
  DelaunayTree *dt = new DelaunayTree(p1,p2,p3);
  // insert points into DelaunayTree
  dt->insert(p4) ;
  dt->insert(&intersect[0]) ;
  dt->insert(&intersect[1]) ;

  // looking for kink points !!! 2005-09-18
  std::list<Mu::Point*> *vertices = geo->giveMyListOfVertices();
  if(vertices->size() != 2)
  {
	 for(std::list<Mu::Point*>::iterator i = vertices->begin() ; i != vertices->end(); i++)
	 {
		if(this->isWithinMe(*i) && (isAligned(*i,&intersect[0],&intersect[1]) == false))
		  dt->insert(*i) ;
	 }
  }

  double x1 = intersect[0].x ; // first intersection point
  double y1 = intersect[0].y ; // first intersection point

  double x2 = intersect[1].x ; // second intersection point
  double y2 = intersect[1].y ; // second intersection point

  vector<DelaunayTriangle *> *v  = dt->getTriangles(); 
  double A = 0.0  ; double delta ;
  // loop on delaunay triangles and compute the accumulated area
  for(size_t i = 0 ; i < v->size() ; i++)
  {
	 if((*v)[i]->isTriangle == true)
	 {
		Mu::Point *center = (*v)[i]->getCenter();
		if(x2 > x1)
		  delta = (x1-center->x)*(y2-center->y)-(x2-center->x)*(y1-center->y);
		else
		  delta = (x2-center->x)*(y1-center->y)-(x1-center->x)*(y2-center->y);

		if(delta > 0.000001) 
		  A  += (*v)[i]->area(); 
	 }
  }
  delete p1 ; delete p2 ; delete p3 ; delete p4 ;
  delete dt ; delete v ;  // Purify,14-10-2005

  return A ;
 
}

