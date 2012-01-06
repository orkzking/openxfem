/************************************************************************************ 

Copyright (C) 2005
Stephane BORDAS, Cyrille DUNANT, Vinh Phu NGUYEN, Quang Tri TRUONG, Ravindra DUDDU

This file is part of the XFEM C++ Library (XFEMLIB) written 
and maintained by above authors.

This program is free software; you can redistribute it and/or modify it.

This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
license text for more details.

Any feedback is welcome. Emails : nvinhphu@gmail.com, ...

*************************************************************************************/

#include "cracktip.h"

#include "flotarry.h"
#include "intarray.h"
#include "flotmtrx.h"
#include "vertex.h"
#include "node.h"
#include "element.h"
#include "material.h"
#include "gausspnt.h"
#include "feinterpol.h"
#include "geometryentity.h"
#include "geometry_2D.h"
#include "auxiliaryfield.h"
#include "auxiliaryfield.cpp"
#include "assert.h"
#include "delaunay.h"
#include "timestep.h"
#include "crackgrowthdirectionlaw.h"
#include "crackgrowthincrementlaw.h"
#include "crackinterior.h"
#include "timestep.h"
#include "timinteg.h"
#include <vector>
#include <list>
#include <valarray>
#include <cmath>
#include <iostream>

#ifndef M_PI 
#define M_PI  3.14159265358979323846
#endif


CrackTip :: CrackTip(int n,Domain* aDomain)
: EnrichmentItem(n,aDomain)
// ****************************************
{
  K_i  = 0.0 ;
  K_ii = 0.0 ;
  K_eq = 0.0 ;
  enrichRadius = 0 ;
  isActive = true ;
  field = this->giveMyFieldForSIF();
  this->crackTypeInitialization();
  matArray = NULL ;
  tipSegment = NULL ;
  stillInTipElement = false ; 
  myAssociatedCrackInterior = NULL ;
  isOnElementEdge = false ;
  alpha = 1000.0 ; // code convention
}

CrackTip :: ~CrackTip()
// ********************
{
  delete matArray ;
  delete tipSegment ;
}

Vertex* CrackTip ::giveTip()
// *************************
// returns the tip of the receiver
// used to check tip elements, compute asymptotic functions, ...
{
  return dynamic_cast<Vertex*>(myGeometry);  
}

void  CrackTip::crackTypeInitialization()
//***************************************
// Read from the input file the type of the crack
{
  char  type[32] ;

  this -> readString("Type",type) ;

  if (! strcmp(type,"HomoElast"))
	 tipID = HomoElast;
  else if (! strcmp(type,"BiMatElast"))
	 tipID = BiMatElast;
  else
	 assert(false);
}

FieldType CrackTip::giveMyFieldForSIF()
// ************************************
// For 2D, returns the field of the Domain, for 3D ...
{
  field = domain->giveFieldType();

  return field;
}

Mu::Circle* CrackTip::makeBall()
// *****************************
// make a ball centered at the TIP with radius r = fac * sqrt(area of tip element)
// corrected version : 2005-08-11
{
  // get the radius
  double rad = this->giveRadiusOfDomainIntegration();

  // define the center 
  double X = this->giveTip()->giveCoordinate(1); 
  double Y = this->giveTip()->giveCoordinate(2); 
  Mu::Point  *center = new Mu::Point(X,Y);

  return new Mu::Circle(rad,center);
}

std::list<Element*> CrackTip::buildIntegrationDomain1()
// ****************************************************
// IntegrationDomain is a set of elements that belong to the ball 
// of radius r = h*sqrt(area_tip_element) centered at the crack tip.
// Remark: This integration domain is no longer used !!!
{
  // find elements belongs to the circle.

  vector<Element*> *interactedElements = this->giveElementsInteractWithMe();
  Element *tipElem = (*interactedElements)[0] ; 
  
  Mu::Circle *ball = this->makeBall();
  std::list<Element*> ret = tipElem->conflicts(ball);

  // after call Element::conflicts(Circle*), remember to reset checked = false for the next time
  for(size_t i = 0 ; i < domain->giveNumberOfElements() ; i++)
	 domain->giveElement(i+1)->clearChecked();

  delete ball ;

  return ret;
} 

std::list<Element*> CrackTip::buildIntegrationDomain2()
// ****************************************************
// IntegrationDomain is a set of elements that intersect with the ball 
// of radius r = h*sqrt(area_tip_element) centered at the crack tip.
{
  // find elements belongs to the circle.

  vector<Element*> *interactedElements = this->giveElementsInteractWithMe();
  Element *tipElem = (*interactedElements)[0] ; 
 
  Mu::Circle *ball = this->makeBall();

  std::list<Element*> temp = tipElem->conflicts(ball);

  // after call Element::conflicts(Circle*), remember to reset checked = false for the next time
  for(size_t i = 0 ; i < domain->giveNumberOfElements() ; i++)
	 domain->giveElement(i+1)->clearChecked();

  // find out elements intersecting with the perimeter of the circle
  // just loop over temp to check, so it is efficient.
  std::list<Element*> ret ;
  for(std::list<Element*> :: iterator i = temp.begin() ; i != temp.end() ; i++)
  {
	 if( (*i)->intersects(ball) )
		ret.push_back(*i);
  }

  delete ball ;

  return ret;
} 


std::valarray<double> CrackTip::computeInteractionIntegral(TimeStep* stepN)
// ************************************************************************
// compute the interaction integral, see Yau et al., 1980 
// Steps :
// 1- detection of the elements on which we integrate
// 2- loop over these elements
// 3- loop over Gauss points
// 4- computation of stress, strain... in local coordinates !!!   ATTENTION 
// 5- computation of the auxilliary fields: AuxStress and AuxEps and AuxGradDisp
// 6- computation of I1 and I2
{
  FloatArray   *sigma,*epsilon,*U ;
  FloatMatrix  *dNdx,*D,*B,*rdu1dx,*rt,*LocalGradDisp ;
  Node         *aNode ; 
  Mu::Point    *point,*coord,*globalCoordGP ;
  GaussPoint   *gp;
  GaussPoint** gaussPointArray;
  IntArray     *nodeArr ;
  double       dV ;
  size_t       numOfNodes ;
  double       dqdx1, dqdx2 ;
  std::valarray<double> I(0.,2) ;         // interaction integral I(1,2) of modes I and II
  Mu::Circle *myBall = this->makeBall();  // circle centered at tip, r = r_d
 
  // Choose domain used for standard XFEM and modified XFEM (crack nucleation)
  std::list<Element*> intDomain = this->buildIntegrationDomain2();

  /*
  // ******************************************************************************
  // ***                       EXPORT TO MATLAB FOR PLOT                        ***
  // ******************************************************************************

  // plot the elements in domain A to check. 
  FILE       *mlbFile;
  char        mlbFileName[64] = "integration" ;    // name of the file
  strcat(mlbFileName,".m");                        // extension .m ( Matlab M file)
  mlbFile = fopen(mlbFileName, "w");

  fprintf (mlbFile,"aroundElements = [ ");
  for(std::list<Element*> :: iterator i = intDomain.begin() ; i != intDomain.end() ; i++)
  {
  fprintf (mlbFile,"%d \n",(*i)->giveNumber());
  }
  fprintf (mlbFile," ];\n ");
  
  // weight functions at nodes ...
  fprintf (mlbFile,"weight = [ ");
  double a = myAssociatedCrackInterior->giveCrackLength();
  std::cout << " Half crack length " << a << endl;
  for(std::list<Element*> :: iterator i = intDomain.begin() ; i != intDomain.end() ; i++)
  {
	 std::valarray<double> q(0.,(*i)->giveNumberOfNodes()) ;  // weight function q 
	 if(this->domain->giveNumberOfInteractedElements()!=0)
	 {
		for(size_t m = 0 ; m < (*i)->giveNumberOfNodes() ; m++)
		{
		  aNode = (*i)->giveNode(m+1);
		  point = aNode->makePoint();
		  if (myBall->in(point))                          // p locates inside the ball
			 q[m] = 1;
		  delete point ;
		}
	 }
	 else
	 {
      double X = this->giveTip()->giveCoordinate(1); 
		double Y = this->giveTip()->giveCoordinate(2); 
		Mu::Point  *tip = new Mu::Point(X,Y);
		for(size_t m = 0 ; m < (*i)->giveNumberOfNodes() ; m++)
		{
		  aNode = (*i)->giveNode(m+1);
		  point = aNode->makePoint();
		  double r = sqrt(dist(*point,*tip));

		  if(r<=0.25*a)
			 q[m] = 1;
		  else if(r<=a && r>0.25*a)
			 q[m] = 4./3.*(1-r/a) ;
		  else
			 q[m] = 0;

		  delete point ;
		}
	 }
	 // export to Matlab 
	 if( (*i)->giveNumberOfNodes() == 3 )
	 fprintf (mlbFile,"%5.1e %5.1e %5.1e \n",q[0],q[1],q[2]);
	 else if( (*i)->giveNumberOfNodes() == 4 )
	 fprintf (mlbFile,"%5.1e %5.1e %5.1e %5.1e \n",q[0],q[1],q[2],q[3]);
	 else
	 std::cout << " reserve anyone who will do other element types !!! " << std::endl ;
  }
  fprintf (mlbFile," ];\n ");

  fprintf (mlbFile," GP = [ "); // Gauss points for J integral

  // --------------------------------------------------------------------------------
  */

  // inclination angle and rotation matrix
  
  FloatMatrix *r = new FloatMatrix(2,2);
  r->at(1,1) =  cos(alpha) ; r->at(1,2) = sin(alpha) ;
  r->at(2,1) = -sin(alpha) ; r->at(2,2) = cos(alpha) ; 

  rt = r -> GiveTransposition() ;  

 // initialization of used auxiliary field
  AuxiliaryField<Material,NullMaterial,PlaneStrain> *auxFieldHomo 
	 = new AuxiliaryField<Material,NullMaterial,PlaneStrain>();

  AuxiliaryField<Material,Material,PlaneStrain> *auxFieldBiMat 
	 = new AuxiliaryField<Material,Material,PlaneStrain>();

  // ******************************************************************************
  // ***               LOOP ON ELEMENTS IN THE INTEGRATION DOMAIN               ***
  // ******************************************************************************
 
  for(std::list<Element*> :: iterator i = intDomain.begin() ; i != intDomain.end() ; i++)
  {
	 D = (*i)->giveConstitutiveMatrix()->GiveCopy();   // compliance matrix
	 U = (*i)->ComputeVectorOf('d',stepN);             // displacement vector 
	 nodeArr  = (*i)->giveNodeArray() ;                // connectivity of element elem
	 gaussPointArray = (*i)->setGaussQuadForJ_Integral(); // gauss points of elem
	 numOfNodes = (*i)->giveNumberOfNodes() ;   // number of nodes of elem

	 // computation of weight function q 
	 // Attention : need to rewrite for quadratic XFEM !!!
	 std::valarray<double> q(0.,numOfNodes) ;          // weight function q 
	 
		for(size_t m = 0 ; m < numOfNodes ; m++)
		{
		  aNode = (*i)->giveNode(m+1);
		  point = aNode->makePoint();
		  if (myBall->in(point))                          // p locates inside the ball
			 q[m] = 1;
		  delete point ;
		}
	 /*else
	 {
      double X = this->giveTip()->giveCoordinate(1); 
		double Y = this->giveTip()->giveCoordinate(2); 
		Mu::Point  *tip = new Mu::Point(X,Y);
		for(size_t m = 0 ; m < numOfNodes ; m++)
		{
		  aNode = (*i)->giveNode(m+1);
		  point = aNode->makePoint();
		  double r = sqrt(dist(*point,*tip));

		  if(r<=0.25*a)
			 q[m] = 1;
		  else if(r<=a && r>0.25*a)
			 q[m] = 4./3.*(1-r/a) ;
		  else
			 q[m] = 0;

		  delete point ;
		}
	 }*/

	 // ******************************************************************************
	 // ***     LOOP ON GAUSS POINTS OF EACH ELEMENT IN THE INTEGRATION DOMAIN     ***
	 // ******************************************************************************

	 for(size_t p = 0 ; p < (*i)->giveNumOfGPtsForJ_Integral() ; p++)
	 {
		gp    = gaussPointArray[p]; 
		coord = gp -> giveCoordinates() ; 
		dNdx  = (*i)->giveFEInterpolation()->evaldNdx(domain,nodeArr,coord);
		dV    = (*i)->computeVolumeAround(gp);

		// COMPUTATION OF DERIVATIVES OF U W.R.T x, dU1/dx
		B = (*i)->ComputeBmatrixAt(gp);
		FloatMatrix  *dU1dx = new FloatMatrix(2,2);// stores derivative of u1 w.r.t x,y 
		for(size_t k = 1 ; k <= 0.5*B->giveNumberOfColumns() ; k++)
		{
		  dU1dx->at(1,1) += B->at(1,2*k-1)*(*U)[2*k-2];   // du1/dx
		  dU1dx->at(1,2) += B->at(2,2*k)  *(*U)[2*k-2];   // du1/dy
		  dU1dx->at(2,1) += B->at(1,2*k-1)*(*U)[2*k-1];   // du2/dx
		  dU1dx->at(2,2) += B->at(2,2*k)  *(*U)[2*k-1];   // du2/dy
		}

		// COMPUTE DERIVATIVES OF WEIGHT FUNCTION Q 

		//if( this->domain->giveNumberOfInteractedElements() != 0 )
		//{
		  double dqdx = 0.0 ; double dqdy = 0.0 ; 
		  for(size_t m = 0 ; m < numOfNodes ; m++)
		  {
			 dqdx += dNdx->at(m+1,1)*q[m];  // derivative of q w.r.t global x
			 dqdy += dNdx->at(m+1,2)*q[m];  // derivative of q w.r.t global y
		  }
		  dqdx1 =  dqdx * cos(alpha) + dqdy * sin(alpha) ;
		  dqdx2 = -dqdx * sin(alpha) + dqdy * cos(alpha) ;
		//}
		/*else              // compute modified weight function(crack nucleation)
		{
		  globalCoordGP = (*i)->giveFEInterpolation()->local2Global(domain,nodeArr,coord);
		  FloatArray * localCoord = this->computeLocalCoordOf(globalCoordGP) ;

		  double xloc = (*localCoord)[0];
		  double yloc = (*localCoord)[1];
		  double R = sqrt( xloc * xloc + yloc * yloc ); 
		  double theta = atan2(yloc,xloc);

		  double dqdr; 
		  if( R >= 0.25*a && R <= a)
		  {
			 std::cout << " Nonzeo value !!! " << endl ;
			 dqdr = -4/(3*a);
		  }
		  else
			 dqdr = 0.0 ;
		  dqdx1 = dqdr * cos(theta);
		  dqdx2 = dqdr * sin(theta);
		}*/

		// COMPUTE THE STRESS OF ELEM AT GAUSS POINT GP
		epsilon = B -> Times(U) ;
		sigma   = D->Times(epsilon);

		// TRANSFORM THESE TERMS OF STATE (1) TO LOCAL CRACK TIP COORDINATE SYSTEM
		// using formula : (rot_matrix)*H*(rot_matrix)' 

		rdu1dx = r -> Times(dU1dx) ;                        
		LocalGradDisp = rdu1dx -> Times(rt) ;           

		FloatArray *LocalSigma = new FloatArray(3);
		double C2T = cos(2.*alpha) ;
		double S2T = sin(2.*alpha) ;

		(*LocalSigma)[0] = 0.5*((*sigma)[0]+(*sigma)[1])+
		  0.5*((*sigma)[0]-(*sigma)[1])*C2T + (*sigma)[2]*S2T ;

		(*LocalSigma)[1] = 0.5*((*sigma)[0]+(*sigma)[1])-
		  0.5*((*sigma)[0]-(*sigma)[1])*C2T - (*sigma)[2]*S2T ;

		(*LocalSigma)[2] = (*sigma)[2]*C2T - 0.5*((*sigma)[0] - (*sigma)[1])*S2T ;

		// COMPUTATION OF COMPONENTS OF THE AUXILIARY FIELD

		FloatArray  AuxEps(3);
		FloatArray  AuxStress(3);
		FloatMatrix AuxGradDisp(2,2) ;   

		globalCoordGP = (*i)->giveFEInterpolation()->local2Global(domain,nodeArr,coord);

		// loop on modes 
		for(size_t u = 0 ; u < 2 ; u++)
		{
		  if( matArray->size() == 1 ) // homogeneous media
		  {
			 if(u == 0)
				auxFieldHomo->ComputeComponentsOfAuxField(this,globalCoordGP,Mode_i,AuxGradDisp,AuxEps,AuxStress); 
			 else
				auxFieldHomo->ComputeComponentsOfAuxField(this,globalCoordGP,Mode_ii,AuxGradDisp,AuxEps,AuxStress); 
		  }
		  else                        // bimaterial media
		  {
			 if(u == 0)
			 {
				 auxFieldBiMat->ComputeComponentsOfAuxField(this,globalCoordGP,Mode_i,AuxGradDisp,AuxEps,AuxStress); 
			 }
			 else
				auxFieldBiMat->ComputeComponentsOfAuxField(this,globalCoordGP,Mode_ii,AuxGradDisp,AuxEps,AuxStress); 
		  }

		  // compute the interaction strain energy W(1,2)= sigma_{ij}*epsilon_{ij}		
		  double W = (*LocalSigma)[0]*AuxEps.at(1) + 2.0*(*LocalSigma)[2]*AuxEps.at(3) + 
			 (*LocalSigma)[1]*AuxEps.at(2) ;

		  // computing terms in the interaction integral ...

		  double I1 =
			 ( (*LocalSigma)[0]*AuxGradDisp.at(1,1)+(*LocalSigma)[2]*AuxGradDisp.at(2,1) ) * dqdx1 +
			 ( (*LocalSigma)[2]*AuxGradDisp.at(1,1)+(*LocalSigma)[1]*AuxGradDisp.at(2,1) ) * dqdx2 ;

		  double I2 =
			 ( AuxStress.at(1) * LocalGradDisp->at(1,1) + AuxStress.at(3) * LocalGradDisp->at(2,1) ) * dqdx1 +
			 ( AuxStress.at(3) * LocalGradDisp->at(1,1) + AuxStress.at(2) * LocalGradDisp->at(2,1) ) * dqdx2 ;

		  I[u] += (I1 + I2 - W*dqdx1) * dV ;	
		}		

		delete dNdx ; delete dU1dx ;  delete epsilon; delete sigma; delete LocalGradDisp;
		delete B ; delete LocalSigma ; delete rdu1dx ;

		// export Gauss points to Matlab for plot ...
		//fprintf (mlbFile,"%5.6f %5.6f  \n",globalCoordGP->x,globalCoordGP->y);

	 }                               // end of loop on Gauss points

	 // delete GPs used for J integral  
	 for (size_t j = 0 ; j < (*i)->giveNumOfGPtsForJ_Integral() ; j++)
		delete gaussPointArray[j] ;
	 delete gaussPointArray ;
	 delete D;
  }                                 // end of loop on elements 
  //fprintf (mlbFile," ];\n ");
  // -----------------------------------------------------------------------------------------------
  // EXPORT MATLAB PLOT COMMANDS 
  /*
  fprintf (mlbFile,"figure;\n ");
  fprintf (mlbFile,"hold on ;\n") ;
  fprintf (mlbFile,"set(gcf,'color','white');\n ");
  fprintf (mlbFile,"plot_mesh(node,element,elementType,'g-');\n") ;
  fprintf (mlbFile,"plot_mesh(node,element(aroundElements,:),elementType,'g-');\n") ;
  fprintf (mlbFile,"r = %6.3e ;\n",myBall->getRadius());   
  fprintf (mlbFile,"theta = -pi:0.1:pi;\n");   
  fprintf (mlbFile,"x = %6.3e + r*cos(theta) ;\n",myBall->getCenter()->x);   
  fprintf (mlbFile,"y = %6.3e + r*sin(theta) ;\n",myBall->getCenter()->y);  

  // plot the crack ...
  fprintf (mlbFile,"for i = 1 : size(crack,1)\n");
  fprintf (mlbFile,"  crFig = plot(crack(:,1:2:size(crack,2)),crack(:,2:2:size(crack,2)),'y-');\n");
  fprintf (mlbFile,"  set(crFig,'LineWidth',1.5)\n");
  fprintf (mlbFile,"end\n");
  fprintf (mlbFile,"plot(x,y);\n") ;
  fprintf (mlbFile,"plot(GP(:,1),GP(:,2),'r*');\n") ;
  fprintf (mlbFile,"for i = 1 : size(crack,1)\n");
  fprintf (mlbFile,"  crFig = plot(crack(:,1:2:size(crack,2)),crack(:,2:2:size(crack,2)),'y-');\n");
  fprintf (mlbFile,"  set(crFig,'LineWidth',1.5)\n");
  fprintf (mlbFile,"end\n");
  fprintf (mlbFile,"axis off;\n") ;

  fprintf (mlbFile,"figure;\n ");
  fprintf (mlbFile,"hold on ;\n") ;
  fprintf (mlbFile,"set(gcf,'color','white');\n ");
  fprintf (mlbFile,"plot_mesh(node,element,elementType,'g-');\n") ;
  fprintf (mlbFile,"plot_field(node,element(aroundElements,:),elementType,weight);\n") ;
  fprintf (mlbFile,"for i = 1 : size(crack,1)\n");
  fprintf (mlbFile,"  crFig = plot(crack(:,1:2:size(crack,2)),crack(:,2:2:size(crack,2)),'y-');\n");
  fprintf (mlbFile,"  set(crFig,'LineWidth',1.5)\n");
  fprintf (mlbFile,"end\n");
  fprintf (mlbFile,"colorbar;\n");
  fprintf (mlbFile,"axis off;\n") ;
  fclose  (mlbFile) ;
  // -----------------------------------------------------------------------------------------------
  */
  delete auxFieldHomo ; delete r ; delete rt ;

  return I ;
}

double CrackTip::computeSIFMultiplierForOneMat()
// *********************************************
// Ki = this->computeSIFMultiplierForOneMat() * Interaction integral
{
  this->giveMatArray();

  double E  = (*matArray)[0]->give('E');     
  double nu = (*matArray)[0]->give('n');

  if (field == PlaneStress)
	 return (E /2.0 );
  else if (field == PlaneStrain)
	 return E/(2.0 - nu*nu - nu*nu);
  else
	 assert(false);
}

double CrackTip::computeSIFMultiplierForBiMat()
// ********************************************
// The same as CrackTip::computeSIFMultiplierForOneMat() but for
// interfacial cracks 
{
  this->giveMatArray();

  double E1  =  (*matArray)[0]->give('E');      
  double nu1 =  (*matArray)[0]->give('n');

  double E2  =  (*matArray)[1]->give('E');      
  double nu2 =  (*matArray)[1]->give('n');

  double E1bar,E2bar; 
  if (field == PlaneStress)
  {
	 E1bar = E1;
	 E2bar = E2;
  }
  else if (field == PlaneStrain)
  {
	 E1bar = E1/(1. - nu1 * nu1);
	 E2bar = E2/(1. - nu2 * nu2);
  }
  else
	 assert(false);

  double Eetoile = (E1bar*E2bar)/(E1bar+E2bar); 

  double shear1 = E1/(1. + nu1 + nu1);
  double shear2 = E2/(1. + nu2 + nu2);

  double K1 = Field<PlaneStrain>::K(nu1);     // the Kolosov coeficient 
  double K2 = Field<PlaneStrain>::K(nu2);     

  double epsilon = 0.5/M_PI*log((shear1 + shear2*K1)/(shear1*K2 + shear2));

  return Eetoile * cosh(M_PI*epsilon) * cosh(M_PI*epsilon);
}
double CrackTip::computeSIFMultiplier()
// *************************************
// K_i = this->computeSIFMultiplier()*this->computeInteractionIntegral()
{
  if (tipID == HomoElast)
	 return this->computeSIFMultiplierForOneMat();
  else if (tipID == BiMatElast)
	 return this->computeSIFMultiplierForBiMat();
  else
	 assert(false);
}

void CrackTip::computeSIFs(TimeStep* stepN)
// *****************************************
// compute the stress intensity factors
// using the formula : K_i = (E/2.0)* I(1,1)(one material)
{
  double A = computeSIFMultiplier();
  std::valarray<double> I = this->computeInteractionIntegral(stepN);

  K_i  = A * I[0] ;
  K_ii = A * I[1] ;
}

void CrackTip::computeK_eq(TimeStep* stepN)
// ****************************************
// compute the equivalent stress intensity factor 
{
  K_eq = sqrt(K_i*K_i + K_ii*K_ii);
}

std::vector<Material*>*  CrackTip::giveMatArray()
// **********************************************
// Read from the input file the material info
{
  if (!matArray)
  {
	 matArray = new std::vector<Material*> ;    // initialization
	 size_t nMats  = this->readInteger("Mat") ; // number of materials

	 for (size_t i = 0 ; i< nMats ; i++)
	 {
		size_t matNumber = this->readInteger("Mat",i+2) ;
		matArray->push_back(domain->giveMaterial(matNumber));
	 }
  }

  return matArray ;
}

Mu::Segment* CrackTip :: giveTipSegment()
// **************************************
{
  if(!tipSegment)
  {
	 Vertex *tip = this->giveTip();
	 PiecewiseLinear *myDead = tip->giveMyParent();// myDead = PiecewiseLinear contains this TIP.
	 tipSegment = myDead->giveSegmentContain(tip); // segment contains the TIP
  }

  return tipSegment ;
}
/*
FloatArray* CrackTip ::computeLocalCoordOf(Mu::Point *p)
// *****************************************************
// computes the coordinates of p in the local crack tip coordinate system.
// Used then to define r and theta . 
{

double x1 = this->giveTipSegment()->first()->x ;
double y1 = this->giveTipSegment()->first()->y ;
double x2 = this->giveTipSegment()->second()->x ;
double y2 = this->giveTipSegment()->second()->y ;

// checking the Point p is above or below the extended segment seg
double delta = (x1 - p->x)*(y2 - p->y) - (x2 - p->x)*(y1 - p->y) ;

double tol = std::numeric_limits<float>::epsilon(); // tol = 1.19209e-007
int sign;
if (delta > tol)
sign =  1;
else if (delta < -tol)
sign = -1;
else
sign = 0 ;

Mu::Point v(x2 - x1, y2 - y1); 
Mu::Point w(p->x - x1, p->y - y1);

double c1 = w*v; 
double c2 = v*v;
double b = c1 / c2;
Mu::Point Pb(x1 + b * v.x, y1 + b * v.y) ;
double yloc = dist(Pb, (*p))*sign ;   // distance from the Point p to the extended crack segment

Mu::Point w1(p->x - x2, p->y - y2);
FloatArray * ret = new FloatArray(2);

double xloc = dist(Pb, *(this->giveTipSegment()->second())) ; // the x coord. in local carck tip coordinate system
double c = w1*v; 
if ( c >= 0 )     // the angle is <=90 
xloc = xloc ;
else
xloc = -1. * xloc ;

(*ret)[0] = xloc; // the x coord. in local crack tip coordinate system
(*ret)[1] = yloc; // the y coord. in local crack tip coordinate system

return ret;
}*/

FloatArray* CrackTip ::computeLocalCoordOf(Mu::Point *p)
// *****************************************************
// computes the coordinates of p in the local crack tip coordinate system.
// Used then to define r and theta . 
// New version. 2005-08-09
{
  this->giveTipSegment();
  double x1 = tipSegment->first()->x ;
  double y1 = tipSegment->first()->y ;
  double x2 = tipSegment->second()->x ;
  double y2 = tipSegment->second()->y ;

  double alpha = atan2(y2-y1,x2-x1); // inclination of local coord. system

  double xloc =  (p->x - x2) * cos(alpha) + (p->y - y2) * sin(alpha) ;
  double yloc = -(p->x - x2) * sin(alpha) + (p->y - y2) * cos(alpha) ;

  FloatArray * ret = new FloatArray(2);
  (*ret)[0] = xloc; // the x coord. in local crack tip coordinate system
  (*ret)[1] = yloc; // the y coord. in local crack tip coordinate system

  return ret;
}

 std::valarray<double>* CrackTip ::computePolarCoordOf(Mu::Point *p)
// *****************************************************************
// computes the polar coordinates of p in the local crack tip coordinate system.
{
  FloatArray * localCoord = this->computeLocalCoordOf(p) ;

  double xloc = (*localCoord)[0];
  double yloc = (*localCoord)[1];

  std::valarray<double> *ret = new std::valarray<double>(2);
  (*ret)[0] = sqrt( xloc * xloc + yloc * yloc ); // r
  (*ret)[1] = atan2(yloc,xloc) ;                 // theta

  delete localCoord;

  return ret;
}

void CrackTip ::printOutputAt(TimeStep *stepN, FILE* strFile, FILE* s01File)
// *************************************************************************
// compute the SIFs and print to the file *.str
{
  double r = this->giveRadiusOfDomainIntegration();

  this->computeSIFs(stepN) ;

  std::cout << " SIF of Mode I  : " << K_i << std::endl ;

  std::cout << " SIF of Mode II : " << K_ii << std::endl ;

  fprintf (strFile,"    Stress Intensity Factors results : \n") ;
  fprintf (strFile,"-----------------------------------------\n") ;
  fprintf (strFile,"     Mode I       Mode II    Radius of int. domain: \n") ;
  fprintf (strFile,"Crack tip %d :",number) ;
  fprintf (strFile,"SIFs  % .4f % .4f % .4f \n",K_i,K_ii,r) ;
} 

double CrackTip::giveRadiusOfDomainIntegration()
// *********************************************
// returns the radius of the domain integration.
// Two cases:
// 1. Normal enrichment scheme(enrDetectorID=1): rd=fac*sqrt{tip-element area}
// where fac is read from input file : domainIntFac
// 2. Fixed enrichment scheme(enrDetectorID=2): rd=radii
// where radii is read from the input file : domainIntRad
// Modified version : 13-10-2005.
{
  double radius ;

  size_t enrDetectorID = this -> readInteger("enrichScheme");
  if(enrDetectorID == 1)
  {
	 double fac = this->read("domainIntFac");
	 vector<Element*> *interactedElements = this->giveElementsInteractWithMe();
	 Element *tipElem = (*interactedElements)[0] ; 
	 double he = tipElem->area();
	 radius = fac*sqrt(he) ;       // rd=fac*sqrt{tip-element area}
  }
  else if(enrDetectorID == 2)
  {
	 radius = this->read("domainIntRad");
  }
  else
	 assert(false) ;

  return radius ;
}


double CrackTip::giveEnrichRadius()
// ********************************
// Reads this valuefrom the data file if it does not exist.
{
  if (enrichRadius == 0)
  {
	 enrichRadius = this->read("enrichRadius");
  }
  return enrichRadius ;
}


void CrackTip :: UpdateMyGeometry()
// ********************************
// updates position of the tip.
{
  // compute the crack advance : growth angle theta and crack increment length delta
  double delta = domain->giveGrowthIncrementLaw()->computeIncrementLength();
  double theta = domain->giveGrowthDirectionLaw()->computeGrowthAngle(this);

  double deltaX = delta * cos(theta + alpha); // in global coordinate system
  double deltaY = delta * sin(theta + alpha); // in global coordinate system

  // coord. of current tip 
  double currentX = this->giveTip()->giveCoordinate(1); 
  double currentY = this->giveTip()->giveCoordinate(2); 

  // coord. of the new tip
  double updatedX = currentX + deltaX ; 
  double updatedY = currentY + deltaY ; 

  // update the crack tip coordinates
  static_cast<Vertex*>(myGeometry)->setCoordinates(updatedX,updatedY);

  // update the tip segment 
  tipSegment->set(currentX,currentY,updatedX,updatedY);

  // check if the updated TIP still locates in tip-element or not
  std::vector<Element*> *elemList = this->giveElementsInteractWithMe(); 
  assert(elemList->size() == 1);  // a Tip just interacts with ONE element !!!
  Element *tipElement = (*elemList)[0];
  if( static_cast<Vertex*>(myGeometry)->interactsWith(tipElement) )
	 stillInTipElement = true ;	 
} 

Mu::Circle* CrackTip::DefineDomainForUpdatedEnrichment()
// *****************************************************
// Set of elements belonging to the ball centered at the old tip
// radius = crack increment length
{
  // coord. of current tip 
  double currentX = this->giveTip()->giveCoordinate(1); 
  double currentY = this->giveTip()->giveCoordinate(2); 

  Mu::Point center(currentX,currentY) ; // center of the circle

  double rad = domain->giveGrowthIncrementLaw()->computeIncrementLength();

  return new Mu::Circle(rad,center) ;
}

std::list<Element*> CrackTip::DefineUpdatedElements()
// *****************************************************
{
  std::vector<Element*> *elemList = this->giveElementsInteractWithMe(); 
  assert(elemList->size() == 1);
  Element *tipElement = (*elemList)[0];

  Mu::Circle* c = this->DefineDomainForUpdatedEnrichment() ;

  std::list<Element*> ret = tipElement->conflicts(c) ;

  delete c ;

  return ret ;
}

void CrackTip::updateEnrichment()
// ******************************
// Key method for crack growth simulation !!!
// 2005-09-02 : start 
// 2005-09-23 : corrected 
{
  std::cout << " UPDATE GEOMETRIES ... " << endl ;

  for(size_t i = 0 ; i < domain->giveNumberOfElements() ; i++)
	 domain->giveElement(i+1)->reinitializeStateOfElement() ;

  // get the old tip element 
  std::vector<Element*> *elemList = this->giveElementsInteractWithMe(); 
  assert(elemList->size() == 1);
  Element *oldTipElement = (*elemList)[0];
 
  // update the crack tip's geometry ...
  this->UpdateMyGeometry() ; 
  
  // ***********************************************************************
  // *****            Update the geo-mesh interaction                   ****
  // ***********************************************************************
  if(stillInTipElement)
  {
	 std::cout << " The tip still locates within old tip element !!! " << endl ;
	 oldTipElement->setStateOfElement() ;
    myAssociatedCrackInterior->UpdateMyGeometry();
	 stillInTipElement = false ; // Stephane's reminder :)
  }
  else
  {
	 std::cout << " The tip left old tip element !!! " << endl ;
	 // erase enrichment items of old tip-element before updating 
	 oldTipElement->eraseListOfEnrichmentItems();
	 // erase enrichment items of old tip-element's node before updating 
	 for(size_t i = 0 ; i < oldTipElement->giveNumberOfNodes() ; i++)
	 {         
		oldTipElement->giveNode(i+1)->eraseListOfEnrichmentItems() ;
	 }		

	 std::list<Element*> elems = this->DefineUpdatedElements();

	 // 1. FIND ELEMENTS INTERSECT WITH THE CRACK ADVANCE

	 std::list<Element*> newSplitElems ;

	 for(std::list<Element*> :: iterator i = elems.begin() ; i != elems.end() ; i++)
	 {
		if((*i)->intersects(tipSegment))
		  newSplitElems.push_back(*i) ;
	 }

	 // 2. FIND NEW TIP-ELEMENT
	 // 2.1. The tip falls within element
	 Element *newTipElem = NULL ;
	 std::list<Element*> :: const_iterator ci = newSplitElems.begin() ;
	 bool found = false ;
	 while( ci != newSplitElems.end() && !found)
	 {
		if(static_cast<Vertex*>(myGeometry)->interactsWith(*ci))
		{
		  newTipElem = *ci ;
		  found = true ;
		}
		ci++ ;
	 }
	 // 2.2 if the circle used to do the update is too small, i.e., 
	 // the crack advance is small, the code maybe did not find out the new tip element yet.
	 if(newTipElem == NULL)
	 {
		std::list<Element*> *supOfOldTipElem = oldTipElement->ComputeNeighboringElements();
		std::list<Element*> :: const_iterator ci = supOfOldTipElem->begin() ;
		while( ci != supOfOldTipElem->end() && !found)
		{
		  if(static_cast<Vertex*>(myGeometry)->interactsWith(*ci))
		  {
			 newTipElem = *ci ;
			 found = true ;
		  }
		  ci++ ;
		}
	 }
	 // 2.3 Tip touches element edge or goes through node !!!
	 // for some cases (new tip touches element edge or goes through a node), 
	 // do not find out new tip element yet. Need to do additional jobs. 2005-09-06.
	 if(newTipElem == NULL)
	 {
		std::cout << " The new tip does not locate within element " << std::endl ; 

		double X = this->giveTip()->giveCoordinate(1); 
		double Y = this->giveTip()->giveCoordinate(2); 

		Mu::Point *p = new Mu::Point(floor(X*1e6+0.5)/1e6,floor(Y*1e6+0.5)/1e6) ;

		std::list<Element*> :: iterator i = newSplitElems.begin() ;
		found = false ;
		while( (!found) && (i != elems.end()) )
		{
		  if( (*i)->IsOnEdge(p) )
		  {
			 newTipElem = *i ;
			 found = true ;
		  }
		  i++ ;
		}

		std::list<Element*> *supportOfTipElem = newTipElem->giveNeighboringElements() ;
		std::list<Element*> :: iterator j = supportOfTipElem->begin() ;
		found = false ;
		while( (!found) && (j != supportOfTipElem->end()) )
		{
		  if( (*j)->IsOnEdge(p) )
		  {
			 newTipElem = *j ;
			 found = true ;
		  }
		  j++ ;
		}
		delete p ;

		this->isOnElementEdge = true ;

		if(newTipElem == NULL)
		  assert(false) ;
	 }

	 // Attention, if crack advance goes throught the node, then Cyrille's code
	 // returns false while it must return true => need additional check .
	 /*
	 if( find(newSplitElems.begin(),newSplitElems.end(),oldTipElement)
	 == newSplitElems.end())
	 newSplitElems.push_back(oldTipElement); */

	 // For multi-material problems, the material of crack tip changed. Need to update this
	 // so that the auxiliary field correctly computed. 2005-09-21

	 if(this->domain->isMultiMaterialDomain())
	 {
		Material *newMat = newTipElem->giveMaterial() ;
		(*matArray)[0] = newMat ;
	 }

	 // ***********************************************************************
	 // *****                Update the enrichment                         ****
	 // ***********************************************************************
    std::cout << " UPDATE ENRICHMENT ... " << endl ;
	 // 1. Enrich this new tip-element with the updated CrackTip
	 newTipElem->isEnrichedWith(this);
	 this->interactedElements->clear() ; // the tip did not interact with old tip-element anymore.
	 this->setListOfInteractedElements(newTipElem);
	 newTipElem->setEnrichmentForMyNodes(this) ; // ATTENTION, FIXED ENRICHMENT AREA !!!

	 // 2. Enrich new split-elements with the tip's associated CrackInterior
	 for(std::list<Element*> :: iterator i = newSplitElems.begin() ; i != newSplitElems.end() ; i++)
	 {
		(*i)->isEnrichedWith(myAssociatedCrackInterior) ;
		myAssociatedCrackInterior->setListOfInteractedElements(*i);
		(*i)->setEnrichmentForMyNodes(myAssociatedCrackInterior) ;
	 }

	 // ***********************************************************************
	 // *****             Resolve conflicts in enrichment                  ****
	 // ***********************************************************************
    std::cout << " RESOLVE CONFLICTS IN NEW ENRICHED NODES ... " << endl ;
	 for(size_t i = 0 ; i < newTipElem->giveNumberOfNodes() ; i++)
	 {         
		newTipElem->giveNode(i+1)->resolveConflictsInEnrichment() ;
	 }		

	 // ***********************************************************************
	 // *****         Resolve linear dependency in enrichment              ****
	 // ***********************************************************************
    std::cout << " RESOLVE LINEAR DEPENDENCY FOR NEW STEP ENRICHED NODES ... " << endl ;
	 Node *aNode ;
	 std::vector<Node*> StepEnrichedNodes ;
	
	 for(std::list<Element*> :: iterator i = newSplitElems.begin() ; i != newSplitElems.end() ; i++)
	 {
		for(size_t j = 0 ; j < (*i)->giveNumberOfNodes() ; j++)
		{
		  aNode = (*i)->giveNode(j+1) ;
		  if( find(StepEnrichedNodes.begin(),StepEnrichedNodes.end(),aNode)
			 == StepEnrichedNodes.end() )
			 StepEnrichedNodes.push_back(aNode);
		}
	 }
	 // Update geometry of associated CrackInterior
	 myAssociatedCrackInterior->UpdateMyGeometry();
    /*
	 for(std::vector<Node*>::iterator k = StepEnrichedNodes.begin();k != StepEnrichedNodes.end();k++)
	 {
		if((*k)->isStepEnriched())
		{
		  (*k)->resolveLinearDependency(myAssociatedCrackInterior); 
		}
	 }*/

	 // Find out set of elements whose stiffness matrices need to be recomputed 
	 // Also mark new enriched nodes so that the code can update DOFs...
	 // Steps :
	 //   1. Find nodes of elements intersecting with the  crack advance
	 //   2. Find the nodal support of these newEnrichedNodes

	 // 2005-09-22
	 if(find(newSplitElems.begin(),newSplitElems.end(),newTipElem)
		== newSplitElems.end())
		newSplitElems.push_back(newTipElem);

	 std::vector<Node*> newEnrichedNodes ;
	 for(std::list<Element*> :: iterator i = newSplitElems.begin() ; i != newSplitElems.end() ; i++)
	 {
		for(size_t j = 0 ; j < (*i)->giveNumberOfNodes() ; j++)
		{
		  aNode = (*i)->giveNode(j+1) ;
		  if( find(newEnrichedNodes.begin(),newEnrichedNodes.end(),aNode)
			 == newEnrichedNodes.end())
			 newEnrichedNodes.push_back(aNode);
		}
	 }

	 map<Node*,vector<Element*> > nodeElemMap = this->domain->giveNodalSupports();
	 std::list<Element*> updatedElements ;

	 for (size_t i = 0 ; i < newEnrichedNodes.size() ; i++)
	 {
		updatedElements.insert(updatedElements.end(),
		  nodeElemMap[newEnrichedNodes[i]].begin(),nodeElemMap[newEnrichedNodes[i]].end()) ;
	 }

	 // removing the redundancies in elemsUpdatedStiffnessMatrice ...
	 updatedElements.sort();
	 updatedElements.unique(); 

	 // 3. Mark these elements so that their stiffness matrices are recomputed
	 for(std::list<Element*> :: iterator i = updatedElements.begin() ; i != updatedElements.end() ; i++)
	 {
		(*i)->setStateOfElement() ;
	 }

	 // 4. Mark new enriched nodes
	 for(size_t i = 0 ; i < newEnrichedNodes.size() ; i++)
	 {
		if(newEnrichedNodes[i]->getIsEnriched()) // after resolve linear dependency, some nodes are no longer enriched
		{
		  newEnrichedNodes[i]->setIsUpdated() ;
		}
	 }
  }
}

void CrackTip::resolveConflictsInNodalEnrichment()
// ***********************************************
// do not enrich a node with both H(x) and branch functions of the SAME CRACK
// REMARK : the code could not solve the case an element cut by one crack and contains
// a tip of other crack !!! For this, just modify the method Node::resolveConflictsInEnrichment()
{
  std::vector<Element*> *elemList = this->giveElementsInteractWithMe(); 
  assert(elemList->size() == 1);
  Element *TipElement = (*elemList)[0];

  for(size_t i = 0 ; i < TipElement->giveNumberOfNodes() ; i++)
  {
	 TipElement->giveNode(i+1)->resolveConflictsInEnrichment() ;
  }
}

void CrackTip::printYourSelf()
// ***************************
// useful for debugging
{
  std::cout << " CrackTip numbered " << this->giveNumber() << endl ;
}

double CrackTip::giveInclinationAngle()
{
  if(alpha == 1000.0)
  {
	 this->computeInclinationAngle();
  }
  return alpha ;
}

void CrackTip::computeInclinationAngle()
{
  	this->giveTipSegment();
	double x1 = tipSegment->first()->x ;
	double y1 = tipSegment->first()->y ;
	double x2 = tipSegment->second()->x ;
	double y2 = tipSegment->second()->y ;

	alpha = atan2(y2-y1,x2-x1) ; 
}