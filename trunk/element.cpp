/************************************************************************************

Copyright (C) 2005
Stephane BORDAS, Cyrille DUNANT, Vinh Phu NGUYEN, Quang Tri TRUONG, Ravindra DUDDU

This file is part of the XFEM C++ Library (OpenXFEM++) written
and maintained by above authors.

This program is free software; you can redistribute it and/or modify it.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
license text for more details.

Any feedback is welcome. Emails : nvinhphu@gmail.com, ...

*************************************************************************************/
//   file ELEMENT.CPP

#include "element.h"
#include "tri_u.h"
#include "tri6.h"
#include "quad_u.h"
#include "plateiso4.h"
#include "tetra4.h"
#include "mitc4.h"
#include "domain.h"
#include "timestep.h"
#include "timinteg.h"
#include "node.h"
#include "dof.h"
#include "material.h"
#include "bodyload.h"
#include "gausspnt.h"
#include "intarray.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "diagmtrx.h"
#include "enrichmentitem.h"
#include "cracktip.h"
#include "crackinterior.h"
#include "materialinterface.h"
#include "enrichmentfunction.h"
#include "standardquadrature.h"
#include "splitgaussquadrature.h"
#include "vertex.h"
#include "feinterpol.h"
#include "linsyst.h"
#include "functors.h"
#include "skyline.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>


#include "planelast.h"

using namespace std;

Element :: Element (int n, Domain* aDomain)
: FEMComponent (n, aDomain)
// Constructor. Creates an element with number n, belonging to aDomain.
{
  material           = 0    ;
  numberOfNodes      = 0    ;
  nodeArray          = NULL ;
  locationArray      = NULL ;
  constitutiveMatrix = NULL ;
  massMatrix         = NULL ;
  stiffnessMatrix    = NULL ;
  bodyLoadArray      = NULL ;
  gaussPointArray    = NULL ;
  neighbors          = NULL ;    // XFEM, NVP 2005
  numberOfGaussPoints      = 0 ; // XFEM, NVP 2005
  numOfGPsForJ_Integral    = 0 ; // XFEM, NVP 2005-08-13
  numberOfIntegrationRules = 0 ; // XFEM, NVP 2005
  quadratureRuleArray= NULL ;    // XFEM, NVP 2005
  standardFEInterpolation  = NULL ;  // XFEM, NVP 2005
  enrichmentFEInterpolation= NULL ;  // XFEM, NVP 2005
  enrichmentItemListOfElem = NULL ;  // XFEM, NVP 2005
  checked                  = false;  // XFEM, NVP 2005
  isUpdated = false ;                // XFEM, NVP 2005-09-03
  isMultiMaterial = false ; // one material element
  //materialIDs =
}


Element :: ~Element ()
// Destructor.
{
  delete nodeArray ;
  delete locationArray ;
  delete massMatrix ;
  delete stiffnessMatrix ;
  delete constitutiveMatrix ;
  if (gaussPointArray)
  {
	 for (size_t i = 0 ; i < numberOfGaussPoints ; i++)
		delete gaussPointArray[i] ;
	 delete gaussPointArray ;
  }

  if (quadratureRuleArray)
  {
	 for (size_t i = 0 ; i < numberOfIntegrationRules ; i++)
		delete quadratureRuleArray[i] ;
	 delete quadratureRuleArray ;
  }

  delete bodyLoadArray ;
  delete enrichmentItemListOfElem ;
  delete standardFEInterpolation ;
  delete enrichmentFEInterpolation ;
  delete neighbors ;
}

FloatMatrix*  Element :: ComputeTangentStiffnessMatrix ()
// Computes numerically the tangent stiffness matrix of the receiver, with
// the mesh subject to the total displacements D.
// Remark: geometrical nonlinarities are not acccounted for.
// MODIFIED TO TAKE ACCOUNT FOR MULTIMATERIALS ELEMENT !!!
{
  GaussPoint	    *gp;
  FloatMatrix      *b,*db,*d;
  double           dV;

  //test ConstantStiffness?
  char nlSolverClassName[32] ;
  NLSolver* nlSolver    = this->domain->giveNLSolver();
  nlSolver -> giveClassName(nlSolverClassName) ;

  if (stiffnessMatrix) {
	 delete stiffnessMatrix;
  }

  stiffnessMatrix = new FloatMatrix();
  Material *mat  = this->giveMaterial();
  for (size_t i = 0 ; i < this->giveNumberOfGaussPoints() ; i++)
  {
	 gp = gaussPointArray[i];
	 b  = this->ComputeBmatrixAt(gp);
	 // compute D matrix
	 if(! strcmp(nlSolverClassName,"ConstantStiffness"))
		d = this->giveConstitutiveMatrix()->GiveCopy();    // NO USED LONGER !!!
	 else
		d  = mat->ComputeConstitutiveMatrix(gp,this);      // USE THIS ONE
    //d->printYourself();
	 dV = this->computeVolumeAround(gp);
	 db = d->Times(b);
	 stiffnessMatrix->plusProduct(b,db,dV);
	 delete d;
	 delete db;
	 delete b;	 // SC-Purify-10.09.97
  }

  stiffnessMatrix->symmetrized();
  return stiffnessMatrix->GiveCopy(); // fix memory leaks !!!
}

FloatArray*  Element :: ComputeInternalForces (FloatArray* dElem)
// Computes the internal force vector of the receiver, with the domain sub-
// ject to the displacements D.
{
  GaussPoint		 *gp;
  FloatMatrix      *b;
  double           dV;

  Material   *mat  = this->giveMaterial();
  FloatArray *f    = new FloatArray();

  for(size_t i = 0 ; i < this->giveNumberOfGaussPoints() ; i++)
  {
	 gp    = gaussPointArray[i];
	 b     = this->ComputeBmatrixAt(gp);
	 mat -> ComputeStress(dElem,this,gp);
	 dV    = this->computeVolumeAround(gp);
	 f->plusProduct(b,gp->giveStressVector(),dV);
	 delete b;
  }
  return f;
}

void  Element :: assembleLhsAt (TimeStep* stepN)
// Assembles the left-hand side (stiffness matrix) of the receiver to
// the linear system' left-hand side, at stepN.
{
  //FloatMatrix* elemLhs ;
  //Skyline*     systLhs ;
  //IntArray*    locArray ;

  //elemLhs  = this -> ComputeLhsAt(stepN) ;
  //systLhs  = domain -> giveNLSolver() -> giveLinearSystem() -> giveLhs() ;
  //locArray = this -> giveLocationArray() ;
  //systLhs -> assemble(elemLhs,locArray) ;

  //delete elemLhs ;
}


void  Element :: assembleRhsAt (TimeStep* stepN)
// Assembles the right-hand side (load vector) of the receiver to
// the linear system' right-hand side, at stepN.
{
  //FloatArray* elemRhs ;
  //FloatArray* systRhs ;
  //IntArray*   locArray ;

  //elemRhs = this -> ComputeRhsAt(stepN) ;
  //if (elemRhs) {
  //   systRhs  = domain -> giveNLSolver() -> giveLinearSystem() -> giveRhs() ;
  //   locArray = this -> giveLocationArray() ;
  //   systRhs -> assemble(elemRhs,locArray) ;
  //   delete elemRhs ;}
}


void  Element :: assembleYourselfAt (TimeStep* stepN)
// Assembles the contributions of the receiver to the linear system, at
// time step stepN. This may, or may not, require assembling the receiver's
// left-hand side.
{
#  ifdef VERBOSE
  printf ("assembling element %d\n",number) ;
#  endif

  //CB - Modified by SC - 25.07.97
  //because we ALWAYS reform the system!
  //if (stepN -> requiresNewLhs())
  //CE - Modified by SC - 25.07.97
  //this -> assembleLhsAt(stepN) ;
  //this -> assembleRhsAt(stepN) ;
}


FloatArray*  Element :: ComputeBcLoadVectorAt (TimeStep* stepN)
// Computes the load vector due to the boundary conditions acting on the
// receiver's nodes, at stepN. Returns NULL if this array contains only
// zeroes.
// Modified by NVP 2005-09-04 for XFEM implementation.
{
  FloatArray *d, *answer ;
  FloatMatrix	*k;

  d = this -> ComputeVectorOfPrescribed('d',stepN) ;

  if(this->domain->isXFEMorFEM() == false && stepN->giveNumber() > 1)
  {
	 FloatArray *previousDPr;
	 previousDPr = this -> ComputeVectorOfPrescribed ('d',domain -> giveTimeIntegrationScheme() -> givePreviousStep()) ;
	 d->minus(previousDPr);
	 delete previousDPr;
  }
  if (d -> containsOnlyZeroes())
  {
	 answer = NULL ;
  } else
  {
	 k = this -> GiveStiffnessMatrix() ;
	 answer = k -> Times(d) -> negated() ;
	 delete k ;
  }

  delete d ;
  return answer ;
}


FloatArray*  Element :: ComputeBodyLoadVectorAt (TimeStep* stepN)
// Computes numerically the load vector of the receiver due to the body
// loads, at stepN.
{
  double      dV ;
  GaussPoint* gp ;
  FloatArray  *answer,*f,*ntf ;
  FloatMatrix *n,*nt ;

  if (this -> giveBodyLoadArray() -> isEmpty())         // no loads
	 return NULL ;

  else {
	 f = this -> ComputeResultingBodyForceAt(stepN) ;
	 if (! f)                                           // nil resultant
		return NULL ;
	 else {
		answer = new FloatArray(0) ;
		for (size_t i = 0 ; this->giveNumberOfGaussPoints() ; i++)
		{
		  gp  = gaussPointArray[i] ;
		  n   = this -> ComputeNmatrixAt(gp) ;
		  dV  = this -> computeVolumeAround(gp) ;
		  nt  = n    -> GiveTransposition() ;
		  ntf = nt   -> Times(f) -> times(dV) ;
		  answer -> add(ntf) ;
		  delete n ;
		  delete nt ;
		  delete ntf ;
		}
		delete f ;
		return answer ;
	 }
  }
}

FloatMatrix*  Element :: ComputeConsistentMassMatrix ()
// Computes numerically the consistent (full) mass matrix of the receiver.
{
  double      density,dV ;
  FloatMatrix *n,*answer ;
  GaussPoint  *gp ;

  answer  = new FloatMatrix() ;
  density = this -> giveMaterial() -> give('d') ;
  for (size_t i = 0 ; this->giveNumberOfGaussPoints() ; i++)
  {
	 gp      = gaussPointArray[i] ;
	 n       = this -> ComputeNmatrixAt(gp) ;
	 dV      = this -> computeVolumeAround(gp) ;
	 answer -> plusProduct(n,n,density*dV) ;
	 delete n ;
  }

  return  answer->symmetrized() ;
}

FloatMatrix*  Element :: computeLhsAt (TimeStep* stepN)
// Computes the contribution of the receiver to the left-hand side of the
// linear system.
{
  TimeIntegrationScheme* scheme ;

  scheme = domain -> giveTimeIntegrationScheme() ;
  if (scheme -> isStatic())
	 return  this -> ComputeStaticLhsAt (stepN) ;
  else if (scheme -> isNewmark())
	 return  this -> ComputeNewmarkLhsAt(stepN) ;
  else
  {
	 printf ("Error : unknown time integration scheme : %c\n",scheme) ;
	 exit(0) ;
	 return NULL ; //SC
  }
}


FloatArray*  Element :: ComputeLoadVectorAt (TimeStep* stepN)
// Computes the load vector of the receiver, at stepN.
{
  FloatArray* loadVector ;
  FloatArray* bodyLoadVector = NULL ;
  FloatArray* bcLoadVector   = NULL ;

  loadVector = new FloatArray(0) ;

  bodyLoadVector = this -> ComputeBodyLoadVectorAt(stepN) ;
  if (bodyLoadVector)
  {
	 loadVector -> add(bodyLoadVector) ;
	 delete bodyLoadVector ;
  }

  // BCLoad vector only at the first iteration
  if (this->domain->giveNLSolver()->giveCurrentIteration() == 1)
  {
	 bcLoadVector = this -> ComputeBcLoadVectorAt(stepN) ;
	 if (bcLoadVector)
	 {
		loadVector -> add(bcLoadVector) ;
		delete bcLoadVector ;
	 }
  }

  if (loadVector -> isNotEmpty())
	 return loadVector ;
  else
  {
	 delete loadVector ;
	 return NULL ;
  }
}

FloatMatrix*  Element :: computeMassMatrix ()
// Returns the lumped mass matrix of the receiver.
{
  FloatMatrix* consistentMatrix ;

  consistentMatrix = this -> ComputeConsistentMassMatrix() ;
  massMatrix       = consistentMatrix -> Lumped() ;
  delete consistentMatrix ;

  return massMatrix ;
}

FloatMatrix*  Element :: ComputeNewmarkLhsAt (TimeStep* stepN)
// Computes the contribution of the receiver to the left-hand side of the
// linear system, using Newmark's formula.
{
  FloatMatrix *m,*k,*lhs ;
  double      beta,dt ;

  if (stepN->giveNumber() == 0) {
	 lhs = this -> GiveStiffnessMatrix() ;
  } else {
	 beta = domain -> giveTimeIntegrationScheme() -> giveBeta() ;
	 if (beta == 0.0)
	 {
		printf ("Error: beta = 0.0 in Newmark \n") ;
		exit(0) ;
	 }
	 else
	 {
		dt  = stepN -> giveTimeIncrement() ;
		m   = this -> giveMassMatrix() -> Times(1.0 / (beta*dt*dt));
		lhs	= this -> GiveStiffnessMatrix();
		lhs->plus(m);
		delete m;
	 }
  }

  return lhs ;
}

FloatArray*  Element :: ComputeNewmarkRhsAt (TimeStep* stepN, FloatArray* dxacc)
// Computes the contribution of the receiver to the right-hand side of the
// linear system, using Newmark's formula.
{
  FloatMatrix *K;
  DiagonalMatrix *M;
  FloatArray  *fExt,*dPred,*rhs,*a,*d,*dElem,*dPrev,*temp;
  double      beta,dt ;

  fExt = this -> ComputeLoadVectorAt(stepN) ;

  if (stepN->giveNumber() == 0)
  {
	 // computes also the true stress state at the intial step, through the internal
	 // forces computation.
	 K = this -> GiveStiffnessMatrix () ;
	 d     = this -> ComputeVectorOf ('d',stepN) ;
	 dElem = dxacc->Extract(this->giveLocationArray());
	 rhs   = K -> Times(d) -> add(fExt) -> add (this->ComputeInternalForces(dElem)->negated());
	 delete d;
	 delete dElem;
	 delete K;
  }
  else
  {
	 dPred = this -> ComputeVectorOf ('D',stepN) ;
	 dPrev = this -> ComputeVectorOf ('d',domain -> giveTimeIntegrationScheme() -> givePreviousStep()) ;
	 double aNorm = dxacc->giveNorm();
	 if (aNorm == 0.)
	 { //means iteration zero (dxacc = 0)
		dElem = dPred->Minus(dPrev);
		rhs   = (this->ComputeInternalForces(dElem)->negated()) -> add(fExt);
		delete dElem;
	 } else
	 {
		M = (DiagonalMatrix*) this -> giveMassMatrix () ;
		dt  = stepN -> giveTimeIncrement() ;
		beta = domain -> giveTimeIntegrationScheme() -> giveBeta() ;
		dElem = dxacc->Extract(this->giveLocationArray());
		a = dElem->Times(1.0 / (beta*dt*dt));
		temp = dPred->Minus(dPrev);
		temp->add(dElem);
		rhs   = (M->Times(a->negated())) -> add(this->ComputeInternalForces(temp)->negated()) -> add(fExt);
		delete a ;
		delete dElem;
		delete temp;
	 }
	 delete dPred;
	 delete dPrev;
  }
  delete fExt ;
  return rhs ;
}

int  Element :: computeNumberOfDofs ()
// Returns the total number of dofs of the receiver's nodes.
{
  int n = 0 ;

  for (size_t i = 0 ; i < numberOfNodes ; i++)
	 n += this -> giveNode(i+1) -> giveNumberOfDofs() ;

  return n ;
}

size_t  Element :: computeNumberOfDisplacementDofs ()
// **************************************************
// Returns the total number of "true" dofs of the receiver's nodes.
// just read from the input file, the dofs of each node
{
  size_t n = 0 ;

  for (size_t i = 0 ; i < numberOfNodes ; i++)
  {
	 n += this -> giveNode(i+1) -> readInteger("nDofs") ;
  }

  return n ;
}


FloatArray*  Element :: ComputeResultingBodyForceAt (TimeStep* stepN)
// Computes at stepN the resulting force due to all body loads that act
// on the receiver. This force is used by the element for computing its
// body load vector.
{
  int         n ;
  BodyLoad*   load ;
  FloatArray  *force,*resultant ;

  resultant = new FloatArray(0) ;

  int nLoads = this -> giveBodyLoadArray() -> giveSize() ;

  for (size_t i = 1 ; i <= nLoads ; i++)
  {
	 n     = bodyLoadArray -> at(i) ;
	 load  = (BodyLoad*) domain->giveLoad(n) ;
	 force = load -> ComputeForceOn(this,stepN) ;
	 resultant -> add(force) ;
	 delete force ;
  }

  if (resultant->giveSize() == 0)
  {
	 delete resultant ;
	 return NULL ;
  }
  else
	 return resultant ;
}


FloatArray*  Element :: computeRhsAt (TimeStep* stepN, FloatArray* dxacc)
// Computes the contribution of the receiver to the right-hand side of the
// linear system.
{
  TimeIntegrationScheme* scheme = domain -> giveTimeIntegrationScheme() ;

  if (scheme -> isStatic())
	 return  this -> ComputeStaticRhsAt (stepN,dxacc) ;
  else if (scheme -> isNewmark())
	 return  this -> ComputeNewmarkRhsAt(stepN,dxacc) ;
  else
  {
	 printf ("Error : unknown time integration scheme : %c\n",scheme) ;
	 assert(false) ;
	 return NULL ; //SC
  }
}


FloatMatrix*  Element :: ComputeStaticLhsAt (TimeStep* stepN)
// Computes the contribution of the receiver to the left-hand side of the
// linear system, in a static analysis.

// Modified by NVP 21-10-2005 for XFEM update
// Only recompute the stiffness matrices for updated elements!!!
{
  //if (stepN->giveNumber() == 1)
  return this -> ComputeTangentStiffnessMatrix() ;
  //else
  //{
  // if(isUpdated)
  ///	return this -> ComputeTangentStiffnessMatrix() ;

  // return stiffnessMatrix->GiveCopy() ;
  //}
}


FloatArray*  Element :: ComputeStaticRhsAt (TimeStep* stepN, FloatArray* dxacc)
// Computes the contribution of the receiver to the right-hand side of the
// linear system, in a static analysis.
// Modified by NVP 2005-09-05 for XFEM ( crack growth simulation part).
{
  FloatArray *answer,*fInternal,*fExternal,*dElem,*dElemTot;

  dElem = dxacc->Extract(this->giveLocationArray());

  //add delta prescribed displacement vector - SC 29.04.99
  if(this->domain->giveNLSolver()->giveCurrentIteration() != 1)// from second iteration on
  {
	 FloatArray *currentDPr,*previousDPr;
	 currentDPr = this -> ComputeVectorOfPrescribed('d',stepN) ;
	 if( (stepN->giveNumber() > 1) && (this->domain->isXFEMorFEM() == false) )
	 {
		previousDPr = this -> ComputeVectorOfPrescribed('d',domain -> giveTimeIntegrationScheme() -> givePreviousStep()) ;
		dElemTot = (dElem->Plus(currentDPr))->Minus(previousDPr);
		delete previousDPr;
	 }
	 else
	 {
		dElemTot = dElem->Plus(currentDPr);
	 }
	 delete currentDPr;
  }
  else // first iteration
  {
	 dElemTot = dElem->Times(1.);
  }

  fInternal = this->ComputeInternalForces(dElemTot);
  fExternal = this->ComputeLoadVectorAt(stepN);

  answer = fExternal->Minus(fInternal);

  delete dElem;
  delete dElemTot;
  delete fExternal;
  delete fInternal;

  return answer;
}


FloatMatrix*  Element :: computeStiffnessMatrix ()
// Computes numerically the stiffness matrix of the receiver.
// NOT USED ANYMORE !!!
{
  double      dV ;
  FloatMatrix *b,*db,*d ;
  GaussPoint  *gp ;
  stiffnessMatrix = new FloatMatrix() ;
  for (size_t i = 0 ; this->giveNumberOfGaussPoints() ; i++)
  {
	 gp = gaussPointArray[i] ;
	 b  = this -> ComputeBmatrixAt(gp) ;
	 d  = this -> giveConstitutiveMatrix() ;
	 dV = this -> computeVolumeAround(gp) ;
	 db = d -> Times(b) ;
	 stiffnessMatrix -> plusProduct(b,db,dV) ;

	 delete b ;
	 delete db ;
  }

  return stiffnessMatrix -> symmetrized() ;
}


FloatArray*  Element :: computeStrainVector (GaussPoint* gp, TimeStep* stepN)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step stepN. The nature of these strains depends
// on the element's type.
{
  FloatMatrix *b ;
  FloatArray  *u,*Epsilon ;

  b       = this -> ComputeBmatrixAt(gp) ;
  u       = this -> ComputeVectorOf('d',stepN) ;
  Epsilon = b -> Times(u) ;

  gp -> letStrainVectorBe(Epsilon) ;    // gp stores Epsilon, not a copy
  delete b ;
  delete u ;
  return Epsilon ;
}
/*
FloatArray*  Element :: computeStressAt(GaussPoint* gp,TimeStep* stepN)
// ********************************************************************
// computes stress at Stress point gp.
{
FloatMatrix *b  = this -> ComputeBmatrixAt(gp) ;
FloatArray  *u  = this -> ComputeVectorOf('d',stepN) ;
FloatArray  *Epsilon = b -> Times(u) ;
FloatMatrix *D = this->giveConstitutiveMatrix();

FloatArray *stress = D->Times(Epsilon);

delete b ;
delete u ;
delete Epsilon ;

return stress ;
}*/

void  Element :: computeStrain (TimeStep* stepN)
// compute strain vector at stepN
{
  GaussPoint* gp ;

  for (size_t i = 1 ; i <= this->giveNumberOfGaussPoints() ; i++)
  {
	 gp = gaussPointArray[i-1] ;
	 this -> computeStrainVector(gp,stepN) ;
  }
}

FloatArray*  Element :: computeStrainIncrement (GaussPoint* gp, FloatArray* dElem)
// Returns the vector containing the strains incr. at point 'gp', with the nodal
// displacements incr. of the receiver given by dElem.
{
  FloatArray* depsilon;
  FloatMatrix* b;

  b = this -> ComputeBmatrixAt(gp);

  depsilon = b -> Times(dElem);

  delete b;

  return depsilon;
}

FloatArray* Element :: ComputeVectorOf (char u, TimeStep* stepN)
// *************************************************************
// Forms the vector containing the values of the unknown 'u' (e.g., the
// displacement) of the dofs of the receiver's nodes.

// Modified by NVP to take into account the enriched DOFs for XFEM. 2005/07/15
// u = [u1 v1 ...un vn | a1 b1 ..an bn ]
{
  Node       *nodeI ;
  int        nDofs ;

  FloatArray *answer = new FloatArray(this->computeNumberOfDofs()) ;

  int k = 0 ;

  for (size_t i = 0 ; i < numberOfNodes ; i++)
  {
	 nodeI = this->giveNode(i+1) ;
	 nDofs = nodeI->giveNumberOfTrueDofs () ;
	 for (size_t j = 1 ; j <= nDofs ; j++)
		answer->at(++k) = nodeI->giveDof(j)->giveUnknown(u,stepN) ;
  }

  if(this->isEnriched() == false) // non-enriched element
	 return answer ;

  // enriched element
  size_t numOfTrueDofs ;
  for (size_t i = 0 ; i < numberOfNodes ; i++)
  {
	 nodeI = this->giveNode(i+1) ;
	 numOfTrueDofs = nodeI->giveNumberOfTrueDofs() ;
	 nDofs = nodeI->giveNumberOfDofs() - numOfTrueDofs ;
	 for (size_t j = 1 ; j <= nDofs ; j++)
		answer->at(++k) = nodeI->giveDof(j+numOfTrueDofs)->giveUnknown(u,stepN) ;
  }

  return answer ;
}

FloatArray* Element :: ComputeVectorOfDisplacement (char u, TimeStep* stepN)
// *************************************************************************
// computes the "true" displacement vector of the receiver
// XFEM implementation.
{
  Node       *nodeI ;
  size_t      nDofs ;

  FloatArray *answer = new FloatArray(this->computeNumberOfDisplacementDofs()) ;
  int k      = 0 ;

  for (size_t i = 0 ; i < numberOfNodes ; i++) {
	 nodeI = this->giveNode(i+1) ;
	 nDofs = nodeI->giveNumberOfTrueDofs() ;
	 for (size_t j = 1 ; j <= nDofs ; j++)
		answer->at(++k) = nodeI->giveDof(j)->giveUnknown(u,stepN) ;
  }

  return answer ;
}

FloatArray* Element :: ComputeVectorOfPrescribed (char u, TimeStep* stepN)
// Forms the vector containing the prescribed values of the unknown 'u'
// (e.g., the prescribed displacement) of the dofs of the receiver's
// nodes. Puts 0 at each free dof.
{
  Node       *nodeI ;
  Dof        *dofJ ;
  FloatArray *answer ;
  int        nDofs ;

  answer = new FloatArray(this->computeNumberOfDofs()) ;
  int k  = 0 ;
  for (size_t i = 0 ; i < numberOfNodes ; i++)
  {
	 nodeI = this->giveNode(i+1) ;
	 nDofs = nodeI->giveNumberOfDofs() ;
	 for (size_t j = 1 ; j <= nDofs ; j++)
	 {
		dofJ = nodeI->giveDof(j) ;
		if (dofJ -> hasBc())
		  answer->at(++k) = dofJ->giveUnknown(u,stepN) ;
		else
		  answer->at(++k) = 0. ;
	 }
  }

  return answer ;
}


IntArray*  Element :: giveBodyLoadArray ()
// Returns the array which contains the number of every body load that act
// on the receiver.
{
  int numberOfLoads ;

  if (! bodyLoadArray)
  {
	 numberOfLoads = this -> readIfHas("bodyLoads") ;
	 bodyLoadArray = new IntArray(numberOfLoads) ;
	 for (size_t i = 1 ; i <= numberOfLoads ; i++)
		bodyLoadArray->at(i) = this->readInteger("bodyLoads",i+1) ;
  }

  return bodyLoadArray ;
}


FloatMatrix*  Element :: giveConstitutiveMatrix ()
// Returns the elasticity matrix {E} of the receiver.
{
  if (! constitutiveMatrix)
	 this -> computeConstitutiveMatrix() ;

  return constitutiveMatrix ;
}

IntArray*  Element :: giveLocationArray ()
// Returns the location array of the receiver.
// Modified by NVP to take into account the presence of enriched Dofs
{
  if (! locationArray)
  {
	 this->computeLocationArray();
  }
  return locationArray ;
}

IntArray*  Element :: computeLocationArray ()
// ******************************************
// Returns the location array of the receiver.
// Modified by NVP to take into account the presence of enriched Dofs
{
  IntArray* nodalStandardArray ;      // standard scatter vector of node
  IntArray* nodalEnrichedArray ;      // enriched scatter vector of node

  locationArray = new IntArray(0) ;   // total scatter vector of element
  for(size_t i = 0 ; i < numberOfNodes ; i++)
  {
	 nodalStandardArray = this->giveNode(i+1)->giveStandardLocationArray() ;
	 locationArray = locationArray -> followedBy(nodalStandardArray) ;
  }
  // non-enriched elements
  if(this->isEnriched() == false)
  {
	 /* // DEBUG ...
	 std::cout << "Dealing with element " << this->giveNumber() << std::endl ;
	 for(size_t i = 0 ; i < locationArray->giveSize() ; i++)
	 std::cout << (*locationArray)[i] << " " ;
	 std::cout << std::endl ; */
	 return locationArray ;
  }

  // enriched elements
  for (size_t i = 0 ; i < numberOfNodes ; i++)
  {
	 nodalEnrichedArray = this -> giveNode(i+1) -> giveEnrichedLocationArray() ;
	 if(nodalEnrichedArray)   // enriched node
		locationArray = locationArray -> followedBy(nodalEnrichedArray) ;
  }

  /*// -----------------    DEBUG ONLY 2005-09-29 -------------------------
  if(this->isEnriched())
  {
  std::cout << " Location array of element " << this->giveNumber() << " is :" << endl ;
  for(size_t i = 0 ; i < locationArray->giveSize() ; i++)
  std::cout << (*locationArray)[i] << " " ;
  std::cout << std::endl ;
  }
  // --------------------------------------------------------------------*/

  return locationArray ;
}

GaussPoint**  Element :: giveGaussPointArray()
// *******************************************
// 2005-09-03 : modify for crack growth problem
// Some elements need changed Gauss Quadrature.
{
  if (gaussPointArray == NULL)
	 this->computeGaussPoints();
  else if(isUpdated)
	 this->computeGaussPoints();

  return gaussPointArray;
}


FloatMatrix*  Element :: giveMassMatrix ()
// Returns the mass matrix of the receiver.
{
  if (! massMatrix)
	 this -> computeMassMatrix() ;
  return massMatrix ;
}

Material*  Element :: giveMaterial ()
// **********************************
// Returns the material of the receiver.
{
  if (! material)
  {
	 material = this -> readInteger("mat") ;
	 //std::cout <<this->giveNumber() << " CUC CUT !!! " ;
  }

  return  domain -> giveMaterial(material) ;
}



Node*  Element :: giveNode (int i)
// Returns the i-th node of the receiver.
{
  int  n ;

  if (! nodeArray)
	 nodeArray = new IntArray(numberOfNodes) ;

  n = nodeArray->at(i) ;
  if (! n) {
	 n = this -> readInteger("nodes",i) ;
	 nodeArray->at(i) = n ;}

  return  domain -> giveNode(n) ;
}

IntArray* Element ::giveNodeArray()
{
  if (! nodeArray)
  {
	 nodeArray = new IntArray(numberOfNodes) ;
	 for (size_t i = 0 ; i < numberOfNodes ; i++)
		(*nodeArray)[i] = this->readInteger("nodes",i+1) ;
  }

  return nodeArray;
}

size_t Element :: giveNumberOfGaussPoints()
// ****************************************
{
  if(numberOfGaussPoints == 0)
	 this->computeGaussPoints();

  return numberOfGaussPoints ;
}

size_t Element :: giveNumOfGPtsForJ_Integral()
// *******************************************
{
  if(numOfGPsForJ_Integral == 0)
	 this->setGaussQuadForJ_Integral();

  return numOfGPsForJ_Integral ;
}

FloatMatrix*  Element :: GiveStiffnessMatrix ()
// ********************************************
// Returns the stiffness matrix of the receiver.
// Modification made by NVP for multistep problems (XFEM).
// 2005-09-03
{
  if (! stiffnessMatrix)
	 return this -> ComputeTangentStiffnessMatrix() ;
  //else if(isUpdated)
  // return this -> ComputeTangentStiffnessMatrix() ;

  return stiffnessMatrix->GiveCopy() ;
}


void  Element :: instanciateYourself ()
// Gets from input file all data of the receiver.
{
  int i ;

#  ifdef VERBOSE
  printf ("instanciating element %d\n",number) ;
#  endif

  material = this -> readInteger("mat") ;
  nodeArray = new IntArray(numberOfNodes) ;
  for (i=1 ; i<=numberOfNodes ; i++)
	 nodeArray->at(i) = this->readInteger("nodes",i) ;
  this -> giveBodyLoadArray() ;
}


void  Element :: printOutputAt (TimeStep* stepN, FILE* strFile, FILE* s01File)
// Performs end-of-step operations.
{
  GaussPoint* gp ;

#  ifdef VERBOSE
  printf ("element %d printing output\n",number) ;
#  endif

  fprintf (strFile,"element %d :\n",number) ;

  for (size_t i  = 0 ; i < this->giveNumberOfGaussPoints() ; i++) {
	 gp = gaussPointArray[i] ;
	 this -> computeStrainVector(gp,stepN) ;
	 //no computation of the stress vector here, it
	 //has already been done during the calculation
	 //of Finternal!!!
	 gp   -> computeStressLevel() ;
	 gp   -> printOutput(strFile) ;
	 gp   -> printBinaryResults(s01File) ;
  }
}


Element*  Element :: typed ()
// Returns a new element, which has the same number than the receiver,
// but is typed (PlaneProblem, or Truss2D,..).
{
  Element* newElement ;
  char     type[32] ;

  this -> readString("class",type) ;
  newElement = this -> ofType(type) ;

  return newElement ;
}


void  Element :: updateYourself()
// ******************************
// Updates the receiver at end of step.
// Modified by NVP for XFEM implementation. 2005-09-05
{
#  ifdef VERBOSE
  printf ("updating element %d\n",number) ;
#  endif

  if(this->domain->isXFEMorFEM() == false)
  {
	 for (size_t i = 0 ; i < numberOfGaussPoints ; i++)
		gaussPointArray[i] -> updateYourself() ;
  }
  else
  {
	 for (size_t i = 0 ; i < numberOfGaussPoints ; i++)
		gaussPointArray[i] -> updateYourselfForXFEM() ;
  }

  delete locationArray ;
  locationArray = NULL ;
}



FloatMatrix*  Element ::ComputeBmatrixAt(GaussPoint *aGausspoint)
// **************************************************************
// Computes the general B matrix  of the receiver, B = [Bu Ba]
// Including non enriched part and enriched parts if element is enriched
// B = [Bu Ba] with
// Bu = [N1,x  0    N2,x  0    N3,x  0
//       0     N1,y 0     N2,y 0     N3,y
//       N1,y  N1,x N2,y  N2,x N3,y  N3,x]
//
// Ba = [(Nbar1(phi-phiAtNode1)),x 0 (Nbar2(phi-phiAtNode2)),x 0 (Nbar3(phi-phiAtNode3)),x 0
//        0(Nbar1(phi-phiAtNode1)),y 0 (Nbar2(phi-phiAtNode2)),y 0 (Nbar3(phi-phiAtNode3)),y
//       (Nbar1(phi-phiAtNode1)),y (Nbar1(phi-phiAtNode1)),x ...]
{
  // computes the standard part of B matrix : Bu
  FloatMatrix *Bu = this->ComputeBuMatrixAt(aGausspoint);

  // non enriched elements
  if (this->isEnriched() == false)
	 return Bu;

  // enriched elements,then computes the enriched part Ba
  FloatMatrix *Ba = new FloatMatrix();
  vector<EnrichmentFunction*> *enrFnVector; // vector of enrichment functions
  FloatArray  *gradPhiGP ;                  // grad of enr. func. at Gauss points
  FloatMatrix *temp ;
  double      N,dNdx,dNdy,phiGP,phiNode,dPhidXGP,dPhidYGP ;

  Mu::Point *Coord = aGausspoint -> giveCoordinates() ; // local coord. of Gauss point

  // Get the shape functions multiplied with the enr. functions...
  FloatArray  *Nbar = this->giveXFEInterpolation()->evalN(Coord);
  FloatMatrix *gradNbar = this->giveXFEInterpolation()->evaldNdx(domain,this->giveNodeArray(),Coord);

  for(size_t i = 0 ; i < numberOfNodes ; i++)
  {
	 Node* nodeI = this->giveNode(i+1) ;
	 if (nodeI->getIsEnriched())          // this node is enriched, then continue ...
	 {
		N = (*Nbar)[i];                    // shape function Ni
		dNdx = gradNbar->at(i+1,1) ;       // derivative of Ni w.r.t  x coord.
		dNdy = gradNbar->at(i+1,2) ;       // derivative of Ni w.r.t  y coord.

		// loop on enrichment items of nodeI
		list<EnrichmentItem*> *enrItemList = nodeI->giveEnrItemListOfNode();
		for (list<EnrichmentItem*>::iterator iter = enrItemList->begin(); iter != enrItemList->end(); ++iter)
		{
		  // get the enrichment funcs of current enr. item
		  enrFnVector = (*iter)->giveEnrFuncVector();

		  //loop on vector of enrichment functions ...
		  for(size_t k = 0 ; k < enrFnVector->size() ; k++ )
		  {
			 EnrichmentFunction* enrFn = (*enrFnVector)[k];

			 // let Enr. function,enrFn, know for which enr. item it is modeling
			 enrFn->findActiveEnrichmentItem(*iter);

			 // value of enrichment function at gauss point
			 phiGP = enrFn->EvaluateYourSelfAt(aGausspoint);

			 // value of enrichment function at node
			 phiNode = enrFn->EvaluateYourSelfAt(this,nodeI);

			 // grad of enrichment function at Gauss point
			 gradPhiGP = enrFn->EvaluateYourGradAt(aGausspoint);
			 dPhidXGP = (*gradPhiGP)[0] ; // derivative of Phi w.r.t  x coord.
			 dPhidYGP = (*gradPhiGP)[1] ; // derivative of Phi w.r.t  y coord.

			 double a = dNdx * (phiGP - phiNode) + N * dPhidXGP ;
			 double b = dNdy * (phiGP - phiNode) + N * dPhidYGP ;

			 temp = new FloatMatrix(4,2);

			 temp->at(1,1) = a   ; temp->at(1,2) = 0.0 ;
			 temp->at(2,1) = 0.0 ; temp->at(2,2) = b   ;
			 temp->at(3,1) = b   ; temp->at(3,2) = a   ;

			 Ba = Ba->FollowedBy(temp);
			 delete temp ;              // Purify, 11-10-05
			 delete gradPhiGP ;         // Purify, 11-10-05

		  }  // end of loop on enr. functions
		}    // end of loop on enr. items
	 }
  }        // end of loop on element nodes

  FloatMatrix *answer = Bu->FollowedBy(Ba);

  delete Ba ;                     // Purify, 11-10-05
  delete Nbar ; delete gradNbar ; // Purify, 11-10-05

  return answer;
}

bool Element:: isEnriched()
// ************************
//  Returns true if element is enriched (at least one node is enriched)
//  and false otherwise
{
  size_t count = 0 ;

  for(size_t i = 0 ; i < numberOfNodes ; i++)
  {
	 if(this->giveNode(i+1)->getIsEnriched())
	 {
		count += 1 ;
		i = numberOfNodes ;
	 }
  }

  return (count != 0)? true:false ;
}

Element*  Element :: ofType (char* aClass)
// ***************************************
// Returns a new element, which has the same number than the receiver,
// but belongs to aClass (Tri3_U, Tetra4, ...).
{
  Element* newElement ;

  if (! strcmp(aClass,"Q4U"))
	 newElement = new Quad4_U(number,domain);
  else if (! strcmp(aClass,"T3U"))
	 newElement = new Tri3_U(number,domain) ;
  else if (! strcmp(aClass,"T6U"))
	 newElement = new Tri6_U(number,domain) ;
  else if (! strcmp(aClass,"PQ4"))
	 newElement = new PlateIsoQ4(number,domain) ;
  else if (! strcmp(aClass,"MITC4"))
	 newElement = new MITC4(number,domain) ;
  else if (! strcmp(aClass,"H4U"))
	 newElement = new Tetra4(number,domain) ;
  else
  {
	 printf ("%s : unknown element type \n",aClass) ;
	 assert(false) ;
  }

  return newElement ;
}

void Element::treatGeoMeshInteraction()
// ************************************
// Check if element interacts with enrichment item or not. If so, insert this element
// into the list of each enrichment item
// Modified at 2005-09-07 to make it more efficient than before.
// 28-12-2005: MATERIAL INTERFACE IMPLEMENTATION, ASSUMING 2 MATERIALS
{
  EnrichmentItem *enrItem;

  for(size_t i = 0 ; i < domain->giveNumberOfEnrichmentItems() ; i++)
  {
	 enrItem = domain->giveEnrichmentItem(i+1);
	 enrItem->treatMeshGeoInteraction(this);
  }
}

void Element :: isEnrichedWith(EnrichmentItem* enrItem)
// *****************************************************
// If element is enriched with enrichment item enrItem, then insert
// enrItem into the list of enrichment items
{
  if(enrichmentItemListOfElem == NULL)
	 enrichmentItemListOfElem = new std::list<EnrichmentItem*>;

  if( find(enrichmentItemListOfElem->begin(),enrichmentItemListOfElem->end(),enrItem)
	 == enrichmentItemListOfElem->end())
	 enrichmentItemListOfElem->push_back(enrItem);
}

std::list<Element*>*  Element :: ComputeNeighboringElements()
// **********************************************************
// Loop on Element's nodes and get the nodal support of these nodes
// and insert into list<Element*> neighbors
// Since there are redundancies, first sort this list and then list.unique()
// Criterion for sort is defined by operator < ( compare the element.number)
{
  map<Node*,vector<Element*> > nodeElemMap = this->domain->giveNodalSupports();
  neighbors = new std::list<Element*> ;
  for (size_t i = 0 ; i < numberOfNodes ; i++)
  {
	 Node *aNode = this->giveNode(i+1);
	 neighbors->insert(neighbors->end(),nodeElemMap[aNode].begin(),nodeElemMap[aNode].end()) ;
  }
  // removing the redundancies in neighbors ...
  neighbors->sort();
  neighbors->unique();
  neighbors->remove(this); // not contain the receiver !!!

  return neighbors ;
}

std::list<Element*>*  Element :: giveNeighboringElements()
// *******************************************************
// returns the neighboring elements of the receiver, if it does not exist yet,
// compute it.
{
  if (neighbors == NULL)
	 neighbors = this->ComputeNeighboringElements();
  return neighbors ;
}

void Element :: printMyNeighbors()
// *******************************
// Print neighbors of the receiver.
// Useful for debugging
{
  neighbors = this->giveNeighboringElements();

  std::cout << " Neighbors of element " << number << " : ";
  for (list<Element*>::iterator it = neighbors->begin(); it != neighbors->end(); ++it)
	 std::cout << (*it)->giveNumber() << " " ;
  std::cout << std::endl;
}

Mu::Point* Element :: giveMyCenter()
// *********************************
// compute the gravity center of the receiver
{
  double sumX = 0.0 ; double sumY = 0.0 ;
  for(size_t i = 0 ; i < numberOfNodes ; i++)
  {
	 Node *nodeI = this -> giveNode(i+1);
	 double xI = nodeI->giveCoordinate(1);
	 double yI = nodeI->giveCoordinate(2);
	 sumX += xI ;
	 sumY += yI ;
  }

  double xc = sumX/numberOfNodes ;
  double yc = sumY/numberOfNodes ;

  return new Mu::Point(xc,yc);
}

bool Element :: in(Mu::Circle *c)
// ******************************
// if at least one node of the receiver belong to the circle c, then
// this element is considered locate inside c.
{
  size_t count = 0 ;
  Mu::Point *p;
  Node *aNode;

  for(size_t i = 0 ; i < numberOfNodes ; i++)
  {
	 aNode = this->giveNode(i+1) ;
	 p = aNode->makePoint();
	 if(c->in(p))
	 {
		count += 1 ;
		delete p;
		i = numberOfNodes ;
	 }
  }
  return (count != 0)? true:false ;
}

std::list<Element*> Element :: conflicts(Mu::Circle *c)
// ****************************************************
// With the help of Cyrille Dunant.
{
  checked = true ;
  std::list<Element*> ret,temp ;
  ret.push_back(this);

  std::list<Element*>* myNeighbors = this->giveNeighboringElements();
  for (std::list<Element*>::iterator it = myNeighbors->begin(); it != myNeighbors->end(); it++)
  {
	 if( ((*it)->checked == false) && ((*it)->in(c)))
	 {
		temp = (*it)->conflicts(c);
		ret.insert(ret.end(),temp.begin(),temp.end());
	 }
  }

  return ret ;
}

bool Element :: intersects(Mu::Circle *c)
// **************************************
// Check the intersection between element and circle c
// used for determining the annular domain for J integral computation.
// 2005-08-10.
{
  std::vector<double> distanceVect ;

  double radius  = c->getRadius() ;
  double centerX = c->getCenter()->x ;
  double centerY = c->getCenter()->y ;

  for(size_t i = 0 ; i < numberOfNodes ; i++)
  {
	 double x = this->giveNode(i+1)->giveCoordinate(1);
	 double y = this->giveNode(i+1)->giveCoordinate(2);

	 double r = sqrt((x-centerX)*(x-centerX)+(y-centerY)*(y-centerY));

	 distanceVect.push_back(r-radius);
  }
  double max = (*std::max_element(distanceVect.begin(),distanceVect.end()));
  double min = (*std::min_element(distanceVect.begin(),distanceVect.end()));

  return ( max*min <= 0 ? true : false ) ;
}

bool Element::intersects(Mu::Segment *seg)
// ***************************************
{
  std::vector<Mu::Point> pts ;
  for(size_t i = 0 ; i < numberOfNodes ; i++)
  {
	 Mu::Point *p = this->giveNode(i+1)->makePoint() ;
	 pts.push_back(*p);
	 delete p ;
  }

  bool ret = false ;

  for(size_t i = 0 ; i < pts.size() ;  i++)
  {
	 Mu::Segment s(pts[i],pts[(i+1)%pts.size()]) ;
	 ret = ret || seg->intersects(&s) ;
  }

  return ret ;
}

bool Element :: IsOnEdge(Mu::Point* testPoint)
// *******************************************
// This method allows the tip touch element edge or element node
// 2005-09-06


{
	size_t i;
	std::vector<Mu::Point> pts ;
	for(i = 0 ; i < numberOfNodes ; i++)
	{
		Mu::Point *p = this->giveNode(i+1)->makePoint() ;
		pts.push_back(*p);
		delete p ;
	}

	bool found = false ;
	i = 0 ;
	while( (i < pts.size()) && (!found) )
	{
		Mu::Segment s(pts[i],pts[(i+1)%pts.size()]) ;
		if(s.on(testPoint))
		{
			std:: cout << " Ah, I found you ! " << std::endl ;
			found = true ;
		}
		i++ ;
  }

  return found ;
}

bool Element :: isWithinMe(Mu::Point *p)
// *************************************
// check if p is within the receiver using orientation test.
// Used in method PartitionMyself() to detect kink points
// 2005-09-06
{
  double const EPSILON = 0.00000001;

  double x0 = p->x ;
  double y0 = p->y ;

  double delta,x1,y1,x2,y2 ;
  size_t count = 0;

  for (size_t i = 1 ; i <= numberOfNodes ; i++)
  {
	 // coordinates of first node
	 x1 = this->giveNode(i)->giveCoordinate(1);
	 y1 = this->giveNode(i)->giveCoordinate(2);

	 // coordinates of second node
	 if(i != numberOfNodes)
	 {
		x2 = this->giveNode(i+1)->giveCoordinate(1);
		y2 = this->giveNode(i+1)->giveCoordinate(2);
	 }
	 else
	 {
		x2 = this->giveNode(1)->giveCoordinate(1);
		y2 = this->giveNode(1)->giveCoordinate(2);
	 }

	 delta = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0) ;

	 if (delta > EPSILON)
		count += 1 ;
  }

  return (count == numberOfNodes);
}

bool Element :: isCoincidentToMyNode(Mu::Point *p)
// ***********************************************
{
  for(size_t i = 0 ; i < numberOfNodes ; i++)
  {
	 Mu::Point *pp = this->giveNode(i+1)->makePoint() ;
	 if( pp == p )
		return true ;
	 delete pp ;
  }

  return false ;
}

void Element :: clearChecked()
// ***************************
// set checked false for the next time checking
{
  checked = false ;
}

void Element :: printMyEnrItems()
// ******************************
// debug only
{
  if(enrichmentItemListOfElem != NULL)
  {
	 std::cout << " Element " << number << " interacted with  " ;
	 for (list<EnrichmentItem*>::iterator it = enrichmentItemListOfElem->begin(); it != enrichmentItemListOfElem->end(); ++it)
	 {
		(*it)->printYourSelf();
	 }
	 std::cout << std::endl ;
  }
}

/*! RB-SB-2004-10-29
takes a string and appends it with the values of stress
for the receiver. Goes to the next line to allow further
information to be stored in the string. */
void Element :: exportStressResultsToMatlab(string& theStringStress)
// *****************************************************************
{
  double val1 = 0.0; double val2 = 0.0; double val3 = 0.0; double val4 = 0.0;
  double val5 = 0.0; double vonMises ;
  FloatArray* stressveci;

  for(size_t i = 1 ; i <= this->giveNumberOfGaussPoints(); i++)
  {
	 stressveci = this->giveGaussPointNumber(i)->giveStressVector();
	 vonMises = this->giveGaussPointNumber(i)->computeVonMisesStress();//NVP 2005
	 assert(stressveci->giveSize() == 4);
	 val1 = stressveci->at(1);
	 val2 = stressveci->at(2);
	 val3 = stressveci->at(3);
	 val4 = stressveci->at(4);
	 val5 = vonMises ;            //NVP 2005

	 char val1c[50];
	 char val2c[50];
	 char val3c[50];
	 char val4c[50];
	 char val5c[50];              //NVP 2005

	 _gcvt( val1, 17, val1c );
	 _gcvt( val2, 17, val2c );
	 _gcvt( val3, 17, val3c );
	 _gcvt( val4, 17, val4c );
	 _gcvt( val5, 17, val5c );    //NVP 2005

	 string space(" ");
	 string newline("\n");
	 string valString;

	 valString += val1c;
	 valString += space;
	 valString += val2c;
	 valString += space;
	 valString += val3c;
	 valString += space;
	 valString += val4c;
	 valString += space;
	 valString += val5c;
	 valString += newline;
	 theStringStress += valString;
  }
}

/*! RB-SB-2004-10-29
takes a string and appends it with the values of strain
for all of the receiver's Gauss points. Goes to the next line to allow further
information to be stored in the string. Results are written as
for GAUSS POINT 1 : on row 1 : sigmaxx sigmayy sigmaxy sigmazz
...
for GAUSS POINT N : on row N : sigmaxx sigmayy sigmaxy sigmazz
N is the number of Gauss points of the receiver.
*/
void Element :: exportStrainResultsToMatlab(string& theStringStrain)
// *****************************************************************
{
  double val1 = 0.0; double val2 = 0.0; double val3 = 0.0; double val4 = 0.0;//values
  FloatArray* strainveci;//strain vector
  int ngp = this->giveNumberOfGaussPoints();//number of Gauss points
  GaussPoint* gp;//Gauss point
  TimeStep* stepN = this->domain->giveTimeIntegrationScheme()->giveCurrentStep();//current step

  for (size_t i = 1 ; i <= ngp ; i++)
  {
	 gp = this->giveGaussPointNumber(i);  //get Gauss point number i
	 this->computeStrainVector(gp,stepN); //compute the strain vector to make sure its value exists
	 strainveci = this->giveGaussPointNumber(i)->giveStrainVector();//get the strain vector
	 assert(strainveci->giveSize() == 4); //make sure it's ok
	 val1 = strainveci->at(1);            //get the values of the strains
	 val2 = strainveci->at(2);
	 val3 = strainveci->at(3);
	 val4 = strainveci->at(4);
	 char val1c[50];                     //for transformation into chars and strings
	 char val2c[50];
	 char val3c[50];
	 char val4c[50];

	 _gcvt( val1, 17, val1c );//transform the values val1... into chars val1c taking 14 digits
	 _gcvt( val2, 17, val2c );
	 _gcvt( val3, 17, val3c );
	 _gcvt( val4, 17, val4c );

	 string space(" ");//define some useful strings for output
	 string newline("\n");
	 string valString;
	 valString += val1c;//concatenate the strings together
	 valString +=space;
	 valString += val2c;
	 valString +=space;
	 valString += val3c;
	 valString +=space;
	 valString += val4c;
	 valString +=newline;
	 //final string composed of the four values of the strains at the Gauss point
	 theStringStrain+=valString;//the final, modified, string used for output
  }
}

/*! RB-SB-2004-10-29
Returns the global coordinates of the Gauss points of
the receiver.
*/
void Element :: exportGaussPointsToMatlab(string& theString)
{

  GaussPoint* gp;  //Gauss point
  FloatMatrix* N;  //array of shape functions
  FloatArray* XI;  //array of nodal coordinates
  FloatArray* globalCoords; //array of global coords.
  //double  weight = 0. ; // total weight of Gauss points, for debug only. NVP 2005-07-28

  char value[30];
  string space(" ");
  string newline("\n");
  for (size_t k = 1 ; k <= this->giveNumberOfGaussPoints() ; k++)
  {
	 gp = this->giveGaussPointNumber(k);
	 //weight += gp->giveWeight();
	 N  = this->ComputeNmatrixAt(gp);
	 XI = this->ComputeGlobalNodalCoordinates();
	 globalCoords = N->Times(XI);
	 _gcvt(globalCoords->at(1),17,value);//transforms the double param1 in char param3 with 14 digits
	 theString+=value;
	 theString+=space;
	 _gcvt(globalCoords->at(2),17,value);
	 theString+=value;
	 theString+=space;
	 theString+=newline;
	 delete N;
	 delete XI;
	 delete globalCoords;
  }

  // _gcvt(weight,3,value);
  //theString+=value;
  //theString+=space;

    //Changed by M. Forton 5/11/11 - Uncomment to fix
    //theString += newline;
}

/*! 		  RB-SB-2004-10-29
Returns the array of nodal coordinates for an element.
The values are stored as follows:
[xNode1
yNode1
xNode2
yNode2
...
xNode_numnode
yNode_numnode]
*/
FloatArray* Element :: ComputeGlobalNodalCoordinates()
{
  size_t n = this->giveNumberOfNodes();
  FloatArray* result = new FloatArray(2*n);
  Node* node;

  for (int k = 1 ; k <= n ; k++)
  {
	 node = this->giveNode(k);
	 result->at(2*k-1) = node->giveCoordinate(1);
	 result->at(2*k)   = node->giveCoordinate(2);
  }

  return result;
}

GaussPoint* Element :: giveGaussPointNumber(int i)
// ************************************************
{
  if (gaussPointArray) {
	 return gaussPointArray[i-1];
  }
  else {
	 printf("Sorry, cannot give Gauss point number %d of element %d \n",i,number);
	 printf("The Gauss Point Array for element %d was never created \n",number);
	 return NULL;
	 exit(0);
  }
}
/*
void Element :: exportStressPointsToMatlab(string& theStringStress,TimeStep* stepN)
// ********************************************************************************
// NVP - 2005-07-19
{
double val1 = 0.0; double val2 = 0.0; double val3 = 0.0; double val4 = 0.0;
FloatArray *stressveci;

this->computeStressPoints();

for (size_t i = 0 ; i < this->giveNumberOfStressPoints(); i++)
{
stressveci = this->computeStressAt(gpForStressPlot[i],stepN);

double II = stressveci->computeInvariantJ2();
double equiStress = sqrt(3.0 * II) ; // equivalent von Mises stress.

val1 = stressveci->at(1);
val2 = stressveci->at(2);
val3 = stressveci->at(3);
val4 = stressveci->at(4);
double val5 = equiStress ;

char val1c[50];
char val2c[50];
char val3c[50];
char val4c[50];
char val5c[50];

_gcvt( val1, 17, val1c );
_gcvt( val2, 17, val2c );
_gcvt( val3, 17, val3c );
_gcvt( val4, 17, val4c );
_gcvt( val5, 17, val5c );

string space(" ");
string valString;

valString += val1c;
valString += space;
valString += val2c;
valString += space;
valString += val3c;
valString += space;
valString += val4c;
valString += space;
valString += val5c;
valString += space;

theStringStress += valString;
}
string newline("\n");
theStringStress += newline ;
}*/

void	Element::reinitializeStiffnessMatrix()
{
  delete stiffnessMatrix;
  stiffnessMatrix = NULL;
}

void Element::computeNodalLevelSets(EnrichmentItem* enrItem)
{
  GeometryEntity *geo = enrItem->giveMyGeo();
  for(size_t i = 0 ; i < numberOfNodes ; i++)
  {
	 Node *nodeI = this->giveNode(i+1);
	 Mu::Point *p = nodeI->makePoint();
	 double ls = geo->computeSignedDistanceOfPoint(p);
	 nodeI->setLevelSets(enrItem,ls);
  }
}

void Element::updateMaterialID()
{
  if(material == 0)
	 material = 2;
  else
	 material += 1 ;
}

void Element :: resolveConflictsInEnrItems()
// ****************************************
// check if this  list contains both a CrackTip and a CrackInterior
// then remove the CrackInterior from the list.
// functor IsType<CrackTip,EnrichmentItem>() defined generically in file "functors.h"
{
  list<EnrichmentItem*> ::iterator iter1,iter2;

  iter1=find_if(enrichmentItemListOfElem->begin(),enrichmentItemListOfElem->end(),IsType<MaterialInterface,EnrichmentItem>());
  iter2=find_if(enrichmentItemListOfElem->begin(),enrichmentItemListOfElem->end(),IsType<CrackInterior,EnrichmentItem>());

  if (iter1 != enrichmentItemListOfElem->end() && iter2 != enrichmentItemListOfElem->end())
	 enrichmentItemListOfElem->remove(*iter2);
}

