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


#ifndef _ANSPLATE_H_
#define _ANSPLATE_H_

#include "plateiso.h"
#include "enrichmentitem.h"
#include <vector>



//! An assumed natural strain elements for bending plate
/*!
* This class is the super class of an ANS elements. Avoid shear locking and distorted meshes
* are the major advantages of this kind of element.
* Each node has three DOFs: tranverse displacement w and two rotations theta_x and theta_y
* The shear strain and bending strain are interpolated separately.
* The more well known name is the MITC-n elements, where MITC stands for 
* Mixed Interpolation of Tensorial Components. Or Bathe-Divorkin elements.
*
* Inputs: 
*  - number of nodes and normal isoparametric shape functions as usual
*  - number and location of sampling points
*  - shape functions used for interpolating the natural shear strain 

* Remarks: for enriched elements (cracked plate), attention with the integration scheme!!!
* Start implemented: 1-12-2005 by Nguyen Vinh Phu, some portions are taken from work
* of Mourad Belgasmia.
*/

class ANSPlate : public PlateIso
{
public :

  ANSPlate (int,Domain*) ;              //!< constructor
  ~ANSPlate ()  {;}                     //!< destructor

  /** @name  Numerical integration of stiffness matrix
  *  Functions to perform the numerical integration of stiffness matrix
  */
  //@{

  FloatMatrix*       ComputeBSmatrixAt (GaussPoint*) ; // shearing stiffness matrix
  /*!
   * Assumed strain shape functions, used to interpolate the natural shear strain
   * Two rows corresponding to two shear strains along Ksi and Eta directions
	*/
  virtual FloatMatrix*  ComputeAssumedStrainShapeFunction(GaussPoint*) = 0 ;  
  /*!
   * Discretized natural shear strains matrix Bns
	*/
  FloatMatrix*       ComputeBnsMatrixAt(Mu::Point*) ;  
  //@}


  /** @name  Stiffness matrix
  *  Functions to compute the total stiffness matrix from its contributions 
  */
  //@{
  FloatMatrix*       computeShearingStiffnessMatrix();    // [Ks]
  //@}

  /*!
  Setting the sampling points used for interpolation the natural shear strain
  */
  virtual void       setttingSamplingPoints() = 0 ;

private:
  std::vector<Mu::Point*> samplingPoints; //<! points whose shear strain values used to
  // to interpolate
} ;

#endif // _ANSPLATE_H_
