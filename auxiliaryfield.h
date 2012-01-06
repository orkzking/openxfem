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
   
	Any feedback is welcome. Emails : nvphu80@yahoo.com,stephane.bordas@epfl.ch

*************************************************************************************/

#ifndef _AUXILIARYFIELD_H_
#define _AUXILIARYFIELD_H_

#include "flotmtrx.h"
#include "flotarry.h"
#include "nullmaterial.h"
#include "geometry_base.h"
#include <typeinfo>
#include <cmath>

class CrackTip;class NullMaterial;

typedef enum
{
	Mode_i,
	Mode_ii
} ModeType ;

typedef enum
{
	PlaneStrain,
	PlaneStress
} FieldType ;


//! A template class for auxiliary fields
/*!
	This class is a tempalte class used to compute the interaction integrals.
	Each auxiliary field is attached closely with Materials and 1 CrackTip
	TASK :
	   - compute the derivatives of the displacement field 
		- compute the strain field
		- compute the stress field
 */


template<class M1, class M2, const FieldType field=PlaneStrain>
class AuxiliaryField 
{
public:

  AuxiliaryField(){};             //!< default constructor
  virtual ~AuxiliaryField(){};    //!< destructor

  /*!
   Compute the derivatives of the displacement field, the strain vector and the stress
	vector of the auxiliary field. 
	Two auxiliary fields are implemented including one of homogeneous and 
	bi-material media.
   */
  void ComputeComponentsOfAuxField(CrackTip* tip, Mu::Point* p,ModeType mode, FloatMatrix& AuxGradDisp,FloatArray& AuxEps,FloatArray& AuxStress);
  /*!
   Compute the derivatives of the displacement field of homogeneous material.
	@param : 
	   tip is the crack tip under consideration
	   p is the point where the auxiliary field is evaluated
		mode : deformed modes (mode I or mode II )
	@return : 
	   AuxGradDisp is the derivatives of the displacement field
      AuxEps    is the strain array.
	   AuxStress is the stress array.
   */
  void ComputeComponentsOfOneMat(CrackTip* tip, Mu::Point* p,ModeType mode, FloatMatrix& AuxGradDisp,FloatArray& AuxEps,FloatArray& AuxStress);
  /*!
   Compute the derivatives of the displacement field of bi-material.
	@param : tip is the crack tip under consideration
	p is the point where the auxiliary field is evaluated
	@return : AuxGradDisp is the derivatives of the displacement field
	AuxEps    is the strain array.
	AuxStress is the stress array.
   */
  void ComputeComponentsOfBiMat(CrackTip* tip, Mu::Point* p,ModeType mode, FloatMatrix& AuxGradDisp,FloatArray& AuxEps,FloatArray& AuxStress);

protected:

  M1 *material1;
  M2 *material2;

  size_t giveNumberOfMaterials() const
  {
	 return (size_t)(typeid(M2) != typeid(NullMaterial))+1 ;
  }


};    // end of class AuxiliaryField 


// template class Field used to compute the Kolosov coefficient

// See Thinking in C++Bruce Eckel vol2 
// Compile time selection, Template idioms

template<FieldType field> class Field{};

template<> class Field<PlaneStrain>
{
public:
  static double K(const double nu)
  {
	 return 3.0 - 4.0*nu; 
  }
};

template<> class Field<PlaneStress>
{
public:
  static double K(const double nu)
  {
	 return (3.0 - nu)/(1.0 + nu);
  }
}; 


#endif // _AUXILIARYFIELD_H_
