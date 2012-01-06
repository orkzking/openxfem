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
// file SPLITGAUSSQUADRATURE.CPP

#include "splitgaussquadrature.h"

#include "flotarry.h"
#include "flotmtrx.h"
#include "delaunay.h"
#include "element.h"
#include "fei2dtrilin.h"
#include "standardquadrature.h"
#include <ctype.h>
#include <stdlib.h>
#include <vector>
#include <iostream>


void SplitGaussLegendreQuadrature::setUpIntegrationPoints(Element* elem,size_t nPoints)
{
  this->SetUpPointsOn2dDomain(elem,nPoints);
  // this->SetUpPointsOn3dDomain(elem); // WILL BE IMPROVED LATER !!!
}

void SplitGaussLegendreQuadrature::SetUpPointsOn2dDomain(Element* elem,size_t nPoints)
// ************************************************************************************
// Computes the Gauss points for partitioned element
// using the sub-triangles got from the method Element::PartitionMySelf()
{
  // get the subtriangles which already defined in local coordinate system
  
  std::vector<DelaunayTriangle *> * tri = elem->PartitionMySelf();
  
  /*// for debugging ...
  std::cout << " This element is partitioned into " << tri->size() << " sub-triangles " << std::endl;
  for(size_t i = 0 ; i < tri->size() ; i++)
  {
	 (*tri)[i]->print() ;
	 std::cout << (*tri)[i]->area() << endl ;
  } */

  // get Gauss points for each sub-triangles
  StandardGaussLegendreQuadrature* standardQuad = new StandardGaussLegendreQuadrature;
  standardQuad->setUpIntegrationPoints(::TRIANGLE,nPoints);

  std::vector<double>*     weight;
  std::vector<Mu::Point*>* coordArray;  

  weight = standardQuad->giveWeightArray();
  coordArray = standardQuad->giveIntegrationPointVector();	 

  // initialization ...
  weightArray = new std::vector<double> ; 
  integrationPointVector = new std::vector<Mu::Point*>   ;
  FEI2dTriLin* FEInterpol = new FEI2dTriLin(1,2); 

  FloatMatrix* coord ;
  // ********************************************************************************
  // ******                     LOOP ON SUB-TRIANGLES                          ******
  // ********************************************************************************
  for(size_t j = 0 ; j < tri->size() ; j++ )      
  {  
	 if( ((*tri)[j]->isTriangle) && ((*tri)[j]->area() > 1e-4 ) ) // new !!! 2005-09-13
	 { 
		coord = new FloatMatrix(3,2);
		coord->at(1,1) = (*tri)[j]->first->x;   // coord. of first vertex
		coord->at(1,2) = (*tri)[j]->first->y;

		coord->at(2,1) = (*tri)[j]->second->x;  // coord. of second vertex
		coord->at(2,2) = (*tri)[j]->second->y;

		coord->at(3,1) = (*tri)[j]->third->x;   // coord. of third vertex
		coord->at(3,2) = (*tri)[j]->third->y;

		double area = (*tri)[j]->area();        // area of subtriangle
       
		// ********************************************************************************
		// ******            LOOP ON GAUSS POINTS OF EACH SUB-TRIANGLE               ******
		// ********************************************************************************
		for (size_t i = 0 ; i < coordArray->size() ; i++) 
		{

#ifdef NVP5JULY05
		  std::cout << "   " << " - the " << i+1 << " Gauss point : " << std::endl;
		  std::cout << "        " ;
		  (*coordArray)[i]->print();
        std::cout << std::endl;
#endif    

		  // transform this Gauss point to the global coord : \sum{N_i*x_i}
		  Mu::Point* newIntPoint = new Mu::Point();
		  FloatArray *temp ;
		  temp = FEInterpol->evalN((*coordArray)[i]);

		  for (size_t j = 0 ; j < 3 ; j++)
		  {
			 double x = coord->at(j+1,1) ;
			 double y = coord->at(j+1,2) ;
			 newIntPoint->x += (*temp)[j]*x;
			 newIntPoint->y += (*temp)[j]*y;
		  }
		  delete temp; // Purify, 11-10-05

		  double newWeight = 2.0 * weight->at(i) * area;       // !!!  2.0 *

		  integrationPointVector->push_back(newIntPoint);
		  weightArray->push_back(newWeight);

#ifdef NVP5JULY05
		  for (size_t i = 0 ; i < integrationPointVector->size() ; i++)
		  {		
			 std::cout << "   " << " - the new " << i+1 << " Gauss point : " << std::endl;
			 std::cout << "        " ;
			 (*integrationPointVector)[i]->print();
			 std::cout << std::endl;
		  }		  
#endif
		}    // end of loop on GPs
		delete coord ;                        // Purify, 11-10-05
	 }
  }        // end of loop on sub-triangles

  delete weight ; delete coordArray ;       // Purify, 11-10-05
  delete FEInterpol ;                       // Purify, 11-10-05
  delete tri ;
}