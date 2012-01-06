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

Any feedback is welcome. Emails : nvphu80@yahoo.com, ...

*************************************************************************************/


#include "standardquadrature.h"
#include "flotarry.h"
#include "mathfem.h"
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>


void StandardGaussLegendreQuadrature :: SetUpPointsOnLine (size_t nPoints)
// creates array of nPoints Gauss Integration Points 
// ( don't confuse with GaussPoint - elem is only the container where to 
//   store corrdinates and weights)
{
  integrationPointVector = new std::vector<Mu::Point*>;   
  weightArray = new std::vector<double>; 

  switch (nPoints) 
  {

	 // 	case 1:
	 //   
	 // // 		coord = new FloatArray(1);
	 // // 		(*coord)[0] = 0.0;
	 // // 		weight = 2.0;
	 // // 
	 // // 		integrationPointVector->push_back(coord);
	 // // 		weightArray->push_back(weight) ;
	 // 		//! Unsupported !!
	 // 		assert(false) ;
	 // 		break;

  case 2:
	 {
		numberOfIntegrationPoints = 2 ;

		weightArray->push_back(1.0) ;
		weightArray->push_back(1.0) ;

		Mu::Point * coord = new Mu::Point(0.577350269189626, -0.577350269189626) ;
		integrationPointVector->push_back(coord) ;	
		break;
	 }

  default:
	 std::cout << "SetUpPointsOnLine: unsupported number of IPs (" << nPoints << ")" << std::endl;
  }
}

void StandardGaussLegendreQuadrature :: SetUpPointsOnTriangle (size_t nPoints)
// creates array of nPoints Gauss Integration Points 
// ( don't confuse with GaussPoint - elem is only the container where to 
//   store corrdinates and weights)
{
  integrationPointVector = new std::vector<Mu::Point*>;   
  weightArray = new std::vector<double>; 

  switch (nPoints) 
  {
  case 1:
	 {
		numberOfIntegrationPoints = 1 ;

		weightArray->push_back(0.5) ;

		Mu::Point * coord1  = new Mu::Point(1/3. , 1/3.);
		integrationPointVector->push_back(coord1) ;
		break;
	 }

  case 3: 
	 {

		numberOfIntegrationPoints = 3 ;

		weightArray->push_back(1/6.) ;
		weightArray->push_back(1/6.) ;
		weightArray->push_back(1/6.) ;

		Mu::Point *coord0 = new Mu::Point(1/6., 1/6.) ;
		Mu::Point *coord1 = new Mu::Point(2/3., 1/6.) ;
		Mu::Point *coord2 = new Mu::Point(1/6., 2/3.) ;

		integrationPointVector->push_back(coord0) ;
		integrationPointVector->push_back(coord1) ;
		integrationPointVector->push_back(coord2) ;


		break;
	 }

  case 7:
	 {
		numberOfIntegrationPoints = 7 ;

		Mu::Point *coord0 = new Mu::Point(0.1012865073235, 0.1012865073235) ;
		Mu::Point *coord1 = new Mu::Point(0.7974269853531, 0.1012865073235) ;
		Mu::Point *coord2 = new Mu::Point(0.1012865073235, 0.7974269853531) ;
		Mu::Point *coord3 = new Mu::Point(0.4701420641051, 0.0597158717898) ;
		Mu::Point *coord4 = new Mu::Point(0.4701420641051,0.4701420641051) ;
		Mu::Point *coord5 = new Mu::Point(0.0597158717898, 0.4701420641051) ;
		Mu::Point *coord6 = new Mu::Point(0.3333333333333, 0.3333333333333) ;

		integrationPointVector->push_back(coord0) ;
		weightArray->push_back(0.1259391805448) ;
		integrationPointVector->push_back(coord1) ;
		weightArray->push_back(0.1259391805448) ;
		integrationPointVector->push_back(coord2) ;
		weightArray->push_back(0.1259391805448) ;
		integrationPointVector->push_back(coord3) ;
		weightArray->push_back(0.1323941527885) ;
		integrationPointVector->push_back(coord4) ;
		weightArray->push_back(0.1323941527885) ;
		integrationPointVector->push_back(coord5) ;
		weightArray->push_back(0.1323941527885) ;
		integrationPointVector->push_back(coord6) ;
		weightArray->push_back(0.2250000000000) ;

		for(size_t i = 0 ; i < weightArray->size() ; i++)
		  (*weightArray)[i] *= 0.5 ;

		break ;
	 }

  case 13:
	 {
		numberOfIntegrationPoints = 13 ;

		Mu::Point *coord0 = new Mu::Point(0.0651301029022, 0.0651301029022) ;
		Mu::Point *coord1 = new Mu::Point( 0.8697397941956, 0.0651301029022) ;
		Mu::Point *coord2 = new Mu::Point( 0.0651301029022, 0.8697397941956) ;
		Mu::Point *coord3 = new Mu::Point(0.3128654960049, 0.0486903154253) ;
		Mu::Point *coord4 = new Mu::Point(0.6384441885698, 0.3128654960049) ;
		Mu::Point *coord5 = new Mu::Point(0.0486903154253, 0.6384441885698) ;
		Mu::Point *coord6 = new Mu::Point(0.6384441885698, 0.0486903154253) ;
		Mu::Point *coord7 = new Mu::Point(0.3128654960049, 0.6384441885698) ;
		Mu::Point *coord8 = new Mu::Point(0.0486903154253, 0.3128654960049) ;
		Mu::Point *coord9 = new Mu::Point(0.2603459660790, 0.2603459660790) ;
		Mu::Point *coord10 = new Mu::Point( 0.4793080678419, 0.2603459660790) ;
		Mu::Point *coord11 = new Mu::Point(0.2603459660790, 0.4793080678419) ;
		Mu::Point *coord12 = new Mu::Point(0.3333333333333, 0.3333333333333) ;

		integrationPointVector->push_back(coord0) ;
		weightArray->push_back(0.0533472356088) ;
		integrationPointVector->push_back(coord1) ;
		weightArray->push_back(0.0533472356088) ;
		integrationPointVector->push_back(coord2) ;
		weightArray->push_back(0.0533472356088) ;
		integrationPointVector->push_back(coord3) ;
		weightArray->push_back(0.0771137608903) ;
		integrationPointVector->push_back(coord4) ;
		weightArray->push_back(0.0771137608903) ;
		integrationPointVector->push_back(coord5) ;
		weightArray->push_back(0.0771137608903) ;
		integrationPointVector->push_back(coord6) ;
		weightArray->push_back(0.0771137608903) ;
		integrationPointVector->push_back(coord7) ;
		weightArray->push_back(0.0771137608903) ;
		integrationPointVector->push_back(coord8) ;
		weightArray->push_back(0.0771137608903) ;
		integrationPointVector->push_back(coord9) ;
		weightArray->push_back(0.1756152576332) ;
		integrationPointVector->push_back(coord10) ;
		weightArray->push_back(0.1756152576332) ;
		integrationPointVector->push_back(coord11) ;
		weightArray->push_back(0.1756152576332) ;
		integrationPointVector->push_back(coord12) ;
		weightArray->push_back(-0.1495700444677) ;

		for(size_t i = 0 ; i < weightArray->size() ; i++)
		  (*weightArray)[i] *= 0.5 ;

		break;
	 }
  default:
	 printf ("SetUpPointsOnTriangle: unsupported number of IPs (%d)", nPoints);
  }
}

void StandardGaussLegendreQuadrature :: SetUpPointsOnSquare(size_t nPoints) 
// creates array of nPoints Gauss Integration Points 
// ( don't confuse with GaussPoint - elem is only the container where to 
//   store corrdinates and weights)
{
  integrationPointVector = new std::vector<Mu::Point*>;   
  weightArray = new std::vector<double>; 

  switch (nPoints) 
  {

  case 1:
	 {
		numberOfIntegrationPoints = 1 ;

		weightArray->push_back(4.0) ;

		Mu::Point  *coord = new Mu::Point(0,0) ;
		integrationPointVector->push_back(coord) ;

		break;
	 }
  case 4:
	 {
		numberOfIntegrationPoints = 4 ;

		Mu::Point  *coord0 = new Mu::Point(-0.577350269189626,-0.577350269189626) ;
		Mu::Point  *coord1 = new Mu::Point(0.577350269189626,-0.577350269189626) ;
		Mu::Point  *coord2 = new Mu::Point(0.577350269189626,0.577350269189626) ;
		Mu::Point  *coord3 = new Mu::Point(-0.577350269189626,0.577350269189626) ;


		integrationPointVector->push_back(coord0) ;
		weightArray->push_back(1.0) ;
		integrationPointVector->push_back(coord1) ;
		weightArray->push_back(1.0) ;
		integrationPointVector->push_back(coord2) ;
		weightArray->push_back(1.0) ;
		integrationPointVector->push_back(coord3) ;
		weightArray->push_back(1.0) ;

		break;
	 }
	 
  case 9:  // code 2005-08-23 for Q4 XFEM implementation.
	 {
		numberOfIntegrationPoints = 9 ;

		std::valarray<double> c(3);
		std::valarray<double> w(3);

		c[0] = -0.774596669241483;
		c[1] =  0.0;
		c[2] =  0.774596669241483;

		w[0] =  0.555555555555555;
		w[1] =  0.888888888888888;
		w[2] =  0.555555555555555;

		for (size_t i = 0 ; i < 3 ; i++)
		{
		  for (size_t j = 0 ; j < 3 ; j++) 
		  {
			 Mu::Point  *coord = new Mu::Point(c[i],c[j]) ;
			 integrationPointVector->push_back(coord) ;
			 weightArray->push_back(w[i] * w[j]);
		  }
		}

		break;
	 }

	 case 16:  // 25-10-2005
	 {
		numberOfIntegrationPoints = 16 ;

		std::valarray<double> c(4);
		std::valarray<double> w(4);

		c[0] =  0.861134311594053;
		c[1] = -0.861134311594053;
		c[2] =  0.339981043584856;
		c[3] = -0.339981043584856;

		w[0] =  0.347854845137454;
		w[1] =  0.347854845137454; 
		w[2] =  0.652145154862546;
		w[3] =  0.652145154862546;  

		for (size_t i = 0 ; i < 4 ; i++)
		{
		  for (size_t j = 0 ; j < 4 ; j++) 
		  {
			 Mu::Point  *coord = new Mu::Point(c[i],c[j]) ;
			 integrationPointVector->push_back(coord) ;
			 weightArray->push_back(w[i] * w[j]);
		  }
		}

		break;
	 }

	 case 25:  // 25-10-2005, nucleation crack implementation!!!
	 {
		numberOfIntegrationPoints = 25 ;

		std::valarray<double> c(5);
		std::valarray<double> w(5);

		c[0] =  0.906179845938664;
		c[1] = -0.906179845938664;
		c[2] =  0.538469310105683;
		c[3] = -0.538469310105683;
		c[4] =  0.000000000000000;

		w[0] =  0.236926885056189;
		w[1] =  0.236926885056189;
		w[2] =  0.478628670499366;
		w[3] =  0.478628670499366;  
		w[4] =  0.568888888888889;  

		for (size_t i = 0 ; i < 5 ; i++)
		{
		  for (size_t j = 0 ; j < 5 ; j++) 
		  {
			 Mu::Point  *coord = new Mu::Point(c[i],c[j]) ;
			 integrationPointVector->push_back(coord) ;
			 weightArray->push_back(w[i] * w[j]);
		  }
		}

		break;
	 }
	 case 36:  // 25-10-2005, nucleation crack implementation!!!
	 {
		numberOfIntegrationPoints = 36 ;

		std::valarray<double> c(6);
		std::valarray<double> w(6);

		c[0] =  0.932469514203152;
		c[1] = -0.932469514203152;
		c[2] =  0.661209386466265;
		c[3] = -0.661209386466265;
		c[4] =  0.238619186003152;
		c[5] = -0.238619186003152;

		w[0] =  0.171324492379170;
		w[1] =  0.171324492379170;
		w[2] =  0.360761573048139;
		w[3] =  0.360761573048139;   
		w[4] =  0.467913934572691; 
		w[5] =  0.467913934572691;

		for (size_t i = 0 ; i < 6 ; i++)
		{
		  for (size_t j = 0 ; j < 6 ; j++) 
		  {
			 Mu::Point  *coord = new Mu::Point(c[i],c[j]) ;
			 integrationPointVector->push_back(coord) ;
			 weightArray->push_back(w[i] * w[j]);
		  }
		}

		break;
	 }
	 case 49:  // 25-10-2005, nucleation crack implementation!!!
	 {
		numberOfIntegrationPoints = 49 ;

		std::valarray<double> c(7);
		std::valarray<double> w(7);

		c[0] =  0.949107912342759;
		c[1] = -0.949107912342759;
		c[2] =  0.741531185599394;
		c[3] = -0.741531185599394;
		c[4] =  0.405845151377397;
		c[5] = -0.405845151377397;
		c[6] =  0.000000000000000;

		w[0] =  0.129484966168870;
		w[1] =  0.129484966168870;
		w[2] =  0.279705391489277;
		w[3] =  0.279705391489277; 
		w[4] =  0.381830050505119;
		w[5] =  0.381830050505119;
		w[6] =  0.417959183673469;

		for (size_t i = 0 ; i < 7 ; i++)
		{
		  for (size_t j = 0 ; j < 7 ; j++) 
		  {
			 Mu::Point  *coord = new Mu::Point(c[i],c[j]) ;
			 integrationPointVector->push_back(coord) ;
			 weightArray->push_back(w[i] * w[j]);
		  }
		}

		break;
	 }
	  case 56:  // 25-10-2005, nucleation crack implementation!!!
	 {
		numberOfIntegrationPoints = 56 ;

		std::valarray<double> c(8);
		std::valarray<double> w(8);

		c[0] =  0.960289856497536;
		c[1] = -0.960289856497536;
		c[2] =  0.796666477413627;
		c[3] = -0.796666477413627;
		c[4] =  0.525532409916329;
		c[5] = -0.525532409916329;
		c[6] =  0.183434642495650;
		c[7] = -0.183434642495650;

		w[0] = 0.101228536290376;
		w[1] = 0.101228536290376;
		w[2] = 0.222381034453374;
		w[3] = 0.222381034453374;
		w[4] = 0.313706645877887;
		w[5] = 0.313706645877887;
		w[6] = 0.362683783378362;
		w[7] = 0.362683783378362;

		for (size_t i = 0 ; i < 8 ; i++)
		{
		  for (size_t j = 0 ; j < 8 ; j++) 
		  {
			 Mu::Point  *coord = new Mu::Point(c[i],c[j]) ;
			 integrationPointVector->push_back(coord) ;
			 weightArray->push_back(w[i] * w[j]);
		  }
		}

		break;
	 }
	 
  default:
	 std::cout << "SetUpPointsOnSquare: unsupported number of IPs (" << nPoints << ")" << std::endl ;
  }
}

void StandardGaussLegendreQuadrature :: SetUpPointsOnTetrahedra(size_t nPoints)
// creates array of nPoints Gauss Integration Points 
// ( don't confuse with GaussPoint - elem is only the container where to 
//   store corrdinates and weights)
{
  integrationPointVector = new std::vector<Mu::Point*>;   
  weightArray = new std::vector<double>; 

  switch (nPoints) 
  {
  case 1:
	 {
		numberOfIntegrationPoints = 1 ;
		weightArray->push_back(1.0) ;
		Mu::Point * coord = new Mu::Point(0.25,0.25,0.25) ;
		integrationPointVector->push_back(coord) ;	
		break;
  
	 }
  case 4:
	 {
		numberOfIntegrationPoints = 4 ;

		weightArray->push_back(0.25) ;
		weightArray->push_back(0.25) ;
		weightArray->push_back(0.25) ;
		weightArray->push_back(0.25) ;

		Mu::Point * point1 = new Mu::Point(0.58541020,0.13819660,0.13819660) ;
		Mu::Point * point2 = new Mu::Point(0.13819660,0.58541020,0.13819660) ;
		Mu::Point * point3 = new Mu::Point(0.13819660,0.13819660,0.58541020) ;
		Mu::Point * point4 = new Mu::Point(0.13819660,0.13819660,0.13819660) ;

		integrationPointVector->push_back(point1) ;	
		integrationPointVector->push_back(point2) ;	
		integrationPointVector->push_back(point3) ;	
		integrationPointVector->push_back(point4) ;	
		break;
	 }

  default:
	 std::cout << "SetUpPointsOnTetrahedra: unsupported number of IPs (" << nPoints << ")" << std::endl;
  }
}

int
StandardGaussLegendreQuadrature::getRequiredNumberOfIntegrationPoints (IntDomain dType, int approxOrder) 
{
  int requiredNIP ;

  if (approxOrder <0) return 0;

  if (dType == LINE){
	 requiredNIP = (approxOrder+1)/2;
	 if (requiredNIP > 4) return -1;
	 return requiredNIP;}

  if (dType == ::TRIANGLE){
	 if (approxOrder <= 1) return 1;
	 if (approxOrder <= 3) return 4;
	 if (approxOrder <= 5) return 7;
	 return -1;}

  if (dType == SQUARE){
	 requiredNIP = std::max((approxOrder+1)/2, 2); 
	 requiredNIP *= requiredNIP;
	 if (requiredNIP > 16) return -1;
	 return requiredNIP;}
  else
	 printf ("GaussIntegrationRule::setUpIntegrationPoints - unknown integrationDomain");
  return -1;
}
