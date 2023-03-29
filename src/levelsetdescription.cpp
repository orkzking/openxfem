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

#include "levelsetdescription.h"

#include "element.h"
#include "node.h"
#include "geometry_base.h"
#include "geometryentity.h"
#include "flotarry.h"
#include <math.h>
#include <vector>
#include <algorithm>
    

bool LevelSetDescription::interactsWith(Element *elem, GeometryEntity* geo)
// ************************************************************************
// computes the LS(signed distance) at nodes of elem, if at least there are
// two LS with different signs => elem cut by geo, see Stolarska 2001.
{
	double dMax,dMin;
	Mu::Point* coord = new Mu::Point();
	std::vector<double> distVector;


	for (size_t i=0 ; i< elem->giveNumberOfNodes() ; i++)
	{
		coord->x = elem->giveNode(i+1)->giveCoordinate(1);
		coord->x = elem->giveNode(i+1)->giveCoordinate(2);
		distVector.push_back(geo->computeSignedDistanceOfPoint(coord)) ;
    }

	 typedef std::vector<double>::iterator vec_iter ;

	 vec_iter start,end;
	 
	 start = distVector.begin(); 
	 end = distVector.end();

	 dMax = *(max_element(start,end));
	 dMin = *(min_element(start,end));

	 return (dMax*dMin<0)? true : false ;
}

