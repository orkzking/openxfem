//   file STANDARDDESCRIPTION.CPP

#include "standarddescription.h"

#include "element.h"
#include "geometry_base.h"
#include "node.h"
#include "tri_u.h"
#include "geometryentity.h"
#include "piecewiselinear.h"
#include <iostream>
#include <typeinfo>

/*
bool StandardDescription::interactsWith(Element *elem, GeometryEntity* geo)
//*************************************************************************
// current implementation just deals with Tri_U
{
  Tri_U *tri = static_cast<Tri_U*>(elem);
  Mu::Triangle *t = tri->makeTriangle();

  return geo->intersects(t);
} 
*/

bool StandardDescription::interactsWith(Element *elem, GeometryEntity* geo)
//*************************************************************************
// New implementation. 2005-08-23
// works for T3 and Q4 elements.
// \todo code for T6 elements or high order elements !!!
{
  std::vector<Mu::Point> pts ;
  for(size_t i = 0 ; i < elem->giveNumberOfNodes() ; i++)
  {
	 Mu::Point *p = elem->giveNode(i+1)->makePoint() ;
	 pts.push_back(*p);
	 delete p ;
  }

  bool ret = false ;	

  for(size_t i = 0 ; i < pts.size() ;  i++)
  {
	 Mu::Segment s(pts[i],pts[(i+1)%pts.size()]) ;
	 ret = ret || geo->intersects(&s) ;
	 if(ret) 
		i = pts.size() ; // no need to test more 
  }

  return ret ;
  
} 

/*
bool StandardDescription::interactsWith(Element *elem, GeometryEntity* geo)
{
  FloatArray *coord;

  size_t Above = 0; 
  size_t Below = 0;

  for (size_t i = 1 ; i <= elem->giveNumberOfNodes() ; i++)
  {
	 coord = elem->giveNode(i)->giveCoordinates();
	 Mu::Point *P = new Mu::Point((*coord)[0],(*coord)[1]);
	 int test = geo->givePositionComparedTo(P);
	 switch(test)
	 {
	 case  1 :
		Above += 1; 
		break ; 
	 case -1 :
		Below += 1;
		break ; 
	 case  0 :
		std::cout<< "Node is on the discontinuity" << std::endl;
	 }
  }
  
  return ( (Above != elem->giveNumberOfNodes()) && (Below != elem->giveNumberOfNodes()) );
} 
*/



