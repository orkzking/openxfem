//
// C++ Implementation: delaunay
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "delaunay.h"
#include <algorithm>

#define DEBUG
#undef DEBUG

DelaunayTreeItem::DelaunayTreeItem( DelaunayTreeItem * father,  const Point * c)
{
	this->father = father;
	m_c  = c ;
	dead = false ;
	visited = false ;
	
	isPlane =false ;
	isTriangle =false ;
}
	
DelaunayTreeItem::~DelaunayTreeItem()
{
// 	for(size_t i = 0 ; i < neighbour.size() ; i++)
// 	{
// 		neighbour[i]->removeNeighbour(this) ;
// 	}
}
	
const Point * DelaunayTreeItem::killer() const 
{
	return m_k ;
}

const Point * DelaunayTreeItem::creator() const 
{
	return m_c ;
}

void  DelaunayTreeItem::setCreator(const Point * p)  
{
	m_c = p ;
}

void Star::updateNeighbourhood()
{
	std::vector<DelaunayTreeItem *> sons ;
	std::vector<DelaunayTreeItem *> neighbouringPlanes ;
	
	for(size_t i = 0 ; i < treeitem.size() ;i++)
	{
		sons.insert(sons.end() , treeitem[i]->son.begin() , treeitem[i]->son.end()) ;
		for(size_t j = 0 ; j < treeitem[i]->neighbour.size() ; j++)
		{
			if (treeitem[i]->neighbour[j]->isPlane)
				neighbouringPlanes.push_back(treeitem[i]->neighbour[j]) ;
		}
	}
	
	for(size_t i = 0 ; i < this->size() ;i++)
	{
		for(size_t j = 0 ; j < sons.size() ;j++)
		{
			for(size_t k = 0 ; k < sons.size() ;k++)
			{
				if(sons[j]->isVertex(edge[i]) && sons[k]->isVertex(edge[i])&& sons[j]->isAlive() && sons[k]->isAlive())
				{
						makeNeighbours(sons[j], sons[k]) ;
				}
			}
			
			for(size_t k = 0 ; k < neighbouringPlanes.size() ; k++)
			{
				if (neighbouringPlanes[k]->isNeighbour(sons[j]))
					makeNeighbours(neighbouringPlanes[k], sons[j]) ;
			}
		}
	}
}

void DelaunayTree::addSharedNodes(size_t nodes_per_side)
{
	for(size_t i = 0 ; i < this->tree.size() ; i++)
	{
		if(this->tree[i]->isAlive() && this->tree[i]->isTriangle)
		{
			for(size_t j = 0 ; j < this->tree[i]->neighbour.size() ; j++)
			{
				shareNodes(this->tree[i]->neighbour[j], this->tree[i], nodes_per_side) ;
			}
		}
	}
	
	for(size_t i = 0 ; i < this->tree.size() ; i++)
	{
		if(this->tree[i]->isAlive() && this->tree[i]->isTriangle)
		{
			this->global_counter += static_cast<DelaunayTriangle *>(this->tree[i])->insertSharedNodes(this->global_counter) ;
		}
	}
}

size_t DelaunayTriangle::insertSharedNodes(size_t start_id)
{
	assert(first->id > -1) ;
	assert(second->id > -1) ;
	assert(third->id > -1) ;
	size_t nodes_per_side = sharedPoints[0].second.size() ;
	this->Triangle::boundingPoints->resize(sharedPoints.size()*nodes_per_side + 3) ;
	size_t index = 0 ;
	(*this->Triangle::boundingPoints)[0] = first ; 
	(*this->Triangle::boundingPoints)[nodes_per_side+1] = second ; 
	(*this->Triangle::boundingPoints)[nodes_per_side*2+2] = third ; 
	
	std::vector<Point> placeholders ;
	placeholders.push_back((*first)) ;
	
// 	for(size_t i = 0 ; i < nodes_per_side ; i++)
// 	{
		placeholders.push_back((*second) * 0.5 + (*first) * 0.5) ;
// 	}
	placeholders.push_back((*second)) ;
// 	for(size_t i = 0 ; i < nodes_per_side ; i++)
// 	{
		placeholders.push_back((*third) * 0.5 + (*second) * 0.5) ;
// 	}
	placeholders.push_back((*third)) ;
// 	for(size_t i = 0 ; i < nodes_per_side ; i++)
// 	{
		placeholders.push_back((*first) * 0.5 + (*third) * 0.5) ;
// 	}
	
	for(size_t i = 0 ; i < sharedPoints.size() ; i++)
	{
		for(size_t j = 0 ; j <sharedPoints[i].second.size() ; j++)
		{
			bool in = false ;
			for(size_t k = 0 ; k < placeholders.size() ; k++)
			{
				if((*sharedPoints[i].second[j]) == placeholders[k])
				{
					in = true ;
					(*this->boundingPoints)[k] = sharedPoints[i].second[j] ;
					if (sharedPoints[i].second[j]->id == -1)
						sharedPoints[i].second[j]->id = start_id + index++ ;
					break ;
				}
			}
			assert(in) ;
		}
	}
	for(size_t i  = 0 ; i < this->Triangle::getBoundingPoints()->size() ; i++)
		assert(this->Triangle::getBoundingPoint(i)->id > -1) ;
	return index ;
}

size_t DelaunayTree::numPoints() const
{
	return this->global_counter ;
}

void shareNodes(DelaunayTreeItem *t0,DelaunayTreeItem *t1, size_t nodes_per_side)
{
	
	std::pair<Point*, Point*> seg ;
	if(t0->isTriangle)
		seg = t0->commonEdge(t1) ;
	else if(t1->isTriangle)
		seg = t1->commonEdge(t0) ;
	else
		return ;
	
	std::vector<Point*> shared ;
	
	for(size_t i = 0 ; i < t0->sharedPoints.size() ; i++)
		if(t0->sharedPoints[i].first == t1)
			return ;
	
// 	for(size_t i = 0 ; i < nodes_per_side ; i++)
// 	{
		shared.push_back(new Point( (*seg.first) * 0.5 + (*seg.second) * 0.5)) ;
// 	}
	(*shared.rbegin())->id = -1 ;
	
	t0->sharedPoints.push_back(std::pair<DelaunayTreeItem *, std::vector<Point*> >(t1, shared)) ;
	t1->sharedPoints.push_back(std::pair<DelaunayTreeItem *, std::vector<Point*> >(t0, shared)) ;

}

std::pair<std::vector<DelaunayTriangle *>, std::vector<DelaunayTreeItem *> > DelaunayTreeItem::conflicts( const Geometry *g)
{
	std::pair<std::vector<DelaunayTriangle *>, 
		std::vector<DelaunayTreeItem *> > ret ;

	if(visited)
	{
		return ret ;
	}
	visited = true ;
	
	ret.second.push_back(this) ;
	
	
	if(!g->in(first) && !g->in(second) && ( !g->in(third) && isTriangle ) || ( g->in(third) && isPlane ) )
	{
		return ret ;
	}
	
	for (size_t i  = 0 ;  i < stepson.size() ; i++)
	{
		if( !stepson[i]->visited) 
		{
			std::pair<std::vector<DelaunayTriangle *>, std::vector<DelaunayTreeItem *> > temp  =  stepson[i]->conflicts(g) ;
			ret.first.insert(ret.first.end(), temp.first.begin(), temp.first.end()) ;
			ret.second.insert(ret.second.end(), temp.second.begin(), temp.second.end()) ;
		}
	}
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		if( !son[i]->visited )
		{
			std::pair<std::vector<DelaunayTriangle *>, std::vector<DelaunayTreeItem *> >  temp  = son[i]->conflicts(g) ;
			ret.first.insert(ret.first.end(), temp.first.begin(), temp.first.end()) ;
			ret.second.insert(ret.second.end(), temp.second.begin(), temp.second.end()) ;
		}
	}
	
	if(isAlive() && isTriangle)
	{
		ret.first.push_back(dynamic_cast<DelaunayTriangle *>(this)) ;
	}
	
	for (size_t i  = 0 ;  i < neighbour.size() ; i++)
	{
		if( !neighbour[i]->visited )
		{
			std::pair<std::vector<DelaunayTriangle *>, std::vector<DelaunayTreeItem *> > temp  = neighbour[i]->conflicts(g) ;
			ret.first.insert(ret.first.end(), temp.first.begin(), temp.first.end()) ;
			ret.second.insert(ret.second.end(), temp.second.begin(), temp.second.end()) ;
// 			for (size_t j  = 0 ;  j < neighbour[i]->neighbour.size() ; j++)
// 			{
// 				if( !neighbour[i]->neighbour[j]->visited )
// 				{
// 					temp  = neighbour[i]->neighbour[j]->conflicts(g) ;
// 					ret->first.insert(ret->first.end(), temp->first.begin(), temp->first.end()) ;
// 					ret->second.insert(ret->second.end(), temp->second.begin(), temp->second.end()) ;
// 					delete temp;
// 				}
// 			}
			
		}
	}
	return ret ;
}

std::pair<std::vector<Point *>, std::vector<DelaunayTreeItem *> > * DelaunayTreeItem::conflicts( const Segment *s)  
{
	std::pair<std::vector<Point *>, 
		std::vector<DelaunayTreeItem *> > * ret = new std::pair<std::vector<Point *>, 
		std::vector<DelaunayTreeItem *> > ;
	std::pair<std::vector<Point *>, std::vector<DelaunayTreeItem *> > * temp ;
	
	if(visited)
		return ret ;
	visited = true ;
	ret->second.push_back(this) ;
	
	if(isTriangle && (isVertex(s->first()) || isVertex(s->second())))
		return ret ;
	
	if(isPlane && ( isVertex(s->first()) || isVertex(s->second() ) ) )
		return ret ;
	
	if(isTriangle)
	{
		std::vector<Point *> intersections ;
		Segment s0((*first), (*second)) ;
		Segment s1((*second), (*third)) ;
		Segment s2((*third), (*first)) ;
		
		if (s->intersects(&s0) )
		{
			Point * i = new Point(s->intersection(&s0)) ;
			intersections.push_back(i) ;
		}
		if (s->intersects(&s1) )
		{
			Point * i = new Point(s->intersection(&s1)) ;
			intersections.push_back(i) ;
		}
		if (s->intersects(&s2) )
		{
			Point * i = new Point(s->intersection(&s2)) ;

			intersections.push_back(i) ;
		}
			
		if(intersections.size() == 2)
		{
			Point *i = new Point( (*intersections[0])*0.5 + (*intersections[1])*0.5 ) ;
			ret->first.push_back(intersections[0]) ;
			ret->first.push_back(i) ;
			ret->first.push_back(intersections[1]) ;
			
		}
		if(intersections.size() == 1)
			ret->first.push_back(intersections[0]) ;
		
	}
	else if(isPlane)
	{
		Segment s0((*first), (*second)) ;
		if(s->intersects(&s0))
		{
			Point * i = new Point(s->intersection(&s0)) ;

			ret->first.push_back(i) ;
		}
	}
	
	if(!isAlive())
	{
		for (size_t i  = 0 ;  i < ret->first.size() ; i++)
			delete ret->first[i] ;
		
		ret->first.clear() ;
	}
	
	for (size_t i  = 0 ;  i < stepson.size() ; i++)
	{
		if( !stepson[i]->visited) 
		{
			temp  =  stepson[i]->conflicts(s) ;
			for(size_t j  = 0 ; j < temp->first.size() ;j++)
			{
				bool unique = true ;
				for(size_t k  = 0 ; k < ret->first.size() && unique == true ;k++)
				{
					if((*temp->first[j]) == (*ret->first[k]))
						unique = false ; 
				}
				
				if(unique)
					ret->first.push_back(temp->first[j]) ;
			}
			ret->second.insert(ret->second.end(), temp->second.begin(), temp->second.end()) ;
			delete temp;
		}
	}
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		if( !son[i]->visited )
		{
			temp  = son[i]->conflicts(s) ;
			for(size_t j  = 0 ; j < temp->first.size() ;j++)
			{
				bool unique = true ;
				for(size_t k  = 0 ; k < ret->first.size() && unique == true ;k++)
				{
					if((*temp->first[j]) == (*ret->first[k]))
						unique = false ; 
				}
				
				if(unique)
					ret->first.push_back(temp->first[j]) ;
			}
			ret->second.insert(ret->second.end(), temp->second.begin(), temp->second.end()) ;
			delete temp;
		}
	}
	
	for (size_t i  = 0 ;  i < neighbour.size() ; i++)
	{
		if( !neighbour[i]->visited )
		{
			temp  = neighbour[i]->conflicts(s) ;
			for(size_t j  = 0 ; j < temp->first.size() ;j++)
			{
				bool unique = true ;
				for(size_t k  = 0 ; k < ret->first.size() && unique == true ;k++)
				{
					if((*temp->first[j]) == (*ret->first[k]))
						unique = false ; 
				}
				
				if(unique)
					ret->first.push_back(temp->first[j]) ;
			}
			ret->second.insert(ret->second.end(), temp->second.begin(), temp->second.end()) ;
			delete temp;
		}
	}
	
	return ret ;
}


std::pair<std::vector<DelaunayTreeItem *>, std::vector<DelaunayTreeItem *> > * DelaunayTreeItem::conflicts(const Point *p)
{
	std::pair<std::vector<DelaunayTreeItem *>, 
		std::vector<DelaunayTreeItem *> > * ret = new std::pair<std::vector<DelaunayTreeItem *>, 
		std::vector<DelaunayTreeItem *> > ;
	
	if(visited)
		return ret ;
	visited = true ;
	ret->second.push_back(this) ;

	if(isTriangle && isVertex(p))
		return ret ;
	
	if(isPlane && (*p == *first || *p == *second) )
		return ret ;
	
	if(!inCircumCircle(p))
		return ret ;
	
	for (size_t i  = 0 ;  i < stepson.size() ; i++)
	{
		if( !stepson[i]->visited) 
		{
			std::pair<std::vector<DelaunayTreeItem *>, std::vector<DelaunayTreeItem *> > * temp  =  stepson[i]->conflicts(p) ;
			ret->first.insert(ret->first.end(), temp->first.begin(), temp->first.end()) ;
			ret->second.insert(ret->second.end(), temp->second.begin(), temp->second.end()) ;
			delete temp;
		}
	}
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		if( !son[i]->visited )
		{
			std::pair<std::vector<DelaunayTreeItem *>, std::vector<DelaunayTreeItem *> > * temp  = son[i]->conflicts(p) ;
			ret->first.insert(ret->first.end(), temp->first.begin(), temp->first.end()) ;
			ret->second.insert(ret->second.end(), temp->second.begin(), temp->second.end()) ;
			delete temp;
		}
	}

	if(isAlive())
		ret->first.push_back(this) ;
	
	for (size_t i  = 0 ;  i < neighbour.size() ; i++)
	{
		if( !neighbour[i]->visited )
		{
			std::pair<std::vector<DelaunayTreeItem *>, std::vector<DelaunayTreeItem *> > * temp  = neighbour[i]->conflicts(p) ;
			ret->first.insert(ret->first.end(), temp->first.begin(), temp->first.end()) ;
			ret->second.insert(ret->second.end(), temp->second.begin(), temp->second.end()) ;
			delete temp;
		}
	}

	return ret ;
}

void DelaunayTreeItem::removeNeighbour(DelaunayTreeItem * t)
{
	if(std::find(neighbour.begin(), neighbour.end(),t) !=  neighbour.end())
	{
		neighbour.erase(std::find(neighbour.begin(), neighbour.end(), t)) ;
		deadneighbour.push_back(t) ;
	}
}
	
void DelaunayTreeItem::addNeighbour(DelaunayTreeItem * t)
{
	assert(t != this) ;
	if(std::find(neighbour.begin(), neighbour.end(), t) != neighbour.end())
	{
		return ;
	}
	
	if(t->isAlive())
		neighbour.push_back(t) ;
}
	

DelaunayTreeItem * DelaunayTreeItem::getNeighbour(size_t i)
{
	return this->neighbour[i] ;
}


void DelaunayTreeItem::erase(const Point * p)
{
	dead = true ;
	m_k  = p ;
}

void DelaunayTreeItem::kill(const Point * p)
{
	dead = true ;
	m_k  = p ;
	
	for(size_t i = 0 ; i < this->neighbour.size() ; i++)
	{
		this->neighbour[i]->removeNeighbour(this) ;
	}
}
	
bool DelaunayTreeItem::isAlive() const
{
	return !this->dead ;
}
	
void DelaunayTreeItem::addStepson(DelaunayTreeItem * s)
{
	if(s == this)
		return ;
	stepson.push_back(s) ;
	s->setStepfather(this) ;
	addNeighbour(s) ;
}

void DelaunayTreeItem::addSon(DelaunayTreeItem * s)
{
	son.push_back(s) ;
}
	
void DelaunayTreeItem::removeSon(DelaunayTreeItem * t)
{
	if(std::find(son.begin(), son.end(),t) !=  son.end())
	{
		son.erase(std::find(son.begin(), son.end(), t)) ;
	}
}

void DelaunayTreeItem::removeStepson(DelaunayTreeItem * t)
{
	if(std::find(stepson.begin(), stepson.end(),t) !=  stepson.end())
	{
		stepson.erase(std::find(stepson.begin(), stepson.end(), t)) ;
	}
}

void DelaunayTreeItem::setStepfather(DelaunayTreeItem * s)
{
	stepfather = s ;
	addNeighbour(s) ;
}
	
void DelaunayTreeItem::clearVisited()
{
	visited = false ;
}

DelaunayTriangle::DelaunayTriangle(DelaunayTreeItem * father,  Point *p0,  Point *p1, Point *p2,  Point * c) : Triangle(p0, p1, p2), DelaunayTreeItem(father, c)
{
	first = getBoundingPoint(0) ;
	second = getBoundingPoint(1) ;
	third = getBoundingPoint(2) ;
	
	assert(in(this->getCenter())) ;
	
	neighbour.clear() ;
	isPlane = false ;
	isTriangle = true ;
	assert(first->id > -1) ;
	assert(second->id > -1) ;
	assert(third->id > -1) ;
}

DelaunayTriangle::DelaunayTriangle() : Triangle(), DelaunayTreeItem(NULL, NULL)
{
	first = getBoundingPoint(0) ;
	second = getBoundingPoint(1) ;
	third = getBoundingPoint(2) ;
	
	assert(in(this->getCenter())) ;
	
	neighbour.clear() ;
	isPlane = false ;
	isTriangle = true ;
}
	
DelaunayTriangle::~DelaunayTriangle()
{
}
	
inline bool DelaunayTriangle::isVertex(const Point * p) const
{
	return ((*p) == (*first) || (*p) == (*second) || (*p) == (*third)) ;
}
	

std::pair<Point*, Point*> DelaunayTriangle::commonEdge(const DelaunayTreeItem * t) 
{
	if(t == this)
		return std::pair<Point*, Point*>(reinterpret_cast<Point*>(NULL), reinterpret_cast<Point*>(NULL)) ;
	
	if(t->isTriangle)
	{
		if(this->isVertex(t->first) && this->isVertex(t->second))
			return std::pair< Point*,  Point*>(t->first , t->second) ;
		
		if(this->isVertex(t->second) && this->isVertex(t->third))
			return std::pair< Point*,  Point*>(t->second , t->third) ;
		
		if(this->isVertex(t->first) && this->isVertex(t->third))
			return std::pair< Point*,  Point*>(t->first , t->third) ;
	}
	else
	{
		if(isAligned((*t->first), (*first), (*second)) && isAligned((*t->second), (*first), (*second)))
			return std::pair< Point*,  Point*>(first, second ) ;
		if(isAligned((*t->first), (*third), (*second)) && isAligned((*t->second), (*third), (*second)))
			return std::pair< Point*,  Point*>(third, second ) ;
		if(isAligned((*t->first), (*third), (*first)) && isAligned((*t->second), (*third), (*first)))
			return std::pair< Point*,  Point*>(first, third ) ;
	}
	
	//assert(false) ;
	return std::pair< Point*,  Point*>(reinterpret_cast<Point*>(NULL), reinterpret_cast<Point*>(NULL)) ;
}

std::pair< Point*,  Point*> DelaunayDemiPlane::commonEdge(const DelaunayTreeItem * t) 
{
	if(t->isTriangle)
	{
		Point test = (*first)*0.5 + (*second)*0.5 ;
		
		if(isAligned(test, (*t->first), (*t->second)))
		{
#ifdef DEBUG
			std::cout << "\t -->(" <<  t->first->x << ", " << t->first->y << ") (" ;
			std::cout << t->second->x << ", " << t->second->y << ")" << std::endl;
#endif
			return std::pair< Point*,  Point*>(t->first, t->second ) ;
		}
		if(isAligned(test, (*t->third), (*t->second)))
		{
#ifdef DEBUG
			std::cout << "\t -->(" <<  t->third->x << ", " << t->third->y << ") (" ;
			std::cout << t->second->x << ", " << t->second->y << ")" << std::endl;
#endif
			return std::pair< Point*,  Point*>(t->third, t->second ) ;
		}
		if(isAligned(test, (*t->first), (*t->third)))
		{
#ifdef DEBUG
			std::cout << "\t -->(" <<  t->third->x << ", " << t->third->y << ") (" ;
			std::cout << t->first->x << ", " << t->first->y << ")" << std::endl;
#endif
			return std::pair< Point*,  Point*>( t->first, t->third ) ;
		}
	}
	else
	{
		if((t->first == first && t->second == second) ||
		   (t->first == second &&  t->second == first))
			return std::pair< Point*,  Point*>(t->first , t->second) ;
	}
	
	return std::pair< Point*,  Point*>(reinterpret_cast<Point*>(NULL), reinterpret_cast<Point*>(NULL)) ;
}
	

void DelaunayDemiPlane::merge(DelaunayDemiPlane * p)
{
	if(isAlive() && p != this && p->isAlive())
	{

		if( (first == p->first && second == p->second ) || (first == p->second && second == p->first) )
		{

			for(size_t i = 0 ; i <  p->neighbour.size() ; i++)
			{
				for(size_t j = 0 ; j <  neighbour.size() ; j++)
				{
					if(p->neighbour[i]->isNeighbour(neighbour[j]))
						makeNeighbours(neighbour[j], p->neighbour[i]) ;
				}
			}
			p->kill(first) ;
			kill(first) ;
			return ;
		}
		if(isAligned((*first), (*second), (*p->first)) && isAligned((*first), (*second), (*p->second)) &&
		   isAligned((*first), (*p->first) ,(*p->second)) && isAligned((*second), (*p->first) ,(*p->second)))
		{

			for(size_t i = 0 ; i <  p->neighbour.size() ; i++)
			{
				makeNeighbours(this, p->neighbour[i]) ;
			}
			for(size_t i = 0 ; i <  p->son.size() ; i++)
			{
				p->son[i]->father = this ;
				son.push_back(p->son[i]) ;
			}
			for(size_t i = 0 ; i <  p->stepson.size() ; i++)
			{
				p->stepson[i]->stepfather = this ;
				stepson.push_back(p->stepson[i]) ;
				makeNeighbours(this, p->stepson[i]) ;
			}
			p->father->son.push_back(this) ;
			p->son.push_back(this) ;
			p->kill(first) ;
		}
	}
}


std::pair< Point*,  Point*> DelaunayTriangle::nearestEdge(const Point p)
{
	std::map<double, Point> cen ;
	Point c0(((*first) + (*second))/2) ;
	cen[squareDist(c0, p)] = c0 ;
	Point c1(((*third) + (*second))/2) ;
	cen[squareDist(c1, p)] = c1 ;
	Point c2(((*third) + (*first))/2) ;
	cen[squareDist(c2, p)] = c2 ;
	
	if(cen.begin()->second == c0)
		return std::pair< Point*,  Point*>(first, second) ;
	if(cen.begin()->second == c1)
		return std::pair< Point*,  Point*>(third, second) ;
	if(cen.begin()->second == c2)
		return std::pair< Point*,  Point*>(third, first) ;
	
	return std::pair< Point*,  Point*>(first, second) ;
}
	
bool DelaunayTriangle::inCircumCircle(const Point *p) const 
{
	return  this->Triangle::inCircumCircle(p) ;
// 	return  (p->x - circumCenter.x)*(p->x - circumCenter.x) + (p->y - circumCenter.y)*(p->y - circumCenter.y)  < radius - 10*std::numeric_limits<double>::epsilon() ;
}
	
std::vector<DelaunayTreeItem *> * DelaunayTriangle::insert( Point *p,  Star* s)
{
	assert(first->id > -1) ;
	assert(second->id > -1) ;
	assert(third->id > -1) ;
	
#ifdef DEBUG
	std::cout << "inserting isVertex triangle" ; this->print() ;
#endif
	
	std::vector<DelaunayTreeItem *> * ret = new std::vector<DelaunayTreeItem *>(0) ;
	
	if(!isAlive())
		return ret ;
	
	if(isVertex(p))
		return ret ;
	
	if(visited)
		return ret ;
	
	if(!inCircumCircle(p))
		return ret ;
	
	if(isVertex(p))
		return ret ;
	
	assert(neighbour.size() <= 3) ;
	assert(son.size() <= 3) ;
	visited = true ;
	
	if(son.size() == 0)
	{
#ifdef DEBUG
		std::cout << "we are a leaf" << std::endl ;
#endif
		for (size_t i = 0 ; i < neighbour.size() ; i++)
		{
			assert( neighbour[i]->isNeighbour(this) );
			std::pair< Point*,  Point*> pp = this->commonEdge(neighbour[i]) ;

			if (!neighbour[i]->inCircumCircle(p))
			{
#ifdef DEBUG
				std::cout << "no conflict with neighbour" ; neighbour[i]->print() ;
#endif
				if(!isAligned((*p), (*pp.first), (*pp.second) ))
				{
#ifdef DEBUG
					std::cout << "we are not creating a degenerate triangle" << std::flush ;
#endif
					DelaunayTriangle *ss = new DelaunayTriangle(this, p, pp.first, pp.second, p) ;
					son.push_back(ss) ;
#ifdef DEBUG
					ss->print() ;
#endif
					neighbour[i]->addStepson(ss) ;
					ret->push_back(ss) ;
					
				}
				else
				{
#ifdef DEBUG
					std::cout << "degenerate case" << std::endl ;
#endif
				}
			}
#ifdef DEBUG
			else
			{
				std::cout << "there IS a conflict with neighbour " ; neighbour[i]->print() ;
			}
#endif
		}
		this->kill(p) ;
		s->updateNeighbourhood() ;
	}
#ifdef DEBUG
	std::cout << "done inserting" << std::endl ;
#endif
	return ret ;
}
	

bool DelaunayTriangle::isNeighbour(DelaunayTreeItem * t)  
{
	if(t == this)
		return false ;
	
	if(t->isTriangle )
	{
		return (commonEdge(t).first != NULL) ;
	}
	else
	{
		if (isAligned((*t->first),(*first), (*second) ) && isAligned((*t->second), (*first), (*second)))
		    return true ;
		if (isAligned((*t->first),(*third), (*second) ) && isAligned((*t->second), (*third), (*second)))
			return true ;
		if (isAligned((*t->first),(*third), (*first) ) && isAligned((*t->second), (*third), (*first))) 
			return true ;
	}
	
	return false ;
}

void DelaunayTriangle::print() const
{
	for(size_t i = 0 ; i < this->getBoundingPoints()->size() ; i++)
	{
		std::cout << "(" << getBoundingPoint(i)->x << ", " << getBoundingPoint(i)->y << ") " ;
	}
	std::cout <<  ":: "<< isAlive() << std::endl ;
}


// void DelaunayTriangle::displace(std::valarray<double> * eps)
// {
// 	(*eps)[first->id*2]+=(*eps)[first->id*2] ;
// 	(*eps)[first->id*2+1]+=(*eps)[first->id*2+1] ;
// }

DelaunayDemiPlane::DelaunayDemiPlane(DelaunayTreeItem * father,  Point  * _begin,  Point  * _end,  Point  * p,  Point * c) : DelaunayTreeItem(father, c)
{
	second  = _end ;
	first = _begin ;
	third = p ;
	dead = false ;
	vector =(*second)- (*first) ;
	Point pseudonormal = (*third) - (*first);
	direction = (vector.x*pseudonormal.y - vector.y*pseudonormal.x) ;
	isPlane =true ;
	isTriangle = false ;
}

std::pair< Point*,  Point*> DelaunayDemiPlane::nearestEdge(const Point p) 
{
	return std::pair< Point*,  Point*>(first, second) ;
}
	
bool DelaunayDemiPlane::inCircumCircle(const Point *p) const
{
	return ((vector.x*(p->y - first->y) - vector.y*(p->x - first->x)) * direction < -10*std::numeric_limits<double>::epsilon()) ;
}
	
bool DelaunayDemiPlane::isVertex(const Point *p) const
{
	return ( (*p) == (*first) || (*p) == (*second) ) || isAligned((*p), (*first), (*second)) ;
}
	

bool  DelaunayDemiPlane::isNeighbour(DelaunayTreeItem * t)
{
	
	if(t->isTriangle)
	{
		return t->isNeighbour(this) ;
	}
	else
	{
		return (first ==  t->first || first == t->first ||
		        second == t->first || second == t->second) ;
	}
return false ;
}

std::vector<DelaunayTreeItem *> * DelaunayDemiPlane::insert( Point *p, Star *s)
{
#ifdef DEBUG
	std::cout << "inserting isVertex the demi-plane " ; this->print() ;
#endif
	std::vector<DelaunayTreeItem *> * ret = new std::vector<DelaunayTreeItem *>(0) ;
	
	if (visited)
		return ret ;

	if(!isAlive())
		return ret ;
	
	if(!inCircumCircle(p))
		return ret ;
	
	if(isVertex(p))
		return ret ;
		
	visited = true ;
	
	if(son.size() == 0 )
	{
#ifdef DEBUG
 		std::cout << "we are a leaf" << std::endl ;
#endif
		std::vector<Point*> lims ;
		
		for(size_t i = 0 ; i < neighbour.size() ; i++)
		{
			std::pair< Point*,  Point*> pp = neighbour[i]->commonEdge(this) ;
			if(std::find(lims.begin(), lims.end(), pp.first) != lims.end())
				lims.erase(std::find(lims.begin(), lims.end(), pp.first)) ;
			else
				lims.push_back(pp.first) ;
			
			if(std::find(lims.begin(), lims.end(), pp.second) != lims.end())
				lims.erase(std::find(lims.begin(), lims.end(), pp.second)) ;
			else
				lims.push_back(pp.second) ;
			
			if (!neighbour[i]->inCircumCircle(p))
			{
#ifdef DEBUG
				std::cout << "there is no conflict with neighbour " ; neighbour[i]->print() ;
#endif
				assert(neighbour[i]->isNeighbour(this) );
				
				if(!isAligned((*p), (*first), (*second)))
				{
#ifdef DEBUG
					std::cout << "we are not creating a degenerate triangle" << std::flush ;
#endif
					DelaunayTriangle *ss = new DelaunayTriangle(this, p, pp.first, pp.second , p) ;
					ss->visited = true ;
#ifdef DEBUG
					ss->print() ;
#endif
					son.push_back(ss) ;
					neighbour[i]->addStepson(ss) ;
		
					ret->push_back(ss) ;
				}
			}
#ifdef DEBUG
			else if(neighbour[i]->inCircumCircle(p) ) 
			{
				std::cout << "there IS a conflict with neighbour " ; neighbour[i]->print() ;
				//neighbour[i]->insert(p, s) ;
			}
#endif
		}
		this->kill(p) ;

		for(size_t i = 0 ; i < deadneighbour.size() ; i++)
		{
			std::pair< Point*, Point*> pp = deadneighbour[i]->commonEdge(this) ;
			if(std::find(lims.begin(), lims.end(), pp.first) != lims.end())
				lims.erase(std::find(lims.begin(), lims.end(), pp.first)) ;
			else
				lims.push_back(pp.first) ;
			
			if(std::find(lims.begin(), lims.end(), pp.second) != lims.end())
				lims.erase(std::find(lims.begin(), lims.end(), pp.second)) ;
			else
				lims.push_back(pp.second) ;
		}
		if(lims.size() == 0)
		{
			lims.push_back(first) ;
			lims.push_back(second) ;
		}
		
		assert(lims.size() < 3) ;
		
		assert(!isAligned((*lims[0]), (*p), (*lims[1]))) ;
		assert((*lims[0]) != (*lims[1])) ;
		
		DelaunayDemiPlane *p0 = new DelaunayDemiPlane(this, lims[0], p, lims[1], p) ;
#ifdef DEBUG
		std::cout << "creating plane " << std::flush ; p0->print() ;
#endif
		DelaunayDemiPlane *p1 = new DelaunayDemiPlane(this, lims[1], p, lims[0], p) ;
#ifdef DEBUG
		std::cout << "creating plane " << std::flush ; p1->print() ;
#endif
		son.push_back(p0) ;
		son.push_back(p1) ;
		
		ret->push_back(p0) ;
		ret->push_back(p1) ;
		
		for(size_t i = 0 ; i < son.size()-2 ; i++)
		{
			if(son[i]->isNeighbour(p0))
				son[i]->addStepson(p0) ;
			
			if(son[i]->isNeighbour(p1))
				son[i]->addStepson(p1) ;
			
			for(size_t j = i ; j < son.size()-2 ; j++)
				if(son[i]->isNeighbour(son[j]))
					makeNeighbours(son[i], son[j]) ;
		}
		s->updateNeighbourhood() ;
	}
#ifdef DEBUG
	std::cout << "done inserting isVertex the plane" << std::endl ;
#endif
	return ret ;
}

void DelaunayDemiPlane::print() const 
{
	std::cout << "###############(" << first->x << ", " << first->y << ") (" <<
		second->x << ", " << second->y << ")" << "X (" << third->x << ", " << third->y << ") :: " << isAlive()  << std::endl ;
}

void makeNeighbours(DelaunayTreeItem *t0, DelaunayTreeItem *t1 ) 
{
	if(t0 == t1 || !t0->isAlive() || !t1->isAlive())
		return ;
	if(t0->isPlane && t1->isPlane)
	   return ;
	
	t0->addNeighbour(t1) ;
	t1->addNeighbour(t0) ;
} 

void updateNeighbours(std::vector<DelaunayTreeItem *> * t)
{
	for(size_t i = 0 ; i < t->size() ; i++)
	{
		for(size_t j = i ; j < t->size() ; j++)
		{
			if((*t)[i]->isNeighbour((*t)[j]))
			{
				makeNeighbours( (*t)[i],(*t)[j] ) ;
			}
		}
	}
} 

DelaunayRoot::DelaunayRoot(Point * p0, Point * p1, Point * p2) : DelaunayTreeItem(NULL, NULL)
{
	isPlane = false ;
	isTriangle =false ;
	this->father = NULL ;
	DelaunayTriangle *t = new DelaunayTriangle(this, p0, p1, p2, NULL) ;
	DelaunayDemiPlane * pl0 = new DelaunayDemiPlane(this, p0, p1, p2, NULL);
	DelaunayDemiPlane * pl1 = new DelaunayDemiPlane(this, p0, p2, p1, NULL);
	DelaunayDemiPlane * pl2 = new DelaunayDemiPlane(this, p1, p2, p0, NULL);
	
	makeNeighbours(t,pl0 ) ;
	makeNeighbours(t,pl1) ;
	makeNeighbours(t,pl2) ;

	son.push_back(t) ;
	son.push_back(pl0) ;
	son.push_back(pl1) ;
	son.push_back(pl2) ;
	kill(p0) ;
}
	
void DelaunayRoot::print() const
{
	std::cout << "I am root !" << std::endl ;
}


DelaunayTreeItem * DelaunayRoot::getSon(size_t i)
{
	return son[i] ;
}
	
bool DelaunayRoot::isVertex(const Point *p) const
{
	return false ;
}
	
bool DelaunayRoot::inCircumCircle(const Point *p) const 
{
	return true ;
}
	
std::pair< Point*,  Point*> DelaunayRoot::nearestEdge(const Point p)
{
	assert(false) ;
	return std::pair< Point*,  Point*>(reinterpret_cast<Point*>(NULL), reinterpret_cast<Point*>(NULL)) ;
}
	
std::vector<DelaunayTreeItem *> * DelaunayRoot::insert(Point *p, Star *s)
{
	std::vector<DelaunayTreeItem *> * ret = new std::vector<DelaunayTreeItem *>() ;
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		
		std::vector<DelaunayTreeItem *> * temp = son[i]->insert(p, s) ;
		
		for(size_t j = 0 ; j< temp->size() ; j++)
		{
			ret->push_back((*temp)[j]) ;
		}
		
		delete temp ;
	}
	
	updateNeighbours(ret) ;
	return ret ;
}



std::pair<std::vector<DelaunayTriangle *>, std::vector<DelaunayTreeItem *> > DelaunayRoot::conflicts(const Geometry *g)
{
	std::pair<std::vector<DelaunayTriangle *>, std::vector<DelaunayTreeItem *> > ret ;

	visited = true ;
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		std::pair<std::vector<DelaunayTriangle *>, std::vector<DelaunayTreeItem *> > temp  = son[i]->conflicts(g) ;
		ret.first.insert(ret.first.end(),temp.first.begin(), temp.first.end()) ;
		ret.second.insert(ret.second.end(),temp.second.begin(), temp.second.end()) ;
	}
	
	return ret ;
}


std::pair<std::vector<DelaunayTreeItem *>, std::vector<DelaunayTreeItem *> > * DelaunayRoot::conflicts(const Point *p)
{
	std::pair<std::vector<DelaunayTreeItem *>, std::vector<DelaunayTreeItem *> > * ret = new std::pair<std::vector<DelaunayTreeItem *>, std::vector<DelaunayTreeItem *> >  ;

	visited = true ;
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		std::pair<std::vector<DelaunayTreeItem *>, std::vector<DelaunayTreeItem *> > * temp  = son[i]->conflicts(p) ;
		ret->first.insert(ret->first.end(),temp->first.begin(), temp->first.end()) ;
		ret->second.insert(ret->second.end(),temp->second.begin(), temp->second.end()) ;
		delete temp ;
	}
	
	return ret ;
}

Star::Star(std::vector<DelaunayTreeItem *> *t, const Point *p)
{
	treeitem.insert(treeitem.end(), t->begin(), t->end()) ;
	for(size_t i = 0 ; i < t->size() ; i++)
	{
		if((*t)[i]->isTriangle)
		{
			this->edge.push_back((*t)[i]->first) ;
			this->edge.push_back((*t)[i]->second) ;
			this->edge.push_back((*t)[i]->third) ;
		}
		else if((*t)[i]->isPlane)
		{
			this->edge.push_back((*t)[i]->first) ;
			this->edge.push_back((*t)[i]->second) ;
		}
	}
	std::sort(edge.begin(), edge.end()) ;
	edge.erase(unique(edge.begin(), edge.end()), edge.end()) ;
}
	
size_t Star::size()
{
	return edge.size() ;
}
	
const Point * Star::getEdge(size_t i) const
{
	assert(i < edge.size()) ;
	return this->edge[i] ;
}


DelaunayTree::DelaunayTree(Point * p0, Point *p1, Point *p2)
{
	this->global_counter = 3;
	p0->id = 0 ; p1->id = 1 ; p2->id = 2 ;
	DelaunayRoot *root = new DelaunayRoot( p0, p1, p2) ;
	tree.push_back(root) ;
	tree.push_back(root->getSon(0)) ;
	tree.push_back(root->getSon(1)) ;
	tree.push_back(root->getSon(2)) ;
	tree.push_back(root->getSon(3)) ;
	plane.push_back(dynamic_cast<DelaunayDemiPlane *>(root->getSon(1))) ;
	plane.push_back(dynamic_cast<DelaunayDemiPlane *>(root->getSon(2))) ;
	plane.push_back(dynamic_cast<DelaunayDemiPlane *>(root->getSon(3))) ;
}
	
DelaunayTree::~DelaunayTree() 
{ 
	for(size_t i = 0 ;  i < this->tree.size() ; i++)
	{
// 		for(size_t j = 0 ; j < this->tree[i]->size() ; j++)
// 			delete this->tree[i]->getPoint(j) ;
// 		
		delete this->tree[i] ;
	}
} ;

void DelaunayTree::insertIf( Point *p, std::vector<SamplingCriterion *> v, double minScore )
{
	std::vector<DelaunayTreeItem *> * cons = new std::vector<DelaunayTreeItem *>;
	
	cons = this->conflicts(p) ;
	
	
	for(size_t i = 0 ; i < cons->size() ; i++)
	{
		(*cons)[i]->deadneighbour.clear() ;
	}
	
	//We store all the pointers to the affected elements of the tree, and make copies of those
	std::vector<DelaunayTreeItem *> backup = tree ;
	
	
#ifdef DEBUG
	for(size_t i = 0 ; i < cons->size() ; i++)
	{
		std::cout << "we have conflicts with" ;
		
		(*cons)[i]->print() ;
		
		for(size_t j = 0 ; j< (*cons)[i]->neighbour.size() ; j++)
		{
			std::cout << "\t --> " ;
			(*cons)[i]->neighbour[j]->print() ;
		}
	}
#endif
	
	Star * s = new Star(cons, p) ;
	
	std::vector<DelaunayTreeItem *> ret ;
	
	for(size_t i = 0 ; i < cons->size() ; i++)
	{
		std::vector<DelaunayTreeItem *> * temp = (*cons)[i]->insert(p, s) ;
		ret.insert(ret.end(), temp->begin(), temp->end()) ;
		delete temp ;
	}
	
	for(size_t i =  0 ; i < ret.size() ; i++)
	{
		if(ret[i]->isTriangle)
		{
			double score = 0;
			for(size_t j = 0 ; j < v.size() ; j++)
				if (v[j]->meetsCriterion(dynamic_cast<DelaunayTriangle *>(ret[i])))
					score += 1/v.size() ;
			
			if (score < minScore)
			{
				tree = backup ;	
				delete cons ;
				delete s ;
				for(size_t j =  0 ; j < ret.size() ; j++)
				{
					delete ret[j] ;
				}

				for(size_t j = 0 ; j < cons->size() ; j++)
					(*cons)[j]->clearVisited() ;
				
				return ;
			}
			else
			{
				p->id =this->global_counter++ ;
			}
		}
	}
	
	bool weGotPlanes = false ;
	
	for(size_t j = 0 ; j< ret.size() ; j++)
	{
		if(ret[j]->isAlive() && ret[j]->isPlane )
		{
			weGotPlanes = true ;
			
			if(ret[j]->neighbour.size() > 1 || ret[j]->neighbour.size() == 0)
			{
				ret[j]->kill(p) ;
			}
			else
			{
				plane.push_back(dynamic_cast<DelaunayDemiPlane*>(ret[j])) ;
			}
		}
	}
	
	if(weGotPlanes)
	{
		for(size_t k = 0 ; k< plane.size()-1 ; k++)
		{
			for(size_t j = k ; j< plane.size() ; j++)
			{
				plane[j]->merge(plane[k]) ;
			}
		}
	}
	
	for(size_t i = 0 ; i < ret.size() ; i++)
		ret[i]->clearVisited() ;

	for(size_t i = 0 ; i < cons->size() ; i++)
		(*cons)[i]->clearVisited() ;
	
	for(size_t i = 0 ; i < ret.size() ; i++)
	{
		
#ifdef DEBUG
		if(ret[i]->isAlive())
		{
			std::cout << "we the new guys..." ;
			
			ret[i]->print() ;
			
			for(size_t j = 0 ; j< ret[i]->neighbour.size() ; j++)
			{
				std::cout << "\t --> " ;
				ret[i]->neighbour[j]->print() ;
			}
		}
#endif
		assert(ret[i]->isPlane  || ret[i]->neighbour.size() == 3) ;
		
		if(ret[i]->isAlive())
			tree.push_back(ret[i]) ;
	}
	
	std::vector<DelaunayDemiPlane *> * hull = this->getConvexHull() ;
	plane.clear() ;
	plane.insert(plane.end(), hull->begin(), hull->end()) ;
	delete hull ;
	delete cons ;
	delete s ;
}

void DelaunayTree::insert(Segment *s)
{
	std::vector<Point *> * cons = this->conflicts(s);
	
	this->insert(s->first()) ;
	
	for(size_t i = 0 ;  i < cons->size() ; i++)
			this->insert((*cons)[i]) ;
		
	this->insert(s->second()) ;
	
	delete cons ;
	
}

void DelaunayTree::insert(Point *p)
{
	std::vector<DelaunayTreeItem *> * cons = this->conflicts(p) ;
	if(cons->size() == 0)
	{
		return ;
	}
	p->id = this->global_counter++ ;
	
	for(size_t i = 0 ; i < cons->size() ; i++)
	{
		(*cons)[i]->deadneighbour.clear() ;
	}
	
#ifdef DEBUG
	for(size_t i = 0 ; i < cons->size() ; i++)
	{
		std::cout << "we have conflicts with" ;
		
		(*cons)[i]->print() ;
		
		for(size_t j = 0 ; j< (*cons)[i]->neighbour.size() ; j++)
		{
			std::cout << "\t --> " ;
			(*cons)[i]->neighbour[j]->print() ;
		}
	}
#endif
	Star * s = new Star(cons, p) ;
	
	std::vector<DelaunayTreeItem *> ret ;
	
	for(size_t i = 0 ; i < cons->size() ; i++)
	{
		std::vector<DelaunayTreeItem *> * temp = (*cons)[i]->insert(p, s) ;
		ret.insert(ret.end(), temp->begin(), temp->end()) ;
		delete temp ;
	}
	
	bool weGotPlanes = false ;
	
	for(size_t j = 0 ; j< ret.size() ; j++)
	{
		if(ret[j]->isAlive() && ret[j]->isPlane )
		{
			weGotPlanes = true ;
			
			if(ret[j]->neighbour.size() > 1 || ret[j]->neighbour.size() == 0)
			{
				ret[j]->kill(p) ;
			}
			else
			{
				plane.push_back(dynamic_cast<DelaunayDemiPlane*>(ret[j])) ;
			}
		}
	}
	
	if(weGotPlanes)
	{
		for(size_t k = 0 ; k< plane.size()-1 ; k++)
		{
			for(size_t j = k ; j< plane.size() ; j++)
			{
				plane[j]->merge(plane[k]) ;
			}
		}
	}

	for(size_t i = 0 ; i < ret.size() ; i++)
	{
		ret[i]->clearVisited() ;
// 			if(ret[i]->isAlive() && ret[i]->isTriangle && ret[i]->neighbour.size() != 3)
// 				ret[i]->kill(p);	
	}
	for(size_t i = 0 ; i < cons->size() ; i++)
	{
		(*cons)[i]->clearVisited() ;
	}
	for(size_t i = 0 ; i < ret.size() ; i++)
	{
		
#ifdef DEBUG
		if(ret[i]->isAlive())
		{
			std::cout << "we the new guys..." ;
			
			ret[i]->print() ;
			
			for(size_t j = 0 ; j< ret[i]->neighbour.size() ; j++)
			{
				std::cout << "\t --> " ;
				ret[i]->neighbour[j]->print() ;
			}
		}
#endif
		assert(ret[i]->isPlane  || ret[i]->neighbour.size() == 3) ;
		
		if(ret[i]->isAlive())
			tree.push_back(ret[i]) ;
		
#ifdef NDEBUG
		if(ret[i]->isTriangle)
			assert(dynamic_cast<DelaunayTriangle*>(ret[i])->Triangle::area() >0 ) ;
#endif
	}

	std::vector<DelaunayDemiPlane *> * hull = this->getConvexHull() ;
	plane.clear() ;
	plane.insert(plane.end(), hull->begin(), hull->end()) ;
	delete hull ;
	delete cons ;
	delete s ;
}
	

std::vector<DelaunayTreeItem *> * DelaunayTree::slowConflicts( const Point *p) const
{
	std::vector<DelaunayTreeItem *> * ret = new std::vector<DelaunayTreeItem *> ;
	
	for(size_t i = 1 ; i < tree.size() ; i++)
	{
		if(tree[i]->isAlive() && tree[i]->inCircumCircle(p))
			ret->push_back(tree[i]) ;
	}
	
	return ret ;
}

std::vector<Point *> * DelaunayTree::slowConflicts( const Segment *s) const
{
	std::vector<Point *> * ret = new std::vector<Point *> ;
	
	for(size_t i = 1 ; i < tree.size() ; i++)
	{
		if(tree[i]->isAlive())
		{
			if(tree[i]->isTriangle)
			{
				Segment s0((*tree[i]->first), (*tree[i]->second)) ;
				Segment s1((*tree[i]->second), (*tree[i]->third)) ;
				Segment s2((*tree[i]->third), (*tree[i]->first)) ;
			
				if (s->intersects(&s0) )
				{
					Point * i = new Point(s->intersection(&s0)) ;
		
					bool unique = true ;
					for(size_t k  = 0 ; k < ret->size() && unique == true ;k++)
					{
						if((*(*ret)[k]) == (*i))
							unique = false ; 
					}
						
					if(unique)
						ret->push_back(i) ;
					else
						delete i ;
				}
				if (s->intersects(&s1) )
				{
					Point * i = new Point(s->intersection(&s1)) ;
					
					bool unique = true ;
					for(size_t k  = 0 ; k < ret->size() && unique == true ;k++)
					{
						if((*(*ret)[k]) == (*i))
							unique = false ; 
					}
					
					if(unique)
						ret->push_back(i) ;
					else
						delete i ;
				}
				if (s->intersects(&s2) )
				{
					Point * i = new Point(s->intersection(&s2)) ;
					
					bool unique = true ;
					for(size_t k  = 0 ; k < ret->size() && unique == true ;k++)
					{
						if((*(*ret)[k]) == (*i))
							unique = false ; 
					}
					
					if(unique)
						ret->push_back(i) ;
					else
						delete i ;
				}
				
			}
			else
			{
				Segment s0((*tree[i]->first), (*tree[i]->second)) ;
				
				if (s->intersects(&s0) )
				{
					Point * i = new Point(s->intersection(&s0)) ;
					
					bool unique = true ;
					for(size_t k  = 0 ; k < ret->size() && unique == true ;k++)
					{
						if((*(*ret)[k]) == (*i))
							unique = false ; 
					}
					
					if(unique)
						ret->push_back(i) ;
					else
						delete i ;
				}
			}
		}
	}
	
	return ret ;
}

std::vector<DelaunayTriangle *> * DelaunayTree::slowConflicts( const Geometry *g) const
{
	std::vector<DelaunayTriangle *> * ret = new std::vector<DelaunayTriangle *> ;
	
	for(size_t i = 1 ; i < tree.size() ; i++)
	{
		if(tree[i]->isAlive() && 
		   tree[i]->isTriangle && 
		   ( 
		     g->in(tree[i]->first) || 
		     g->in(tree[i]->second) || 
		     g->in(tree[i]->third) 
		   )
		  )
			ret->push_back(dynamic_cast<DelaunayTriangle *>(tree[i])) ;
	}
	
	return ret ;
}

std::vector<DelaunayTreeItem *> * DelaunayTree::conflicts( const Point *p) const
{
	std::pair< std::vector<DelaunayTreeItem *>,std::vector<DelaunayTreeItem *> > * cons = this->tree[0]->conflicts(p) ;
	
	for(size_t i = 0 ; i < plane.size() ; i++)
	{
	
		if(!plane[i]->visited)
		{
			std::pair< std::vector<DelaunayTreeItem *>,std::vector<DelaunayTreeItem *> > * temp = plane[i]->conflicts(p) ;
			
			cons->first.insert(cons->first.end(), temp->first.begin(),temp->first.end()) ;
			cons->second.insert(cons->second.end(), temp->second.begin(),temp->second.end()) ;
			delete temp ;
			
			
			for(size_t j = 0 ; j < plane[i]->neighbour.size() ; j++)
			{
				std::pair< std::vector<DelaunayTreeItem *>, std::vector<DelaunayTreeItem *> > * temp = plane[i]->neighbour[j]->conflicts(p) ;
				
				cons->first.insert(cons->first.end(), temp->first.begin(),temp->first.end()) ;
				cons->second.insert(cons->second.end(), temp->second.begin(),temp->second.end()) ;
				delete temp ;
			}
			
		}
		
	}
	
	for(size_t i = 0 ; i < cons->second.size() ; i++)
	{
		cons->second[i]->clearVisited() ;
	}
	
	std::vector<DelaunayTreeItem *> * ret = new std::vector<DelaunayTreeItem *> ;
	ret->insert(ret->end(), cons->first.begin(), cons->first.end()) ;
	delete cons ;
	if(ret->size() == 0)
	{
		delete ret ;
		std::vector<DelaunayTreeItem *> * ret = slowConflicts(p) ;
		return ret ;
	}
	
	return ret ;
}

std::vector<Point *> * DelaunayTree::conflicts( const Segment *s) const
{
	std::pair< std::vector<Point *>,std::vector<DelaunayTreeItem *> > * cons = new std::pair< std::vector<Point *>,std::vector<DelaunayTreeItem *> >;
	
	cons = this->tree[0]->conflicts(s) ;
	
	for(size_t i = 0 ; i < plane.size() ; i++)
	{
		
		if(!plane[i]->visited)
		{
			std::pair< std::vector<Point *>,std::vector<DelaunayTreeItem *> > * temp = plane[i]->conflicts(s) ;
			
			cons->first.insert(cons->first.end(), temp->first.begin(),temp->first.end()) ;
			cons->second.insert(cons->second.end(), temp->second.begin(),temp->second.end()) ;
			delete temp ;
			
			
			for(size_t j = 0 ; j < plane[i]->neighbour.size() ; j++)
			{
				std::pair< std::vector<Point *>, std::vector<DelaunayTreeItem *> > * temp = plane[i]->neighbour[j]->conflicts(s) ;
				
				cons->first.insert(cons->first.end(), temp->first.begin(),temp->first.end()) ;
				cons->second.insert(cons->second.end(), temp->second.begin(),temp->second.end()) ;
				delete temp ;
			}
			
		}
		
	}
	
	for(size_t i = 0 ; i < cons->second.size() ; i++)
	{
		cons->second[i]->clearVisited() ;
	}
	
	std::vector<Point *> * ret = new std::vector<Point *> ;
	ret->insert(ret->end(), cons->first.begin(), cons->first.end()) ;
	delete cons ;
	if(ret->size() == 0)
	{
		delete ret ;
		//std::cout << "this is wrong... Oh, well, applying slow but sure method." << std::endl ;
		std::vector<Point *> * ret = slowConflicts(s) ;
		return ret ;
	}
	
	return ret ;
}

std::vector<DelaunayTriangle *> * DelaunayTree::conflicts(const Geometry *g) const
{
	std::pair< std::vector<DelaunayTriangle *>,std::vector<DelaunayTreeItem *> > cons = this->tree[0]->conflicts(g) ;
	
	for(size_t i = 0 ; i < plane.size() ; i++)
	{
		
		if(!plane[i]->visited)
		{
			std::pair< std::vector<DelaunayTriangle *>,std::vector<DelaunayTreeItem *> > temp = plane[i]->conflicts(g) ;
			
			cons.first.insert(cons.first.end(), temp.first.begin(),temp.first.end()) ;
			cons.second.insert(cons.second.end(), temp.second.begin(),temp.second.end()) ;
			
			
			for(size_t j = 0 ; j < plane[i]->neighbour.size() ; j++)
			{
				std::pair< std::vector<DelaunayTriangle *>, std::vector<DelaunayTreeItem *> > temp = plane[i]->neighbour[j]->conflicts(g) ;
				
				cons.first.insert(cons.first.end(), temp.first.begin(),temp.first.end()) ;
				cons.second.insert(cons.second.end(), temp.second.begin(),temp.second.end()) ;
			}
			
		}
		
	}
	
	for(size_t i = 0 ; i < cons.second.size() ; i++)
	{
		cons.second[i]->clearVisited() ;
	}
	
	std::vector<DelaunayTriangle *> * ret = new std::vector<DelaunayTriangle *> ;
	ret->insert(ret->end(), cons.first.begin(), cons.first.end()) ;
	
	
	if(ret->size() == 0)
	{
		delete ret ;
		//std::cout << "this is wrong... Oh, well, applying slow but sure method." << std::endl ;
		std::vector<DelaunayTriangle *> * ret = slowConflicts(g) ;
		return ret ;
	}
	
	return ret ;
}

std::vector<DelaunayDemiPlane *> * DelaunayTree::getConvexHull()
{
	std::vector<DelaunayDemiPlane *> * ret = new std::vector<DelaunayDemiPlane *> ;
	
	for(size_t i = 0 ; i < plane.size() ; i++)
	{
		if(plane[i]->isAlive())
			ret->push_back(plane[i]) ;
	}
	
	return ret ;
}

std::vector<DelaunayTriangle *> * DelaunayTree::getTriangles() const
{
	std::vector<DelaunayTriangle *> * ret = new std::vector<DelaunayTriangle *> ;
	
	for(size_t i = 0 ; i < tree.size() ; i++)
	{
		if(tree[i]->isAlive() && tree[i]->isTriangle)
		{
			ret->push_back(dynamic_cast<DelaunayTriangle *>(tree[i])) ;
		}
	}
	
	return ret ;
}

void DelaunayTree::print() const
{
	size_t alive = 0 ;
	std::cout << "we have a total of " << tree.size() << " elements" << std::endl ;
	for(size_t i = 0 ; i < tree.size() ; i++)
	{
		if (tree[i]->isAlive())
		{
			alive++ ;
		}
	}
	std::cout << "of which " << alive << "are alive" << std::endl ;
	#ifdef DEBUG
	for(size_t i = 0 ; i < tree.size() ; i++)
	{
		if (tree[i]->isAlive() && tree[i]->neighbour.size() > 3 )
		{
			tree[i]->print() ;
			for(size_t j = 0 ; j< tree[i]->neighbour.size() ; j++)
			{
				std::cout << "\t --> " ;
				tree[i]->neighbour[j]->print() ;
			}
		}
	}
	#endif
}



