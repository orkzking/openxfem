#include "geometry_2D.h"

#ifndef M_PI 
#define M_PI  3.14159265358979323846
#endif

double round(double x){
  return floor(x+0.5);
}

using namespace Mu ;

Triangle::Triangle() : ConvexGeometry(3)
{
	gType = TRIANGLE ;
	assert(this->size() == 3) ;
	
	(*boundingPoints)[0] = new Point(0, 1) ;
	(*boundingPoints)[1] = new Point(0, 0) ;
	(*boundingPoints)[2] = new Point(1, 0) ;
	
	computeCenter() ;
}

Triangle::Triangle( Point p0,  Point p1,  Point p2) : ConvexGeometry(3)
{
	gType = TRIANGLE ;
	
	assert(this->size() == 3) ;
	
	(*boundingPoints)[0] = new Point(p0) ;
	(*boundingPoints)[1] = new Point(p1) ;
	(*boundingPoints)[2] = new Point(p2) ;
	
	
	if(p0.z == p1.z == p2.z == 0)
	{
		if(!isTrigoOriented())
		{
			std::swap((*boundingPoints)[1], (*boundingPoints)[2]) ;
		}
	}
	
	computeCircumCenter() ;
	computeCenter() ;
	
	if(p0.z == p1.z == p2.z == 0)
	{
		if(!this->in(this->getCenter()))
		{
			assert(false) ;
		}
	}
	
	radius = (squareDist(p1, circumCenter) + squareDist(p0, circumCenter) + squareDist(p2, circumCenter))/3.;
	
}

Triangle::Triangle( Point *p0,  Point *p1,  Point *p2): ConvexGeometry(3)
{
	gType = TRIANGLE ;
	assert(this->size() == 3) ;
	

	(*boundingPoints)[0] = p0 ;
	(*boundingPoints)[1] = p1 ;
	(*boundingPoints)[2] = p2 ;

	
	if(p0->z == p1->z == p2->z == 0)
	{
		if(!isTrigoOriented())
		{
			std::swap((*boundingPoints)[1], (*boundingPoints)[2]) ;
		}
	}
	
	
	computeCircumCenter() ;
	computeCenter() ;
	
	if(p0->z == p1->z == p2->z == 0)
	{
		if(!this->in(this->getCenter()))
		{
			assert(false) ;
		}
		
	}
	
	radius = (squareDist(*p1, circumCenter) + squareDist(*p0, circumCenter) + squareDist(*p2, circumCenter))/3.;
}


void Triangle::computeCenter()
{
	for(size_t i = 0 ; i < this->size() ; i++)
	{
		this->center += *this->getPoint(i) ;
	}
	
	this->center = this->center/this->size() ;
}

double Triangle::getRadius() const
{
	return sqrt(radius) ;
}

Point Triangle::getCircumCenter() const
{
	return this->circumCenter ;
}

void Triangle::computeCircumCenter()
{	
	if (fabs((*boundingPoints)[1]->y-(*boundingPoints)[0]->y) < 20*std::numeric_limits<double>::epsilon()) 
	{
		double m2 = - ((*boundingPoints)[2]->x-(*boundingPoints)[1]->x) / ((*boundingPoints)[2]->y-(*boundingPoints)[1]->y);
		double mx2 = ((*boundingPoints)[1]->x + (*boundingPoints)[2]->x) / 2.0;
		double my2 = ((*boundingPoints)[1]->y + (*boundingPoints)[2]->y) / 2.0;
		double xc = ((*boundingPoints)[1]->x + (*boundingPoints)[0]->x) / 2.0;
		double yc = m2 * (xc - mx2) + my2;
		
		circumCenter = Point(xc, yc) ;
	} 
	else if (fabs((*boundingPoints)[2]->y-(*boundingPoints)[1]->y) < 20*std::numeric_limits<double>::epsilon()) 
	{
		double m1 = - ((*boundingPoints)[1]->x-(*boundingPoints)[0]->x) / ((*boundingPoints)[1]->y-(*boundingPoints)[0]->y);
		double mx1 = ((*boundingPoints)[0]->x + (*boundingPoints)[1]->x) / 2.0;
		double my1 = ((*boundingPoints)[0]->y + (*boundingPoints)[1]->y) / 2.0;
		double xc = ((*boundingPoints)[2]->x + (*boundingPoints)[1]->x) / 2.0;
		double yc = m1 * (xc - mx1) + my1;
		
		circumCenter = Point(xc, yc) ;
	} 
	else 
	{
		double m1 = - ((*boundingPoints)[1]->x-(*boundingPoints)[0]->x) / ((*boundingPoints)[1]->y-(*boundingPoints)[0]->y);
		double m2 = - ((*boundingPoints)[2]->x-(*boundingPoints)[1]->x) / ((*boundingPoints)[2]->y-(*boundingPoints)[1]->y);
		double mx1 = ((*boundingPoints)[0]->x + (*boundingPoints)[1]->x) / 2.0;
		double mx2 = ((*boundingPoints)[1]->x + (*boundingPoints)[2]->x) / 2.0;
		double my1 = ((*boundingPoints)[0]->y + (*boundingPoints)[1]->y) / 2.0;
		double my2 = ((*boundingPoints)[1]->y + (*boundingPoints)[2]->y) / 2.0;
		double xc = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
		double yc = m1 * (xc - mx1) + my1;
		
		circumCenter = Point(xc, yc) ;
	}
}

bool Triangle::inCircumCircle(const Point p) const
{
	return  squareDist(circumCenter, p) < radius - 10*std::numeric_limits<double>::epsilon() ;
}

bool Triangle::inCircumCircle(const Point *p) const
{
	return  squareDist(circumCenter, (*p)) < radius - 10*std::numeric_limits<double>::epsilon() ;
}

bool Triangle::is1D() const
{ 
	return false ; 
} 

double Triangle::area() const
{
	assert(this->boundingPoints->size() == 3) ;
	
	Segment s0((*(*boundingPoints)[0]), (*(*boundingPoints)[1])) ;
	Segment s1((*(*boundingPoints)[0]), (*(*boundingPoints)[2])) ;

	return 0.5*((*s0.vector())^(*s1.vector())).norm() ;
}


void Triangle::project(Point * p) const
{
	Segment s(*getCenter(), *p) ;
	
	for(size_t i = 0 ; i < boundingPoints->size() ; i++)
	{
		if(s.on(getPoint(i)) && !Segment(*getPoint(i), *p).on(getCenter()))
		{
			p->x = getPoint(i)->x ;
			p->y = getPoint(i)->y ;
			return ;
		}
	}
	
	for(size_t i = 0 ; i < boundingPoints->size() ; i++)
	{
		Segment seg(*getPoint(i), *getPoint((i+1)%boundingPoints->size())) ;
		if(s.intersects(&seg))
		{
			Point t = s.intersection(&seg);
			p->x = t.x ;
			p->y = t.y ;
			return ;
		}
	}
	
}

bool Triangle::in(const Point p) 
{
	return ConvexPolygon::in(p) ;
}

bool Triangle::in(const Point *p) const
{
	return ConvexPolygon::in(p) ;
}

bool Triangle::in(const Point p) const
{
	return ConvexPolygon::in(p) ;
}

void Triangle::sampleBoundingSurface(size_t num_points)
{
	assert(num_points%3 == 0) ;
	
	Point v0(*(*this->boundingPoints)[0]) ;
	Point v1(*(*this->boundingPoints)[1]) ;
	Point v2(*(*this->boundingPoints)[2]) ;
	
	this->boundingPoints->resize(num_points) ;
	
	(*this->boundingPoints)[0] = new Point(v0) ;
	
	for(size_t i = 1 ; i < num_points/3 ; i++)
	{
		(*this->boundingPoints)[i] = new Point(v0*3.*i/num_points + v1*(1.-3.*i/num_points)) ;
	}
	
	(*this->boundingPoints)[num_points/3] = new Point(v1) ;
	
	for(size_t i = num_points/3+1 ; i < 2*num_points/3 ; i++)
	{
		(*this->boundingPoints)[i] = new Point(v1*3.*(i-num_points/3)/num_points + v2*(1.-3.*(i-num_points/3)/num_points)) ;
	}
	
	(*this->boundingPoints)[2*num_points/3] = new Point(v2) ;
	
	for(size_t i = 2*num_points/3+1 ; i < num_points ; i++)
	{
		(*this->boundingPoints)[i] = new Point(v2*3.*(i-2*num_points/3)/num_points + v0*(1.-3.*(i-2*num_points/3)/num_points)) ;
	}
}

void Triangle::sampleSurface(size_t num_points)
{
	if(this->size() == 3)
		this->sampleBoundingSurface(num_points) ;
}

Rectangle::Rectangle(double x, double y, double originX, double originY) : ConvexGeometry(4), size_y(y), size_x(x)
{
	gType = RECTANGLE ;
	this->center = Point(originX,originY) ;
	(*boundingPoints)[0] = new Point(originX-0.5*x, originY+0.5*y) ;
	(*boundingPoints)[1] = new Point(originX-0.5*x, originY-0.5*y) ;
	(*boundingPoints)[2] = new Point(originX+0.5*x, originY-0.5*y) ;
	(*boundingPoints)[3] = new Point(originX+0.5*x, originY+0.5*y) ;
} 

Rectangle::Rectangle(double x, double y, Point &center) : ConvexGeometry(4), size_y(y), size_x(x)
{
	gType = RECTANGLE ;
	this->center = center ;
	(*boundingPoints)[0] = new Point(center.x-0.5*x, center.y+0.5*y) ;
	(*boundingPoints)[1] = new Point(center.x-0.5*x, center.y-0.5*y) ;
	(*boundingPoints)[2] = new Point(center.x+0.5*x, center.y-0.5*y) ;
	(*boundingPoints)[3] = new Point(center.x+0.5*x, center.y+0.5*y) ;
}

Rectangle::Rectangle() : ConvexGeometry(4), size_y(2), size_x(2)
{
	gType = RECTANGLE ;
	this->center = Point(0,0) ;
	(*boundingPoints)[0] = new Point(-1, 1) ;
	(*boundingPoints)[1] = new Point(-1, -1) ;
	(*boundingPoints)[2] = new Point(1, -1) ;
	(*boundingPoints)[3] = new Point(1,1) ;
}


void Rectangle::computeCenter()
{
	for(size_t i = 0 ; i < this->size() ; i++)
		this->center += *this->getPoint(i) ;
	
	this->center = this->center/this->size() ;
}

double  Rectangle::getRadius() const
{
	return sqrt(width()*width()*0.25 + height()*height()*0.25) ;
}

bool Rectangle::in(const Point p)
{
	if(p.x < getCenter()->x - 0.5*size_x)
		return false ;
	if(p.x > getCenter()->x + 0.5*size_x)
		return false ;
	if(p.y > getCenter()->y + 0.5*size_y)
		return false ;
	if(p.y < getCenter()->y - 0.5*size_y)
		return false ;
	
	return true ;
}
bool Rectangle::in(const Point * p) const
{
	if(p->x < getCenter()->x - 0.5*size_x)
		return false ;
	if(p->x > getCenter()->x + 0.5*size_x)
		return false ;
	if(p->y > getCenter()->y + 0.5*size_y)
		return false ;
	if(p->y < getCenter()->y - 0.5*size_y)
		return false ;
	
	return true ;
}


double Rectangle::area() const
{
	return this->size_x*this->size_y ;
}

bool Rectangle::is1D() const
{ 
	return false ; 
} 

double Rectangle::width() const
{
	return size_x ;
}

double Rectangle::height() const
{
	return size_y ;
}


void Rectangle::project(Point * p) const
{
	
	Line l(*getCenter(), *p) ;
	
	for(size_t i = 0 ; i < boundingPoints->size() ; i++)
	{
		Segment test(*getBoundingPoint(i), *getBoundingPoint((i+1)%boundingPoints->size())) ;
		
		if(l.intersects(&test))
		{
			Point inter = l.intersection(&test) ;
			p->x = inter.x ;
			p->y = inter.y ;
			return ;
		}
	}
}

void Rectangle::sampleBoundingSurface(size_t num_points)
{
	assert(num_points%4 == 0) ;
	double perimeter = 2*(size_x+size_y) ;
	
	double distanceBetweenPoints = perimeter/num_points ;
	
	this->numberOfPointsAlongX = static_cast<size_t>(std::ceil(size_x/distanceBetweenPoints) + 1);
	double distanceBetweenPointsAlongX = size_x/(this->numberOfPointsAlongX-1) ;
	
	this->numberOfPointsAlongY = static_cast<size_t>(std::ceil(size_y/distanceBetweenPoints) + 1);
	double distanceBetweenPointsAlongY = size_y/(this->numberOfPointsAlongY-1) ;
	
	num_points = ((numberOfPointsAlongX)*2 + (numberOfPointsAlongY)*2 - 4) ;
	
	boundingPoints->resize(num_points) ;
	
	for (size_t i = 0 ; i < numberOfPointsAlongY; i++)
	{
		(*boundingPoints)[i] = new Point(center.x-0.5*size_x, center.y + 0.5*size_y - i*distanceBetweenPointsAlongY) ;
	}
	for (size_t i = 1 ; i < numberOfPointsAlongX ; i++)
	{
		(*boundingPoints)[numberOfPointsAlongY+i-1] = new Point( center.x-0.5*size_x+i*distanceBetweenPointsAlongX, 
			getCenter()->y-0.5*size_y);
	}
	for (size_t i = 1 ; i < numberOfPointsAlongY ; i++)
	{
		(*boundingPoints)[numberOfPointsAlongX+numberOfPointsAlongY+i-2] = new Point(center.x+0.5*size_x,
			center.y-0.5*size_y+i*distanceBetweenPointsAlongY);
	}
	for (size_t i = 1 ; i < numberOfPointsAlongX-1 ; i++)
	{
		assert(2*numberOfPointsAlongY+numberOfPointsAlongX+i-3< num_points) ;
		(*boundingPoints)[2*numberOfPointsAlongY+numberOfPointsAlongX+i-3] = new Point(
			center.x + 0.5*size_x - i*distanceBetweenPointsAlongX ,
			center.y + 0.5*size_y) ;
		
	}
}

void Rectangle::sampleSurface(size_t num_points)
{
	if(this->size() == 4)
	{
		this->Rectangle::sampleBoundingSurface(num_points) ;
	}
	
	size_t nip = static_cast<size_t>((numberOfPointsAlongX-2)*(numberOfPointsAlongY-2)) ;
	
	inPoints->resize(nip) ;
	
	if(nip > 0)
	{
		for(size_t i = 0 ; i < this->numberOfPointsAlongX-2 ; i++)
		{
			for(size_t j = 0 ; j < this->numberOfPointsAlongY-2 ; j++)
			{	
				(*inPoints)[i*(numberOfPointsAlongY-2)+j] = new Point((*this->boundingPoints)[numberOfPointsAlongY+i]->x, (*this->boundingPoints)[j+1]->y) ;
			}
		}
	}
}

Circle::Circle(double r, double originX, double originY) : ConvexGeometry(0)
{
	gType = CIRCLE ;
	this->center = Point(originX, originY) ;
	this->radius = r ;
}

Circle::Circle(double r, const Point *center) : ConvexGeometry(0)
{
	gType = CIRCLE ;
	this->center = Point(*center) ;
	this->radius = r ; 
}

Circle::Circle(double r, const Point center) : ConvexGeometry(0)
{
	gType = CIRCLE ;
	this->center = center ;
	this->radius = r ; 
}

void Circle::computeCenter()
{
}


bool Circle::is1D() const
{
	return false ;
} 

void Circle::project(Point * p) const
{
	double hypothenuse = dist(*p, center) ;
	double oppose = std::abs(p->y - center.y) ;
	double adjacent = std::abs(p->x - center.x) ;
	
	if(p->x >= center.x && p->y >= center.y)
	{
		p->x = radius * adjacent/hypothenuse + center.x ;
		p->y = radius * oppose/hypothenuse + center.y ;
	}
	else if (p->x < center.x && p->y >= center.y)
	{
		p->x = - radius * adjacent/hypothenuse + center.x ;
		p->y = radius * oppose/hypothenuse + center.y ;
	}
	else if (p->x < center.x && p->y < center.y)
	{
		p->x = - radius * adjacent/hypothenuse + center.x ;
		p->y = - radius * oppose/hypothenuse + center.y ;
	}
	else
	{
		p->x = radius * adjacent/hypothenuse + center.x ;
		p->y = - radius * oppose/hypothenuse + center.y ;
	}
	
	return ;
}

void Circle::sampleBoundingSurface(size_t num_points)
{
	this->boundingPoints->resize(num_points) ;
	
	double angle = 2*M_PI/ (num_points) ;
	
	for (size_t i = 0 ; i< num_points ; i++)
	{
		(*this->boundingPoints)[i] = new Point(radius*cos((double)i*angle) + getCenter()->x, radius*sin((double)i*angle) + getCenter()->y);
// 		std::cout << "x = " << (*this->boundingPoints)[i]->x << ", y = " << (*this->boundingPoints)[i]->y << std::endl ;
	}
}

void Circle::sampleSurface(size_t num_points)
{
	assert(!sampled) ;
	
	if(boundingPoints->size() == 0)
		this->sampleBoundingSurface(num_points) ;
	sampled = true ;
	
	size_t numberOfRings = static_cast<size_t>((double)num_points/(2 * M_PI )) ;
	if(numberOfRings > 0)
		numberOfRings-- ;
	assert(numberOfRings >= 0) ;
	double angle = 2*M_PI/ (num_points) ;
	double offset = 0 ;
	
	//std::cout << "we have " << numberOfRings<< " rings" << std::endl ;
	
	std::vector<Point*> temp ;
	
	for (size_t i = 0 ; i< numberOfRings ; ++i)
	{
		double r = radius*(1. - (double)(pow(i, 1.1) + 1)/(pow(numberOfRings, 1.1)+1)) ;
		//std::cout << "radius is " << r << std::endl ;
		
		for (size_t j = 0 ; j< num_points ; ++j)
		{
			temp.push_back(new Point(r*cos((double)(j+0.5*(i))*angle+offset) + getCenter()->x, r*sin((double)(j+0.5*(i))*angle) + getCenter()->y));
		}
		
// 		num_points/=1.1 ;
// 		
// 		angle = 2*M_PI/ (num_points) ;
// 		
// 		offset = 0.5*(2*M_PI/ (num_points) -  2*M_PI/ (num_points*1.1)) ;
	}
	
	inPoints->resize(temp.size() + 1) ;
	(*inPoints)[0] = new Point(center) ;
	std::copy(temp.begin(), temp.end(),&(*inPoints)[1]) ;
	//std::cout << "we have " << num_points << " sample points" << std::endl ;
	
}

bool Circle::in(const Point v) 
{
	return (squareDist(center, v) <= radius*radius + std::numeric_limits<double>::epsilon()) ;
}

bool Circle::in(const Point * v) const
{
	return (squareDist(center, (*v)) <= radius*radius + std::numeric_limits<double>::epsilon()) ;
}

double Circle::getRadius() const
{
	return radius ;
}

double Circle::area() const
{
	return M_PI*this->radius*this->radius ;
}


SegmentedLine::SegmentedLine(PointSet * points) : NonConvexGeometry(1)
{
	gType = SEGMENTED_LINE ;
	this->center = Point() ;
	boundingPoints->resize(points->size()) ;
	std::copy(points->begin(), points->end(), &(*boundingPoints)[0]) ;
}

void SegmentedLine::computeCenter()
{
}

double SegmentedLine::getRadius() const
{
	//! \todo make it do something
	return 0 ;
}

Point * SegmentedLine::getHead() const 
{
	return (*this->begin()) ;
}

Point * SegmentedLine::getTail() const 
{
	return (*this->end()-1) ;
}

void SegmentedLine::sampleBoundingSurface(size_t num_points) 
{ 
}

void SegmentedLine::project(Point *p) const
{
	std::map<double, Segment> m ;
	for(size_t i = 0 ; i < boundingPoints->size()-1 ; i++)
	{
		Segment s(*getBoundingPoint(i), *getBoundingPoint(i+1)) ;
		m[squareDist(p,s.midPoint())] = s ;
		
	}
	Line l(*m.rbegin()->second.first(), *m.rbegin()->second.vector()) ;
	
	Point proj = l.projection(p) ;
	p->x = proj.x ;
	p->y = proj.y ;
}

void SegmentedLine::sampleSurface(size_t num_points) 
{ 
	
} 

bool SegmentedLine::in(const Point v) 
{
	for (size_t i = 0 ; i < boundingPoints->size() ;i++)
		if((*(*boundingPoints)[i]) == v)
			return true ;
	
	return false ;
}

bool SegmentedLine::is1D() const
{ 
	return true ; 
} 
