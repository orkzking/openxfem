#include "geometry_base.h"

using namespace Mu ;

Point::Point() 
{
	this->x = 0 ; 
	this->y = 0 ;
	this->z = 0 ;
	this->id = -1 ;
}

double Point::angle() const
{
	assert(z==0) ;
	return atan2(y,x) ;
}

Point::Point(double x, double y)
{
	//(*(dynamic_cast< valarray<double> *>(this))) ;
	this->x = x ; 
	this->y = y ;
	this->z = 0 ;
	this->id = -1 ;
}

Point::Point(double x, double y, double z)
{
	//(*(dynamic_cast< valarray<double> *>(this))) ;
	this->x = x ; 
	this->y = y ;
	this->z = z ;
	this->id = -1 ;
}

void Point::print() const
{
	std::cout << " (" << x << ", " << y << ") " ;
}

double Point::norm() const
{
	return sqrt(x*x + y*y + z*z) ;
}

double Point::sqNorm() const
{
	return x*x + y*y + z*z;
}
void Point::setX(double v) 
{ 
	x = v ;
}

void Point::setY(double v) 
{ 
	y = v ;
}

void Point::setZ(double v) 
{ 
	z = v ;
}

void Point::set(const Point p)
{
	x = p.x ; 
	y = p.y;
	z = p.z;
}

void Point::set(const Point * p)
{
	x = p->x ; 
	y = p->y;
	z = p->z;
}

void Point::set(double v, double vv)
{
	x = v ; 
	y = vv ;
}

void Point::set(double v, double vv, double vvv)
{
	x = v ; 
	y = vv ;
	z = vvv ;
}

bool Point::operator==(Point p) const
{ 
	return  std::abs(x-p.x) < 2*std::numeric_limits<double>::epsilon() && 
		std::abs(y-p.y) < 2*std::numeric_limits<double>::epsilon() &&  std::abs(z-p.z) < 2*std::numeric_limits<double>::epsilon(); 
}

bool Point::operator!=(Point p) const
{
	return   std::abs(x-p.x) >  2*std::numeric_limits<double>::epsilon() ||
		std::abs(y-p.y) >  2*std::numeric_limits<double>::epsilon() ||  std::fabs(z-p.z) >  2*std::numeric_limits<double>::epsilon(); 
}

Point Point::operator-(Point p) const
{
	Point ret((*this)) ;
	ret.x -= p.x; 
	ret.y -= p.y ; 
	ret.z -= p.z ; 
	return ret ;
}

Point Point::operator-(Vector p) const
{
	Point ret((*this)) ;
	ret.x -= p[0] ; 
	ret.y -= p[1] ; 
	ret.z -= p[2] ; 
	return ret ; 
}

Point Point::operator+(Point p) const
{
	Point ret((*this)) ;
	ret.x += p.x; 
	ret.y += p.y ; 
	ret.z += p.z ; 
	return ret ;
}

void Point::operator+=(Point p)
{
	x += p.x; 
	y += p.y ; 
	z += p.z ; 
}

Point Point::operator+(Vector p) const
{
	Point ret((*this)) ;
	ret.x += p[0] ; 
	ret.y += p[1] ; 
	ret.z += p[2] ; 
	return ret ; 
}

Point Point::operator/(double p) const 
{
	Point ret((*this)) ;
	ret.x /= p ; 
	ret.y /= p ; 
	ret.z /= p ; 
	return ret ; 
}


bool Point::operator <(Point p) const 
{
	return (y < p.y) || ( (y == p.y) && (x < p.x) ||  (y == p.y) && (x == p.x) && (z < p.z));
}

bool Point::operator >(Point p) const 
{
	return (y > p.y) || ( (y == p.y) && (x > p.x) ||  (y == p.y) && (x == p.x) && (z > p.z));
}

Point Point::operator*(double p)  const 
{
	Point ret((*this)) ;
	ret.x *= p ; 
	ret.y *= p ; 
	ret.z *= p ; 
	return ret ; 
}


double Point::operator*(Point p) const
{
	return x*p.x + y*p.y + z*p.z;
}

double Point::operator*(Vector p) const
{
	return x*p[0] + y*p[1] + z*p[2]; 
}

Point Point::operator^(Point p) const
{
	Point ret ;
	ret.x = y*p.z - z*p.y ;
	ret.y = z*p.x - x*p.z ;
	ret.z = x*p.y - y*p.x ;
	
	return ret ;
}

Point Point::operator^(Vector p) const
{
	Point ret ;
	ret.x = y*p[2] - z*p[1] ;
	ret.y = z*p[0] - x*p[2] ;
	ret.z = x*p[1] - y*p[0] ;
	
	return ret ;
}


PointSet::PointSet()
{
	this->boundingPoints = new std::valarray<Point *>(0); 
	this->chullEndPos = 0;
}

std::valarray<Point *> * ConvexGeometry::getBoundingPoints() const
{ 
	return this->boundingPoints ; 
}

Point* ConvexGeometry::getBoundingPoint(size_t i) const
{ 
	return (*boundingPoints)[i] ; 
}

Point * ConvexGeometry::getPoint(size_t i) const
{
	if (i < inPoints->size())
	{
		return (*inPoints)[i] ;
	}
	
	return (*boundingPoints)[i-inPoints->size()] ;
}

std::valarray<Point *> * NonConvexGeometry::getBoundingPoints() const
{ 
	return this->boundingPoints ; 
}

Point* NonConvexGeometry::getBoundingPoint(size_t i) const
{ 
	return (*this->boundingPoints)[i] ; 
}

Point* Geometry::getInPoint(size_t i) const
{ 
	return (*this->inPoints)[i] ; 
}

double & Point::operator[](size_t i)
{
	return (*(&x+i)) ;
}
double Point::operator[](size_t i) const
{
	return (*(&x+i)) ;
}

std::valarray<Point* > * Geometry::getInPoints() const
{
	return this->inPoints ;
}


ConvexGeometry::ConvexGeometry(): ConvexPolygon((size_t)0)
{
}

ConvexGeometry::ConvexGeometry(size_t s) : ConvexPolygon(s)
{
}

const Point* Geometry::getCenter() const
{
	return &center;
}

Point * Geometry::getCenter() 
{
	return &center;
}

PointSet::PointSet(size_t npoints) 
{
	this->boundingPoints = new std::valarray<Point *>(npoints) ; 
	for(size_t i = 0 ; i < npoints ;i++)
		(*this->boundingPoints)[i] = NULL ;
	this->chullEndPos = 0;
} ;

size_t ConvexGeometry::size() const
{
	return 	boundingPoints->size() + inPoints->size() ;
}

size_t NonConvexGeometry::size() const
{
	return 	boundingPoints->size() + inPoints->size() ;
}

double PointSet::x(size_t i) 
{
	return (*boundingPoints)[i]->x ; 
}

double PointSet::y(size_t i) 
{ 
	return (*boundingPoints)[i]->y ; 
}

double PointSet::z(size_t i) 
{ 
	return (*boundingPoints)[i]->z ; 
}

void PointSet::setX(size_t i, double vv) 
{ 
	(*boundingPoints)[i]->Point::setX(vv) ; 
}

void PointSet::setY(size_t i, double vv) 
{ 
	(*boundingPoints)[i]->Point::setY(vv) ; 
}

void PointSet::setZ(size_t i, double vv) 
{ 
	(*boundingPoints)[i]->Point::setZ(vv) ; 
}

void PointSet::set(size_t i, Point * p) 
{ 
	delete (*boundingPoints)[i] ;
	(*boundingPoints)[i] = p; 
}

void PointSet::set(size_t i, double x, double y) 
{ 
	int id = -1 ;
	if((*boundingPoints)[i] != NULL)
		id = (*boundingPoints)[i]->id ;
	
	delete (*boundingPoints)[i] ;
	
	(*boundingPoints)[i] = new Point(x,y) ; 
	(*boundingPoints)[i]->id = id ;
}

void PointSet::set(size_t i, double x, double y, double z) 
{ 
	int id = -1 ;
	if((*boundingPoints)[i] != NULL)
		id = (*boundingPoints)[i]->id ;
	
	delete (*boundingPoints)[i] ;
	
	(*boundingPoints)[i] = new Point(x,y,z) ; 
	(*boundingPoints)[i]->id = id ;
}

Point * PointSet::operator[](size_t i) 
{ 
	return (*boundingPoints)[i] ; 
}

Point * PointSet::operator[](size_t i) const
{ 
	return (*boundingPoints)[i] ; 
}

Point * PointSet::getPoint(size_t i) const
{ 
	return (*boundingPoints)[i] ; 
}

Point * PointSet::getPoint(size_t i)
{ 
	return (*boundingPoints)[i] ; 
}

PointSet::iterator PointSet::begin() const
{
	return &(*boundingPoints)[0] ; 
}

PointSet::iterator PointSet::end() const
{
	return  &(*boundingPoints)[boundingPoints->size()] ;
}

ConvexPolygon * PointSet::convexHull()
{
		//!have we allready done that ?
	
	if(this->chullEndPos != 0 )
	{
		ConvexPolygon *hull = new ConvexPolygon(chullEndPos) ;
		std::copy(begin(), end(), hull->begin()) ;
		return hull ;
	}
	
	ConvexPolygon *ret = new ConvexPolygon(this) ;
	chullEndPos = ret->size() ;
	
	return ret ;
	
}


bool PointSet::in(const Point p)  
{
	ConvexPolygon * hull = convexHull() ;
	bool ret = hull->in(p) ;
	delete hull ;
	return ret ;
}

size_t PointSet::size() const
{ 
	return  boundingPoints->size() ;
}

Point * PointSet::computeCenter()
{
	Point * ret = new Point();
	for (size_t i = 0 ; i < boundingPoints->size() ; i++)
	{
		ret->x += (*boundingPoints)[i]->x/boundingPoints->size() ;
		ret->y += (*boundingPoints)[i]->y/boundingPoints->size() ;
	}
	return ret ;
}

void PointSet::removePoint(size_t index)
{
	std::valarray<Point *> n(size()-1) ;
	//std::copy((*this)[0], (*this)[index-1], (*n)[0]) ;
	for (size_t i = 0 ; i < index ; i++)
	{
		n[i] = (*boundingPoints)[i] ;
	}
	for (size_t i = index+1 ; i < size() ; i++)
	{
		n[i] = (*boundingPoints)[i] ;
	}
	
	boundingPoints->resize(n.size()) ;
	(*boundingPoints) = n ;
}


Geometry::Geometry(): gType(NULL_GEOMETRY)
{
	sampled = false ;
	inPoints  = new std::valarray<Point * >(0) ;
}

Geometry::Geometry(size_t numPoints): gType(NULL_GEOMETRY)
{
	inPoints  = new std::valarray<Point * >(numPoints) ;
}

GeometryType Geometry::getGeometryType() const
{
	return gType ;
}

// ConvexPolygon::ConvexPolygon()
// {
// 	std::cout << "calling default ConvexPolygon constructor" << std::endl ;
// }

ConvexPolygon::ConvexPolygon(size_t npoints) : PointSet(npoints)
{
}

ConvexPolygon::ConvexPolygon(PointSet * po) : PointSet(po->size())
{
	std::vector<Point*> *points = new std::vector<Point*>(po->size()) ;
	std::copy(po->begin(), po->end(), points->begin()) ;
	
	Point *pivot( (*points)[0] ) ;
	for(size_t i = 1 ; i < points->size() ; i++)
	{
		if (pivot->y < (*points)[i]->y)
			pivot = (*points)[i] ;
	}
	
	//! we build a map ordered by the angle to the pivot.
	
	std::map< double, Point *>  pointSet ;
	
	for(size_t i = 0 ; i < points->size() ; i++)
	{
		if ((*(*points)[i]) != (*pivot))
		{
			double angle = atan2(pivot->y-(*points)[i]->y, pivot->x-(*points)[i]->x) ;
			
			if(pointSet.find(angle) != pointSet.end())
			{
				if(squareDist((*pointSet[angle]), (*pivot)) < squareDist((*(*points)[i]), (*pivot)))
				{
					pointSet[angle] = (*points)[i] ;
				}
			}
			else
				pointSet[angle] = (*points)[i] ;
		}
	}
	
	pointSet[-10] = pivot ;
	
	//!we build the convex hull : we are on the convex hull if and only if we only turned left.
	
	std::vector<Point *> temphull ;
	std::vector<Point *> orderedset ;
	
	for(std::map<double, Point*>::const_iterator i = pointSet.begin() ; i!= pointSet.end() ; ++i)
	{
		orderedset.push_back(i->second) ;
	}
	
	temphull.push_back(orderedset[0] ) ;
	temphull.push_back(orderedset[1] ) ;
	
	for(std::vector< Point *>::iterator i = orderedset.begin()+2 ; i != orderedset.end() ; ++i)
	{
		for(size_t j = 0 ; j < temphull.size() ; j++)
			
		//! this is a usual cross product of two vectors...
			if(  ((*i)->y - (*temphull.rbegin())->y)*((*temphull.rbegin())->x - temphull[temphull.size()-2]->x ) - 
			     ((*i)->x - (*temphull.rbegin())->x)*((*temphull.rbegin())->y - temphull[temphull.size()-2]->y ) > 
			     std::numeric_limits<double>::epsilon() )
			{
				temphull.push_back(*i) ;
			}	
		else
		{
			while( !(((*i)->y - (*temphull.rbegin())->y)*((*temphull.rbegin())->x - temphull[temphull.size()-2]->x ) - 
			         ((*i)->x - (*temphull.rbegin())->x)*((*temphull.rbegin())->y -temphull[temphull.size()-2]->y ) > 
			         std::numeric_limits<double>::epsilon()))
			{
				temphull.pop_back();
			}
			temphull.push_back(*i) ;
		}
		
		
	}
	
	delete points ;
	std::copy(temphull.begin(), temphull.end(),this->begin()) ;
}


bool ConvexPolygon::in(const Point p)  
{
	
	bool in = false ;
	
	for (size_t i = 0, j  =  boundingPoints->size()-1; i <  boundingPoints->size(); j = i++)
	{
		if (((((*boundingPoints)[i]->y <= p.y ) && (p.y<(*boundingPoints)[j]->y)) || (((*boundingPoints)[j]->y <= p.y) && (p.y<(*boundingPoints)[i]->y))) &&
		    (p.x < ((*boundingPoints)[j]->x - (*boundingPoints)[i]->x) * (p.y - (*boundingPoints)[i]->y) / ((*boundingPoints)[j]->y - (*boundingPoints)[i]->y) + (*boundingPoints)[i]->x))
			in = !in;
	}
	
	return in ;
	
}

bool ConvexPolygon::in(const Point *p)  const
{
	
	bool in = false ;
	
	for (size_t i = 0, j  =  boundingPoints->size()-1; i <  boundingPoints->size(); j = i++)
	{
		if (((((*boundingPoints)[i]->y <= p->y ) && (p->y < (*boundingPoints)[j]->y)) || 
		     (((*boundingPoints)[j]->y <= p->y) && (p->y<(*boundingPoints)[i]->y))) &&
		    (p->x < ((*boundingPoints)[j]->x - (*boundingPoints)[i]->x) * (p->y - (*boundingPoints)[i]->y) / ((*boundingPoints)[j]->y - (*boundingPoints)[i]->y) + (*boundingPoints)[i]->x))
			in = !in;
	}
	
	return in ;
}

bool ConvexPolygon::in(const Point p)  const
{
	bool in = false ;
	
	for (size_t i = 0, j  =  boundingPoints->size()-1; i <  boundingPoints->size(); j = i++)
	{
		if (((((*boundingPoints)[i]->y <= p.y ) && (p.y<(*boundingPoints)[j]->y)) || (((*boundingPoints)[j]->y <= p.y) && (p.y<(*boundingPoints)[i]->y))) &&
		    (p.x < ((*boundingPoints)[j]->x - (*boundingPoints)[i]->x) * (p.y - (*boundingPoints)[i]->y) / ((*boundingPoints)[j]->y - (*boundingPoints)[i]->y) + (*boundingPoints)[i]->x))
			in = !in;
	}
	
	return in ;
}

bool ConvexPolygon::isTrigoOriented()  const
{
	
	for (size_t i = 0 ;  i <  boundingPoints->size()-1; i++)
	{
		Point v_0 = *(*boundingPoints)[(i+1)%boundingPoints->size()] - *(*boundingPoints)[(i)%boundingPoints->size()] ;
		Point v_1 = *(*boundingPoints)[(i+2)%boundingPoints->size()] - *(*boundingPoints)[(i+1)%boundingPoints->size()] ;
		if((v_0^v_1).z < 0 )
			return false ;
	}
	
	return true ;
	
}

NonConvexGeometry::NonConvexGeometry() : PointSet(1)
{
	orderedSet.resize(1) ;
	orderedSet[0] = (*this->inPoints)[0] ;
}

NonConvexGeometry::NonConvexGeometry(size_t numPoints) : PointSet(numPoints)
{
	orderedSet.resize(numPoints) ;
	for(size_t i = 0 ; i < numPoints ; i++)
		orderedSet[i] = (*this->inPoints)[i] ;
}

NonConvexGeometry::NonConvexGeometry(PointSet * p)
{
	this->boundingPoints = new std::valarray<Point *>(p->size()) ;
	std::copy(p->begin(), p->end(), begin()) ;//&boundingPoints[0]) ;
}

Point * NonConvexGeometry::getPoint(size_t i) const
{
	if (i < inPoints->size())
		return (*inPoints)[i] ;
	return (*boundingPoints)[i-inPoints->size()] ;
}

Line::Line(Point origin, Point vector)
{
	p = origin ;
	v = vector ;
}

bool Line::intersects(const Line *l) const
{
	return v.x * l->vector()->y - v.y * l->vector()->x != 0 ;
}

bool Line::intersects(const Segment *s) const
{
	
	if (v.x * s->vector()->y - v.y * s->vector()->x == 0)
		return false ;
	
	Matrix m(2,2) ;
	Vector vv(2) ;
	
	m[0][0] = v.x ; m[0][1] = -s->vector()->x ;
	m[1][0] = v.y ; m[1][1] = -s->vector()->y ;
	
	vv[0] = s->first()->x - p.x ; vv[1] = s->first()->y - p.y ; 
	
	Matrix m_inv = inverse2x2Matrix(m) ;
	
	Vector fac = m_inv * vv ;
	
	return fac[0] < 1 && fac[0] > 0 ;
}

bool Line::on(const Point *m) const
{
	return isAligned((*m), p, (p+v)) ;
}

const Point * Line::vector() const
{
	return dynamic_cast<const Point *>(&v) ;
}

const Point * Line::origin() const
{
	return &p ;
}

Point Line::intersection(const Line *l) const
{
	double t = 0;
	if(v.x != 0 && v.y != 0)
		t = ((p.y - l->origin()->y) + v.y/v.x * (p.x - l->origin()->x)) / (v.x/v.y - l->vector()->y) ;
	else if (v.x == 0)
		t = ( p.x - l->origin()->x ) / l->vector()->x ;
	else if (v.y == 0)
		t = ( p.y - l->origin()->y ) / l->vector()->y ;
	
	return (*l->origin()) + (*l->vector())*t ;
}

Point Line::intersection(const Segment *s) const
{
	double t = 0;
	if(v.x != 0 && v.y != 0)
		t = ((p.y - s->first()->y) + v.y/v.x * (p.x - s->first()->x)) / (v.x/v.y - s->vector()->y) ;
	else if (v.x == 0)
		t = ( p.x - s->first()->x ) / s->vector()->x ;
	else if (v.y == 0)
		t = ( p.y - s->first()->y ) / s->vector()->y ;
	
	return (*s->first()) + (*s->vector())*t ;
}


bool Line::intersects(const Geometry *g) const
{
	switch(g->getGeometryType())
	{
	case CIRCLE:
		{
			return squareDist(projection(g->getCenter()), (*g->getCenter())) < g->getRadius()*g->getRadius() ;
		}
	default:
		return false ;
	}
}


std::vector<Point> Line::intersection(const Geometry * g) const
{
	switch(g->getGeometryType())
	{
	case CIRCLE:
		{
			double a = v.sqNorm() ;
			double b = p.x*v.x + p.y*v.y ;
			double c = p.sqNorm()-g->getRadius()*g->getRadius() ;
			double delta = b*b - 4*a*c ;
			
			if(delta == 0)
			{
				std::vector<Point> ret ;
				ret.push_back(p+v*(-b/(2.*a))) ;
				return ret ;
			}
			else if (delta > 0)
			{
				std::vector<Point> ret ;
				ret.push_back(p+v*((-b + sqrt(delta))/(2.*a))) ;
				ret.push_back(p+v*((-b - sqrt(delta))/(2.*a))) ;
				return ret ;
			}

			return std::vector<Point>(0) ;
		}
	default:
		return std::vector<Point>() ;
	}
}

Point Line::projection(const Point *m ) const
{
	Line l((*m), Point(-v.y, v.x)) ;
	return l.intersection(this) ;
}


Segment::Segment(const Point p0, const Point p1)
{
	f = p0 ;
	s = p1 ;
	mid = p0*0.5 + p1*0.5;
	vec = p1-p0 ;
}


Segment::Segment()
{
	f = Point(0,0) ;
	s = Point(0,0) ;
	mid = f*0.5 + s*0.5;
	vec = f-s ;
}
Segment::~Segment()
{
}


void Segment::print() const 
{
	std::cout << "[ (" << f.x << ", " << f.y << ") ; (" << s.x << ", " << s.y << ") ]" << std::endl ;
}

bool Segment::intersects(const Line *l) const
{
	double t = 0;
	double t0 = 0;
	
	if(vec.x != 0 && vec.y != 0)
	{
		if(l->vector()->y != 0 && l->vector()->x != 0)
		{
			t = ( l->origin()->y + (f.x  - l->origin()->x) / l->vector()->x * l->vector()->y) /
				(vec.y + vec.x/ l->vector()->x * l->vector()->y) ;
			//t0 = (f.x + t*vec.x - l->origin()->x) / l->vector()->x ;
		}
		else if (l->vector()->y == 0)
		{
			t = (l->origin()->y - f.y) / vec.y ;
			//t0 = (f.x + t*vec.x - l->origin()->x) / l->vector()->x ;
		}
		else if (l->vector()->x == 0)
		{
			t = (l->origin()->x - f.x) / vec.x ;
			//t0 = (f.y + t*vec.y - l->origin()->y) / l->vector()->y ;
		}
	}
	else if (vec.x == 0)
	{
		t0 = ( f.x - l->origin()->x ) / l->vector()->x ;
		t = (l->origin()->y + t0* l->origin()->y - f.y)/vec.y ;
	}
	else if (vec.y == 0)
	{
		t0 =( f.y - l->origin()->y ) / l->vector()->y ;
		t = (l->origin()->x + t0* l->vector()->x - f.x)/vec.x ;
	}
	return t > 0 && t<1 ;
	
	
	
	if (vec.x * l->vector()->y - vec.y * l->vector()->x == 0)
		return false ;
	
	return std::abs((l->vector()->y*vec.x - l->vector()->x*vec.y)/(l->vector()->norm()*vec.norm())) < 1 ;
}

bool Segment::intersects(const Geometry *g) const
{
	switch(g->getGeometryType())
	{
	case CIRCLE: // corrected code by Truong Quang Tri
		{
		   double d1 = squareDist(g->getCenter(),&f);
			double d2 = squareDist(g->getCenter(),&s);
		   double chk= (d1 - g->getRadius()*g->getRadius()) * (d2 - g->getRadius()*g->getRadius());
			return (chk <= 0) ? true : false ;
		}
	case TRIANGLE:
		{
			bool ret = false ;	
			
			for(size_t i = 0 ; i <  g->getBoundingPoints()->size() ;  i++)
			{
				Segment s(*g->getBoundingPoint(i), *g->getBoundingPoint((i+1)%g->getBoundingPoints()->size())) ;
				ret = ret || s.intersects(this) ;
			}
			
			return ret ;
		}
	case RECTANGLE:
		{
			Segment s0(*g->getBoundingPoint(0), *g->getBoundingPoint(1)) ;
			Segment s1(*g->getBoundingPoint(1), *g->getBoundingPoint(2)) ;
			Segment s2(*g->getBoundingPoint(2), *g->getBoundingPoint(3)) ;
			Segment s3(*g->getBoundingPoint(3), *g->getBoundingPoint(0)) ;
			
			return s0.intersects(this) || s1.intersects(this) || s2.intersects(this) || s3.intersects(this);
		}
	default:
		return false ;
	}
}


std::vector<Point> Segment::intersection(const Geometry *g) const
{
	switch(g->getGeometryType())
	{
	case TRIANGLE:
		{
			std::vector<Point> ret ;
			
			for(size_t i = 0 ; i <  g->getBoundingPoints()->size() ;  i++)
			{
				Segment s(*g->getBoundingPoint(i), *g->getBoundingPoint((i+1)%g->getBoundingPoints()->size())) ;
				if(s.intersects(this))
					ret.push_back( s.intersection(this)) ;
			}
			
			return ret ;
		}
	case RECTANGLE:
		{
			std::vector<Point> ret ;
			
			Segment s0(*g->getBoundingPoint(0), *g->getBoundingPoint(1)) ;
			
			if(s0.intersects(this))
				ret.push_back( s0.intersection(this)) ;
			
			Segment s1(*g->getBoundingPoint(1), *g->getBoundingPoint(2)) ;
			
			if(s1.intersects(this))
				ret.push_back( s1.intersection(this)) ;
			
			Segment s2(*g->getBoundingPoint(2), *g->getBoundingPoint(3)) ;
			
			if(s2.intersects(this))
				ret.push_back( s2.intersection(this)) ;
			
			Segment s3(*g->getBoundingPoint(3), *g->getBoundingPoint(0)) ;
			
			if(s3.intersects(this))
				ret.push_back( s3.intersection(this)) ;
			
			return ret ;
		}
	case CIRCLE:
		{
			double a = vec.sqNorm() ;
			double b = f.x*vec.x + f.y*vec.y ;
			double c = f.sqNorm()-g->getRadius()*g->getRadius() ;
			double delta = b*b - 4*a*c ;
			
			if(delta == 0)
			{
				std::vector<Point> ret ;
				ret.push_back(f+vec*(-b/(2.*a))) ;
				return ret ;
			}
			else if (delta > 0)
			{
				std::vector<Point> ret ;
				ret.push_back(f+vec*(-b + sqrt(delta)/(2.*a))) ;
				ret.push_back(f+vec*(-b - sqrt(delta)/(2.*a))) ;
				return ret ;
			}
			else
			{
				return std::vector<Point>(0) ;
			}
		}
	default:
		return std::vector<Point>(0) ;
	}
}


bool Segment::intersects(const Segment *s) const
{
	if (vec.x * s->vector()->y - vec.y * s->vector()->x == 0)
		return false ;
	
	Matrix m(2,2) ;
	Vector v(2) ;
	
	m[0][0] = vec.x ; m[0][1] = -s->vector()->x ;
	m[1][0] = vec.y ; m[1][1] = -s->vector()->y ;
	
	v[0] = s->first()->x - f.x ; v[1] = s->first()->y - f.y ; 
	
	Matrix m_inv = inverse2x2Matrix(m) ;
	
	Vector fac = m_inv * v ;
	
	return fac[0] < 1 && fac[0] > 0 && fac[1] < 1 && fac[1] > 0;
	
}

bool Segment::on(const Point *p) const
{
	//if(vec.x == vec.y == 0) can not use with MSVC
	//	return false ;
   if(vec.x == 0 && vec.y == 0)
		return false ;
	
	if(!isAligned(p, &this->f, &this->s))
		return false ;
	else if(std::abs(vec.x) > std::abs(vec.y))
	{
		double t = (p->x - f.x ) / vec.x ;
		return t > -2.*std::numeric_limits<double>::epsilon() && t<1 + 2.*std::numeric_limits<double>::epsilon() ;
	}
	else
	{
		double t = (p->y - f.y ) / vec.y ;
		return t > -2.*std::numeric_limits<double>::epsilon() && t<1 + 2.*std::numeric_limits<double>::epsilon() ;
	}
}

void Segment::setFirst(Point p) 
{
	f = p ;
	mid = f*0.5 + s*0.5;
	vec = s-f ;
}

void Segment::setFirst(double x, double y)
{
	f.set(x, y) ;
	mid = f*0.5 + s*0.5;
	vec = s-f ;
}

void Segment::setSecond(Point p)
{
	s = p ;
	mid = f*0.5 + s*0.5;
	vec = s-f ;
}

void Segment::setSecond(double x, double y)
{
	s.set(x, y) ;
	mid = f*0.5 + s*0.5;
	vec = s-f ;
}

void Segment::set(Point p0, Point p1)
{
	f.set(p0) ;
	s.set(p1) ;
	mid = f*0.5 + s*0.5;
	vec = s-f ;
}

void Segment::set(double x0, double y0, double x1, double y1)
{
	f.set(x0, y0) ;
	s.set(x1, y1) ;
	mid = f*0.5 + s*0.5;
	vec = s-f ;
}

const Point * Segment::first() const
{
	return &f ;
}

const Point * Segment::second() const
{
	return &s ;
}

Point * Segment::first() 
{
	return &f ;
}

Point * Segment::second() 
{
	return &s ;
}


const Point * Segment::midPoint() const
{
	return &mid ;
}

const Point * Segment::vector() const
{
	return &vec ;
}


Point Segment::intersection(const Line *l) const
{
	Matrix m(2,2) ;
	Vector v(2) ;
	
	m[0][0] = vec.x ; m[0][1] = -l->vector()->x ;
	m[1][0] = vec.y ; m[1][1] = -l->vector()->y ;
	
	v[0] = l->origin()->x - f.x ; v[1] = l->origin()->y - f.y ; 
	
	Matrix m_inv = inverse2x2Matrix(m) ;
	
	Vector fac = m_inv * v ;
	
	return f + vec*fac[0];
}

Point Segment::intersection(const Segment *s) const
{
	
	Matrix m(2,2) ;
	Vector v(2) ;
	
	m[0][0] = vec.x ; m[0][1] = -s->vector()->x ;
	m[1][0] = vec.y ; m[1][1] = -s->vector()->y ;
	
	v[0] = s->first()->x - f.x ; v[1] = s->first()->y - f.y ; 
	
	Matrix m_inv = inverse2x2Matrix(m) ;
	
	Vector fac = m_inv * v ;
	
	return f + vec*fac[0];
}

bool isInTriangle(const Point test, const Point p0, const Point p1, const Point p2) 
{
	return isOnTheSameSide( test, p0, p1, p2) && isOnTheSameSide(test, p1, p0, p2) && isOnTheSameSide(test, p2, p1, p2) ;
}

bool isOnTheSameSide(const Point test, const Point witness, const Point f0, const Point f1) 
{
	return (((f1.x-f0.x)*(test.y-f0.y) - (f1.y-f0.y)*(test.x-f0.x)) *
	        ((f1.x-f0.x)*(witness.y-f0.y) - (f1.y-f0.y)*(witness.x-f0.x)) > -2*std::numeric_limits<double>::epsilon()) ;
}

double dist(const Point v1, const Point v2)
{
	return sqrt ((v2.x-v1.x)*(v2.x-v1.x)+(v2.y-v1.y)*(v2.y-v1.y)) ;
}

double squareDist(const Point v1, const Point v2)
{
	return (v2.x-v1.x)*(v2.x-v1.x)+(v2.y-v1.y)*(v2.y-v1.y) ;
}

double squareDist(const Point *v1, const Point *v2)
{
	return (v2->x-v1->x)*(v2->x-v1->x)+(v2->y-v1->y)*(v2->y-v1->y) ;
}

ConvexPolygon* convexHull(const std::vector<Point *> * points)
{
	Point *pivot = (*points)[0]  ;
	for(size_t i = 1 ; i < points->size() ; i++)
	{
		if (pivot->y < (*points)[i]->y)
			pivot = (*points)[i] ;
	}
	std::cout << "pivot = " << pivot->x << ", " << pivot->y << std::endl ;
	
	//!then we build a map ordered by the angle to the pivot.
	
	std::map< double, Point* >  pointSet ;
	
	for(size_t i = 0 ; i < points->size() ; i++)
	{
		if ((*(*points)[i]) != (*pivot))
		{
			double angle = atan2(pivot->y-(*points)[i]->y, pivot->x-(*points)[i]->x) ;
			
			if(pointSet.find(angle) != pointSet.end())
			{
				std::cout << "already there..." << std::endl ;
				if(squareDist((*pointSet[angle]), (*pivot)) < squareDist((*(*points)[i]), (*pivot)))
				{
					pointSet[angle] = (*points)[i] ;
				}
			}
			else
				pointSet[angle] = (*points)[i] ;
		}
	}
	
	pointSet[-10] = pivot ;
	
	//!we build the convex hull : we are on the convex hull if and only if we only turned left.
	
	std::vector<Point*> temphull ;
	std::vector<Point*> orderedset ;
	
	for(std::map<double, Point*>::const_iterator i = pointSet.begin() ; i!= pointSet.end() ; ++i)
	{
		std::cout << "point " <<  i->second->x << ", " <<  i->second->y << std::endl ;
		orderedset.push_back(i->second) ;
	}
	
	temphull.push_back(orderedset[0] ) ;
	temphull.push_back(orderedset[1] ) ;
	
	for(std::vector< Point *>::iterator i = orderedset.begin()+2 ; i != orderedset.end() ; ++i)
	{
		for(size_t j = 0 ; j < temphull.size() ; j++)
			std::cout << "(" << temphull[j]->x << ", " << temphull[j]->y << ")" << std::endl ;
		
		//! this is a usual cross product of two vectors...
		if(  ((*i)->y - (*temphull.rbegin())->y)*((*temphull.rbegin())->x - temphull[temphull.size()-2]->x ) - 
		     ((*i)->x - (*temphull.rbegin())->x)*((*temphull.rbegin())->y - temphull[temphull.size()-2]->y ) > 
		     std::numeric_limits<double>::epsilon() )
		{
			temphull.push_back(*i) ;
			std::cout << "new point in hull = " <<  (*i)->x << ", " <<  (*i)->y << std::endl ;
		}	
		else
		{
			while( !(((*i)->y - (*temphull.rbegin())->y)*((*temphull.rbegin())->x - temphull[temphull.size()-2]->x ) - 
			         ((*i)->x - (*temphull.rbegin())->x)*((*temphull.rbegin())->y -temphull[temphull.size()-2]->y ) > 
			         std::numeric_limits<double>::epsilon()))
			{
				std::cout << "out of the hull = " <<  (*temphull.rbegin())->x << ", " <<  (*temphull.rbegin())->y << std::endl ;
				temphull.pop_back();
			}
			temphull.push_back(*i) ;
			std::cout << "new point in hull = " <<  (*i)->x << ", " <<  (*i)->y << std::endl ;
		}
		
		
	}
	
	ConvexPolygon *hull = new ConvexPolygon(temphull.size()) ;
	
	std::copy(temphull.begin(), temphull.end(),hull->begin()) ;
	
	for(size_t i = 0 ; i < hull->size() ; i++)
		std::cout << "(" << (*hull)[i]->x << ", " << (*hull)[i]->y << ")" << std::endl ;
	
	return hull ;
} 

