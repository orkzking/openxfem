#ifndef __GEOMETRY_BASE_H_
#define __GEOMETRY_BASE_H_



#include <valarray>
#include <vector>
#include <map>
#include <limits>
#include <iostream>
#include "matrixops.h"


namespace Mu
{
typedef enum
{
	NULL_GEOMETRY,
	CIRCLE,
	TRIANGLE,
	RECTANGLE,
	CONVEX_POLYGON,
	SEGMENTED_LINE,
	TETRAHEDRON
} GeometryType ;

class Segment ;
class ConvexPolygon ;

class Point
{
	
public:
	double x ;
	double y ;
	double z ;
	int id ;
	Point();
	Point(double x, double y) ;
	Point(double x, double y, double z) ;
	
	void setX(double v) ;
	void setY(double v) ;
	void setZ(double v) ;
	void set(double v, double vv) ;
	void set(double v, double vv, double vvv) ;
	void set(const Point p) ;
	void set(const Point * p) ;
	
	
	/** Robust point comparator.
	 * 
	 * @param p \c Point to compare to.
	 * @return is within 2*std::numeric_limits\<double\>::epsilon() of p.
	 */
	bool operator==(Point p) const ;
	
	/** Robust point comparator.
	 * 
	 * @param p \c Point to compare to.
	 * @return is outside of 2*std::numeric_limits\<double\>::epsilon() of p.
	 */
	bool operator!=(Point p) const ;
	
	/** Operator defining order in the plane.
	 * 
	 * @param p \c Point to compare to 
	 * @return this is bottommost then leftmost compared to p.
	 */
	bool operator <(Point p) const ;
	
	/** Operator defining order in the plane.
	 * 
	 * @param p \c Point to compare to 
	 * @return this is topmost then rightmost compared to p.
	 */
	bool operator >(Point p) const ;
	
	Point operator-(Point p) const ;
	Point operator-(Vector p) const ;
	Point operator+(Point p) const ;
	Point operator+(Vector p) const ;
	Point operator/(double p) const ;
	Point operator*(double p) const ;
	double operator*(Point p) const ;
	double operator*(Vector p) const ;
	Point operator^(Point p) const ;
	Point operator^(Vector p) const ;
	
	void operator+=(Point p) ;
	
	/** Returns the norm of the Vector. The norm is simply
	 * \f$ \sqrt{x^2 + y^2} \f$.
	 * 
	 * @return the norm of the vector, 
	 */
	double norm() const;
	double sqNorm() const;
	void print() const  ;
	double angle() const ;
	double & operator[](size_t i) ;
	double operator[](size_t i) const ;
	
	
} ;

class Geometry 
{
protected:
	
	std::valarray<Point *> * inPoints ;
	bool sampled ;
	
	Point center ;
	
	virtual void computeCenter() = 0;
	
	GeometryType gType ;
	
public:
	
	Geometry() ;
	Geometry(size_t numPoints) ;
	virtual ~Geometry() { delete inPoints ;}
	
	virtual std::valarray<Point *> * getBoundingPoints() const = 0;
	virtual Point* getBoundingPoint(size_t i) const = 0;
	virtual std::valarray<Point *> * getInPoints() const ;
	virtual Point* getInPoint(size_t i) const ;
	virtual Point * getPoint(size_t i) const = 0;
	virtual const Point * getCenter() const ;
	virtual Point * getCenter() ;
	virtual void project(Point *) const = 0;
	
	GeometryType getGeometryType() const ;
	
	virtual double getRadius() const = 0;
	
	virtual void sampleBoundingSurface(size_t num_points) = 0 ;
	virtual void sampleSurface(size_t num_points) = 0 ;
	virtual bool in(const Point p) = 0;
	
	virtual bool in(const Point *p) const { return false ;}
	
	virtual size_t size() const = 0 ;
	virtual double area() const = 0;
	virtual size_t sides() const { return 3 ; }
	
	virtual bool is1D() const = 0 ;
} ;

class Line
{
protected:
	Point p ;
	Point v ;
public:
	Line(Point origin, Point vector) ;
	
	bool intersects(const Line *l) const;
	bool intersects(const Segment *s) const;
	bool intersects(const Geometry *g) const;
	bool on(const Point *p) const;
	
	std::vector<Point> intersection(const Geometry * g) const ;
	Point intersection(const Line *l) const;
	Point intersection(const Segment *l) const;
	
	const Point * vector() const ;
	const Point * origin() const ;
	
	Point projection(const Point *p ) const ;
	
};

class Segment
{
protected:
	Point f ;
	Point s ;
	Point mid ;
	Point vec ;
	
public:
	Segment(const Point p0, const  Point p1) ;
	Segment() ;
	virtual ~Segment() ;
	
	bool intersects(const Line *l) const;
	bool intersects(const Segment *s) const;
	bool intersects(const Geometry *g) const;
	bool on(const Point *p) const;
	
	void setFirst(Point p) ;
	void setFirst(double x, double y) ;
	void setSecond(Point p) ;
	void setSecond(double x, double y) ;
	void set(Point p0, Point p1) ;
	void set(double x0, double y0, double x1, double y1) ;
	
	void print() const ;
	
	const Point * first() const;
	const Point * second() const;
	
	/** first point of the segment.*/
	Point * first() ;
	/** second point of the segment.*/
	Point * second() ;
	
	/** midpoint of the segment. is recalculated if the endpoints change.*/
	const Point * midPoint() const;
	const Point * vector() const ;
	
	Point intersection(const Line *l) const;
	Point intersection(const Segment *l) const;
	std::vector<Point> intersection(const Geometry * g) const;
} ;


class PointSet
{
protected:
	std::valarray<Point *> * boundingPoints  ;
	size_t chullEndPos ;
	
public:
	PointSet() ;
	PointSet(size_t npoints) ;
	
	virtual ~PointSet() { delete boundingPoints ;}
	double x(size_t i);
	double y(size_t i) ;
	double z(size_t i) ;
	
	void setX(size_t i, double v);
	void setY(size_t i, double v) ;
	void setZ(size_t i, double v) ;
	void set(size_t i, Point *p);
	void set(size_t i, double x, double y) ;
	void set(size_t i, double x, double y, double z) ;
	
	Point * operator[](size_t i) ;
	Point * operator[](size_t i) const;
	Point * getPoint(size_t i) const;
	Point * getPoint(size_t i) ;
	
	virtual bool in(const Point p)  ;
	virtual size_t size() const;
	
	void removePoint(size_t index) ;
	Point * computeCenter() ;
	
	ConvexPolygon * convexHull() ;
	
	typedef Point** iterator;
	
	iterator begin() const ;
	iterator end() const ;
} ;

class NonConvexGeometry : public PointSet, public Geometry
{
protected:
	std::valarray<Point *> orderedSet ;
	std::vector<size_t> stopPos ;
	
public:
	NonConvexGeometry() ;
	NonConvexGeometry(size_t numPoints) ;
	NonConvexGeometry(PointSet * p) ;
	virtual ~NonConvexGeometry() { } ;
	
	virtual std::valarray<Point * > * getBoundingPoints() const ;
	virtual Point* getBoundingPoint(size_t i) const ;
	virtual Point * getPoint(size_t i) const ;
	virtual size_t size() const ;
	virtual double area() const = 0;
	
	virtual void project(Point *) const = 0;
	
} ;

class ConvexPolygon : public PointSet
{
public:
	ConvexPolygon(size_t npoints) ;
	ConvexPolygon(PointSet * po) ;
	virtual ~ConvexPolygon() { } ;
	
	virtual bool in(const Point p) ;
	virtual bool in(const Point p) const;
	virtual bool in(const Point *p)  const ;
	virtual bool isTrigoOriented()  const ;
} ;

class ConvexGeometry :  public ConvexPolygon, public Geometry
{
public:
	ConvexGeometry() ;
	ConvexGeometry(size_t numPoints) ;
	virtual ~ConvexGeometry() { } ;
	
	virtual std::valarray<Point *> * getBoundingPoints() const ;
	virtual Point* getBoundingPoint(size_t i) const ;
	
	virtual void sampleBoundingSurface(size_t num_points) = 0 ;
	virtual void sampleSurface(size_t num_points) = 0 ;
	virtual Point * getPoint(size_t i) const ;
	virtual size_t size() const ;
	virtual double area() const = 0;
	
	virtual void project(Point *) const = 0;
	
} ;


} ;


/** Check the alignment of three points.
 * 
 * @param test first point.
 * @param f0 second point.
 * @param f1 third point.
 * @return true if the points are aligned.
 */
inline bool isAligned(const Mu::Point test, const Mu::Point f0, const Mu::Point f1)  
{
	return ( std::abs((f1.x-test.x)*(f0.y-test.y) - (f0.x-test.x)*(f1.y-test.y)) < 10*std::numeric_limits<double>::epsilon()) ;
} ;

/** Check the alignment of three points.
 * 
 * @param test first point.
 * @param f0 second point.
 * @param f1 third point.
 * @return true if the points are aligned.
 */
inline bool isAligned(const Mu::Point *test, const Mu::Point *f0, const Mu::Point *f1)  
{
	return ( std::abs((f1->x-test->x)*(f0->y-test->y) - (f0->x-test->x)*(f1->y-test->y)) < 10*std::numeric_limits<double>::epsilon()) ;
} ;


/** Test if a point is in a triangle defined by three points.
 * 
 * @param test point to test.
 * @param p0 vertex 0.
 * @param p1 vertex 1.
 * @param p2 vertex 2.
 * @return true if test is in.
 */
bool isInTriangle(const Mu::Point test, const Mu::Point p0, const Mu::Point p1, const Mu::Point p2)  ;

/** Test wether two points lie on the same demi-plane.
 * 
 * @param test first point to test.
 * @param witness second point to test.
 * @param f0 first point defining the plane boundary.
 * @param f1 second point defining the plane boundary.
 * @return true if both points are on the same side of the demi-plane.
 */
bool isOnTheSameSide(const Mu::Point test, const Mu::Point witness, const Mu::Point f0, const Mu::Point f1)  ;

//bool isAligned(const Point test, const Point f0, const Point f1)  ;

/**Return the distance between two points
 * 
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ \sqrt{(x_0-x_1)^2 + (y_0-y_1)^2} \f$
 */
double dist(const Mu::Point v1, const Mu::Point v2) ;

/**Return the sqare distance between two points
 * 
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ (x_0-x_1)^2 + (y_0-y_1)^2 \f$
 */
double squareDist(const Mu::Point v1, const Mu::Point v2) ;

/**Return the sqare distance between two points
 * 
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ (x_0-x_1)^2 + (y_0-y_1)^2 \f$
 */
double squareDist(const Mu::Point *v1, const Mu::Point *v2) ;


/** Return the convex hull of a set of points.
 * 
 * @param points 
 * @return a convex polygon (all boundary points anti-cockwise-ordered).
 */
Mu::ConvexPolygon* convexHull(const std::vector<Mu::Point> * points) ;

#endif
