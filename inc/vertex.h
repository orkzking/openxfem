//   ********************************
//   ***       CLASS VERTEX       ***
//   ********************************

 
#ifndef _VERTEX_H_
#define _VERTEX_H_

#include "geometryentity.h"
#include "element.h"
#include "piecewiselinear.h"
#include <iostream>
//class Element;


//! Vertex class 
/*!
 This class implements a point in space.\n
 */


class Vertex : public GeometryEntity
{

friend class PiecewiseLinear; // so that these classes can access private data member
friend class CrackTip;        // myParent of this class 
public:

	Vertex(int,Domain*); //!< Constructor
	~Vertex();           //!< Destructor
	
   /** @name Get and Give coordinates 
	 *  get the coord. from the input file and return this coord.
	 */
   //@{
	void           getCoordinates () ;    
	FloatArray*    giveCoordinates();
	double         giveCoordinate(int) ;
	//@}

	/*!
	 Update the coordinates of the receiver. Used to update the crack tip.
	 */
    void  setCoordinates(double x, double y) ; 
	/*!
	 Check if a Vertex is in an element or not
	 Override the virtual function of base class GeomeryEntity since Vertex has
	 special treatment.
	 */
	 bool  interactsWith(Element*); 
	/*!
	 Check the intersection of the receiver with a triangle
	 More precisely, check if the segment containing the receiver intersects with tri or not
	 */
	bool   intersects(const Mu::Triangle* tri);
	bool   intersects(const Mu::Segment* s);
	/*!
	 Compute the intersection of the receiver with a triangle
	 @param tri the triangle 
	 @return intersection points stored in std::vector<>
	 */
	std::vector<Mu::Point> intersection(const Mu::Triangle* tri);
	std::vector<Mu::Point> intersection(const Mu::Segment * s);
	PiecewiseLinear*       giveMyParent(){return myParent;} 

	/*!
    Export geometry to Matlab file for plotting.
	 \todo need to implement this ?
    */
	void  exportToMatlab(std::string&){;}
	/*! 
	 Print the coordinate of the receiver on the screen. 
	 Debug only. 2005-09-06
	 */
	void  print() ;

private:
	FloatArray*         coordinates; //!< coord. of the Vertex
	PiecewiseLinear*    myParent;    //!< parent of Vertex is a PiecewiseLinear, need to be accessed by external : PiecewiseLinear
	
} ;


#endif // _VERTEX_H_
