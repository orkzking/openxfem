//   file PIECEWIS.CPP
 
#include "piecewis.h"
#include <math.h>


double  PiecewiseLinFunction :: at (double time)
   // Returns the value of the receiver at time 'time'.
{
	int    index;
	double answer;

	if (! numberOfPoints)
		this -> getPoints() ;

	if (time > dates[numberOfPoints-1]) return values[numberOfPoints-1];
	if (time < dates[0]) return values[0];

	for (int i=0 ; i<numberOfPoints ; i++) {
		if (time - dates[i] >= 0.) {
			index = i;
		}
	}

	answer = values[index] + ((values[index+1] - values[index]) * ((time - dates[index]) / (dates[index+1] - dates[index])));
	return answer;
}


void  PiecewiseLinFunction :: getPoints ()
   // Reads in the data file the date and the value of every point of the
   // receiver.
{
   int i ;

   numberOfPoints = this->readInteger("nPoints") ;

   dates = new double [numberOfPoints] ;
   for (i=0 ; i<numberOfPoints ; i++)
      dates[i] = this -> read("t",i+1) ;

   values = new double [numberOfPoints] ;
   for (i=0 ; i<numberOfPoints ; i++)
      values[i] = this -> read("f(t)",i+1) ;
}


