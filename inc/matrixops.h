// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __MATRIXOPS_H__
#define __MATRIXOPS_H__

#include <valarray>
#include "sliceiters.h"
#include <iostream>
#include <numeric>

#include <assert.h>

typedef std::valarray<double> Vector ;


namespace Mu
{

struct MtM ;

class Matrix 
{
	std::valarray< double > *v;
	size_t r, c;

public:
	Matrix(size_t x, size_t y);
	Matrix()
	{
		r = 2 ;
		c = 2 ;
		v = new std::valarray< double >(0.,4) ;
	}
	Matrix(const Matrix&) ;
	virtual ~Matrix() {delete v; }
	
	Matrix &operator=(const Matrix &m);
	
	size_t size() const {return r*c ;}
	size_t numCols() const {return c ;}
	size_t numRows() const {return r ;}
	
	Slice_iter< double > column(size_t i )
	{
		return Slice_iter< double >(v, std::slice(i, r, c)) ;
	}
	Cslice_iter< double > column(size_t i ) const 
	{
		return Cslice_iter< double >(v, std::slice(i, r, c)) ;
	}
	
	Slice_iter< double > row(size_t i )
	{
		return Slice_iter< double >(v, std::slice(i*c, c, 1)) ;
	}
	Cslice_iter< double > row(size_t i ) const
	{
		return Cslice_iter< double >(v, std::slice(i*c, c, 1)) ;
	}
	
	Matrix transpose() const ;

	double& operator()(size_t x, size_t y) ;
	double operator()(size_t x, size_t y) const;
	
	Slice_iter< double > operator()(size_t i) {return row(i) ;}
	Cslice_iter< double > operator()(size_t i) const {return row(i) ;}
	
	Slice_iter< double > operator[](size_t i) {return row(i) ;}
	Cslice_iter< double > operator[](size_t i) const {return row(i) ;}
	
	Matrix& operator *=(double);
	Matrix operator *(double) const;
	Matrix operator /(double) const;
	Matrix& operator /=(double) ;
	Matrix& operator *=(const Matrix &m);
// 	const Matrix& operator *(const Matrix &m) const ;
	Matrix& operator +=(const Matrix &m) ;
	Matrix operator +(const Matrix &m) const;
	Matrix& operator -=(const Matrix &m) ;
	Matrix operator -(const Matrix &m) const;
	Matrix& operator =(const MtM& m) ;
	
	bool operator ==(const Matrix &m) ;
	bool operator !=(const Matrix &m) ;
	
	std::valarray< double > &array() {return *v ;}
	std::valarray< double > array() const {return *v ;}
} ;


struct MtV
{
	const Matrix &m;
	const Vector &v;
	
	MtV(const Matrix &mm, const Vector &vv) : m(mm), v(vv) { }
	
	operator const Vector();
} ;


struct MtM
{
	const Matrix &first;
	const Matrix &second;
	
	MtM(const Matrix &mm, const Matrix &mmm) : first(mm), second(mmm) { }
	
	operator const Matrix() const;
} ;

}

inline Mu::MtV operator*(const Mu::Matrix& mm, const Vector& v)
{
	return Mu::MtV(mm, v) ;
} ;

inline Mu::MtM operator*(const Mu::Matrix& mm, const Mu::Matrix& mmm)
{
	return Mu::MtM(mm, mmm) ;
} ;

inline const Mu::Matrix matrix_multiply(const Mu::Matrix &m0, const Mu::Matrix &m1 )
{
	assert(m0.numCols() == m1.numRows()) ;
	
	Mu::Matrix ret(m0.numRows(), m1.numCols()) ;
	
	for(size_t i = 0 ; i < m0.numRows() ; i++)
	{
		for(size_t j = 0 ; j < m1.numCols() ; j++)
		{
			const Mu::Cslice_iter<double>& ri = m0.row(i) ;
			const Mu::Cslice_iter<double>& cj = m1.column(j) ;
			ret[i][j] = std::inner_product(&ri[0], &ri[m0.numCols()], cj, double(0) ) ;
		}
	}
	return ret ;
}

inline const Vector matrix_vector_multiply(const Mu::Matrix &m, const Vector &v )
{
	assert(m.numRows() == v.size()) ;
	
	Vector ret(0.,v.size()) ;
	
	for(size_t i = 0 ; i < m.numRows() ; i++)
	{
		const Mu::Cslice_iter<double>& ri = m.row(i) ;
      for(size_t j = 0 ; j < v.size() ; j++)
		  //ret[i] = std::inner_product(ri, ri.end(), &v[0], double(0) ) ;
		  ret[i] += ri[j] * v[j] ;
	}
	return ret ;
}

inline const Vector operator*(const Vector &v , const Mu::Matrix &m )
{
	assert(m.numCols() == v.size()) ;
	
	Vector ret(0.,v.size()) ;
	
	for(size_t i = 0 ; i < m.numCols() ; i++)
	{

		const Mu::Cslice_iter<double>& ri = m.column(i) ;
		for(size_t j = 0 ; j < v.size() ; j++)
		  //ret[i] = std::inner_product(ri, ri.end(), &v[0], double(0) ) ;
		  ret[i] += ri[j] * v[j] ;
	}
	return ret ;
}

//clever 2x2 Matrix inversion. Thanks the numerical cookbook :)
inline Mu::Matrix inverse2x2Matrix(const Mu::Matrix s)
{

	if(s[0][0] == 0 || s[1][1] == 0)
	{
		Mu::Matrix swap(2,2) ; swap[0][0] = 0 ; swap[0][1] = 1 ; swap[1][0] = 1 ; swap[1][1] = 0 ;
		Mu::Matrix s_ = s*swap ;
	
		Mu::Matrix ret(2,2) ;
		double r1 = 1./s_[0][0] ;
		double r2 = s_[1][0] * r1 ;
		double r3 = r1* s_[0][1] ;
		double r4 = s_[1][0] * r3 ;
		double r5 = r4 - s_[1][1] ;
		double r6 = 1./r5 ;
		ret[0][1] = r3*r6 ;
		ret[1][0] = r6*r2 ;
		double r7 = r3*ret[1][0] ;
		ret[0][0] = r1-r7 ;
		ret[1][1] = -r6 ;
		return swap*ret ;
	}
	
	Mu::Matrix ret(2,2) ;
	double r1 = 1./s[0][0] ;
	double r2 = s[1][0] * r1 ;
	double r3 = r1* s[0][1] ;
	double r4 = s[1][0] * r3 ;
	double r5 = r4 - s[1][1] ;
	double r6 = 1./r5 ;
	ret[0][1] = r3*r6 ;
	ret[1][0] = r6*r2 ;
	double r7 = r3*ret[1][0] ;
	ret[0][0] = r1-r7 ;
	ret[1][1] = -r6 ;
	
	return ret ;
} ;

inline Mu::Matrix inverse3x3Matrix(const Mu::Matrix m)
{
	Mu::Matrix ret(3,3) ;
	double r11 = m[1][1]*m[2][2]-m[1][2]*m[2][1] ;
	double r21 = m[1][2]*m[2][0]-m[1][0]*m[2][2] ;
	double r31 = m[1][0]*m[2][1]-m[1][1]*m[2][0] ;
	double det = m[0][0]*(r11) + m[0][1]*(r21) + m[0][2]*(r31) ;
	
	ret[0][0] = r11 ; 
	ret[0][1] = m[0][2]*m[2][1] - m[0][1]*m[2][2] ;
	ret[0][2] = m[0][1]*m[1][2] - m[0][2]*m[1][1] ;
	
	ret[1][0] = r21 ; 
	ret[1][1] = m[0][0]*m[2][2] - m[0][2]*m[2][0] ;
	ret[1][2] = m[0][2]*m[1][0] - m[1][1]*m[1][2] ;
	
	ret[2][0] = r31 ; 
	ret[2][1] = m[0][1]*m[2][0] - m[1][1]*m[2][1] ;
	ret[2][2] = m[1][1]*m[2][2] - m[0][1]*m[1][0] ;
	
	return ret/det ;
} ;



# endif  // __MATRIXOPS_H__
