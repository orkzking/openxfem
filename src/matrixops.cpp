#include "matrixops.h"

using namespace Mu ;

Matrix::Matrix(size_t x, size_t y)
{
	r=x ;
	c=y ;
	v = new std::valarray<double>(0.,x*y) ;
}

double& Matrix::operator()(size_t x, size_t y)
{
	return row(x)[y] ;
}

double Matrix::operator()(size_t x, size_t y) const
{
	return row(x)[y] ;
}

Matrix Matrix::transpose() const
{
	Matrix ret(c,r) ;
	for(size_t i = 0 ; i < r ; i++)
	{
		for(size_t j = 0 ; j < c ; j++)
		{
			ret[j][i] = (*this)[i][j] ;
		}
	}
	
	return ret ;
}

Matrix &Matrix::operator*=(double d)
{
	(*v) *=d ;
	return *this ;
}

Matrix Matrix::operator*(double d) const
{
	Matrix ret(*this) ;
	ret*=d ;
		
	return ret ;
}

Matrix Matrix::operator/(double d) const
{
	Matrix ret(*this) ;
	ret /=d ;
	return ret ;
}

Matrix &Matrix::operator/=(double d)
{
	(*v) /=d ;
	return *this ;
}

Matrix & Matrix::operator =(const MtM& m)
{
	(*this) = (Matrix) m ;
	
	return *this ;
}

Matrix &Matrix::operator*=(const Matrix &m)
{
	assert(m.numRows() == this->numCols()) ;
	
	Matrix ret =matrix_multiply((*this), m) ;
	
	(*v) = ret.array() ;
	r = ret.numRows() ;
	c = ret.numCols() ;
	
	return *this ;
}


Matrix &Matrix::operator +=(const Matrix &m)
{
	(*v) += (m.array()) ;
	return *this ;
}

Matrix Matrix::operator +(const  Matrix &m) const
{
	Matrix ret(*this) ;
	ret += m ; ;
	return ret ;
}

Matrix &Matrix::operator -=(const Matrix &m)
{
	(*v) -= (m.array()) ;
	return *this ;
}

Matrix Matrix::operator -(const  Matrix &m) const
{
	Matrix ret(*this) ;
	ret -= m ; ;
	return ret ;
}

bool Matrix::operator ==(const Matrix &m)
{
	if(v->size() != m.array().size())
		return false ;
	
	for(size_t i = 0 ; i < v->size() ; i++)
		if((*v)[i] != m.array()[i])
			return false ;
	
	return true ;
}

bool Matrix::operator !=(const Matrix &m)
{
	if(v->size() != m.array().size())
		return true ;
	else
		for(size_t i = 0 ; i < v->size() ; i++)
			if((*v)[i] != m.array()[i])
				return true ;
	
	return false ;
}

Matrix::Matrix(const Matrix& m) : r(m.numRows()), c( m.numCols())
{
	v = new std::valarray<double>(m.size()) ;
	(*v) = m.array() ;
}

Matrix &Matrix::operator =(const Matrix &m)
{
	delete v ;
	v = new std::valarray<double>(m.size()) ;
	(*v) = m.array() ;
	r = m.numRows() ;
	c = m.numCols() ;
	return *this ;
}

MtM::operator const Matrix() const
{
	return matrix_multiply(first, second) ;
}


MtV::operator const Vector()
{
	return matrix_vector_multiply(m,v) ;
}

