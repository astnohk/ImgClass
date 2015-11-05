#include <cmath>
#include "Vector.h"



template <class T>
VECTOR_2D<T>::VECTOR_2D(void)
{
	x = 0;
	y = 0;
}

template <class T>
VECTOR_2D<T>::VECTOR_2D(const T& init_x, const T& init_y)
{
	x = init_x;
	y = init_y;
}

template <class T>
void
VECTOR_2D<T>::reset(void)
{
	x = 0;
	y = 0;
}

template <class T>
void
VECTOR_2D<T>::reset(const T& init_x, const T& init_y)
{
	x = init_x;
	y = init_y;
}


template <class T>
VECTOR_2D<T> &
VECTOR_2D<T>::operator+=(const VECTOR_2D<T>& vector)
{
	if (this != &vector) {
		this->x += vector.x;
		this->y += vector.y;
	}
	return *this;
}

template <class T>
VECTOR_2D<T> &
VECTOR_2D<T>::operator-=(const VECTOR_2D& vector)
{
	if (this != &vector) {
		this->x -= vector.x;
		this->y -= vector.y;
	}
	return *this;
}

template <class T>
VECTOR_2D<T> &
VECTOR_2D<T>::operator*=(const T& value) // Scalar multiplication
{
	this->x *= value;
	this->y *= value;
	return *this;
}


template <class T>
bool
VECTOR_2D<T>::operator==(const VECTOR_2D<T>& vector)
{
	if (this->x == vector.x
	    && this->y == vector.y) {
		return true;
	} else {
		return false;
	}
}

template <class T>
bool
VECTOR_2D<T>::operator!=(const VECTOR_2D<T>& vector)
{
	if (this->x == vector.x
	    && this->y == vector.y) {
		return false;
	} else {
		return true;
	}
}




//
// ----- Global operators -----
//
template<class Type>
VECTOR_2D<Type> &
operator+(const VECTOR_2D<Type> vector)
{
	return vector;
}

template<class Type>
VECTOR_2D<Type> &
operator-(const VECTOR_2D<Type> vector)
{
	vector.x = -vector.x;
	vector.y = -vector.y;
	return vector;
}


template <class Type>
VECTOR_2D<Type> &
operator+(const VECTOR_2D<Type> lvector, const VECTOR_2D<Type>& rvector)
{
	lvector.x += rvector.x;
	lvector.y += rvector.y;
	return lvector;
}

template <class Type>
VECTOR_2D<Type> &
operator-(const VECTOR_2D<Type> lvector, const VECTOR_2D<Type> &rvector)
{
	lvector.x = lvector.x - rvector.x;
	lvector.y = lvector.y - rvector.y;
	return lvector;
}


// Here "*" means product of vectors
template <class Type>
Type &
operator*(const VECTOR_2D<Type>& lvector, const VECTOR_2D<Type>& rvector)
{
	Type tmp;
	tmp = lvector.x * rvector.x + lvector.y * rvector.y;
	return tmp;
}

// Here "*" means scaling by scalar
template <class Type>
Type &
operator*(const VECTOR_2D<Type> lvector, const Type& rvalue)
{
	lvector.x = lvector.x * rvalue;
	lvector.y = lvector.y * rvalue;
	return lvector;
}

// Here "*" means scaling by scalar
template <class Type>
Type &
operator*(const Type& lvalue, const VECTOR_2D<Type> rvector)
{
	rvector.x = rvector.x * lvalue;
	rvector.y = rvector.y * lvalue;
	return rvector;
}


namespace Vector_2D {
	template<class Type>
	double
	norm(const VECTOR_2D<Type> &vector)
	{
		return sqrt((double)vector.x * vector.x + vector.y * vector.y);
	}

	template<class Type>
	double
	arg(const VECTOR_2D<Type> &vector)
	{
		return atan2(vector.y, vector.x);
	}
}

