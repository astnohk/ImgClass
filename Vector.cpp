#include "Vector.h"



template <class T>
VECTOR_2D<T>::VECTOR_2D(void)
{
	x = T();
	y = T();
}

template <class T>
VECTOR_2D::VECTOR_2D(const T& init_x, const T& init_y)
{
	x = init_x;
	y = init_y;
}

template <class T>
void
VECTOR_2D::reset(void)
{
	x = T();
	y = T();
}

template <class T>
void
VECTOR_2D::reset(const T& init_x, const T& init_y)
{
	x = init_x;
	y = init_y;
}


template <class T>
VECTOR_2D<T> &
VECTOR_2D::operator+=(const VECTOR_2D<T>& vector)
{
	if (this != &vector) {
		this->x += vector.x;
		this->y += vector.y;
	}
	return *this;
}

template <class T>
VECTOR_2D<T> &
VECTOR_2D::operator-=(const VECTOR_2D& vector)
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
	if (this != &vector) {
		this->x *= value;
		this->y *= value;
	}
	return *this;
}


template <class T>
bool
VECTOR_2D<T>::operator==(const T& vector)
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
VECTOR_2D<T>::operator!=(const T& vector)
{
	if (this->x == vector.x
	    && this->y == vector.y) {
		return false;
	} else {
		return true;
	}
}




// Global operators
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
operator-(const VECTOR_2D lvector, const VECTOR_2D &rvector)
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

