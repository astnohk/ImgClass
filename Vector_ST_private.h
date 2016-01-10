#include <cmath>


template <class T>
Vector_ST<T>::Vector_ST(void)
{
	x = T();
	y = T();
	t = 0;
}

template <class T>
Vector_ST<T>::Vector_ST(const Vector_ST<T>& value)
{
	x = value.x;
	y = value.y;
	t = value.t;
}

template <class T>
Vector_ST<T>::Vector_ST(const T& init_x, const T& init_y, const int init_t)
{
	x = init_x;
	y = init_y;
	t = init_t;
}

template <class T>
void
Vector_ST<T>::reset(void)
{
	x = T();
	y = T();
	t = 0;
}

template <class T>
void
Vector_ST<T>::reset(const T& init_x, const T& init_y, const int init_t)
{
	x = init_x;
	y = init_y;
	t = init_t;
}




// ----- Operators -----
template <class T>
Vector_ST<T> &
Vector_ST<T>::operator=(const Vector_ST<T>& vector)
{
	x = vector.x;
	y = vector.y;
	t = vector.t;
	return *this;
}


template <class T>
Vector_ST<T> &
Vector_ST<T>::operator+=(const Vector_ST<T>& vector)
{
	if (this != &vector) {
		this->x += vector.x;
		this->y += vector.y;
	}
	return *this;
}

template <class T>
Vector_ST<T> &
Vector_ST<T>::operator-=(const Vector_ST& vector)
{
	if (this != &vector) {
		this->x -= vector.x;
		this->y -= vector.y;
	}
	return *this;
}

template <class T>
Vector_ST<T> &
Vector_ST<T>::operator*=(const T& value) // Scalar multiplication
{
	this->x *= value;
	this->y *= value;
	return *this;
}

template <class T>
Vector_ST<T> &
Vector_ST<T>::operator/=(const T& value) // Scalar multiplication
{
	this->x /= value;
	this->y /= value;
	return *this;
}


template <class T>
bool
Vector_ST<T>::operator==(const Vector_ST<T>& vector)
{
	if (this->x == vector.x
	    && this->y == vector.y
	    && this->t == vector.t) {
		return true;
	} else {
		return false;
	}
}

template <class T>
bool
Vector_ST<T>::operator!=(const Vector_ST<T>& vector)
{
	if (this->x == vector.x
	    && this->y == vector.y
	    && this->t == vector.t) {
		return false;
	} else {
		return true;
	}
}




//
// ----- Global operators -----
//
template <class Type>
Vector_ST<Type>
operator+(const Vector_ST<Type> vector)
{
	return vector;
}

template <class Type>
Vector_ST<Type>
operator-(Vector_ST<Type> vector)
{
	vector.x = -vector.x;
	vector.y = -vector.y;
	return vector;
}


template <class Type>
Vector_ST<Type>
operator+(Vector_ST<Type> lvector, const Vector_ST<Type>& rvector)
{
	lvector.x += rvector.x;
	lvector.y += rvector.y;
	return lvector;
}

template <class Type>
Vector_ST<Type>
operator-(Vector_ST<Type> lvector, const Vector_ST<Type>& rvector)
{
	lvector.x = lvector.x - rvector.x;
	lvector.y = lvector.y - rvector.y;
	return lvector;
}


// Here "*" means product of vectors
template <class Type>
Type
operator*(const Vector_ST<Type>& lvector, const Vector_ST<Type>& rvector)
{
	Type tmp;
	tmp = lvector.x * rvector.x + lvector.y * rvector.y;
	return tmp;
}

// Here "*" means scaling by scalar
template <class Type, class Tval>
Vector_ST<Type>
operator*(Vector_ST<Type> lvector, const Tval& rvalue)
{
	lvector.x = lvector.x * rvalue;
	lvector.y = lvector.y * rvalue;
	return lvector;
}

// Here "*" means scaling by scalar
template <class Type, class Tval>
Vector_ST<Type>
operator*(const Tval& lvalue, Vector_ST<Type> rvector)
{
	rvector.x = rvector.x * lvalue;
	rvector.y = rvector.y * lvalue;
	return rvector;
}

// Here "/" means scaling by scalar
template <class Type, class Tval>
Vector_ST<Type>
operator/(Vector_ST<Type> lvector, const Tval& rvalue)
{
	lvector.x = lvector.x / rvalue;
	lvector.y = lvector.y / rvalue;
	return lvector;
}


template <class Type>
double
norm_squared(const Vector_ST<Type> &value)
{
	return value.x * value.x + value.y * value.y;
}

template <class Type>
double
norm(const Vector_ST<Type> &value)
{
	return sqrt(value.x * value.x + value.y * value.y);
}

template <class Type>
double
arg(const Vector_ST<Type> &vector)
{
	return atan2(vector.y, vector.x);
}

template <class Type>
Vector_ST<Type>
floor(Vector_ST<Type> vector)
{
	vector.x = floor(vector.x);
	vector.y = floor(vector.y);
	return vector;
}

template <class Type>
Vector_ST<Type>
round(Vector_ST<Type> vector)
{
	vector.x = round(vector.x);
	vector.y = round(vector.y);
	return vector;
}

template <class Type>
Vector_ST<Type>
ceil(Vector_ST<Type> vector)
{
	vector.x = ceil(vector.x);
	vector.y = ceil(vector.y);
	return vector;
}


// Stream
template <class Type>
std::ostream &
operator<<(std::ostream& os, const Vector_ST<Type>& rvalue)
{
	os << "[x:" << rvalue.x << " y:" << rvalue.y << " t:" << rvalue.t << "]";
	return os;
}

