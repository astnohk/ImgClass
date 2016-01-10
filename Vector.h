#ifndef LIB_ImgClass_Vector
#define LIB_ImgClass_Vector

#include <ostream>

template <class T>
struct VECTOR_2D
{
	public:
		T x;
		T y;

		// Constructors
		VECTOR_2D(void);
		VECTOR_2D(const VECTOR_2D<T>& value);
		VECTOR_2D(const T& init_x, const T& init_y);
		void reset(void);
		void reset(const T& init_x, const T& init_y);

		// Operators
		template<class RT> operator VECTOR_2D<RT>() const;

		VECTOR_2D<T>& operator=(const VECTOR_2D<T>& vector);
		template<class Tval> VECTOR_2D<T>& operator=(const Tval& value);

		VECTOR_2D<T>& operator+=(const VECTOR_2D<T>& vector);

		VECTOR_2D<T>& operator-=(const VECTOR_2D<T>& vector);

		VECTOR_2D<T>& operator*=(const VECTOR_2D<T>& vector);
		VECTOR_2D<T>& operator*=(const T& value);

		VECTOR_2D<T>& operator/=(const T& value);

		bool operator==(const VECTOR_2D<T>& vector);
		bool operator!=(const VECTOR_2D<T>& vector);

		// Friend operators
		template<class Type> friend VECTOR_2D<Type> operator+(VECTOR_2D<Type> vector);
		template<class Type> friend VECTOR_2D<Type> operator-(VECTOR_2D<Type> vector);
                                            
		template<class Type> friend VECTOR_2D<Type> operator+(VECTOR_2D<Type> lvector, const VECTOR_2D<Type>& rvector);
		template<class Type> friend VECTOR_2D<Type> operator-(VECTOR_2D<Type> lvector, const VECTOR_2D<Type>& rvector);
                                            
		template<class Type> friend Type operator*(const VECTOR_2D<Type>& lvector, const VECTOR_2D<Type>& rvector);
		template<class Type, class Tval> friend VECTOR_2D<Type> operator*(VECTOR_2D<Type> lvector, const Tval& rvalue);
		template<class Type, class Tval> friend VECTOR_2D<Type> operator*(const Tval& lvalue, VECTOR_2D<Type> rvector);

		template<class Type, class Tval> friend VECTOR_2D<Type> operator/(VECTOR_2D<Type> lvector, const Tval& rvalue);

		template<class Type> friend double norm(const VECTOR_2D<Type> &vector);
		template<class Type> friend double arg(const VECTOR_2D<Type> &vector);
		template<class Type> friend VECTOR_2D<Type> floor(const VECTOR_2D<Type> &vector);
		template<class Type> friend VECTOR_2D<Type> round(const VECTOR_2D<Type> &vector);
		template<class Type> friend VECTOR_2D<Type> ceil(const VECTOR_2D<Type> &vector);
};

// Global operator overloading
template<class Type> VECTOR_2D<Type> operator+(VECTOR_2D<Type> vector);
template<class Type> VECTOR_2D<Type> operator-(VECTOR_2D<Type> vector);

template<class Type> VECTOR_2D<Type> operator+(VECTOR_2D<Type> lvector, const VECTOR_2D<Type>& rvector);
template<class Type> VECTOR_2D<Type> operator-(VECTOR_2D<Type> lvector, const VECTOR_2D<Type>& rvector);

template<class Type> Type operator*(const VECTOR_2D<Type>& lvector, const VECTOR_2D<Type>& rvector);
template<class Type, class Tval> VECTOR_2D<Type> operator*(VECTOR_2D<Type> lvector, const Tval& rvalue);
template<class Type, class Tval> VECTOR_2D<Type> operator*(const Tval& lvalue, VECTOR_2D<Type> rvector);

template<class Type, class Tval> VECTOR_2D<Type> operator/(VECTOR_2D<Type> lvector, const Tval& rvalue);

// Arithmetic
template<class Type> double norm_squared(const VECTOR_2D<Type>& vector);
template<class Type> double norm(const VECTOR_2D<Type>& vector);
template<class Type> double arg(const VECTOR_2D<Type>& vector);
template<class Type> VECTOR_2D<Type> floor(const VECTOR_2D<Type> &vector);
template<class Type> VECTOR_2D<Type> round(const VECTOR_2D<Type> &vector);
template<class Type> VECTOR_2D<Type> ceil(const VECTOR_2D<Type> &vector);

// Stream
template<class Type> std::ostream& operator<<(std::ostream& os, const VECTOR_2D<Type>& rvalue);




template <class T>
struct Vector_ST
{
	public:
		T x;
		T y;
		T t;

		// Constructors
		Vector_ST(void);
		Vector_ST(const Vector_ST<T>& value);
		Vector_ST(const T& init_x, const T& init_y, const int init_t);
		void clear(void);
		void reset(void);
		void reset(const T& init_x, const T& init_y, const int init_t);

		// Operators
		Vector_ST<T>& operator=(const Vector_ST<T>& vector);

		Vector_ST<T>& operator+=(const Vector_ST<T>& vector);

		Vector_ST<T>& operator-=(const Vector_ST<T>& vector);

		Vector_ST<T>& operator*=(const Vector_ST<T>& vector);
		Vector_ST<T>& operator*=(const T& value);

		Vector_ST<T>& operator/=(const T& value);

		bool operator==(const Vector_ST<T>& vector);
		bool operator!=(const Vector_ST<T>& vector);

		// Friend operators
		template<class Type> friend Vector_ST<Type> operator+(Vector_ST<Type> vector);
		template<class Type> friend Vector_ST<Type> operator-(Vector_ST<Type> vector);
                                            
		template<class Type> friend Vector_ST<Type> operator+(Vector_ST<Type> lvector, const Vector_ST<Type>& rvector);
		template<class Type> friend Vector_ST<Type> operator-(Vector_ST<Type> lvector, const Vector_ST<Type>& rvector);
                                            
		template<class Type> friend Type operator*(const Vector_ST<Type>& lvector, const Vector_ST<Type>& rvector);
		template<class Type, class Tval> friend Vector_ST<Type> operator*(Vector_ST<Type> lvector, const Tval& rvalue);
		template<class Type, class Tval> friend Vector_ST<Type> operator*(const Tval& lvalue, Vector_ST<Type> rvector);

		template<class Type, class Tval> friend Vector_ST<Type> operator/(Vector_ST<Type> lvector, const Tval& rvalue);

		template<class Type> friend double norm(const Vector_ST<Type> &vector);
		template<class Type> friend double arg(const Vector_ST<Type> &vector);
		template<class Type> friend Vector_ST<Type> floor(const Vector_ST<Type> &vector);
		template<class Type> friend Vector_ST<Type> round(const Vector_ST<Type> &vector);
		template<class Type> friend Vector_ST<Type> ceil(const Vector_ST<Type> &vector);
};


// Global operator overloading
template<class Type> Vector_ST<Type> operator+(Vector_ST<Type> vector);
template<class Type> Vector_ST<Type> operator-(Vector_ST<Type> vector);

template<class Type> Vector_ST<Type> operator+(Vector_ST<Type> lvector, const Vector_ST<Type>& rvector);
template<class Type> Vector_ST<Type> operator-(Vector_ST<Type> lvector, const Vector_ST<Type>& rvector);

template<class Type> Type operator*(const Vector_ST<Type>& lvector, const Vector_ST<Type>& rvector);
template<class Type, class Tval> Vector_ST<Type> operator*(Vector_ST<Type> lvector, const Tval& rvalue);
template<class Type, class Tval> Vector_ST<Type> operator*(const Tval& lvalue, Vector_ST<Type> rvector);

template<class Type, class Tval> Vector_ST<Type> operator/(Vector_ST<Type> lvector, const Tval& rvalue);

// Arithmetic
template<class Type> double norm_squared(const Vector_ST<Type>& vector);
template<class Type> double norm(const Vector_ST<Type>& vector);
template<class Type> double arg(const Vector_ST<Type>& vector);
template<class Type> Vector_ST<Type> floor(const Vector_ST<Type> &vector);
template<class Type> Vector_ST<Type> round(const Vector_ST<Type> &vector);
template<class Type> Vector_ST<Type> ceil(const Vector_ST<Type> &vector);

// Stream
template<class Type> std::ostream& operator<<(std::ostream& os, const Vector_ST<Type>& rvalue);

#include "Vector_private.h"
#include "Vector_ST_private.h"

#endif

