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
template<class Type> double norm(const VECTOR_2D<Type>& vector);
template<class Type> double arg(const VECTOR_2D<Type>& vector);
template<class Type> VECTOR_2D<Type> floor(const VECTOR_2D<Type> &vector);
template<class Type> VECTOR_2D<Type> round(const VECTOR_2D<Type> &vector);
template<class Type> VECTOR_2D<Type> ceil(const VECTOR_2D<Type> &vector);

// Stream
template<class Type> std::ostream& operator<<(std::ostream& os, const VECTOR_2D<Type>& rvalue);

#include "Vector_private.h"

#endif

