#ifndef LIB_Vector
#define LIB_Vector

template <class T>
struct VECTOR_2D
{
	public:
		T x;
		T y;

		// Constructors
		VECTOR_2D(void);
		VECTOR_2D(const T& init_x, const T& init_y);
		void reset(void);
		void reset(const T& init_x, const T& init_y);

		// Operators
		VECTOR_2D<T>& operator+=(const VECTOR_2D<T>& vector);
		VECTOR_2D<T>& operator-=(const VECTOR_2D<T>& vector);
		VECTOR_2D<T>& operator*=(const VECTOR_2D<T>& vector);
		VECTOR_2D<T>& operator*=(const T& value);

		bool operator==(const VECTOR_2D<T>& vector);
		bool operator!=(const VECTOR_2D<T>& vector);

		// Friend operators
		template<class Type> friend VECTOR_2D<Type>& operator+(const VECTOR_2D<Type> vector);
		template<class Type> friend VECTOR_2D<Type>& operator-(const VECTOR_2D<Type> vector);
                                            
		template<class Type> friend VECTOR_2D<Type>& operator+(const VECTOR_2D<Type> lvector, const VECTOR_2D<Type>& rvector);
		template<class Type> friend VECTOR_2D<Type>& operator-(const VECTOR_2D<Type> lvector, const VECTOR_2D<Type>& rvector);
                                            
		template<class Type> friend Type& operator*(const VECTOR_2D<Type> lvector, const VECTOR_2D<Type>& rvector);
		template<class Type> friend Type& operator*(const VECTOR_2D<Type> lvector, const Type& rvalue);
		template<class Type> friend Type& operator*(const Type& lvalue, const VECTOR_2D<Type> rvector);

		template<class Type> friend Type& norm(const VECTOR_2D<Type> &vector);
};


// Global operator overloading
template<class Type> VECTOR_2D<Type>& operator+(const VECTOR_2D<Type> vector);
template<class Type> VECTOR_2D<Type>& operator-(const VECTOR_2D<Type> vector);

template<class Type> VECTOR_2D<Type>& operator+(const VECTOR_2D<Type> lvector, const VECTOR_2D<Type>& rvector);
template<class Type> VECTOR_2D<Type>& operator-(const VECTOR_2D<Type> lvector, const VECTOR_2D<Type>& rvector);

template<class Type> Type& operator*(const VECTOR_2D<Type> lvector, const VECTOR_2D<Type>& rvector);
template<class Type> Type& operator*(const VECTOR_2D<Type> lvector, const Type& rvalue);
template<class Type> Type& operator*(const Type& lvalue, const VECTOR_2D<Type> rvector);

template<class Type> double norm(const VECTOR_2D<Type> &vector);
template<class Type> double arg(const VECTOR_2D<Type> &vector);

#include "Vector_private.h"

#endif

