#ifndef LIB_ImgClass_RGB
#define LIB_ImgClass_RGB

namespace ImgClass {
	template <class T>
	class RGB
	{
		public:
			T R;
			T G;
			T B;

			// Constructor
			RGB(void);
			RGB(const T& red, const T& green, const T& blue);
			RGB(const RGB<T>& color); // Copy constructor
			explicit RGB(const T& value);

			// Operators
			template<class ConvertType> operator RGB<ConvertType>() const;
			template<class ConvertType> operator ConvertType() const; // return intensity with YUV coefficient

			template<class RT> RGB<T>& operator=(const RGB<RT>& rvalue);
			RGB<T>& operator=(const T& rvalue);

			RGB<T>& operator+=(const RGB<T>& rcolor);
			RGB<T>& operator+=(const T& rvalue);

			RGB<T>& operator-=(const RGB<T>& rcolor);
			RGB<T>& operator-=(const T& rvalue);

			RGB<T>& operator*=(const RGB<T>& rcolor);
			RGB<T>& operator*=(const T& rvalue);

			RGB<T>& operator/=(const RGB<T>& rcolor);
			RGB<T>& operator/=(const T& rvalue);

			// friend Global Operators
			// Arithmetic
			template<class Type> friend const RGB<Type> operator+(RGB<Type> color);
			template<class Type> friend const RGB<Type> operator-(RGB<Type> color);

			template<class Type> friend const RGB<Type> operator+(RGB<Type> lcolor, const RGB<Type>& rcolor);
			template<class Type> friend const RGB<Type> operator+(RGB<Type> lcolor, const Type& rvalue);
			template<class Type> friend const RGB<Type> operator+(const Type& lvalue, RGB<Type> rcolor);

			template<class Type> friend const RGB<Type> operator-(RGB<Type> lcolor, const RGB<Type>& rcolor);
			template<class Type> friend const RGB<Type> operator-(RGB<Type> lcolor, const Type& rvalue);
			template<class Type> friend const RGB<Type> operator-(const Type& lvalue, RGB<Type> rcolor);

			template<class Type> friend const RGB<Type> operator*(RGB<Type> lcolor, const RGB<Type>& rcolor);
			template<class Type> friend const RGB<Type> operator*(RGB<Type> lcolor, const Type& rvalue);
			template<class Type> friend const RGB<Type> operator*(const Type& lvalue, RGB<Type> rcolor);

			template<class Type> friend const RGB<Type> operator/(RGB<Type> lcolor, const RGB<Type>& rcolor);
			template<class Type> friend const RGB<Type> operator/(RGB<Type> lcolor, const Type& rvalue);
			template<class Type> friend const RGB<Type> operator/(const Type& lvalue, RGB<Type> rcolor);

			// Comparator
			template<class Type> friend bool operator==(const RGB<Type>& lcolor, const RGB<Type>& rcolor);
			template<class Type> friend bool operator!=(const RGB<Type>& lcolor, const RGB<Type>& rcolor);
	};
}

// Global Operators
// Arithmetic
template<class Type> const ImgClass::RGB<Type>& operator+(ImgClass::RGB<Type> color);
template<class Type> const ImgClass::RGB<Type>& operator-(ImgClass::RGB<Type> color);

template<class Type> const ImgClass::RGB<Type>& operator+(ImgClass::RGB<Type> lcolor, const ImgClass::RGB<Type>& rcolor);
template<class Type> const ImgClass::RGB<Type>& operator+(ImgClass::RGB<Type> lcolor, const Type& rvalue);
template<class Type> const ImgClass::RGB<Type>& operator+(const Type& lvalue, ImgClass::RGB<Type> rcolor);

template<class Type> const ImgClass::RGB<Type>& operator-(ImgClass::RGB<Type> lcolor, const ImgClass::RGB<Type>& rcolor);
template<class Type> const ImgClass::RGB<Type>& operator-(ImgClass::RGB<Type> lcolor, const Type& rvalue);
template<class Type> const ImgClass::RGB<Type>& operator-(const Type& lvalue, ImgClass::RGB<Type> rcolor);

template<class Type> const ImgClass::RGB<Type>& operator*(ImgClass::RGB<Type> lcolor, const ImgClass::RGB<Type>& rcolor);
template<class Type> const ImgClass::RGB<Type>& operator*(ImgClass::RGB<Type> lcolor, const Type& rvalue);
template<class Type> const ImgClass::RGB<Type>& operator*(const Type& lvalue, ImgClass::RGB<Type> rcolor);

template<class Type> const ImgClass::RGB<Type>& operator/(ImgClass::RGB<Type> lcolor, const ImgClass::RGB<Type>& rcolor);
template<class Type> const ImgClass::RGB<Type>& operator/(ImgClass::RGB<Type> lcolor, const Type& rvalue);
template<class Type> const ImgClass::RGB<Type>& operator/(const Type& lvalue, ImgClass::RGB<Type> rcolor);

// Comparator
template<class Type> bool operator==(const ImgClass::RGB<Type>& lcolor, const ImgClass::RGB<Type>& rcolor);
template<class Type> bool operator!=(const ImgClass::RGB<Type>& lcolor, const ImgClass::RGB<Type>& rcolor);

#include "RGB_private.h"

#endif

