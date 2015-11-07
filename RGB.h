#ifndef nullptr
#define nullptr 0
#endif



#ifndef LIB_ImgClass_RGB
#define LIB_ImgClass_RGB

namespace ImgClass {
	template <class T>
	class RGB
	{
		private:
			T _red;
			T _green;
			T _blue;
		public:
			RGB(void);
			RGB(const T& red, const T& green, const T& blue);
			template<class RT> RGB(const RT& value);

			// Data access
			T R(void);
			T G(void);
			T B(void);

			// Operators
			template<class RT> RGB<T>& operator=(const RGB<RT>& value);
			template<class RT> RGB<T>& operator=(const RT& value);

			RGB<T>& operator+=(const RGB<T>& color);
			RGB<T>& operator+=(const T& value);

			RGB<T>& operator-=(const RGB<T>& color);
			RGB<T>& operator-=(const T& value);

			RGB<T>& operator*=(const RGB<T>& color);
			RGB<T>& operator*=(const T& value);

			RGB<T>& operator/=(const RGB<T>& color);
			RGB<T>& operator/=(const T& value);

			// friend Global Operators
			// Arithmetic
			template<class Type> friend RGB<Type>& operator+(RGB<Type> color);
			template<class Type> friend RGB<Type>& operator-(RGB<Type> color);

			template<class Type> friend RGB<Type>& operator+(RGB<Type> lcolor, const RGB<Type>& rcolor);
			template<class Type> friend RGB<Type>& operator+(RGB<Type> lcolor, const Type& rvalue);
			template<class Type> friend RGB<Type>& operator+(const Type& lvalue, RGB<Type> rcolor);

			template<class Type> friend RGB<Type>& operator-(RGB<Type> lcolor, const RGB<Type>& rcolor);
			template<class Type> friend RGB<Type>& operator-(RGB<Type> lcolor, const Type& rvalue);
			template<class Type> friend RGB<Type>& operator-(const Type& lvalue, RGB<Type> rcolor);

			template<class Type> friend RGB<Type>& operator*(RGB<Type> lcolor, const RGB<Type>& rcolor);
			template<class Type> friend RGB<Type>& operator*(RGB<Type> lcolor, const Type& rvalue);
			template<class Type> friend RGB<Type>& operator*(const Type& lvalue, RGB<Type> rcolor);

			template<class Type> friend RGB<Type>& operator/(RGB<Type> lcolor, const RGB<Type>& rcolor);
			template<class Type> friend RGB<Type>& operator/(RGB<Type> lcolor, const Type& rvalue);
			template<class Type> friend RGB<Type>& operator/(const Type& lvalue, RGB<Type> rcolor);

			// Comparator
			template<class Type> friend bool operator==(const RGB<Type>& lcolor, const RGB<Type>& rcolor);
			template<class Type> friend bool operator!=(const RGB<Type>& lcolor, const RGB<Type>& rcolor);
	};
}

// Global Operators
// Arithmetic
template<class Type> ImgClass::RGB<Type>& operator+(ImgClass::RGB<Type> color);
template<class Type> ImgClass::RGB<Type>& operator-(ImgClass::RGB<Type> color);

template<class Type> ImgClass::RGB<Type>& operator+(ImgClass::RGB<Type> lcolor, const ImgClass::RGB<Type>& rcolor);
template<class Type> ImgClass::RGB<Type>& operator+(ImgClass::RGB<Type> lcolor, const Type& rvalue);
template<class Type> ImgClass::RGB<Type>& operator+(const Type& lvalue, ImgClass::RGB<Type> rcolor);

template<class Type> ImgClass::RGB<Type>& operator-(ImgClass::RGB<Type> lcolor, const ImgClass::RGB<Type>& rcolor);
template<class Type> ImgClass::RGB<Type>& operator-(ImgClass::RGB<Type> lcolor, const Type& rvalue);
template<class Type> ImgClass::RGB<Type>& operator-(const Type& lvalue, ImgClass::RGB<Type> rcolor);

template<class Type> ImgClass::RGB<Type>& operator*(ImgClass::RGB<Type> lcolor, const ImgClass::RGB<Type>& rcolor);
template<class Type> ImgClass::RGB<Type>& operator*(ImgClass::RGB<Type> lcolor, const Type& rvalue);
template<class Type> ImgClass::RGB<Type>& operator*(const Type& lvalue, ImgClass::RGB<Type> rcolor);

template<class Type> ImgClass::RGB<Type>& operator/(ImgClass::RGB<Type> lcolor, const ImgClass::RGB<Type>& rcolor);
template<class Type> ImgClass::RGB<Type>& operator/(ImgClass::RGB<Type> lcolor, const Type& rvalue);
template<class Type> ImgClass::RGB<Type>& operator/(const Type& lvalue, ImgClass::RGB<Type> rcolor);

// Comparator
template<class Type> bool operator==(const ImgClass::RGB<Type>& lcolor, const ImgClass::RGB<Type>& rcolor);
template<class Type> bool operator!=(const ImgClass::RGB<Type>& lcolor, const ImgClass::RGB<Type>& rcolor);

#include "RGB_private.h"

#endif

