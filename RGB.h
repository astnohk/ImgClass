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

			// Operators
			explicit operator T() const; // return intensity
			template<class ConvertType> operator RGB<ConvertType>() const;

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

			// Non substituting
			const RGB<T> operator+(const RGB<T>& rcolor) const;
			const RGB<T> operator+(const T& rvalue) const;

			const RGB<T> operator-(const RGB<T>& rcolor) const;
			const RGB<T> operator-(const T& rvalue) const;
	};
}

#include "RGB_private.h"

#endif

