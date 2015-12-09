#ifndef LIB_ImgClass_Lab
#define LIB_ImgClass_Lab

#include "RGB.h"

/* L*a*b* color space
 * constructor will map RGB to L*a*b* color space
 */
namespace ImgClass {
	class Lab
	{
		public:
			double L;
			double a;
			double b;

			// Constructor
			Lab(void);
			Lab(const double& _L, const double& _a, const double& _b);
			Lab(const RGB<double>& value); // Copy constructor
			explicit Lab(const double& value);

			// Operators
			Lab& operator=(const Lab& value);
			Lab& operator=(const double& value);
			Lab& operator=(const RGB<double>& value);

			Lab& operator+=(const Lab& color);
			Lab& operator+=(const double& value);

			Lab& operator-=(const Lab& color);
			Lab& operator-=(const double& value);

			Lab& operator*=(const Lab& color);
			Lab& operator*=(const double& value);

			Lab& operator/=(const Lab& color);
			Lab& operator/=(const double& value);

			// Arithmetic non-substituting operator
			Lab operator+(const Lab& rcolor) const;
			Lab operator+(const double& rvalue) const;

			Lab operator-(const Lab& rcolor) const;
			Lab operator-(const double& rvalue) const;

			Lab operator*(const Lab& rcolor) const;
			Lab operator*(const double& rvalue) const;

			Lab operator/(const Lab& rcolor) const;
			Lab operator/(const double& rvalue) const;

			// friend Global Operators
			// Arithmetic
			friend const Lab operator+(Lab color);
			friend const Lab operator-(Lab color);

			// Comparator
			friend bool operator==(const Lab& lcolor, const Lab& rcolor);
			friend bool operator!=(const Lab& lcolor, const Lab& rcolor);

			// Norm
			friend double norm(const ImgClass::Lab& color);
	};
}

// Global Operators
// Arithmetic
const ImgClass::Lab operator+(ImgClass::Lab color);
const ImgClass::Lab operator-(ImgClass::Lab color);

// Comparator
bool operator==(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor);
bool operator!=(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor);

// Norm
double norm_squared(const ImgClass::Lab& color);

#endif

