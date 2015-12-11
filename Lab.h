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
			Lab(const Lab& color); // Copy constructor
			Lab(const RGB<double>& color);

			// Operators
			explicit operator double() const;

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

