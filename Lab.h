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
			Lab(const RGB& color);

			Lab& set(const RGB& color);
			// Operators
			explicit operator double() const;

			Lab& operator=(const Lab& value);
			Lab& operator=(const double& value);

			Lab& operator+=(const Lab& color);
			Lab& operator+=(const double& value);

			Lab& operator-=(const Lab& color);
			Lab& operator-=(const double& value);

			Lab& operator*=(const Lab& color);
			Lab& operator*=(const double& value);

			Lab& operator/=(const Lab& color);
			Lab& operator/=(const double& value);

		protected:
			// Converter
			double f(const double t);
	};
}

// Global Operators
// Arithmetic
const ImgClass::Lab operator+(ImgClass::Lab color);
const ImgClass::Lab operator-(ImgClass::Lab color);

// Comparator
bool operator==(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor);
bool operator==(const ImgClass::Lab& lcolor, const double& rvalue);
bool operator==(const double& lvalue, const ImgClass::Lab& rcolor);

bool operator!=(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor);
bool operator!=(const ImgClass::Lab& lcolor, const double& rvalue);
bool operator!=(const double& lvalue, const ImgClass::Lab& rcolor);

bool operator<(const ImgClass::Lab& lcolor, const double& rvalue);
bool operator<(const double& lvalue, const ImgClass::Lab& rcolor);

bool operator>(const ImgClass::Lab& lcolor, const double& rvalue);
bool operator>(const double& lvalue, const ImgClass::Lab& rcolor);

// Arithmetic non-substituting operator
ImgClass::Lab operator+(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor);
ImgClass::Lab operator+(const ImgClass::Lab& lcolor, const double& rvalue);
ImgClass::Lab operator+(const double& lvalue, const ImgClass::Lab& rcolor);

ImgClass::Lab operator-(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor);
ImgClass::Lab operator-(const ImgClass::Lab& lcolor, const double& rvalue);
ImgClass::Lab operator-(const double& lvalue, const ImgClass::Lab& rcolor);

ImgClass::Lab operator*(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor);
ImgClass::Lab operator*(const ImgClass::Lab& lcolor, const double& rvalue);
ImgClass::Lab operator*(const double& lvalue, const ImgClass::Lab& rcolor);

ImgClass::Lab operator/(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor);
ImgClass::Lab operator/(const ImgClass::Lab& lcolor, const double& rvalue);

// Product
double inner_prod(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor);

// Norm
const ImgClass::Lab abs(const ImgClass::Lab& color);
double norm_squared(const ImgClass::Lab& color);
double norm(const ImgClass::Lab& color);

#endif

