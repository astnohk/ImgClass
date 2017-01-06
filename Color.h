#ifndef LIB_ImgClass_COLOR
#define LIB_ImgClass_COLOR

#include <ostream>

namespace ImgClass {
	class RGB;
	class HSV;
	class Lab;
}
// --------------------------------


namespace ImgClass {
	class RGB
	{
		public:

		double R;
		double G;
		double B;

		// Constructor
		RGB(void);
		RGB(const double& red, const double& green, const double& blue);
		RGB(const RGB& rgb); // Copy constructor
		RGB(const HSV& hsv); // Conversion
		RGB(const Lab& lab); // Conversion

		RGB& set(const double& red, const double& green, const double& blue);
		RGB& gamma(const double& gamma_val);

		// Operators
		explicit operator double() const; // return intensity

		RGB& operator=(const RGB& rvalue);
		RGB& operator=(const double& rvalue);

		RGB& operator+=(const RGB& rcolor);
		RGB& operator+=(const double& rvalue);

		RGB& operator-=(const RGB& rcolor);
		RGB& operator-=(const double& rvalue);

		RGB& operator*=(const RGB& rcolor);
		RGB& operator*=(const double& rvalue);

		RGB& operator/=(const RGB& rcolor);
		RGB& operator/=(const double& rvalue);
	};
}


const ImgClass::RGB operator+(ImgClass::RGB rcolor);
const ImgClass::RGB operator-(ImgClass::RGB rcolor);

const ImgClass::RGB operator+(const ImgClass::RGB& lcolor, const ImgClass::RGB& rcolor);
const ImgClass::RGB operator-(const ImgClass::RGB& lcolor, const ImgClass::RGB& rcolor);

const ImgClass::RGB operator*(const ImgClass::RGB& lcolor, const ImgClass::RGB& rcolor);
const ImgClass::RGB operator*(const ImgClass::RGB& lcolor, const double& rvalue);
const ImgClass::RGB operator*(const double& lvalue, const ImgClass::RGB& rcolor);

const ImgClass::RGB operator/(const ImgClass::RGB& lcolor, const ImgClass::RGB& rcolor);
const ImgClass::RGB operator/(const ImgClass::RGB& lcolor, const double& rvalue);

// Comparator
bool operator==(const ImgClass::RGB& lcolor, const ImgClass::RGB& rcolor);
bool operator!=(const ImgClass::RGB& lcolor, const ImgClass::RGB& rcolor);

// Product
double inner_prod(const ImgClass::RGB& lcolor, const ImgClass::RGB& rcolor);

// Norm
double norm_squared(const ImgClass::RGB& color);
double norm(const ImgClass::RGB& color);

// Saturation
ImgClass::RGB saturate(const ImgClass::RGB& value, const double& min, const double& max);

// Quantization
ImgClass::RGB color_quantize(const ImgClass::RGB &value, const double &max = 255.0);

// Stream
std::ostream& operator<<(std::ostream& os, const ImgClass::RGB& rcolor);




/******** HSV Color Space ********/
namespace ImgClass {
	class HSV
	{
		public:

		double H;
		double S;
		double V;

		// Constructor
		HSV(void);
		HSV(const double& hue, const double& saturation, const double& value);
		HSV(const HSV& hsv); // Copy constructor
		HSV(const RGB& rgb);

		HSV& set(const double& hue, const double& saturation, const double& value);
		HSV& set(const RGB& rgb);
		HSV& set_H(const double val);
		HSV& set_S(const double val);
		HSV& set_V(const double val);

		// Operators
		HSV& operator=(const HSV& rval);
		HSV& operator+=(const HSV& rval);
		HSV& operator-=(const HSV& rval);
	};
}


const ImgClass::HSV operator+(ImgClass::HSV rval);
const ImgClass::HSV operator-(ImgClass::HSV rval);

const ImgClass::HSV operator+(const ImgClass::HSV& lval, const ImgClass::HSV& rval);
const ImgClass::HSV operator-(const ImgClass::HSV& lval, const ImgClass::HSV& rval);

bool operator==(const ImgClass::HSV& lval, const ImgClass::HSV& rval);
bool operator!=(const ImgClass::HSV& lval, const ImgClass::HSV& rval);

std::ostream& operator<<(std::ostream& os, const ImgClass::HSV& rval);




/******** L*a*b* color space ********/
// constructor will map RGB to L*a*b* color space
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
bool operator!=(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor);

bool operator<(const ImgClass::Lab& lcolor, const double& rvalue);
bool operator<(const double& lvalue, const ImgClass::Lab& rcolor);

bool operator>(const ImgClass::Lab& lcolor, const double& rvalue);
bool operator>(const double& lvalue, const ImgClass::Lab& rcolor);

// Arithmetic non-substituting operator
ImgClass::Lab operator+(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor);

ImgClass::Lab operator-(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor);

ImgClass::Lab operator*(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor);
ImgClass::Lab operator*(const ImgClass::Lab& lcolor, const double& rvalue);
ImgClass::Lab operator*(const double& lvalue, const ImgClass::Lab& rcolor);

ImgClass::Lab operator/(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor);
ImgClass::Lab operator/(const ImgClass::Lab& lcolor, const double& rvalue);

// Norm
const ImgClass::Lab abs(const ImgClass::Lab& color);
double norm_squared(const ImgClass::Lab& color);
double norm(const ImgClass::Lab& color);
#ifndef NORM_DOUBLE
#define NORM_DOUBLE
double norm(const double& value);
double norm_squared(const double& value);
#endif

// Product
double inner_prod(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor);
#ifndef INNER_PROD_DOUBLE
#define INNER_PROD_DOUBLE
double inner_prod(const double& lvalue, const double& rvalue);
#endif

// Saturation
ImgClass::Lab saturate(const ImgClass::Lab& value, const double& min, const double& max);

// Quantization
ImgClass::Lab color_quantize(const ImgClass::Lab& value);

// Stream
std::ostream& operator<<(std::ostream& os, const ImgClass::Lab& rcolor);

#endif

