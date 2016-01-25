#ifndef LIB_ImgClass_RGB
#define LIB_ImgClass_RGB

#include <ostream>

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
			RGB(const RGB& color); // Copy constructor

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

// Stream
std::ostream& operator<<(std::ostream& os, const ImgClass::RGB& rcolor);

#endif

