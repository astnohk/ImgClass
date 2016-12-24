#ifndef LIB_ImgClass_HSV
#define LIB_ImgClass_HSV

#include <ostream>
#include "RGB.h"

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

		// Converters
		RGB get_RGB(void);
	};
}


const ImgClass::HSV operator+(ImgClass::HSV rval);
const ImgClass::HSV operator-(ImgClass::HSV rval);

const ImgClass::HSV operator+(const ImgClass::HSV& lval, const ImgClass::HSV& rval);
const ImgClass::HSV operator-(const ImgClass::HSV& lval, const ImgClass::HSV& rval);

bool operator==(const ImgClass::HSV& lval, const ImgClass::HSV& rval);
bool operator!=(const ImgClass::HSV& lval, const ImgClass::HSV& rval);

std::ostream& operator<<(std::ostream& os, const ImgClass::HSV& rval);

#endif

