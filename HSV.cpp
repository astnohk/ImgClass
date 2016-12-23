#include <cfloat>
#include <cmath>
#include <stdexcept>

#include "HSV.h"


namespace ImgClass {
	HSV::HSV(void)
	{
		H = 0;
		S = 0;
		V = 0;
	}

	HSV::HSV(const double& hue, const double& saturation, const double& value)
	{
		H = hue;
		S = saturation;
		V = value;
	}

	HSV::HSV(const HSV& hsv)
	{
		H = hsv.H;
		S = hsv.S;
		V = hsv.V;
	}

	HSV::HSV(const RGB& rgb)
	{
		double max = rgb.R;
		double min = rgb.R;
		if (max < rgb.G) {
			max = rgb.G;
		}
		if (max < rgb.B) {
			max = rgb.B;
		}
		if (min > rgb.G) {
			min = rgb.G;
		}
		if (min > rgb.B) {
			min = rgb.B;
		}
		H = fmod((rgb.G - rgb.R) / (max - min) + 1.0 / 6.0, 1.0);
		if (H < 0) {
			H += 1.0;
		}
		S = max - min;
		V = max;
	}

	HSV &
	HSV::set(const double& hue, const double& saturation, const double& value)
	{
		H = hue;
		S = saturation;
		V = value;
		return *this;
	}

	HSV &
	HSV::set(const RGB& rgb)
	{
		double max = rgb.R;
		double min = rgb.R;
		if (max < rgb.G) {
			max = rgb.G;
		}
		if (max < rgb.B) {
			max = rgb.B;
		}
		if (min > rgb.G) {
			min = rgb.G;
		}
		if (min > rgb.B) {
			min = rgb.B;
		}
		H = fmod((rgb.G - rgb.R) / (max - min) + 1.0 / 6.0, 1.0);
		if (H < 0) {
			H += 1.0;
		}
		S = max - min;
		V = max;
		return *this;
	}


	HSV &
	HSV::operator=(const RGB& rgb)
	{
		double max = rgb.R;
		double min = rgb.R;
		if (max < rgb.G) {
			max = rgb.G;
		}
		if (max < rgb.B) {
			max = rgb.B;
		}
		if (min > rgb.G) {
			min = rgb.G;
		}
		if (min > rgb.B) {
			min = rgb.B;
		}
		H = fmod((rgb.G - rgb.R) / (max - min) + 1.0 / 6.0, 1.0);
		if (H < 0) {
			H += 1.0;
		}
		S = max - min;
		V = max;
		return *this;
	}

	HSV &
	HSV::operator=(const HSV& rval)
	{
		H = rval.H;
		S = rval.S;
		V = rval.V;
		return *this;
	}

	HSV &
	HSV::operator+=(const HSV& rval)
	{
		H += rval.H;
		S += rval.S;
		V += rval.V;
		return *this;
	}

	HSV &
	HSV::operator-=(const HSV& rval)
	{
		H -= rval.H;
		S -= rval.S;
		V -= rval.V;
		return *this;
	}
}

const ImgClass::HSV
operator+(ImgClass::HSV rval)
{
	return rval;
}

const ImgClass::HSV
operator-(ImgClass::HSV rval)
{
	rval.H = -rval.H;
	rval.S = -rval.S;
	rval.V = -rval.V;
	return rval;
}

const ImgClass::HSV
operator+(const ImgClass::HSV& lval, const ImgClass::HSV& rval)
{
	ImgClass::HSV color;
	color.H = lval.H + rval.H;
	color.S = lval.S + rval.S;
	color.V = lval.V + rval.V;
	return color;
}

const ImgClass::HSV
operator-(const ImgClass::HSV& lval, const ImgClass::HSV& rval)
{
	ImgClass::HSV color;
	color.H = lval.H - rval.H;
	color.S = lval.S - rval.S;
	color.V = lval.V - rval.V;
	return color;
}


bool
operator==(const ImgClass::HSV& lval, const ImgClass::HSV& rval)
{
	if (fabs(lval.H - rval.H) <= DBL_EPSILON
	    && fabs(lval.S - rval.S) <= DBL_EPSILON
	    && fabs(lval.V - rval.V) <= DBL_EPSILON) {
		return true;
	} else {
		return false;
	}
}

bool
operator!=(const ImgClass::HSV& lval, const ImgClass::HSV& rval)
{
	if (fabs(lval.H - rval.H) > DBL_EPSILON
	    || fabs(lval.S - rval.S) > DBL_EPSILON
	    || fabs(lval.V - rval.V) > DBL_EPSILON) {
		return true;
	} else {
		return false;
	}
}


// Stream
std::ostream &
operator<<(std::ostream& os, const ImgClass::HSV& rval)
{
	os << "[H:" << rval.H << " S:" << rval.S << " V:" << rval.V << ":";
	return os;
}

