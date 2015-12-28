/*
 * by default use CIELAB-CIEXYZ conversions
 */
#include <cfloat>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "Lab.h"




// d65
static const double X_n = 0.95045;
static const double Y_n = 1.0;
static const double Z_n = 1.088917;

// d50
//static const double X_n = 0.9642;
//static const double Y_n = 1.0;
//static const double Z_n = 0.8249;

// flat
//static const double X_n = 1.0;
//static const double Y_n = 1.0;
//static const double Z_n = 1.0;




namespace ImgClass {
	Lab::Lab(void)
	{
		L = 0;
		a = 0;
		b = 0;
	}

	Lab::Lab(const double& _L, const double& _a, const double& _b)
	{
		L = _L;
		a = _a;
		b = _b;
	}

	Lab::Lab(const Lab& color)
	{
		L = color.L;
		a = color.a;
		b = color.b;
	}

	Lab::Lab(const RGB& color)
	{
		RGB linear_sRGB(color);
		linear_sRGB.gamma(2.2); // Convert sRGB to linear sRGB
		double X = 0.4124564 * linear_sRGB.R + 0.3575761 * linear_sRGB.G + 0.1804375 * linear_sRGB.B;
		double Y = 0.2126729 * linear_sRGB.R + 0.7151522 * linear_sRGB.G + 0.0721750 * linear_sRGB.B;
		double Z = 0.0193339 * linear_sRGB.R + 0.1191920 * linear_sRGB.G + 0.9503041 * linear_sRGB.B;
		L = 116.0 * this->f(Y / Y_n) - 16.0;
		a = 500.0 * (this->f(X / X_n) - this->f(Y / Y_n));
		b = 200.0 * (this->f(Y / Y_n) - this->f(Z / Z_n));
	}

	double
	Lab::f(const double t)
	{
		if (t > pow(6.0 / 29.0, 3.0)) {
			return pow(t, 1.0 / 3.0);
		} else {
			return pow(29.0 / 6.0, 2.0) / 3.0 * t + 4.0 / 29.0;
		}
	}




	Lab &
	Lab::set(const RGB& color)
	{
		Lab tmp(color); // Convert with constructor
		L = tmp.L;
		a = tmp.a;
		b = tmp.b;
		return *this;
	}

	// Operators
	Lab::operator double() const
	{
		return L;
	}


	Lab &
	Lab::operator=(const Lab& color)
	{
		L = color.L;
		a = color.a;
		b = color.b;
		return *this;
	}

	Lab &
	Lab::operator=(const double& value)
	{
		L = value;
		a = value;
		b = value;
		return *this;
	}

	Lab &
	Lab::operator+=(const Lab& color)
	{
		L += color.L;
		a += color.a;
		b += color.b;
		return *this;
	}

	Lab &
	Lab::operator-=(const Lab& color)
	{
		L -= color.L;
		a -= color.a;
		b -= color.b;
		return *this;
	}

	Lab &
	Lab::operator*=(const Lab& color)
	{
		L *= color.L;
		a *= color.a;
		b *= color.b;
		return *this;
	}

	Lab &
	Lab::operator*=(const double& value)
	{
		L *= value;
		a *= value;
		b *= value;
		return *this;
	}

	Lab &
	Lab::operator/=(const Lab& color)
	{
		L /= color.L;
		a /= color.a;
		b /= color.b;
		return *this;
	}

	Lab &
	Lab::operator/=(const double& value)
	{
		L /= value;
		a /= value;
		b /= value;
		return *this;
	}
}


// ----- Global Operators -----

// Arithmetic

const ImgClass::Lab
operator+(ImgClass::Lab color)
{
	return color;
}

const ImgClass::Lab
operator-(ImgClass::Lab color)
{
	color.L = -color.L;
	color.a = -color.a;
	color.b = -color.b;
	return color;
}




// Comparator
bool
operator==(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor)
{
	if (fabs(lcolor.L - rcolor.L) <= DBL_EPSILON
	    && fabs(lcolor.a - rcolor.a) <= DBL_EPSILON
	    && fabs(lcolor.b - rcolor.b) <= DBL_EPSILON) {
		return true;
	} else {
		return false;
	}
}

bool
operator!=(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor)
{
	if (fabs(lcolor.L - rcolor.L) > DBL_EPSILON
	    || fabs(lcolor.a - rcolor.a) > DBL_EPSILON
	    || fabs(lcolor.b - rcolor.b) > DBL_EPSILON) {
		return true;
	} else {
		return false;
	}
}


bool
operator<(const ImgClass::Lab& lcolor, const double& rvalue)
{
	if (lcolor.L * lcolor.L + lcolor.a * lcolor.a + lcolor.b * lcolor.b
	    < rvalue * rvalue) {
		return true;
	} else {
		return false;
	}
}

bool
operator<(const double& lvalue, const ImgClass::Lab& rcolor)
{
	if (lvalue * lvalue
	    <rcolor.L * rcolor.L + rcolor.a * rcolor.a + rcolor.b * rcolor.b) {
		return true;
	} else {
		return false;
	}
}

bool
operator>(const ImgClass::Lab& lcolor, const double& rvalue)
{
	if (lcolor.L * lcolor.L + lcolor.a * lcolor.a + lcolor.b * lcolor.b
	    > rvalue * rvalue) {
		return true;
	} else {
		return false;
	}
}

bool
operator>(const double& lvalue, const ImgClass::Lab& rcolor)
{
	if (lvalue * lvalue
	    > rcolor.L * rcolor.L + rcolor.a * rcolor.a + rcolor.b * rcolor.b) {
		return true;
	} else {
		return false;
	}
}



// Arithmetic operators
ImgClass::Lab operator+(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor)
{
	ImgClass::Lab color;

	color.L = lcolor.L + rcolor.L;
	color.a = lcolor.a + rcolor.a;
	color.b = lcolor.b + rcolor.b;
	return color;
}

ImgClass::Lab operator+(const ImgClass::Lab& lcolor, const double& rvalue)
{
	ImgClass::Lab color;

	color.L = lcolor.L + rvalue;
	color.a = lcolor.a + rvalue;
	color.b = lcolor.b + rvalue;
	return color;
}

ImgClass::Lab operator+(const double& lvalue, const ImgClass::Lab& rcolor)
{
	ImgClass::Lab color;

	color.L = lvalue + rcolor.L;
	color.a = lvalue + rcolor.a;
	color.b = lvalue + rcolor.b;
	return color;
}


ImgClass::Lab operator-(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor)
{
	ImgClass::Lab color;

	color.L = lcolor.L - rcolor.L;
	color.a = lcolor.a - rcolor.a;
	color.b = lcolor.b - rcolor.b;
	return color;
}

ImgClass::Lab operator-(const ImgClass::Lab& lcolor, const double& rvalue)
{
	ImgClass::Lab color;

	color.L = lcolor.L - rvalue;
	color.a = lcolor.a - rvalue;
	color.b = lcolor.b - rvalue;
	return color;
}

ImgClass::Lab operator-(const double& lvalue, const ImgClass::Lab& rcolor)
{
	ImgClass::Lab color;

	color.L = lvalue - rcolor.L;
	color.a = lvalue - rcolor.a;
	color.b = lvalue - rcolor.b;
	return color;
}


ImgClass::Lab operator*(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor)
{
	ImgClass::Lab color;

	color.L = lcolor.L * rcolor.L;
	color.a = lcolor.a * rcolor.a;
	color.b = lcolor.b * rcolor.b;
	return color;
}

ImgClass::Lab operator*(const ImgClass::Lab& lcolor, const double& rvalue)
{
	ImgClass::Lab color;

	color.L = lcolor.L * rvalue;
	color.a = lcolor.a * rvalue;
	color.b = lcolor.b * rvalue;
	return color;
}

ImgClass::Lab operator*(const double& lvalue, const ImgClass::Lab& rcolor)
{
	ImgClass::Lab color;

	color.L = lvalue * rcolor.L;
	color.a = lvalue * rcolor.a;
	color.b = lvalue * rcolor.b;
	return color;
}


ImgClass::Lab operator/(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor)
{
	ImgClass::Lab color;

	color.L = lcolor.L / rcolor.L;
	color.a = lcolor.a / rcolor.a;
	color.b = lcolor.b / rcolor.b;
	return color;
}

ImgClass::Lab operator/(const ImgClass::Lab& lcolor, const double& rvalue)
{
	ImgClass::Lab color;

	color.L = lcolor.L / rvalue;
	color.a = lcolor.a / rvalue;
	color.b = lcolor.b / rvalue;
	return color;
}




const ImgClass::Lab
abs(const ImgClass::Lab& color)
{
	ImgClass::Lab ret;

	ret.L = color.L;
	ret.a = color.a;
	ret.b = color.b;
	return ret;
}


// Product
double
inner_prod(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor)
{
	return lcolor.L * rcolor.L
	    + lcolor.a * rcolor.a
	    + lcolor.b * rcolor.b;
}


// Norm
double
norm_squared(const ImgClass::Lab& color)
{
	return color.L * color.L
	    + color.a * color.a
	    + color.b * color.b;
}

double
norm(const ImgClass::Lab& color)
{
	return sqrt(color.L * color.L
	    + color.a * color.a
	    + color.b * color.b);
}


// Stream
std::ostream &
operator<<(std::ostream& os, const ImgClass::Lab& rcolor)
{
	os << "[L*:" << rcolor.L << " a*:" << rcolor.a << " b*:" << rcolor.b << "]";
	return os;
}

