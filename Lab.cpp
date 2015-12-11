/*
 * by default use CIELAB-CIEXYZ conversions
 */

#include <cmath>
#include <stdexcept>

#include "Lab.h"




static const double X_n = 95.047;
static const double Y_n = 100.000;
static const double Z_n = 108.883;




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

	Lab::Lab(const RGB<double>& color)
	{
		double X = (0.49 * color.R + 0.31 * color.G + 0.20 * color.G) / 0.17697;
		double Y = (0.17697 * color.R + 0.81240 * color.G + 0.01063 * color.G) / 0.17697;
		double Z = (0.01 * color.G + 0.99 * color.G) / 0.17697;
		double t0 = Y / Y_n;
		double t1 = 0.0;
		L = 116.0 * (t0 > pow(6.0 / 29.0, 3.0) ? pow(t0, 1.0 / 3.0) : pow(29.0 / 6.0, 2.0) / 3.0 * t0 + 4.0 / 29.0);
		t0 = X / X_n;
		t1 = Y / Y_n;
		a = 500.0 * (
		    (t0 > pow(6.0 / 29.0, 3.0) ? pow(t0, 1.0 / 3.0) : pow(29.0 / 6.0, 2.0) / 3.0 * t0 + 4.0 / 29.0)
		    - (t1 > pow(6.0 / 29.0, 3.0) ? pow(t1, 1.0 / 3.0) : pow(29.0 / 6.0, 2.0) / 3.0 * t1 + 4.0 / 29.0));
		t0 = Y / Y_n;
		t1 = Z / Z_n;
		b = 200.0 * (
		    (t0 > pow(6.0 / 29.0, 3.0) ? pow(t0, 1.0 / 3.0) : pow(29.0 / 6.0, 2.0) / 3.0 * t0 + 4.0 / 29.0)
		    - (t1 > pow(6.0 / 29.0, 3.0) ? pow(t1, 1.0 / 3.0) : pow(29.0 / 6.0, 2.0) / 3.0 * t1 + 4.0 / 29.0));
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
	Lab::operator=(const RGB<double>& value)
	{
		double X = (0.49 * value.R + 0.31 * value.G + 0.20 * value.G) / 0.17697;
		double Y = (0.17697 * value.R + 0.81240 * value.G + 0.01063 * value.G) / 0.17697;
		double Z = (0.01 * value.G + 0.99 * value.G) / 0.17697;

		double t0 = Y / Y_n;
		double t1 = 0.0;
		L = 116.0 * (t0 > pow(6.0 / 29.0, 3.0) ? pow(t0, 1.0 / 3.0) : pow(29.0 / 6.0, 2.0) / 3.0 * t0 + 4.0 / 29.0);

		t0 = X / X_n;
		t1 = Y / Y_n;
		a = 500.0 * (
		    (t0 > pow(6.0 / 29.0, 3.0) ? pow(t0, 1.0 / 3.0) : pow(29.0 / 6.0, 2.0) / 3.0 * t0 + 4.0 / 29.0)
		    - (t1 > pow(6.0 / 29.0, 3.0) ? pow(t1, 1.0 / 3.0) : pow(29.0 / 6.0, 2.0) / 3.0 * t1 + 4.0 / 29.0));

		t0 = Y / Y_n;
		t1 = Z / Z_n;
		b = 200.0 * (
		    (t0 > pow(6.0 / 29.0, 3.0) ? pow(t0, 1.0 / 3.0) : pow(29.0 / 6.0, 2.0) / 3.0 * t0 + 4.0 / 29.0)
		    - (t1 > pow(6.0 / 29.0, 3.0) ? pow(t1, 1.0 / 3.0) : pow(29.0 / 6.0, 2.0) / 3.0 * t1 + 4.0 / 29.0));

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
	if (lcolor.L == rcolor.L
	    && lcolor.a == rcolor.a
	    && lcolor.b == rcolor.b) {
		return true;
	} else {
		return false;
	}
}

bool
operator==(const ImgClass::Lab& lcolor, const double& rvalue)
{
	if (lcolor.L == rvalue
	    && lcolor.a == rvalue
	    && lcolor.b == rvalue) {
		return true;
	} else {
		return false;
	}
}

bool
operator==(const double& lvalue, const ImgClass::Lab& rcolor)
{
	if (lvalue == rcolor.L
	    && lvalue == rcolor.a
	    && lvalue == rcolor.b) {
		return true;
	} else {
		return false;
	}
}


bool
operator!=(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor)
{
	if (lcolor.L != rcolor.L
	    || lcolor.a != rcolor.a
	    || lcolor.b != rcolor.b) {
		return true;
	} else {
		return false;
	}
}

bool
operator!=(const ImgClass::Lab& lcolor, const double& rvalue)
{
	if (lcolor.L != rvalue
	    || lcolor.a != rvalue
	    || lcolor.b != rvalue) {
		return true;
	} else {
		return false;
	}
}

bool
operator!=(const double& lvalue, const ImgClass::Lab& rcolor)
{
	if (lvalue != rcolor.L
	    || lvalue != rcolor.a
	    || lvalue != rcolor.b) {
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

double
fabs(const ImgClass::Lab& color)
{
	return sqrt(color.L * color.L + color.a * color.a + color.b * color.b);
}

// Norm
double
norm_squared(const ImgClass::Lab& color)
{
	return color.L * color.L
	    + color.a * color.a
	    + color.b * color.b;
}

