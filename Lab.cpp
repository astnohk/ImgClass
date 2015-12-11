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



	// ----- Not substitute -----
	Lab
	Lab::operator+(const Lab& rcolor) const
	{
		Lab ret;

		ret.L = L + rcolor.L;
		ret.a = a + rcolor.a;
		ret.b = b + rcolor.b;
		return ret;
	}

	Lab
	Lab::operator+(const double& rvalue) const
	{
		Lab ret;

		ret.L = L + rvalue;
		ret.a = a + rvalue;
		ret.b = b + rvalue;
		return ret;
	}


	Lab
	Lab::operator-(const Lab& rcolor) const
	{
		Lab ret;

		ret.L = L - rcolor.L;
		ret.a = a - rcolor.a;
		ret.b = b - rcolor.b;
		return ret;
	}

	Lab
	Lab::operator-(const double& rvalue) const
	{
		Lab ret;

		ret.L = L - rvalue;
		ret.a = a - rvalue;
		ret.b = b - rvalue;
		return ret;
	}


	Lab
	Lab::operator*(const Lab& rcolor) const
	{
		Lab ret;

		ret.L = L * rcolor.L;
		ret.a = a * rcolor.a;
		ret.b = b * rcolor.b;
		return ret;
	}

	Lab
	Lab::operator*(const double& rvalue) const
	{
		Lab ret;

		ret.L = L * rvalue;
		ret.a = a * rvalue;
		ret.b = b * rvalue;
		return ret;
	}


	Lab
	Lab::operator/(const Lab& rcolor) const
	{
		Lab ret;

		ret.L = L / rcolor.L;
		ret.a = a / rcolor.a;
		ret.b = b / rcolor.b;
		return ret;
	}

	Lab
	Lab::operator/(const double& rvalue) const
	{
		Lab ret;

		ret.L = L / rvalue;
		ret.a = a / rvalue;
		ret.b = b / rvalue;
		return ret;
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


/* Can be overloaded in C++11 or later
const ImgClass::Lab
operator+(ImgClass::Lab lcolor, const ImgClass::Lab& rcolor)
{
	lcolor.L += rcolor.L;
	lcolor.a += rcolor.a;
	lcolor.b += rcolor.b;
	return lcolor;
}


const ImgClass::Lab
operator-(ImgClass::Lab lcolor, const ImgClass::Lab& rcolor)
{
	lcolor.L -= rcolor.L;
	lcolor.a -= rcolor.a;
	lcolor.b -= rcolor.b;
	return lcolor;
}


const ImgClass::Lab
operator*(ImgClass::Lab lcolor, const ImgClass::Lab& rcolor)
{
	lcolor.L *= rcolor.L;
	lcolor.a *= rcolor.a;
	lcolor.b *= rcolor.b;
	return lcolor;
}


const ImgClass::Lab
operator/(ImgClass::Lab lcolor, const ImgClass::Lab& rcolor)
{
	lcolor.L /= rcolor.L;
	lcolor.a /= rcolor.a;
	lcolor.b /= rcolor.b;
	return lcolor;
}
*/


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


// Norm
double
norm_squared(const ImgClass::Lab& color)
{
	return color.L * color.L
	    + color.a * color.a
	    + color.b * color.b;
}

