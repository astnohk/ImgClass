#include <cfloat>
#include <cmath>
#include <stdexcept>

#include "RGB.h"


namespace ImgClass {
	RGB::RGB(void)
	{
		R = 0;
		G = 0;
		B = 0;
	}

	RGB::RGB(const double& red, const double& green, const double& blue)
	{
		R = red;
		G = green;
		B = blue;
	}

	RGB::RGB(const RGB& color)
	{
		R = color.R;
		G = color.G;
		B = color.B;
	}


	RGB &
	RGB::set(const double& red, const double& green, const double& blue)
	{
		R = red;
		G = green;
		B = blue;
		return *this;
	}


	RGB &
	RGB::gamma(const double& gamma_val)
	{
		R = pow(R, gamma_val);
		G = pow(G, gamma_val);
		B = pow(B, gamma_val);
		return *this;
	}




	// Operators
	RGB::operator double() const // return intensity
	{
		const double yum_y_red = 0.299;
		const double yum_y_green = 0.587;
		const double yum_y_blue = 0.114;

		return yum_y_red * R + yum_y_green * G + yum_y_blue * B;
	}


	RGB &
	RGB::operator=(const RGB& rcolor)
	{
		R = rcolor.R;
		G = rcolor.G;
		B = rcolor.B;
		return *this;
	}

	RGB &
	RGB::operator=(const double& rvalue)
	{
		R = rvalue;
		G = rvalue;
		B = rvalue;
		return *this;
	}


	RGB &
	RGB::operator+=(const RGB& rcolor)
	{
		R += rcolor.R;
		G += rcolor.G;
		B += rcolor.B;
		return *this;
	}

	RGB &
	RGB::operator+=(const double& rvalue)
	{
		R += rvalue;
		G += rvalue;
		B += rvalue;
		return *this;
	}

	RGB &
	RGB::operator-=(const RGB& rcolor)
	{
		R -= rcolor.R;
		G -= rcolor.G;
		B -= rcolor.B;
		return *this;
	}

	RGB &
	RGB::operator-=(const double& rvalue)
	{
		R -= rvalue;
		G -= rvalue;
		B -= rvalue;
		return *this;
	}

	RGB &
	RGB::operator*=(const RGB& rcolor)
	{
		R *= rcolor.R;
		G *= rcolor.G;
		B *= rcolor.B;
		return *this;
	}

	RGB &
	RGB::operator*=(const double& rvalue)
	{
		R *= rvalue;
		G *= rvalue;
		B *= rvalue;
		return *this;
	}

	RGB &
	RGB::operator/=(const RGB& rcolor)
	{
		R /= rcolor.R;
		G /= rcolor.G;
		B /= rcolor.B;
		return *this;
	}

	RGB &
	RGB::operator/=(const double& rvalue)
	{
		R /= rvalue;
		G /= rvalue;
		B /= rvalue;
		return *this;
	}
}




// ----- Global Operators -----

// Arithmetic

const ImgClass::RGB
operator+(ImgClass::RGB rcolor)
{
	return rcolor;
}

const ImgClass::RGB
operator-(ImgClass::RGB rcolor)
{
	rcolor.R = -rcolor.R;
	rcolor.G = -rcolor.G;
	rcolor.B = -rcolor.B;
	return rcolor;
}


const ImgClass::RGB
operator+(const ImgClass::RGB& lcolor, const ImgClass::RGB& rcolor)
{
	ImgClass::RGB color;
	color.R = lcolor.R + rcolor.R;
	color.G = lcolor.G + rcolor.G;
	color.B = lcolor.B + rcolor.B;
	return color;
}

const ImgClass::RGB
operator-(const ImgClass::RGB& lcolor, const ImgClass::RGB& rcolor)
{
	ImgClass::RGB color;
	color.R = lcolor.R - rcolor.R;
	color.G = lcolor.G - rcolor.G;
	color.B = lcolor.B - rcolor.B;
	return color;
}


const ImgClass::RGB
operator*(const ImgClass::RGB& lcolor, const ImgClass::RGB& rcolor)
{
	ImgClass::RGB color;
	color.R = lcolor.R * rcolor.R;
	color.G = lcolor.G * rcolor.G;
	color.B = lcolor.B * rcolor.B;
	return color;
}

const ImgClass::RGB
operator*(const ImgClass::RGB& lcolor, const double& rvalue)
{
	ImgClass::RGB color;
	color.R = lcolor.R * rvalue;
	color.G = lcolor.G * rvalue;
	color.B = lcolor.B * rvalue;
	return color;
}

const ImgClass::RGB
operator*(const double& lvalue, const ImgClass::RGB& rcolor)
{
	ImgClass::RGB color;
	color.R = lvalue * rcolor.R;
	color.G = lvalue * rcolor.G;
	color.B = lvalue * rcolor.B;
	return color;
}

const ImgClass::RGB
operator/(const ImgClass::RGB& lcolor, const ImgClass::RGB& rcolor)
{
	ImgClass::RGB color;
	color.R = lcolor.R / rcolor.R;
	color.G = lcolor.G / rcolor.G;
	color.B = lcolor.B / rcolor.B;
	return color;
}

const ImgClass::RGB
operator/(const ImgClass::RGB& lcolor, const double& rvalue)
{
	ImgClass::RGB color;
	color.R = lcolor.R / rvalue;
	color.G = lcolor.G / rvalue;
	color.B = lcolor.B / rvalue;
	return color;
}


// Comparator
bool
operator==(const ImgClass::RGB& lcolor, const ImgClass::RGB& rcolor)
{
	if (fabs(lcolor.R - rcolor.R) <= DBL_EPSILON
	    && fabs(lcolor.G - rcolor.G) <= DBL_EPSILON
	    && fabs(lcolor.B - rcolor.B) <= DBL_EPSILON) {
		return true;
	} else {
		return false;
	}
}

bool
operator!=(const ImgClass::RGB& lcolor, const ImgClass::RGB& rcolor)
{
	if (fabs(lcolor.R - rcolor.R) > DBL_EPSILON
	    || fabs(lcolor.G - rcolor.G) > DBL_EPSILON
	    || fabs(lcolor.B - rcolor.B) > DBL_EPSILON) {
		return true;
	} else {
		return false;
	}
}


// Product
double
inner_prod(const ImgClass::RGB& lcolor, const ImgClass::RGB& rcolor)
{
	return lcolor.R * rcolor.R
	    + lcolor.G * rcolor.G
	    + lcolor.B * rcolor.B;
}


// Norm
double
norm_squared(const ImgClass::RGB& color)
{
	return color.R * color.R
	    + color.G * color.G
	    + color.B * color.B;
}

double
norm(const ImgClass::RGB& color)
{
	return sqrt(color.R * color.R
	    + color.G * color.G
	    + color.B * color.B);
}


// Stream
std::ostream &
operator<<(std::ostream& os, const ImgClass::RGB& rcolor)
{
	os << "[R:" << rcolor.R << " G:" << rcolor.G << " B:" << rcolor.B << "]";
	return os;
}

