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

	RGB::RGB(const double& red, const double& green, const double& blue, const double& gamma_val)
	{
		if (gamma_val > 0.0) {
			R = pow(red, gamma_val);
			G = pow(green, gamma_val);
			B = pow(blue, gamma_val);
		} else {
			R = red;
			G = green;
			B = blue;
		}
	}

	RGB::RGB(const RGB& color)
	{
		R = color.R;
		G = color.G;
		B = color.B;
	}


	RGB &
	RGB::set(const double& red, const double& green, const double& blue, const double& gamma_val)
	{
		if (gamma_val > 0.0) {
			R = pow(red, gamma_val);
			G = pow(green, gamma_val);
			B = pow(blue, gamma_val);
		} else {
			R = red;
			G = green;
			B = blue;
		}
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


	const RGB
	RGB::operator+(const RGB& rcolor) const
	{
		RGB color;
		color.R += rcolor.R;
		color.G += rcolor.G;
		color.B += rcolor.B;
		return color;
	}

	const RGB
	RGB::operator+(const double& rvalue) const
	{
		RGB color;
		color.R += rvalue;
		color.G += rvalue;
		color.B += rvalue;
		return color;
	}


	const RGB
	RGB::operator-(const RGB& rcolor) const
	{
		RGB color;
		color.R -= rcolor.R;
		color.G -= rcolor.G;
		color.B -= rcolor.B;
		return color;
	}

	const RGB
	RGB::operator-(const double& rvalue) const
	{
		RGB color;
		color.R -= rvalue;
		color.G -= rvalue;
		color.B -= rvalue;
		return color;
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

const ImgClass::RGB
operator/(const double& lvalue, const ImgClass::RGB& rcolor)
{
	ImgClass::RGB color;
	color.R = lvalue / rcolor.R;
	color.G = lvalue / rcolor.G;
	color.B = lvalue / rcolor.B;
	return color;
}


// Comparator
bool
operator==(const ImgClass::RGB& lcolor, const ImgClass::RGB& rcolor)
{
	if (lcolor.R == rcolor.R
	    && lcolor.G == rcolor.G
	    && lcolor.B == rcolor.B) {
		return true;
	} else {
		return false;
	}
}

bool
operator!=(const ImgClass::RGB& lcolor, const ImgClass::RGB& rcolor)
{
	if (lcolor.R != rcolor.R
	    || lcolor.G != rcolor.G
	    || lcolor.B != rcolor.B) {
		return true;
	} else {
		return false;
	}
}

