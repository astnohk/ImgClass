#include <cfloat>
#include <cmath>
#include <stdexcept>
#include <iostream>

#include "Color.h"


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

	RGB::RGB(const RGB& rgb)
	{
		R = rgb.R;
		G = rgb.G;
		B = rgb.B;
	}

	RGB::RGB(const HSV& hsv)
	{
		double C = hsv.V * hsv.S;
		double X = C * (1.0 - fabs(fmod(hsv.H * 6.0, 2.0) - 1.0));

		R = hsv.V - C;
		G = hsv.V - C;
		B = hsv.V - C;
		switch (int(floor(hsv.H * 6.0))) {
			case 0: *this += RGB(C, X, 0); break;
			case 1: *this += RGB(X, C, 0); break;
			case 2: *this += RGB(0, C, X); break;
			case 3: *this += RGB(0, X, C); break;
			case 4: *this += RGB(X, 0, C); break;
			case 5: *this += RGB(C, 0, X);
		}
	}

	RGB::RGB(const Lab& lab)
	{
		auto f_inv = [](double t_inv) -> double {
			if (t_inv > 6.0 / 29.0) {
				return pow(t_inv, 3.0);
			} else {
				return (116.0 * t_inv - 16.0) * 27.0 / 24389.0;
			}
		};
		double f_y = (lab.L + 16.0) / 116.0;
		double f_x = f_y + lab.a / 500.0;
		double f_z = f_y - lab.b / 200.0;
		double Y;
		if (lab.L > 216.0 / 27.0) {
			Y = Y_n * pow(f_y, 3.0);
		} else {
			Y = Y_n * lab.L * 27.0 / 24389.0;
		}
		double X = X_n * f_inv(f_x);
		double Z = Z_n * f_inv(f_z);
		R = 3.2404542 * X - 1.5371385 * Y - 0.4985314 * Z;
		G = -0.969266 * X + 1.8760108 * Y + 0.0415560 * Z;
		B = 0.0556434 * X - 0.2040259 * Y + 1.0572252 * Z;
		this->gamma(1.0 / 2.2); // Convert linear sRGB to sRGB
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




// Saturation
ImgClass::RGB
saturate(const ImgClass::RGB& value, const double& min, const double& max)
{
	ImgClass::RGB ret(value);
	auto lambda = [&min, &max](const double& val) -> double {
		if (val < min) {
			return min;
		} else if (val > max) {
			return max;
		} else {
			return val;
		}
	};
	ret.R = lambda(ret.R);
	ret.G = lambda(ret.G);
	ret.B = lambda(ret.B);
	return ret;
}


// Quantization
ImgClass::RGB
color_quantize(const ImgClass::RGB &value, const double &max)
{
	ImgClass::RGB ret;
	ret.R = round(max * value.R);
	ret.G = round(max * value.G);
	ret.B = round(max * value.B);
	return ret;
}


// Stream
std::ostream &
operator<<(std::ostream& os, const ImgClass::RGB& rcolor)
{
	os << "[R:" << rcolor.R << " G:" << rcolor.G << " B:" << rcolor.B << "]";
	return os;
}

