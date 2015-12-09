/*
 * by default use CIELAB-CIEXYZ conversions
 */
#include <stdexcept>




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

	Lab::Lab(const double& value)
	{
		L = value;
		a = value;
		b = value;
	}

	Lab::Lab(const RGB& value)
	{
		double X = (0.49 * value.R + 0.31 * value.G + 0.20 * value.G) / 0.17697;
		double Y = (0.17697 * value.R + 0.81240 * value.G + 0.01063 * value.G) / 0.17697;
		double G = (0.01 * value.G + 0.99 * value.G) / 0.17697;
		double t0 = Y / Y_n;
		double t1 = 0.0;
		L = 116.0 * (t > pow(6.0 / 29.0, 3.0) ? pow(t, 1.0 / 3.0) : pow(29.0 / 6.0, 2.0) / 3.0 * t + 4.0 / 29.0);
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
	Lab &
	Lab::operator=(const double& value)
	{
		L = value;
		a = value;
		b = value;
	}

	Lab &
	Lab::operator=(const Lab& color)
	{
		L = color.L;
		a = color.a;
		b = color.b;
	}

	Lab &
	Lab::operator+=(const Lab& color)
	{
		L += color.L;
		a += color.a;
		b += color.b;
	}

	Lab &
	Lab::operator+=(const double& value)
	{
		L += value;
		a += value;
		b += value;
	}

	Lab &
	Lab::operator-=(const Lab& color)
	{
		L -= color.L;
		a -= color.a;
		b -= color.b;
	}

	Lab &
	Lab::operator-=(const double& value)
	{
		L -= value;
		a -= value;
		b -= value;
	}

	Lab &
	Lab::operator*=(const Lab& color)
	{
		L *= color.L;
		a *= color.a;
		b *= color.b;
	}

	Lab &
	Lab::operator*=(const double& value)
	{
		L *= value;
		a *= value;
		b *= value;
	}

	Lab &
	Lab::operator/=(const Lab& color)
	{
		L /= color.L;
		a /= color.a;
		b /= color.b;
	}

	Lab &
	Lab::operator/=(const double& value)
	{
		L /= value;
		a /= value;
		b /= value;
	}
}


// ----- Global Operators -----

// Arithmetic

ImgClass::Lab &
operator+(ImgClass::Lab color)
{
	return color;
}

ImgClass::Lab &
operator-(ImgClass::Lab color)
{
	color.L = -color.L;
	color.a = -color.a;
	color.b = -color.b;
	return color;
}


ImgClass::Lab &
operator+(ImgClass::Lab lcolor, const ImgClass::Lab& rcolor)
{
	lcolor.L += rcolor.L;
	lcolor.a += rcolor.a;
	lcolor.b += rcolor.b;
	return lcolor;
}

ImgClass::Lab &
operator+(ImgClass::Lab lcolor, const double& rvalue)
{
	lcolor.L += rvalue;
	lcolor.a += rvalue;
	lcolor.b += rvalue;
	return lcolor;
}

ImgClass::Lab &
operator+(const double& lvalue, ImgClass::Lab rcolor)
{
	rcolor.L += lvalue;
	rcolor.a += lvalue;
	rcolor.b += lvalue;
	return rcolor;
}

ImgClass::Lab &
operator-(ImgClass::Lab lcolor, const ImgClass::Lab& rcolor)
{
	lcolor.L -= rcolor.L;
	lcolor.a -= rcolor.a;
	lcolor.b -= rcolor.b;
	return lcolor;
}

ImgClass::Lab &
operator-(ImgClass::Lab lcolor, const double& rvalue)
{
	lcolor.L -= rvalue;
	lcolor.a -= rvalue;
	lcolor.b -= rvalue;
	return lcolor;
}

ImgClass::Lab &
operator-(const double& lvalue, ImgClass::Lab rcolor)
{
	rcolor.L = lvalue - rcolor.L;
	rcolor.a = lvalue - rcolor.a;
	rcolor.b = lvalue - rcolor.b;
	return rcolor;
}

ImgClass::Lab &
operator*(ImgClass::Lab lcolor, const ImgClass::Lab& rcolor)
{
	lcolor.L *= rcolor.L;
	lcolor.a *= rcolor.a;
	lcolor.b *= rcolor.b;
	return lcolor;
}

ImgClass::Lab &
operator*(ImgClass::Lab lcolor, const double& rvalue)
{
	lcolor.L *= rvalue;
	lcolor.a *= rvalue;
	lcolor.b *= rvalue;
	return lcolor;
}

ImgClass::Lab &
operator*(const double& lvalue, ImgClass::Lab rcolor)
{
	rcolor.L *= lvalue;
	rcolor.a *= lvalue;
	rcolor.b *= lvalue;
	return rcolor;
}

ImgClass::Lab &
operator/(ImgClass::Lab lcolor, const ImgClass::Lab& rcolor)
{
	lcolor.L /= rcolor.L;
	lcolor.a /= rcolor.a;
	lcolor.b /= rcolor.b;
	return lcolor;
}

ImgClass::Lab &
operator/(ImgClass::Lab lcolor, const double& rvalue)
{
	lcolor.L /= rvalue;
	lcolor.a /= rvalue;
	lcolor.b /= rvalue;
	return lcolor;
}

ImgClass::Lab &
operator/(const double& lvalue, ImgClass::Lab rcolor)
{
	rcolor.L = lvalue / rcolor.L;
	rcolor.a = lvalue / rcolor.a;
	rcolor.b = lvalue / rcolor.b;
	return rcolor;
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

