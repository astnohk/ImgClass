#include <stdexcept>


namespace ImgClass {
	template <class T>
	RGB<T>::RGB(void)
	{
		R = 0;
		G = 0;
		B = 0;
	}

	template <class T>
	RGB<T>::RGB(const T& red, const T& green, const T& blue)
	{
		R = red;
		G = green;
		B = blue;
	}

	template <class T>
	RGB<T>::RGB(const RGB<T>& color)
	{
		R = color.R;
		G = color.G;
		B = color.B;
	}

	template <class T>
	RGB<T>::RGB(const T& value)
	{
		R = value;
		G = value;
		B = value;
	}




	// Operators
	template <class T>
	template <class ConvertType>
	RGB<T>::operator RGB<ConvertType>() const
	{
		RGB<ConvertType> color;

		color.R = ConvertType(color.R);
		color.G = ConvertType(color.G);
		color.B = ConvertType(color.B);

		return color;
	}

	template <class T>
	template <class ConvertType>
	RGB<T>::operator ConvertType() const
	{
		const double yuv_y_red = 0.299;
		const double yuv_y_green = 0.587;
		const double yuv_y_blue = 0.114;

		return ConvertType(yuv_y_red * R
		    + yuv_y_green * G
		    + yuv_y_blue * B);
	}


	template <class T>
	template <class RT>
	RGB<T> &
	RGB<T>::operator=(const RGB<RT>& rcolor)
	{
		R = T(rcolor.R);
		G = T(rcolor.G);
		B = T(rcolor.B);
		return *this;
	}

	template <class T>
	RGB<T> &
	RGB<T>::operator=(const T& rvalue)
	{
		R = rvalue;
		G = rvalue;
		B = rvalue;
		return *this;
	}


	template <class T>
	RGB<T> &
	RGB<T>::operator+=(const RGB<T>& rcolor)
	{
		R += rcolor.R;
		G += rcolor.G;
		B += rcolor.B;
		return *this;
	}

	template <class T>
	RGB<T> &
	RGB<T>::operator+=(const T& rvalue)
	{
		R += rvalue;
		G += rvalue;
		B += rvalue;
		return *this;
	}

	template <class T>
	RGB<T> &
	RGB<T>::operator-=(const RGB<T>& rcolor)
	{
		R -= rcolor.R;
		G -= rcolor.G;
		B -= rcolor.B;
		return *this;
	}

	template <class T>
	RGB<T> &
	RGB<T>::operator-=(const T& rvalue)
	{
		R -= rvalue;
		G -= rvalue;
		B -= rvalue;
		return *this;
	}

	template <class T>
	RGB<T> &
	RGB<T>::operator*=(const RGB<T>& rcolor)
	{
		R *= rcolor.R;
		G *= rcolor.G;
		B *= rcolor.B;
		return *this;
	}

	template <class T>
	RGB<T> &
	RGB<T>::operator*=(const T& rvalue)
	{
		R *= rvalue;
		G *= rvalue;
		B *= rvalue;
		return *this;
	}

	template <class T>
	RGB<T> &
	RGB<T>::operator/=(const RGB<T>& rcolor)
	{
		R /= rcolor.R;
		G /= rcolor.G;
		B /= rcolor.B;
		return *this;
	}

	template <class T>
	RGB<T> &
	RGB<T>::operator/=(const T& rvalue)
	{
		R /= rvalue;
		G /= rvalue;
		B /= rvalue;
		return *this;
	}
}


// ----- Global Operators -----

// Arithmetic

template <class Type>
const ImgClass::RGB<Type>
operator+(ImgClass::RGB<Type> color)
{
	return color;
}

template <class Type>
const ImgClass::RGB<Type>
operator-(ImgClass::RGB<Type> color)
{
	color.R = -color.R;
	color.G = -color.G;
	color.B = -color.B;
	return color;
}


template <class Type>
const ImgClass::RGB<Type>
operator+(ImgClass::RGB<Type> lcolor, const ImgClass::RGB<Type>& rcolor)
{
	lcolor.R += rcolor.R;
	lcolor.G += rcolor.G;
	lcolor.B += rcolor.B;
	return lcolor;
}

template <class Type>
const ImgClass::RGB<Type>
operator+(ImgClass::RGB<Type> lcolor, const Type& rvalue)
{
	lcolor.R += rvalue;
	lcolor.G += rvalue;
	lcolor.B += rvalue;
	return lcolor;
}

template <class Type>
const ImgClass::RGB<Type>
operator+(const Type& lvalue, ImgClass::RGB<Type> rcolor)
{
	rcolor.R += lvalue;
	rcolor.G += lvalue;
	rcolor.B += lvalue;
	return rcolor;
}

template <class Type>
const ImgClass::RGB<Type>
operator-(ImgClass::RGB<Type> lcolor, const ImgClass::RGB<Type>& rcolor)
{
	lcolor.R -= rcolor.R;
	lcolor.G -= rcolor.G;
	lcolor.B -= rcolor.B;
	return lcolor;
}

template <class Type>
const ImgClass::RGB<Type>
operator-(ImgClass::RGB<Type> lcolor, const Type& rvalue)
{
	lcolor.R -= rvalue;
	lcolor.G -= rvalue;
	lcolor.B -= rvalue;
	return lcolor;
}

template <class Type>
const ImgClass::RGB<Type>
operator-(const Type& lvalue, ImgClass::RGB<Type> rcolor)
{
	rcolor.R = lvalue - rcolor.R;
	rcolor.G = lvalue - rcolor.G;
	rcolor.B = lvalue - rcolor.B;
	return rcolor;
}

template <class Type>
const ImgClass::RGB<Type>
operator*(ImgClass::RGB<Type> lcolor, const ImgClass::RGB<Type>& rcolor)
{
	lcolor.R *= rcolor.R;
	lcolor.G *= rcolor.G;
	lcolor.B *= rcolor.B;
	return lcolor;
}

template <class Type>
const ImgClass::RGB<Type>
operator*(ImgClass::RGB<Type> lcolor, const Type& rvalue)
{
	lcolor.R *= rvalue;
	lcolor.G *= rvalue;
	lcolor.B *= rvalue;
	return lcolor;
}

template <class Type>
const ImgClass::RGB<Type>
operator*(const Type& lvalue, ImgClass::RGB<Type> rcolor)
{
	rcolor.R *= lvalue;
	rcolor.G *= lvalue;
	rcolor.B *= lvalue;
	return rcolor;
}

template <class Type>
const ImgClass::RGB<Type>
operator/(ImgClass::RGB<Type> lcolor, const ImgClass::RGB<Type>& rcolor)
{
	lcolor.R /= rcolor.R;
	lcolor.G /= rcolor.G;
	lcolor.B /= rcolor.B;
	return lcolor;
}

template <class Type>
const ImgClass::RGB<Type>
operator/(ImgClass::RGB<Type> lcolor, const Type& rvalue)
{
	lcolor.R /= rvalue;
	lcolor.G /= rvalue;
	lcolor.B /= rvalue;
	return lcolor;
}

template <class Type>
const ImgClass::RGB<Type>
operator/(const Type& lvalue, ImgClass::RGB<Type> rcolor)
{
	rcolor.R = lvalue / rcolor.R;
	rcolor.G = lvalue / rcolor.G;
	rcolor.B = lvalue / rcolor.B;
	return rcolor;
}


// Comparator
template <class Type>
bool
operator==(const ImgClass::RGB<Type>& lcolor, const ImgClass::RGB<Type>& rcolor)
{
	if (lcolor.R == rcolor.R
	    && lcolor.G == rcolor.G
	    && lcolor.B == rcolor.B) {
		return true;
	} else {
		return false;
	}
}

template <class Type>
bool
operator!=(const ImgClass::RGB<Type>& lcolor, const ImgClass::RGB<Type>& rcolor)
{
	if (lcolor.R != rcolor.R
	    || lcolor.G != rcolor.G
	    || lcolor.B != rcolor.B) {
		return true;
	} else {
		return false;
	}
}

