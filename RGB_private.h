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




	// Operators
	template <class T>
	RGB<T>::operator T() const // return intensity
	{
		const double yum_y_red = 0.299;
		const double yum_y_green = 0.587;
		const double yum_y_blue = 0.114;

		return yum_y_red * R + yum_y_green * G + yum_y_blue * B;
	}

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


	template <class T>
	const RGB<T>
	RGB<T>::operator+(const RGB<T>& rcolor) const
	{
		RGB<T> color;
		color.R += rcolor.R;
		color.G += rcolor.G;
		color.B += rcolor.B;
		return color;
	}

	template <class T>
	const RGB<T>
	RGB<T>::operator+(const T& rvalue) const
	{
		RGB<T> color;
		color.R += rvalue;
		color.G += rvalue;
		color.B += rvalue;
		return color;
	}


	template <class T>
	const RGB<T>
	RGB<T>::operator-(const RGB<T>& rcolor) const
	{
		RGB<T> color;
		color.R -= rcolor.R;
		color.G -= rcolor.G;
		color.B -= rcolor.B;
		return color;
	}

	template <class T>
	const RGB<T>
	RGB<T>::operator-(const T& rvalue) const
	{
		RGB<T> color;
		color.R -= rvalue;
		color.G -= rvalue;
		color.B -= rvalue;
		return color;
	}
}




// ----- Global Operators -----

// Arithmetic

template <class Type>
const ImgClass::RGB<Type>
operator+(ImgClass::RGB<Type> rcolor)
{
	return rcolor;
}

template <class Type>
const ImgClass::RGB<Type>
operator-(ImgClass::RGB<Type> rcolor)
{
	rcolor.R = -rcolor.R;
	rcolor.G = -rcolor.G;
	rcolor.B = -rcolor.B;
	return rcolor;
}


template <class Type>
const ImgClass::RGB<Type>
operator*(const ImgClass::RGB<Type>& lcolor, const ImgClass::RGB<Type>& rcolor)
{
	ImgClass::RGB<Type> color;
	color.R = lcolor.R * rcolor.R;
	color.G = lcolor.G * rcolor.G;
	color.B = lcolor.B * rcolor.B;
	return color;
}

template <class Type>
const ImgClass::RGB<Type>
operator*(const ImgClass::RGB<Type>& lcolor, const Type& rvalue)
{
	ImgClass::RGB<Type> color;
	color.R = lcolor.R * rvalue;
	color.G = lcolor.G * rvalue;
	color.B = lcolor.B * rvalue;
	return color;
}

template <class Type>
const ImgClass::RGB<Type>
operator*(const Type& lvalue, const ImgClass::RGB<Type>& rcolor)
{
	ImgClass::RGB<Type> color;
	color.R = lvalue * rcolor.R;
	color.G = lvalue * rcolor.G;
	color.B = lvalue * rcolor.B;
	return color;
}

template <class Type>
const ImgClass::RGB<Type>
operator/(const ImgClass::RGB<Type>& lcolor, const ImgClass::RGB<Type>& rcolor)
{
	ImgClass::RGB<Type> color;
	color.R = lcolor.R / rcolor.R;
	color.G = lcolor.G / rcolor.G;
	color.B = lcolor.B / rcolor.B;
	return color;
}

template <class Type>
const ImgClass::RGB<Type>
operator/(const ImgClass::RGB<Type>& lcolor, const Type& rvalue)
{
	ImgClass::RGB<Type> color;
	color.R = lcolor.R / rvalue;
	color.G = lcolor.G / rvalue;
	color.B = lcolor.B / rvalue;
	return color;
}

template <class Type>
const ImgClass::RGB<Type>
operator/(const Type& lvalue, const ImgClass::RGB<Type>& rcolor)
{
	ImgClass::RGB<Type> color;
	color.R = lvalue / rcolor.R;
	color.G = lvalue / rcolor.G;
	color.B = lvalue / rcolor.B;
	return color;
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

