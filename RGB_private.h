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
	template <class RT>
	RGB<T>::RGB(const RT& value)
	{
		R = value;
		G = value;
		B = value;
	}




	// Operators
	template <class T>
	template <class RT>
	RGB<T> &
	RGB<T>::operator=(const RT& value)
	{
		R = value;
		G = value;
		B = value;
	}

	template <class T>
	template <class RT>
	RGB<T> &
	RGB<T>::operator=(const RGB<RT>& color)
	{
		R = color.R;
		G = color.G;
		B = color.B;
	}

	template <class T>
	RGB<T> &
	RGB<T>::operator+=(const RGB<T>& color)
	{
		R += color.R;
		G += color.G;
		B += color.B;
	}

	template <class T>
	RGB<T> &
	RGB<T>::operator+=(const T& value)
	{
		R += value;
		G += value;
		B += value;
	}

	template <class T>
	RGB<T> &
	RGB<T>::operator-=(const RGB<T>& color)
	{
		R -= color.R;
		G -= color.G;
		B -= color.B;
	}

	template <class T>
	RGB<T> &
	RGB<T>::operator-=(const T& value)
	{
		R -= value;
		G -= value;
		B -= value;
	}

	template <class T>
	RGB<T> &
	RGB<T>::operator*=(const RGB<T>& color)
	{
		R *= color.R;
		G *= color.G;
		B *= color.B;
	}

	template <class T>
	RGB<T> &
	RGB<T>::operator*=(const T& value)
	{
		R *= value;
		G *= value;
		B *= value;
	}

	template <class T>
	RGB<T> &
	RGB<T>::operator/=(const RGB<T>& color)
	{
		R /= color.R;
		G /= color.G;
		B /= color.B;
	}

	template <class T>
	RGB<T> &
	RGB<T>::operator/=(const T& value)
	{
		R /= value;
		G /= value;
		B /= value;
	}
}


// ----- Global Operators -----

// Arithmetic

template <class Type>
ImgClass::RGB<Type> &
operator+(ImgClass::RGB<Type> color)
{
	return color;
}

template <class Type>
ImgClass::RGB<Type> &
operator-(ImgClass::RGB<Type> color)
{
	color.R = -color.R;
	color.G = -color.G;
	color.B = -color.B;
	return color;
}


template <class Type>
ImgClass::RGB<Type> &
operator+(ImgClass::RGB<Type> lcolor, const ImgClass::RGB<Type>& rcolor)
{
	lcolor.R += rcolor.R;
	lcolor.G += rcolor.G;
	lcolor.B += rcolor.B;
	return lcolor;
}

template <class Type>
ImgClass::RGB<Type> &
operator+(ImgClass::RGB<Type> lcolor, const Type& rvalue)
{
	lcolor.R += rvalue;
	lcolor.G += rvalue;
	lcolor.B += rvalue;
	return lcolor;
}

template <class Type>
ImgClass::RGB<Type> &
operator+(const Type& lvalue, ImgClass::RGB<Type> rcolor)
{
	rcolor.R += lvalue;
	rcolor.G += lvalue;
	rcolor.B += lvalue;
	return rcolor;
}

template <class Type>
ImgClass::RGB<Type> &
operator-(ImgClass::RGB<Type> lcolor, const ImgClass::RGB<Type>& rcolor)
{
	lcolor.R -= rcolor.R;
	lcolor.G -= rcolor.G;
	lcolor.B -= rcolor.B;
	return lcolor;
}

template <class Type>
ImgClass::RGB<Type> &
operator-(ImgClass::RGB<Type> lcolor, const Type& rvalue)
{
	lcolor.R -= rvalue;
	lcolor.G -= rvalue;
	lcolor.B -= rvalue;
	return lcolor;
}

template <class Type>
ImgClass::RGB<Type> &
operator-(const Type& lvalue, ImgClass::RGB<Type> rcolor)
{
	rcolor.R = lvalue - rcolor.R;
	rcolor.G = lvalue - rcolor.G;
	rcolor.B = lvalue - rcolor.B;
	return rcolor;
}

template <class Type>
ImgClass::RGB<Type> &
operator*(ImgClass::RGB<Type> lcolor, const ImgClass::RGB<Type>& rcolor)
{
	lcolor.R *= rcolor.R;
	lcolor.G *= rcolor.G;
	lcolor.B *= rcolor.B;
	return lcolor;
}

template <class Type>
ImgClass::RGB<Type> &
operator*(ImgClass::RGB<Type> lcolor, const Type& rvalue)
{
	lcolor.R *= rvalue;
	lcolor.G *= rvalue;
	lcolor.B *= rvalue;
	return lcolor;
}

template <class Type>
ImgClass::RGB<Type> &
operator*(const Type& lvalue, ImgClass::RGB<Type> rcolor)
{
	rcolor.R *= lvalue;
	rcolor.G *= lvalue;
	rcolor.B *= lvalue;
	return rcolor;
}

template <class Type>
ImgClass::RGB<Type> &
operator/(ImgClass::RGB<Type> lcolor, const ImgClass::RGB<Type>& rcolor)
{
	lcolor.R /= rcolor.R;
	lcolor.G /= rcolor.G;
	lcolor.B /= rcolor.B;
	return lcolor;
}

template <class Type>
ImgClass::RGB<Type> &
operator/(ImgClass::RGB<Type> lcolor, const Type& rvalue)
{
	lcolor.R /= rvalue;
	lcolor.G /= rvalue;
	lcolor.B /= rvalue;
	return lcolor;
}

template <class Type>
ImgClass::RGB<Type> &
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

