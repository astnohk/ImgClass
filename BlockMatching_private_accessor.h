#include <algorithm>
#include <cmath>
#include <iostream>
#include <new>
#include <stdexcept>




// ----- Get data -----
template <class T>
int
BlockMatching<T>::width(void) const
{
	return _width;
}

template <class T>
int
BlockMatching<T>::height(void) const
{
	return _height;
}

template <class T>
int
BlockMatching<T>::block_size(void) const
{
	return _block_size;
}

template <class T>
int
BlockMatching<T>::vector_width(void) const
{
	return _cells_width;
}

template <class T>
int
BlockMatching<T>::vector_height(void) const
{
	return _cells_height;
}

template <class T>
bool
BlockMatching<T>::isNULL(void)
{
	if (_width <= 0 || _height <= 0) {
		return true;
	} else {
		return false;
	}
}

template <class T>
bool
BlockMatching<T>::isForward(void)
{
	return _MV_forward;
}




// ----- Get reference -----
template <class T>
ImgVector<VECTOR_2D<double> > &
BlockMatching<T>::data(void)
{
	if (_motion_vector.isNULL()) {
		this->block_matching();
	}
	return _motion_vector;
}

template <class T>
VECTOR_2D<double> &
BlockMatching<T>::operator[](int n)
{
	if (_motion_vector.isNULL()) {
		this->block_matching();
	}
	return _motion_vector[n];
}

template <class T>
VECTOR_2D<double>&
BlockMatching<T>::at(int x, int y)
{
	if (_motion_vector.isNULL()) {
		this->block_matching();
	}
	return _motion_vector.at(x, y);
}




// ----- Get Vector Field data -----
template <class T>
const VECTOR_2D<double>
BlockMatching<T>::get(int x, int y)
{
	if (_motion_vector.isNULL()) {
		this->block_matching();
	}
	return _motion_vector.at(x, y);
}

