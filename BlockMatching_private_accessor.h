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
BlockMatching<T>::vector_field_width(void) const
{
	return _cells_width;
}

template <class T>
int
BlockMatching<T>::vector_field_height(void) const
{
	return _cells_height;
}

template <class T>
bool
BlockMatching<T>::isNULL(void) const
{
	if (_width <= 0 || _height <= 0) {
		return true;
	} else {
		return false;
	}
}




// ----- Get reference -----
template <class T>
ImgVector<Vector_ST<double> > &
BlockMatching<T>::ref_motion_vector_time(void)
{
	if (_motion_vector_time.isNULL()) {
		this->block_matching();
	}
	return _motion_vector_time;
}

template <class T>
ImgVector<VECTOR_2D<double> > &
BlockMatching<T>::ref_motion_vector_prev(void)
{
	if (_motion_vector_prev.isNULL()) {
		this->block_matching();
	}
	return _motion_vector_prev;
}

template <class T>
ImgVector<VECTOR_2D<double> > &
BlockMatching<T>::ref_motion_vector_next(void)
{
	if (_motion_vector_next.isNULL()) {
		this->block_matching();
	}
	return _motion_vector_next;
}


template <class T>
VECTOR_2D<double> &
BlockMatching<T>::operator[](int n)
{
	if (_motion_vector_prev.isNULL()) {
		this->block_matching();
	}
	return _motion_vector_prev[n];
}


template <class T>
Vector_ST<double>&
BlockMatching<T>::at(int x, int y)
{
	if (_motion_vector_time.isNULL()) {
		this->block_matching();
	}
	assert(0 <= x && x < _width
	    && 0 <= y && y < _height);
	return _motion_vector_time.at(int(floor(x / _block_size)), int(floor(y / _block_size)));
}

template <class T>
VECTOR_2D<double>&
BlockMatching<T>::at_prev(int x, int y)
{
	if (_motion_vector_prev.isNULL()) {
		this->block_matching();
	}
	assert(0 <= x && x < _width
	    && 0 <= y && y < _height);
	return _motion_vector_prev.at(int(floor(x / _block_size)), int(floor(y / _block_size)));
}

template <class T>
VECTOR_2D<double>&
BlockMatching<T>::at_next(int x, int y)
{
	if (_motion_vector_next.isNULL()) {
		this->block_matching();
	}
	assert(0 <= x && x < _width
	    && 0 <= y && y < _height);
	return _motion_vector_next.at(int(floor(x / _block_size)), int(floor(y / _block_size)));
}

template <class T>
Vector_ST<double>&
BlockMatching<T>::at_block(int x, int y)
{
	if (_motion_vector_time.isNULL()) {
		this->block_matching();
	}
	assert(0 <= x && x < _cells_width
	    && 0 <= y && y < _cells_height);
	return _motion_vector_time.at(x, y);
}

template <class T>
VECTOR_2D<double>&
BlockMatching<T>::at_block_prev(int x, int y)
{
	if (_motion_vector_prev.isNULL()) {
		this->block_matching();
	}
	assert(0 <= x && x < _cells_width
	    && 0 <= y && y < _cells_height);
	return _motion_vector_prev.at(x, y);
}

template <class T>
VECTOR_2D<double>&
BlockMatching<T>::at_block_next(int x, int y)
{
	if (_motion_vector_next.isNULL()) {
		this->block_matching();
	}
	assert(0 <= x && x < _cells_width
	    && 0 <= y && y < _cells_height);
	return _motion_vector_next.at(x, y);
}




// ----- Get Vector Field data -----
template <class T>
const Vector_ST<double>
BlockMatching<T>::get(int x, int y)
{
	if (_motion_vector_time.isNULL()) {
		this->block_matching();
	}
	assert(0 <= x && x < _width
	    && 0 <= y && y < _height);
	return _motion_vector_time.get(int(floor(x / _block_size)), int(floor(y / _block_size)));
}

template <class T>
const VECTOR_2D<double>
BlockMatching<T>::get_prev(int x, int y)
{
	if (_motion_vector_prev.isNULL()) {
		this->block_matching();
	}
	assert(0 <= x && x < _width
	    && 0 <= y && y < _height);
	return _motion_vector_prev.get(int(floor(x / _block_size)), int(floor(y / _block_size)));
}

template <class T>
const VECTOR_2D<double>
BlockMatching<T>::get_next(int x, int y)
{
	if (_motion_vector_next.isNULL()) {
		this->block_matching();
	}
	assert(0 <= x && x < _width
	    && 0 <= y && y < _height);
	return _motion_vector_next.get(int(floor(x / _block_size)), int(floor(y / _block_size)));
}

template <class T>
const Vector_ST<double>
BlockMatching<T>::get_block(int x, int y)
{
	if (_motion_vector_time.isNULL()) {
		this->block_matching();
	}
	assert(0 <= x && x < _cells_width
	    && 0 <= y && y < _cells_height);
	return _motion_vector_time.get(x, y);
}

template <class T>
const VECTOR_2D<double>
BlockMatching<T>::get_block_prev(int x, int y)
{
	if (_motion_vector_prev.isNULL()) {
		this->block_matching();
	}
	assert(0 <= x && x < _cells_width
	    && 0 <= y && y < _cells_height);
	return _motion_vector_prev.get(x, y);
}

template <class T>
const VECTOR_2D<double>
BlockMatching<T>::get_block_next(int x, int y)
{
	if (_motion_vector_next.isNULL()) {
		this->block_matching();
	}
	assert(0 <= x && x < _cells_width
	    && 0 <= y && y < _cells_height);
	return _motion_vector_next.get(x, y);
}

