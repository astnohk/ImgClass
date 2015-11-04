#include <cmath>
#include <cstdio>
#include <new>
#include <stdexcept>

#include "BlockMatching.h"



template <class T>
BlockMatching<T>::BlockMatching(void)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0;
}

template <class T>
BlockMatching<T>::BlockMatching(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const int BlockSize)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0;
	if (img_prev.isNULL() || img_next.isNULL()) {
		throw std::invalid_argument("img_prev, img_next");
	} else if (img_prev.width() != img_next.width() || img_prev.height() != img_next.height()) {
		throw std::invalid_argument("width or height of (img_prev, img_next)");
	} else if (BlockSize < 0) {
		throw std::out_of_range("BlockSize");
	}

	_width = img_prev.width();
	_height = img_prev.height();
	_block_size = BlockSize;
	_image_prev.copy(image_prev);
	_image_next.copy(image_next);
	_motion_vector.reset(0, 0);

	_cells_width = (int)ceil(_width / BlockSize);
	_cells_height = (int)ceil(_height / BlockSize);
}

template <class T>
BlockMatching<T>::~BlockMatching(void)
{
}

template <class T>
void
BlockMatching<T>::reset(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const int BlockSize)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0;
	if (img_prev.isNULL() || img_next.isNULL()) {
		throw std::invalid_argument("img_prev, img_next");
	} else if (img_prev.width() != img_next.width() || img_prev.height() != img_next.height()) {
		throw std::invalid_argument("width or height of (img_prev, img_next)");
	} else if (BlockSize < 0) {
		throw std::out_of_range("BlockSize");
	}

	_width = img_prev.width();
	_height = img_prev.height();
	_block_size = BlockSize;
	_image_prev.copy(image_prev);
	_image_next.copy(image_next);
	_motion_vector.reset(0, 0);

	_cells_width = (int)ceil(_width / BlockSize);
	_cells_height = (int)ceil(_height / BlockSize);
}


template <class T>
void
BlockMatching<T>::block_matching(void)
{
	if (this->isNULL()) {
		return;
	} else if (_motion_vector.isNULL() == false) {
		return;
	}

	// Initialize
	_motion_vector.reset(_width, _height);

	for (int y_c = 0; y_c < _cells_height; y_c++) {
		for (int x_c = 0; x_c < _cells_width; x_c++) {
			_motion_vector.ref(x_c, y_c) = max_crosscorr(x_c * _block_size, y_c * _block_size);
		}
	}
}

template <class T>
VECTOR_2D<double> &
BlockMatching<T>::max_crosscorr(const int x_prev, const int y_prev)
{
	VECTOR_2D<double> vector;
	T max = 0;

	for (int y = 0; y < _cells_height; y++) {
		for (int x = 0; x < _cells_width; x++) {
			T tmp = this->MAD(x_prev, y_prev, x, y, _block_size, _block_size);
			if (tmp > max) {
				max = tmp;
				vector.x = x - x_prev;
				vector.y = y - y_prev;
			}
		}
	}
	return vector;
}

template <class T>
T
BlockMatching<T>::MAD(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height)
{
	T mad = 0;

	for (int m = 0; m < block_heihgt; m++) {
		for (int n = 0; n < block_width; n++) {
			for (int y = 0; y < block_height; y++) {
				for (int x = 0; x < block_width; x++) {
					mad = mad + abs(image_next.get(x_prev + x, y_prev + y) - image_prev.get(x_next + n, y_next + m));
				}
			}
		}
	}
	return mad / (block_width * block_height);
}


// Check state
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
		true;
	} else {
		false;
	}
}

// Get reference
template <class T>
ImgVector<VECTOR_2D<double> > &
BlockMatching<T>::data(void)
{
	return _motion_vector;
}

template <class T>
VECTOR_2D<double> &
BlockMatching<T>::operator[](int n)
{
	return _motion_vector[];
}

template <class T>
VECTOR_2D<double>&
BlockMatching<T>::ref(int x, int y)
{
	return _motion_vector.ref(x, y);
}

// Get data
template <class T>
VECTOR_2D<double>
BlockMatching<T>::get(int x, int y) const
{
	return _motion_vector.get(x, y);
}

