#include <algorithm>
#include <cmath>
#include <cstdio>
#include <new>
#include <stdexcept>

#if defined(_OPENMP)
#include <omp.h>
#endif



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
	if (image_prev.isNULL()) {
		throw std::invalid_argument("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : ImgVector<T>& image_prev");
	} else if (image_next.isNULL()) {
		throw std::invalid_argument("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : ImgVector<T>& image_next");
	} else if (image_prev.width() != image_next.width()) {
		throw std::invalid_argument("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : width of image_prev and image_next not match");
	} else if (image_prev.height() != image_next.height()) {
		throw std::invalid_argument("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : height of image_prev and image_next not match");
	} else if (BlockSize < 0) {
		throw std::out_of_range("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : BlockSize");
	}

	_width = image_prev.width();
	_height = image_prev.height();
	_block_size = BlockSize;
	_image_prev.copy(image_prev);
	_image_next.copy(image_next);
	_motion_vector.reset(0, 0);

	_cells_width = (int)ceil((double)_width / BlockSize);
	_cells_height = (int)ceil((double)_height / BlockSize);
}

template <class T>
BlockMatching<T>::BlockMatching(const ImgVector<T>* image_prev, const ImgVector<T>* image_next, const int BlockSize)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0;
	if (image_prev == nullptr) {
		throw std::invalid_argument("const ImgVector<T>* image_prev");
	} else if (image_next == nullptr) {
		throw std::invalid_argument("const ImgVector<T>* image_next");
	} if (image_prev->isNULL()) {
		throw std::invalid_argument("const ImgVector<T>* image_prev");
	} if (image_next->isNULL()) {
		throw std::invalid_argument("const ImgVector<T>* image_next");
	} else if (image_prev.width() != image_next.width() || image_prev.height() != image_next.height()) {
		throw std::invalid_argument("width or height of image_prev and image_next not match");
	} else if (BlockSize < 0) {
		throw std::out_of_range("BlockSize");
	}

	_width = image_prev.width();
	_height = image_prev.height();
	_block_size = BlockSize;
	_image_prev.copy(image_prev);
	_image_next.copy(image_next);
	_motion_vector.reset(0, 0);

	_cells_width = (int)ceil((double)_width / BlockSize);
	_cells_height = (int)ceil((double)_height / BlockSize);
}

template <class T>
BlockMatching<T>::BlockMatching(const BlockMatching& copy)
{
	_width = copy._width;
	_height = copy._height;
	_block_size = copy._block_size;
	_cells_width = copy._cells_width;
	_cells_height = copy._cells_height;

	_image_prev.copy(copy._image_prev);
	_image_next.copy(copy._image_next);
	_motion_vector.copy(copy._motion_vector);
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
	if (image_prev.isNULL()) {
		throw std::invalid_argument("const ImgVector<T>* image_prev");
	} else if (image_next.isNULL()) {
		throw std::invalid_argument("const ImgVector<T>* image_next");
	} else if (image_prev.width() != image_next.width() || image_prev.height() != image_next.height()) {
		throw std::invalid_argument("width or height of image_prev and image_next not match");
	} else if (BlockSize < 0) {
		throw std::out_of_range("BlockSize");
	}

	_width = image_prev.width();
	_height = image_prev.height();
	_block_size = BlockSize;
	_image_prev.copy(image_prev);
	_image_next.copy(image_next);
	_motion_vector.reset(0, 0);

	_cells_width = (int)ceil((double)_width / BlockSize);
	_cells_height = (int)ceil((double)_height / BlockSize);
}

template <class T>
void
BlockMatching<T>::reset(const ImgVector<T>* image_prev, const ImgVector<T>* image_next, const int BlockSize)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0; if (image_prev == nullptr) {
		throw std::invalid_argument("void BlockMatching<T>::reset(const ImgVector<T>*, const ImgVector<T>*, const int) : const ImgVector<T>* image_prev");
	} else if (image_next == nullptr) {
		throw std::invalid_argument("void BlockMatching<T>::reset(const ImgVector<T>*, const ImgVector<T>*, const int) : const ImgVector<T>* image_next");
	} else if (image_prev->isNULL()) {
		throw std::invalid_argument("void BlockMatching<T>::reset(const ImgVector<T>*, const ImgVector<T>*, const int) : const ImgVector<T>* image_prev");
	} else if (image_next->isNULL()) {
		throw std::invalid_argument("void BlockMatching<T>::reset(const ImgVector<T>*, const ImgVector<T>*, const int) : const ImgVector<T>* image_next");
	} else if (image_prev->width() != image_next->width()) {
		throw std::invalid_argument("void BlockMatching<T>::reset(const ImgVector<T>*, const ImgVector<T>*, const int) : the width of image_prev and image_next not match");
	} else if (image_prev->height() != image_next->height()) {
		throw std::invalid_argument("void BlockMatching<T>::reset(const ImgVector<T>*, const ImgVector<T>*, const int) : the height of image_prev and image_next not match");
	} else if (BlockSize < 0) {
		throw std::out_of_range("void BlockMatching<T>::reset(const ImgVector<T>*, const ImgVector<T>*, const int) : BlockSize");
	}

	_width = image_prev->width();
	_height = image_prev->height();
	_block_size = BlockSize;
	_image_prev.copy(image_prev);
	_image_next.copy(image_next);
	_motion_vector.reset(0, 0);

	_cells_width = (int)ceil((double)_width / BlockSize);
	_cells_height = (int)ceil((double)_height / BlockSize);
}


template <class T>
void
BlockMatching<T>::block_matching(const int search_range)
{
	if (this->isNULL()) {
		return;
	} else if (_motion_vector.isNULL() == false) {
		return;
	}

	// Initialize
	_motion_vector.reset(_cells_width, _cells_height);

	int x_c, y_c;
#pragma omp parallel for private(x_c)
	for (y_c = 0; y_c < _cells_height; y_c++) {
		for (x_c = 0; x_c < _cells_width; x_c++) {
			_motion_vector.ref(x_c, y_c) = max_crosscorr(x_c * _block_size, y_c * _block_size, search_range);
		}
	}
}

template <class T>
VECTOR_2D<double>
BlockMatching<T>::max_crosscorr(const int x_prev, const int y_prev, const int search_range)
{
	VECTOR_2D<double> vector(.0, .0);
	T min;
	int x_start, x_end;
	int y_start, y_end;

	if (search_range < 0) {
		x_start = 1 - _block_size;
		x_end = _width - 1;
		y_start = 1 - _block_size;
		y_end = _height - 1;
	} else {
		x_start = std::max(x_prev - search_range / 2, 1 - _block_size);
		x_end = x_prev + search_range / 2;
		y_start = std::max(y_prev - search_range / 2, 1 - _block_size);
		y_end = y_prev + search_range / 2;
	}
	// Initialize temporal minimum of SAD
	min = this->SAD(x_prev, y_prev, x_start, y_start, _block_size, _block_size);
	// Search infimum of SAD
	for (int y = y_start; y <= y_end; y++) {
		for (int x = x_start; x <= x_end; x++) {
			T sad = this->SAD(x_prev, y_prev, x, y, _block_size, _block_size);

			if (sad < min) {
				min = sad;
				vector.x = (double)x - x_prev;
				vector.y = (double)y - y_prev;
			} else if (fabs(sad - min) < 1E-6) {
				if (Vector_2D::norm(vector) > (double)sqrt((x - x_prev) * (x - x_prev) + (y - y_prev) * (y - y_prev))) {
					min = sad;
					vector.x = (double)x - x_prev;
					vector.y = (double)y - y_prev;
				}
			}
		}
	}
	return vector;
}

template <class T>
T
BlockMatching<T>::SAD(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height)
{
	T sad = 0;

	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			sad = sad + fabs((double)_image_next.get_zeropad(x_next + x, y_next + y) - _image_prev.get_zeropad(x_prev + x, y_prev + y));
		}
	}
	return sad;
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

// Get reference
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
BlockMatching<T>::ref(int x, int y)
{
	if (_motion_vector.isNULL()) {
		this->block_matching();
	}
	return _motion_vector.ref(x, y);
}

// Get data
template <class T>
VECTOR_2D<double>
BlockMatching<T>::get(int x, int y)
{
	if (_motion_vector.isNULL()) {
		this->block_matching();
	}
	return _motion_vector.get(x, y);
}

