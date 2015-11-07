#include <algorithm>
#include <cmath>
#include <iostream>
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
	_MV_forward = false;
}

template <class T>
BlockMatching<T>::BlockMatching(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const int BlockSize)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0;
	_MV_forward = false;
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
	_MV_forward = false;
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
	_MV_forward = copy._MV_forward;

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
	_MV_forward = false;
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
	_cells_height = 0;
	_MV_forward = false;
	if (image_prev == nullptr) {
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




// ----- Algorithm -----

template <class T>
void
BlockMatching<T>::block_matching(const int search_range, bool normalize)
{
	T (BlockMatching<T>::*SAD_func)(int, int, int, int, int, int);

	if (this->isNULL()) {
		return;
	}
	if (normalize) {
		SAD_func = &BlockMatching<T>::NSAD;
	} else {
		SAD_func = &BlockMatching<T>::SAD;
	}
	// Initialize
	_motion_vector.reset(_cells_width, _cells_height);
	// Compute Motion Vectors
#pragma omp parallel for
	for (int Y_b = 0; Y_b < _cells_height; Y_b++) {
		int y_b = Y_b * _block_size;
		for (int X_b = 0; X_b < _cells_width; X_b++) {
			int x_b = X_b * _block_size;
			VECTOR_2D<double> MV(.0, .0);
			T SAD_min;
			int x_start, x_end;
			int y_start, y_end;
			// Compute start and end coordinates
			if (search_range < 0) {
				x_start = 1 - _block_size;
				x_end = _width - 1;
				y_start = 1 - _block_size;
				y_end = _height - 1;
			} else {
				x_start = std::max(x_b - (search_range / 2), 1 - _block_size);
				x_end = std::min(x_b + search_range / 2, _width - 1);
				y_start = std::max(y_b - (search_range / 2), 1 - _block_size);
				y_end = std::min(y_b + search_range / 2, _height - 1);
			}
			SAD_min = (this->*SAD_func)(
			    x_start, y_start,
			    x_b, y_b,
			    _block_size, _block_size);
			for (int y = y_start; y <= y_end; y++) {
				for (int x = x_start; x <= x_end; x++) {
					VECTOR_2D<double> v_tmp((double)x - x_b, (double) y - y_b);
					T SAD_tmp = (this->*SAD_func)(
					    x, y,
					    x_b, y_b,
					    _block_size, _block_size);
					if (SAD_tmp < SAD_min) {
						SAD_min = SAD_tmp;
						MV = v_tmp;
					} else if (fabs(SAD_tmp - SAD_min) < 1E-6) {
						if (Vector_2D::norm(MV) > Vector_2D::norm(v_tmp)) {
							SAD_min = SAD_tmp;
							MV = v_tmp;
						}
					}
				}
			}
			_motion_vector.ref(X_b, Y_b) = MV;
		}
	}
	_MV_forward = false;
}

template <class T>
void
BlockMatching<T>::block_matching_forward(const int search_range, bool normalize)
{
	T (BlockMatching<T>::*SAD_func)(int, int, int, int, int, int);

	if (this->isNULL()) {
		return;
	}
	if (normalize) {
		SAD_func = &BlockMatching<T>::NSAD;
	} else {
		SAD_func = &BlockMatching<T>::SAD;
	}
	// Initialize
	_motion_vector.reset(_cells_width, _cells_height);
	// Compute Motion Vectors
#pragma omp parallel for
	for (int Y_b = 0; Y_b < _cells_height; Y_b++) {
		int y_b = Y_b * _block_size;
		for (int X_b = 0; X_b < _cells_width; X_b++) {
			int x_b = X_b * _block_size;
			VECTOR_2D<double> MV(.0, .0);
			T SAD_min;
			int x_start, x_end;
			int y_start, y_end;
			// Compute start and end coordinates
			if (search_range < 0) {
				x_start = 1 - _block_size;
				x_end = _width - 1;
				y_start = 1 - _block_size;
				y_end = _height - 1;
			} else {
				x_start = std::max(x_b - (search_range / 2), 1 - _block_size);
				x_end = std::min(x_b + search_range / 2, _width - 1);
				y_start = std::max(y_b - (search_range / 2), 1 - _block_size);
				y_end = std::min(y_b + search_range / 2, _height - 1);
			}
			SAD_min = (this->*SAD_func)(
			    x_b, y_b,
			    x_start, y_start,
			    _block_size, _block_size);
			for (int y = y_start; y <= y_end; y++) {
				for (int x = x_start; x <= x_end; x++) {
					VECTOR_2D<double> v_tmp((double)x - x_b, (double)y - y_b);
					T SAD_tmp = (this->*SAD_func)(
					    x_b, y_b,
					    x, y,
					    _block_size, _block_size);
					if (SAD_tmp < SAD_min) {
						SAD_min = SAD_tmp;
						MV = v_tmp;
					} else if (fabs(SAD_tmp - SAD_min) < 1E-6) {
						if (Vector_2D::norm(MV) > Vector_2D::norm(v_tmp)) {
							SAD_min = SAD_tmp;
							MV = v_tmp;
						}
					}
				}
			}
			_motion_vector.ref(X_b, Y_b) = MV;
		}
	}
	_MV_forward = true;
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

template <class T>
T
BlockMatching<T>::NSAD(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height)
{
	T sad = 0;
	ImgVector<T>* norm_image_prev = nullptr;
	ImgVector<T>* norm_image_next = nullptr;

	norm_image_prev = _image_prev.crop(x_prev, y_prev, x_prev + block_width, y_prev + block_height);
	norm_image_next = _image_next.crop(x_next, y_next, x_next + block_width, y_next + block_height);
	// Contrast stretching for normalization
	T max = std::max(_image_prev.max(), _image_next.max());
	norm_image_prev->contrast_stretching(0, max);
	norm_image_next->contrast_stretching(0, max);
	// Compute SAD
	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			sad = sad + fabs((double)norm_image_next->get_zeropad(x, y) - norm_image_prev->get_zeropad(x, y));
		}
	}
	delete norm_image_prev;
	delete norm_image_next;
	return sad;
}

