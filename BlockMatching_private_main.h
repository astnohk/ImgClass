#include <algorithm>
#include <cmath>
#include <iostream>
#include <new>
#include <stdexcept>

#if defined(_OPENMP)
#include <omp.h>
#endif




// ----- Algorithm -----

template <class T>
void
BlockMatching<T>::block_matching(const int search_range)
{
	T (BlockMatching<T>::*Corr_func)(int, int, int, int, int, int);

	if (this->isNULL()) {
		return;
	}
	// Initialize
	_motion_vector.reset(_cells_width, _cells_height);
	// Compute Motion Vectors
	Corr_func = &BlockMatching<T>::SAD;
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
			SAD_min = (this->*Corr_func)(
			    x_start, y_start,
			    x_b, y_b,
			    _block_size, _block_size);
			for (int y = y_start; y <= y_end; y++) {
				for (int x = x_start; x <= x_end; x++) {
					VECTOR_2D<double> v_tmp((double)x - x_b, (double) y - y_b);
					T SAD_tmp = (this->*Corr_func)(
					    x, y,
					    x_b, y_b,
					    _block_size, _block_size);
					if (SAD_tmp < SAD_min) {
						SAD_min = SAD_tmp;
						MV = v_tmp;
					} else if (fabs(SAD_tmp - SAD_min) < 1E-6) {
						if (Vector_2D::norm(MV) >= Vector_2D::norm(v_tmp)) {
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
BlockMatching<T>::block_matching_forward(const int search_range)
{
	T (BlockMatching<T>::*Corr_func)(int, int, int, int, int, int);

	if (this->isNULL()) {
		return;
	}
	// Initialize
	_motion_vector.reset(_cells_width, _cells_height);
	// Compute Motion Vectors
	Corr_func = &BlockMatching<T>::SAD;
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
			SAD_min = fabs((this->*Corr_func)(
			    x_b, y_b,
			    x_start, y_start,
			    _block_size, _block_size));
			for (int y = y_start; y <= y_end; y++) {
				for (int x = x_start; x <= x_end; x++) {
					VECTOR_2D<double> v_tmp((double)x - x_b, (double)y - y_b);
					T SAD_tmp = fabs((this->*Corr_func)(
					    x_b, y_b,
					    x, y,
					    _block_size, _block_size));
					if (SAD_tmp < SAD_min) {
						SAD_min = SAD_tmp;
						MV = v_tmp;
					} else if (fabs(SAD_tmp - SAD_min) < 1E-6) {
						if (Vector_2D::norm(MV) >= Vector_2D::norm(v_tmp)) {
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
void
BlockMatching<T>::block_matching_SAD_ZNCC(const int search_range)
{
	const double coeff_SAD = 1.0;
	const double coeff_ZNCC = 10.0;
	double (BlockMatching<T>::*SAD_func)(int, int, int, int, int, int) = &BlockMatching<T>::MAD;
	double (BlockMatching<T>::*NCC_func)(int, int, int, int, int, int) = &BlockMatching<T>::ZNCC;

	if (this->isNULL()) {
		return;
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
			T ZNCC_max;
			T E_min;
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
			ZNCC_max = (this->*NCC_func)(
			    x_start, y_start,
			    x_b, y_b,
			    _block_size, _block_size);
			E_min = coeff_SAD * SAD_min + coeff_ZNCC * (1.0 - ZNCC_max);
			for (int y = y_start; y <= y_end; y++) {
				for (int x = x_start; x <= x_end; x++) {
					VECTOR_2D<double> v_tmp((double)x - x_b, (double) y - y_b);
					T SAD_min = (this->*SAD_func)(
					    x, y,
					    x_b, y_b,
					    _block_size, _block_size);
					T ZNCC_max = (this->*NCC_func)(
					    x, y,
					    x_b, y_b,
					    _block_size, _block_size);
					T E_tmp = coeff_SAD * SAD_min + coeff_ZNCC * (1.0 - ZNCC_max);
					if (E_tmp < E_min) {
						E_min = E_tmp;
						MV = v_tmp;
					} else if (fabs(E_tmp - E_min) < 1E-6) {
						if (Vector_2D::norm(MV) >= Vector_2D::norm(v_tmp)) {
							E_min = E_tmp;
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
BlockMatching<T>::MAD(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height)
{
	T sad = 0;

	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			sad = sad + fabs((double)_image_next.get_zeropad(x_next + x, y_next + y) - _image_prev.get_zeropad(x_prev + x, y_prev + y));
		}
	}
	return sad / (double)(block_width * block_height);
}

template <class T>
T
BlockMatching<T>::NCC(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height)
{
	T prod = 0;
	T norm_prev = 0;
	T norm_next = 0;

	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			prod += _image_next.get_zeropad(x_next + x, y_next + y) * _image_prev.get_zeropad(x_prev + x, y_prev + y);
			norm_prev += _image_prev.get_zeropad(x_prev + x, y_prev + y) * _image_prev.get_zeropad(x_prev + x, y_prev + y);
			norm_next += _image_next.get_zeropad(x_next + x, y_next + y) * _image_next.get_zeropad(x_next + x, y_next + y);
		}
	}
	return prod / sqrt(norm_prev) / sqrt(norm_next);
}

template <class T>
T
BlockMatching<T>::ZNCC(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height)
{
	double N = block_width * block_height;
	T sum_prev = 0;
	T sum_next = 0;
	T sum_sq_prev = 0;
	T sum_sq_next = 0;
	T sum_sq_prev_next = 0;

	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			sum_prev += _image_prev.get_zeropad(x_prev + x, y_prev + y);
			sum_next += _image_next.get_zeropad(x_next + x, y_next + y);
			sum_sq_prev += _image_prev.get_zeropad(x_prev + x, y_prev + y) * _image_prev.get_zeropad(x_prev + x, y_prev + y);
			sum_sq_next += _image_next.get_zeropad(x_next + x, y_next + y) * _image_next.get_zeropad(x_next + x, y_next + y);
			sum_sq_prev_next += _image_prev.get_zeropad(x_prev + x, y_prev + y) * _image_next.get_zeropad(x_next + x, y_next + y);
		}
	}
	return (N * sum_sq_prev_next - sum_prev * sum_next) / sqrt(N * sum_sq_prev - sum_prev * sum_prev) / sqrt(N * sum_sq_next - sum_next * sum_next);
}

