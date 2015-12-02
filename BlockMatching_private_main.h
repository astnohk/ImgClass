#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
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
	const double coeff_SAD = 1.0;
	const double coeff_ZNCC = 10.0;
	double (BlockMatching<T>::*SAD_func)(int, int, int, int, int, int) = &BlockMatching<T>::MAD;
	double (BlockMatching<T>::*NCC_func)(int, int, int, int, int, int) = &BlockMatching<T>::ZNCC;

	ImgVector<bool> estimated(_cells_width, _cells_height, false);
	std::list<VECTOR_2D<int> > flat_blocks;
	std::list<VECTOR_2D<int> >::iterator itr_flat_blocks;

	if (this->isNULL()) {
		std::cerr << "void BlockMatching<T>::block_matching(const int) : _block_size < 0" << std::endl;
		throw std::logic_error("void BlockMatching<T>::block_matching(const int) : this is NULL");
	} else if (_block_size < 0) {
		std::cerr << "void BlockMatching<T>::block_matching(const int) : _block_size < 0" << std::endl;
		throw std::logic_error("void BlockMatching<T>::block_matching(const int) : _block_size < 0");
	}
	// Initialize
	_motion_vector.reset(_cells_width, _cells_height);
	ImgVector<VECTOR_2D<double> >* grad_prev = this->grad_prev(0, 0, _image_prev.width(), _image_prev.height());
	ImgVector<double> grad_prev_x(grad_prev->width(), grad_prev->height());
	ImgVector<double> grad_prev_y(grad_prev->width(), grad_prev->height());
	for (int i = 0; i < grad_prev->size(); i++) {
		grad_prev_x[i] = grad_prev->get(i).x;
		grad_prev_y[i] = grad_prev->get(i).y;
	}
	delete grad_prev;
	grad_prev = nullptr;
	T grad_prev_max_min_x = grad_prev_x.max() - grad_prev_x.min();
	T grad_prev_max_min_y = grad_prev_y.max() - grad_prev_y.min();
	// Compute Motion Vectors
#pragma omp parallel for
	for (int Y_b = 0; Y_b < _cells_height; Y_b++) {
		int y_b = Y_b * _block_size;
		for (int X_b = 0; X_b < _cells_width; X_b++) {
			int x_b = X_b * _block_size;
			// Check the block can be estimated
			ImgVector<VECTOR_2D<double> >* grad_prev_local = this->grad_prev(x_b, y_b, _block_size, _block_size);
			ImgVector<double> grad_prev_local_x(_block_size, _block_size);
			ImgVector<double> grad_prev_local_y(_block_size, _block_size);
			for (int i = 0; i < grad_prev_local->size(); i++) {
				grad_prev_local_x[i] = grad_prev_local->get(i).x;
				grad_prev_local_y[i] = grad_prev_local->get(i).y;
			}
			delete grad_prev_local;
			grad_prev_local = nullptr;
			if ((grad_prev_local_x.max() - grad_prev_local_x.min()) / grad_prev_max_min_x < 0.1
			    || (grad_prev_local_y.max() - grad_prev_local_y.min()) / grad_prev_max_min_y < 0.1) {
				if (estimated.get_zeropad(X_b - 1, Y_b) || estimated.get_zeropad(X_b, Y_b - 1)) {
					flat_blocks.push_front(VECTOR_2D<int>(X_b, Y_b));
				} else {
					flat_blocks.push_back(VECTOR_2D<int>(X_b, Y_b));
				}
				_motion_vector.at(X_b, Y_b) = 0;
				continue;
			}
			estimated.at(X_b, Y_b) = true;

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
			_motion_vector.at(X_b, Y_b) = MV;
		}
	}
	// Interpolate Motion Vector on skipped blocks
	vector_interpolation(flat_blocks, &estimated);
}


template <class T>
void
BlockMatching<T>::vector_interpolation(const std::list<VECTOR_2D<int> >& flat_blocks, ImgVector<bool>* estimated)
{
	for (std::list<VECTOR_2D<int> >::const_iterator itr_flat_blocks = flat_blocks.begin();
	    itr_flat_blocks != flat_blocks.end();
	    ++itr_flat_blocks) {
		int X_b = itr_flat_blocks->x;
		int Y_b = itr_flat_blocks->y;
		int x_b = X_b * _block_size;
		int y_b = Y_b * _block_size;
		T MAD_min = _image_prev.max();
		T MAD_tmp;

		// Compare 4 adjacent
		if (estimated->get_mirror(X_b, Y_b - 1)
		    && (MAD_tmp = this->MAD(x_b, y_b, x_b, y_b - _block_size, _block_size, _block_size, _image_prev, _image_prev)) < MAD_min) {
			_motion_vector.at(X_b, Y_b) = _motion_vector.get_mirror(X_b, Y_b - 1);
			MAD_min = MAD_tmp;
		}
		if (estimated->get_mirror(X_b - 1, Y_b)
		    && (MAD_tmp = this->MAD(x_b, y_b, x_b - _block_size, y_b, _block_size, _block_size, _image_prev, _image_prev)) < MAD_min) {
			_motion_vector.at(X_b, Y_b) = _motion_vector.get_mirror(X_b - 1, Y_b);
			MAD_min = MAD_tmp;
		}
		if (estimated->get_mirror(X_b, Y_b + 1)
		    && (MAD_tmp = this->MAD(x_b, y_b, x_b, y_b + _block_size, _block_size, _block_size, _image_prev, _image_prev)) < MAD_min) {
			_motion_vector.at(X_b, Y_b) = _motion_vector.get_mirror(X_b, Y_b + 1);
			MAD_min = MAD_tmp;
		}
		if (estimated->get_mirror(X_b + 1, Y_b)
		    && (MAD_tmp = this->MAD(x_b, y_b, x_b + _block_size, y_b, _block_size, _block_size, _image_prev, _image_prev)) < MAD_min) {
			_motion_vector.at(X_b, Y_b) = _motion_vector.get_mirror(X_b + 1, Y_b);
			MAD_min = MAD_tmp;
		}
		estimated->at(X_b, Y_b) = true;
	}
}

template <class T>
ImgVector<VECTOR_2D<double> > *
BlockMatching<T>::grad_prev(const int top_left_x, const int top_left_y, const int crop_width, const int crop_height)
{
	ImgVector<VECTOR_2D<double> >* gradients = new ImgVector<VECTOR_2D<double> >(crop_width, crop_height);

	for (int y = 0; y < crop_height; y++) {
		for (int x = 0; x < crop_width; x++) {
			gradients->at(x, y).x =
			    _image_prev.get_mirror(top_left_x + x + 1, top_left_y + y)
			    - _image_prev.get_mirror(top_left_x + x, top_left_y + y);
			gradients->at(x, y).y =
			    _image_prev.get_mirror(top_left_x + x, top_left_y + y + 1)
			    - _image_prev.get_mirror(top_left_x + x, top_left_y + y);
		}
	}
	return gradients;
}




template <class T>
void
BlockMatching<T>::block_matching_subset(const ImgVector<int>* region_map, const int search_range)
{
	const double coeff_SAD = 1.0;
	const double coeff_ZNCC = 10.0;
	double (BlockMatching<T>::*SAD_func)(int, int, int, int, int, int) = &BlockMatching<T>::MAD;
	double (BlockMatching<T>::*NCC_func)(int, int, int, int, int, int) = &BlockMatching<T>::ZNCC;

	std::list<VECTOR_2D<int> > regions;

	if (this->isNULL()) {
		std::cerr << "void BlockMatching<T>::block_matching_region(const ImgVector<int>*, const int) : this is NULL" << std::endl;
		throw std::logic_error("void BlockMatching<T>::block_matching_region(const ImgVector<int>*, const int) : this is NULL");
	} else if (region_map == nullptr) {
		std::cerr << "void BlockMatching<T>::block_matching_region(const ImgVector<int>*, const int)" << std::endl;
		throw std::invalid_argument("const ImgVector<int>* region_map");
	}
	// Initialize
	_block_size = -1; // Do not use block_size
	// MV are expanded on entire image pixels
	_cells_width = _width;
	_cells_height = _height;
	_motion_vector.reset(_width, _height);
	// Make connected region list
	std::list<VECTOR_2D<int> > regions = get_connected_region_list(region_map);
	// Compute Motion Vectors
#pragma omp parallel for schedule(dynamic)
	for (int n = 0; n < ; n++) {
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
			_motion_vector.at(X_b, Y_b) = MV;
		}
	}
}




template <class T>
const T
BlockMatching<T>::MAD(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height, const ImgVector<T>& img_prev, const ImgVector<T>& img_next)
{
	T sad = 0;

	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			sad = sad + fabs((double)img_next.get_zeropad(x_next + x, y_next + y) - img_prev.get_zeropad(x_prev + x, y_prev + y));
		}
	}
	return sad / double(block_width * block_height);
}

template <class T>
const T
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
const T
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
const T
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
	return prod / (sqrt(norm_prev) * sqrt(norm_next) + 1E-10);
}

template <class T>
const T
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
	return (N * sum_sq_prev_next - sum_prev * sum_next) / (sqrt(N * sum_sq_prev - sum_prev * sum_prev) * sqrt(N * sum_sq_next - sum_next * sum_next) + 1E-10);
}




template <class T>
const T
BlockMatching<T>::MAD(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height, const ImgVector<bool>& mask)
{
	double N = .0;
	T sad = 0;

	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			if (mask.get(x, y)) {
				N += 1.0;
				sad = sad + fabs((double)_image_next.get_zeropad(x_next + x, y_next + y) - _image_prev.get_zeropad(x_prev + x, y_prev + y));
			}
		}
	}
	return sad / (double)N;
}

template <class T>
const T
BlockMatching<T>::ZNCC(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height, const ImgVector<bool>& mask)
{
	double N = block_width * block_height;
	T sum_prev = 0;
	T sum_next = 0;
	T sum_sq_prev = 0;
	T sum_sq_next = 0;
	T sum_sq_prev_next = 0;

	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			if (mask.get(x, y)) {
				sum_prev += _image_prev.get_zeropad(x_prev + x, y_prev + y);
				sum_next += _image_next.get_zeropad(x_next + x, y_next + y);
				sum_sq_prev += _image_prev.get_zeropad(x_prev + x, y_prev + y) * _image_prev.get_zeropad(x_prev + x, y_prev + y);
				sum_sq_next += _image_next.get_zeropad(x_next + x, y_next + y) * _image_next.get_zeropad(x_next + x, y_next + y);
				sum_sq_prev_next += _image_prev.get_zeropad(x_prev + x, y_prev + y) * _image_next.get_zeropad(x_next + x, y_next + y);
			}
		}
	}
	return (N * sum_sq_prev_next - sum_prev * sum_next) / sqrt(N * sum_sq_prev - sum_prev * sum_prev) / sqrt(N * sum_sq_next - sum_next * sum_next);
}

