#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
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
	if (_connected_regions_next.size() > 0) {
		block_matching_arbitrary_shaped(search_range);
	} else if (_width == _cells_width && _height == _cells_height) {
		block_matching_dense_lattice(search_range);
	} else {
		block_matching_lattice(search_range);
	}
}




template <class T>
void
BlockMatching<T>::block_matching_lattice(const int search_range)
{
	const double coeff_SAD = 1.0;
	const double coeff_ZNCC = 10.0;
	double (BlockMatching<T>::*SAD_func)(const int, const int, const int, const int, const int, const int) = &BlockMatching<T>::MAD;
	double (BlockMatching<T>::*NCC_func)(const int, const int, const int, const int, const int, const int) = &BlockMatching<T>::ZNCC;

	ImgVector<bool> estimated(_cells_width, _cells_height, false);
	std::list<VECTOR_2D<int> > flat_blocks;
	std::list<VECTOR_2D<int> >::iterator itr_flat_blocks;

	if (this->isNULL()) {
		std::cerr << "void BlockMatching<T>::block_matching_lattice(const int) : _block_size < 0" << std::endl;
		throw std::logic_error("void BlockMatching<T>::block_matching(const int) : this is NULL");
	} else if (_block_size < 0) {
		std::cerr << "void BlockMatching<T>::block_matching_lattice(const int) : _block_size < 0" << std::endl;
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
	double grad_prev_max_min_x = grad_prev_x.max() - grad_prev_x.min();
	double grad_prev_max_min_y = grad_prev_y.max() - grad_prev_y.min();
	// Compute Motion Vectors
#ifdef _OPENMP
#pragma omp parallel for
#endif
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
			double SAD = (this->*SAD_func)(
			    x_start, y_start,
			    x_b, y_b,
			    _block_size, _block_size);
			double ZNCC = (this->*NCC_func)(
			    x_start, y_start,
			    x_b, y_b,
			    _block_size, _block_size);
			double E_min = coeff_SAD * SAD + coeff_ZNCC * (1.0 - ZNCC);
			for (int y = y_start; y <= y_end; y++) {
				for (int x = x_start; x <= x_end; x++) {
					VECTOR_2D<double> v_tmp(double(x - x_b), double(y - y_b));
					SAD = (this->*SAD_func)(
					    x, y,
					    x_b, y_b,
					    _block_size, _block_size);
					ZNCC = (this->*NCC_func)(
					    x, y,
					    x_b, y_b,
					    _block_size, _block_size);
					double E_tmp = coeff_SAD * SAD + coeff_ZNCC * (1.0 - ZNCC);
					if (E_tmp < E_min) {
						E_min = E_tmp;
						MV = v_tmp;
					} else if (fabs(E_tmp - E_min) < 1E-6) {
						if (norm(MV) >= norm(v_tmp)) {
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
BlockMatching<T>::block_matching_dense_lattice(const int search_range)
{
	const double coeff_SAD = 1.0;
	const double coeff_ZNCC = 8.0;
	double (BlockMatching<T>::*SAD_func)(const int, const int, const int, const int, const int, const int) = &BlockMatching<T>::MAD_centered;
	double (BlockMatching<T>::*NCC_func)(const int, const int, const int, const int, const int, const int) = &BlockMatching<T>::ZNCC;

	if (this->isNULL()) {
		std::cerr << "void BlockMatching<T>::block_matching_dense_lattice(const int) : _block_size < 0" << std::endl;
		throw std::logic_error("void BlockMatching<T>::block_matching_dense_lattice(const int) : this is NULL");
	} else if (_block_size < 0) {
		std::cerr << "void BlockMatching<T>::block_matching_dense_lattice(const int) : _block_size < 0" << std::endl;
		throw std::logic_error("void BlockMatching<T>::block_matching_dense_lattice(const int) : _block_size < 0");
	}
	// Initialize
	_motion_vector.reset(_width, _height);
	// Compute Motion Vectors
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int y_c = 0; y_c < _cells_height; y_c++) {
		int y_b = y_c - _block_size / 2;
		for (int x_c = 0; x_c < _cells_width; x_c++) {
			int x_b = x_c - _block_size / 2;

			VECTOR_2D<double> MV(.0, .0);
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
			double SAD = (this->*SAD_func)(
			    x_start, y_start,
			    x_b, y_b,
			    _block_size, _block_size);
			double ZNCC = (this->*NCC_func)(
			    x_start, y_start,
			    x_b, y_b,
			    _block_size, _block_size);
			double E_min = coeff_SAD * SAD + coeff_ZNCC * (1.0 - ZNCC);
			for (int y = y_start; y <= y_end; y++) {
				for (int x = x_start; x <= x_end; x++) {
					if (norm(_image_prev.get_zeropad(x, y) - _image_next.get(x_c, y_c)) < 0.07) {
						VECTOR_2D<double> v_tmp(double(x - x_c), double(y - y_c));
						SAD = (this->*SAD_func)(
						    x, y,
						    x_b, y_b,
						    _block_size, _block_size);
						ZNCC = (this->*NCC_func)(
						    x, y,
						    x_b, y_b,
						    _block_size, _block_size);
						double E_tmp = coeff_SAD * SAD + coeff_ZNCC * (1.0 - ZNCC);
						if (E_tmp < E_min) {
							E_min = E_tmp;
							MV = v_tmp;
						} else if (fabs(E_tmp - E_min) < 1E-6) {
							if (norm(MV) >= norm(v_tmp)) {
								E_min = E_tmp;
								MV = v_tmp;
							}
						}
					}
				}
			}
			_motion_vector.at(x_c, y_c) = MV;
		}
	}
}


template <class T>
void
BlockMatching<T>::block_matching_arbitrary_shaped(const int search_range)
{
	const int x_start = -search_range / 2;
	const int x_end = search_range / 2;
	const int y_start = -search_range / 2;
	const int y_end = search_range / 2;

	const double coeff_SAD = 1.0;
	const double coeff_ZNCC = 1.0;
	double (BlockMatching<T>::*SAD_func)(const int, const int, const std::list<VECTOR_2D<int> >&) = &BlockMatching<T>::MAD_region;
	double (BlockMatching<T>::*NCC_func)(const int, const int, const std::list<VECTOR_2D<int> >&) = &BlockMatching<T>::ZNCC_region;
	//double (BlockMatching<T>::*SAD_func)(const int, const int, const std::list<VECTOR_2D<int> >&) = &BlockMatching<T>::MAD_region_nearest_intensity;
	//double (BlockMatching<T>::*NCC_func)(const int, const int, const std::list<VECTOR_2D<int> >&) = &BlockMatching<T>::ZNCC_region_nearest_intensity;

	if (this->isNULL()) {
		std::cerr << "void BlockMatching<T>::block_matching_region(const ImgVector<int>*, const int) : this is NULL" << std::endl;
		throw std::logic_error("void BlockMatching<T>::block_matching_region(const ImgVector<int>*, const int) : this is NULL");
	}
	// MV are expanded on entire image pixels
	_motion_vector.reset(_width, _height);
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_BLOCKMATCHING)
	unsigned int finished_regions = 0;
	unsigned int progress = .0;
	printf(" Block Matching :   0.0%%\x1b[1A\n");
#endif
	// Compute Motion Vectors
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (unsigned int n = 0; n < _connected_regions_next.size(); n++) {
		VECTOR_2D<double> MV(.0, .0);

		// Compute initial value
		double SAD = (this->*SAD_func)(x_start, y_start, _connected_regions_next[n]);
		double ZNCC = (this->*NCC_func)(x_start, y_start, _connected_regions_next[n]);
		double E_min = coeff_SAD * SAD + coeff_ZNCC * (1.0 - ZNCC);
		// Search minimum value
		for (int y_diff = y_start; y_diff <= y_end; y_diff++) {
			for (int x_diff = x_start; x_diff <= x_end; x_diff++) {
				std::vector<VECTOR_2D<int> > region;
//				region = get_nearest_color_region(_connected_regions_next[n], x_diff, y_diff);

				VECTOR_2D<double> v_tmp(x_diff, y_diff);
				SAD = (this->*SAD_func)(
				    x_diff, y_diff,
				    _connected_regions_next[n]);
				ZNCC = (this->*NCC_func)(
				    x_diff, y_diff,
				    _connected_regions_next[n]);
				double E_tmp = coeff_SAD * SAD + coeff_ZNCC * (1.0 - ZNCC);
				if (E_tmp < E_min) {
					E_min = E_tmp;
					MV = v_tmp;
				} else if (fabs(E_tmp - E_min) < 1E-6) {
					if (norm(MV) >= norm(v_tmp)) {
						E_min = E_tmp;
						MV = v_tmp;
					}
				}
			}
		}
		for (std::list<VECTOR_2D<int> >::iterator ite = _connected_regions_next[n].begin();
		    ite != _connected_regions_next[n].end();
		    ++ite) {
			_motion_vector.at(ite->x, ite->y) = MV;
		}
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_BLOCKMATCHING)
#ifdef _OPENMP
#pragma omp critical
#endif
		{
			double ratio = double(++finished_regions) / _connected_regions_next.size();
			if (round(ratio * 1000.0) > progress) {
				progress = static_cast<unsigned int>(round(ratio * 1000.0)); // Take account of Over-Run
				printf("\r Block Matching : %5.1f%%\x1b[1A\n", progress * 0.1);
			}
		}
#endif
	}
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_BLOCKMATCHING)
	printf("\n Block Matching : Finished\n");
#endif
}




// ----- Interpolation -----
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
		double MAD_min = 10.0;
		double MAD_tmp;

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




// ----- Gradient -----
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

template <>
ImgVector<VECTOR_2D<double> > *
BlockMatching<ImgClass::Lab>::grad_prev(const int top_left_x, const int top_left_y, const int crop_width, const int crop_height);




// ----- Others -----
template <class T>
std::vector<VECTOR_2D<int> >
BlockMatching<T>::get_nearest_color_region(const std::list<VECTOR_2D<int> >& connected_region, const int x_diff, const int y_diff)
{
	std::vector<VECTOR_2D<int> > region;
	std::list<int> overlapping_regions;

	// Collect neighborhood regions
	for (std::list<VECTOR_2D<int> >::const_iterator ite_region = connected_region.begin();
	    ite_region != connected_region.end();
	    ++ite_region) {
		VECTOR_2D<int> r(ite_region->x + x_diff, ite_region->y + y_diff);
		if (r.x < 0 && _width <= r.x
		    && r.y < 0 && _height <= r.y) {
			continue;
		}
		int num = _region_map_prev.get(r.x, r.y);
		if (overlapping_regions.size() == 0) {
			overlapping_regions.push_back(num);
		} else {
			// Get overlapped region IDs
			std::list<int>::iterator ite = overlapping_regions.begin();
			for (;
			    ite != overlapping_regions.end();
			    ++ite) {
				if (*ite == num) {
					break;
				}
			}
			if (ite != overlapping_regions.end()) {
				overlapping_regions.push_back(num);
			}
		}
	}
	// Select neighborhood regions which has near colors to the region and covers connected_region greater than 9%
	T color_next = _color_quantized_next.get(connected_region.front().x, connected_region.front().y);
	double min = DBL_MAX;
	for (std::list<int>::const_iterator ite = overlapping_regions.begin();
	    ite != overlapping_regions.end();
	    ++ite) {
		VECTOR_2D<int> r = _connected_regions_prev[*ite].front();
		if (norm_squared(color_next - _color_quantized_prev.get(r.x, r.y)) < min) {
			min = norm_squared(color_next - _color_quantized_prev.get(r.x, r.y));
		}
	}
	return region;
}




// ----- Correlation -----
template <class T>
double
BlockMatching<T>::MAD(const int x_l, const int y_l, const int x_r, const int y_r, const int block_width, const int block_height, const ImgVector<T>& limage, const ImgVector<T>& rimage)
{
	double sad = 0;

	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			sad = sad + fabs(double(limage.get_zeropad(x_l+ x, y_l+ y) - rimage.get_zeropad(x_r + x, y_r + y)));
		}
	}
	return sad / double(block_width * block_height);
}

template <>
double
BlockMatching<ImgClass::RGB>::MAD(const int x_l, const int y_l, const int x_r, const int y_r, const int block_width, const int block_height, const ImgVector<ImgClass::RGB>& limage, const ImgVector<ImgClass::RGB>& rimage);

template <>
double
BlockMatching<ImgClass::Lab>::MAD(const int x_l, const int y_l, const int x_r, const int y_r, const int block_width, const int block_height, const ImgVector<ImgClass::Lab>& limage, const ImgVector<ImgClass::Lab>& rimage);




template <class T>
double
BlockMatching<T>::MAD_centered(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height)
{
	T sad = 0;
	T center_color;
	double N = 0.0;

	center_color = _image_next.get_zeropad(x_next + block_width / 2, y_next + block_height / 2);
	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			double coeff = 1.0 - fabs(double(center_color - _image_prev.get_zeropad(x_prev + x, y_prev + y)));
			N += coeff;
			sad += coeff * fabs(double(_image_next.get_zeropad(x_next + x, y_next + y) - _image_prev.get_zeropad(x_prev + x, y_prev + y)));
		}
	}
	return sad / N;
}

template <>
double
BlockMatching<ImgClass::RGB>::MAD_centered(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height);

template <>
double
BlockMatching<ImgClass::Lab>::MAD_centered(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height);




template <class T>
double
BlockMatching<T>::MAD(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height)
{
	double sad = 0;

	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			sad = sad + fabs(double(_image_next.get_zeropad(x_next + x, y_next + y) - _image_prev.get_zeropad(x_prev + x, y_prev + y)));
		}
	}
	return sad / double(block_width * block_height);
}

template <>
double
BlockMatching<ImgClass::RGB>::MAD(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height);

template <>
double
BlockMatching<ImgClass::Lab>::MAD(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height);




template <class T>
double
BlockMatching<T>::ZNCC(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height)
{
	double N = block_width * block_height;
	double sum_prev = 0;
	double sum_next = 0;
	double sum_sq_prev = 0;
	double sum_sq_next = 0;
	double sum_sq_prev_next = 0;

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

template <>
double
BlockMatching<ImgClass::RGB>::ZNCC(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height);

template <>
double
BlockMatching<ImgClass::Lab>::ZNCC(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height);




template <class T>
double
BlockMatching<T>::MAD_region(const int x_diff_prev, const int y_diff_prev, const std::list<VECTOR_2D<int> >& region)
{
	double N = .0;
	T sad = 0;

	for (std::list<VECTOR_2D<int> >::const_iterator ite = region.begin();
	    ite != region.end();
	    ++ite) {
		N += 1.0;
		sad = sad
		    + fabs(double(_image_next.get_zeropad(ite->x, ite->y)
		    - _image_prev.get_zeropad(ite->x + x_diff_prev, ite->y + y_diff_prev)));
	}
	return sad / N;
}

template <>
double
BlockMatching<ImgClass::RGB>::MAD_region(const int x_diff_prev, const int y_diff_prev, const std::list<VECTOR_2D<int> >& region);

template <>
double
BlockMatching<ImgClass::Lab>::MAD_region(const int x_diff_prev, const int y_diff_prev, const std::list<VECTOR_2D<int> >& region);




template <class T>
double
BlockMatching<T>::ZNCC_region(const int x_diff_prev, const int y_diff_prev, const std::list<VECTOR_2D<int> >& region)
{
	double N = 0;
	T sum_prev = 0;
	T sum_next = 0;
	T sum_sq_prev = 0;
	T sum_sq_next = 0;
	T sum_sq_prev_next = 0;

	for (std::list<VECTOR_2D<int> >::const_iterator ite = region.begin();
	    ite != region.end();
	    ++ite) {
		N += 1.0;
		// Previous frame
		sum_prev += _image_prev.get_zeropad(ite->x + x_diff_prev, ite->y + y_diff_prev);
		sum_sq_prev += _image_prev.get_zeropad(ite->x + x_diff_prev, ite->y + y_diff_prev)
		    * _image_prev.get_zeropad(ite->x + x_diff_prev, ite->y + y_diff_prev);
		// Next frame
		sum_next += _image_next.get_zeropad(ite->x, ite->y);
		sum_sq_next += _image_next.get_zeropad(ite->x, ite->y)
		    * _image_next.get_zeropad(ite->x, ite->y);
		// Co-frame
		sum_sq_prev_next += _image_prev.get_zeropad(ite->x + x_diff_prev, ite->y + y_diff_prev)
		    * _image_next.get_zeropad(ite->x, ite->y);
	}
	// Calculate Covariance
	return (N * sum_sq_prev_next - sum_prev * sum_next)
	    / (sqrt(N * sum_sq_prev - sum_prev * sum_prev) * sqrt(N * sum_sq_next - sum_next * sum_next));
}

template <>
double
BlockMatching<ImgClass::RGB>::ZNCC_region(const int x_diff_prev, const int y_diff_prev, const std::list<VECTOR_2D<int> >& region);

template <>
double
BlockMatching<ImgClass::Lab>::ZNCC_region(const int x_diff_prev, const int y_diff_prev, const std::list<VECTOR_2D<int> >& region);




/* Regions concerning correlation method
 *
 * Compute correlation only with the regions which have near intensity.
 */
template <class T>
double
BlockMatching<T>::MAD_region_nearest_intensity(const int x_diff_prev, const int y_diff_prev, const std::list<VECTOR_2D<int> >& region)
{
	double N = .0;
	double sad = 0;

	for (std::list<VECTOR_2D<int> >::const_iterator ite = region.begin();
	    ite != region.end();
	    ++ite) {
		VECTOR_2D<int> r(ite->x + x_diff_prev, ite->y + y_diff_prev);
		double coeff = fabs(
		    _color_quantized_next.get(ite->x, ite->y)
		    - _color_quantized_prev.get_zeropad(r.x, r.y));
		N += coeff;
		sad += coeff * fabs(
		    double(_image_next.get(ite->x, ite->y)
		    - _image_prev.get_zeropad(r.x, r.y)));
	}
	return sad / N;
}

template <>
double
BlockMatching<ImgClass::RGB>::MAD_region_nearest_intensity(const int x_diff_prev, const int y_diff_prev, const std::list<VECTOR_2D<int> >& region);

template <>
double
BlockMatching<ImgClass::Lab>::MAD_region_nearest_intensity(const int x_diff_prev, const int y_diff_prev, const std::list<VECTOR_2D<int> >& region);




template <class T>
double
BlockMatching<T>::ZNCC_region_nearest_intensity(const int x_diff_prev, const int y_diff_prev, const std::list<VECTOR_2D<int> >& region)
{
	double N = 0;
	double sum_prev = 0;
	double sum_next = 0;
	double sum_sq_prev = 0;
	double sum_sq_next = 0;
	double sum_sq_prev_next = 0;

	for (std::list<VECTOR_2D<int> >::const_iterator ite = region.begin();
	    ite != region.end();
	    ++ite) {
		VECTOR_2D<int> r(ite->x + x_diff_prev, ite->y + y_diff_prev);
		double coeff = 1.0 - fabs(
		    _color_quantized_next.get(ite->x, ite->y)
		    - _color_quantized_prev.get_zeropad(r.x, r.y));
		N += coeff;
		// Previous frame
		sum_prev += coeff * _image_prev.get_zeropad(r.x, r.y);
		sum_sq_prev += coeff
		    * _image_prev.get_zeropad(r.x, r.y)
		    * _image_prev.get_zeropad(r.x, r.y);
		// Next frame
		sum_next += coeff * _image_next.get_zeropad(ite->x, ite->y);
		sum_sq_next += coeff
		    * _image_next.get_zeropad(ite->x, ite->y)
		    * _image_next.get_zeropad(ite->x, ite->y);
		// Co-frame
		sum_sq_prev_next += coeff
		    * _image_prev.get_zeropad(r.x, r.y)
		    * _image_next.get_zeropad(ite->x, ite->y);
	}
	// Calculate Covariance
	return (N * sum_sq_prev_next - sum_prev * sum_next)
	    / (sqrt(N * sum_sq_prev - sum_prev * sum_prev) * sqrt(N * sum_sq_next - sum_next * sum_next));
}

template <>
double
BlockMatching<ImgClass::RGB>::ZNCC_region_nearest_intensity(const int x_diff_prev, const int y_diff_prev, const std::list<VECTOR_2D<int> >& region);

template <>
double
BlockMatching<ImgClass::Lab>::ZNCC_region_nearest_intensity(const int x_diff_prev, const int y_diff_prev, const std::list<VECTOR_2D<int> >& region);

