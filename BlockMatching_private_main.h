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
BlockMatching<T>::block_matching(const int search_range, const double coeff_MAD, const double coeff_ZNCC)
{
	if (_connected_regions_current.size() > 0) {
		block_matching_arbitrary_shaped(search_range, coeff_MAD, coeff_ZNCC);
	} else {
		block_matching_lattice(search_range, coeff_MAD, coeff_ZNCC);
	}
}




template <class T>
void
BlockMatching<T>::block_matching_lattice(const int search_range, const double coeff_MAD, const double coeff_ZNCC)
{
	double (BlockMatching<T>::*MAD_func)(const ImgVector<T>&, const ImgVector<T>&, const int, const int, const int, const int) = &BlockMatching<T>::MAD;
	double (BlockMatching<T>::*NCC_func)(const ImgVector<T>&, const ImgVector<T>&, const int, const int, const int, const int) = &BlockMatching<T>::ZNCC;

	if (this->isNULL()) {
		std::cerr << "void BlockMatching<T>::block_matching_lattice(const int) : _block_size < 0" << std::endl;
		throw std::logic_error("void BlockMatching<T>::block_matching(const int) : this is NULL");
	} else if (_block_size < 0) {
		std::cerr << "void BlockMatching<T>::block_matching_lattice(const int) : _block_size < 0" << std::endl;
		throw std::logic_error("void BlockMatching<T>::block_matching(const int) : _block_size < 0");
	}
	// Initialize
	// Reset Motion Vector
	_motion_vector_time.reset(_cells_width, _cells_height);
	_motion_vector_prev.reset(_cells_width, _cells_height);
	if (_image_next.isNULL() == false) {
		_motion_vector_next.reset(_cells_width, _cells_height);
	}
	// Set reference_images
	std::vector<ImgVector<T> *> reference_images;
	std::vector<ImgVector<VECTOR_2D<double> > *> motion_vectors;
	reference_images.push_back(&_image_prev);
	motion_vectors.push_back(&_motion_vector_prev);
	if (_image_next.isNULL() == false) {
		reference_images.push_back(&_image_next);
		motion_vectors.push_back(&_motion_vector_next);
	}
	// Compute Motion Vectors for previous and next frame
	for (size_t ref = 0; ref < reference_images.size(); ref++) {
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_BLOCKMATCHING)
		unsigned int finished = 0;
		unsigned int progress = .0;
		printf(" Block Matching :   0.0%%\x1b[1A\n");
#endif
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int Y_b = 0; Y_b < _cells_height; Y_b++) {
			int y_b = Y_b * _block_size;
			for (int X_b = 0; X_b < _cells_width; X_b++) {
				int x_b = X_b * _block_size;
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
				double E_min = DBL_MAX;
				VECTOR_2D<double> MV(.0, .0);
				for (int y = y_start; y <= y_end; y++) {
					for (int x = x_start; x <= x_end; x++) {
						VECTOR_2D<double> v_tmp(double(x - x_b), double(y - y_b));
						double MAD = (this->*MAD_func)(
						    *(reference_images[ref]), _image_current,
						    x, y, x_b, y_b);
						double ZNCC = (this->*NCC_func)(
						    *(reference_images[ref]), _image_current,
						    x, y, x_b, y_b);
						double E_tmp = coeff_MAD * MAD + coeff_ZNCC * (1.0 - ZNCC);
						if (E_tmp < E_min) {
							E_min = E_tmp;
							MV = v_tmp;
						} else if (fabs(E_tmp - E_min) < 1E-6
						    && norm_squared(MV) >= norm_squared(v_tmp)) {
							E_min = E_tmp;
							MV = v_tmp;
						}
					}
				}
				if (_subpixel_scale > 1) { // Sub-pixel scale search of infimum
					VECTOR_2D<double> MV_subpel(.0, .0);
					double MAD_min = DBL_MAX;
					for (int y = -_subpixel_scale + 1; y < _subpixel_scale; y++) {
						for (int x = -_subpixel_scale + 1; x < _subpixel_scale; x++) {
							VECTOR_2D<double> v_tmp(double(x) / double(_subpixel_scale), double(y) / double(_subpixel_scale));
							double MAD = MAD_cubic(
							    *(reference_images[ref]), _image_current,
							    double(x_b) + MV.x + v_tmp.x,
							    double(y_b) + MV.y + v_tmp.y,
							    double(x_b) + MV.x,
							    double(y_b) + MV.y);
							if (MAD < MAD_min) {
								MAD_min = MAD;
								MV_subpel = v_tmp;
							} else if (fabs(MAD - MAD_min) < 1E-6
							    && norm_squared(MV_subpel) >= norm_squared(v_tmp)) {
								MAD_min = MAD;
								MV_subpel = v_tmp;
							}
						}
					}
					std::cout << MV_subpel << std::endl;
					MV += MV_subpel;
				}
				motion_vectors[ref]->at(X_b, Y_b) = MV;
			}
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_BLOCKMATCHING)
			double ratio = double(++finished) / _cells_height;
			if (round(ratio * 1000.0) > progress) {
				progress = static_cast<unsigned int>(round(ratio * 1000.0)); // Take account of Over-Run
				printf("\r Block Matching : %5.1f%%\x1b[1A\n", progress * 0.1);
			}
#endif
		}
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_BLOCKMATCHING)
		printf("\n");
#endif
	}
	// Output
	if (_image_next.isNULL() == false) { // Use bi-directional motion estimation
		for (int Y_b = 0; Y_b < _cells_height; Y_b++) {
			for (int X_b = 0; X_b < _cells_width; X_b++) {
				VECTOR_2D<double> mv_p = _motion_vector_prev.get(X_b, Y_b);
				VECTOR_2D<double> mv_n = _motion_vector_next.get(X_b, Y_b);
				double diff_prev = 0.0;
				double diff_next = 0.0;
				for (int dy = 0; dy < _block_size; dy++) {
					int y = _block_size * Y_b + dy;
					for (int dx = 0; dx < _block_size; dx++) {
						int x = _block_size * X_b + dx;
						diff_prev += norm(_image_prev.get_zeropad_cubic(x + mv_p.x, y + mv_p.y) - _image_current.get(x, y));
						diff_next += norm(_image_next.get_zeropad_cubic(x + mv_n.x, y + mv_n.y) - _image_current.get(x, y));
					}
				}
				if (diff_prev <= diff_next) {
					_motion_vector_time.at(X_b, Y_b).x = _motion_vector_prev.at(X_b, Y_b).x;
					_motion_vector_time.at(X_b, Y_b).y = _motion_vector_prev.at(X_b, Y_b).y;
					_motion_vector_time.at(X_b, Y_b).t = -1;
				} else {
					_motion_vector_time.at(X_b, Y_b).x = _motion_vector_next.at(X_b, Y_b).x;
					_motion_vector_time.at(X_b, Y_b).y = _motion_vector_next.at(X_b, Y_b).y;
					_motion_vector_time.at(X_b, Y_b).t = 1;
				}
			}
		}
	} else { // Use forward motion estimation
		for (size_t n = 0; n < _motion_vector_prev.size(); n++) {
			_motion_vector_time[n].x = _motion_vector_prev[n].x;
			_motion_vector_time[n].y = _motion_vector_prev[n].y;
			_motion_vector_time[n].t = -1;
		}
	}
}


template <class T>
void
BlockMatching<T>::block_matching_arbitrary_shaped(const int search_range, const double coeff_MAD, const double coeff_ZNCC)
{
	double (BlockMatching<T>::*MAD_func)(const ImgVector<T>&, const ImgVector<T>&, const int, const int, const std::list<VECTOR_2D<int> >&) = &BlockMatching<T>::MAD_region;
	double (BlockMatching<T>::*NCC_func)(const ImgVector<T>&, const ImgVector<T>&, const int, const int, const std::list<VECTOR_2D<int> >&) = &BlockMatching<T>::ZNCC_region;
	//double (BlockMatching<T>::*MAD_func)(const ImgVector<T>&, const ImgVector<T>&, const int, const int, const std::list<VECTOR_2D<int> >&) = &BlockMatching<T>::MAD_region_nearest_intensity;
	//double (BlockMatching<T>::*NCC_func)(const ImgVector<T>&, const ImgVector<T>&, const int, const int, const std::list<VECTOR_2D<int> >&) = &BlockMatching<T>::ZNCC_region_nearest_intensity;

	if (this->isNULL()) {
		std::cerr << "void BlockMatching<T>::block_matching_region(const ImgVector<int>*, const int) : this is NULL" << std::endl;
		throw std::logic_error("void BlockMatching<T>::block_matching_region(const ImgVector<int>*, const int) : this is NULL");
	}
	// MV are expanded on entire image pixels
	_motion_vector_time.reset(_width, _height);
	_motion_vector_prev.reset(_width, _height);
	if (_image_next.isNULL() == false) {
		_motion_vector_next.reset(_width, _height);
	}
	// Set reference_images
	std::vector<ImgVector<T> *> reference_images;
	std::vector<ImgVector<VECTOR_2D<double> > *> motion_vectors;
	reference_images.push_back(&_image_prev);
	motion_vectors.push_back(&_motion_vector_prev);
	if (_image_next.isNULL() == false) {
		reference_images.push_back(&_image_next);
		motion_vectors.push_back(&_motion_vector_next);
	}
	// Compute Motion Vectors
	for (size_t ref = 0; ref < reference_images.size(); ref++) {
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_BLOCKMATCHING)
		unsigned int finished_regions = 0;
		unsigned int progress = .0;
		printf(" Block Matching :   0.0%%\x1b[1A\n");
#endif
		for (unsigned int n = 0; n < _connected_regions_current.size(); n++) {
			VECTOR_2D<double> MV(.0, .0);
			double E_min = DBL_MAX;
			// Search minimum value
			{
				int x;
#ifdef _OPENMP
#pragma omp parallel for
#endif
				for (x = 0; x < search_range * search_range; x++) {
					int y_diff = x / search_range - search_range / 2;
					int x_diff = x % search_range - search_range / 2;
					VECTOR_2D<double> v_tmp(x_diff, y_diff);
					double MAD = (this->*MAD_func)(
					    *(reference_images[ref]), _image_current,
					    x_diff, y_diff,
					    _connected_regions_current[n]);
					double ZNCC = (this->*NCC_func)(
					    *(reference_images[ref]), _image_current,
					    x_diff, y_diff,
					    _connected_regions_current[n]);
					double E_tmp = coeff_MAD * MAD + coeff_ZNCC * (1.0 - ZNCC);
#ifdef _OPENMP
#pragma omp critical
#endif
					{
						if (E_tmp < E_min) {
							E_min = E_tmp;
							MV = v_tmp;
						} else if (fabs(E_tmp - E_min) < DBL_EPSILON) {
							if (norm_squared(v_tmp) < norm_squared(MV)) {
								E_min = E_tmp;
								MV = v_tmp;
							}
						}
					}
				}
			}
			if (_subpixel_scale > 1) { // Sub-pixel scale search of infimum
				VECTOR_2D<double> MV_subpel(.0, .0);
				double MAD_min = DBL_MAX;
				for (int y = -_subpixel_scale + 1; y < _subpixel_scale; y++) {
					for (int x = -_subpixel_scale + 1; x < _subpixel_scale; x++) {
						VECTOR_2D<double> v_tmp(
						    double(x) / double(_subpixel_scale),
						    double(y) / double(_subpixel_scale));
						double MAD = MAD_region_cubic(
						    *(reference_images[ref]), _image_current,
						    MV.x + v_tmp.x, MV.y + v_tmp.y,
						    _connected_regions_current[n]);
						if (MAD < MAD_min) {
							MAD_min = MAD;
							MV_subpel = v_tmp;
						} else if (fabs(MAD - MAD_min) < 1E-6
						    && norm_squared(MV_subpel) >= norm_squared(v_tmp)) {
							MAD_min = MAD;
							MV_subpel = v_tmp;
						}
					}
				}
				MV += MV_subpel;
			}
			for (std::list<VECTOR_2D<int> >::iterator ite = _connected_regions_current[n].begin();
			    ite != _connected_regions_current[n].end();
			    ++ite) {
				motion_vectors[ref]->at(ite->x, ite->y) = MV;
			}
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_BLOCKMATCHING)
			double ratio = double(++finished_regions) / _connected_regions_current.size();
			if (round(ratio * 1000.0) > progress) {
				progress = static_cast<unsigned int>(round(ratio * 1000.0)); // Take account of Over-Run
				printf("\r Block Matching : %5.1f%%\x1b[1A\n", progress * 0.1);
			}
#endif
		}
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_BLOCKMATCHING)
		printf("\n");
#endif
	}
	for (unsigned int n = 0; n < _connected_regions_current.size(); n++) {
		for (std::list<VECTOR_2D<int> >::const_iterator ite = _connected_regions_current[n].begin();
		    ite != _connected_regions_current[n].end();
		    ++ite) {
			if (_image_next.isNULL()) { // Use forward motion estimation
				_motion_vector_time.at(ite->x, ite->y).x = _motion_vector_prev.get(ite->x, ite->y).x;
				_motion_vector_time.at(ite->x, ite->y).y = _motion_vector_prev.get(ite->x, ite->y).y;
				_motion_vector_time.at(ite->x, ite->y).t = -1;
			} else { // Use bi-directional motion estimation
				VECTOR_2D<double> mv_p = _motion_vector_prev.get(ite->x, ite->y);
				VECTOR_2D<double> mv_n = _motion_vector_next.get(ite->x, ite->y);
				double diff_prev = norm(_image_prev.get_zeropad_cubic(ite->x + mv_p.x, ite->y + mv_p.y) - _image_current.get(ite->x, ite->y));
				double diff_next = norm(_image_next.get_zeropad_cubic(ite->x + mv_n.x, ite->y + mv_n.y) - _image_current.get(ite->x, ite->y));
				if (diff_prev <= diff_next) { // Use forward motion vector
					_motion_vector_time.at(ite->x, ite->y).x = _motion_vector_prev.get(ite->x, ite->y).x;
					_motion_vector_time.at(ite->x, ite->y).y = _motion_vector_prev.get(ite->x, ite->y).y;
					_motion_vector_time.at(ite->x, ite->y).t = -1;
				} else {
					_motion_vector_time.at(ite->x, ite->y).x = _motion_vector_next.get(ite->x, ite->y).x;
					_motion_vector_time.at(ite->x, ite->y).y = _motion_vector_next.get(ite->x, ite->y).y;
					_motion_vector_time.at(ite->x, ite->y).t = 1;
				}
			}
		}
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
		    && (MAD_tmp = this->MAD(_image_prev, _image_prev, x_b, y_b, x_b, y_b - _block_size)) < MAD_min) {
			_motion_vector_prev.at(X_b, Y_b) = _motion_vector_prev.get_mirror(X_b, Y_b - 1);
			MAD_min = MAD_tmp;
		}
		if (estimated->get_mirror(X_b - 1, Y_b)
		    && (MAD_tmp = this->MAD(_image_prev, _image_prev, x_b, y_b, x_b - _block_size, y_b)) < MAD_min) {
			_motion_vector_prev.at(X_b, Y_b) = _motion_vector_prev.get_mirror(X_b - 1, Y_b);
			MAD_min = MAD_tmp;
		}
		if (estimated->get_mirror(X_b, Y_b + 1)
		    && (MAD_tmp = this->MAD(_image_prev, _image_prev, x_b, y_b, x_b, y_b + _block_size)) < MAD_min) {
			_motion_vector_prev.at(X_b, Y_b) = _motion_vector_prev.get_mirror(X_b, Y_b + 1);
			MAD_min = MAD_tmp;
		}
		if (estimated->get_mirror(X_b + 1, Y_b)
		    && (MAD_tmp = this->MAD(_image_prev, _image_prev, x_b, y_b, x_b + _block_size, y_b)) < MAD_min) {
			_motion_vector_prev.at(X_b, Y_b) = _motion_vector_prev.get_mirror(X_b + 1, Y_b);
			MAD_min = MAD_tmp;
		}
		estimated->at(X_b, Y_b) = true;
	}
}




// ----- Correlation -----
template <class T>
double
BlockMatching<T>::MAD(const ImgVector<T>& reference, const ImgVector<T>& interest, const int x_ref, const int y_ref, const int x_int, const int y_int)
{
	double sad = 0;

	for (int y = 0; y < _block_size; y++) {
		for (int x = 0; x < _block_size; x++) {
			sad += norm(
			    reference.get_zeropad(x_ref + x, y_ref + y)
			    - interest.get_zeropad(x_int + x, y_int + y));
		}
	}
	return sad / double(_block_size * _block_size);
}

template <class T>
double
BlockMatching<T>::MAD_cubic(const ImgVector<T>& reference, const ImgVector<T>& interest, const double x_ref, const double y_ref, const double x_int, const double y_int)
{
	double sad = 0;

	for (int y = 0; y < _block_size; y++) {
		for (int x = 0; x < _block_size; x++) {
			sad += norm(
			    reference.get_mirror_cubic(x_ref + x, y_ref + y)
			    - interest.get_mirror_cubic(x_int + x, y_int + y));
		}
	}
	return sad / double(_block_size * _block_size);
}




template <class T>
double
BlockMatching<T>::ZNCC(const ImgVector<T>& reference, const ImgVector<T>& interest, const int x_ref, const int y_ref, const int x_int, const int y_int)
{
	double N = _block_size * _block_size;
	T sum_reference = T();
	T sum_interest = T();
	double sum_sq_reference = 0;
	double sum_sq_interest = 0;
	double sum_sq_reference_interest = 0;

	for (int y = 0; y < _block_size; y++) {
		for (int x = 0; x < _block_size; x++) {
			sum_reference += reference.get_zeropad(x_ref + x, y_ref + y);
			sum_interest += interest.get_zeropad(x_int + x, y_int + y);
			sum_sq_reference += inner_prod(
			    reference.get_zeropad(x_ref + x, y_ref + y),
			    reference.get_zeropad(x_ref + x, y_ref + y));
			sum_sq_interest += inner_prod(
			    interest.get_zeropad(x_int + x, y_int + y),
			    interest.get_zeropad(x_int + x, y_int + y));
			sum_sq_reference_interest += inner_prod(
			    reference.get_zeropad(x_ref + x, y_ref + y),
			    interest.get_zeropad(x_int + x, y_int + y));
		}
	}
	return (N * sum_sq_reference_interest - inner_prod(sum_reference, sum_interest))
	    / (sqrt((N * sum_sq_reference - inner_prod(sum_reference, sum_reference))
	    * (N * sum_sq_interest - inner_prod(sum_interest, sum_interest)))
	    + DBL_EPSILON);
}




template <class T>
double
BlockMatching<T>::MAD_region(const ImgVector<T>& reference, const ImgVector<T>& interest, const int x_diff, const int y_diff, const std::list<VECTOR_2D<int> >& region_interest)
{
	double N = .0;
	double sad = .0;

	for (std::list<VECTOR_2D<int> >::const_iterator ite = region_interest.begin();
	    ite != region_interest.end();
	    ++ite) {
		N += 1.0;
		sad += norm(
		    interest.get_zeropad(ite->x, ite->y)
		    - reference.get_zeropad(ite->x + x_diff, ite->y + y_diff));
	}
	return sad / N;
}

template <class T>
double
BlockMatching<T>::MAD_region_cubic(const ImgVector<T>& reference, const ImgVector<T>& interest, const double x_diff, const double y_diff, const std::list<VECTOR_2D<int> >& region_interest)
{
	double N = .0;
	double sad = .0;

	for (std::list<VECTOR_2D<int> >::const_iterator ite = region_interest.begin();
	    ite != region_interest.end();
	    ++ite) {
		N += 1.0;
		sad += norm(
		    interest.get_zeropad(ite->x, ite->y)
		    - reference.get_mirror_cubic(double(ite->x) + x_diff, double(ite->y) + y_diff));
	}
	return sad / N;
}




template <class T>
double
BlockMatching<T>::ZNCC_region(const ImgVector<T>& reference, const ImgVector<T>& interest, const int x_diff, const int y_diff, const std::list<VECTOR_2D<int> >& region_interest)
{
	double N = .0;
	T sum_reference = T();
	T sum_interest = T();
	double sum_sq_reference = .0;
	double sum_sq_interest = .0;
	double sum_sq_reference_interest = .0;

	for (std::list<VECTOR_2D<int> >::const_iterator ite = region_interest.begin();
	    ite != region_interest.end();
	    ++ite) {
		VECTOR_2D<int> r(ite->x + x_diff, ite->y + y_diff);
		N += 1.0;
		// Previous frame
		sum_reference += reference.get_zeropad(r.x, r.y);
		sum_sq_reference += inner_prod(
		    reference.get_zeropad(r.x, r.y)
		    , reference.get_zeropad(r.x, r.y));
		// Next frame
		sum_interest += interest.get_zeropad(ite->x, ite->y);
		sum_sq_interest += inner_prod(
		    interest.get_zeropad(ite->x, ite->y)
		    , interest.get_zeropad(ite->x, ite->y));
		// Co-frame
		sum_sq_reference_interest += inner_prod(
		    reference.get_zeropad(r.x, r.y)
		    , interest.get_zeropad(ite->x, ite->y));
	}
	// Calculate Covariance
	return (N * sum_sq_reference_interest - inner_prod(sum_reference, sum_interest))
	    / (sqrt((N * sum_sq_reference - inner_prod(sum_reference, sum_reference))
	    * (N * sum_sq_interest - inner_prod(sum_interest, sum_interest)))
	    + DBL_EPSILON);
}




/* Regions concerning correlation method
 *
 * Compute correlation only with the regions which have near intensity.
 */
template <class T>
double
BlockMatching<T>::MAD_region_nearest_intensity(const int x_diff, const int y_diff, const std::list<VECTOR_2D<int> >& region_current)
{
	double N = .0;
	double sad = 0;

	for (std::list<VECTOR_2D<int> >::const_iterator ite = region_current.begin();
	    ite != region_current.end();
	    ++ite) {
		VECTOR_2D<int> r(ite->x + x_diff, ite->y + y_diff);

		double coeff = fabs(
		    _color_quantized_current.get(ite->x, ite->y)
		    - _color_quantized_prev.get_zeropad(r.x, r.y));
		N += coeff;
		sad += coeff * fabs(
		    double(_image_current.get(ite->x, ite->y)
		    - _image_prev.get_zeropad(r.x, r.y)));
	}
	return sad / N;
}

template <>
double
BlockMatching<ImgClass::RGB>::MAD_region_nearest_intensity(const int x_diff, const int y_diff, const std::list<VECTOR_2D<int> >& region_current);

template <>
double
BlockMatching<ImgClass::Lab>::MAD_region_nearest_intensity(const int x_diff, const int y_diff, const std::list<VECTOR_2D<int> >& region_current);




template <class T>
double
BlockMatching<T>::ZNCC_region_nearest_intensity(const int x_diff, const int y_diff, const std::list<VECTOR_2D<int> >& region_current)
{
	double N = 0;
	double sum_prev = 0;
	double sum_current = 0;
	double sum_sq_prev = 0;
	double sum_sq_current = 0;
	double sum_sq_prev_current = 0;

	for (std::list<VECTOR_2D<int> >::const_iterator ite = region_current.begin();
	    ite != region_current.end();
	    ++ite) {
		VECTOR_2D<int> r(ite->x + x_diff, ite->y + y_diff);

		double coeff = 1.0 - fabs(
		    _color_quantized_current.get(ite->x, ite->y)
		    - _color_quantized_prev.get_zeropad(r.x, r.y));
		N += coeff;
		// Previous frame
		sum_prev += coeff * _image_prev.get_zeropad(r.x, r.y);
		sum_sq_prev += coeff
		    * _image_prev.get_zeropad(r.x, r.y)
		    * _image_prev.get_zeropad(r.x, r.y);
		// Next frame
		sum_current += coeff * _image_current.get_zeropad(ite->x, ite->y);
		sum_sq_current += coeff
		    * _image_current.get_zeropad(ite->x, ite->y)
		    * _image_current.get_zeropad(ite->x, ite->y);
		// Co-frame
		sum_sq_prev_current += coeff
		    * _image_prev.get_zeropad(r.x, r.y)
		    * _image_current.get_zeropad(ite->x, ite->y);
	}
	// Calculate Covariance
	return (N * sum_sq_prev_current - sum_prev * sum_current) /
	    (sqrt((N * sum_sq_prev - sum_prev * sum_prev)
	    * (N * sum_sq_current - sum_current * sum_current))
	    + 1.0E-10);
}

template <>
double
BlockMatching<ImgClass::RGB>::ZNCC_region_nearest_intensity(const int x_diff, const int y_diff, const std::list<VECTOR_2D<int> >& region_current);

template <>
double
BlockMatching<ImgClass::Lab>::ZNCC_region_nearest_intensity(const int x_diff, const int y_diff, const std::list<VECTOR_2D<int> >& region_current);

