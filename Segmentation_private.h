#include <algorithm>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <list>

#include <cstdio>
#include <fstream>
#include <string>
#include <iostream>

#define SQUARE(a) ((a) * (a))



// ----- Constructor -----
template <class T>
Segmentation<T>::Segmentation(void)
{
	_size = 0;
	_width = 0;
	_height = 0;
	_kernel_spatial = 10.0;
	_kernel_intensity = 0.1;
}

template <class T>
Segmentation<T>::Segmentation(const ImgVector<T>& image, const double kernel_spatial_radius, const double kernel_intensity_radius)
{
	_image.copy(image);
	_size = _image.size();
	_width = _image.width();
	_height = _image.height();
	_kernel_spatial = kernel_spatial_radius;
	_kernel_intensity = kernel_intensity_radius;
	if (_kernel_spatial <= 0.0) {
		_kernel_spatial = 1.0;
	}
	if (_kernel_intensity <= 0.0) {
		_kernel_intensity = 1.0;
	}

	_color_quantized_image.reset(_width, _height);
	_shift_vector.reset(_width, _height);
	_segmentation_map.reset(_width, _height);

	// Initial Segmentation
	Segmentation_MeanShift();
}

template <class T>
Segmentation<T>::Segmentation(const Segmentation<T>& segmentation) // Copy constructor
{
	_size = segmentation._size;
	_width = segmentation._width;
	_height = segmentation._height;
	_kernel_spatial = segmentation._kernel_spatial;
	_kernel_intensity = segmentation._kernel_intensity;

	_image.copy(segmentation._image);
	_color_quantized_image.copy(segmentation._color_quantized_image);
	_shift_vector.copy(segmentation._shift_vector);
	_segmentation_map.copy(segmentation._segmentation_map);
	_regions.assign(segmentation._regions.begin(), segmentation._regions.end());
}


template <class T>
Segmentation<T> &
Segmentation<T>::reset(const ImgVector<T>& image, const double kernel_spatial_radius, const double kernel_intensity_radius)
{
	_image.copy(image);
	_size = _image.size();
	_width = _image.width();
	_height = _image.height();
	_kernel_spatial = kernel_spatial_radius;
	_kernel_intensity = kernel_intensity_radius;
	if (_kernel_spatial <= 0.0) {
		_kernel_spatial = 1.0;
	}
	if (_kernel_intensity <= 0.0) {
		_kernel_intensity = 1.0;
	}

	_color_quantized_image.reset(_width, _height);
	_segmentation_map.reset(_width, _height);
	_regions.clear();
	_shift_vector.reset(_width, _height);

	// Initial Segmentation
	Segmentation_MeanShift();
	return *this;
}


template <class T>
Segmentation<T> &
Segmentation<T>::copy(const Segmentation<T>& segmentation)
{
	_size = segmentation._size;
	_width = segmentation._width;
	_height = segmentation._height;
	_kernel_spatial = segmentation._kernel_spatial;
	_kernel_intensity = segmentation._kernel_intensity;

	_image.copy(segmentation._image);
	_color_quantized_image.copy(segmentation._color_quantized_image);
	_shift_vector.copy(segmentation._shift_vector);
	_segmentation_map.copy(segmentation._segmentation_map);
	_regions.assign(segmentation._regions.begin(), segmentation._regions.end());
	return *this;
}


// ----- Destructor -----
template <class T>
Segmentation<T>::~Segmentation(void)
{
}




// ----- Setter ------
template <class T>
void
Segmentation<T>::set_kernel(const double kernel_spatial_radius, const double kernel_intensity_radius)
{
	_kernel_spatial = kernel_spatial_radius;
	_kernel_intensity = kernel_intensity_radius;
}


template <class T>
Segmentation<T> &
Segmentation<T>::operator=(const Segmentation<T>& rvalue)
{
	_size = rvalue._size;
	_width = rvalue._width;
	_height = rvalue._height;
	_kernel_spatial = rvalue._kernel_spatial;
	_kernel_intensity = rvalue._kernel_intensity;

	_image.copy(rvalue._image);
	_color_quantized_image.copy(rvalue._color_quantized_image);
	_shift_vector.copy(rvalue._shift_vector);
	_segmentation_map.copy(rvalue._segmentation_map);
	_regions.assign(rvalue._regions.begin(), rvalue._regions.end());
	return *this;
}




// ----- Data -----
template <class T>
const ImgVector<T> &
Segmentation<T>::ref_color_quantized_image(void) const
{
	return _color_quantized_image;
}

template <class T>
const ImgVector<size_t> &
Segmentation<T>::ref_vector_converge_map(void) const
{
	return _vector_converge_map;
}

template <class T>
const ImgVector<size_t> &
Segmentation<T>::ref_segmentation_map(void) const
{
	return _segmentation_map;
}

template <class T>
const ImgVector<VECTOR_2D<double> > &
Segmentation<T>::ref_shift_vector(void) const
{
	return _shift_vector;
}

template <class T>
const std::vector<std::list<VECTOR_2D<int> > > &
Segmentation<T>::ref_regions(void) const
{
	return _regions;
}


// ----- Accessor -----
template <class T>
int
Segmentation<T>::width(void) const
{
	return _width;
}

template <class T>
int
Segmentation<T>::height(void) const
{
	return _height;
}

template <class T>
size_t
Segmentation<T>::size(void) const
{
	return _size;
}


template <class T>
size_t & 
Segmentation<T>::operator[](size_t n)
{
	return _segmentation_map[n];
}

template <class T>
size_t & 
Segmentation<T>::at(size_t n)
{
	assert(0 <= n && n < _size);
	return _segmentation_map[n];
}

template <class T>
size_t & 
Segmentation<T>::at(int x, int y)
{
	return _segmentation_map.at(x, y);
}

template <class T>
size_t & 
Segmentation<T>::at_repeat(int x, int y)
{
	return _segmentation_map.at_repeat(x, y);
}

template <class T>
size_t & 
Segmentation<T>::at_mirror(int x, int y)
{
	return _segmentation_map.at_mirror(x, y);
}


template <class T>
size_t
Segmentation<T>::get(size_t n) const
{
	return _segmentation_map.get(n);
}

template <class T>
size_t
Segmentation<T>::get(int x, int y) const
{
	return _segmentation_map.get(x, y);
}

template <class T>
size_t
Segmentation<T>::get_zeropad(int x, int y) const
{
	return _segmentation_map.get_zeropad(x, y);
}

template <class T>
size_t
Segmentation<T>::get_repeat(int x, int y) const
{
	return _segmentation_map.get_repeat(x, y);
}

template <class T>
size_t
Segmentation<T>::get_mirror(int x, int y) const
{
	return _segmentation_map.get_mirror(x, y);
}




// ----- Mean Shift -----
template <class T>
void
Segmentation<T>::Segmentation_MeanShift(const int Iter_Max, const unsigned int Min_Number_of_Pixels, const int Search_Range)
{
	const double Decreased_Gray_Max = 255.0;
	const VECTOR_2D<int> adjacent[8] = {
	    VECTOR_2D<int>(-1, -1), VECTOR_2D<int>(0, -1), VECTOR_2D<int>(1, -1),
	    VECTOR_2D<int>(-1, 0), VECTOR_2D<int>(1, 0),
	    VECTOR_2D<int>(-1, 1), VECTOR_2D<int>(0, 1), VECTOR_2D<int>(1, 1)};
	std::vector<VECTOR_2D<int> > pel_list;
	std::list<std::list<VECTOR_2D<int> > > regions_list(0);
	std::vector<std::list<VECTOR_2D<int> > > regions_vector(0);

	if (_width <= 0 || _height <= 0) {
		return;
	}
	// Make pixel list
	pel_list.reserve(SQUARE(static_cast<unsigned int>(ceil(2.0 * _kernel_spatial)) + 1));
	{
		for (int m = -int(ceil(_kernel_spatial)); m <= int(ceil(_kernel_spatial)); m++) {
			for (int n = -int(ceil(_kernel_spatial)); n <= int(ceil(_kernel_spatial)); n++) {
				if (n * n + m * m <= SQUARE(_kernel_spatial)) {
					pel_list.push_back(VECTOR_2D<int>(n, m));
				}
			}
		}
	}
	// Initialize vector converge map
	ImgVector<bool> vector_converge_mask(_width, _height, false);
	_vector_converge_map.reset(_width, _height);
	// Compute Mean Shift vector
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_SEGMENTATION)
	unsigned int finished = 0;
	unsigned int progress = .0;
	printf(" Mean-Shift method :   0.0%%\x1b[1A\n");
#endif
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			_shift_vector.at(x, y) = MeanShift(x, y, pel_list, Iter_Max);
			VECTOR_2D<int> r(
			    static_cast<int>(round(_shift_vector.get(x, y).x)),
			    static_cast<int>(round(_shift_vector.get(x, y).y)));
			// Check if the shift vector do NOT point on out of bounds
			if (r.x < 0) {
				_shift_vector.at(x, y).x = 0;
				r.x = 0;
			} else if (_width <= r.x) {
				_shift_vector.at(x, y).x = _width - 1;
				r.x = _width - 1;
			}
			if (r.y < 0) {
				_shift_vector.at(x, y).y = 0;
				r.y = 0;
			} else if (_height <= r.y) {
				_shift_vector.at(x, y).y = _height - 1;
				r.y = _height - 1;
			}
			// Check converge point
			vector_converge_mask.at(r.x, r.y) = true;
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_SEGMENTATION)
#ifdef _OPENMP
#pragma omp critical
#endif
		{
			double ratio = double(++finished) / _image.size();
			if (round(ratio * 1000.0) > progress) {
				progress = static_cast<unsigned int>(round(ratio * 1000.0)); // Take account of Over-Run
				printf("\r Mean-Shift method : %5.1f%%\x1b[1A\n", progress * 0.1);
			}
		}
#endif
		}
	}
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_SEGMENTATION)
	printf("\n Mean-Shift method : Finished\n");
#endif
	// Collect converged points
	{
		size_t num = 1;
		for (int y = 0; y < _height; y++) {
			for (int x = 0; x < _width; x++) {
				if (vector_converge_mask.get(x, y)) {
					VECTOR_2D<int> r(x, y);
					regions_list.push_back(std::list<VECTOR_2D<int> >(1, r));
					for (std::list<VECTOR_2D<int> >::iterator ite = regions_list.back().begin();
					    ite != regions_list.back().end();
					    ++ite) {
						for (int k = 0; k < 8; k++) { // Search 8-adjacent
							r.x = ite->x + adjacent[k].x;
							r.y = ite->y + adjacent[k].y;
							if (vector_converge_mask.get_zeropad(r.x, r.y)) {
								vector_converge_mask.at(r.x, r.y) = false; // eliminate collected point from map
								regions_list.back().push_back(r);
							}
						}
					}
					for (std::list<VECTOR_2D<int> >::iterator ite = regions_list.back().begin();
					    ite != regions_list.back().end();
					    ++ite) {
						_vector_converge_map.at(ite->x, ite->y) = num;
					}
					num++;
				}
			}
		}
	}
	regions_list.clear();
	// Pick all pixels related to limit points and Make the list of regions
	for (int y = 0; y < _shift_vector.height(); y++) {
		for (int x = 0; x < _shift_vector.width(); x++) {
			VECTOR_2D<int> r(int(round(_shift_vector.get(x, y).x)), int(round(_shift_vector.get(x, y).y)));
			size_t n_region = _vector_converge_map.get_zeropad(r.x, r.y);
			_segmentation_map.at(x, y) = n_region;
		}
	}
	// Collect connected regions from Mean-Shift filtered image
	unsigned int num_region = collect_regions_in_segmentation_map(&regions_vector);
	// Set _segmentation_map by _regions No.
	for (unsigned int n = 0; n < regions_vector.size(); n++) {
		for (std::list<VECTOR_2D<int> >::iterator ite = regions_vector[n].begin();
		    ite != regions_vector[n].end();
		    ++ite) {
			_segmentation_map.at(ite->x, ite->y) = n;
		}
	}
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_SEGMENTATION)
	printf(" Mean-Shift method : Eliminate small regions\n");
#endif
	// Eliminate small connected regions
	unsigned int num_small_region = small_region_eliminate(&regions_vector, Min_Number_of_Pixels, Search_Range);
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_SEGMENTATION)
	printf(" Mean-Shift method : the number of regions %u -> %u\n", num_region, num_region - num_small_region);
#endif
	// Copy connected regions list
	num_region -= num_small_region;
	_regions.resize(num_region);
	{
		std::vector<std::vector<VECTOR_2D<int> > >::iterator ite = _regions.begin();
		for (unsigned int n = 0; n < regions_vector.size(); n++) {
			if (regions_vector[n].size() > 0) {
				ite->assign(regions_vector[n].begin(), regions_vector[n].end());
				++ite;
			}
		}
	}
	// Reset _segmentation_map by _regions No.
	for (unsigned int n = 0; n < _regions.size(); n++) {
		for (unsigned int i = 0; i < _regions[n].size(); i++) {
			_segmentation_map.at(_regions[n][i].x, _regions[n][i].y) = n;
		}
	}
	// Make color-quantized image
	for (unsigned int n = 0; n < _regions.size(); n++) {
		T sum_color = T();
		for (unsigned int k = 0; k < _regions[n].size(); k++) {
			sum_color += _image.get(_regions[n][k].x, _regions[n][k].y);
		}
		T quantized_color = sum_color * Decreased_Gray_Max / double(_regions[n].size());
		for (unsigned int k = 0; k < _regions[n].size(); k++) {
			_color_quantized_image.at(_regions[n][k].x, _regions[n][k].y) += quantized_color;
		}
	}
}

template <class T>
unsigned int
Segmentation<T>::collect_regions_in_segmentation_map(std::vector<std::list<VECTOR_2D<int> > >* vector)
{
	const VECTOR_2D<int> adjacent[8] = {
	    VECTOR_2D<int>(-1, -1), VECTOR_2D<int>(0, -1), VECTOR_2D<int>(1, -1),
	    VECTOR_2D<int>(-1, 0), VECTOR_2D<int>(1, 0),
	    VECTOR_2D<int>(-1, 1), VECTOR_2D<int>(0, 1), VECTOR_2D<int>(1, 1)};
	std::list<std::list<VECTOR_2D<int> > > list(0);
	ImgVector<bool> collected(_width, _height, false);

	unsigned int num_region = 0;
	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			if (collected.get(x, y)) {
				continue;
			}
			VECTOR_2D<int> r(x, y);
			collected.at(x, y) = true;
			size_t N = _segmentation_map.get(r.x, r.y);
			list.push_back(std::list<VECTOR_2D<int> >(1, r));
			num_region++;
			// Search connected regions with 8-adjacent
			for (std::list<VECTOR_2D<int> >::iterator ite = list.back().begin();
			    ite != list.back().end();
			    ++ite) {
				for (int k = 0; k < 8; k++) {
					r.x = ite->x + adjacent[k].x;
					r.y = ite->y + adjacent[k].y;
					if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height
					    && collected.get(r.x, r.y) == false
					    && _segmentation_map.get(r.x, r.y) == N) {
						collected.at(r.x, r.y) = true;
						list.back().push_back(r);
					}
				}
			}
		}
	}
	vector->resize(num_region);
	{
		unsigned int N = 0;
		for (std::list<std::list<VECTOR_2D<int> > >::iterator ite = list.begin();
		    ite != list.end();
		    ++ite) {
			vector->at(N).assign(ite->begin(), ite->end());
			N++;
		}
	}
	return num_region;
}

template <class T>
unsigned int
Segmentation<T>::small_region_eliminate(std::vector<std::list<VECTOR_2D<int> > >* regions_vector, const unsigned int Min_Number_of_Pixels, const int search_range)
{
	std::vector<bool> small_regions(regions_vector->size(), false);
	unsigned int num_small_region = 0;

	// Check small regions
	for (unsigned int n = 0; n < regions_vector->size(); n++) {
		if (regions_vector->at(n).size() < Min_Number_of_Pixels
		    && regions_vector->at(n).size() > 0) {
			small_regions[n] = true;
		}
	}
	// Concatenate small regions
	for (size_t n = 0; n < regions_vector->size(); n++) {
		if (small_regions[n] == false) {
			continue;
		}
		num_small_region++; // Count the number of small regions here because some regions may become larger on this routine
		// Search nearest neighbor of the small region
		T center_color = _image.get(regions_vector->at(n).begin()->x, regions_vector->at(n).begin()->y);
		size_t concatenate_target = n;
		double min = DBL_MAX;
		bool check = false;
		for (std::list<VECTOR_2D<int> >::iterator ite = regions_vector->at(n).begin();
		    ite != regions_vector->at(n).end();
		    ++ite) {
			for (int y = -search_range; y <= search_range; y++) {
				for (int x = -search_range; x <= search_range; x++) {
					VECTOR_2D<int> r(ite->x + x, ite->y + y);
					if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height
					    && _segmentation_map.get(r.x, r.y) != n) {
						double diff = this->distance(center_color, _image.get(r.x, r.y));
						double dist = norm_squared(_shift_vector.get(r.x, r.y) - _shift_vector.get(ite->x, ite->y));
						if (diff + dist < min) {
							check = true;
							min = diff + dist;
							concatenate_target = _segmentation_map.get(r.x, r.y);
						}
					}
				}
			}
		}
		if (check) { // Concatenate small region to the neighborhood region
			// Update _segmentation_map
			VECTOR_2D<int> r = regions_vector->at(concatenate_target).front();
			for (std::list<VECTOR_2D<int> >::iterator ite = regions_vector->at(n).begin();
			    ite != regions_vector->at(n).end();
			    ++ite) {
				_segmentation_map.at(ite->x, ite->y) = concatenate_target;
				_color_quantized_image.at(ite->x, ite->y) = _color_quantized_image.get(r.x, r.y);
			}
			// splice the list
			regions_vector->at(concatenate_target).splice(regions_vector->at(concatenate_target).end(), regions_vector->at(n));
			// Uncheck the region which has became larger
			if (small_regions[concatenate_target]
			    && regions_vector->at(n).size() >= Min_Number_of_Pixels) {
				small_regions[concatenate_target] = false;
			}
		}
	}
	return num_small_region;
}


// Distance in any space
template <class T>
double
Segmentation<T>::distance(const T& lvalue, const T& rvalue)
{
	return fabs(lvalue - rvalue);
}

template <>
double
Segmentation<ImgClass::RGB>::distance(const ImgClass::RGB& lvalue, const ImgClass::RGB& rvalue)
{
	return sqrt(SQUARE(lvalue.R - rvalue.R)
	    + SQUARE(lvalue.G - rvalue.G)
	    + SQUARE(lvalue.B - rvalue.B));
}

template <>
double
Segmentation<ImgClass::Lab>::distance(const ImgClass::Lab& lvalue, const ImgClass::Lab& rvalue)
{
	return sqrt(SQUARE(lvalue.L - rvalue.L)
	    + SQUARE(lvalue.a - rvalue.a)
	    + SQUARE(lvalue.b - rvalue.b));
}


/*
 * std::vector<double> kernel has kernel radius for each dimensions.
 * The values it needs are below:
 *	_kernel_spatial : the spatial radius of mean shift kernel
 *	_kernel_intensity : the intensity threshold of mean shift kernel
 */
template <class T>
const VECTOR_2D<double>
Segmentation<T>::MeanShift(const int x, const int y, std::vector<VECTOR_2D<int> >& pel_list, int Iter_Max)
{
	const double radius_spatial_squared = SQUARE(_kernel_spatial);
	const double radius_intensity_squared = SQUARE(_kernel_intensity);
	const double displacement_min = SQUARE(0.01);
	// Initialize
	VECTOR_2D<double> u(static_cast<double>(x), static_cast<double>(y));
	double intensity = _image.get(x, y);
	// Iterate until it converge
	for (int i = 0; i < Iter_Max; i++) {
		double N = 0.0;
		double sum_intensity_diff = 0.0;
		VECTOR_2D<double> sum_d(0.0, 0.0);
		for (unsigned int n = 0; n < pel_list.size(); n++) {
			VECTOR_2D<int> r(
			    static_cast<int>(round(u.x) + pel_list[n].x),
			    static_cast<int>(round(u.y) + pel_list[n].y));
			if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height) {
				double intensity_diff = _image.get(r.x, r.y) - intensity;
				VECTOR_2D<double> d(r.x - u.x, r.y - u.y);
				double ratio_intensity = SQUARE(intensity_diff) / radius_intensity_squared;
				double ratio_spatial = norm_squared(d) / radius_spatial_squared;
				if (ratio_intensity <= 1.0 && ratio_spatial <= 1.0) {
					double coeff = 1.0 - (ratio_intensity * ratio_spatial);
					N += coeff;
					sum_intensity_diff += intensity_diff * coeff;
					sum_d += d * coeff;
				}
			}
		}
		intensity += sum_intensity_diff / N;
		VECTOR_2D<double> displacement(sum_d.x / N, sum_d.y / N);
		u += displacement;
		if (norm_squared(displacement) < displacement_min) {
			break;
		}
	}
	return u;
}

/*
 * std::vector<double> kernel has kernel radius for each dimensions.
 * The values it needs are below:
 *	_kernel_spatial : the spatial radius of mean shift kernel
 *	_kernel_intensity : the norm threshold of mean shift kernel in L*a*b* space
 */
template <> // Specialized for ImgClass::RGB<double>
const VECTOR_2D<double>
Segmentation<ImgClass::RGB>::MeanShift(const int x, const int y, std::vector<VECTOR_2D<int> >& pel_list, int Iter_Max)
{
	const double radius_spatial_squared = SQUARE(_kernel_spatial);
	const double radius_intensity_squared = SQUARE(_kernel_intensity);
	const double displacement_min = SQUARE(0.01);
	// Initialize
	VECTOR_2D<double> u(static_cast<double>(x), static_cast<double>(y));
	ImgClass::RGB center(_image.get(x, y));
	// Iterate until it converge
	for (int i = 0; i < Iter_Max; i++) {
		double N = 0.0;
		ImgClass::RGB sum_diff(0.0, 0.0, 0.0);
		VECTOR_2D<double> sum_d(0.0, 0.0);
		for (unsigned int n = 0; n < pel_list.size(); n++) {
			VECTOR_2D<int> r(
			    static_cast<int>(round(u.x) + pel_list[n].x),
			    static_cast<int>(round(u.y) + pel_list[n].y));
			if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height) {
				ImgClass::RGB diff(_image.get(r.x, r.y) - center);
				VECTOR_2D<double> d(r.x - u.x, r.y - u.y);
				double ratio_intensity = norm_squared(diff) / radius_intensity_squared;
				double ratio_spatial = norm_squared(d) / radius_spatial_squared;
				if (ratio_intensity <= 1.0 && ratio_spatial <= 1.0) {
					double coeff = 1.0 - (ratio_intensity * ratio_spatial);
					N += coeff;
					sum_diff += diff * coeff;
					sum_d += d * coeff;
				}
			}
		}
		center += sum_diff / N;
		VECTOR_2D<double> displacement(sum_d.x / N, sum_d.y / N);
		u += displacement;
		if (norm_squared(displacement) < displacement_min) {
			break;
		}
	}
	return u;
}

/*
 * std::vector<double> kernel has kernel radius for each dimensions.
 * The values it needs are below:
 *	_kernel_spatial : the spatial radius of mean shift kernel
 *	_kernel_intensity : the norm threshold of mean shift kernel in L*a*b* space
 */
template <> // Specialized for ImgClass::Lab
const VECTOR_2D<double>
Segmentation<ImgClass::Lab>::MeanShift(const int x, const int y, std::vector<VECTOR_2D<int> >& pel_list, int Iter_Max)
{
	const double radius_spatial_squared = SQUARE(_kernel_spatial);
	const double radius_intensity_squared = SQUARE(100.0 * _kernel_intensity);
	const double displacement_min = SQUARE(0.01);

	// Initialize
	VECTOR_2D<double> u(static_cast<double>(x), static_cast<double>(y));
	ImgClass::Lab center = _image.get(x, y);
	// Iterate until it converge
	for (int i = 0; i < Iter_Max; i++) {
		double N = 0.0;
		ImgClass::Lab sum_diff(0.0, 0.0, 0.0);
		VECTOR_2D<double> sum_d(0.0, 0.0);
		for (unsigned int n = 0; n < pel_list.size(); n++) {
			VECTOR_2D<int> r(
			    static_cast<int>(round(u.x) + pel_list[n].x),
			    static_cast<int>(round(u.y) + pel_list[n].y));
			if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height) {
				ImgClass::Lab diff(_image.get(r.x, r.y) - center);
				diff.L /= 4.0; // Difference of Lighting is not so important in segmentation
				VECTOR_2D<double> d(r.x - u.x, r.y - u.y);
				double ratio_intensity = norm_squared(diff) / radius_intensity_squared;
				double ratio_spatial = norm_squared(d) / radius_spatial_squared;
				if (ratio_intensity <= 1.0 && ratio_spatial <= 1.0) {
					double coeff = 1.0 - (ratio_intensity * ratio_spatial);
					N += coeff;
					sum_diff += diff * coeff;
					sum_d += d * coeff;
				}
			}
		}
		center += sum_diff / N;
		VECTOR_2D<double> displacement(sum_d.x / N, sum_d.y / N);
		u += displacement;
		if (norm_squared(displacement) < displacement_min) {
			break;
		}
	}
	return u;
}

