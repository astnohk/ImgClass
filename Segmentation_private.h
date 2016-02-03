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
	_shift_vector_spatial.reset(_width, _height);
	_shift_vector_color.reset(_width, _height);
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
	_shift_vector_spatial.copy(segmentation._shift_vector_spatial);
	_shift_vector_color.copy(segmentation._shift_vector_color);
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
	_shift_vector_spatial.reset(_width, _height);
	_shift_vector_color.reset(_width, _height);

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
	_shift_vector_spatial.copy(segmentation._shift_vector_spatial);
	_shift_vector_color.copy(segmentation._shift_vector_color);
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
	_shift_vector_spatial.copy(rvalue._shift_vector_spatial);
	_shift_vector_color.copy(rvalue._shift_vector_color);
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
Segmentation<T>::ref_segmentation_map(void) const
{
	return _segmentation_map;
}

template <class T>
const ImgVector<VECTOR_2D<double> > &
Segmentation<T>::ref_shift_vector_spatial(void) const
{
	return _shift_vector_spatial;
}

template <class T>
const ImgVector<T> &
Segmentation<T>::ref_shift_vector_color(void) const
{
	return _shift_vector_color;
}

template <class T>
const std::vector<std::vector<VECTOR_2D<int> > > &
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
Segmentation<T>::Segmentation_MeanShift(const int Iter_Max, const size_t Min_Number_of_Pixels, const int Search_Range)
{
	const VECTOR_2D<int> adjacent[8] = {
	    VECTOR_2D<int>(-1, -1), VECTOR_2D<int>(0, -1), VECTOR_2D<int>(1, -1),
	    VECTOR_2D<int>(-1, 0), VECTOR_2D<int>(1, 0),
	    VECTOR_2D<int>(-1, 1), VECTOR_2D<int>(0, 1), VECTOR_2D<int>(1, 1)};
	const double Decreased_Gray_Max = 255.0;
	std::vector<VECTOR_2D<int> > pel_list;

	if (_width <= 0 || _height <= 0) {
		return;
	}
	// Make pixel list
	pel_list.reserve(SQUARE(size_t(ceil(2.0 * _kernel_spatial)) + 1));
	{
		for (int m = -int(ceil(_kernel_spatial)); m <= int(ceil(_kernel_spatial)); m++) {
			for (int n = -int(ceil(_kernel_spatial)); n <= int(ceil(_kernel_spatial)); n++) {
				if (n * n + m * m <= SQUARE(_kernel_spatial)) {
					pel_list.push_back(VECTOR_2D<int>(n, m));
				}
			}
		}
	}
	// Compute Mean Shift vector
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_SEGMENTATION)
	size_t finished = 0;
	size_t progress = .0;
	printf(" Mean-Shift method :   0.0%%\x1b[1A\n");

#endif
	// Compute Mean Shift
	ImgVector<std::list<VECTOR_2D<int> > > vector_converge_list_map(_width, _height);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			ImgClass::Segmentation::tuple<T> tmp = MeanShift(x, y, pel_list, Iter_Max);
			auto lambda = [](double x, double max) -> double {
				return x >= 0 ? x < max ? x : max - 1 : 0; };
			// Saturation
			tmp.spatial.x = lambda(tmp.spatial.x, _width);
			tmp.spatial.y = lambda(tmp.spatial.y, _height);
			// Set vector
			_shift_vector_spatial.at(x, y) = tmp.spatial;
			_shift_vector_color.at(x, y) = tmp.color;
			// Assign start point to converge list
			VECTOR_2D<int> shift(int(tmp.spatial.x), int(tmp.spatial.y));
#pragma omp critical
			{
				vector_converge_list_map.at(shift.x, shift.y).push_back(VECTOR_2D<int>(x, y));
			}
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_SEGMENTATION)
#ifdef _OPENMP
#pragma omp critical
#endif
			{
				double ratio = double(++finished) / _image.size();
				if (round(ratio * 1000.0) > progress) {
					progress = size_t(round(ratio * 1000.0)); // Take account of Over-Run
					printf("\r Mean-Shift method : %5.1f%%\x1b[1A\n", progress * 0.1);
				}
			}
#endif
		}
	}
	// Concatenate the list of connected regions
	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			if (vector_converge_list_map.at(x, y).size() > 0) {
				std::list<VECTOR_2D<int> > tmp_list;
				tmp_list.push_back(VECTOR_2D<int>(x, y));
				for (auto ite = tmp_list.begin(); ite != tmp_list.end(); ++ite) {
					for (size_t i = 0; i < 8; i++) {
						VECTOR_2D<int> r(ite->x + adjacent[i].x, ite->y + adjacent[i].y);
						if (0 <= r.x && r.x < _width
						    && 0 <= r.y && r.y < _height
						    && vector_converge_list_map.at(r.x, r.y).size() > 0) {
							tmp_list.push_back(r);
							vector_converge_list_map.at(x, y).splice(
							    vector_converge_list_map.at(x, y).end(),
							    vector_converge_list_map.at(r.x, r.y));
						}
					}
				}
			}
		}
	}
	std::list<std::list<VECTOR_2D<int> > > regions_list; // start_point, shift_vector_spatial, shift_vector_color
	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			// Search converge point
			if (vector_converge_list_map.at(x, y).size() > 0) {
				std::list<std::list<VECTOR_2D<int> > > tmp_list;
				for (VECTOR_2D<int>& start : vector_converge_list_map.at(x, y)) {
					bool found = false;
					const T& color = _shift_vector_color.get(start.x, start.y);
					for (std::list<VECTOR_2D<int> >& region : tmp_list) {
						VECTOR_2D<int> r_tmp = region.front();
						if (normalized_distance(color, _shift_vector_color.get(r_tmp.x, r_tmp.y)) < 0.2) {
							region.push_back(start);
							found = true;
							break;
						}
					}
					if (found == false) {
						tmp_list.push_back(std::list<VECTOR_2D<int> >(1, start));
					}
				}
				regions_list.splice(regions_list.end(), tmp_list);
			}
		}
	}
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_SEGMENTATION)
	printf("\n Mean-Shift method : Finished\n");
#endif
	// Make region map
	_segmentation_map.resize(_width, _height);
	{
		size_t num = 1;
		for (const std::list<VECTOR_2D<int> >& region : regions_list) {
			for (const VECTOR_2D<int>& r : region) {
				_segmentation_map.at(r.x, r.y) = num;
			}
			num++;
		}
	}
	regions_list.clear();
	// Collect connected regions from Mean-Shift filtered image
	std::vector<std::list<VECTOR_2D<int> > > regions_vector;
	size_t num_region = collect_regions_in_segmentation_map(&regions_vector);
	// Set _segmentation_map by _regions No.
	for (size_t n = 0; n < regions_vector.size(); n++) {
		for (const VECTOR_2D<int>& r : regions_vector[n]) {
			_segmentation_map.at(r.x, r.y) = n;
		}
	}
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_SEGMENTATION)
	printf(" Mean-Shift method : Eliminate small regions\n");
#endif
	// Eliminate small connected regions
	size_t num_small_region = small_region_eliminate(&regions_vector, Min_Number_of_Pixels, Search_Range);
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_SEGMENTATION)
	std::cout << " Mean-Shift method : The number of regions " << num_region << " -> " << num_region - num_small_region << "\n";
#endif
	// Copy connected regions list
	num_region -= num_small_region;
	_regions.resize(num_region);
	{
		std::vector<std::vector<VECTOR_2D<int> > >::iterator ite = _regions.begin();
		for (size_t n = 0; n < regions_vector.size(); n++) {
			if (regions_vector[n].size() > 0) {
				ite->assign(regions_vector[n].begin(), regions_vector[n].end());
				++ite;
			}
		}
	}
	// Reset _segmentation_map by _regions No.
	for (size_t n = 0; n < _regions.size(); n++) {
		for (size_t i = 0; i < _regions[n].size(); i++) {
			_segmentation_map.at(_regions[n][i].x, _regions[n][i].y) = n;
		}
	}
	// Make color-quantized image
	for (size_t n = 0; n < _regions.size(); n++) {
		T sum_color = T();
		for (size_t k = 0; k < _regions[n].size(); k++) {
			sum_color += _image.get(_regions[n][k].x, _regions[n][k].y);
		}
		T quantized_color = sum_color * Decreased_Gray_Max / double(_regions[n].size());
		for (size_t k = 0; k < _regions[n].size(); k++) {
			_color_quantized_image.at(_regions[n][k].x, _regions[n][k].y) += quantized_color;
		}
	}
}


template <class T>
size_t
Segmentation<T>::collect_regions_in_segmentation_map(std::vector<std::list<VECTOR_2D<int> > >* vector)
{
	const VECTOR_2D<int> adjacent[8] = {
	    VECTOR_2D<int>(-1, -1), VECTOR_2D<int>(0, -1), VECTOR_2D<int>(1, -1),
	    VECTOR_2D<int>(-1, 0), VECTOR_2D<int>(1, 0),
	    VECTOR_2D<int>(-1, 1), VECTOR_2D<int>(0, 1), VECTOR_2D<int>(1, 1)};
	std::list<std::list<VECTOR_2D<int> > > regions(0);
	ImgVector<bool> collected(_width, _height, false);

	size_t num_region = 0;
	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			if (collected.get(x, y)) {
				continue;
			}
			VECTOR_2D<int> r(x, y);
			collected.at(x, y) = true;
			size_t N = _segmentation_map.get(r.x, r.y);
			regions.push_back(std::list<VECTOR_2D<int> >(1, r));
			num_region++;
			// Search connected regions with 8-adjacent
			for (const VECTOR_2D<int>& element : regions.back()) {
				for (int k = 0; k < 8; k++) {
					r.x = element.x + adjacent[k].x;
					r.y = element.y + adjacent[k].y;
					if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height
					    && collected.get(r.x, r.y) == false
					    && _segmentation_map.get(r.x, r.y) == N) {
						collected.at(r.x, r.y) = true;
						regions.back().push_back(r);
					}
				}
			}
		}
	}
	vector->resize(num_region);
	{
		size_t N = 0;
		for (const std::list<VECTOR_2D<int> >& region : regions) {
			vector->at(N).assign(region.begin(), region.end());
			N++;
		}
	}
	return num_region;
}


template <class T>
size_t
Segmentation<T>::small_region_eliminate(std::vector<std::list<VECTOR_2D<int> > >* regions_vector, const size_t Min_Number_of_Pixels, const int search_range)
{
	std::vector<bool> small_regions(regions_vector->size(), false);
	size_t num_small_region = 0;

	// Check small regions
	for (size_t n = 0; n < regions_vector->size(); n++) {
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
		for (const VECTOR_2D<int>& element : regions_vector->at(n)) {
			for (int y = -search_range; y <= search_range; y++) {
				for (int x = -search_range; x <= search_range; x++) {
					VECTOR_2D<int> r(element.x + x, element.y + y);
					if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height
					    && _segmentation_map.get(r.x, r.y) != n) {
						double diff = this->distance(center_color, _image.get(r.x, r.y));
						double dist = norm_squared(_shift_vector_spatial.get(r.x, r.y) - _shift_vector_spatial.get(element.x, element.y));
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
			for (const VECTOR_2D<int>& element : regions_vector->at(n)) {
				_segmentation_map.at(element.x, element.y) = concatenate_target;
				_color_quantized_image.at(element.x, element.y) = _color_quantized_image.get(r.x, r.y);
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
Segmentation<ImgClass::RGB>::distance(const ImgClass::RGB& lvalue, const ImgClass::RGB& rvalue);

template <>
double
Segmentation<ImgClass::Lab>::distance(const ImgClass::Lab& lvalue, const ImgClass::Lab& rvalue);


/*
 * std::vector<double> kernel has kernel radius for each dimensions.
 * The values it needs are below:
 *	_kernel_spatial : the spatial radius of mean shift kernel
 *	_kernel_intensity : the intensity threshold of mean shift kernel
 */
template <class T>
const ImgClass::Segmentation::tuple<T>
Segmentation<T>::MeanShift(const int x, const int y, std::vector<VECTOR_2D<int> >& pel_list, int Iter_Max)
{
	const double radius_spatial_squared = SQUARE(_kernel_spatial);
	const double radius_intensity_squared = SQUARE(_kernel_intensity);
	const double displacement_min = SQUARE(0.01);
	// Initialize
	ImgClass::Segmentation::tuple<T> tuple;
	tuple.spatial = VECTOR_2D<double>(static_cast<double>(x), static_cast<double>(y));
	tuple.color = _image.get(x, y);
	// Iterate until it converge
	for (int i = 0; i < Iter_Max; i++) {
		double N = 0.0;
		double sum_intensity_diff = 0.0;
		VECTOR_2D<double> sum_d(0.0, 0.0);
		for (size_t n = 0; n < pel_list.size(); n++) {
			VECTOR_2D<int> r(
			    static_cast<int>(round(tuple.spatial.x) + pel_list[n].x),
			    static_cast<int>(round(tuple.spatial.y) + pel_list[n].y));
			if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height) {
				double intensity_diff = _image.get(r.x, r.y) - tuple.color;
				VECTOR_2D<double> d(r.x - tuple.spatial.x, r.y - tuple.spatial.y);
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
		tuple.color += sum_intensity_diff / N;
		VECTOR_2D<double> displacement(sum_d.x / N, sum_d.y / N);
		tuple.spatial += displacement;
		if (norm_squared(sum_intensity_diff / N) * norm_squared(displacement) < displacement_min) {
			break;
		}
	}
	return tuple;
}

/*
 * std::vector<double> kernel has kernel radius for each dimensions.
 * The values it needs are below:
 *	_kernel_spatial : the spatial radius of mean shift kernel
 *	_kernel_intensity : the norm threshold of mean shift kernel in L*a*b* space
 */
template <> // Specialized for ImgClass::RGB<double>
const ImgClass::Segmentation::tuple<ImgClass::RGB>
Segmentation<ImgClass::RGB>::MeanShift(const int x, const int y, std::vector<VECTOR_2D<int> >& pel_list, int Iter_Max);

/*
 * std::vector<double> kernel has kernel radius for each dimensions.
 * The values it needs are below:
 *	_kernel_spatial : the spatial radius of mean shift kernel
 *	_kernel_intensity : the norm threshold of mean shift kernel in L*a*b* space
 */
template <> // Specialized for ImgClass::Lab
const ImgClass::Segmentation::tuple<ImgClass::Lab>
Segmentation<ImgClass::Lab>::MeanShift(const int x, const int y, std::vector<VECTOR_2D<int> >& pel_list, int Iter_Max);

