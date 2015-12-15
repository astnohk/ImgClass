#include <algorithm>
#include <cassert>
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
	_width = 0;
	_height = 0;
	_kernel_spatial = 10.0;
	_kernel_intensity = 0.1;
}

template <class T>
Segmentation<T>::Segmentation(const ImgVector<T>& image, const double kernel_spatial_radius, const double kernel_intensity_radius)
{
	_image.copy(image);
	if (_image.max() > 1.0) {
		// Normalize
		_image.contrast_stretching(0.0, 1.0);
	}
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

	_decrease_color_image.reset(_width, _height);
	_shift_vector.reset(_width, _height);
	_segments_map.reset(_width, _height);

	// Initial Segmentation
	Segmentation_MeanShift();
}

template <class T>
Segmentation<T>::Segmentation(const Segmentation<T>& segments) // Copy constructor
{
	_width = segments._width;
	_height = segments._height;
	_kernel_spatial = segments._kernel_spatial;
	_kernel_intensity = segments._kernel_intensity;

	_image.copy(segments._image);
	_decrease_color_image.copy(segments._decrease_color_image);
	_shift_vector.copy(segments._shift_vector);
	_segments_map.copy(segments._segments_map);
	_regions.assign(segments._regions.begin(), segments._regions.end());
}


template <class T>
Segmentation<T> &
Segmentation<T>::reset(const ImgVector<T>& image, const double kernel_spatial_radius, const double kernel_intensity_radius)
{
	_image.copy(image);
	if (_image.max() > 1.0) {
		// Normalize
		_image.contrast_stretching(0.0, 1.0);
	}
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

	_decrease_color_image.reset(_width, _height);
	_segments_map.reset(_width, _height);
	_regions.clear();
	_shift_vector.reset(_width, _height);

	// Initial Segmentation
	Segmentation_MeanShift();
	return *this;
}

// For L*a*b* color
template <>
Segmentation<ImgClass::RGB> &
Segmentation<ImgClass::RGB>::reset(const ImgVector<ImgClass::RGB>& image, const double kernel_spatial_radius, const double kernel_intensity_radius)
{
	_image.copy(image);
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

	_decrease_color_image.reset(_width, _height);
	_segments_map.reset(_width, _height);
	_regions.clear();
	_shift_vector.reset(_width, _height);

	// Initial Segmentation
	Segmentation_MeanShift();
	return *this;
}

// For L*a*b* color
template <>
Segmentation<ImgClass::Lab> &
Segmentation<ImgClass::Lab>::reset(const ImgVector<ImgClass::Lab>& image, const double kernel_spatial_radius, const double kernel_intensity_radius)
{
	_image.copy(image);
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

	_decrease_color_image.reset(_width, _height);
	_segments_map.reset(_width, _height);
	_regions.clear();
	_shift_vector.reset(_width, _height);

	// Initial Segmentation
	Segmentation_MeanShift();
	return *this;
}


template <class T>
Segmentation<T> &
Segmentation<T>::copy(const Segmentation<T>& segments)
{
	_width = segments._width;
	_height = segments._height;
	_kernel_spatial = segments._kernel_spatial;
	_kernel_intensity = segments._kernel_intensity;

	_image.copy(segments._image);
	_decrease_color_image.copy(segments._decrease_color_image);
	_shift_vector.copy(segments._shift_vector);
	_segments_map.copy(segments._segments_map);
	_regions.assign(segments._regions.begin(), segments._regions.end());
	return *this;
}


// ----- Destructor -----
template <class T>
Segmentation<T>::~Segmentation(void)
{
}




// ----- Data -----
template <class T>
const ImgVector<int> &
Segmentation<T>::ref_decrease_color_image(void) const
{
	return _decrease_color_image;
}

template <class T>
const ImgVector<int> &
Segmentation<T>::ref_segments_map(void) const
{
	return _segments_map;
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
int & 
Segmentation<T>::operator[](int n)
{
	return _segments_map[n];
}

template <class T>
int & 
Segmentation<T>::at(int n)
{
	assert(0 <= n && n < _width * _height);
	return _segments_map[n];
}

template <class T>
int & 
Segmentation<T>::at(int x, int y)
{
	return _segments_map.at(x, y);
}

template <class T>
int & 
Segmentation<T>::at_repeat(int x, int y)
{
	return _segments_map.at_repeat(x, y);
}

template <class T>
int & 
Segmentation<T>::at_mirror(int x, int y)
{
	return _segments_map.at_mirror(x, y);
}


template <class T>
int
Segmentation<T>::get(int n) const
{
	return _segments_map.get(n);
}

template <class T>
int
Segmentation<T>::get(int x, int y) const
{
	return _segments_map.get(x, y);
}

template <class T>
int
Segmentation<T>::get_zeropad(int x, int y) const
{
	return _segments_map.get_zeropad(x, y);
}

template <class T>
int
Segmentation<T>::get_repeat(int x, int y) const
{
	return _segments_map.get_repeat(x, y);
}

template <class T>
int
Segmentation<T>::get_mirror(int x, int y) const
{
	return _segments_map.get_mirror(x, y);
}




// ----- Mean Shift -----

template <class T>
void
Segmentation<T>::Segmentation_MeanShift(const int Iter_Max, const unsigned int Min_Number_of_Pixels)
{
	const double Decreased_Gray_Max = 255;
	const VECTOR_2D<int> adjacent[4] = {(VECTOR_2D<int>){1, 0}, (VECTOR_2D<int>){0, 1}, (VECTOR_2D<int>){-1, 0}, (VECTOR_2D<int>){0, -1}};
	ImgVector<int> vector_converge_map(_width, _height);
	std::vector<VECTOR_2D<int> > pel_list((2 * _kernel_spatial + 1) * (2 * _kernel_spatial + 1));
	std::list<std::list<VECTOR_2D<int> > > regions_list(0);
	std::vector<std::list<VECTOR_2D<int> > > regions_vector(0);
	unsigned int num_region = 0;
	unsigned int num_small_region = 0;

	if (_width <= 0 || _height <= 0) {
		return;
	}
	// Make pixel list
	unsigned int num = 0;
	for (int m = -_kernel_spatial; m <= _kernel_spatial; m++) {
		for (int n = -_kernel_spatial; n <= _kernel_spatial; n++) {
			if (n * n + m * m <= SQUARE(_kernel_spatial)) {
				pel_list[num] = (VECTOR_2D<int>){n, m};
				num++;
			}
		}
	}
	pel_list.resize(num);
	// Compute Mean Shift vector
#pragma omp parallel for schedule(dynamic)
	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			_shift_vector.at(x, y) = MeanShift(x, y, pel_list, Iter_Max);
			VECTOR_2D<int> r((int)round(_shift_vector.get(x, y).x), (int)round(_shift_vector.get(x, y).y));
			if (r.x < 0) {
				r.x = 0;
			} else if (_width <= r.x) {
				r.x = _width - 1;
			}
			if (r.y < 0) {
				r.y = 0;
			} else if (_height <= r.y) {
				r.y = _height - 1;
			}
			vector_converge_map.at(r.x, r.y) = -1;
		}
	}
	// Collect connected region
	num_region = 1; // To consider the out of image as region #0
	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			if (vector_converge_map.get(x, y) < 0) {
				VECTOR_2D<int> r(x, y);
				regions_list.push_back(std::list<VECTOR_2D<int> >(1, r));
				for (std::list<VECTOR_2D<int> >::iterator ite = regions_list.back().begin();
				    ite != regions_list.back().end();
				    ++ite) {
					for (int k = 0; k < 4; k++) { // Search 4-adjacent
						r = (VECTOR_2D<int>){ite->x + adjacent[k].x, ite->y + adjacent[k].y};
						if (vector_converge_map.get_zeropad(r.x, r.y) < 0) {
							vector_converge_map.at(r.x, r.y) = 0; // eliminate collected point from map
							regions_list.back().push_back(r);
						}
					}
				}
				for (std::list<VECTOR_2D<int> >::iterator ite = regions_list.back().begin();
				    ite != regions_list.back().end();
				    ++ite) {
					vector_converge_map.at(ite->x, ite->y) = num_region;
				}
				num_region++;
			}
		}
	}
	// Make the list of regions
	regions_list.clear();
	regions_vector.resize(num_region);
	for (int y = 0; y < _shift_vector.height(); y++) {
		for (int x = 0; x < _shift_vector.width(); x++) {
			VECTOR_2D<int> r((int)round(_shift_vector.get(x, y).x), (int)round(_shift_vector.get(x, y).y));
			int n_region = vector_converge_map.get_zeropad(r.x, r.y);
			regions_vector[n_region].push_back(VECTOR_2D<int>(x, y));
		}
	}
	// Set _segments_map by _regions No.
	for (unsigned int n = 0; n < regions_vector.size(); n++) {
		for (std::list<VECTOR_2D<int> >::iterator ite = regions_vector[n].begin();
		    ite != regions_vector[n].end();
		    ++ite) {
			_segments_map.at(ite->x, ite->y) = n;
		}
	}
	// Eliminate small connected regions
	num_small_region = small_region_eliminate(&regions_vector, Min_Number_of_Pixels);
	// Make decreased color image
	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			double mean = mean_kernel(
			    (int)round(_shift_vector.at(x, y).x),
			    (int)round(_shift_vector.at(x, y).y),
			    pel_list);
			_decrease_color_image.at(x, y) = 1 + round(Decreased_Gray_Max * mean);
		}
	}
	// ----- Output -----
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
	// Set _segments_map by _regions No.
	for (unsigned int n = 0; n < _regions.size(); n++) {
		for (unsigned int i = 0; i < _regions[n].size(); i++) {
			_segments_map.at(_regions[n][i].x, _regions[n][i].y) = n;
		}
	}
}

template <class T>
unsigned int
Segmentation<T>::small_region_eliminate(std::vector<std::list<VECTOR_2D<int> > >* regions_vector, const unsigned int Min_Number_of_Pixels)
{
	const VECTOR_2D<int> adjacent[4] = {(VECTOR_2D<int>){1, 0}, (VECTOR_2D<int>){0, 1}, (VECTOR_2D<int>){-1, 0}, (VECTOR_2D<int>){0, -1}};
	std::vector<bool> small_regions(regions_vector->size(), false);
	int num_small_region = 0;

	// Check small regions
	for (unsigned int n = 0; n < regions_vector->size(); n++) {
		if (regions_vector->at(n).size() < Min_Number_of_Pixels) {
			for (std::list<VECTOR_2D<int> >::iterator ite = regions_vector->at(n).begin();
			    ite != regions_vector->at(n).end();
			    ++ite) {
				small_regions[n] = true;
			}
		}
	}
	// Concatenate small regions
	for (int n = 0; n < small_regions.size(); n++) {
		if (small_regions[n] == false) {
			continue;
		}
		T center_color = _image.get(regions_vector->at(n).begin()->x, regions_vector->at(n).begin()->y);
		int concatenate_target = n;
		double min = 1000.0;
		bool check = false;
		for (std::list<VECTOR_2D<int> >::iterator ite = regions_vector->at(n).begin();
		    ite != regions_vector->at(n).end();
		    ++ite) {
			for (int k = 0; k < 4; k++) {
				VECTOR_2D<int> r(ite->x + adjacent[k].x, ite->y + adjacent[k].y);
				double dist = 0.0;
				if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height
				    && _segments_map.get(r.x, r.y) != n
				    && (dist = this->distance(center_color, _image.get(r.x, r.y))) < min) {
					check = true;
					min = dist;
					concatenate_target = _segments_map.get(r.x, r.y);
				}
			}
		}
		if (check) {
			regions_vector->at(concatenate_target).splice(regions_vector->at(concatenate_target).end(), regions_vector->at(n));
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
	VECTOR_2D<double> u;
	double intensity;

	// Initialize
	u.x = x;
	u.y = y;
	intensity = _image.get(x, y);
	// Iterate until it converge
	for (int i = 0; i < Iter_Max; i++) {
		double N = 0.0;
		double sum_intensity_diff = 0.0;
		VECTOR_2D<double> sum_d(0.0, 0.0);
		VECTOR_2D<double> d_tmp;

		for (unsigned int n = 0; n < pel_list.size(); n++) {
			VECTOR_2D<double> r(u.x + pel_list[n].x, u.y + pel_list[n].y);

			if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height) {
				double intensity_diff = _image.get(r.x, r.y) - intensity;

				if (fabs(intensity_diff) <= _kernel_intensity) {
					double coeff = 1.0 - (
					    SQUARE(intensity_diff / _kernel_intensity)
					    * ((SQUARE(pel_list[n].x) + SQUARE(pel_list[n].y)) / SQUARE(_kernel_spatial))
					    );
					N += coeff;
					sum_intensity_diff += intensity_diff * coeff;
					sum_d.x += pel_list[n].x * coeff;
					sum_d.y += pel_list[n].y * coeff;
				}
			}
		}
		intensity += sum_intensity_diff / N;
		d_tmp.x = std::min(sum_d.x / N, 1.0);
		d_tmp.y = std::min(sum_d.y / N, 1.0);
		if (Vector_2D::norm(d_tmp) < 0.01) {
			u += d_tmp;
			break;
		}
		u += d_tmp;
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
	// Initialize
	VECTOR_2D<double> u(x, y);
	ImgClass::RGB center(_image.get(x, y));
	// Iterate until it converge
	for (int i = 0; i < Iter_Max; i++) {
		double N = 0.0;
		ImgClass::RGB sum_diff(0.0, 0.0, 0.0);
		VECTOR_2D<double> sum_d(0.0, 0.0);
		VECTOR_2D<double> d_tmp;

		for (unsigned int n = 0; n < pel_list.size(); n++) {
			VECTOR_2D<double> r(u.x + pel_list[n].x, u.y + pel_list[n].y);

			if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height) {
				ImgClass::RGB diff(_image.get(r.x, r.y) - center);

				if (norm_squared(diff) <= _kernel_intensity * _kernel_intensity) {
					double coeff = 1.0 - (
					    (norm_squared(diff) / SQUARE(_kernel_intensity))
					    * ((SQUARE(pel_list[n].x) + SQUARE(pel_list[n].y)) / SQUARE(_kernel_spatial))
					    );
					N += coeff;
					sum_diff += diff * coeff;
					sum_d += VECTOR_2D<double>(pel_list[n]) * coeff;
				}
			}
		}
		center += sum_diff / N;
		d_tmp.x = std::min(sum_d.x / N, 1.0);
		d_tmp.y = std::min(sum_d.y / N, 1.0);
		if (Vector_2D::norm(d_tmp) < 0.01) {
			u += d_tmp;
			break;
		}
		u += d_tmp;
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
	// Initialize
	VECTOR_2D<double> u(x, y);
	ImgClass::Lab center = _image.get(x, y);
	// Iterate until it converge
	for (int i = 0; i < Iter_Max; i++) {
		double N = 0.0;
		ImgClass::Lab sum_diff(0.0, 0.0, 0.0);
		VECTOR_2D<double> sum_d(0.0, 0.0);
		VECTOR_2D<double> d_tmp;

		for (unsigned int n = 0; n < pel_list.size(); n++) {
			VECTOR_2D<double> r(u.x + pel_list[n].x, u.y + pel_list[n].y);

			if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height) {
				ImgClass::Lab diff(_image.get(r.x, r.y) - center);

				if (norm_squared(diff) <= _kernel_intensity * _kernel_intensity) {
					double coeff = 1.0 - (
					    (norm_squared(diff) / SQUARE(_kernel_intensity))
					    * ((SQUARE(pel_list[n].x) + SQUARE(pel_list[n].y)) / SQUARE(_kernel_spatial))
					    );
					N += coeff;
					sum_diff += diff * coeff;
					sum_d += VECTOR_2D<double>(pel_list[n]) * coeff;
				}
			}
		}
		center += sum_diff / N;
		d_tmp.x = std::min(sum_d.x / N, 1.0);
		d_tmp.y = std::min(sum_d.y / N, 1.0);
		if (Vector_2D::norm(d_tmp) < 0.01) {
			u += d_tmp;
			break;
		}
		u += d_tmp;
	}
	return u;
}


template <class T>
double
Segmentation<T>::mean_kernel(const int x, const int y, std::vector<VECTOR_2D<int> >& pel_list)
{
	double N = 0.0;
	double mean = 0.0;
	for (unsigned int n = 0; n < pel_list.size(); n++) {
		VECTOR_2D<double> r(x + pel_list[n].x, y + pel_list[n].y);

		if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height) {
			double coeff = 1.0
			    - (SQUARE(pel_list[n].x) + SQUARE(pel_list[n].y)) / SQUARE(_kernel_spatial);
			N += coeff;
			mean += double(_image.get(r.x, r.y)) * coeff;
		}
	}
	return mean / N;
}

