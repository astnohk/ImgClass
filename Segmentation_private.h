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
Segmentation<T>::ref_segments_map(void) const
{
	return _segments_map;
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
	const double Decreased_Gray_Max = 63; // 1 / 4 quantize
	const VECTOR_2D<int> adjacent[4] = {(VECTOR_2D<int>){1, 0}, (VECTOR_2D<int>){0, 1}, (VECTOR_2D<int>){-1, 0}, (VECTOR_2D<int>){0, -1}};
	std::vector<VECTOR_2D<int> > pel_list((2 * _kernel_spatial + 1) * (2 * _kernel_spatial + 1));
	std::list<std::list<VECTOR_2D<int> > > regions_list(0);
	std::vector<std::list<VECTOR_2D<int> > > regions_vector(0);
	ImgVector<int> decrease_color(_width, _height);
	unsigned int num_region;

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
	// Mean Shift method
#pragma omp parallel for schedule(dynamic)
	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			_shift_vector.at(x, y) = MeanShift(x, y, pel_list, Iter_Max);
			double mean = mean_kernel(
			    (int)round(_shift_vector.at(x, y).x),
			    (int)round(_shift_vector.at(x, y).y),
			    pel_list);
			decrease_color.at(x, y) = 1 + round(Decreased_Gray_Max * mean);
		}
	}
	// Collect connected region
	num_region = 0;
	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			if (decrease_color.get(x, y) > 0) {
				int color = decrease_color.get(x, y);
				VECTOR_2D<int> r(x, y);
				regions_list.push_back(std::list<VECTOR_2D<int> >(1, r));
				decrease_color.at(x, y) *= -1; // negative value indicates already checked
				for (std::list<VECTOR_2D<int> >::iterator ite = regions_list.back().begin();
				    ite != regions_list.back().end();
				    ++ite) {
					for (int k = 0; k < 4; k++) { // Search 4-adjacent
						if (decrease_color.get_zeropad(ite->x + adjacent[k].x, ite->y + adjacent[k].y) == color) {
							r = (VECTOR_2D<int>){ite->x + adjacent[k].x, ite->y + adjacent[k].y};
							decrease_color.at(r.x, r.y) *= -1; // negative value indicates already checked
							regions_list.back().push_back(r);
						}
					}
				}
				for (std::list<VECTOR_2D<int> >::iterator ite = regions_list.back().begin();
				    ite != regions_list.back().end();
				    ++ite) {
					_segments_map.at(ite->x, ite->y) = num_region;
				}
				num_region++;
			}
		}
	}
	regions_vector.assign(regions_list.begin(), regions_list.end());
	// Check small regions
	std::list<std::vector<int> > small_regions;
	for (unsigned int n = 0; n < regions_vector.size(); n++) {
		if (regions_vector[n].size() < Min_Number_of_Pixels) {
			for (std::list<VECTOR_2D<int> >::iterator ite = regions_vector[n].begin();
			    ite != regions_vector[n].end();
			    ++ite) {
				_segments_map.at(ite->x, ite->y) = -1;
				small_regions.push_back(std::vector<int>(3));
				small_regions.back()[0] = n;
				small_regions.back()[1] = ite->x;
				small_regions.back()[2] = ite->y;
			}
			num_region--;
			regions_vector[n].clear();
		}
	}
	// Eliminate small regions
	decrease_color *= -1; // negative all pixels to restore original intensity
	unsigned int small_regions_size_prev = 0;
	while (small_regions.size() > 0) { // Interpolate small region with of its neighborhood
		if (small_regions.size() == small_regions_size_prev) {
			break;
		}
		small_regions_size_prev = small_regions.size();
		for (std::list<std::vector<int> >::iterator ite = small_regions.begin();
		    ite != small_regions.end();
		    ) {
			int color = decrease_color.get(ite->at(1), ite->at(2));
			int min = Decreased_Gray_Max;
			bool check = false;
			for (int k = 0; k < 4; k++) {
				VECTOR_2D<int> r(ite->at(1) + adjacent[k].x, ite->at(2) + adjacent[k].y);
				if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height
				    && _segments_map.get_zeropad(r.x, r.y) >= 0
				    && _segments_map.get_zeropad(r.x, r.y) != ite->at(0)
				    && abs(color - decrease_color.get_zeropad(r.x, r.y)) < min) {
					check = true;
					min = fabs(color - decrease_color.get(r.x, r.y));
					decrease_color.at(ite->at(1), ite->at(2)) = decrease_color.get(r.x, r.y);
					_segments_map.at(ite->at(1), ite->at(2)) = _segments_map.at(r.x, r.y);
				}
			}
			if (check) {
				// remove from the list
				VECTOR_2D<int> v_tmp(ite->at(1), ite->at(2));
				regions_vector[_segments_map.at(ite->at(1), ite->at(2))].push_back(v_tmp);
				ite = small_regions.erase(ite);
			} else {
				// left current pixel and increment the iterator
				++ite;
			}
		}
	}
	if (small_regions.size() > 0) {
		std::cerr << "There are some pixels which are NOT able to be eliminated" << std::endl;
		for (std::list<std::vector<int> >::iterator ite = small_regions.begin();
		    ite != small_regions.end();
		    ++ite) {
			std::cerr << ite->at(0) << " " << ite->at(1) << " " << ite->at(2) << std::endl;
			std::cerr << "    " << decrease_color.get(ite->at(1), ite->at(2)) << std::endl;
		}
	}
	// Copy connected regions list
	_regions.resize(num_region);
	{
		std::vector<std::list<VECTOR_2D<int> > >::iterator ite = _regions.begin();
		for (unsigned int n = 0; n < regions_vector.size(); n++) {
			if (regions_vector[n].size() > 0) {
				ite->assign(regions_vector[n].begin(), regions_vector[n].end());
				++ite;
			}
		}
	}
	// Update _segments_map by _regions num
	for (unsigned int n = 0; n < _regions.size(); n++) {
		for (std::list<VECTOR_2D<int> >::iterator ite = _regions[n].begin();
		    ite != _regions[n].end();
		    ++ite) {
			_segments_map[ite->x, ite->y] = n;
		}
	}
	// Output vectors
	FILE *fp;
	FILE *fp_img;
	fp = fopen("segment_vector.dat", "w");
	fp_img = fopen("segments_color.pgm", "w");
	fprintf(fp, "%d %d\n", _width, _height);
	fprintf(fp_img, "P2\n%d %d\n%d\n", _width, _height, (int)Decreased_Gray_Max);
	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			_shift_vector.at(x, y).x -= x;
			_shift_vector.at(x, y).y -= y;
			fwrite(&(_shift_vector.at(x, y).x), sizeof(double), 1, fp);
			fwrite(&(_shift_vector.at(x, y).y), sizeof(double), 1, fp);
			fprintf(fp_img, "%d ", decrease_color.get(x, y));
		}
	}
	fclose(fp);
	fclose(fp_img);
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
			mean += _image.get(r.x, r.y) * coeff;
		}
	}
	return mean / N;
}


/*
 * std::vector<double> kernel has kernel radius for each dimensions.
 * The values it needs are below:
 *	_kernel_spatial : the spatial radius of mean shift kernel
 *	_kernel_intensity : the norm threshold of mean shift kernel in L*a*b* space
 */
template <> // Specialized for ImgClass::Lab
const VECTOR_2D<double>
Segmentation<ImgClass::Lab>::MeanShift_Lab(const int x, const int y, std::vector<VECTOR_2D<int> >& pel_list, int Iter_Max)
{
	// Initialize
	VECTOR_2D<double> u(x, y);
	ImgClass::Lab center = _image.get(x, y);
	// Iterate until it converge
	for (int i = 0; i < Iter_Max; i++) {
		double N = 0.0;
		ImgClass::Lab sum_diff(0.0);
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

