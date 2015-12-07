#include <algorithm>
#include <cassert>
#include <cmath>
#include <list>

#include <cstdio>
#include <fstream>
#include <string>
#include <iostream>




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
	_segments.reset(_width, _height);

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
	_segments.copy(segments._segments);
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

	_segments.reset(_width, _height);
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
	_segments.copy(segments._segments);
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
Segmentation<T>::ref_segments(void) const
{
	return _segments;
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
	return _segments[n];
}

template <class T>
int & 
Segmentation<T>::at(int n)
{
	assert(0 <= n && n < _width * _height);
	return _segments[n];
}

template <class T>
int & 
Segmentation<T>::at(int x, int y)
{
	return _segments.at(x, y);
}

template <class T>
int & 
Segmentation<T>::at_repeat(int x, int y)
{
	return _segments.at_repeat(x, y);
}

template <class T>
int & 
Segmentation<T>::at_mirror(int x, int y)
{
	return _segments.at_mirror(x, y);
}


template <class T>
int
Segmentation<T>::get(int n) const
{
	return _segments.get(n);
}

template <class T>
int
Segmentation<T>::get(int x, int y) const
{
	return _segments.get(x, y);
}

template <class T>
int
Segmentation<T>::get_zeropad(int x, int y) const
{
	return _segments.get_zeropad(x, y);
}

template <class T>
int
Segmentation<T>::get_repeat(int x, int y) const
{
	return _segments.get_repeat(x, y);
}

template <class T>
int
Segmentation<T>::get_mirror(int x, int y) const
{
	return _segments.get_mirror(x, y);
}




// ----- Mean Shift -----

template <class T>
void
Segmentation<T>::Segmentation_MeanShift(const int Iter_Max, const unsigned int Min_Number_of_Pixels)
{
	const double Decreased_Gray_Max = 63; // 1 / 4 quantize
	const VECTOR_2D<int> adjacent[4] = {(VECTOR_2D<int>){1, 0}, (VECTOR_2D<int>){0, 1}, (VECTOR_2D<int>){-1, 0}, (VECTOR_2D<int>){0, -1}};
	std::vector<VECTOR_2D<int> > pel_list((2 * _kernel_spatial + 1) * (2 * _kernel_spatial + 1));
	std::list<std::list<VECTOR_2D<int> > > regions(0);
	ImgVector<int> decrease_color(_width, _height);

	if (_width <= 0 || _height <= 0) {
		return;
	}
	// Make pixel list
	int num = 0;
	for (int m = -_kernel_spatial; m <= _kernel_spatial; m++) {
		for (int n = -_kernel_spatial; n <= _kernel_spatial; n++) {
			if (n * n + m * m < _kernel_spatial * _kernel_spatial) {
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
			_shift_vector.at(x, y) = MeanShift_Grayscale(x, y, pel_list, Iter_Max);
			decrease_color.at(x, y) = 1 + int(Decreased_Gray_Max * _image.get_zeropad((int)round(_shift_vector.get(x, y).x), (int)round(_shift_vector.get(x, y).y)));
		}
	}
	// Collect connected region
	num = 1;
	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			if (decrease_color.get(x, y) > 0) {
				int color = decrease_color.get(x, y);
				VECTOR_2D<int> r(x, y);
				regions.push_back(std::list<VECTOR_2D<int> >(1, r));
				decrease_color.at(x, y) *= -1; // negative value indicates already checked
				for (std::list<VECTOR_2D<int> >::iterator ite = regions.back().begin();
				    ite != regions.back().end();
				    ++ite) {
					for (int k = 0; k < 4; k++) { // Search 4-adjacent
						if (decrease_color.get_zeropad(ite->x + adjacent[k].x, ite->y + adjacent[k].y) == color) {
							r = (VECTOR_2D<int>){ite->x + adjacent[k].x, ite->y + adjacent[k].y};
							decrease_color.at(r.x, r.y) *= -1; // negative value indicates already checked
							regions.back().push_back(r);
						}
					}
				}
				for (std::list<VECTOR_2D<int> >::iterator ite = regions.back().begin();
				    ite != regions.back().end();
				    ++ite) {
					_segments.at(ite->x, ite->y) = num;
				}
				num++;
			}
		}
	}
	// Check small regions
	std::list<VECTOR_2D<int> > small_region_coordinate;
	for (std::list<std::list<VECTOR_2D<int> > >::iterator ite_region = regions.begin();
	    ite_region != regions.end();
	    num++, ++ite_region) {
		if (ite_region->size() < Min_Number_of_Pixels) {
			for (std::list<VECTOR_2D<int> >::iterator ite = ite_region->begin();
			    ite != ite_region->end();
			    ++ite) {
				_segments.at(ite->x, ite->y) = 0;
				small_region_coordinate.push_back(*ite);
			}
		}
	}
	/*
	// Eliminate small regions
	decrease_color *= -1; // negative all pixels to restore original intensity
	num = 1;
	while (small_region_coordinate.size() > 0) { // Interpolate small region with of its neighborhood
		for (std::list<VECTOR_2D<int> >::iterator ite = small_region_coordinate.begin();
		    ite != small_region_coordinate.end();
		    num++, ++ite) {
			int color = decrease_color.get(ite->x, ite->y);
			int min = Decreased_Gray_Max;
			bool check = false;
			for (int k = 0; k < 4; k++) {
				VECTOR_2D<int> r(ite->x + adjacent[k].x, ite->y + adjacent[k].y);
				if (_segments.get_zeropad(r.x, r.y) > 0
				    && _segments.get_zeropad(r.x, r.y) != num
				    && abs(color - decrease_color.get_zeropad(r.x, r.y)) < min) {
					check = true;
					min = fabs(color - decrease_color.get(r.x, r.y));
					decrease_color.at(ite->x, ite->y) = decrease_color.get(r.x, r.y);
					_segments.at(ite->x, ite->y) = _segments.get(r.x, r.y);
				}
			}
			if (check) {
				// remove from the list
				ite = small_region_coordinate.erase(ite);
			} else {
				// left current pixel and increment iterator
				++ite;
			}
		}
	}*/
	// Output vectors
	FILE *fp;
	FILE *fp_img;
	fp = fopen("segment_vector.dat", "w");
	fp_img = fopen("segments_color.pgm", "w");
	fprintf(fp, "%d %d\n", _width, _height);
	fprintf(fp_img, "P2\n%d %d\n%d\n", _width, _height, 255);
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
 *	kernel_spatial : the spatial radius of mean shift kernel
 *	kernel_intensity : the intensity threshold of mean shift kernel
 */
template <class T>
const VECTOR_2D<double>
Segmentation<T>::MeanShift_Grayscale(const int x, const int y, std::vector<VECTOR_2D<int> >& pel_list, int Iter_Max)
{
	VECTOR_2D<double> u;
	T intensity;

	// Initialize
	u.x = x;
	u.y = y;
	intensity = _image.get(x, y);
	// Iterate until it converge
	for (int i = 0; i < Iter_Max; i++) {
		double N = 0.0;
		T sum_intensity_diff = 0.0;
		VECTOR_2D<double> sum_d(0.0, 0.0);
		VECTOR_2D<double> d_tmp;

		for (unsigned int n = 0; n < pel_list.size(); n++) {
			VECTOR_2D<double> r(u.x + pel_list[n].x, u.y + pel_list[n].y);

			if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height) {
				double intensity_diff = _image.get(r.x, r.y) - intensity;

				if (fabs(intensity_diff) <= _kernel_intensity) {
					double coeff = (fabs(intensity_diff / _kernel_intensity) - 1.0) * (fabs(intensity_diff / _kernel_intensity) - 1.0);
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
		if (fabs(_image.get(u.x + d_tmp.x, u.y + d_tmp.y) - intensity) > _kernel_intensity / 2) {
			break;
		} else if (Vector_2D::norm(d_tmp) < 0.1) {
			u += d_tmp;
			break;
		}
		u += d_tmp;
	}
	return u;
}

