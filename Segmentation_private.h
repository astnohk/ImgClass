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
Segmentation<T>::Segmentation(const Segmentation<T>& segments) // Copy constructor
{
	_width = segments._width;
	_height = segments._height;
	_kernel_spatial = segments._kernel_spatial;
	_kernel_intensity = segments._kernel_intensity;

	_image.copy(segments._image);
	_segments.copy(segments._segments);
}

template <class T>
Segmentation<T>::Segmentation(const ImgVector<T>* image, const double kernel_spatial_radius, const double kernel_intensity_radius)
{
	_image.copy(image);
	if (_image.max() > 1.0) {
		// Normalize
		_image.contrast_stretching(0.0, 1.0);
	}
	_width = _image.width();
	_height = _image.height();
	_segments.reset(_width, _height);
	_kernel_spatial = kernel_spatial_radius;
	_kernel_intensity = kernel_intensity_radius;

	// Initial Segmentation
	Segmentation_MeanShift();
}


template <class T>
Segmentation<T> &
Segmentation<T>::reset(const ImgVector<T>* image, const double kernel_spatial_radius, const double kernel_intensity_radius)
{
	_image.copy(image);
	if (_image.max() > 1.0) {
		// Normalize
		_image.contrast_stretching(0.0, 1.0);
	}
	_width = _image.width();
	_height = _image.height();
	_segments.reset(_width, _height);
	_kernel_spatial = kernel_spatial_radius;
	_kernel_intensity = kernel_intensity_radius;

	// Initial Segmentation
	Segmentation_MeanShift();
	return *this;
}

template <class T>
Segmentation<T> &
Segmentation<T>::copy(const Segmentation<T>* segments)
{
	_width = segments->_width;
	_height = segments->_height;
	_kernel_spatial = segments->_kernel_spatial;
	_kernel_intensity = segments->_kernel_intensity;

	_image.copy(segments->_image);
	_segments.copy(segments->_segments);
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
Segmentation<T>::Segmentation_MeanShift(int Iter_Max)
{
	ImgVector<VECTOR_2D<double> > shift_vector(_width, _height);

	if (_width <= 0 || _height <= 0) {
		return;
	}
	for (int i = 0; i < _width * _height; i++) {
		_segments[i] = i;
	}
	FILE *fp;
	fp = fopen("segment_vector.dat", "w");
	fprintf(fp, "%d %d\n", _width, _height);
	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			shift_vector.at(x, y) = MeanShift_Grayscale(x, y, Iter_Max);
			_segments.at(x, y) = _segments.at(round(shift_vector.at(x, y).x), round(shift_vector.at(x, y).y));
			shift_vector.at(x, y).x -= x;
			shift_vector.at(x, y).y -= y;
			fwrite(&(shift_vector.get(x, y).x), sizeof(double), 1, fp);
			fwrite(&(shift_vector.get(x, y).y), sizeof(double), 1, fp);
		}
	}
	fclose(fp);
}

/*
 * std::vector<double> kernel has kernel radius for each dimensions.
 * The values it needs are below:
 *	kernel_spatial : the spatial radius of mean shift kernel
 *	kernel_intensity : the intensity threshold of mean shift kernel
 */
template <class T>
VECTOR_2D<double>
Segmentation<T>::MeanShift_Grayscale(const int x, const int y, int Iter_Max)
{
	VECTOR_2D<double> u;
	T intensity;
	std::list<VECTOR_2D<int> > pel_list;

	// Make pixel list
	for (int m = -_kernel_spatial; m <= _kernel_spatial; m++) {
		for (int n = -_kernel_spatial; n <= _kernel_spatial; n++) {
			if (sqrt(n * n + m * m) < _kernel_spatial) {
				VECTOR_2D<int> r(n, m);
				pel_list.push_back(r);
			}
		}
	}
	// Initialize
	u.x = x;
	u.y = y;
	intensity = _image.get(x, y);
	// Iterate until it converge
	for (int i = 0; i < Iter_Max; i++) {
		int N = 0;
		T sum_intensity = 0.0;
		VECTOR_2D<double> sum_r(0.0, 0.0);
		VECTOR_2D<double> u_tmp;

		for (std::list<VECTOR_2D<int> >::iterator ite = pel_list.begin(); ite != pel_list.end(); ++ite) {
			VECTOR_2D<double> r(u.x + ite->x, u.y + ite->y);

			if (x == 0 && y == 0)
				printf("intensity : %f - %f\n", intensity, _image.get_zeropad(r.x, r.y));

			if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height
			    && fabs(intensity - _image.get(r.x, r.y)) <= _kernel_intensity) {
				++N;
				sum_intensity += _image.get(r.x, r.y);
				sum_r += r;
			}
		}
		intensity = sum_intensity / double(N);
		u_tmp.x = sum_r.x / double(N);
		u_tmp.y = sum_r.y / double(N);
		if (Vector_2D::norm(u - u_tmp) < 0.01) {
			u = u_tmp;
			break;
		}
		u = u_tmp;
	}
	return u;
}

