#include <cassert>
#include <cmath>
#include <list>




template <class T>
Segmentation::Segmentation(void)
{
	_width = 0;
	_height = 0;
	_kernel_spatial = 10.0;
	_kernel_intensity = 0.1;
}

template <class T>
Segmentation::Segmentation(const ImgVector<T>& image, const double kernel_spatial_radius, const double kernel_intensity_radius)
{
	_image.copy(image);
	if (_image.max() > 1.0) {
		// Normalize
		_image.contrast_stretching(0.0, 1.0);
	}
	_width = image.width();
	_height = image.height();
	_segments.reset(_width, _height);
	_kernel_spatial = kernel_spatial_radius;
	_kernel_intensity = kernel_intensity_radius;

	// Initial Segmentation
	Segmentation_MeanShift();
}

template <class T>
~Segmentation::Segmentation(void)
{
}




// ----- Accessor -----
int width(void) const
{
	return _width;
}

int height(void) const
{
	return _height;
}


int& operator[](int n)
{
	return _segments[n];
}

int& at(int n)
{
	assert(0 <= n && n < _width * _height);
	return _segments[n];
}

int& at(int x, int y)
{
	return _segments.at(x, y);
}

int& at_repeat(int x, int y)
{
	return _segments.at_repeat(x, y);
}

int& at_mirror(int x, int y)
{
	return _segments.at_mirror(x, y);
}


int
get(int n) const
{
	return _segments.get(n);
}

int
get(int x, int y) const
{
	return _segments.get(x, y);
}

int
get_zeropad(int x, int y) const
{
	return _segments.get_zeropad(x, y);
}

int
get_repeat(int x, int y) const
{
	return _segments.get_repeat(x, y);
}

int
get_mirror(int x, int y) const
{
	return _segments.get_mirror(x, y);
}




// ----- Mean Shift -----

template <class T>
void
Segmentation::Segmentation_MeanShift(int Iter_Max)
{
	ImgVector<VECTOR_2D<double> > shift_vector(_width, _height);

	if (_width <= 0 || _height <= 0) {
		return;
	}
	for (int i = 0; i < _width * _height; i++) {
		_segments[i] = i;
	}
	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			shift_vector.at(x, y) = MeanShift_Grayscale(x, y, Iter_Max);
			_segments.at(x, y) = _segments.at(round(shift_vector.at(x, y).x), round(shift_vector.at(x, y).y));
		}
	}
}


/*
 * std::vector<double> kernel has kernel radius for each dimensions.
 * The values it needs are below:
 *	kernel_spatial : the spatial radius of mean shift kernel
 *	kernel_intensity : the intensity threshold of mean shift kernel
 */
template <class T>
VECTOR_2D<double>
Segmentation::MeanShift_Grayscale(const int x, const int y, int Iter_Max)
{
	VECTOR_2D<double> u;
	T intensity;
	std::list<VECTOR_2D<int> > pel_list;

	// Make pixel list
	for (int m = -_kernel_spatial; m <= _kernel_spatial; m++) {
		for (int n = -_kernel_spatial; n <= _kernel_spatial; n++) {
			if (sqrt(n * n + m * m) < _kernel_spatial) {
				pel_list.push_back(VECTOR_2D<int>({n, m}));
			}
		}
	}
	// Initialize
	u = VECTOR_2D({x, y});
	intensity = _image.get(x, y);
	// Iterate until it converge
	for (int i = 0; i < Iter_Max; i++) {
		int N = 0;
		T sum_intensity = 0.0;
		VECTOR_2D<double> sum_r = 0;
		VECTOR_2D<double> u_tmp;

		for (std::list<VECTOR_2D<int> >::iterator ite = pel_list.begin(); ite != pel_list.end(); ++it) {
			VECTOR_2D<double> r(u.x + ite.x, u.y + ite.y);

			if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height
			    && fabs(intensity - _image.get(r.x, r.y)) <= kernel_intensity) {
				++N;
				sum_intensity += _image.get(r.x, r.y);
				sum_r += r;
			}
		}
		intensity = sum_intensity / double(N);
		u_tmp.x = sum_r.x / double(N);
		u_tmp.y = sum_r.y / double(N);
		if (Vector_2D::norm(u - u_tmp) < 0.5) {
			u = u_tmp;
			break;
		}
		u = u_tmp;
	}
	return u;
}

