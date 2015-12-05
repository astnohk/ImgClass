#include <cassert>
#include <cmath>
#include <iostream>
#include <new>
#include <stdexcept>




template <class T>
ImgVector<T>::ImgVector(void)
{
	_data = nullptr;
	_width = 0;
	_height = 0;
}


// Copy Constructor
template <class T>
ImgVector<T>::ImgVector(const ImgVector<T> &target)
{
	_data = nullptr;
	_width = 0;
	_height = 0;
	if (target._width > 0 && target._height > 0) {
		try {
			_data = new T[target._width * target._height]();
		}
		catch (const std::bad_alloc &bad) {
			std::cerr << "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory" << std::endl;
			_data = nullptr;
			throw;
		}
		_width = target._width;
		_height = target._height;
		for (int i = 0; i < _width * _height; i++) {
			_data[i] = target._data[i];
		}
	}
}

template <class T>
ImgVector<T>::ImgVector(int W, int H, const T& value)
{
	_data = nullptr;
	_width = 0;
	_height = 0;
	if (W > 0 && H > 0) {
		try {
			_data = new T[W * H]();
		}
		catch (const std::bad_alloc &bad) {
			std::cerr << "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory" << std::endl;
			_data = nullptr;
			throw;
		}
		_width = W;
		_height = H;
		for (int i = 0; i < _width * _height; i++) {
			_data[i] = value;
		}
	}
}

template <class T>
ImgVector<T>::ImgVector(int W, int H, const T* array)
{
	_data = nullptr;
	_width = 0;
	_height = 0;
	if (W > 0 && H > 0) {
		try {
			_data = new T[W * H]();
		}
		catch (const std::bad_alloc &bad) {
			std::cerr << "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory" << std::endl;
			_data = nullptr;
			throw;
		}
		_width = W;
		_height = H;
		if (array != nullptr) {
			for (int i = 0; i < _width * _height; i++) {
				_data[i] = array[i];
			}
		}
	}
}


template <class T>
ImgVector<T>::~ImgVector(void)
{
	delete[] _data;
}


template <class T>
void
ImgVector<T>::clear(void)
{
	delete[] _data;
	_data = nullptr;
	_width = 0;
	_height = 0;
}


template <class T>
void
ImgVector<T>::reset(int W, int H, const T& value)
{
	delete[] _data;
	_data = nullptr;
	_width = 0;
	_height = 0;
	if (W > 0 && H > 0) {
		try {
			_data = new T[W * H]();
		}
		catch (const std::bad_alloc &bad) {
			std::cerr << "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory" << std::endl;
			_data = nullptr;
			throw;
		}
		_width = W;
		_height = H;
		for (int i = 0; i < _width * _height; i++) {
			_data[i] = value;
		}
	}
}

template <class T>
void
ImgVector<T>::reset(int W, int H, const T* array)
{
	delete[] _data;
	_data = nullptr;
	_width = 0;
	_height = 0;
	if (W > 0 && H > 0) {
		try {
			_data = new T[W * H]();
		}
		catch (const std::bad_alloc &bad) {
			std::cerr << "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory" << std::endl;
			_data = nullptr;
			throw;
		}
		_width = W;
		_height = H;
		if (array != nullptr) {
			for (int i = 0; i < _width * _height; i++) {
				_data[i] = array[i];
			}
		}
	}
}


template <class T>
ImgVector<T> &
ImgVector<T>::copy(const ImgVector<T>& vector)
{
	if (this != &vector
	    && vector._width > 0 && vector._height > 0) {
		T *tmp_data = nullptr;
		try {
			tmp_data = new T[vector._width * vector._height]();
		}
		catch (const std::bad_alloc &bad) {
			std::cerr << "ImgVector::copy(const ImgVector<T>&) : Cannot Allocate Memory" << std::endl;
			throw;
			return *this;
		}
		_width = vector._width;
		_height = vector._height;
		delete[] _data;
		_data = tmp_data;
		for (int i = 0; i < _width * _height; i++) {
			_data[i] = vector._data[i];
		}
	}
	return *this;
}

template <class T>
ImgVector<T> &
ImgVector<T>::operator=(const ImgVector<T>& vector)
{
	if (this != &vector
	    && vector._width > 0 && vector._height > 0) {
		T *tmp_data = nullptr;
		try {
			tmp_data = new T[vector._width * vector._height]();
		}
		catch (const std::bad_alloc &bad) {
			std::cerr << "ImgVector::operator=(ImgVector<T>&) : Cannot Allocate Memory" << std::endl;
			throw;
			return *this;
		}
		_width = vector._width;
		_height = vector._height;
		delete[] _data;
		_data = tmp_data;
		for (int i = 0; i < _width * _height; i++) {
			_data[i] = vector._data[i];
		}
	}
	return *this;
}




template <class T>
int
ImgVector<T>::width(void) const
{
	return _width;
}


template <class T>
int
ImgVector<T>::height(void) const
{
	return _height;
}


template <class T>
int
ImgVector<T>::size(void) const
{
	return _width * _height;
}


template <class T>
bool
ImgVector<T>::isNULL(void) const
{
	if (_data == nullptr) {
		return true;
	} else {
		return false;
	}
}




template <class T>
T *
ImgVector<T>::data(void) const
{
	return _data;
}


template <class T>
T &
ImgVector<T>::operator[](int n)
{
	return _data[n];
}

template <class T>
T &
ImgVector<T>::at(int x, int y)
{
	assert(0 <= x && x < _width && 0 <= y && y < _height);
	return _data[_width * y + x];
}

template <class T>
T &
ImgVector<T>::at_repeat(int x, int y)
{
	int x_repeat, y_repeat;

	if (x >= 0) {
		x_repeat = x % _width;
	} else {
		x_repeat = _width - ((int)std::abs((double)x + 1.0) % _width);
	}
	if (y >= 0) {
		y_repeat = y % _height;
	} else {
		y_repeat = _height - ((int)std::abs((double)y + 1.0) % _height);
	}
	return _data[_width * y_repeat + x_repeat];
}

template <class T>
T &
ImgVector<T>::at_mirror(int x, int y)
{
	int x_mirror, y_mirror;

	if (x < 0) {
		x = -x - 1; // Mirroring over negative has offset
	}
	if (y < 0) {
		y = -y - 1;
	}
	x_mirror = (int)round(_width - 0.5 - std::fabs(_width - 0.5 - (x % (2 * _width))));
	y_mirror = (int)round(_height - 0.5 - std::fabs(_height - 0.5 - (y % (2 * _height))));
	return _data[_width * y_mirror + x_mirror];
}


template <class T>
const T
ImgVector<T>::get(int n) const
{
	assert(0 <= n && n < _width * _height);
	return _data[n];
}

template <class T>
const T
ImgVector<T>::get(int x, int y) const
{
	assert(0 <= x && x < _width && 0 <= y && y < _height);
	return _data[_width * y + x];
}

template <class T>
const T
ImgVector<T>::get_zeropad(int x, int y) const
{
	T zero = T();

	if (x < 0 || _width <= x || y < 0 || _height <= y) {
		return zero;
	} else {
		return _data[_width * y + x];
	}
}

template <class T>
const T
ImgVector<T>::get_repeat(int x, int y) const
{
	int x_repeat, y_repeat;

	if (x >= 0) {
		x_repeat = x % _width;
	} else {
		x_repeat = _width - ((int)std::abs((double)x + 1.0) % _width);
	}
	if (y >= 0) {
		y_repeat = y % _height;
	} else {
		y_repeat = _height - ((int)std::abs((double)y + 1.0) % _height);
	}
	return _data[_width * y_repeat + x_repeat];
}

template <class T>
const T
ImgVector<T>::get_mirror(int x, int y) const
{
	int x_mirror, y_mirror;

	if (x < 0) {
		x = -x - 1; // Mirroring over negative has offset
	}
	if (y < 0) {
		y = -y - 1;
	}
	x_mirror = (int)round(_width - 0.5 - std::fabs(_width - 0.5 - (x % (2 * _width))));
	y_mirror = (int)round(_height - 0.5 - std::fabs(_height - 0.5 - (y % (2 * _height))));
	return _data[_width * y_mirror + x_mirror];
}




// Get continuous function interpolated by bicubic
template <class T>
const T
ImgVector<T>::get_zeropad(double x, double y, double B, double C) const
{
	double bicubic_x[4];
	double bicubic_y[4];
	T value = T();

	for (int n = 0; n < 4; n++) {
		bicubic_x[n] = ImgVector<double>::cubic(n - 1.0 - (x - floor(x)), B, C);
		bicubic_y[n] = ImgVector<double>::cubic(n - 1.0 - (y - floor(y)), B, C);
	}
	for (int m = 0; m < 4; m++) {
		for (int n = 0; n < 4; n++) {
			value = value
			    + this->get_zeropad((int)floor(x), (int)floor(y)) * bicubic_x[n] * bicubic_y[m];
		}
	}
	return value;
}

template <class T>
const T
ImgVector<T>::get_repeat(double x, double y, double B, double C) const
{
	double bicubic_x[4];
	double bicubic_y[4];
	T value = T();

	for (int n = 0; n < 4; n++) {
		bicubic_x[n] = ImgVector<double>::cubic(n - 1.0 - (x - floor(x)), B, C);
		bicubic_y[n] = ImgVector<double>::cubic(n - 1.0 - (y - floor(y)), B, C);
	}
	for (int m = 0; m < 4; m++) {
		for (int n = 0; n < 4; n++) {
			value = value
			    + this->get_repeat((int)floor(x), (int)floor(y)) * bicubic_x[n] * bicubic_y[m];
		}
	}
	return value;
}

template <class T>
const T
ImgVector<T>::get_mirror(double x, double y, double B, double C) const
{
	double bicubic_x[4];
	double bicubic_y[4];
	T value = T();

	for (int n = 0; n < 4; n++) {
		bicubic_x[n] = ImgVector<double>::cubic(n - 1.0 - (x - floor(x)), B, C);
		bicubic_y[n] = ImgVector<double>::cubic(n - 1.0 - (y - floor(y)), B, C);
	}
	for (int m = 0; m < 4; m++) {
		for (int n = 0; n < 4; n++) {
			value = value
			    + this->get_mirror((int)floor(x), (int)floor(y)) * bicubic_x[n] * bicubic_y[m];
		}
	}
	return value;
}




template <class T>
const T
ImgVector<T>::min(void) const
{
	if (_width <= 0 || _height <= 0) {
		throw std::logic_error("T ImgVector<T>::min(void) : vector is empty");
	}
	T min = _data[0];
	for (int i = 1; i < _width * _height; i++) {
		if (_data[i] < min) {
			min = _data[i];
		}
	}
	return min;
}

template <class T>
const T
ImgVector<T>::max(void) const
{
	if (_width <= 0 || _height <= 0) {
		throw std::logic_error("T ImgVector<T>::max(void) : vector is empty");
	}
	T max = _data[0];
	for (int i = 1; i < _width * _height; i++) {
		if (_data[i] > max) {
			max = _data[i];
		}
	}
	return max;
}


template <class T>
const T
ImgVector<T>::min(int top_left_x, int top_left_y, int crop_width, int crop_height) const
{
	if (_width <= 0 || _height <= 0) {
		throw std::logic_error("T ImgVector<T>::min(void) : vector is empty");
	} else if (crop_width <= 0) {
		throw std::invalid_argument("T ImgVector<T>::min(int, int, int, int) : crop_width <= 0");
	} else if (crop_height <= 0) {
		throw std::invalid_argument("T ImgVector<T>::min(int, int, int, int) : crop_height<= 0");
	}
	T min = _data[_width * top_left_y + top_left_x];
	for (int y = 0; y < crop_height; y++) {
		for (int x = 0; x < crop_width; x++) {
			if (0 <= x + top_left_x && x + top_left_x < _width
			    && 0 <= y + top_left_y && y + top_left_y < _height
			    && _data[_width * (y + top_left_y) + x + top_left_x] < min) {
				min = _data[_width * (y + top_left_y) + x + top_left_x];
			}
		}
	}
	return min;
}

template <class T>
const T
ImgVector<T>::max(int top_left_x, int top_left_y, int crop_width, int crop_height) const
{
	if (_width <= 0 || _height <= 0) {
		throw std::logic_error("T ImgVector<T>::max(void) : vector is empty");
	} else if (crop_width <= 0) {
		throw std::invalid_argument("T ImgVector<T>::max(int, int, int, int) : crop_width <= 0");
	} else if (crop_height <= 0) {
		throw std::invalid_argument("T ImgVector<T>::max(int, int, int, int) : crop_height <= 0");
	}
	T max = _data[_width * top_left_y + top_left_x];
	for (int y = 0; y < crop_height; y++) {
		for (int x = 0; x < crop_width; x++) {
			if (0 <= x + top_left_x && x + top_left_x < _width
			    && 0 <= y + top_left_y && y + top_left_y < _height
			    && _data[_width * (y + top_left_y) + x + top_left_x] > max) {
				max = _data[_width * (y + top_left_y) + x + top_left_x];
			}
		}
	}
	return max;
}




// ----- Statistics -----
template <class T>
const T
ImgVector<T>::variance() const
{
	double N = double(_width * _height);
	T sum_squared = .0;
	T sum = .0;

	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			sum += _data[_width * y + x];
			sum_squared += _data[_width * y + x] * _data[_width * y + x];
		}
	}
	return sum_squared / N - sum * sum / (N * N);
}

template <class T>
const T
ImgVector<T>::variance(const int top_left_x, const int top_left_y, const int crop_width, const int crop_height) const
{
	double N = .0;
	T sum_squared = .0;
	T sum = .0;

	for (int y = top_left_y; y < top_left_y + crop_height; y++) {
		for (int x = top_left_x; x < top_left_x + crop_width; x++) {
			if (0 <= x && x < _width
			    && 0 <= y && y < _height) {
				N += 1.0;
				sum += _data[_width * y + x];
				sum_squared += _data[_width * y + x] * _data[_width * y + x];
			}
		}
	}
	return sum_squared / N - sum * sum / (N * N);
}




// ----- Image Operation -----
template <class T>
ImgVector<T> *
ImgVector<T>::crop(int top_left_x, int top_left_y, int crop_width, int crop_height) const
{
	ImgVector<T>* tmp = nullptr;

	if (crop_width <= 0) {
		throw std::invalid_argument("ImgVector<T>* ImgVector<T>::crop(int, int, int, int) : crop_width <= 0");
	} else if (crop_height <= 0) {
		throw std::invalid_argument("ImgVector<T>* ImgVector<T>::crop(int, int, int, int) : crop_height <= 0");
	}
	try {
		tmp = new ImgVector<T>;
	}
	catch (const std::bad_alloc &bad) {
		std::cerr << "ImgVector<T>* ImgVector<T>::crop(int, int, int, int) : Cannot allocate memory" << std::endl;
		throw;
	}
	// Initialize
	tmp->reset(crop_width, crop_height);
	// Crop
	for (int y = 0; y < crop_width; y++) {
		for (int x = 0; x < crop_height; x++) {
			if (0 <= top_left_y + y && top_left_y + y < _height
			    && 0 <= top_left_x + x && top_left_x + x < _width) {
				tmp->at(x, y) = _data[_width * (y + top_left_y) + x + top_left_x];
			} else {
				tmp->at(x, y) = 0;
			}
		}
	}
	return tmp;
}




template <class T>
void
ImgVector<T>::contrast_stretching(const T& Min, const T& Max)
{
	if (Min > Max) {
		throw std::invalid_argument("void ImgVector<T>::contrast_stretching(const T&, const T&) : min > max");
	}
	// Initialize
	T min_tmp = _data[0];
	T max_tmp = _data[0];

	for (int i = 1; i < _width * _height; i++) {
		if (_data[i] < min_tmp) {
			min_tmp = _data[i];
		} else if (_data[i] > max_tmp) {
			max_tmp = _data[i];
		}
	}
	// Stretching
	for (int i = 0; i < _width * _height; i++) {
		_data[i] = Min + (_data[i] - min_tmp) * Max / (max_tmp - min_tmp);
	}
}




template <class T>
void
ImgVector<T>::resize_zerohold(int W, int H)
{
	T *resized = nullptr;
	T additive_identity = T();
	double scale_x = .0;
	double scale_y = .0;
	int area_x;
	int area_y;
	int m, n;
	int x, y;
	T sum;

	if (W <= 0) {
		throw std::out_of_range("ImgVector<T>::resize_zerohold(int, int) : int W");
	}else if (H <= 0) {
		throw std::out_of_range("ImgVector<T>::resize_zerohold(int, int) : int H");
	}
	scale_x = (double)W / _width;
	scale_y = (double)H / _height;
	try {
		resized = new T[W * H]();
	}
	catch (const std::bad_alloc &bad) {
		std::cerr << "ImgVector<T>::resize_zerohold(int, int) : Cannot allocate memory" << std::endl;
		throw;
	}
	area_x = ceil((double)_width / W);
	area_y = ceil((double)_height / H);
	for (y = 0; y < H; y++) {
		for (x = 0; x < W; x++) {
			sum = additive_identity;
			for (m = 0; m < area_y; m++) {
				for (n = 0; n < area_x; n++) {
					sum += this->get((int)floor(x / scale_x) + n, (int)floor(y / scale_y) + m);
				}
			}
			resized[W * y + x] = sum / (area_x * area_y);
		}
	}
	delete[] _data;
	_data = resized;
	_width = W;
	_height = H;
}


/*
    bool ImgVector<T>::resize_bicubic(int W, int H, double min, double max, T (*Nearest_Integer_Method)(double &d), double A)
    int W, int H : width and height of resized image
    double min : minimum value of saturated value
    double max : maximum value of saturated value
    T (*Nearest_Integer_Method)(double &d) : round method (e.g. floor(), round(), etc.)
    A : cubic method's parameter (default A = -0.5 which correspond to Hermite)
*/
template <class T>
void
ImgVector<T>::resize_bicubic(int W, int H, double min, double max, T (*Nearest_Integer_Method)(double &d), double B, double C)
{
	T *resized = nullptr;
	double *conv = nullptr;
	ImgVector<double> Tmp;
	double scale_x, scale_y;
	double scale_conv;
	int L, L_center;
	double dx, dy;
	int x, y;
	int m, n;
	double sum;

	if (W <= 0) {
		throw std::out_of_range("ImgVector<T>::resize_bicubic(int, int, double, double, T (*)(double &d), double, double) :int W");
	} else if (H <= 0) {
		throw std::out_of_range("ImgVector<T>::resize_bicubic(int, int, double, double, T (*)(double &d), double, double) :int H");
	}
	scale_x = (double)W / _width;
	scale_y = (double)H / _height;
	Tmp.reset(W, _height);
	// The length of cubic convolution coefficient
	scale_conv = 1.0;
	if (scale_x < 1.0 || scale_y < 1.0) {
		scale_conv = ceil(1.0 / (scale_x < scale_y ? scale_x : scale_y));
	}
	try {
		resized = new T[W * H]();
		conv = new double[(int)scale_conv * 4]();
	}
	catch (const std::bad_alloc &bad) {
		std::cerr << "ImgVector<double>::resize_bicubic(int, int) error : Cannot allocate memory" << std::endl;
		delete[] resized;
		delete[] conv;
		throw;
	}
	// Horizontal convolution
	if (scale_x >= 1.0) {
		scale_conv = 1.0;
		L = 4;
		L_center = floor((L - 1.0) / 2);
	} else {
		scale_conv = 1.0 / scale_x;
		L = 4 * (int)ceil(scale_conv);
		L_center = floor((L - 1.0) / 2);
	}
	for (x = 0; x < W; x++) {
		if (scale_x >= 1.0) {
			dx = (x - (scale_x - 1.0) / 2.0) / scale_x;
			for (n = 0; n < L; n++) {
				conv[n] = ImgVector<T>::cubic((double)(n - L_center) - (dx - floor(dx)), B, C);
			}
		} else {
			dx = x / scale_x + (1.0 / scale_x - 1.0) / 2.0;
			for (n = 0; n < L; n++) {
				conv[n] = ImgVector<T>::cubic(((double)(n - L_center) - (dx - floor(dx))) * scale_x, B, C) / scale_conv;
			}
		}
		for (y = 0; y < _height; y++) {
			sum = 0.0;
			for (n = 0; n < L; n++) {
				sum += conv[n] * this->get_mirror((int)floor(dx) + n - L_center, y);
			}
			Tmp.at(x, y) = sum;
		}
	}
	// Vertical convolution
	if (scale_y >= 1.0) {
		scale_conv = 1.0;
		L = 4;
		L_center = floor((L - 1.0) / 2);
	} else {
		scale_conv = 1.0 / scale_y;
		L = 4 * (int)ceil(scale_conv);
		L_center = floor((L - 1.0) / 2);
	}
	for (y = 0; y < H; y++) {
		if (scale_y >= 1.0) {
			dy = (y - (scale_y - 1.0) / 2.0) / scale_y;
			for (m = 0; m < L; m++) {
				conv[m] = ImgVector<T>::cubic((double)(m - L_center) - (dy - floor(dy)), B, C);
			}
		} else {
			dy = y / scale_y + (1.0 / scale_y - 1.0) / 2.0;
			for (m = 0; m < L; m++) {
				conv[m] = ImgVector<T>::cubic(((double)(m - L_center) - (dy - floor(dy))) / scale_conv, B, C) / scale_conv;
			}
		}
		for (x = 0; x < W; x++) {
			sum = 0.0;
			for (m = 0; m < L; m++) {
				sum += conv[m] * Tmp.get_mirror(x, (int)floor(dy) + m - L_center);
			}
			if (min != max) {
				sum = sum >= min ? sum <= max ? sum : max : min;
			}
			if (Nearest_Integer_Method != nullptr) {
				resized[W * y + x] = Nearest_Integer_Method(sum);
			} else {
				resized[W * y + x] = sum;
			}
		}
	}
	delete[] conv;
	delete[] _data;
	_data = resized;
	_width = W;
	_height = H;
	return true;
}


template <class T>
double
ImgVector<T>::cubic(double x, double B, double C) const
{
	double x_abs = fabs(x);

	if (x_abs <= 1.0) {
		return ((2.0 - 1.5 * B - C) * x_abs + (-3.0 + 2.0 * B + C)) * x_abs * x_abs + 1.0 - B / 3.0;
	} else if (x_abs < 2.0) {
		return (((-B / 6.0 - C) * x_abs + B + 5.0 * C) * x_abs - 2.0 * B - 8.0 * C) * x_abs + 8.0 / 6.0 * B + 4.0 * C;
	} else {
		return 0.0;
	}
}


template <class T>
void
ImgVector<T>::map(T (*func)(T &value))
{
	if (func != nullptr) {
		for (int i = 0 ; i < _width * _height; i++) {
			_data[i] = func(_data[i]);
		}
	}
}




// ----- Arithmetic Operators -----
template <class T>
ImgVector<T> &
ImgVector<T>::operator+=(const T& val)
{
	for (int i = 0; i < _width * _height; i++) {
		_data[i] += val;
	}
	return *this;
}

template <class T>
ImgVector<T> &
ImgVector<T>::operator-=(const T& val)
{
	for (int i = 0; i < _width * _height; i++) {
		_data[i] -= val;
	}
	return *this;
}

template <class T>
ImgVector<T> &
ImgVector<T>::operator*=(const T& val)
{
	for (int i = 0; i < _width * _height; i++) {
		_data[i] *= val;
	}
	return *this;
}

template <class T>
ImgVector<T> &
ImgVector<T>::operator/=(const T& val)
{
	for (int i = 0; i < _width * _height; i++) {
		_data[i] /= val;
	}
	return *this;
}

