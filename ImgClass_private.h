#include <cassert>
#include <cmath>
#include <iostream>
#include <new>
#include <stdexcept>




template <class T>
ImgVector<T>::ImgVector(void)
{
	_data = nullptr;
	_reserved_size = 0;
	_width = 0;
	_height = 0;
}


template <class T>
ImgVector<T>::ImgVector(const int Width, const int Height, const T& value)
{
	_data = nullptr;
	_reserved_size = 0;
	_width = 0;
	_height = 0;
	if (Width > 0 && Height > 0) {
		_reserved_size = size_t(Width) * size_t(Height);
		try {
			_data = new T[_reserved_size];
		}
		catch (const std::bad_alloc& bad) {
			std::cerr << bad.what() << std::endl
			    << "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory" << std::endl;
			_data = nullptr;
			throw;
		}
		_width = Width;
		_height = Height;
		for (size_t n = 0; n < _reserved_size; n++) {
			_data[n] = value;
		}
	}
}

template <class T>
ImgVector<T>::ImgVector(const int Width, const int Height, const T* array)
{
	_data = nullptr;
	_reserved_size = 0;
	_width = 0;
	_height = 0;
	if (Width > 0 && Height > 0) {
		_reserved_size = size_t(Width) * size_t(Height);
		try {
			_data = new T[_reserved_size];
		}
		catch (const std::bad_alloc& bad) {
			std::cerr << bad.what() << std::endl
			    << "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory" << std::endl;
			_data = nullptr;
			throw;
		}
		_width = Width;
		_height = Height;
		if (array != nullptr) {
			for (size_t n = 0; n < _reserved_size; n++) {
				_data[n] = array[n];
			}
		}
	}
}


template <class T>
ImgVector<T>::ImgVector(const ImgVector<T>& copy)
{
	_data = nullptr;
	_reserved_size = 0;
	_width = 0;
	_height = 0;
	if (copy._width > 0 && copy._height > 0) {
		_reserved_size = copy.size();
		try {
			_data = new T[_reserved_size];
		}
		catch (const std::bad_alloc& bad) {
			std::cerr << bad.what() << std::endl
			    << "ImgVector::ImgVector(const ImgVector<T>&) : Cannot Allocate Memory" << std::endl;
			_data = nullptr;
			throw;
		}
		_width = copy._width;
		_height = copy._height;
		for (size_t n = 0; n < _reserved_size; n++) {
			_data[n] = copy._data[n];
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
ImgVector<T>::reserve(const int Width, const int Height)
{
	size_t new_size = size_t(Width) * size_t(Height);
	if (_reserved_size < new_size) {
		T* new_data = nullptr;
		try {
			new_data = new T[new_size];
		}
		catch (const std::bad_alloc& bad) {
			std::cerr << bad.what() << std::endl
			    << "ImgVector<T>::reserve(const int, const int) : Cannot allocate Memory" << std::endl;
			throw;
		}
		_reserved_size = new_size;
		for (size_t y = 0; y < size_t(_height); y++) {
			for (size_t x = 0; x < size_t(_width); x++) {
				new_data[size_t(_width) * y + x] = _data[size_t(_width) * y + x];
			}
		}
		delete[] _data;
		_data = new_data;
	}
}


template <class T>
void
ImgVector<T>::clear(void)
{
	delete[] _data;
	_data = nullptr;
	_reserved_size = 0;
	_width = 0;
	_height = 0;
}


template <class T>
void
ImgVector<T>::reset(const int Width, const int Height, const T& value)
{
	if (Width > 0 && Height > 0) {
		size_t new_size = size_t(Width) * size_t(Height);
		if (_reserved_size < new_size) {
			T* new_data = nullptr;
			try {
				new_data = new T[new_size];
			}
			catch (const std::bad_alloc& bad) {
				std::cerr << bad.what() << std::endl
				    << "ImgVector::reset(const int, const int, const T&) : Cannot Allocate Memory" << std::endl;
				throw;
			}
			delete[] _data;
			_data = new_data;
			_reserved_size = new_size;
		}
		_width = Width;
		_height = Height;
		for (size_t n = 0; n < new_size; n++) {
			_data[n] = value;
		}
	}
}

template <class T>
void
ImgVector<T>::reset(const int Width, const int Height, const T* array)
{
	if (Width > 0 && Height > 0) {
		size_t new_size = size_t(Width) * size_t(Height);
		if (_reserved_size < new_size) {
			T* new_data = nullptr;
			try {
				new_data = new T[new_size];
			}
			catch (const std::bad_alloc& bad) {
				std::cerr << bad.what() << std::endl
				    << "ImgVector::reset(const int, const int, const T*) : Cannot Allocate Memory" << std::endl;
				throw;
			}
			_data = new_data;
			_reserved_size = new_size;
		}
		_width = Width;
		_height = Height;
		if (array != nullptr) {
			for (size_t n = 0; n < new_size; n++) {
				_data[n] = array[n];
			}
		}
	}
}


template <class T>
void
ImgVector<T>::resize(const int Width, const int Height, const T& value)
{
	if (Width > 0 && Height > 0) {
		size_t new_size = size_t(Width) * size_t(Height);
		if (_reserved_size < new_size) { // New size is greater than previous size
			T *new_data = nullptr;
			try {
				new_data = new T[new_size];
			}
			catch (const std::bad_alloc& bad) {
				std::cerr << bad.what() << std::endl
				    << "ImgVector::resize(const int, const int, const T&) : Cannot Allocate Memory" << std::endl;
				throw;
			}
			for (size_t y = 0; y < size_t(Height); y++) {
				for (size_t x = 0; x < size_t(Width); x++) {
					if (y < size_t(_height) && x < size_t(_width)) {
						new_data[size_t(Width) * y + x] =
						    _data[size_t(_width) * y + x];
					} else {
						new_data[size_t(Width) * y + x] = value;
					}
				}
			}
			delete[] _data;
			_data = new_data;
			_reserved_size = new_size;
		} else if (Width > _width) { // New size is less than or equal to previous but new_width > previous_width
			for (int y = Height - 1; y >= 0; y--) {
				for (int x = Width - 1; x >= 0; x--) {
					if (y < _height && x < _width) {
						_data[size_t(Width) * size_t(y) + size_t(x)] =
						    _data[size_t(_width) * size_t(y) + size_t(x)];
					} else {
						_data[size_t(Width) * size_t(y) + size_t(x)] = value;
					}
				}
			}
		} else {
			for (size_t y = 0; y < size_t(Height); y++) {
				for (size_t x = 0; x < size_t(Width); x++) {
					if (y < size_t(_height) && x < size_t(_width)) {
						_data[size_t(Width) * y + x] = _data[size_t(_width) * y + x];
					} else {
						_data[size_t(Width) * y + x] = value;
					}
				}
			}
		}
	} else {
		delete[] _data;
		_data = nullptr;
	}
	_width = Width;
	_height = Height;
}


template <class T>
ImgVector<T> &
ImgVector<T>::copy(const ImgVector<T>& vector)
{
	if (this != &vector
	    && vector._width > 0 && vector._height > 0) {
		size_t new_size = vector.size();
		if (_reserved_size < new_size) {
			T *new_data = nullptr;
			try {
				new_data = new T[new_size];
			}
			catch (const std::bad_alloc& bad) {
				std::cerr << bad.what() << std::endl
				    << "ImgVector::copy(const ImgVector<T>&) : Cannot Allocate Memory" << std::endl;
				throw;
			}
			delete[] _data;
			_data = new_data;
			_reserved_size = new_size;
		}
		_width = vector._width;
		_height = vector._height;
		for (size_t n = 0; n < new_size; n++) {
			_data[n] = vector._data[n];
		}
	}
	return *this;
}

template <class T>
template <class RT>
ImgVector<T> &
ImgVector<T>::cast_copy(const ImgVector<RT>& vector)
{
	if (vector.width() > 0 && vector.height() > 0) {
		size_t new_size = vector.size();
		if (_reserved_size < new_size) {
			T *new_data = nullptr;
			try {
				new_data = new T[new_size];
			}
			catch (const std::bad_alloc& bad) {
				std::cerr << bad.what() << std::endl
				    << "ImgVector::operator=(ImgVector<T>&) : Cannot Allocate Memory" << std::endl;
				throw;
			}
			delete[] _data;
			_data = new_data;
			_reserved_size = new_size;
		}
		_width = vector.width();
		_height = vector.height();
		for (size_t n = 0; n < new_size; n++) {
			_data[n] = T(vector.get(n)); // copy with cast
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
		size_t new_size = vector.size();
		if (_reserved_size < new_size) {
			T *new_data = nullptr;
			try {
				new_data = new T[new_size];
			}
			catch (const std::bad_alloc& bad) {
				std::cerr << bad.what() << std::endl
				    << "ImgVector::operator=(ImgVector<T>&) : Cannot Allocate Memory" << std::endl;
				throw;
			}
			delete[] _data;
			_data = new_data;
			_reserved_size = new_size;
		}
		_width = vector._width;
		_height = vector._height;
		for (size_t n = 0; n < new_size; n++) {
			_data[n] = vector.get(n); // copy with cast
		}
	}
	return *this;
}




// ----- Accessors -----
template <class T>
size_t
ImgVector<T>::reserved_size(void) const
{
	return _reserved_size;
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
size_t
ImgVector<T>::size(void) const
{
	return size_t(_width) * size_t(_height);
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
ImgVector<T>::operator[](const size_t n)
{
	return _data[n];
}

template <class T>
const T &
ImgVector<T>::operator[](const size_t n) const
{
	return _data[n];
}


template <class T>
T &
ImgVector<T>::at(const size_t n)
{
	assert(0 <= n && n < size_t(_width) * size_t(_height));
	return _data[n];
}

template <class T>
const T &
ImgVector<T>::at(const size_t n) const
{
	assert(0 <= n && n < size_t(_width) * size_t(_height));
	return _data[n];
}


template <class T>
T &
ImgVector<T>::at(const int x, const int y)
{
	assert(0 <= x && x < _width && 0 <= y && y < _height);
	return _data[size_t(_width) * size_t(y) + size_t(x)];
}

template <class T>
const T &
ImgVector<T>::at(const int x, const int y) const
{
	assert(0 <= x && x < _width && 0 <= y && y < _height);
	return _data[size_t(_width) * size_t(y) + size_t(x)];
}


template <class T>
T &
ImgVector<T>::at_repeat(const int x, const int y)
{
	size_t x_repeat, y_repeat;

	if (x >= 0) {
		x_repeat = static_cast<size_t>(x % _width);
	} else {
		x_repeat = static_cast<size_t>(_width - (int(std::abs(double(x) + 1.0)) % _width));
	}
	if (y >= 0) {
		y_repeat = static_cast<size_t>(y % _height);
	} else {
		y_repeat = static_cast<size_t>(_height - (int(std::abs(double(y) + 1.0)) % _height));
	}
	return _data[size_t(_width) * y_repeat + x_repeat];
}

template <class T>
const T &
ImgVector<T>::at_repeat(const int x, const int y) const
{
	size_t x_repeat, y_repeat;

	if (x >= 0) {
		x_repeat = static_cast<size_t>(x % _width);
	} else {
		x_repeat = static_cast<size_t>(_width - (int(std::abs(double(x) + 1.0)) % _width));
	}
	if (y >= 0) {
		y_repeat = static_cast<size_t>(y % _height);
	} else {
		y_repeat = static_cast<size_t>(_height - (int(std::abs(double(y) + 1.0)) % _height));
	}
	return _data[size_t(_width) * y_repeat + x_repeat];
}


template <class T>
T &
ImgVector<T>::at_mirror(const int x, const int y)
{
	int x_mirror = x;
	int y_mirror = y;

	if (x_mirror < 0) {
		x_mirror = -x_mirror - 1; // should be set the offset when Mirroring over negative
	}
	if (y_mirror < 0) {
		y_mirror = -y_mirror - 1;
	}
	x_mirror = int(round(_width - 0.5 - std::fabs(_width - 0.5 - (x_mirror % (2 * _width)))));
	y_mirror = int(round(_height - 0.5 - std::fabs(_height - 0.5 - (y_mirror % (2 * _height)))));
	return _data[size_t(_width) * size_t(y_mirror) + size_t(x_mirror)];
}

template <class T>
const T &
ImgVector<T>::at_mirror(const int x, const int y) const
{
	int x_mirror = x;
	int y_mirror = y;

	if (x_mirror < 0) {
		x_mirror = -x_mirror - 1; // should be set the offset when Mirroring over negative
	}
	if (y_mirror < 0) {
		y_mirror = -y_mirror - 1;
	}
	x_mirror = int(round(_width - 0.5 - std::fabs(_width - 0.5 - (x_mirror % (2 * _width)))));
	y_mirror = int(round(_height - 0.5 - std::fabs(_height - 0.5 - (y_mirror % (2 * _height)))));
	return _data[size_t(_width) * size_t(y_mirror) + size_t(x_mirror)];
}


template <class T>
const T
ImgVector<T>::get(const size_t n) const
{
	assert(0 <= n && n < size_t(_width) * size_t(_height));
	return _data[n];
}

template <class T>
const T
ImgVector<T>::get(const int x, const int y) const
{
	assert(0 <= x && x < _width && 0 <= y && y < _height);
	return _data[size_t(_width) * size_t(y) + size_t(x)];
}

template <class T>
const T
ImgVector<T>::get_zeropad(const int x, const int y) const
{
	T zero = T();

	if (x < 0 || _width <= x || y < 0 || _height <= y) {
		return zero;
	} else {
		return _data[size_t(_width) * size_t(y) + size_t(x)];
	}
}

template <class T>
const T
ImgVector<T>::get_repeat(const int x, const int y) const
{
	size_t x_repeat, y_repeat;

	if (x >= 0) {
		x_repeat = static_cast<size_t>(x % _width);
	} else {
		x_repeat = static_cast<size_t>(_width - (int(std::abs(double(x) + 1.0)) % _width));
	}
	if (y >= 0) {
		y_repeat = static_cast<size_t>(y % _height);
	} else {
		y_repeat = static_cast<size_t>(_height - (int(std::abs(double(y) + 1.0)) % _height));
	}
	return _data[size_t(_width) * y_repeat + x_repeat];
}

template <class T>
const T
ImgVector<T>::get_mirror(const int x, const int y) const
{
	int x_mirror = x;
	int y_mirror = y;

	if (x_mirror < 0) {
		x_mirror = -x_mirror - 1; // should be set the offset when Mirroring over negative
	}
	if (y_mirror < 0) {
		y_mirror = -y_mirror - 1;
	}
	x_mirror = int(round(_width - 0.5 - std::fabs(_width - 0.5 - (x_mirror % (2 * _width)))));
	y_mirror = int(round(_height - 0.5 - std::fabs(_height - 0.5 - (y_mirror % (2 * _height)))));
	return _data[size_t(_width) * size_t(y_mirror) + size_t(x_mirror)];
}




// Get continuous function interpolated by bicubic
template <class T>
const T
ImgVector<T>::get_zeropad_cubic(const double& x, const double& y, const double& B, const double& C) const
{
	double bicubic_x[4];
	double bicubic_y[4];
	T value = T();

	if (fabs(x - floor(x)) < DBL_EPSILON
	    && fabs(y - floor(y)) < DBL_EPSILON) {
		value = this->get_zeropad(int(x), int(y));
	} else {
		for (int n = 0; n < 4; n++) {
			bicubic_x[n] = this->cubic(n - 1.0 - (x - floor(x)), B, C);
			bicubic_y[n] = this->cubic(n - 1.0 - (y - floor(y)), B, C);
		}
		for (int m = 0; m < 4; m++) {
			for (int n = 0; n < 4; n++) {
				value += this->get_zeropad(int(floor(x)) + n - 1, int(floor(y)) + m - 1)
				    * bicubic_x[n] * bicubic_y[m];
			}
		}
	}
	return value;
}

template <class T>
const T
ImgVector<T>::get_repeat_cubic(const double& x, const double& y, const double& B, const double& C) const
{
	double bicubic_x[4];
	double bicubic_y[4];
	T value = T();

	if (fabs(x - floor(x)) < DBL_EPSILON
	    && fabs(y - floor(y)) < DBL_EPSILON) {
		value = this->get_zeropad(int(x), int(y));
	} else {
		for (int n = 0; n < 4; n++) {
			bicubic_x[n] = this->cubic(n - 1.0 - (x - floor(x)), B, C);
			bicubic_y[n] = this->cubic(n - 1.0 - (y - floor(y)), B, C);
		}
		for (int m = 0; m < 4; m++) {
			for (int n = 0; n < 4; n++) {
				value += this->get_repeat(int(floor(x)) + n - 1, int(floor(y)) + m - 1)
				    * bicubic_x[n] * bicubic_y[m];
			}
		}
	}
	return value;
}

template <class T>
const T
ImgVector<T>::get_mirror_cubic(const double& x, const double& y, const double& B, const double& C) const
{
	double bicubic_x[4];
	double bicubic_y[4];
	T value = T();

	if (fabs(x - floor(x)) < DBL_EPSILON
	    && fabs(y - floor(y)) < DBL_EPSILON) {
		value = this->get_zeropad(int(x), int(y));
	} else {
		for (int n = 0; n < 4; n++) {
			bicubic_x[n] = this->cubic(n - 1.0 - (x - floor(x)), B, C);
			bicubic_y[n] = this->cubic(n - 1.0 - (y - floor(y)), B, C);
		}
		for (int m = 0; m < 4; m++) {
			for (int n = 0; n < 4; n++) {
				value += this->get_mirror(int(floor(x)) + n - 1, int(floor(y)) + m - 1)
				    * bicubic_x[n] * bicubic_y[m];
			}
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
	for (size_t i = 1; i < size_t(_width) * size_t(_height); i++) {
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
	for (size_t i = 1; i < size_t(_width) * size_t(_height); i++) {
		if (_data[i] > max) {
			max = _data[i];
		}
	}
	return max;
}


template <class T>
const T
ImgVector<T>::min(const int top_left_x, const int top_left_y, const int crop_width, const int crop_height) const
{
	if (_width <= 0 || _height <= 0) {
		throw std::logic_error("T ImgVector<T>::min(void) : vector is empty");
	} else if (crop_width <= 0) {
		throw std::invalid_argument("T ImgVector<T>::min(int, int, int, int) : crop_width <= 0");
	} else if (crop_height <= 0) {
		throw std::invalid_argument("T ImgVector<T>::min(int, int, int, int) : crop_height<= 0");
	}
	T min = _data[size_t(_width) * size_t(top_left_y) + size_t(top_left_x)];
	for (int y = 0; y < crop_height; y++) {
		for (int x = 0; x < crop_width; x++) {
			if (0 <= x + top_left_x && x + top_left_x < _width
			    && 0 <= y + top_left_y && y + top_left_y < _height
			    && _data[size_t(_width) * size_t(y + top_left_y) + size_t(x + top_left_x)] < min) {
				min = _data[size_t(_width) * size_t(y + top_left_y) + size_t(x + top_left_x)];
			}
		}
	}
	return min;
}

template <class T>
const T
ImgVector<T>::max(const int top_left_x, const int top_left_y, const int crop_width, const int crop_height) const
{
	if (_width <= 0 || _height <= 0) {
		throw std::logic_error("T ImgVector<T>::max(void) : vector is empty");
	} else if (crop_width <= 0) {
		throw std::invalid_argument("T ImgVector<T>::max(int, int, int, int) : crop_width <= 0");
	} else if (crop_height <= 0) {
		throw std::invalid_argument("T ImgVector<T>::max(int, int, int, int) : crop_height <= 0");
	}
	T max = _data[size_t(_width) * size_t(top_left_y) + size_t(top_left_x)];
	for (int y = 0; y < crop_height; y++) {
		for (int x = 0; x < crop_width; x++) {
			if (0 <= x + top_left_x && x + top_left_x < _width
			    && 0 <= y + top_left_y && y + top_left_y < _height
			    && _data[size_t(_width) * size_t(y + top_left_y) + size_t(x + top_left_x)] > max) {
				max = _data[size_t(_width) * size_t(y + top_left_y) + size_t(x + top_left_x)];
			}
		}
	}
	return max;
}




// ----- Statistics -----
template <class T>
const T
ImgVector<T>::variance(void) const
{
	double N = double(_width) * double(_height);
	T sum_squared = T();
	T sum = T();

	for (size_t y = 0; y < size_t(_height); y++) {
		for (size_t x = 0; x < size_t(_width); x++) {
			sum += _data[size_t(_width) * y + x];
			sum_squared += _data[size_t(_width) * y + x] * _data[size_t(_width) * y + x];
		}
	}
	return sum_squared / N - sum * sum / (N * N);
}

template <class T>
const T
ImgVector<T>::variance(const int top_left_x, const int top_left_y, const int crop_width, const int crop_height) const
{
	double N = .0;
	T sum_squared = T();
	T sum = T();

	for (int y = top_left_y; y < top_left_y + crop_height; y++) {
		for (int x = top_left_x; x < top_left_x + crop_width; x++) {
			if (0 <= x && x < _width
			    && 0 <= y && y < _height) {
				N += 1.0;
				sum += _data[size_t(_width) * size_t(y) + size_t(x)];
				sum_squared += _data[size_t(_width) * size_t(y) + size_t(x)]
				    * _data[size_t(_width) * size_t(y) + size_t(x)];
			}
		}
	}
	return sum_squared / N - sum * sum / (N * N);
}




// ----- Image Operation -----
template <class T>
ImgVector<T> *
ImgVector<T>::crop(const int top_left_x, const int top_left_y, const int crop_width, const int crop_height) const
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
	catch (const std::bad_alloc& bad) {
		std::cerr << bad.what() << std::endl
		    << "ImgVector<T>* ImgVector<T>::crop(int, int, int, int) : Cannot allocate memory" << std::endl;
		throw;
	}
	// Initialize
	tmp->reset(crop_width, crop_height);
	// Crop
	for (int y = 0; y < crop_width; y++) {
		for (int x = 0; x < crop_height; x++) {
			if (0 <= top_left_y + y && top_left_y + y < _height
			    && 0 <= top_left_x + x && top_left_x + x < _width) {
				tmp->at(x, y) = _data[size_t(_width) * size_t(y + top_left_y) + size_t(x + top_left_x)];
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

	for (size_t i = 1; i < size_t(_width) * size_t(_height); i++) {
		if (_data[i] < min_tmp) {
			min_tmp = _data[i];
		} else if (_data[i] > max_tmp) {
			max_tmp = _data[i];
		}
	}
	// Stretching
	for (size_t i = 0; i < size_t(_width) * size_t(_height); i++) {
		_data[i] = Min + (_data[i] - min_tmp) * Max / (max_tmp - min_tmp);
	}
}




template <class T>
void
ImgVector<T>::resample_zerohold(const int Width, const int Height)
{
	T *resized = nullptr;
	T additive_identity = T();
	double scale_x = .0;
	double scale_y = .0;
	int area_x;
	int area_y;
	T sum;

	if (Width <= 0) {
		throw std::out_of_range("ImgVector<T>::resample_zerohold(const int, const int) : int Width");
	}else if (Height <= 0) {
		throw std::out_of_range("ImgVector<T>::resample_zerohold(const int, const int) : int Height");
	}
	scale_x = double(Width) / _width;
	scale_y = double(Height) / _height;
	try {
		resized = new T[size_t(Width) * size_t(Height)];
	}
	catch (const std::bad_alloc &bad) {
		std::cerr << bad.what() << std::endl
		    << "ImgVector<T>::resample_zerohold(const int, const int) : Cannot allocate memory" << std::endl;
		throw;
	}
	area_x = ceil(double(_width) / Width);
	area_y = ceil(double(_height) / Height);
	for (size_t y = 0; y < size_t(Height); y++) {
		for (size_t x = 0; x < size_t(Width); x++) {
			sum = additive_identity;
			for (int m = 0; m < area_y; m++) {
				for (int n = 0; n < area_x; n++) {
					sum += this->get(int(floor(x / scale_x)) + n, int(floor(y / scale_y)) + m);
				}
			}
			resized[size_t(Width) * y + x] = sum / (area_x * area_y);
		}
	}
	delete[] _data;
	_data = resized;
	_width = Width;
	_height = Height;
}


/*
    bool ImgVector<T>::resize_bicubic(int Width, int Height, double min, double max, T (*Nearest_Integer_Method)(double &d), double A)
    int Width, int Height : width and height of resized image
    double min : minimum value of saturated value
    double max : maximum value of saturated value
    T (*Nearest_Integer_Method)(double &d) : round method (e.g. floor(), round(), etc.)
    A : cubic method's parameter (default A = -0.5 which correspond to Hermite)
*/
template <class T>
void
ImgVector<T>::resample_bicubic(const int Width, const int Height, const bool saturate, const double min, const double max, T (*Nearest_Integer_Method)(double &d), const double B, const double C)
{
	T *resized = nullptr;
	double *conv = nullptr;
	ImgVector<T> Tmp;
	double scale_x, scale_y;
	double scale_conv;
	int L, L_center;
	double dx, dy;
	double sum;

	if (Width <= 0) {
		throw std::out_of_range("ImgVector<T>::resample_bicubic(const int, const int, const double, const double, T (*)(double &d), const double, const double) :int Width");
	} else if (Height <= 0) {
		throw std::out_of_range("ImgVector<T>::resample_bicubic(const int, const int, const double, const double, T (*)(double &d), const double, const double) :int Height");
	}
	scale_x = double(Width) / _width;
	scale_y = double(Height) / _height;
	Tmp.reset(Width, _height);
	// The length of cubic convolution coefficient
	scale_conv = 1.0;
	if (scale_x < 1.0 || scale_y < 1.0) {
		scale_conv = ceil(1.0 / (scale_x < scale_y ? scale_x : scale_y));
	}
	try {
		resized = new T[size_t(Width) * size_t(Height)];
		conv = new double[int(scale_conv * 4)];
	}
	catch (const std::bad_alloc& bad) {
		std::cerr << bad.what() << std::endl
		    << "ImgVector<double>::resample_bicubic(const int, const int, const double, const double, T (*)(double &d), const double, const double) error : Cannot allocate memory" << std::endl;
		delete[] resized;
		delete[] conv;
		throw;
	}
	// Horizontal convolution
	if (scale_x >= 1.0) {
		scale_conv = 1.0;
		L = 4;
		L_center = int(floor((L - 1.0) / 2.0));
	} else {
		scale_conv = 1.0 / scale_x;
		L = 4 * int(ceil(scale_conv));
		L_center = int(floor((L - 1.0) / 2.0));
	}
	for (int x = 0; x < Width; x++) {
		if (scale_x >= 1.0) {
			dx = (x - (scale_x - 1.0) / 2.0) / scale_x;
			for (int n = 0; n < L; n++) {
				conv[n] = ImgVector<T>::cubic(double(n - L_center) - (dx - floor(dx)), B, C);
			}
		} else {
			dx = x / scale_x + (1.0 / scale_x - 1.0) / 2.0;
			for (int n = 0; n < L; n++) {
				conv[n] = ImgVector<T>::cubic((double(n - L_center) - (dx - floor(dx))) * scale_x, B, C) / scale_conv;
			}
		}
		for (int y = 0; y < _height; y++) {
			sum = 0.0;
			for (int n = 0; n < L; n++) {
				sum += conv[n] * this->get_mirror(int(floor(dx)) + n - L_center, y);
			}
			Tmp.at(x, y) = sum;
		}
	}
	// Vertical convolution
	if (scale_y >= 1.0) {
		scale_conv = 1.0;
		L = 4;
		L_center = int(floor((L - 1.0) / 2.0));
	} else {
		scale_conv = 1.0 / scale_y;
		L = 4 * int(ceil(scale_conv));
		L_center = int(floor((L - 1.0) / 2.0));
	}
	for (int y = 0; y < Height; y++) {
		if (scale_y >= 1.0) {
			dy = (y - (scale_y - 1.0) / 2.0) / scale_y;
			for (int m = 0; m < L; m++) {
				conv[m] = ImgVector<T>::cubic(double(m - L_center) - (dy - floor(dy)), B, C);
			}
		} else {
			dy = y / scale_y + (1.0 / scale_y - 1.0) / 2.0;
			for (int m = 0; m < L; m++) {
				conv[m] = ImgVector<T>::cubic((double(m - L_center) - (dy - floor(dy))) / scale_conv, B, C) / scale_conv;
			}
		}
		for (int x = 0; x < Width; x++) {
			sum = 0.0;
			for (int m = 0; m < L; m++) {
				sum += conv[m] * Tmp.get_mirror(x, int(floor(dy)) + m - L_center);
			}
			if (saturate) {
				sum = sum >= min ? sum <= max ? sum : max : min;
			}
			if (Nearest_Integer_Method != nullptr) {
				resized[size_t(Width) * size_t(y) + size_t(x)] = Nearest_Integer_Method(sum);
			} else {
				resized[size_t(Width) * size_t(y) + size_t(x)] = sum;
			}
		}
	}
	delete[] conv;
	delete[] _data;
	_data = resized;
	_width = Width;
	_height = Height;
}


template <class T>
double
ImgVector<T>::cubic(const double x, const double B, const double C) const
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
		for (size_t i = 0 ; i < size_t(_width) * size_t(_height); i++) {
			_data[i] = func(_data[i]);
		}
	}
}




// ----- Arithmetic Operators -----
template <class T>
template <class RT>
ImgVector<T> &
ImgVector<T>::operator+=(const RT& val)
{
	for (size_t i = 0; i < size_t(_width) * size_t(_height); i++) {
		_data[i] += val;
	}
	return *this;
}

template <class T>
template <class RT>
ImgVector<T> &
ImgVector<T>::operator-=(const RT& val)
{
	for (size_t i = 0; i < size_t(_width) * size_t(_height); i++) {
		_data[i] -= val;
	}
	return *this;
}

template <class T>
template <class RT>
ImgVector<T> &
ImgVector<T>::operator*=(const RT& val)
{
	for (size_t i = 0; i < size_t(_width) * size_t(_height); i++) {
		_data[i] *= val;
	}
	return *this;
}

template <class T>
template <class RT>
ImgVector<T> &
ImgVector<T>::operator/=(const RT& val)
{
	for (size_t i = 0; i < size_t(_width) * size_t(_height); i++) {
		_data[i] /= val;
	}
	return *this;
}




template <class T>
template <class RT>
ImgVector<T> &
ImgVector<T>::operator+=(const ImgVector<RT>& rvector)
{
	if (_width != rvector._width
	    || _height != rvector._height) {
		std::cerr << "ImgVector<T>& ImgVector<T>::operator+=(const ImgVector<T>&) : Size of const ImgVector<T>& rvalue is not match" << std::endl;
		throw std::invalid_argument("Size of const ImgVector<T>& rvalue is not match");
	}
	for (size_t i = 0; i < size_t(_width) * size_t(_height); i++) {
		_data[i] += rvector._data[i];
	}
	return *this;
}

template <class T>
template <class RT>
ImgVector<T> &
ImgVector<T>::operator-=(const ImgVector<RT>& rvector)
{
	if (_width != rvector._width
	    || _height != rvector._height) {
		std::cerr << "ImgVector<T>& ImgVector<T>::operator-=(const ImgVector<T>&) : Size of const ImgVector<T>& rvalue is not match" << std::endl;
		throw std::invalid_argument("Size of const ImgVector<T>& rvalue is not match");
	}
	for (size_t i = 0; i < size_t(_width) * size_t(_height); i++) {
		_data[i] -= rvector._data[i];
	}
	return *this;
}

template <class T>
template <class RT>
ImgVector<T> &
ImgVector<T>::operator*=(const ImgVector<RT>& rvector)
{
	if (_width != rvector._width
	    || _height != rvector._height) {
		std::cerr << "ImgVector<T>& ImgVector<T>::operator*=(const ImgVector<T>&) : Size of const ImgVector<T>& rvalue is not match" << std::endl;
		throw std::invalid_argument("Size of const ImgVector<T>& rvalue is not match");
	}
	for (size_t i = 0; i < size_t(_width) * size_t(_height); i++) {
		_data[i] *= rvector._data[i];
	}
	return *this;
}

template <class T>
template <class RT>
ImgVector<T> &
ImgVector<T>::operator/=(const ImgVector<RT>& rvector)
{
	if (_width != rvector._width
	    || _height != rvector._height) {
		std::cerr << "ImgVector<T>& ImgVector<T>::operator/=(const ImgVector<T>&) : Size of const ImgVector<T>& rvalue is not match" << std::endl;
		throw std::invalid_argument("Size of const ImgVector<T>& rvalue is not match");
	}
	for (size_t i = 0; i < size_t(_width) * size_t(_height); i++) {
		_data[i] /= rvector._data[i];
	}
	return *this;
}



// Stream
template <class Type>
std::ostream &
operator<<(std::ostream& os, const ImgVector<Type>& rvalue)
{
#ifdef _CXXABI_H
	char *name = abi::__cxa_demangle(typeid(Type).name(), 0, 0, 0);
	os << "[ size: " << rvalue.width() << "x" << rvalue.height()
	    << ", type: " << name << " ]" << std::endl;
	std::free(name);
#else
	os << "[ size: " << rvalue.width() << "x" << rvalue.height()
	    << ", type: " << typeid(Type).name() << " ]" << std::endl;
#endif
	return os;
}

