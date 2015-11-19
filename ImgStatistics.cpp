#include <cmath>
#include <iostream>
#include <new>

#include "ImgStatistics.h"



ImgStatistics::ImgStatistics(void)
{
	_width = 0;
	_height = 0;
	_data = nullptr;
}

ImgStatistics::ImgStatistics(const ImgStatistics &copy)
{
	_width = 0;
	_height = 0;
	_data = nullptr;
	try {
		_data = new double[copy._width * copy._height];
	}
	catch (const std::bad_alloc &bad) {
		std::cerr << "ImgStatistics::ImgStatistics(const ImgStatistics &) error : memory allocation" << std::endl;
		_data = nullptr;
		return;
	}
	_width = copy._width;
	_height = copy._height;
	for (int i = 0; i < _width * _height; i++) {
		_data[i] = copy._data[i];
	}
}

ImgStatistics::ImgStatistics(int W, int H, double *Img)
{
	_width = 0;
	_height = 0;
	_data = nullptr;
	if (W > 0 && H > 0) {
		try {
			_data = new double[W * H]();
		}
		catch (const std::bad_alloc &bad) {
			std::cerr << "ImgStatistics::ImgStatistics(int, int, double *) error : memory allocation" << std::endl;
			return;
		}
		_width = W;
		_height = H;
		if (Img != nullptr) {
			for (int i = 0; i < W * H; i++) {
				_data[i] = Img[i];
			}
		}
	}
}

ImgStatistics::~ImgStatistics(void)
{
	delete _data;
}

void
ImgStatistics::set(int W, int H, double *Img)
{
	_width = 0;
	_height = 0;
	delete[] _data;
	_data = nullptr;
	if (W > 0 && H > 0) {
		try {
			_data = new double[W * H]();
		}
		catch (const std::bad_alloc &bad) {
			std::cerr << "ImgStatistics::set(int, int, double *) error : memory allocation" << std::endl;
			return;
		}
		_width = W;
		_height = H;
		if (Img != nullptr) {
			for (int i = 0; i < W * H; i++) {
				_data[i] = Img[i];
			}
		}
	}
}

ImgStatistics &
ImgStatistics::copy(const ImgStatistics &copy)
{
	if (this != &copy) {
		double *tmp_data = nullptr;
		try {
			tmp_data = new double[copy._width * copy._height];
		}
		catch (const std::bad_alloc &bad) {
			std::cerr << "ImgStatistics::ImgStatistics(const ImgStatistics &) error : memory allocation" << std::endl;
			return *this;
		}
		_width = copy._width;
		_height = copy._height;
		delete[] _data;
		_data = tmp_data;
		for (int i = 0; i < _width * _height; i++) {
			_data[i] = copy._data[i];
		}
	}
	return *this;
}

ImgStatistics &
ImgStatistics::operator=(const ImgStatistics &copy)
{
	if (this != &copy) {
		double *tmp_data = nullptr;
		try {
			tmp_data = new double[copy._width * copy._height];
		}
		catch (const std::bad_alloc &bad) {
			std::cerr << "ImgStatistics::ImgStatistics(const ImgStatistics &) error : memory allocation" << std::endl;
			return *this;
		}
		_width = copy._width;
		_height = copy._height;
		delete[] _data;
		_data = tmp_data;
		for (int i = 0; i < _width * _height; i++) {
			_data[i] = copy._data[i];
		}
	}
	return *this;
}

double &
ImgStatistics::image(int x, int y)
{
	return _data[_width * y + x];
}

int
ImgStatistics::width(void)
{
	return _width;
}

int
ImgStatistics::height(void)
{
	return _height;
}

double
ImgStatistics::mean(void)
{
	double sum = 0.0;

	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			sum += _data[_width * y + x];
		}
	}
	return sum / (_width * _height);
}

double
ImgStatistics::mean(int center_x, int center_y, int window_width, int window_height)
{
	double sum = 0.0;

	for (int y = 0; y < window_height; y++) {
		int y_tmp = center_y + y - (window_height - 1) / 2;
		if (y_tmp < 0 || _height <= y_tmp) {
			continue;
		}
		for (int x = 0; x < window_width; x++) {
			int x_tmp = center_x + x - (window_width - 1) / 2;
			if (x_tmp < 0 || _width <= x_tmp) {
				continue;
			}
			sum += _data[_width * y_tmp + x_tmp];
		}
	}
	return sum / (window_width * window_height);
}

double
ImgStatistics::variance(void)
{
	double sum = 0.0;
	double mu;

	mu = this->mean();
	for (int n = 0; n < _width * _height; n++) {
		sum += (_data[n] - mu) * (_data[n] - mu);
	}
	return sum;
}

double
ImgStatistics::std_deviation(void)
{
	double sum = 0.0;
	double mu;

	mu = this->mean();
	for (int n = 0; n < _width * _height; n++) {
		sum += (_data[n] - mu) * (_data[n] - mu);
	}
	return sqrt(sum);
}

double
ImgStatistics::variance(int center_x, int center_y, int window_width, int window_height)
{
	double sum = 0.0;
	double mu;
	int x_tmp, y_tmp;

	mu = this->mean(center_x, center_y, window_width, window_height);
	for (int y = 0; y < window_height; y++) {
		y_tmp = center_y + y - (window_height - 1) / 2;
		if (y_tmp < 0 || _height <= y_tmp) {
			continue;
		}
		for (int x = 0; x < window_width; x++) {
			x_tmp = center_x + x - (window_width - 1) / 2;
			if (x_tmp < 0 || _width <= x_tmp) {
				continue;
			}
			sum += (_data[_width * y_tmp + x_tmp] - mu) * (_data[_width * y_tmp + x_tmp] - mu);
		}
	}
	return sum;
}

double
ImgStatistics::std_deviation(int center_x, int center_y, int window_width, int window_height)
{
	double sum = 0.0;
	double mu;
	int x_tmp, y_tmp;

	mu = this->mean(center_x, center_y, window_width, window_height);
	for (int y = 0; y < window_height; y++) {
		y_tmp = center_y + y - (window_height - 1) / 2;
		if (y_tmp < 0 || _height <= y_tmp) {
			continue;
		}
		for (int x = 0; x < window_width; x++) {
			x_tmp = center_x + x - (window_width - 1) / 2;
			if (x_tmp < 0 || _width <= x_tmp) {
				continue;
			}
			sum += (_data[_width * y_tmp + x_tmp] - mu) * (_data[_width * y_tmp + x_tmp] - mu);
		}
	}
	return sqrt(sum);
}




Histogram::Histogram(void)
{
	_bins = 0;
	_hist = nullptr;
}

Histogram::Histogram(const Histogram &copy) // copy constructor
{
	_bins = 0;
	_hist = nullptr;
	double *tmp_hist = nullptr;
	try {
		tmp_hist = new double[copy._bins];
	}
	catch (const std::bad_alloc &bad) {
		std::cerr << "Histogram::Histogram(const Histogram &) error : Cannot allocate memory" << std::endl;
		throw;
	}
	for (int i = 0; i < copy._bins; i++) {
		tmp_hist[i] = copy._hist[i];
	}
	_bins = copy._bins;
	_hist = tmp_hist;
}

Histogram::Histogram(int init_bins)
{
	_bins = 0;
	_hist = nullptr;
	double *tmp_hist = nullptr;
	try {
		tmp_hist = new double[init_bins];
	}
	catch (const std::bad_alloc &bad) {
		std::cerr << "Histogram::copy(const Histogram &) error : Cannot allocate memory" << std::endl;
		throw;
	}
	for (int i = 0; i < init_bins; i++) {
		tmp_hist[i] = .0;
	}
	_bins = init_bins;
	_hist = tmp_hist;
}

Histogram::~Histogram(void)
{
	delete[] _hist;
}

void
Histogram::free(void)
{
	_bins = 0;
	delete[] _hist;
	_hist = nullptr;
}

Histogram &
Histogram::copy(const Histogram &copy)
{
	if (this != &copy) {
		double *tmp_hist = nullptr;
		try {
			tmp_hist = new double[copy._bins];
		}
		catch (const std::bad_alloc &bad) {
			std::cerr << "Histogram::copy(const Histogram &) error : Cannot allocate memory" << std::endl;
			throw;
			return *this;
		}
		for (int i = 0; i < copy._bins; i++) {
			tmp_hist[i] = copy._hist[i];
		}
		delete[] _hist;
		_hist = tmp_hist;
		_bins = copy._bins;
	}
	return *this;
}

Histogram &
Histogram::reset(int init_bins)
{
	if (init_bins >= 0) {
		double *tmp_hist = nullptr;
		try {
			tmp_hist = new double[init_bins];
		}
		catch (const std::bad_alloc &bad) {
			std::cerr << "Histogram::reset(int) error : Cannot allocate memory" << std::endl;
			throw;
			return *this;
		}
		for (int i = 0; i < init_bins; i++) {
			tmp_hist[i] = .0;
		}
		delete[] _hist;
		_hist = tmp_hist;
		_bins = init_bins;
	}
	return *this;
}

const double *
Histogram::data(void) const
{
	return _hist;
}

int
Histogram::bins(void) const
{
	return _bins;
}

double
Histogram::get(int bin) const
{
	if (bin < 0 || _bins <= bin) {
		return 0;
	}
	return _hist[bin];
}

bool
Histogram::add(int bin, double val)
{
	if (bin < 0 || _bins <= bin) {
		return false;
	}
	_hist[bin] += val;
	return true;
}

