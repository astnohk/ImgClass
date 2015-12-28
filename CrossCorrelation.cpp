#include <algorithm>
#include <cmath>
#include <iostream>
#include <new>

#include "CrossCorrelation.h"



CrossCorrelation::CrossCorrelation(void)
{
}

CrossCorrelation::CrossCorrelation(const CrossCorrelation& copy) // copy constructor
{
	_width = copy._width;
	_height = copy._height;
	_img0.copy(copy._img0);
	_img1.copy(copy._img1);
}

CrossCorrelation::CrossCorrelation(const ImgStatistics& img0, const ImgStatistics& img1)
{
	if (img0.width() != img1.width()
	    || img0.height() != img1.height()) {
		return;
	} else {
		_img0.copy(img0);
		_img1.copy(img1);
		_width = img0.width();
		_height = img0.height();
	}
}


CrossCorrelation::~CrossCorrelation(void)
{
}


CrossCorrelation &
CrossCorrelation::copy(const CrossCorrelation& copy)
{
	if (this != &copy) {
		_img0.copy(copy._img0);
		_img1.copy(copy._img1);
		_width = copy._width;
		_height = copy._height;
	}
	return *this;
}

CrossCorrelation &
CrossCorrelation::operator=(const CrossCorrelation& copy)
{
	if (this != &copy) {
		_img0.copy(copy._img0);
		_img1.copy(copy._img1);
		_width = copy._width;
		_height = copy._height;
	}
	return *this;
}




int
CrossCorrelation::width(void)
{
	return _width;
}

int
CrossCorrelation::height(void)
{
	return _height;
}




double
CrossCorrelation::NCC(const int x, const int y, const int window_width, const int window_height)
{
	double std_deviation;
	double sum = 0.0;

	std_deviation = _img0.std_deviation(x, y, window_width, window_height) * _img1.std_deviation(x, y, window_width, window_height);
	for (int m = 0; m < window_height; m++) {
		for (int n = 0; n < window_width; n++) {
			sum += (_img0.image(n, m) - _img0.mean(n, m, window_width, window_height))
			    * (_img0.image(n, m) - _img0.mean(n, m, window_width, window_height));
		}
	}
	return sum / std_deviation;
}

double
CrossCorrelation::TruncatedNCC(const int x, const int y, const int window_width, const int window_height)
{
	return std::min(1.0, 1.0 - this->NCC(x, y, window_width, window_height));
}

