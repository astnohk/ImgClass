#include <algorithm>
#include <cmath>
#include <iostream>
#include <new>

#include "CrossCorrelation.h"



CrossCorrelation::CrossCorrelation(void)
{
}

CrossCorrelation::CrossCorrelation(const CrossCorrelation &copy) // copy constructor
{
	_width = copy._width;
	_height = copy._height;
	_img0.copy(copy._img0);
	_img1.copy(copy._img1);
}

CrossCorrelation::CrossCorrelation(ImgStatistics &img0, ImgStatistics &img1)
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
CrossCorrelation::copy(const CrossCorrelation &copy)
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
CrossCorrelation::operator=(const CrossCorrelation &copy)
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
CrossCorrelation::NCC(int x, int y, int window_width, int window_height)
{
	double std_deviation;
	double sum = 0.0;

	std_deviation = _img0.std_deviation(x, y, window_width, window_height) * _img1.std_deviation(x, y, window_width, window_height);
	for (int y = 0; y < window_height; y++) {
		for (int x = 0; x < window_width; x++) {
			sum += (_img0.image(x, y) - _img0.mean(x, y, window_width, window_height))
			    * (_img0.image(x, y) - _img0.mean(x, y, window_width, window_height));
		}
	}
	return sum / std_deviation;
}

double
CrossCorrelation::TruncatedNCC(int x, int y, int window_width, int window_height)
{
	return std::min(1.0, 1.0 - this->NCC(x, y, window_width, window_height));
}

