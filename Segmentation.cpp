#include "Segmentation.h"




template <>
double
Segmentation<ImgClass::RGB>::distance(const ImgClass::RGB& lvalue, const ImgClass::RGB& rvalue)
{
	return sqrt(
	    SQUARE(lvalue.R - rvalue.R)
	    + SQUARE(lvalue.G - rvalue.G)
	    + SQUARE(lvalue.B - rvalue.B));
}

template <>
double
Segmentation<ImgClass::Lab>::distance(const ImgClass::Lab& lvalue, const ImgClass::Lab& rvalue)
{
	return sqrt(
	    SQUARE(lvalue.L - rvalue.L)
	    + SQUARE(lvalue.a - rvalue.a)
	    + SQUARE(lvalue.b - rvalue.b));
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
	const double radius_spatial_squared = SQUARE(_kernel_spatial);
	const double radius_intensity_squared = SQUARE(_kernel_intensity);
	const double displacement_min = SQUARE(0.01);
	// Initialize
	VECTOR_2D<double> u(static_cast<double>(x), static_cast<double>(y));
	ImgClass::RGB center(_image.get(x, y));
	// Iterate until it converge
	for (int i = 0; i < Iter_Max; i++) {
		double N = 0.0;
		ImgClass::RGB sum_diff(0.0, 0.0, 0.0);
		VECTOR_2D<double> sum_d(0.0, 0.0);
		for (unsigned int n = 0; n < pel_list.size(); n++) {
			VECTOR_2D<int> r(
			    static_cast<int>(round(u.x) + pel_list[n].x),
			    static_cast<int>(round(u.y) + pel_list[n].y));
			if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height) {
				ImgClass::RGB diff(_image.get(r.x, r.y) - center);
				VECTOR_2D<double> d(r.x - u.x, r.y - u.y);
				double ratio_intensity = norm_squared(diff) / radius_intensity_squared;
				double ratio_spatial = norm_squared(d) / radius_spatial_squared;
				if (ratio_intensity <= 1.0 && ratio_spatial <= 1.0) {
					double coeff = 1.0 - (ratio_intensity * ratio_spatial);
					N += coeff;
					sum_diff += diff * coeff;
					sum_d += d * coeff;
				}
			}
		}
		center += sum_diff / N;
		VECTOR_2D<double> displacement(sum_d.x / N, sum_d.y / N);
		u += displacement;
		if (norm_squared(displacement) < displacement_min) {
			break;
		}
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
	const double radius_spatial_squared = SQUARE(_kernel_spatial);
	const double radius_intensity_squared = SQUARE(100.0 * _kernel_intensity);
	const double displacement_min = SQUARE(0.01);

	// Initialize
	VECTOR_2D<double> u(static_cast<double>(x), static_cast<double>(y));
	ImgClass::Lab center = _image.get(x, y);
	// Iterate until it converge
	for (int i = 0; i < Iter_Max; i++) {
		double N = 0.0;
		ImgClass::Lab sum_diff(0.0, 0.0, 0.0);
		VECTOR_2D<double> sum_d(0.0, 0.0);
		for (unsigned int n = 0; n < pel_list.size(); n++) {
			VECTOR_2D<int> r(
			    static_cast<int>(round(u.x) + pel_list[n].x),
			    static_cast<int>(round(u.y) + pel_list[n].y));
			if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height) {
				ImgClass::Lab diff(_image.get(r.x, r.y) - center);
				diff.L /= 4.0; // Difference of Lighting is not so important in segmentation
				VECTOR_2D<double> d(r.x - u.x, r.y - u.y);
				double ratio_intensity = norm_squared(diff) / radius_intensity_squared;
				double ratio_spatial = norm_squared(d) / radius_spatial_squared;
				if (ratio_intensity <= 1.0 && ratio_spatial <= 1.0) {
					double coeff = 1.0 - (ratio_intensity * ratio_spatial);
					N += coeff;
					sum_diff += diff * coeff;
					sum_d += d * coeff;
				}
			}
		}
		center += sum_diff / N;
		VECTOR_2D<double> displacement(sum_d.x / N, sum_d.y / N);
		u += displacement;
		if (norm_squared(displacement) < displacement_min) {
			break;
		}
	}
	return u;
}

