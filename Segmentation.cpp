#include "Segmentation.h"
#include "Color.h"




namespace ImgClass {
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
	Segmentation<ImgClass::RGB>::normalized_distance(const ImgClass::RGB& lvalue, const ImgClass::RGB& rvalue)
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

	template <>
	double
	Segmentation<ImgClass::Lab>::normalized_distance(const ImgClass::Lab& lvalue, const ImgClass::Lab& rvalue)
	{
		return sqrt(
		    SQUARE(lvalue.L - rvalue.L)
		    + SQUARE(lvalue.a - rvalue.a)
		    + SQUARE(lvalue.b - rvalue.b)) / 100.0;
	}




	/*
	 * std::vector<double> kernel has kernel radius for each dimensions.
	 * The values it needs are below:
	 *	_kernel_spatial : the spatial radius of mean shift kernel
	 *	_kernel_intensity : the norm threshold of mean shift kernel in L*a*b* space
	 */
	template <> // Specialized for ImgClass::RGB<double>
	const Segmentation<ImgClass::RGB>::tuple
	Segmentation<ImgClass::RGB>::MeanShift(const int x, const int y, std::vector<VECTOR_2D<int> >& pel_list, int Iter_Max)
	{
		const double radius_spatial_squared = SQUARE(_kernel_spatial);
		const double radius_intensity_squared = SQUARE(_kernel_intensity);
		const double displacement_color_min = SQUARE(0.0001);
		const double displacement_spatial_min = SQUARE(0.01);
		Segmentation<ImgClass::RGB>::tuple tuple;
		// Initialize
		tuple.spatial = VECTOR_2D<double>(static_cast<double>(x), static_cast<double>(y));
		tuple.color = _image.get(x, y);
		// Iterate until it converge
		ImgClass::RGB sum_diff(0.0, 0.0, 0.0);
		VECTOR_2D<double> sum_d(0.0, 0.0);
		VECTOR_2D<int> r(0, 0);
		VECTOR_2D<double> d(0.0, 0.0);
		for (int i = 0; i < Iter_Max; i++) {
			double N = 0.0;
			sum_diff.R = 0; sum_diff.G = 0; sum_diff.B = 0;
			sum_d.x = 0.0; sum_d.y = 0.0;
			for (size_t n = 0; n < pel_list.size(); n++) {
				r.x = static_cast<int>(round(tuple.spatial.x) + pel_list[n].x);
				r.x = static_cast<int>(round(tuple.spatial.y) + pel_list[n].y);
				if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height) {
					ImgClass::RGB diff(_image.get(r.x, r.y) - tuple.color);
					d.x = r.x - tuple.spatial.x;
					d.y = r.y - tuple.spatial.y;
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
			tuple.color += sum_diff / N;
			VECTOR_2D<double> displacement(sum_d.x / N, sum_d.y / N);
			tuple.spatial += displacement;
			if (norm_squared(sum_diff / N) < displacement_color_min && norm_squared(displacement) < displacement_spatial_min) {
				break;
			}
		}
		return tuple;
	}

	/*
	 * std::vector<double> kernel has kernel radius for each dimensions.
	 * The values it needs are below:
	 *	_kernel_spatial : the spatial radius of mean shift kernel
	 *	_kernel_intensity : the norm threshold of mean shift kernel in L*a*b* space
	 */
	template <> // Specialized for ImgClass::Lab
	const Segmentation<ImgClass::Lab>::tuple
	Segmentation<ImgClass::Lab>::MeanShift(const int x, const int y, std::vector<VECTOR_2D<int> >& pel_list, int Iter_Max)
	{
		const double radius_spatial_squared = SQUARE(_kernel_spatial);
		const double radius_intensity_squared = SQUARE(100.0 * _kernel_intensity);
		const double displacement_color_min = SQUARE(0.0001);
		const double displacement_spatial_min = SQUARE(0.01);
		Segmentation<ImgClass::Lab>::tuple tuple;

		// Initialize
		tuple.spatial = VECTOR_2D<double>(static_cast<double>(x), static_cast<double>(y));
		tuple.color = _image.get(x, y);
		// Iterate until it converge
		ImgClass::Lab sum_diff(0.0, 0.0, 0.0);
		VECTOR_2D<double> sum_d(0.0, 0.0);
		VECTOR_2D<int> r(0, 0);
		VECTOR_2D<double> d(0.0, 0.0);
		for (int i = 0; i < Iter_Max; i++) {
			double N = 0.0;
			sum_diff.L = 0.0; sum_diff.a = 0.0; sum_diff.b = 0.0;
			sum_d.x = 0.0; sum_d.y = 0.0;
			for (size_t n = 0; n < pel_list.size(); n++) {
				r.x = static_cast<int>(round(tuple.spatial.x) + pel_list[n].x);
				r.y = static_cast<int>(round(tuple.spatial.y) + pel_list[n].y);
				if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height) {
					ImgClass::Lab diff(_image.get(r.x, r.y) - tuple.color);
					//diff.L *= 0.8; // Difference of Lighting is not so important in segmentation
					d.x = r.x - tuple.spatial.x;
					d.y = r.y - tuple.spatial.y;
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
			tuple.color += sum_diff / N;
			VECTOR_2D<double> displacement(sum_d.x / N, sum_d.y / N);
			tuple.spatial += displacement;
			if (norm_squared(sum_diff / N) < displacement_color_min && norm_squared(displacement) < displacement_spatial_min) {
				break;
			}
		}
		return tuple;
	}
}




// For polymorphism
double
inner_prod(const double& lvalue, const double& rvalue)
{
	return lvalue * rvalue;
}


double norm(const double& value)
{
	return sqrt(value * value);
}

double norm_squared(const double& value)
{
	return value * value;
}

