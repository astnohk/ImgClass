#include "BlockMatching.h"
#include "Lab.h"
#include "RGB.h"




// ----- Global function -----
double
norm_squared(const double& value)
{
	return value * value;
}

double
norm(const double& value)
{
	return fabs(value);
}

double
inner_prod(const double& lvalue, const double& rvalue)
{
	return lvalue * rvalue;
}




// ----- Specialize -----
template <>
void
BlockMatching<ImgClass::RGB>::image_normalizer(void)
{
	// Previous
	double max_int = .0;
	for (size_t i = 0; i < _image_prev.size(); i++) {
		if (norm(_image_prev[i]) > max_int) {
			max_int = norm(_image_prev[i]);
		}
	}
	if (max_int > 1.0) {
		_image_prev /= max_int;
	}
	// Current
	max_int = .0;
	for (size_t i = 0; i < _image_current.size(); i++) {
		if (norm(_image_current[i]) > max_int) {
			max_int = norm(_image_current[i]);
		}
	}
	if (max_int > 1.0) {
		_image_current /= max_int;
	}
	// Next
	max_int = .0;
	for (size_t i = 0; i < _image_next.size(); i++) {
		if (norm(_image_next[i]) > max_int) {
			max_int = norm(_image_next[i]);
		}
	}
	if (max_int > 1.0) {
		_image_next /= max_int;
	}
}

template <>
void
BlockMatching<ImgClass::Lab>::image_normalizer(void)
{
	// Previous
	double max_int = .0;
	for (size_t i = 0; i < _image_prev.size(); i++) {
		if (norm(_image_prev[i]) > max_int) {
			max_int = norm(_image_prev[i]);
		}
	}
	if (max_int > 1.0) {
		_image_prev /= max_int;
	}
	// Current
	max_int = .0;
	for (size_t i = 0; i < _image_current.size(); i++) {
		if (norm(_image_current[i]) > max_int) {
			max_int = norm(_image_current[i]);
		}
	}
	if (max_int > 1.0) {
		_image_current /= max_int;
	}
	// Next
	max_int = .0;
	for (size_t i = 0; i < _image_next.size(); i++) {
		if (norm(_image_next[i]) > max_int) {
			max_int = norm(_image_next[i]);
		}
	}
	if (max_int > 1.0) {
		_image_next /= max_int;
	}
}




// ---------- Correlation functions ----------
// ----- region nearest intensity -----
template <>
double
BlockMatching<ImgClass::RGB>::MAD_region_nearest_intensity(const int x_diff_prev, const int y_diff_prev, const std::list<VECTOR_2D<int> >& region)
{
	double N = .0;
	double sad = .0;

	std::list<VECTOR_2D<int> >::const_iterator ite = region.begin();
	ImgClass::RGB color(_color_quantized_current.get(ite->x, ite->y));
	for (;
	    ite != region.end();
	    ++ite) {
		VECTOR_2D<int> r(ite->x + x_diff_prev, ite->y + y_diff_prev);
		double coeff = norm(color - _color_quantized_prev.get_zeropad(r.x, r.y));
		N += coeff;
		sad += coeff * norm(
		    _image_current.get_zeropad(ite->x, ite->y)
		    - _image_prev.get_zeropad(r.x, r.y));
	}
	return sad / N;
}

template <>
double
BlockMatching<ImgClass::Lab>::MAD_region_nearest_intensity(const int x_diff_prev, const int y_diff_prev, const std::list<VECTOR_2D<int> >& region)
{
	double N = .0;
	double sad = .0;
	
	std::list<VECTOR_2D<int> >::const_iterator ite = region.begin();
	ImgClass::Lab color(_color_quantized_current.get(ite->x, ite->y));
	for (;
	    ite != region.end();
	    ++ite) {
		VECTOR_2D<int> r(ite->x + x_diff_prev, ite->y + y_diff_prev);
		double coeff = norm(color - _color_quantized_prev.get_zeropad(r.x, r.y));
		N += coeff;
		sad += coeff * norm(
		    _image_current.get_zeropad(ite->x, ite->y)
		    - _image_prev.get_zeropad(r.x, r.y));
	}
	return sad / N;
}


template <>
double
BlockMatching<ImgClass::RGB>::ZNCC_region_nearest_intensity(const int x_diff_prev, const int y_diff_prev, const std::list<VECTOR_2D<int> >& region)
{
	double N = .0;
	ImgClass::RGB sum_prev(.0, .0, .0);
	ImgClass::RGB sum_current(.0, .0, .0);
	double sum_sq_prev = .0;
	double sum_sq_current = .0;
	double sum_sq_prev_current = .0;

	std::list<VECTOR_2D<int> >::const_iterator ite = region.begin();
	ImgClass::RGB color(_color_quantized_current.get(ite->x, ite->y));
	for (;
	    ite != region.end();
	    ++ite) {
		VECTOR_2D<int> r(ite->x + x_diff_prev, ite->y + y_diff_prev);
		double coeff = 1.0 - norm(color - _color_quantized_prev.get_zeropad(r.x, r.y));
		N += coeff;
		ImgClass::RGB color_prev(_image_prev.get_zeropad(r.x, r.y));
		ImgClass::RGB color_current(_image_current.get_zeropad(ite->x, ite->y));
		// Previous frame
		sum_prev += coeff * color_prev;
		sum_sq_prev += coeff * inner_prod(color_prev, color_prev);
		// Next frame
		sum_current += coeff * color_current;
		sum_sq_current += coeff * inner_prod(color_current, color_current);
		// Co-frame
		sum_sq_prev_current += coeff * inner_prod(color_prev, color_current);
	}
	// Calculate Covariance
	return (N * sum_sq_prev_current - inner_prod(sum_prev, sum_current)) /
	    (sqrt((N * sum_sq_prev - inner_prod(sum_prev, sum_prev))
	    * (N * sum_sq_current - inner_prod(sum_current, sum_current)))
	    + DBL_EPSILON);
}

template <>
double
BlockMatching<ImgClass::Lab>::ZNCC_region_nearest_intensity(const int x_diff_prev, const int y_diff_prev, const std::list<VECTOR_2D<int> >& region)
{
	double N = .0;
	ImgClass::Lab sum_prev(.0, .0, .0);
	ImgClass::Lab sum_current(.0, .0, .0);
	double sum_sq_prev = .0;
	double sum_sq_current = .0;
	double sum_sq_prev_current = .0;

	std::list<VECTOR_2D<int> >::const_iterator ite = region.begin();
	ImgClass::Lab quantized_color = _color_quantized_current.get(ite->x, ite->y);
	for (;
	    ite != region.end();
	    ++ite) {
		VECTOR_2D<int> r(ite->x + x_diff_prev, ite->y + y_diff_prev);
		double coeff = 1.0 - norm(quantized_color - _color_quantized_prev.get_zeropad(r.x, r.y));
		N += coeff;
		ImgClass::Lab color_prev(_image_prev.get_zeropad(r.x, r.y));
		ImgClass::Lab color_current(_image_current.get_zeropad(ite->x, ite->y));
		// Previous frame
		sum_prev += coeff * color_prev;
		sum_sq_prev += coeff * inner_prod(color_prev, color_prev);
		// Next frame
		sum_current += coeff * color_current;
		sum_sq_current += coeff * inner_prod(color_current, color_current);
		// Co-frame
		sum_sq_prev_current += coeff * inner_prod(color_prev, color_current);
	}
	// Calculate Covariance
	return (N * sum_sq_prev_current - inner_prod(sum_prev, sum_current)) /
	    (sqrt((N * sum_sq_prev - inner_prod(sum_prev, sum_prev))
	    * (N * sum_sq_current - inner_prod(sum_current, sum_current)))
	    + DBL_EPSILON);
}

