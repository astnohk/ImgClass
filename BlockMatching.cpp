#include "BlockMatching.h"




// ----- Miscellaneous -----
template <>
void
BlockMatching<ImgClass::RGB>::image_normalizer(void)
{
	// Previous
	double max_int = .0;
	for (int i = 0; i < _width * _height; i++) {
		if (norm(_image_prev[i]) > max_int) {
			max_int = norm(_image_prev[i]);
		}
	}
	if (max_int > 1.0) {
		_image_prev /= max_int;
	}
	// Current
	max_int = .0;
	for (int i = 0; i < _width * _height; i++) {
		if (norm(_image_current[i]) > max_int) {
			max_int = norm(_image_current[i]);
		}
	}
	if (max_int > 1.0) {
		_image_current /= max_int;
	}
	// Next
	max_int = .0;
	for (int i = 0; i < _image_next.size(); i++) {
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
	for (int i = 0; i < _width * _height; i++) {
		if (norm(_image_prev[i]) > max_int) {
			max_int = norm(_image_prev[i]);
		}
	}
	if (max_int > 1.0) {
		_image_prev /= max_int;
	}
	// Current
	max_int = .0;
	for (int i = 0; i < _width * _height; i++) {
		if (norm(_image_current[i]) > max_int) {
			max_int = norm(_image_current[i]);
		}
	}
	if (max_int > 1.0) {
		_image_current /= max_int;
	}
	// Next
	max_int = .0;
	for (int i = 0; i < _image_next.size(); i++) {
		if (norm(_image_next[i]) > max_int) {
			max_int = norm(_image_next[i]);
		}
	}
	if (max_int > 1.0) {
		_image_next /= max_int;
	}
}




template <>
ImgVector<VECTOR_2D<double> > *
BlockMatching<ImgClass::Lab>::grad_image(const ImgVector<ImgClass::Lab>& image, const int top_left_x, const int top_left_y, const int crop_width, const int crop_height)
{
	ImgVector<VECTOR_2D<double> >* gradients = new ImgVector<VECTOR_2D<double> >(crop_width, crop_height);

	for (int y = 0; y < crop_height; y++) {
		for (int x = 0; x < crop_width; x++) {
			gradients->at(x, y).x =
			    image.get_mirror(top_left_x + x + 1, top_left_y + y).L
			    - image.get_mirror(top_left_x + x, top_left_y + y).L;
			gradients->at(x, y).y =
			    image.get_mirror(top_left_x + x, top_left_y + y + 1).L
			    - image.get_mirror(top_left_x + x, top_left_y + y).L;
		}
	}
	return gradients;
}




// ----- Correlation functions -----
template <>
double
BlockMatching<ImgClass::RGB>::MAD(const int x_l, const int y_l, const int x_r, const int y_r, const int block_width, const int block_height, const ImgVector<ImgClass::RGB>& limage, const ImgVector<ImgClass::RGB>& rimage)
{
	double sad = .0;

	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			sad = sad + norm(limage.get_zeropad(x_l+ x, y_l+ y) - rimage.get_zeropad(x_r + x, y_r + y));
		}
	}
	return sad / double(block_width * block_height);
}

template <>
double
BlockMatching<ImgClass::Lab>::MAD(const int x_l, const int y_l, const int x_r, const int y_r, const int block_width, const int block_height, const ImgVector<ImgClass::Lab>& limage, const ImgVector<ImgClass::Lab>& rimage)
{
	double sad = .0;

	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			sad = sad + norm(limage.get_zeropad(x_l+ x, y_l+ y) - rimage.get_zeropad(x_r + x, y_r + y));
		}
	}
	return sad / double(block_width * block_height);
}




template <>
double
BlockMatching<ImgClass::RGB>::MAD(const int x_prev, const int y_prev, const int x_current, const int y_current, const int block_width, const int block_height)
{
	double sad = .0;

	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			ImgClass::RGB color_prev(_image_prev.get_zeropad(x_prev + x, y_prev + y));
			ImgClass::RGB color_current(_image_current.get_zeropad(x_current + x, y_current + y));
			sad = sad + norm(color_current - color_prev);
		}
	}
	return sad / double(block_width * block_height);
}

template <>
double
BlockMatching<ImgClass::Lab>::MAD(const int x_prev, const int y_prev, const int x_current, const int y_current, const int block_width, const int block_height)
{
	double sad = .0;

	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			ImgClass::Lab color_prev(_image_prev.get_zeropad(x_prev + x, y_prev + y));
			ImgClass::Lab color_current(_image_current.get_zeropad(x_current + x, y_current + y));
			sad = sad + norm(color_current - color_prev);
		}
	}
	return sad / double(block_width * block_height);
}




template <>
double
BlockMatching<ImgClass::RGB>::ZNCC(const int x_prev, const int y_prev, const int x_current, const int y_current, const int block_width, const int block_height)
{
	double N = block_width * block_height;
	ImgClass::RGB sum_prev(.0, .0, .0);
	ImgClass::RGB sum_current(.0, .0, .0);
	double sum_sq_prev = .0;
	double sum_sq_current = .0;
	double sum_sq_prev_current = .0;

	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			sum_prev += _image_prev.get_zeropad(x_prev + x, y_prev + y);
			sum_current += _image_current.get_zeropad(x_current + x, y_current + y);
			sum_sq_prev += inner_prod(_image_prev.get_zeropad(x_prev + x, y_prev + y), _image_prev.get_zeropad(x_prev + x, y_prev + y));
			sum_sq_current += inner_prod(_image_current.get_zeropad(x_current + x, y_current + y), _image_current.get_zeropad(x_current + x, y_current + y));
			sum_sq_prev_current += inner_prod(_image_prev.get_zeropad(x_prev + x, y_prev + y), _image_current.get_zeropad(x_current + x, y_current + y));
		}
	}
	return (N * sum_sq_prev_current - inner_prod(sum_prev, sum_current)) / (sqrt(N * sum_sq_prev - inner_prod(sum_prev, sum_prev)) * sqrt(N * sum_sq_current - inner_prod(sum_current, sum_current)) + 1.0E-10);
}

template <>
double
BlockMatching<ImgClass::Lab>::ZNCC(const int x_prev, const int y_prev, const int x_current, const int y_current, const int block_width, const int block_height)
{
	double N = block_width * block_height;
	ImgClass::Lab sum_prev(.0, .0, .0);
	ImgClass::Lab sum_current(.0, .0, .0);
	double sum_sq_prev = .0;
	double sum_sq_current = .0;
	double sum_sq_prev_current = .0;

	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			sum_prev += _image_prev.get_zeropad(x_prev + x, y_prev + y);
			sum_current += _image_current.get_zeropad(x_current + x, y_current + y);
			sum_sq_prev += inner_prod(_image_prev.get_zeropad(x_prev + x, y_prev + y), _image_prev.get_zeropad(x_prev + x, y_prev + y));
			sum_sq_current += inner_prod(_image_current.get_zeropad(x_current + x, y_current + y), _image_current.get_zeropad(x_current + x, y_current + y));
			sum_sq_prev_current += inner_prod(_image_prev.get_zeropad(x_prev + x, y_prev + y), _image_current.get_zeropad(x_current + x, y_current + y));
		}
	}
	return (N * sum_sq_prev_current - inner_prod(sum_prev, sum_current)) / (sqrt(N * sum_sq_prev - inner_prod(sum_prev, sum_prev)) * sqrt(N * sum_sq_current - inner_prod(sum_current, sum_current)) + 1.0E-10);
}




// ---------- region ----------
template <>
double
BlockMatching<ImgClass::RGB>::MAD_region(const ImgVector<ImgClass::RGB>& reference, const ImgVector<ImgClass::RGB>& current, const int x_diff, const int y_diff, const std::list<VECTOR_2D<int> >& region_current)
{
	double N = .0;
	double sad = .0;

	for (std::list<VECTOR_2D<int> >::const_iterator ite = region_current.begin();
	    ite != region_current.end();
	    ++ite) {
		N += 1.0;
		ImgClass::RGB color_reference(reference.get_zeropad(ite->x + x_diff, ite->y + y_diff));
		ImgClass::RGB color_current(current.get_zeropad(ite->x, ite->y));
		sad += norm(color_current - color_reference);
	}
	return sad / N;
}

template <>
double
BlockMatching<ImgClass::Lab>::MAD_region(const ImgVector<ImgClass::Lab>& reference, const ImgVector<ImgClass::Lab>& current, const int x_diff, const int y_diff, const std::list<VECTOR_2D<int> >& region_current)
{
	double N = .0;
	double sad = .0;

	for (std::list<VECTOR_2D<int> >::const_iterator ite = region_current.begin();
	    ite != region_current.end();
	    ++ite) {
		N += 1.0;
		ImgClass::Lab color_reference(reference.get_zeropad(ite->x + x_diff, ite->y + y_diff));
		ImgClass::Lab color_current(current.get_zeropad(ite->x, ite->y));
		sad += norm(color_current - color_reference);
	}
	return sad / N;
}


template <>
double
BlockMatching<ImgClass::RGB>::ZNCC_region(const ImgVector<ImgClass::RGB>& reference, const ImgVector<ImgClass::RGB>& current, const int x_diff, const int y_diff, const std::list<VECTOR_2D<int> >& region_current)
{
	double N = .0;
	ImgClass::RGB sum_reference(.0, .0, .0);
	ImgClass::RGB sum_current(.0, .0, .0);
	double sum_sq_reference = .0;
	double sum_sq_current = .0;
	double sum_sq_reference_current = .0;

	for (std::list<VECTOR_2D<int> >::const_iterator ite = region_current.begin();
	    ite != region_current.end();
	    ++ite) {
		VECTOR_2D<int> r(ite->x + x_diff, ite->y + y_diff);

		N += 1.0;
		// Previous frame
		sum_reference += reference.get_zeropad(r.x, r.y);
		sum_sq_reference += inner_prod(
		    reference.get_zeropad(r.x, r.y)
		    , reference.get_zeropad(r.x, r.y));
		// Next frame
		sum_current += current.get_zeropad(ite->x, ite->y);
		sum_sq_current += inner_prod(
		    current.get_zeropad(ite->x, ite->y)
		    , current.get_zeropad(ite->x, ite->y));
		// Co-frame
		sum_sq_reference_current += inner_prod(
		    reference.get_zeropad(r.x, r.y)
		    , current.get_zeropad(ite->x, ite->y));
	}
	// Calculate Covariance
	return (N * sum_sq_reference_current - inner_prod(sum_reference, sum_current))
	    / (sqrt(N * sum_sq_reference - inner_prod(sum_reference, sum_reference)) * sqrt(N * sum_sq_current - inner_prod(sum_current, sum_current)) + 1.0E-10);
}

template <>
double
BlockMatching<ImgClass::Lab>::ZNCC_region(const ImgVector<ImgClass::Lab>& reference, const ImgVector<ImgClass::Lab>& current, const int x_diff, const int y_diff, const std::list<VECTOR_2D<int> >& region_current)
{
	double N = .0;
	ImgClass::Lab sum_reference(.0, .0, .0);
	ImgClass::Lab sum_current(.0, .0, .0);
	double sum_sq_reference = .0;
	double sum_sq_current = .0;
	double sum_sq_reference_current = .0;

	for (std::list<VECTOR_2D<int> >::const_iterator ite = region_current.begin();
	    ite != region_current.end();
	    ++ite) {
		VECTOR_2D<int> r(ite->x + x_diff, ite->y + y_diff);

		N += 1.0;
		// Previous frame
		sum_reference += reference.get_zeropad(r.x, r.y);
		sum_sq_reference += inner_prod(
		    reference.get_zeropad(r.x, r.y)
		    , reference.get_zeropad(r.x, r.y));
		// Next frame
		sum_current += current.get_zeropad(ite->x, ite->y);
		sum_sq_current += inner_prod(
		    current.get_zeropad(ite->x, ite->y)
		    , current.get_zeropad(ite->x, ite->y));
		// Co-frame
		sum_sq_reference_current += inner_prod(
		    reference.get_zeropad(r.x, r.y)
		    , current.get_zeropad(ite->x, ite->y));
	}
	// Calculate Covariance
	return (N * sum_sq_reference_current - inner_prod(sum_reference, sum_current))
	    / (sqrt(N * sum_sq_reference - inner_prod(sum_reference, sum_reference)) * sqrt(N * sum_sq_current - inner_prod(sum_current, sum_current)) + 1.0E-10);
}




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
	return (N * sum_sq_prev_current - inner_prod(sum_prev, sum_current))
	    / (sqrt(N * sum_sq_prev - inner_prod(sum_prev, sum_prev)) * sqrt(N * sum_sq_current - inner_prod(sum_current, sum_current)) + 1.0E-10);
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
	return (N * sum_sq_prev_current - inner_prod(sum_prev, sum_current))
	    / (sqrt(N * sum_sq_prev - inner_prod(sum_prev, sum_prev)) * sqrt(N * sum_sq_current - inner_prod(sum_current, sum_current)) + 1.0E-10);
}




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

