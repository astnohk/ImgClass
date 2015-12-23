#include "BlockMatching.h"




// ----- Miscellaneous -----
template <>
void
BlockMatching<ImgClass::RGB>::image_normalizer(void)
{
	double max_int = 0.0;
	for (int i = 0; i < _width * _height; i++) {
		if (norm(_image_prev[i]) > max_int) {
			max_int = norm(_image_prev[i]);
		}
	}
	if (max_int > 1.0) {
		_image_prev /= max_int;
		_image_next /= max_int;
	}
}

template <>
void
BlockMatching<ImgClass::Lab>::image_normalizer(void)
{
	double max_int = 0.0;
	for (int i = 0; i < _width * _height; i++) {
		if (norm(_image_prev[i]) > max_int) {
			max_int = norm(_image_prev[i]);
		}
	}
	if (max_int > 1.0) {
		_image_prev /= max_int;
		_image_next /= max_int;
	}
}




template <>
ImgVector<VECTOR_2D<double> > *
BlockMatching<ImgClass::Lab>::grad_prev(const int top_left_x, const int top_left_y, const int crop_width, const int crop_height)
{
	ImgVector<VECTOR_2D<double> >* gradients = new ImgVector<VECTOR_2D<double> >(crop_width, crop_height);

	for (int y = 0; y < crop_height; y++) {
		for (int x = 0; x < crop_width; x++) {
			gradients->at(x, y).x =
			    _image_prev.get_mirror(top_left_x + x + 1, top_left_y + y).L
			    - _image_prev.get_mirror(top_left_x + x, top_left_y + y).L;
			gradients->at(x, y).y =
			    _image_prev.get_mirror(top_left_x + x, top_left_y + y + 1).L
			    - _image_prev.get_mirror(top_left_x + x, top_left_y + y).L;
		}
	}
	return gradients;
}




// ----- Correlation functions -----
template <>
double
BlockMatching<ImgClass::RGB>::MAD(const int x_l, const int y_l, const int x_r, const int y_r, const int block_width, const int block_height, const ImgVector<ImgClass::RGB>& limage, const ImgVector<ImgClass::RGB>& rimage)
{
	double sad = 0;

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
	double sad = 0;

	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			sad = sad + norm(limage.get_zeropad(x_l+ x, y_l+ y) - rimage.get_zeropad(x_r + x, y_r + y));
		}
	}
	return sad / double(block_width * block_height);
}




template <>
double
BlockMatching<ImgClass::RGB>::MAD_centered(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height)
{
	double sad = 0;
	double N = 0.0;

	ImgClass::RGB center_color = _image_next.get_zeropad(x_next + block_width / 2, y_next + block_height / 2);
	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			double coeff = 1.0 - norm(center_color - _image_prev.get_zeropad(x_prev + x, y_prev + y));
			N += coeff;
			sad += coeff * norm(_image_next.get_zeropad(x_next + x, y_next + y) - _image_prev.get_zeropad(x_prev + x, y_prev + y));
		}
	}
	return sad / N;
}

template <>
double
BlockMatching<ImgClass::Lab>::MAD_centered(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height)
{
	double sad = 0;
	double N = 0.0;

	ImgClass::Lab center_color = _image_next.get_zeropad(x_next + block_width / 2, y_next + block_height / 2);
	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			double coeff = 1.0 - norm(center_color - _image_prev.get_zeropad(x_prev + x, y_prev + y));
			N += coeff;
			sad += coeff * norm(_image_next.get_zeropad(x_next + x, y_next + y) - _image_prev.get_zeropad(x_prev + x, y_prev + y));
		}
	}
	return sad / N;
}




template <>
double
BlockMatching<ImgClass::RGB>::MAD(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height)
{
	double sad = 0;

	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			ImgClass::RGB color_prev(_image_prev.get_zeropad(x_prev + x, y_prev + y));
			ImgClass::RGB color_next(_image_next.get_zeropad(x_next + x, y_next + y));
			sad = sad + norm(color_next - color_prev);
		}
	}
	return sad / double(block_width * block_height);
}

template <>
double
BlockMatching<ImgClass::Lab>::MAD(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height)
{
	double sad = 0;

	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			ImgClass::Lab color_prev(_image_prev.get_zeropad(x_prev + x, y_prev + y));
			ImgClass::Lab color_next(_image_next.get_zeropad(x_next + x, y_next + y));
			sad = sad + norm(color_next - color_prev);
		}
	}
	return sad / double(block_width * block_height);
}




template <>
double
BlockMatching<ImgClass::RGB>::ZNCC(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height)
{
	double N = block_width * block_height;
	ImgClass::RGB sum_prev(.0, .0, .0);
	ImgClass::RGB sum_next(.0, .0, .0);
	double sum_sq_prev = 0;
	double sum_sq_next = 0;
	double sum_sq_prev_next = 0;

	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			sum_prev += _image_prev.get_zeropad(x_prev + x, y_prev + y);
			sum_next += _image_next.get_zeropad(x_next + x, y_next + y);
			sum_sq_prev += inner_prod(_image_prev.get_zeropad(x_prev + x, y_prev + y), _image_prev.get_zeropad(x_prev + x, y_prev + y));
			sum_sq_next += inner_prod(_image_next.get_zeropad(x_next + x, y_next + y), _image_next.get_zeropad(x_next + x, y_next + y));
			sum_sq_prev_next += inner_prod(_image_prev.get_zeropad(x_prev + x, y_prev + y), _image_next.get_zeropad(x_next + x, y_next + y));
		}
	}
	return (N * sum_sq_prev_next - inner_prod(sum_prev, sum_next)) / (sqrt(N * sum_sq_prev - inner_prod(sum_prev, sum_prev)) * sqrt(N * sum_sq_next - inner_prod(sum_next, sum_next)) + 1E-10);
}

template <>
double
BlockMatching<ImgClass::Lab>::ZNCC(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height)
{
	double N = block_width * block_height;
	ImgClass::Lab sum_prev(.0, .0, .0);
	ImgClass::Lab sum_next(.0, .0, .0);
	double sum_sq_prev = 0;
	double sum_sq_next = 0;
	double sum_sq_prev_next = 0;

	for (int y = 0; y < block_height; y++) {
		for (int x = 0; x < block_width; x++) {
			sum_prev += _image_prev.get_zeropad(x_prev + x, y_prev + y);
			sum_next += _image_next.get_zeropad(x_next + x, y_next + y);
			sum_sq_prev += inner_prod(_image_prev.get_zeropad(x_prev + x, y_prev + y), _image_prev.get_zeropad(x_prev + x, y_prev + y));
			sum_sq_next += inner_prod(_image_next.get_zeropad(x_next + x, y_next + y), _image_next.get_zeropad(x_next + x, y_next + y));
			sum_sq_prev_next += inner_prod(_image_prev.get_zeropad(x_prev + x, y_prev + y), _image_next.get_zeropad(x_next + x, y_next + y));
		}
	}
	return (N * sum_sq_prev_next - inner_prod(sum_prev, sum_next)) / (sqrt(N * sum_sq_prev - inner_prod(sum_prev, sum_prev)) * sqrt(N * sum_sq_next - inner_prod(sum_next, sum_next)) + 1E-10);
}




template <>
double
BlockMatching<ImgClass::RGB>::MAD_region(const int x_diff_prev, const int y_diff_prev, const std::list<VECTOR_2D<int> >& region)
{
	double N = .0;
	double sad = 0;

	for (std::list<VECTOR_2D<int> >::const_iterator ite = region.begin();
	    ite != region.end();
	    ++ite) {
		N += 1.0;
		ImgClass::RGB color_prev(_image_prev.get_zeropad(ite->x + x_diff_prev, ite->y + y_diff_prev));
		ImgClass::RGB color_next(_image_next.get_zeropad(ite->x, ite->y));
		sad = sad + norm(color_next - color_prev);
	}
	return sad / N;
}

template <>
double
BlockMatching<ImgClass::Lab>::MAD_region(const int x_diff_prev, const int y_diff_prev, const std::list<VECTOR_2D<int> >& region)
{
	double N = .0;
	double sad = 0;

	for (std::list<VECTOR_2D<int> >::const_iterator ite = region.begin();
	    ite != region.end();
	    ++ite) {
		N += 1.0;
		ImgClass::Lab color_prev(_image_prev.get_zeropad(ite->x + x_diff_prev, ite->y + y_diff_prev));
		ImgClass::Lab color_next(_image_next.get_zeropad(ite->x, ite->y));
		sad = sad + norm(color_next - color_prev);
	}
	return sad / N;
}


template <>
double
BlockMatching<ImgClass::RGB>::ZNCC_region(const int x_diff_prev, const int y_diff_prev, const std::list<VECTOR_2D<int> >& region)
{
	double N = 0;
	ImgClass::RGB sum_prev(.0, .0, .0);
	ImgClass::RGB sum_next(.0, .0, .0);
	double sum_sq_prev = 0;
	double sum_sq_next = 0;
	double sum_sq_prev_next = 0;

	for (std::list<VECTOR_2D<int> >::const_iterator ite = region.begin();
	    ite != region.end();
	    ++ite) {
		N += 1.0;
		// Previous frame
		sum_prev += _image_prev.get_zeropad(ite->x + x_diff_prev, ite->y + y_diff_prev);
		sum_sq_prev += inner_prod(
		    _image_prev.get_zeropad(ite->x + x_diff_prev, ite->y + y_diff_prev)
		    , _image_prev.get_zeropad(ite->x + x_diff_prev, ite->y + y_diff_prev));
		// Next frame
		sum_next += _image_next.get_zeropad(ite->x, ite->y);
		sum_sq_next += inner_prod(
		    _image_next.get_zeropad(ite->x, ite->y)
		    , _image_next.get_zeropad(ite->x, ite->y));
		// Co-frame
		sum_sq_prev_next += inner_prod(
		    _image_prev.get_zeropad(ite->x + x_diff_prev, ite->y + y_diff_prev)
		    , _image_next.get_zeropad(ite->x, ite->y));
	}
	// Calculate Covariance
	return (N * sum_sq_prev_next - inner_prod(sum_prev, sum_next))
	    / (sqrt(N * sum_sq_prev - inner_prod(sum_prev, sum_prev)) * sqrt(N * sum_sq_next - inner_prod(sum_next, sum_next)) + 1.0E-10);
}

template <>
double
BlockMatching<ImgClass::Lab>::ZNCC_region(const int x_diff_prev, const int y_diff_prev, const std::list<VECTOR_2D<int> >& region)
{
	double N = 0;
	ImgClass::Lab sum_prev(.0, .0, .0);
	ImgClass::Lab sum_next(.0, .0, .0);
	double sum_sq_prev = 0;
	double sum_sq_next = 0;
	double sum_sq_prev_next = 0;

	for (std::list<VECTOR_2D<int> >::const_iterator ite = region.begin();
	    ite != region.end();
	    ++ite) {
		N += 1.0;
		// Previous frame
		sum_prev += _image_prev.get_zeropad(ite->x + x_diff_prev, ite->y + y_diff_prev);
		sum_sq_prev += inner_prod(
		    _image_prev.get_zeropad(ite->x + x_diff_prev, ite->y + y_diff_prev)
		    , _image_prev.get_zeropad(ite->x + x_diff_prev, ite->y + y_diff_prev));
		// Next frame
		sum_next += _image_next.get_zeropad(ite->x, ite->y);
		sum_sq_next += inner_prod(
		    _image_next.get_zeropad(ite->x, ite->y)
		    , _image_next.get_zeropad(ite->x, ite->y));
		// Co-frame
		sum_sq_prev_next += inner_prod(
		    _image_prev.get_zeropad(ite->x + x_diff_prev, ite->y + y_diff_prev)
		    , _image_next.get_zeropad(ite->x, ite->y));
	}
	// Calculate Covariance
	return (N * sum_sq_prev_next - inner_prod(sum_prev, sum_next))
	    / (sqrt(N * sum_sq_prev - inner_prod(sum_prev, sum_prev)) * sqrt(N * sum_sq_next - inner_prod(sum_next, sum_next)) + 1.0E-10);
}




// ----- Global function -----
double
norm(const double& value)
{
	return fabs(value);
}

