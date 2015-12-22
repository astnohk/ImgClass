#include "BlockMatching.h"
#include "Lab.h"
#include "RGB.h"




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
	    / (sqrt(N * sum_sq_prev - inner_prod(sum_prev, sum_prev)) * sqrt(N * sum_sq_next - inner_prod(sum_next, sum_next)));
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
	    / (sqrt(N * sum_sq_prev - inner_prod(sum_prev, sum_prev)) * sqrt(N * sum_sq_next - inner_prod(sum_next, sum_next)));
}

