#include "BlockMatching.h"




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

