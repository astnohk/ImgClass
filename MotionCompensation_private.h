#include <cassert>
#include <cstdio>
#include <new>



template <class T>
MotionCompensation<T>::MotionCompensation(void)
{
	_width = 0;
	_height = 0;
}

template <class T>
MotionCompensation<T>::MotionCompensation(const MotionCompensation &copy) // copy constructor
{
	_width = copy._width;
	_height = copy._height;

	_image_prev.copy(copy._image_prev);
	_image_current.copy(copy._image_current);
	_image_next.copy(copy._image_next);

	_vector_prev.copy(copy._vector_prev);
	_vector_next.copy(copy._vector_next);

	_image_compensated.copy(copy._image_compensated);
}


template <class T>
MotionCompensation<T>::MotionCompensation(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<VECTOR_2D<double> >& vector_prev)
{
	_width = 0;
	_height = 0;
	if (image_prev.isNULL()) {
		throw std::invalid_argument("MotionCompensation<T>::MotionCompensation(const ImgVector<double>&, const ImgVector<double>&, const VECTOR_2D<double>&) : image_prev is empty");
	} else if (image_current.isNULL()) {
		throw std::invalid_argument("MotionCompensation<T>::MotionCompensation(const ImgVector<double>&, const ImgVector<double>&, const VECTOR_2D<double>&) : image_current is empty");
	} else if (vector_prev.isNULL()) {
		throw std::invalid_argument("MotionCompensation<T>::MotionCompensation(const ImgVector<double>&, const ImgVector<double>&, const VECTOR_2D<double>&) : vector_prev is empty");
	}
	_width = image_prev.width();
	_height = image_prev.height();
	_image_prev.copy(image_prev);
	_image_current.copy(image_current);
	if (vector_prev.width() == _width && vector_prev.height() == _height) {
		_vector_prev.copy(vector_prev);
	} else {
		// Projection of small vector field to the scaled plane which has same range of images
		_vector_prev.reset(_width, _height);
		for (int y = 0; y < _height; y++) {
			int Y = int(floor(y * vector_prev.height() / _height));
			for (int x = 0; x < _width; x++) {
				int X = int(floor(x * vector_prev.width() / _width));
				_vector_prev.at(x, y) = vector_prev.get(X, Y);
			}
		}
	}
	_image_compensated.reset(_width, _height);
	this->create_image_compensated();
}

template <class T>
MotionCompensation<T>::MotionCompensation(const ImgVector<T> &image_prev, const ImgVector<T> &image_current, const ImgVector<T> &image_next, const ImgVector<VECTOR_2D<double> > &vector_prev, const ImgVector<VECTOR_2D<double> > &vector_next)
{
	_width = 0;
	_height = 0;
	if (image_prev.isNULL()) {
		throw std::invalid_argument("MotionCompensation<T>::MotionCompensation(const ImgVector<double>&, const ImgVector<double>&, const VECTOR_2D<double>&) : image_prev is empty");
	} else if (image_current.isNULL()) {
		throw std::invalid_argument("MotionCompensation<T>::MotionCompensation(const ImgVector<double>&, const ImgVector<double>&, const VECTOR_2D<double>&) : image_current is empty");
	} else if (image_next.isNULL()) {
		throw std::invalid_argument("MotionCompensation<T>::MotionCompensation(const ImgVector<double>&, const ImgVector<double>&, const VECTOR_2D<double>&) : image_next is empty");
	} else if (vector_prev.isNULL()) {
		throw std::invalid_argument("MotionCompensation<T>::MotionCompensation(const ImgVector<double>&, const ImgVector<double>&, const VECTOR_2D<double>&) : vector_prev is empty");
	} else if (vector_next.isNULL()) {
		throw std::invalid_argument("MotionCompensation<T>::MotionCompensation(const ImgVector<double>&, const ImgVector<double>&, const VECTOR_2D<double>&) : vector_next is empty");
	}
	_width = image_prev.width();
	_height = image_prev.height();

	_image_prev.copy(image_prev);
	_image_current.copy(image_current);
	_image_next.copy(image_next);
	if (vector_prev.width() == _width && vector_prev.height() == _height) {
		_vector_prev.copy(vector_prev);
		_vector_next.copy(vector_next);
	} else {
		// Projection of small vector field to the scaled plane which has same range of images
		_vector_prev.reset(_width, _height);
		_vector_next.reset(_width, _height);
		for (int y = 0; y < _height; y++) {
			int Y = int(floor(y * vector_prev.height() / _height));
			for (int x = 0; x < _width; x++) {
				int X = int(floor(x * vector_prev.width() / _width));
				_vector_prev.at(x, y) = vector_prev.get(X, Y);
				_vector_next.at(x, y) = vector_next.get(X, Y);
			}
		}
	}
	_image_compensated.reset(_width, _height);
	this->create_image_compensated();
}


template <class T>
MotionCompensation<T>::~MotionCompensation(void) // Destructor
{
}


template <class T>
MotionCompensation<T> &
MotionCompensation<T>::copy(const MotionCompensation &copy)
{
	_width = copy._width;
	_height = copy._height;
	_image_prev.copy(copy._image_prev);
	_image_current.copy(copy._image_current);
	_image_next.copy(copy._image_next);
	_vector_prev.copy(copy._vector_prev);
	_vector_next.copy(copy._vector_next);

	_image_compensated.copy(copy._image_compensated);
	return *this;
}


template <class T>
MotionCompensation<T> &
MotionCompensation<T>::set(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<VECTOR_2D<double> >& vector_prev)
{
	if (image_prev.isNULL()) {
		throw std::invalid_argument("MotionCompensation& MotionCompensation<T>::set(const ImgVector<double>&, const ImgVector<double>&, const ImgVector<VECTOR_2D<double> >&) : image_prev is empty");
	} else if (image_current.isNULL()) {
		throw std::invalid_argument("MotionCompensation& MotionCompensation<T>::set(const ImgVector<double>&, const ImgVector<double>&, const ImgVector<VECTOR_2D<double> >&) : image_current is empty");
	} else if (vector_prev.isNULL()) {
		throw std::invalid_argument("MotionCompensation& MotionCompensation<T>::set(const ImgVector<double>&, const ImgVector<double>&, const ImgVector<VECTOR_2D<double> >&) : vector_prev is empty");
	}

	_width = image_prev.width();
	_height = image_prev.height();

	_image_prev.copy(image_prev);
	_image_current.copy(image_current);
	_image_next.clear();

	_vector_next.clear();
	if (vector_prev.width() == _width && vector_prev.height() == _height) {
		_vector_prev.copy(vector_prev);
	} else {
		// Projection of small vector field to the scaled plane which has same range of images
		_vector_prev.reset(_width, _height);
		for (int y = 0; y < _height; y++) {
			int Y = int(floor(y * vector_prev.height() / _height));
			for (int x = 0; x < _width; x++) {
				int X = int(floor(x * vector_prev.width() / _width));
				_vector_prev.at(x, y) = vector_prev.get(X, Y);
			}
		}
	}
	_image_compensated.reset(_width, _height);
	this->create_image_compensated();
	return *this;
}

template <class T>
MotionCompensation<T> &
MotionCompensation<T>::set(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<T>& image_next, const ImgVector<VECTOR_2D<double> >& vector_prev, const ImgVector<VECTOR_2D<double> >& vector_next)
{
	if (image_prev.isNULL()) {
		throw std::invalid_argument("MotionCompensation& MotionCompensation<T>::set(const ImgVector<double>&, const ImgVector<double>&, const ImgVector<VECTOR_2D<double> >&) : image_prev is empty");
	} else if (image_current.isNULL()) {
		throw std::invalid_argument("MotionCompensation& MotionCompensation<T>::set(const ImgVector<double>&, const ImgVector<double>&, const ImgVector<VECTOR_2D<double> >&) : image_current is empty");
	} else if (image_next.isNULL()) {
		throw std::invalid_argument("MotionCompensation& MotionCompensation<T>::set(const ImgVector<double>&, const ImgVector<double>&, const ImgVector<VECTOR_2D<double> >&) : image_next is empty");
	} else if (vector_prev.isNULL()) {
		throw std::invalid_argument("MotionCompensation& MotionCompensation<T>::set(const ImgVector<double>&, const ImgVector<double>&, const ImgVector<VECTOR_2D<double> >&) : vector_prev is empty");
	} else if (vector_next.isNULL()) {
		throw std::invalid_argument("MotionCompensation& MotionCompensation<T>::set(const ImgVector<double>&, const ImgVector<double>&, const ImgVector<VECTOR_2D<double> >&) : vector_next is empty");
	}

	_width = image_prev.width();
	_height = image_prev.height();

	_image_prev.copy(image_prev);
	_image_current.copy(image_current);
	_image_next.copy(image_next);
	if (vector_prev.width() == _width && vector_prev.height() == _height) {
		_vector_prev.copy(vector_prev);
		_vector_next.copy(vector_next);
	} else {
		// Projection of small vector field to the scaled plane which has same range of images
		_vector_prev.reset(_width, _height);
		_vector_next.reset(_width, _height);
		for (int y = 0; y < _height; y++) {
			int Y = int(floor(y * vector_prev.height() / _height));
			for (int x = 0; x < _width; x++) {
				int X = int(floor(x * vector_prev.width() / _width));
				_vector_prev.at(x, y) = vector_prev.get(X, Y);
				_vector_next.at(x, y) = vector_next.get(X, Y);
			}
		}
	}
	_image_compensated.reset(_width, _height);
	this->create_image_compensated();
	return *this;
}




template <class T>
int
MotionCompensation<T>::width(void) const
{
	return _width;
}

template <class T>
int
MotionCompensation<T>::height(void) const
{
	return _height;
}

// Reference
template <class T>
const ImgVector<VECTOR_2D<double> > &
MotionCompensation<T>::ref_vector_prev(void) const
{
	return _vector_prev;
}

template <class T>
const ImgVector<VECTOR_2D<double> > &
MotionCompensation<T>::ref_vector_next(void) const
{
	return _vector_next;
}

template <class T>
const ImgVector<T> &
MotionCompensation<T>::ref_image_compensated(void) const
{
	return _image_compensated;
}


template <class T>
T &
MotionCompensation<T>::at_image_compensated(int x, int y) const
{
	return _image_compensated[_width * y + x];
}

template <class T>
T &
MotionCompensation<T>::operator[](int n) const // Get reference to the compensated image[n]
{
	return _image_compensated[n];
}


// Get
template <class T>
T
MotionCompensation<T>::get_image_prev(int n) const
{
	return _image_prev.get(n);
}

template <class T>
T
MotionCompensation<T>::get_image_prev(int x, int y) const
{
	return _image_prev.get(x, y);
}

template <class T>
T
MotionCompensation<T>::get_image_current(int n) const
{
	return _image_current.get(n);
}

template <class T>
T
MotionCompensation<T>::get_image_current(int x, int y) const
{
	return _image_current.get(x, y);
}

template <class T>
T
MotionCompensation<T>::get_image_next(int n) const
{
	return _image_next.get(n);
}

template <class T>
T
MotionCompensation<T>::get_image_next(int x, int y) const
{
	return _image_next.get(x, y);
}


template <class T>
const VECTOR_2D<double>
MotionCompensation<T>::get_vector_prev(int n) const
{
	return _vector_prev.get(n);
}

template <class T>
const VECTOR_2D<double>
MotionCompensation<T>::get_vector_prev(int x, int y) const
{
	return _vector_prev.get(x, y);
}

template <class T>
const VECTOR_2D<double>
MotionCompensation<T>::get_vector_next(int n) const
{
	return _vector_next.get(n);
}

template <class T>
const VECTOR_2D<double>
MotionCompensation<T>::get_vector_next(int x, int y) const
{
	return _vector_next.get(x, y);
}

template <class T>
T
MotionCompensation<T>::get_image_compensated(int n) const
{
	assert(0 <= n && n < _width * _height);
	return _image_compensated.get(n);
}

template <class T>
T
MotionCompensation<T>::get_image_compensated(int x, int y) const
{
	assert(0 <= x && x < _width && 0 <= y && y < _height);
	return _image_compensated.get(x, y);
}




/*
  void MotionCompensation<T>::create_image_masked_compensated(ImgVector<bool> *mask)

  Make compensated image.
  If the mask is specified then it will only compensate the pixel of masked.
  ImgVector<bool> *mask should mean mask of compensated image.
  If mask(x, y) == true then the pixel would be compensated
  and if mask(x, y) == false then the pixel hold the original (image_current) intensity.
*/
template <class T>
void
MotionCompensation<T>::create_image_compensated(const ImgVector<bool>* mask)
{
	if (mask == nullptr) {
		_image_compensated.reset(_width, _height);
		for (int y = 0; y < _height; y++) {
			for (int x = 0; x < _width; x++) {
				VECTOR_2D<int> r(int(round(x + _vector_prev.get(x, y).x)), int(round(y + _vector_prev.get(x, y).y)));
				_image_compensated.at(x, y) = _image_prev.get_zeropad(r.x, r.y);
			}
		}
	} else {
		_image_compensated.copy(_image_current); // Initialize with Original image (image_current)
		for (int y = 0; y < _height; y++) {
			for (int x = 0; x < _width; x++) {
				if (mask->get(x, y) == false) {
					continue;
				}
				VECTOR_2D<int> v(int(round(_vector_prev.get(x, y).x)), int(round(_vector_prev.get(x, y).y)));
				_image_compensated.at(x, y) = _image_prev.get_zeropad(x + v.x, y + v.y);
			}
		}
	}
}


/*
  void MotionCompensation<T>::create_image_masked_compensated(ImgVector<bool> *mask)

  Make estimated (compensated) image.
  if estimate_frame = 2 then it estimate the image of next of next frame.
*/
template <class T>
void
MotionCompensation<T>::create_image_estimated(const double estimate_time, const ImgVector<bool>* mask)
{
	if (mask == nullptr) {
		_image_compensated.reset(_width, _height);
		for (int y = 0; y < _height; y++) {
			for (int x = 0; x < _width; x++) {
				VECTOR_2D<double> v = estimate_time * _vector_prev.get(x, y);
				_image_compensated.at(x, y) = _image_prev.get_zeropad(x + v.x, y + v.y);
			}
		}
	} else {
		_image_compensated.copy(_image_current); // Initialize with Original image (image_current)
		for (int y = 0; y < _height; y++) {
			for (int x = 0; x < _width; x++) {
				if (mask->get(x, y) == false) {
					continue;
				}
				VECTOR_2D<double> v = estimate_time * _vector_prev.get(x, y);
				_image_compensated.at(x, y) = _image_prev.get_zeropad(x + v.x, y + v.y);
			}
		}
	}
}

