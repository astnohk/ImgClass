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
	_vector_time.copy(copy._vector_time);

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
MotionCompensation<T>::MotionCompensation(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<Vector_ST<double> >& vector_prev)
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
	_vector_prev.reset(_width, _height);
	if (vector_prev.width() == _width && vector_prev.height() == _height) {
		for (size_t n = 0; n < vector_prev.size(); n++) {
			_vector_prev[n] = vector_prev[n];
		}
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
MotionCompensation<T>::MotionCompensation(const ImgVector<T> &image_prev, const ImgVector<T> &image_current, const ImgVector<T> &image_next, const std::vector<ImgVector<Vector_ST<double> > >& vectors)
{
	_width = 0;
	_height = 0;
	if (image_prev.isNULL()) {
		throw std::invalid_argument("MotionCompensation<T>::MotionCompensation(const ImgVector<double>&, const ImgVector<double>&, const VECTOR_2D<double>&) : image_prev is empty");
	} else if (image_current.isNULL()) {
		throw std::invalid_argument("MotionCompensation<T>::MotionCompensation(const ImgVector<double>&, const ImgVector<double>&, const VECTOR_2D<double>&) : image_current is empty");
	} else if (image_next.isNULL()) {
		throw std::invalid_argument("MotionCompensation<T>::MotionCompensation(const ImgVector<double>&, const ImgVector<double>&, const VECTOR_2D<double>&) : image_next is empty");
	} else if (vectors.size() < 2) {
		throw std::invalid_argument("MotionCompensation& MotionCompensation<T>::set(const ImgVector<double>&, const ImgVector<double>&, const ImgVector<VECTOR_2D<double> >&) : vectors is too small (vectors.size() < 2)");
	} else if (vectors[0].isNULL()) {
		throw std::invalid_argument("MotionCompensation& MotionCompensation<T>::set(const ImgVector<double>&, const ImgVector<double>&, const ImgVector<VECTOR_2D<double> >&) : vectors[0] is empty");
	} else if (vectors[1].isNULL()) {
		throw std::invalid_argument("MotionCompensation& MotionCompensation<T>::set(const ImgVector<double>&, const ImgVector<double>&, const ImgVector<VECTOR_2D<double> >&) : vectors[1] is empty");
	}
	_width = image_prev.width();
	_height = image_prev.height();

	_image_prev.copy(image_prev);
	_image_current.copy(image_current);
	_image_next.copy(image_next);

	_vector_prev.reset(_width, _height);
	_vector_next.reset(_width, _height);
	if (vectors[0].width() == _width && vectors[0].height() == _height) {
		for (size_t n = 0; n < vectors[0].size(); n++) {
			if (vectors[0][0].t < 0) {
				_vector_prev[n] = vectors[0][n];
				_vector_next[n] = vectors[1][n];
			} else {
				_vector_prev[n] = vectors[1][n];
				_vector_next[n] = vectors[0][n];
			}
		}
	} else {
		// Projection of small vector field to the scaled plane which has same range of images
		for (int y = 0; y < _height; y++) {
			int Y = int(floor(y * vectors[0].height() / _height));
			for (int x = 0; x < _width; x++) {
				int X = int(floor(x * vectors[0].width() / _width));
				if (vectors[0][0].t < 0) {
					_vector_prev.at(x, y) = vectors[0].get(X, Y);
					_vector_next.at(x, y) = vectors[1].get(X, Y);
				} else {
					_vector_prev.at(x, y) = vectors[1].get(X, Y);
					_vector_next.at(x, y) = vectors[0].get(X, Y);
				}
			}
		}
	}
	_image_compensated.reset(_width, _height);
	this->create_image_compensated();
}

template <class T>
MotionCompensation<T>::MotionCompensation(const ImgVector<T> &image_prev, const ImgVector<T> &image_current, const ImgVector<T> &image_next, const ImgVector<Vector_ST<double> > &vector_time)
{
	_width = 0;
	_height = 0;
	if (image_prev.isNULL()) {
		throw std::invalid_argument("MotionCompensation<T>::MotionCompensation(const ImgVector<double>&, const ImgVector<double>&, const VECTOR_2D<double>&) : image_prev is empty");
	} else if (image_current.isNULL()) {
		throw std::invalid_argument("MotionCompensation<T>::MotionCompensation(const ImgVector<double>&, const ImgVector<double>&, const VECTOR_2D<double>&) : image_current is empty");
	} else if (image_next.isNULL()) {
		throw std::invalid_argument("MotionCompensation<T>::MotionCompensation(const ImgVector<double>&, const ImgVector<double>&, const VECTOR_2D<double>&) : image_next is empty");
	} else if (vector_time.isNULL()) {
		throw std::invalid_argument("MotionCompensation<T>::MotionCompensation(const ImgVector<double>&, const ImgVector<double>&, const VECTOR_2D<double>&) : vector_time is empty");
	}
	_width = image_prev.width();
	_height = image_prev.height();

	_image_prev.copy(image_prev);
	_image_current.copy(image_current);
	_image_next.copy(image_next);
	if (vector_time.width() == _width && vector_time.height() == _height) {
		_vector_time.copy(vector_time);
	} else {
		// Projection of small vector field to the scaled plane which has same range of images
		_vector_time.reset(_width, _height);
		for (int y = 0; y < _height; y++) {
			int Y = int(floor(y * vector_time.height() / _height));
			for (int x = 0; x < _width; x++) {
				int X = int(floor(x * vector_time.width() / _width));
				_vector_time.at(x, y) = vector_time.get(X, Y);
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
	_vector_time.clear();
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
MotionCompensation<T>::set(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<Vector_ST<double> >& vector_prev)
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
	_vector_time.clear();
	if (vector_prev.width() == _width && vector_prev.height() == _height) {
		for (size_t n = 0; n < vector_prev.size(); n++) {
			_vector_prev[n] = vector_prev[n];
		}
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
MotionCompensation<T>::set(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<T>& image_next, const std::vector<ImgVector<Vector_ST<double> > >& vectors)
{
	if (image_prev.isNULL()) {
		throw std::invalid_argument("MotionCompensation& MotionCompensation<T>::set(const ImgVector<double>&, const ImgVector<double>&, const ImgVector<VECTOR_2D<double> >&) : image_prev is empty");
	} else if (image_current.isNULL()) {
		throw std::invalid_argument("MotionCompensation& MotionCompensation<T>::set(const ImgVector<double>&, const ImgVector<double>&, const ImgVector<VECTOR_2D<double> >&) : image_current is empty");
	} else if (image_next.isNULL()) {
		throw std::invalid_argument("MotionCompensation& MotionCompensation<T>::set(const ImgVector<double>&, const ImgVector<double>&, const ImgVector<VECTOR_2D<double> >&) : image_next is empty");
	} else if (vectors.size() < 2) {
		throw std::invalid_argument("MotionCompensation& MotionCompensation<T>::set(const ImgVector<double>&, const ImgVector<double>&, const ImgVector<VECTOR_2D<double> >&) : vectors is too small (vectors.size() < 2)");
	} else if (vectors[0].isNULL()) {
		throw std::invalid_argument("MotionCompensation& MotionCompensation<T>::set(const ImgVector<double>&, const ImgVector<double>&, const ImgVector<VECTOR_2D<double> >&) : vectors[0] is empty");
	} else if (vectors[1].isNULL()) {
		throw std::invalid_argument("MotionCompensation& MotionCompensation<T>::set(const ImgVector<double>&, const ImgVector<double>&, const ImgVector<VECTOR_2D<double> >&) : vectors[1] is empty");
	}

	_width = image_prev.width();
	_height = image_prev.height();

	_image_prev.copy(image_prev);
	_image_current.copy(image_current);
	_image_next.copy(image_next);

	_vector_prev.reset(_width, _height);
	_vector_next.reset(_width, _height);
	_vector_time.clear();
	if (vectors[0].width() == _width && vectors[0].height() == _height) {
		for (size_t n = 0; n < vectors[0].size(); n++) {
			if (vectors[0][0].t < 0) {
				_vector_prev[n] = vectors[0][n];
				_vector_next[n] = vectors[1][n];
			} else {
				_vector_prev[n] = vectors[1][n];
				_vector_next[n] = vectors[0][n];
			}
		}
	} else {
		// Projection of small vector field to the scaled plane which has same range of images
		for (int y = 0; y < _height; y++) {
			int Y = int(floor(y * vectors[0].height() / _height));
			for (int x = 0; x < _width; x++) {
				int X = int(floor(x * vectors[0].width() / _width));
				if (vectors[0][0].t < 0) {
					_vector_prev.at(x, y) = vectors[0].get(X, Y);
					_vector_next.at(x, y) = vectors[1].get(X, Y);
				} else {
					_vector_prev.at(x, y) = vectors[1].get(X, Y);
					_vector_next.at(x, y) = vectors[0].get(X, Y);
				}
			}
		}
	}
	_image_compensated.reset(_width, _height);
	this->create_image_compensated();
	return *this;
}

template <class T>
MotionCompensation<T> &
MotionCompensation<T>::set(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<T>& image_next, const ImgVector<Vector_ST<double> >& vector_time)
{
	if (image_prev.isNULL()) {
		throw std::invalid_argument("MotionCompensation& MotionCompensation<T>::set(const ImgVector<double>&, const ImgVector<double>&, const ImgVector<VECTOR_2D<double> >&) : image_prev is empty");
	} else if (image_current.isNULL()) {
		throw std::invalid_argument("MotionCompensation& MotionCompensation<T>::set(const ImgVector<double>&, const ImgVector<double>&, const ImgVector<VECTOR_2D<double> >&) : image_current is empty");
	} else if (image_next.isNULL()) {
		throw std::invalid_argument("MotionCompensation& MotionCompensation<T>::set(const ImgVector<double>&, const ImgVector<double>&, const ImgVector<VECTOR_2D<double> >&) : image_next is empty");
	} else if (vector_time.isNULL()) {
		throw std::invalid_argument("MotionCompensation& MotionCompensation<T>::set(const ImgVector<double>&, const ImgVector<double>&, const ImgVector<VECTOR_2D<double> >&) : vector_time is empty");
	}

	_width = image_prev.width();
	_height = image_prev.height();

	_image_prev.copy(image_prev);
	_image_current.copy(image_current);
	_image_next.copy(image_next);
	_vector_prev.clear();
	_vector_next.clear();
	if (vector_time.width() == _width && vector_time.height() == _height) {
		_vector_time.copy(vector_time);
	} else {
		// Projection of small vector field to the scaled plane which has same range of images
		_vector_time.reset(_width, _height);
		for (int y = 0; y < _height; y++) {
			int Y = int(floor(y * vector_time.height() / _height));
			for (int x = 0; x < _width; x++) {
				int X = int(floor(x * vector_time.width() / _width));
				_vector_time.at(x, y) = vector_time.get(X, Y);
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
const ImgVector<Vector_ST<double> > &
MotionCompensation<T>::ref_vector_time(void) const
{
	return _vector_time;
}

template <class T>
const ImgVector<T> &
MotionCompensation<T>::ref_image_compensated(void) const
{
	return _image_compensated;
}


template <class T>
T &
MotionCompensation<T>::at_image_compensated(const int x, const int y) const
{
	return _image_compensated[size_t(_width) * size_t(y) + size_t(x)];
}

template <class T>
T &
MotionCompensation<T>::operator[](const size_t& n) const // Get reference to the compensated image[n]
{
	return _image_compensated[n];
}


// Get
template <class T>
T
MotionCompensation<T>::get_image_prev(const size_t& n) const
{
	return _image_prev.get(n);
}

template <class T>
T
MotionCompensation<T>::get_image_prev(const int x, const int y) const
{
	return _image_prev.get(x, y);
}

template <class T>
T
MotionCompensation<T>::get_image_current(const size_t& n) const
{
	return _image_current.get(n);
}

template <class T>
T
MotionCompensation<T>::get_image_current(const int x, const int y) const
{
	return _image_current.get(x, y);
}

template <class T>
T
MotionCompensation<T>::get_image_next(const size_t& n) const
{
	return _image_next.get(n);
}

template <class T>
T
MotionCompensation<T>::get_image_next(const int x, const int y) const
{
	return _image_next.get(x, y);
}


template <class T>
const VECTOR_2D<double>
MotionCompensation<T>::get_vector_prev(const size_t& n) const
{
	return _vector_prev.get(n);
}

template <class T>
const VECTOR_2D<double>
MotionCompensation<T>::get_vector_prev(const int x, const int y) const
{
	return _vector_prev.get(x, y);
}

template <class T>
const VECTOR_2D<double>
MotionCompensation<T>::get_vector_next(const size_t& n) const
{
	return _vector_next.get(n);
}

template <class T>
const VECTOR_2D<double>
MotionCompensation<T>::get_vector_next(const int x, const int y) const
{
	return _vector_next.get(x, y);
}

template <class T>
const Vector_ST<double>
MotionCompensation<T>::get_vector_time(const size_t& n) const
{
	return _vector_time.get(n);
}

template <class T>
const Vector_ST<double>
MotionCompensation<T>::get_vector_time(const int x, const int y) const
{
	return _vector_time.get(x, y);
}

template <class T>
T
MotionCompensation<T>::get_image_compensated(const size_t& n) const
{
	assert(0 <= n && n < size_t(_width) * size_t(_height));
	return _image_compensated.get(n);
}

template <class T>
T
MotionCompensation<T>::get_image_compensated(const int x, const int y) const
{
	assert(0 <= x && x < _width && 0 <= y && y < _height);
	return _image_compensated.get(x, y);
}




/*
 * void MotionCompensation<T>::create_image_masked_compensated(void)
 * Make compensated image.
*/
template <class T>
void
MotionCompensation<T>::create_image_compensated(void)
{
	VECTOR_2D<int> r;
	VECTOR_2D<int> r_next;
	_image_compensated.reset(_width, _height);
	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			if (_vector_time.isNULL()) {
				if (_vector_next.isNULL()) { // Forward motion compensation
					r.x = int(round(x + _vector_prev.get(x, y).x));
					r.y = int(round(y + _vector_prev.get(x, y).y));
					_image_compensated.at(x, y) = _image_prev.get_zeropad(r.x, r.y);
				} else { // Mean bi-directional motion compensation
					r.x = int(round(x + _vector_prev.get(x, y).x));
					r.y = int(round(y + _vector_prev.get(x, y).y));
					r_next.x = int(round(x + _vector_next.get(x, y).x));
					r_next.y = int(round(y + _vector_next.get(x, y).y));
					_image_compensated.at(x, y) =
					    (_image_prev.get_zeropad(r.x, r.y) + _image_prev.get_zeropad(r.x, r.y)) * 0.5;
				}
			} else {
				r.x = int(round(x + _vector_time.get(x, y).x));
				r.y = int(round(y + _vector_time.get(x, y).y));
				if (_vector_time.get(x, y).t < 0) {
					_image_compensated.at(x, y) = _image_prev.get_zeropad(r.x, r.y);
				} else {
#ifdef IMG_CLASS_BLOCKMATCHING_OCCLUSION_ZEROPAD
					_image_compensated.at(x, y) = 0;
#else
					_image_compensated.at(x, y) = _image_next.get_zeropad(r.x, r.y);
#endif
				}
			}
		}
	}
}


/*
 * void MotionCompensation<T>::create_image_masked_compensated(const double estimate_time)
 * Make estimated (compensated) image.
 * if estimate_frame = 2 then it estimate the image of next of next frame.
*/
template <class T>
void
MotionCompensation<T>::create_image_estimated(const double estimate_time)
{
	if (_vector_time.isNULL() == false) {
		std::cerr << "Estimate by motion compensation is only for forward motion compensation" << std::endl;
		throw std::invalid_argument("_vector_prev is NULL");
	}
	VECTOR_2D<int> r;
	VECTOR_2D<int> r_est;
	_image_compensated.reset(_width, _height);
	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			r.x = int(round(x + _vector_prev.get(x, y).x));
			r.y = int(round(y + _vector_prev.get(x, y).y));
			r_est = r * (1.0 - estimate_time);
			if (0 <= r_est.x && r_est.x < _width
			    && 0 <= r_est.y && r_est.y < _height) {
				_image_compensated.at(r_est.x, r_est.y) = _image_prev.get_zeropad(r.x, r.y);
			}
		}
	}
}

