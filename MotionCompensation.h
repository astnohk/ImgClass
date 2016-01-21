#ifndef LIB_ImgClass_MotionCompensation
#define LIB_ImgClass_MotionCompensation

#include <vector>

#include "ImgClass.h"
#include "Vector.h"

template <class T>
class MotionCompensation
{
	private:
		int _width;
		int _height;
		ImgVector<T> _image_prev;
		ImgVector<T> _image_current;
		ImgVector<T> _image_next;
		ImgVector<VECTOR_2D<double> > _vector_prev;
		ImgVector<VECTOR_2D<double> > _vector_next;
		ImgVector<Vector_ST<double> > _vector_time;
		ImgVector<T> _image_compensated;

	public:
		MotionCompensation(void);
		explicit MotionCompensation(const MotionCompensation& copy); // copy constructor
		MotionCompensation(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<VECTOR_2D<double> >& vector_prev);
		MotionCompensation(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<Vector_ST<double> >& vector_prev);
		MotionCompensation(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<T>& image_next, const std::vector<ImgVector<Vector_ST<double> > >& vectors);
		MotionCompensation(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<T>& image_next, const ImgVector<Vector_ST<double> >& vector_time);
		~MotionCompensation(void);

		MotionCompensation& copy(const MotionCompensation& copy);
		MotionCompensation& set(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<VECTOR_2D<double> >& vector_prev);
		MotionCompensation& set(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<Vector_ST<double> >& vector_prev);
		MotionCompensation& set(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<T>& image_next, const std::vector<ImgVector<Vector_ST<double> > >& vectors);
		MotionCompensation& set(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<T>& image_next, const ImgVector<Vector_ST<double> >& vector_time);

		// Accessor
		// * parameter
		int width(void) const;
		int height(void) const;

		// * reference
		const ImgVector<VECTOR_2D<double> >& ref_vector_prev(void) const;
		const ImgVector<VECTOR_2D<double> >& ref_vector_next(void) const;
		const ImgVector<Vector_ST<double> >& ref_vector_time(void) const;

		const ImgVector<T>& ref_image_compensated(void) const;
		T& at_image_compensated(const int x, const int y) const;
		T& operator[](const size_t& n) const; // Get reference to compensated_image

		// * original image intensity
		T get_image_prev(const size_t& n) const;
		T get_image_prev(const int x, const int y) const;
		T get_image_current(const size_t& n) const;
		T get_image_current(const int x, const int y) const;
		T get_image_next(const size_t& n) const;
		T get_image_next(const int x, const int y) const;
		// * vector
		// returns with qualifier const to avoid to mistake get_vector() for at_vector()
		const VECTOR_2D<double> get_vector_prev(const size_t& n) const;
		const VECTOR_2D<double> get_vector_prev(const int x, const int y) const;
		const VECTOR_2D<double> get_vector_next(const size_t& n) const;
		const VECTOR_2D<double> get_vector_next(const int x, const int y) const;
		const Vector_ST<double> get_vector_time(const size_t& n) const;
		const Vector_ST<double> get_vector_time(const int x, const int y) const;
		// * compensated image
		T get_image_compensated(const size_t& n) const;
		T get_image_compensated(const int x, const int y) const;

		// Motion compensation methods
		void create_image_compensated(void);
		void create_image_estimated(const double estimate_time);
};

#include "MotionCompensation_private.h"

#endif

