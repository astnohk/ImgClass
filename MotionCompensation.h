#ifndef LIB_ImgClass_MotionCompensation
#define LIB_ImgClass_MotionCompensation

#include "ImgClass.h"
#include "Vector.h"

template <class T>
class MotionCompensation
{
	private:
		bool motion_compensated;
		int _width;
		int _height;
		ImgVector<T> _image_prev;
		ImgVector<T> _image_current;
		ImgVector<T> _image_next;
		ImgVector<VECTOR_2D<double> > _vector_prev;
		ImgVector<VECTOR_2D<double> > _vector_next;

		ImgVector<T> _image_compensated;

	public:
		MotionCompensation(void);
		explicit MotionCompensation(const MotionCompensation& copy); // copy constructor
		MotionCompensation(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<VECTOR_2D<double> >& vector_prev);
		MotionCompensation(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<T>& image_next, const ImgVector<VECTOR_2D<double> >& vector_prev, const ImgVector<VECTOR_2D<double> >& vector_next);
		~MotionCompensation(void);

		MotionCompensation& copy(const MotionCompensation& copy);
		MotionCompensation& set(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<VECTOR_2D<double> >& vector_prev);
		MotionCompensation& set(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<T>& image_next, const ImgVector<VECTOR_2D<double> >& vector_prev, const ImgVector<VECTOR_2D<double> >& vector_next);

		// Accessor
		// * parameter
		int width(void) const;
		int height(void) const;

		// * reference
		const ImgVector<VECTOR_2D<double> >& ref_vector_prev(void);
		const ImgVector<VECTOR_2D<double> >& ref_vector_next(void);

		const ImgVector<T>& ref_image_compensated(void) const;
		T& at_image_compensated(int x, int y) const;
		T& operator[](int n) const; // Get reference to compensated_image

		// * original image intensity
		T get_image_prev(int n) const;
		T get_image_prev(int x, int y) const;
		T get_image_current(int n) const;
		T get_image_current(int x, int y) const;
		T get_image_next(int n) const;
		T get_image_next(int x, int y) const;
		// * vector
		// returns with qualifier const to avoid to mistake get_vector() for at_vector()
		const VECTOR_2D<double> get_vector_prev(int n) const;
		const VECTOR_2D<double> get_vector_prev(int x, int y) const;
		const VECTOR_2D<double> get_vector_next(int n) const;
		const VECTOR_2D<double> get_vector_next(int x, int y) const;
		// * compensated image
		T get_image_compensated(int n);
		T get_image_compensated(int x, int y);

		// Motion compensation methods
		void create_image_compensated(const ImgVector<bool>* mask = nullptr);
		void create_image_compensated_forward(const ImgVector<bool>* mask = nullptr);
		void create_image_estimated(const double estimate_time, const ImgVector<bool>* mask = nullptr);
};

#include "MotionCompensation_private.h"

#endif

