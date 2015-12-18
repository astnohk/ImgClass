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
		ImgVector<T> _image_next;
		ImgVector<T> _image_compensated;
		ImgVector<VECTOR_2D<double> > _vector;
	public:
		MotionCompensation(void);
		explicit MotionCompensation(const MotionCompensation& copy); // copy constructor
		MotionCompensation(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const ImgVector<VECTOR_2D<double> >& vector);
		~MotionCompensation(void);

		MotionCompensation& copy(const MotionCompensation& copy);
		MotionCompensation& set(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const ImgVector<VECTOR_2D<double> >& vector);

		// Accessor
		// * parameter
		int width(void) const;
		int height(void) const;

		// * reference
		const ImgVector<VECTOR_2D<double> >& ref_vector(void);
		const ImgVector<T>& ref_image_compensated(void);

		T& at_image_compensated(int x, int y);
		T& operator[](int n); // Get reference to compensated_image

		// * original image intensity
		T get_image_prev(int n) const;
		T get_image_prev(int x, int y) const;
		T get_image_next(int n) const;
		T get_image_next(int x, int y) const;
		// * vector
		// returns with qualifier const to avoid to mistake get_vector() for at_vector()
		const VECTOR_2D<double> get_vector(int n) const;
		const VECTOR_2D<double> get_vector(int x, int y) const;
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

