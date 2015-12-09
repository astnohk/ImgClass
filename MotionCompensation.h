#ifndef LIB_ImgClass_MotionCompensation
#define LIB_ImgClass_MotionCompensation

#include "ImgClass.h"
#include "Vector.h"

class MotionCompensation
{
	private:
		bool motion_compensated;
		int _width;
		int _height;
		ImgVector<double> _image_prev;
		ImgVector<double> _image_next;
		ImgVector<double> _image_compensated;
		ImgVector<VECTOR_2D<double> > _vector;
	public:
		MotionCompensation(void);
		explicit MotionCompensation(const MotionCompensation &copy); // copy constructor
		MotionCompensation(int width, int height, const double *image_prev, const double *image_next, const VECTOR_2D<double> *vector);
		MotionCompensation(int width, int height, const double *image_prev, const double *image_next, int v_width, int v_height, const VECTOR_2D<double> *vector);
		MotionCompensation(const ImgVector<double> &image_prev, const ImgVector<double> &image_next, const ImgVector<VECTOR_2D<double> > &vector);
		MotionCompensation(const ImgVector<double> *image_prev, const ImgVector<double> *image_next, const ImgVector<VECTOR_2D<double> > *vector);
		~MotionCompensation(void);

		MotionCompensation& copy(const MotionCompensation &copy);
		MotionCompensation& set(int width, int height, const double *image_prev, const double *image_next, const VECTOR_2D<double> *vector);
		MotionCompensation& set(int width, int height, const double *image_prev, const double *image_next, int v_width, int v_height, const VECTOR_2D<double> *vector);
		MotionCompensation& set(const ImgVector<double> &image_prev, const ImgVector<double> &image_next, const ImgVector<VECTOR_2D<double> > &vector);
		MotionCompensation& set(const ImgVector<double> *image_prev, const ImgVector<double> *image_next, const ImgVector<VECTOR_2D<double> > *vector);

		// Accessor
		// * parameter
		int width(void) const;
		int height(void) const;

		// * reference
		const ImgVector<VECTOR_2D<double> >& ref_vector(void);
		const ImgVector<double>& ref_image_compensated(void);

		double& at_image_compensated(int x, int y);
		double& operator[](int n); // Get reference to compensated_image

		// * original image intensity
		double get_image_prev(int n) const;
		double get_image_prev(int x, int y) const;
		double get_image_next(int n) const;
		double get_image_next(int x, int y) const;
		// * vector
		// returns with qualifier const to avoid to mistake get_vector() for at_vector()
		const VECTOR_2D<double> get_vector(int n) const;
		const VECTOR_2D<double> get_vector(int x, int y) const;
		// * compensated image
		double get_image_compensated(int n);
		double get_image_compensated(int x, int y);

		// Motion compensation methods
		void create_image_compensated(ImgVector<bool> *mask = nullptr);
		void create_image_compensated_forward(ImgVector<bool> *mask = nullptr);
		void create_image_estimated(double estimate_frame, ImgVector<bool> *mask = nullptr);
};

#endif

