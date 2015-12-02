#include "ImgClass.h"
#include "Vector.h"

#if defined(_OPENMP)
#include <omp.h>
#endif



#ifndef LIB_ImgClass_Segmentation
#define LIB_ImgClass_Segmentation

template <class T>
class Segmentation
{
	private:
		int _width;
		int _height;
		double _kernel_spatial;
		double _kernel_intensity;
		ImgVector<T> _image;
		ImgVector<VECTOR_2D<double> > _shift_vector;
		ImgVector<int> _segments;

	public:
		// Constructor
		Segmentation(void);
		explicit Segmentation(const Segmentation<T>& segments); // Copy constructor
		Segmentation(const ImgVector<T>* image, const double kernel_spatial_radius = 16.0, const double kernel_intensity_radius = 0.2);
		Segmentation<T>& reset(const ImgVector<T>* image, const double kernel_spatial_radius = 16.0, const double kernel_intensity_radius = 0.2);
		Segmentation<T>& copy(const Segmentation<T>* segments);
		Segmentation<T>& copy(const Segmentation<T>& segments);
		// Destructor
		~Segmentation(void);

		// Accessor
		int width(void) const;
		int height(void) const;

		const ImgVector<int>& ref_segments(void) const;

		int& operator[](int n);
		int& at(int n);
		int& at(int x, int y);
		int& at_repeat(int x, int y);
		int& at_mirror(int x, int y);

		int get(int n) const;
		int get(int x, int y) const;
		int get_zeropad(int x, int y) const;
		int get_repeat(int x, int y) const;
		int get_mirror(int x, int y) const;

		// Mean Shift segmentation
		void Segmentation_MeanShift(int Iter_Max = 64);
		const VECTOR_2D<double> MeanShift_Grayscale(const int x, const int y, int Iter_Max);
};

#include "Segmentation_private.h"

#endif

