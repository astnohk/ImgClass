#include "ImgClass.h"
#include "Vector.h"



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
		ImgVector<int> _segments;

	public:
		Segmentation(void);
		~Segmentation(void);

		// Accessor
		int width(void) const;
		int height(void) const;

		int& operator[](int n);
		int& at(int x, int y);
		int& at_repeat(int x, int y);
		int& at_mirror(int x, int y);

		int get(int n) const;
		int get(int x, int y) const;
		int get_zeropad(int x, int y) const;
		int get_repeat(int x, int y) const;
		int get_mirror(int x, int y) const;

		// Mean Shift segmentation
		void Segmentation_MeanShift_Grayscale(int Iter_Max = 64);
		VECTOR_2D<double> MeanShift_Grayscale(const int x, const int y, int Iter_Max = 64);
}

#include "Segmentation_private.h"

#endif

