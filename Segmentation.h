#ifndef LIB_ImgClass_Segmentation
#define LIB_ImgClass_Segmentation

#include "ImgClass.h"
#include "Lab.h"
#include "RGB.h"
#include "Vector.h"

#if defined(_OPENMP)
#include <omp.h>
#endif


template <class T>
class Segmentation
{
	private:
		int _width;
		int _height;
		double _kernel_spatial;
		double _kernel_intensity;
		ImgVector<T> _image;
		ImgVector<int> _decrease_color_image;
		ImgVector<VECTOR_2D<double> > _shift_vector;
		ImgVector<int> _segments_map;
		std::vector<std::vector<VECTOR_2D<int> > > _regions;

	public:
		// Constructor
		Segmentation(void);
		Segmentation(const ImgVector<T>& image, const double kernel_spatial_radius = 8.0, const double kernel_intensity_radius = 12.0 / 255.0);
		explicit Segmentation(const Segmentation<T>& segments); // Copy constructor

		Segmentation<T>& reset(const ImgVector<T>& image, const double kernel_spatial_radius = 8.0, const double kernel_intensity_radius = 12.0 / 255.0);

		Segmentation<T>& copy(const Segmentation<T>& segments);

		// Destructor
		~Segmentation(void);

		// Accessor
		int width(void) const;
		int height(void) const;

		const ImgVector<int>& ref_decrease_color_image(void) const;
		const ImgVector<int>& ref_segments_map(void) const;
		const ImgVector<VECTOR_2D<double> >& ref_shift_vector(void) const;
		const std::vector<std::list<VECTOR_2D<int> > >& ref_regions(void) const;

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
		void Segmentation_MeanShift(const int Iter_Max = 64, const unsigned int Min_Number_of_Pixels = 25, const int Search_Range = 4);

	protected:
		const VECTOR_2D<double> MeanShift(const int x, const int y, std::vector<VECTOR_2D<int> >& pel_list, int Iter_Max);
		unsigned int small_region_eliminate(std::vector<std::list<VECTOR_2D<int> > >* regions_vector, const unsigned int Min_Number_of_Pixels, const int search_range);
		double mean_kernel(const int x, const int y, std::vector<VECTOR_2D<int> >& pel_list);

		double distance(const T& lcolor, const T& rcolor); // Calculate distance depends on each color space
};

#include "Segmentation_private.h"

#endif

