#ifndef LIB_ImgClass_Segmentation
#define LIB_ImgClass_Segmentation

#include <list>
#include <vector>

#include "Color.h"
#include "Vector.h"
#include "ImgClass.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace ImgClass {
	class RGB;
	class Lab;
	namespace Segmentation {
		template <class T>
		struct tuple
		{
			VECTOR_2D<double> spatial;
			T color;
		};
	}
}

template <class T>
class Segmentation
{
	private:
		int _width;
		int _height;
		size_t _size;
		double _kernel_spatial;
		double _kernel_intensity;
		ImgVector<T> _image;
		ImgVector<T> _color_quantized_image;
		ImgVector<VECTOR_2D<double> > _shift_vector_spatial;
		ImgVector<T> _shift_vector_color;
		ImgVector<size_t> _segmentation_map;
		std::vector<std::vector<VECTOR_2D<int> > > _regions;

	public:
		// Constructor
		Segmentation(void);
		Segmentation(const ImgVector<T>& image, const double kernel_spatial_radius = 16.0, const double kernel_intensity_radius = 10.0 / 255.0);
		explicit Segmentation(const Segmentation<T>& segmentation); // Copy constructor

		Segmentation<T>& reset(const ImgVector<T>& image, const double kernel_spatial_radius = 16.0, const double kernel_intensity_radius = 10.0 / 255.0);

		Segmentation<T>& copy(const Segmentation<T>& segmentation);

		// Destructor
		~Segmentation(void);

		// Setter
		void set_kernel(const double kernel_spatial_radius, const double kernel_intensity_radius);

		Segmentation<T>& operator=(const Segmentation<T>& rvalue);

		// Accessor
		int width(void) const;
		int height(void) const;
		size_t size(void) const;

		const ImgVector<T>& ref_color_quantized_image(void) const;
		const ImgVector<size_t>& ref_segmentation_map(void) const;
		const ImgVector<VECTOR_2D<double> >& ref_shift_vector_spatial(void) const;
		const ImgVector<T>& ref_shift_vector_color(void) const;
		const std::vector<std::vector<VECTOR_2D<int> > >& ref_regions(void) const;

		size_t& operator[](size_t n);
		size_t& at(size_t n);
		size_t& at(int x, int y);
		size_t& at_repeat(int x, int y);
		size_t& at_mirror(int x, int y);

		size_t get(size_t n) const;
		size_t get(int x, int y) const;
		size_t get_zeropad(int x, int y) const;
		size_t get_repeat(int x, int y) const;
		size_t get_mirror(int x, int y) const;

		// Mean Shift segmentation
		void Segmentation_MeanShift(const int Iter_Max = 128, const size_t Min_Number_of_Pixels = 10);

	protected:
		const ImgClass::Segmentation::tuple<T> MeanShift(const int x, const int y, std::vector<VECTOR_2D<int> >& pel_list, int Iter_Max);
		size_t collect_regions_in_segmentation_map(std::list<std::list<VECTOR_2D<int> > >* regions_list);
		size_t small_region_eliminate(std::list<std::list<VECTOR_2D<int> > >* regions_list, const size_t Min_Number_of_Pixels);

		double distance(const T& lcolor, const T& rcolor); // Calculate distance depends on each color space
		double normalized_distance(const T& lcolor, const T& rcolor); // Calculate distance depends on each color space
};

#include "Segmentation_private.h"

#endif

