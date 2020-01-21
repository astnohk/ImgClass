#ifndef LIB_ImgClass_Segmentation
#define LIB_ImgClass_Segmentation

#include <list>
#include <vector>

#include "Vector.h"
#include "ImgClass.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace ImgClass {
	class RGB;
	class Lab;

	template <class T>
	class Segmentation
	{
		private:

		struct tuple
		{
			VECTOR_2D<double> spatial;
			T color;
		};
		struct Region
		{
			std::list<VECTOR_2D<int> > points;
			T color;
		};

		int _width;
		int _height;
		size_t _size;
		size_t _min_pixels;
		double _kernel_spatial;
		double _kernel_intensity;
		ImgVector<T> _image;
		ImgVector<T> _color_quantized_image;
		ImgVector<std::list<VECTOR_2D<int> > > _vector_converge_list_map;
		ImgVector<VECTOR_2D<double> > _shift_vector_spatial;
		ImgVector<T> _shift_vector_color;
		ImgVector<size_t> _segmentation_map;
		std::vector<std::vector<VECTOR_2D<int> > > _regions;


		public:

		// Constructor
		Segmentation(void);
		Segmentation(const ImgVector<T>& image, const double &kernel_spatial_radius = 16.0, const double &kernel_intensity_radius = 10.0 / 255.0, const size_t &min_number_of_pixels = 4);
		Segmentation(const Segmentation<T>& segmentation); // Copy constructor

		Segmentation<T>& reset(const ImgVector<T> &image, const int IterMax, const double &kernel_spatial_radius = 16.0, const double &kernel_intensity_radius = 10.0 / 255.0, const size_t &min_number_of_pixels = 4);

		Segmentation<T>& copy(const Segmentation<T> &segmentation);

		// Destructor
		~Segmentation(void);

		// Setter
		void set_kernel(const double &kernel_spatial_radius, const double &kernel_intensity_radius);
		void set_min_pixels(const size_t &min_number_of_pixels);

		Segmentation<T>& operator=(const Segmentation<T>& rvalue);

		// Accessor
		int width(void) const;
		int height(void) const;
		size_t size(void) const;

		const ImgVector<T>& ref_color_quantized_image(void) const;
		const ImgVector<size_t>& ref_segmentation_map(void) const;
		const ImgVector<std::list<VECTOR_2D<int> > >& ref_vector_converge_list_map(void) const;
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
		void Segmentation_MeanShift(const int Iter_Max = 32);


		protected:

		const ImgClass::Segmentation<T>::tuple MeanShift(const int x, const int y, std::vector<VECTOR_2D<int> >& pel_list, int Iter_Max);
		size_t collect_regions_in_segmentation_map(std::list<std::list<VECTOR_2D<int> > >* regions_list);
		size_t small_region_concatenator(std::list<std::list<VECTOR_2D<int> > >* regions_list);
		void small_region_eliminator(std::list<std::list<VECTOR_2D<int> > >* regions_list);

		double distance(const T& lcolor, const T& rcolor); // Calculate distance depends on each color space
		double normalized_distance(const T& lcolor, const T& rcolor); // Calculate distance depends on each color space
	};
}


double inner_prod(const double& lvalue, const double& rvalue);

double norm(const double& value);
double norm_squared(const double& value);


#include "Segmentation_private.h"

#endif

