#ifndef LIB_ImgClass_BlockMatching
#define LIB_ImgClass_BlockMatching

/* Macro for compatibility where the C++11 is not supported
#ifndef nullptr
#define nullptr 0
#endif
*/

#include <list>
#include <vector>

#include "Vector.h"
#include "ImgClass.h"

namespace ImgClass {
	class RGB;
	class Lab;
}

template <class T>
class BlockMatching
{
	private:
		int _width;
		int _height;
		int _block_size;
		int _cells_width;
		int _cells_height;
		int _subpixel_scale; // Full-pel : 1, Half-pel : 2, Quarter-pel : 4
		ImgVector<T> _image_prev;
		ImgVector<T> _image_current; // Base image for motion estimation
		ImgVector<T> _image_next;
		ImgVector<size_t> _region_map_prev;
		ImgVector<size_t> _region_map_current;
		ImgVector<size_t> _region_map_next;
		ImgVector<T> _color_quantized_prev;
		ImgVector<T> _color_quantized_current;
		ImgVector<T> _color_quantized_next;
		ImgVector<Vector_ST<double> > _motion_vector_time;
		ImgVector<VECTOR_2D<double> > _motion_vector_prev;
		ImgVector<VECTOR_2D<double> > _motion_vector_next;
		// For arbitrary shaped block matching
		std::vector<std::list<VECTOR_2D<int> > > _connected_regions_prev;
		std::vector<std::list<VECTOR_2D<int> > > _connected_regions_current;
		std::vector<std::list<VECTOR_2D<int> > > _connected_regions_next;

	public:
		// Constructors
		BlockMatching(void);
		BlockMatching(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const int BlockSize, const int Subpixel_Scale = 1);
		BlockMatching(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<T>& image_next, const int BlockSize, const int Subpixel_Scale = 1);
		BlockMatching(const ImgVector<T>& image_prev, const ImgVector<size_t>& region_prev, const ImgVector<T>& image_current, const ImgVector<size_t>& region_current, const int Subpixel_Scale = 1);
		BlockMatching(const ImgVector<T>& image_prev, const ImgVector<size_t>& region_prev, const ImgVector<T>& image_current, const ImgVector<size_t>& region_current, const ImgVector<T>& image_next, const ImgVector<size_t>& region_next, const int Subpixel_Scale = 1);
		// Copy constructor
		explicit BlockMatching(const BlockMatching& copy);
		// Destructor
		virtual ~BlockMatching(void);

		// Resetter
		void reset(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const int BlockSize, const int Subpixel_Scale = 1);
		void reset(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<T>& image_next, const int BlockSize, const int Subpixel_Scale = 1);
		void reset(const ImgVector<T>& image_prev, const ImgVector<size_t>& region_prev, const ImgVector<T>& image_current, const ImgVector<size_t>& region_current, const int Subpixel_Scale = 1);
		void reset(const ImgVector<T>& image_prev, const ImgVector<size_t>& region_prev, const ImgVector<T>& image_current, const ImgVector<size_t>& region_current, const ImgVector<T>& image_next, const ImgVector<size_t>& region_next, const int Subpixel_Scale = 1);

		// Get state
		int width(void) const;
		int height(void) const;
		int block_size(void) const;
		int vector_field_width(void) const;
		int vector_field_height(void) const;
		bool isNULL(void) const;

		// Get reference
		ImgVector<Vector_ST<double> >& ref_motion_vector_time(void);
		ImgVector<VECTOR_2D<double> >& ref_motion_vector_prev(void);
		ImgVector<VECTOR_2D<double> >& ref_motion_vector_next(void);
		VECTOR_2D<double>& operator[](int n);
		Vector_ST<double>& at(int x, int y);
		VECTOR_2D<double>& at_prev(int x, int y);
		VECTOR_2D<double>& at_next(int x, int y);
		Vector_ST<double>& at_block(int x, int y);
		VECTOR_2D<double>& at_block_prev(int x, int y);
		VECTOR_2D<double>& at_block_next(int x, int y);

		// Get data
		// returns const to avoid to mistake get() for at()
		const Vector_ST<double> get(int x, int y); // NOT const method because it will make new motion vector when it haven't done block matching
		const VECTOR_2D<double> get_prev(int x, int y); // NOT const method because it will make new motion vector when it haven't done block matching
		const VECTOR_2D<double> get_next(int x, int y); // NOT const method because it will make new motion vector when it haven't done block matching
		const Vector_ST<double> get_block(int x, int y); // NOT const method because it will make new motion vector when it haven't done block matching
		const VECTOR_2D<double> get_block_prev(int x, int y); // NOT const method because it will make new motion vector when it haven't done block matching
		const VECTOR_2D<double> get_block_next(int x, int y); // NOT const method because it will make new motion vector when it haven't done block matching

		// Block Matching methods
		// Search in the range of [-floor(search_range / 2), floor(search_range / 2)]
		void block_matching(const int search_range = 41, const double coeff_MAD = 1.0, const double coeff_ZNCC = 0.0);

	protected:
		void image_normalizer(void);
		// Extract connected region from region_map
		void get_connected_region_list(std::vector<std::list<VECTOR_2D<int> > >* connected_regions, const ImgVector<size_t>& region_map);
		void get_color_quantized_image(ImgVector<T>* decreased_color_image, const ImgVector<T>& image, const std::vector<std::list<VECTOR_2D<int> > >& connected_regions);

		// Main method of block_matching
		void block_matching_lattice(const int search_range, const double coeff_MAD, const double coeff_ZNCC);
		void block_matching_arbitrary_shaped(const int search_range, const double coeff_MAD, const double coeff_ZNCC);
		// Interpolate skipped Motion Vectors
		void vector_interpolation(const std::list<VECTOR_2D<int> >& flat_blocks, ImgVector<bool>* estimated);

		// Correlation function
		double MAD(const ImgVector<T>& reference, const ImgVector<T>& interest, const int x_ref, const int y_ref, const int x_int, const int y_int);
		double MAD_cubic(const ImgVector<T>& reference, const ImgVector<T>& interest, const double x_ref, const double y_ref, const double x_int, const double y_int);
		double ZNCC(const ImgVector<T>& reference, const ImgVector<T>& interest, const int x_ref, const int y_ref, const int x_int, const int y_int);
		// Arbitrary shaped correlation function
		double MAD_region(const ImgVector<T>& reference, const ImgVector<T>& interest, const int x_diff, const int y_diff, const std::list<VECTOR_2D<int> >& region_interest);
		double MAD_region_cubic(const ImgVector<T>& reference, const ImgVector<T>& interest, const double x_diff, const double y_diff, const std::list<VECTOR_2D<int> >& region_interest);
		double ZNCC_region(const ImgVector<T>& reference, const ImgVector<T>& interest, const int x_diff, const int y_diff, const std::list<VECTOR_2D<int> >& region_interest);
		// Arbitrary shaped correlation function with nearest intensity restricted
		double MAD_region_nearest_intensity(const int x_diff, const int y_diff, const std::list<VECTOR_2D<int> >& region_current);
		double ZNCC_region_nearest_intensity(const int x_diff, const int y_diff, const std::list<VECTOR_2D<int> >& region_current);
};

double norm_squared(const double& value);
double norm(const double& value);
double inner_prod(const double& lvalue, const double& rvalue);

#include "BlockMatching_private_initializer.h"
#include "BlockMatching_private_accessor.h"
#include "BlockMatching_private_main.h"

#endif

