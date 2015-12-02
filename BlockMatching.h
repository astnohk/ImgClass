#ifndef nullptr
#define nullptr 0
#endif


#ifndef LIB_ImgClass_BlockMatching
#define LIB_ImgClass_BlockMatching

#include <list>
#include <vector>

#include "ImgClass.h"
#include "Vector.h"


template <class T>
class BlockMatching
{
	private:
		int _width;
		int _height;
		int _block_size;
		int _cells_width;
		int _cells_height;
		ImgVector<T> _image_prev;
		ImgVector<T> _image_next;
		ImgVector<VECTOR_2D<double> > _motion_vector;
		// For arbitrary shaped block matching
		std::vector<std::list<VECTOR_2D<int> > > _connected_regions;

	public:
		// Constructors
		BlockMatching(void);
		BlockMatching(const BlockMatching& copy);
		BlockMatching(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const int BlockSize);
		BlockMatching(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const ImgVector<int>& region_map);
		virtual ~BlockMatching(void);

		void reset(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const int BlockSize);
		void reset(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const ImgVector<int>& region_map);

		void get_connected_region_list(ImgVector<int> region_map); // get copy of region_map to modify

		// Get state
		int width(void) const;
		int height(void) const;
		int block_size(void) const;
		int vector_width(void) const;
		int vector_height(void) const;
		bool isNULL(void);

		// Get reference
		ImgVector<VECTOR_2D<double> >& data(void);
		VECTOR_2D<double>& operator[](int n);
		VECTOR_2D<double>& at(int x, int y);

		// Get data
		// returns const to avoid to mistake get() for at()
		const VECTOR_2D<double> get(int x, int y); // NOT const method because it will make new motion vector when it haven't done block matching

		T MAD(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height, const ImgVector<T>& img_prev, const ImgVector<T>& img_next);
		// Block Matching methods
		// Search in the range of [-floor(search_range / 2), floor(search_range / 2)]
		void block_matching(const int search_range = 41);
		void block_matching_subset(const ImgVector<int>* region_map, const int search_range = 41);
		// Interpolate skipped Motion Vectors
		void vector_interpolation(const std::list<VECTOR_2D<int> >& flat_blocks, ImgVector<bool>* estimated);

		ImgVector<VECTOR_2D<double> >* grad_prev(const int top_left_x, const int top_left_y, const int crop_width, const int crop_height);

		const T SAD(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height);
		const T MAD(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height);
		const T NCC(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height);
		const T ZNCC(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height);
		// Masked one
		const T MAD(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height, const ImgVector<bool>& mask);
		const T ZNCC(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height, const ImgVector<bool>& mask);
};

#include "BlockMatching_private_initializer.h"
#include "BlockMatching_private_accessor.h"
#include "BlockMatching_private_main.h"

#endif

