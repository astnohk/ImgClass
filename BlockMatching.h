#ifndef nullptr
#define nullptr 0
#endif


#ifndef LIB_ImgClass_BlockMatching
#define LIB_ImgClass_BlockMatching

#include "ImgClass.h"
#include "Vector.h"


template <class T>
class BlockMatching
{
	private:
		ImgVector<T> _image_prev;
		ImgVector<T> _image_next;
		ImgVector<VECTOR_2D<double> > _motion_vector;
		int _width;
		int _height;
		int _block_size;
		int _cells_width;
		int _cells_height;
		bool _MV_forward; // Default value is 'false'
		// If it is 'true' it's mean MV displays flow of image
		// while MV points the pixel on the previous frame from where the next frame's pixel coming when it is 'false'

	public:
		// Constructors
		BlockMatching(void);
		BlockMatching(const BlockMatching& copy);
		BlockMatching(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const int BlockSize);
		BlockMatching(const ImgVector<T>* image_prev, const ImgVector<T>* image_next, const int BlockSize);
		virtual ~BlockMatching(void);

		void reset(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const int BlockSize);
		void reset(const ImgVector<T>* image_prev, const ImgVector<T>* image_next, const int BlockSize);

		// Get state
		int width(void) const;
		int height(void) const;
		int block_size(void) const;
		int vector_width(void) const;
		int vector_height(void) const;
		bool isNULL(void);
		bool isForward(void);

		// Get reference
		ImgVector<VECTOR_2D<double> >& data(void);
		VECTOR_2D<double>& operator[](int n);
		VECTOR_2D<double>& ref(int x, int y);

		// Get data
		VECTOR_2D<double> get(int x, int y); // NOT const because it will make new motion vector when it didn't do block matching

		T MAD(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height, const ImgVector<T>& img_prev, const ImgVector<T>& img_next);
		// Block Matching methods
		// Search in the range of [-floor(search_range / 2), floor(search_range / 2)]
		void block_matching(const int search_range = 41);
		ImgVector<VECTOR_2D<double> >* grad_prev(const int top_left_x, const int top_left_y, const int crop_width, const int crop_height);
		void block_matching_forward(const int search_range = 41); // search_range < 0 then do full search

		void block_matching_subset(const int search_range = 41);

		T SAD(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height);
		T MAD(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height);
		T NCC(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height);
		T ZNCC(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height);
		// Masked one
		T MAD(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height, const ImgVector<bool>& mask);
		T ZNCC(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height, const ImgVector<bool>& mask);
};

#include "BlockMatching_private_initializer.h"
#include "BlockMatching_private_accessor.h"
#include "BlockMatching_private_main.h"

#endif

