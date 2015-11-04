#ifndef LIB_BlockMatching
#define LIB_BlockMatching

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

	public:
		// Constructors
		BlockMatching(void);
		BlockMatching(const BlockMatching& copy);
		BlockMatching(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const int BlockSize);
		BlockMatching(const ImgVector<T>* image_prev, const ImgVector<T>* image_next, const int BlockSize);
		virtual ~BlockMatching(void);

		void reset(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const int BlockSize);
		void reset(const ImgVector<T>* image_prev, const ImgVector<T>* image_next, const int BlockSize);

		// Main functions
		void block_matching(void);
		VECTOR_2D<double> max_crosscorr(const int x_prev, const int y_prev);
		T MAD(const int x_prev, const int y_prev, const int x_next, const int y_next, const int block_width, const int block_height);

		// Get state
		int width(void) const;
		int height(void) const;
		int vector_width(void) const;
		int vector_height(void) const;
		bool isNULL(void);
		// Get reference
		ImgVector<VECTOR_2D<double> >& data(void);
		VECTOR_2D<double>& operator[](int n);
		VECTOR_2D<double>& ref(int x, int y);
		// Get data
		VECTOR_2D<double> get(int x, int y) const;
};

#include "BlockMatching_private.h"

#endif

