#include <algorithm>
#include <cmath>
#include <iostream>
#include <new>
#include <stdexcept>




template <class T>
BlockMatching<T>::BlockMatching(void)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0;
	_MV_forward = false;
}

template <class T>
BlockMatching<T>::BlockMatching(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const int BlockSize)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0;
	_MV_forward = false;
	if (image_prev.isNULL()) {
		throw std::invalid_argument("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : ImgVector<T>& image_prev");
	} else if (image_next.isNULL()) {
		throw std::invalid_argument("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : ImgVector<T>& image_next");
	} else if (image_prev.width() != image_next.width()) {
		throw std::invalid_argument("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : width of image_prev and image_next not match");
	} else if (image_prev.height() != image_next.height()) {
		throw std::invalid_argument("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : height of image_prev and image_next not match");
	} else if (BlockSize < 0) {
		throw std::out_of_range("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : BlockSize");
	}

	_width = image_prev.width();
	_height = image_prev.height();
	_block_size = BlockSize;
	_image_prev.copy(image_prev);
	_image_next.copy(image_next);
	_motion_vector.reset(0, 0);

	_cells_width = (int)ceil((double)_width / BlockSize);
	_cells_height = (int)ceil((double)_height / BlockSize);
}

template <class T>
BlockMatching<T>::BlockMatching(const ImgVector<T>* image_prev, const ImgVector<T>* image_next, const int BlockSize)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0;
	_MV_forward = false;
	if (image_prev == nullptr) {
		throw std::invalid_argument("const ImgVector<T>* image_prev");
	} else if (image_next == nullptr) {
		throw std::invalid_argument("const ImgVector<T>* image_next");
	} if (image_prev->isNULL()) {
		throw std::invalid_argument("const ImgVector<T>* image_prev");
	} if (image_next->isNULL()) {
		throw std::invalid_argument("const ImgVector<T>* image_next");
	} else if (image_prev.width() != image_next.width() || image_prev.height() != image_next.height()) {
		throw std::invalid_argument("width or height of image_prev and image_next not match");
	} else if (BlockSize < 0) {
		throw std::out_of_range("BlockSize");
	}

	_width = image_prev.width();
	_height = image_prev.height();
	_block_size = BlockSize;
	_image_prev.copy(image_prev);
	_image_next.copy(image_next);
	_motion_vector.reset(0, 0);

	_cells_width = (int)ceil((double)_width / BlockSize);
	_cells_height = (int)ceil((double)_height / BlockSize);
}

template <class T>
BlockMatching<T>::BlockMatching(const BlockMatching& copy)
{
	_width = copy._width;
	_height = copy._height;
	_block_size = copy._block_size;
	_cells_width = copy._cells_width;
	_cells_height = copy._cells_height;
	_MV_forward = copy._MV_forward;

	_image_prev.copy(copy._image_prev);
	_image_next.copy(copy._image_next);
	_motion_vector.copy(copy._motion_vector);
}

template <class T>
BlockMatching<T>::~BlockMatching(void)
{
}

template <class T>
void
BlockMatching<T>::reset(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const int BlockSize)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0;
	_MV_forward = false;
	if (image_prev.isNULL()) {
		throw std::invalid_argument("const ImgVector<T>* image_prev");
	} else if (image_next.isNULL()) {
		throw std::invalid_argument("const ImgVector<T>* image_next");
	} else if (image_prev.width() != image_next.width() || image_prev.height() != image_next.height()) {
		throw std::invalid_argument("width or height of image_prev and image_next not match");
	} else if (BlockSize < 0) {
		throw std::out_of_range("BlockSize");
	}

	_width = image_prev.width();
	_height = image_prev.height();
	_block_size = BlockSize;
	_image_prev.copy(image_prev);
	_image_next.copy(image_next);
	_motion_vector.reset(0, 0);

	_cells_width = (int)ceil((double)_width / BlockSize);
	_cells_height = (int)ceil((double)_height / BlockSize);
}

template <class T>
void
BlockMatching<T>::reset(const ImgVector<T>* image_prev, const ImgVector<T>* image_next, const int BlockSize)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0;
	_MV_forward = false;
	if (image_prev == nullptr) {
		throw std::invalid_argument("void BlockMatching<T>::reset(const ImgVector<T>*, const ImgVector<T>*, const int) : const ImgVector<T>* image_prev");
	} else if (image_next == nullptr) {
		throw std::invalid_argument("void BlockMatching<T>::reset(const ImgVector<T>*, const ImgVector<T>*, const int) : const ImgVector<T>* image_next");
	} else if (image_prev->isNULL()) {
		throw std::invalid_argument("void BlockMatching<T>::reset(const ImgVector<T>*, const ImgVector<T>*, const int) : const ImgVector<T>* image_prev");
	} else if (image_next->isNULL()) {
		throw std::invalid_argument("void BlockMatching<T>::reset(const ImgVector<T>*, const ImgVector<T>*, const int) : const ImgVector<T>* image_next");
	} else if (image_prev->width() != image_next->width()) {
		throw std::invalid_argument("void BlockMatching<T>::reset(const ImgVector<T>*, const ImgVector<T>*, const int) : the width of image_prev and image_next not match");
	} else if (image_prev->height() != image_next->height()) {
		throw std::invalid_argument("void BlockMatching<T>::reset(const ImgVector<T>*, const ImgVector<T>*, const int) : the height of image_prev and image_next not match");
	} else if (BlockSize < 0) {
		throw std::out_of_range("void BlockMatching<T>::reset(const ImgVector<T>*, const ImgVector<T>*, const int) : BlockSize");
	}

	_width = image_prev->width();
	_height = image_prev->height();
	_block_size = BlockSize;
	_image_prev.copy(image_prev);
	_image_next.copy(image_next);
	_motion_vector.reset(0, 0);

	_cells_width = (int)ceil((double)_width / BlockSize);
	_cells_height = (int)ceil((double)_height / BlockSize);
}

