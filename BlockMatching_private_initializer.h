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
}




template <class T>
BlockMatching<T>::BlockMatching(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const int BlockSize)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0;
	if (image_prev.isNULL()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : ImgVector<T>& image_prev" << std::endl;
		throw std::invalid_argument("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : ImgVector<T>& image_prev");
	} else if (image_next.isNULL()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : ImgVector<T>& image_next" << std::endl;
		throw std::invalid_argument("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : ImgVector<T>& image_next");
	} else if (image_prev.width() != image_next.width()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : width of image_prev and image_next not match" << std::endl;
		throw std::invalid_argument("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : width of image_prev and image_next not match");
	} else if (image_prev.height() != image_next.height()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : height of image_prev and image_next not match" << std::endl;
		throw std::invalid_argument("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : height of image_prev and image_next not match");
	} else if (BlockSize < 0) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : BlockSize" << std::endl;
		throw std::out_of_range("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : BlockSize");
	}

	_width = image_prev.width();
	_height = image_prev.height();
	_block_size = BlockSize;
	_image_prev.copy(image_prev);
	_image_next.copy(image_next);
	_motion_vector.clear();
	_connected_regions.clear();

	_cells_width = (int)ceil((double)_width / BlockSize);
	_cells_height = (int)ceil((double)_height / BlockSize);
}

template <class T>
BlockMatching<T>::BlockMatching(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const ImgVector<int>& region_map)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0;
	if (image_prev.isNULL()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : const ImgVector<T>& image_prev" << std::endl;
		throw std::invalid_argument("const ImgVector<T>& image_prev");
	} else if (image_next.isNULL()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : const ImgVector<T>& image_next" << std::endl;
		throw std::invalid_argument("const ImgVector<T>& image_next");
	} else if (image_prev.width() != image_next.width() || image_prev.height() != image_next.height()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : width or height of image_prev and image_next not match" << std::endl;
		throw std::invalid_argument("width or height of image_prev and image_next not match");
	} else if (region_map.isNULL()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : const ImgVector<int>& region_map" << std::endl;
		throw std::invalid_argument("const ImgVector<int>& region_map");
	} else if (region_map.width() != image_next.width() || region_map.height() != image_next.height()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : width or height of region_map and image_next not match" << std::endl;
		throw std::invalid_argument("width or height of region_map and image_next not match");
	}

	_width = image_prev.width();
	_height = image_prev.height();
	_block_size = 1;
	_cells_width = _width;
	_cells_height = _height;
	_image_prev.copy(image_prev);
	_image_next.copy(image_next);
	_motion_vector.clear();
	// Extract connected regions from region_map
	get_connected_region_list(region_map);
}




template <class T>
BlockMatching<T>::BlockMatching(const BlockMatching& copy)
{
	_width = copy._width;
	_height = copy._height;
	_block_size = copy._block_size;
	_cells_width = copy._cells_width;
	_cells_height = copy._cells_height;

	_image_prev.copy(copy._image_prev);
	_image_next.copy(copy._image_next);
	_motion_vector.copy(copy._motion_vector);
	_connected_regions.assign(copy._connected_regions.begin(), copy._connected_regions.end());
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
	if (image_prev.isNULL()) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const int) : const ImgVector<T>& image_prev" << std::endl;
		throw std::invalid_argument("const ImgVector<T>& image_prev");
	} else if (image_prev.width() != image_next.width() || image_prev.height() != image_next.height()) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const int) : const ImgVector<T>& image_prev, const ImgVector<T>& image_prev" << std::endl;
		throw std::invalid_argument("width or height of image_prev and image_next not match");
	} else if (BlockSize < 0) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const int) : const int BlockSize" << std::endl;
		throw std::out_of_range("BlockSize");
	}

	_width = image_prev.width();
	_height = image_prev.height();
	_block_size = BlockSize;
	_image_prev.copy(image_prev);
	_image_next.copy(image_next);
	_motion_vector.clear();
	_connected_regions.clear();

	_cells_width = (int)ceil((double)_width / BlockSize);
	_cells_height = (int)ceil((double)_height / BlockSize);
}

template <class T>
void
BlockMatching<T>::reset(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const ImgVector<int>& region_map)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0;
	if (image_prev.isNULL()) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const ImgVector<int>&) : const ImgVector<T>& image_prev" << std::endl;
		throw std::invalid_argument("const ImgVector<T>& image_prev");
	} else if (image_prev.width() != image_next.width() || image_prev.height() != image_next.height()) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const ImgVector<int>&) : const ImgVector<T>& image_prev, const ImgVector<T>& image_next" << std::endl;
		throw std::invalid_argument("width or height of image_prev and image_next not match");
	} else if (region_map.width() != image_next.width() || region_map.height() != image_next.height()) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const ImgVector<int>&) : const ImgVector<int>& region_map" << std::endl;
		throw std::invalid_argument("width or height of region_map not match with image_next");
	}

	_width = image_prev.width();
	_height = image_prev.height();
	_block_size = 1;
	_cells_width = _width;
	_cells_height = _height;
	_image_prev.copy(image_prev);
	_image_next.copy(image_next);
	_motion_vector.clear();
	// Extract connected regions from region_map
	get_connected_region_list(region_map);
}




/* Get connected region list
 *
 * If region_map[i] < 0 then it means the pixel is neglected.
 * So usually region_map include only the integer n > 0.
 */
template <class T>
void
BlockMatching<T>::get_connected_region_list(ImgVector<int> region_map) // get copy of region_map
{
	std::list<std::list<VECTOR_2D<int> > > tmp_list;
	// Clear the vector
	_connected_regions.clear();

	for (int y = 0; y < region_map.height(); y++) {
		for (int x = 0; x < region_map.width(); x++) {
			if (region_map.get(x, y) > 0) {
				int num = region_map.get(x, y);
				region_map.at(x, y) = 0;
				tmp_list.push_back(std::list<VECTOR_2D<int> >(0)); // Add new region pixel list
				VECTOR_2D<int> r(x, y);
				tmp_list.back().push_back(r); // Add first element
				for (std::list<VECTOR_2D<int> >::const_iterator ite = tmp_list.back().begin();
				    ite != tmp_list.back().end();
				    ++ite) {
					if (region_map.get_zeropad(ite->x + 1, ite->y) == num) {
						region_map.at(ite->x + 1, ite->y) = 0;
						r.x = ite->x + 1;
						r.y = ite->y;
						tmp_list.back().push_back(r);
					}
					if (region_map.get_zeropad(ite->x, ite->y + 1) == num) {
						region_map.at(ite->x, ite->y + 1) = 0;
						r.x = ite->x;
						r.y = ite->y + 1;
						tmp_list.back().push_back(r);
					}
					if (region_map.get_zeropad(ite->x - 1, ite->y) == num) {
						region_map.at(ite->x - 1, ite->y) = 0;
						r.x = ite->x - 1;
						r.y = ite->y;
						tmp_list.back().push_back(r);
					}
					if (region_map.get_zeropad(ite->x, ite->y + 1) == num) {
						region_map.at(ite->x, ite->y + 1) = 0;
						r.x = ite->x;
						r.y = ite->y + 1;
						tmp_list.back().push_back(r);
					}
				}
			}
		}
	}
	_connected_regions.resize(tmp_list.size());
	// Copy extracted connected region to std::vector _connected_regions
	std::list<std::list<VECTOR_2D<int> > >::iterator ite = tmp_list.begin();
	for (unsigned int n = 0; n < _connected_regions.size(); ++ite, n++) {
		_connected_regions[n].assign(ite->begin(), ite->end());
	}
}

