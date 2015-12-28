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
BlockMatching<T>::BlockMatching(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const int BlockSize, const bool dense)
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
	_image_prev.copy(image_prev);
	_image_next.copy(image_next);

	// Normalize the image
	image_normalizer();

	_block_size = BlockSize;
	if (dense) {
		_cells_width = _width;
		_cells_height = _height;
	} else {
		_cells_width = int(ceil(double(_width) / double(BlockSize)));
		_cells_height = int(ceil(double(_height) / double(BlockSize)));
	}
}

template <class T>
BlockMatching<T>::BlockMatching(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const ImgVector<int>& region_map_prev, const ImgVector<int>& region_map_next)
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
	} else if (region_map_prev.isNULL()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : const ImgVector<int>& region_map_prev" << std::endl;
		throw std::invalid_argument("const ImgVector<int>& region_map_prev");
	} else if (region_map_next.isNULL()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : const ImgVector<int>& region_map_next" << std::endl;
		throw std::invalid_argument("const ImgVector<int>& region_map_prev");
	} else if (region_map_prev.width() != image_prev.width() || region_map_prev.height() != image_prev.height()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : width or height of region_map_prev and image_prev not match" << std::endl;
		throw std::invalid_argument("width or height of region_map_prev and image_prev not match");
	} else if (region_map_next.width() != image_next.width() || region_map_next.height() != image_next.height()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : width or height of region_map_next and image_next not match" << std::endl;
		throw std::invalid_argument("width or height of region_map_next and image_next not match");
	}

	_width = image_prev.width();
	_height = image_prev.height();
	_block_size = 1;
	_cells_width = _width;
	_cells_height = _height;
	_image_prev.copy(image_prev);
	_image_next.copy(image_next);
	_region_map_prev.copy(region_map_prev);
	_region_map_next.copy(region_map_next);
	// Normalize the image
	image_normalizer();
	// Extract connected regions from region_map
	get_connected_region_list(&_connected_regions_prev, region_map_prev);
	get_connected_region_list(&_connected_regions_next, region_map_next);
	// Get decreased image
	get_color_quantized_image(&_color_quantized_prev, _image_prev, _connected_regions_prev);
	get_color_quantized_image(&_color_quantized_next, _image_next, _connected_regions_next);
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
	_region_map_prev.copy(copy._region_map_prev);
	_region_map_next.copy(copy._region_map_next);
	_motion_vector.copy(copy._motion_vector);
	_connected_regions_prev.assign(copy._connected_regions_prev.begin(), copy._connected_regions_prev.end());
	_connected_regions_next.assign(copy._connected_regions_next.begin(), copy._connected_regions_next.end());
}




template <class T>
BlockMatching<T>::~BlockMatching(void)
{
}




template <class T>
void
BlockMatching<T>::reset(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const int BlockSize, const bool dense)
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
	_region_map_prev.clear();
	_region_map_next.clear();
	_connected_regions_prev.clear();
	_connected_regions_next.clear();
	_motion_vector.clear();

	// Normalize the image
	image_normalizer();

	if (dense) {
		_cells_width = _width;
		_cells_height = _height;
	} else {
		_cells_width = int(ceil(double(_width) / double(BlockSize)));
		_cells_height = int(ceil(double(_height) / double(BlockSize)));
	}
}

template <class T>
void
BlockMatching<T>::reset(const ImgVector<T>& image_prev, const ImgVector<T>& image_next, const ImgVector<int>& region_map_prev, const ImgVector<int>& region_map_next)
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
	} else if (region_map_prev.width() != image_prev.width() || region_map_prev.height() != image_prev.height()) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const ImgVector<int>&) : const ImgVector<int>& region_map_prev" << std::endl;
		throw std::invalid_argument("width or height of region_map_prev not match with image_prev");
	} else if (region_map_next.width() != image_next.width() || region_map_next.height() != image_next.height()) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const ImgVector<int>&) : const ImgVector<int>& region_map_next" << std::endl;
		throw std::invalid_argument("width or height of region_map_next not match with image_next");
	}

	_width = image_prev.width();
	_height = image_prev.height();
	_block_size = 1;
	_cells_width = _width;
	_cells_height = _height;
	_image_prev.copy(image_prev);
	_image_next.copy(image_next);
	_region_map_prev.copy(region_map_prev);
	_region_map_next.copy(region_map_next);
	_motion_vector.clear();
	// Normalize the image
	image_normalizer();
	// Extract connected regions from region_map
	get_connected_region_list(&_connected_regions_prev, region_map_prev);
	get_connected_region_list(&_connected_regions_next, region_map_next);
	// Get decreased image
	get_color_quantized_image(&_color_quantized_prev, _image_prev, _connected_regions_prev);
	get_color_quantized_image(&_color_quantized_next, _image_next, _connected_regions_next);
}




/* Get connected region list
 *
 * If region_map[i] < 0 then it means the pixel is neglected.
 * So usually region_map include only the integer n > 0.
 */
template <class T>
void
BlockMatching<T>::get_connected_region_list(std::vector<std::list<VECTOR_2D<int> > >* connected_regions, const ImgVector<int>& region_map)
{
	std::list<std::list<VECTOR_2D<int> > > tmp_list;
	ImgVector<int> region_map_tmp(region_map);
	// Clear the vector
	connected_regions->clear();

	for (int y = 0; y < region_map_tmp.height(); y++) {
		for (int x = 0; x < region_map_tmp.width(); x++) {
			if (region_map_tmp.get(x, y) > 0) {
				int num = region_map_tmp.get(x, y);
				region_map_tmp.at(x, y) = 0;
				tmp_list.push_back(std::list<VECTOR_2D<int> >(0)); // Add new region pixel list
				VECTOR_2D<int> r(x, y);
				tmp_list.back().push_back(r); // Add first element
				for (std::list<VECTOR_2D<int> >::const_iterator ite = tmp_list.back().begin();
				    ite != tmp_list.back().end();
				    ++ite) {
					if (region_map_tmp.get_zeropad(ite->x + 1, ite->y) == num) {
						region_map_tmp.at(ite->x + 1, ite->y) = 0;
						r.x = ite->x + 1;
						r.y = ite->y;
						tmp_list.back().push_back(r);
					}
					if (region_map_tmp.get_zeropad(ite->x, ite->y + 1) == num) {
						region_map_tmp.at(ite->x, ite->y + 1) = 0;
						r.x = ite->x;
						r.y = ite->y + 1;
						tmp_list.back().push_back(r);
					}
					if (region_map_tmp.get_zeropad(ite->x - 1, ite->y) == num) {
						region_map_tmp.at(ite->x - 1, ite->y) = 0;
						r.x = ite->x - 1;
						r.y = ite->y;
						tmp_list.back().push_back(r);
					}
					if (region_map_tmp.get_zeropad(ite->x, ite->y + 1) == num) {
						region_map_tmp.at(ite->x, ite->y + 1) = 0;
						r.x = ite->x;
						r.y = ite->y + 1;
						tmp_list.back().push_back(r);
					}
				}
			}
		}
	}
	connected_regions->resize(tmp_list.size());
	// Copy extracted connected region to std::vector _connected_regions
	std::list<std::list<VECTOR_2D<int> > >::iterator ite = tmp_list.begin();
	for (unsigned int n = 0; n < connected_regions->size(); ++ite, n++) {
		connected_regions->at(n).assign(ite->begin(), ite->end());
	}
}


// ----- Normalizer -----
template <class T>
void
BlockMatching<T>::image_normalizer(void)
{
	double max_int = std::max(_image_prev.max(), _image_next.max());
	if (max_int > 1.0) {
		_image_prev /= max_int;
		_image_next /= max_int;
	}
}

template <>
void
BlockMatching<ImgClass::RGB>::image_normalizer(void);

template <>
void
BlockMatching<ImgClass::Lab>::image_normalizer(void);



// ----- Decrease Color -----
template <class T>
void
BlockMatching<T>::get_color_quantized_image(ImgVector<T>* decreased_color_image, const ImgVector<T>& image, const std::vector<std::list<VECTOR_2D<int> > >& connected_regions)
{
	decreased_color_image->reset(_width, _height);
	for (unsigned int n = 0; n < connected_regions.size(); n++) {
		T sum_color = T();
		for (std::list<VECTOR_2D<int> >::const_iterator ite = connected_regions[n].begin();
		    ite != connected_regions[n].end();
		    ++ite) {
			sum_color += image.get(ite->x, ite->y);
		}
		for (std::list<VECTOR_2D<int> >::const_iterator ite = connected_regions[n].begin();
		    ite != connected_regions[n].end();
		    ++ite) {
			decreased_color_image->at(ite->x, ite->y) = sum_color / double(connected_regions[n].size());
		}
	}
}

