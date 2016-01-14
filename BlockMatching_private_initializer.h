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
BlockMatching<T>::BlockMatching(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const int BlockSize)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0;
	if (image_prev.isNULL()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : ImgVector<T>& image_prev" << std::endl;
		throw std::invalid_argument("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : ImgVector<T>& image_prev");
	} else if (image_current.isNULL()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : ImgVector<T>& image_current" << std::endl;
		throw std::invalid_argument("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : ImgVector<T>& image_current");
	} else if (image_prev.width() != image_current.width()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : width of image_prev and image_current not match" << std::endl;
		throw std::invalid_argument("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : width of image_prev and image_current not match");
	} else if (image_prev.height() != image_current.height()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : height of image_prev and image_current not match" << std::endl;
		throw std::invalid_argument("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : height of image_prev and image_current not match");
	} else if (BlockSize < 0) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : BlockSize" << std::endl;
		throw std::out_of_range("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : BlockSize");
	}

	_width = image_prev.width();
	_height = image_prev.height();
	_block_size = BlockSize;
	_cells_width = int(ceil(double(_width) / double(_block_size)));
	_cells_height = int(ceil(double(_height) / double(_block_size)));

	_image_prev.copy(image_prev);
	_image_current.copy(image_current);

	// Normalize the image
	image_normalizer();
}

template <class T>
BlockMatching<T>::BlockMatching(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<T>& image_next, const int BlockSize)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0;
	if (image_prev.isNULL()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : ImgVector<T>& image_prev" << std::endl;
		throw std::invalid_argument("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : ImgVector<T>& image_prev");
	} else if (image_current.isNULL()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : ImgVector<T>& image_current" << std::endl;
		throw std::invalid_argument("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : ImgVector<T>& image_current");
	} else if (image_prev.width() != image_current.width()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : width of image_prev and image_current not match" << std::endl;
		throw std::invalid_argument("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : width of image_prev and image_current not match");
	} else if (image_prev.height() != image_current.height()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : height of image_prev and image_current not match" << std::endl;
		throw std::invalid_argument("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : height of image_prev and image_current not match");
	} else if (BlockSize < 0) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : BlockSize" << std::endl;
		throw std::out_of_range("BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : BlockSize");
	}

	_width = image_prev.width();
	_height = image_prev.height();
	_block_size = BlockSize;
	_cells_width = int(ceil(double(_width) / double(_block_size)));
	_cells_height = int(ceil(double(_height) / double(_block_size)));

	_image_prev.copy(image_prev);
	_image_current.copy(image_current);
	_image_next.copy(image_next);

	// Normalize the image
	image_normalizer();
}

template <class T>
BlockMatching<T>::BlockMatching(const ImgVector<T>& image_prev, const ImgVector<size_t>& region_map_prev, const ImgVector<T>& image_current, const ImgVector<size_t>& region_map_current)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0;
	if (image_prev.isNULL()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : const ImgVector<T>& image_prev" << std::endl;
		throw std::invalid_argument("const ImgVector<T>& image_prev");
	} else if (image_current.isNULL()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : const ImgVector<T>& image_current" << std::endl;
		throw std::invalid_argument("const ImgVector<T>& image_current");
	} else if (image_prev.width() != image_current.width() || image_prev.height() != image_current.height()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : width or height of image_prev and image_current not match" << std::endl;
		throw std::invalid_argument("width or height of image_prev and image_current not match");
	} else if (region_map_prev.isNULL()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : const ImgVector<int>& region_map_prev" << std::endl;
		throw std::invalid_argument("const ImgVector<int>& region_map_prev");
	} else if (region_map_current.isNULL()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : const ImgVector<int>& region_map_current" << std::endl;
		throw std::invalid_argument("const ImgVector<int>& region_map_prev");
	} else if (region_map_prev.width() != image_prev.width() || region_map_prev.height() != image_prev.height()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : width or height of region_map_prev and image_prev not match" << std::endl;
		throw std::invalid_argument("width or height of region_map_prev and image_prev not match");
	} else if (region_map_current.width() != image_current.width() || region_map_current.height() != image_current.height()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : width or height of region_map_current and image_current not match" << std::endl;
		throw std::invalid_argument("width or height of region_map_current and image_current not match");
	}

	_width = image_prev.width();
	_height = image_prev.height();
	_block_size = 1;
	_cells_width = _width;
	_cells_height = _height;
	_image_prev.copy(image_prev);
	_image_current.copy(image_current);
	_region_map_prev.copy(region_map_prev);
	_region_map_current.copy(region_map_current);
	// Normalize the image
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_BLOCKMATCHING)
	std::cout << " Block Matching : Normalize the input images" << std::endl;
#endif
	image_normalizer();
	// Extract connected regions from region_map
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_BLOCKMATCHING)
	std::cout << " Block Matching : Collect connected region from region map" << std::endl;
#endif
	get_connected_region_list(&_connected_regions_prev, region_map_prev);
	get_connected_region_list(&_connected_regions_current, region_map_current);
	// Get decreased image
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_BLOCKMATCHING)
	std::cout << " Block Matching : Get color quantized image" << std::endl;
#endif
	get_color_quantized_image(&_color_quantized_prev, _image_prev, _connected_regions_prev);
	get_color_quantized_image(&_color_quantized_current, _image_current, _connected_regions_current);
}

template <class T>
BlockMatching<T>::BlockMatching(const ImgVector<T>& image_prev, const ImgVector<size_t>& region_map_prev, const ImgVector<T>& image_current, const ImgVector<size_t>& region_map_current, const ImgVector<T>& image_next, const ImgVector<size_t>& region_map_next)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0;
	if (image_prev.isNULL()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : const ImgVector<T>& image_prev" << std::endl;
		throw std::invalid_argument("const ImgVector<T>& image_prev");
	} else if (image_current.isNULL()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : const ImgVector<T>& image_current" << std::endl;
		throw std::invalid_argument("const ImgVector<T>& image_current");
	} else if (image_prev.width() != image_current.width() || image_prev.height() != image_current.height()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : width or height of image_prev and image_current not match" << std::endl;
		throw std::invalid_argument("width or height of image_prev and image_current not match");
	} else if (region_map_prev.isNULL()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : const ImgVector<int>& region_map_prev" << std::endl;
		throw std::invalid_argument("const ImgVector<int>& region_map_prev");
	} else if (region_map_current.isNULL()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : const ImgVector<int>& region_map_current" << std::endl;
		throw std::invalid_argument("const ImgVector<int>& region_map_prev");
	} else if (region_map_prev.width() != image_prev.width() || region_map_prev.height() != image_prev.height()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : width or height of region_map_prev and image_prev not match" << std::endl;
		throw std::invalid_argument("width or height of region_map_prev and image_prev not match");
	} else if (region_map_current.width() != image_current.width() || region_map_current.height() != image_current.height()) {
		std::cerr << "BlockMatching<T>::BlockMatching(const ImgVector<T>&, const ImgVector<T>&, const int) : width or height of region_map_current and image_current not match" << std::endl;
		throw std::invalid_argument("width or height of region_map_current and image_current not match");
	}

	_width = image_prev.width();
	_height = image_prev.height();
	_block_size = 1;
	_cells_width = _width;
	_cells_height = _height;

	_image_prev.copy(image_prev);
	_image_current.copy(image_current);
	_image_next.copy(image_next);
	_region_map_prev.copy(region_map_prev);
	_region_map_current.copy(region_map_current);
	_region_map_next.copy(region_map_next);

	// Normalize the image
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_BLOCKMATCHING)
	std::cout << " Block Matching : Normalize the input images" << std::endl;
#endif
	image_normalizer();
	// Extract connected regions from region_map
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_BLOCKMATCHING)
	std::cout << " Block Matching : Collect connected region from region map" << std::endl;
#endif
	get_connected_region_list(&_connected_regions_prev, region_map_prev);
	get_connected_region_list(&_connected_regions_current, region_map_current);
	get_connected_region_list(&_connected_regions_next, region_map_next);
	// Get decreased image
#if defined(OUTPUT_IMG_CLASS) || defined(OUTPUT_IMG_CLASS_BLOCKMATCHING)
	std::cout << " Block Matching : Get color quantized image" << std::endl;
#endif
	get_color_quantized_image(&_color_quantized_prev, _image_prev, _connected_regions_prev);
	get_color_quantized_image(&_color_quantized_current, _image_current, _connected_regions_current);
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
	_image_current.copy(copy._image_current);
	_image_next.copy(copy._image_next);

	_region_map_prev.copy(copy._region_map_prev);
	_region_map_current.copy(copy._region_map_current);
	_region_map_next.copy(copy._region_map_next);
	_connected_regions_prev.assign(copy._connected_regions_prev.begin(), copy._connected_regions_prev.end());
	_connected_regions_current.assign(copy._connected_regions_current.begin(), copy._connected_regions_current.end());
	_connected_regions_next.assign(copy._connected_regions_next.begin(), copy._connected_regions_next.end());

	_motion_vector_time.copy(copy._motion_vector_time);
	_motion_vector_prev.copy(copy._motion_vector_prev);
	_motion_vector_next.copy(copy._motion_vector_next);
}




template <class T>
BlockMatching<T>::~BlockMatching(void)
{
}




template <class T>
void
BlockMatching<T>::reset(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const int BlockSize)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0;
	if (image_prev.isNULL()) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const int) : const ImgVector<T>& image_prev" << std::endl;
		throw std::invalid_argument("const ImgVector<T>& image_prev");
	} else if (image_prev.width() != image_current.width() || image_prev.height() != image_current.height()) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const int) : const ImgVector<T>& image_prev, const ImgVector<T>& image_prev" << std::endl;
		throw std::invalid_argument("width or height of image_prev and image_current not match");
	} else if (BlockSize < 0) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const int) : const int BlockSize" << std::endl;
		throw std::out_of_range("BlockSize");
	}

	_width = image_prev.width();
	_height = image_prev.height();
	_block_size = BlockSize;
	_cells_width = int(ceil(double(_width) / double(_block_size)));
	_cells_height = int(ceil(double(_height) / double(_block_size)));

	_image_prev.copy(image_prev);
	_image_current.copy(image_current);
	_image_next.clear();

	_region_map_prev.clear();
	_region_map_current.clear();
	_region_map_next.clear();
	_connected_regions_prev.clear();
	_connected_regions_current.clear();
	_connected_regions_next.clear();

	_motion_vector_time.clear();
	_motion_vector_prev.clear();
	_motion_vector_next.clear();

	// Normalize the image
	image_normalizer();
}

template <class T>
void
BlockMatching<T>::reset(const ImgVector<T>& image_prev, const ImgVector<T>& image_current, const ImgVector<T>& image_next, const int BlockSize)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0;
	if (image_prev.isNULL()) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const int) : const ImgVector<T>& image_prev" << std::endl;
		throw std::invalid_argument("const ImgVector<T>& image_prev");
	} else if (image_prev.width() != image_current.width() || image_prev.height() != image_current.height()) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const int) : const ImgVector<T>& image_prev, const ImgVector<T>& image_prev" << std::endl;
		throw std::invalid_argument("width or height of image_prev and image_current not match");
	} else if (BlockSize < 0) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const int) : const int BlockSize" << std::endl;
		throw std::out_of_range("BlockSize");
	}

	_width = image_prev.width();
	_height = image_prev.height();
	_block_size = BlockSize;
	_cells_width = int(ceil(double(_width) / double(_block_size)));
	_cells_height = int(ceil(double(_height) / double(_block_size)));

	_image_prev.copy(image_prev);
	_image_current.copy(image_current);
	_image_next.copy(image_next);

	_region_map_prev.clear();
	_region_map_current.clear();
	_region_map_next.clear();
	_connected_regions_prev.clear();
	_connected_regions_current.clear();
	_connected_regions_next.clear();

	_motion_vector_time.clear();
	_motion_vector_prev.clear();
	_motion_vector_next.clear();

	// Normalize the image
	image_normalizer();
}

template <class T>
void
BlockMatching<T>::reset(const ImgVector<T>& image_prev, const ImgVector<size_t>& region_map_prev, const ImgVector<T>& image_current, const ImgVector<size_t>& region_map_current)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0;
	if (image_prev.isNULL()) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const ImgVector<int>&) : const ImgVector<T>& image_prev" << std::endl;
		throw std::invalid_argument("const ImgVector<T>& image_prev");
	} else if (image_prev.width() != image_current.width() || image_prev.height() != image_current.height()) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const ImgVector<int>&) : const ImgVector<T>& image_prev, const ImgVector<T>& image_current" << std::endl;
		throw std::invalid_argument("width or height of image_prev and image_current not match");
	} else if (region_map_prev.width() != image_prev.width() || region_map_prev.height() != image_prev.height()) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const ImgVector<int>&) : const ImgVector<int>& region_map_prev" << std::endl;
		throw std::invalid_argument("width or height of region_map_prev not match with image_prev");
	} else if (region_map_current.width() != image_current.width() || region_map_current.height() != image_current.height()) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const ImgVector<int>&) : const ImgVector<int>& region_map_current" << std::endl;
		throw std::invalid_argument("width or height of region_map_current not match with image_current");
	}

	_width = image_prev.width();
	_height = image_prev.height();
	_block_size = 1;
	_cells_width = _width;
	_cells_height = _height;

	_image_prev.copy(image_prev);
	_image_current.copy(image_current);
	_image_next.clear();

	_region_map_prev.copy(region_map_prev);
	_region_map_current.copy(region_map_current);
	_region_map_next.clear();

	_motion_vector_time.clear();
	_motion_vector_prev.clear();
	_motion_vector_next.clear();
	// Normalize the image
	image_normalizer();
	// Extract connected regions from region_map
	get_connected_region_list(&_connected_regions_prev, region_map_prev);
	get_connected_region_list(&_connected_regions_current, region_map_current);
	_connected_regions_next.clear();
	// Get decreased image
	get_color_quantized_image(&_color_quantized_prev, _image_prev, _connected_regions_prev);
	get_color_quantized_image(&_color_quantized_current, _image_current, _connected_regions_current);
	_color_quantized_next.clear();
}

template <class T>
void
BlockMatching<T>::reset(const ImgVector<T>& image_prev, const ImgVector<size_t>& region_map_prev, const ImgVector<T>& image_current, const ImgVector<size_t>& region_map_current, const ImgVector<T>& image_next, const ImgVector<size_t>& region_map_next)
{
	_width = 0;
	_height = 0;
	_block_size = 0;
	_cells_width = 0;
	_cells_height = 0;
	if (image_prev.isNULL()) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const ImgVector<int>&) : const ImgVector<T>& image_prev" << std::endl;
		throw std::invalid_argument("const ImgVector<T>& image_prev");
	} else if (image_prev.width() != image_current.width() || image_prev.height() != image_current.height()) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const ImgVector<int>&) : const ImgVector<T>& image_prev, const ImgVector<T>& image_current" << std::endl;
		throw std::invalid_argument("width or height of image_prev and image_current not match");
	} else if (region_map_prev.width() != image_prev.width() || region_map_prev.height() != image_prev.height()) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const ImgVector<int>&) : const ImgVector<int>& region_map_prev" << std::endl;
		throw std::invalid_argument("width or height of region_map_prev not match with image_prev");
	} else if (region_map_current.width() != image_current.width() || region_map_current.height() != image_current.height()) {
		std::cerr << "void BlockMatching<T>::reset(const ImgVector<T>&, const ImgVector<T>&, const ImgVector<int>&) : const ImgVector<int>& region_map_current" << std::endl;
		throw std::invalid_argument("width or height of region_map_current not match with image_current");
	}

	_width = image_prev.width();
	_height = image_prev.height();
	_block_size = 1;
	_cells_width = _width;
	_cells_height = _height;

	_image_prev.copy(image_prev);
	_image_current.copy(image_current);
	_image_next.copy(image_next);

	_region_map_prev.copy(region_map_prev);
	_region_map_current.copy(region_map_current);
	_region_map_next.copy(region_map_next);

	_motion_vector_time.clear();
	_motion_vector_prev.clear();
	_motion_vector_next.clear();
	// Normalize the image
	image_normalizer();
	// Extract connected regions from region_map
	get_connected_region_list(&_connected_regions_prev, region_map_prev);
	get_connected_region_list(&_connected_regions_current, region_map_current);
	get_connected_region_list(&_connected_regions_next, region_map_next);
	// Get decreased image
	get_color_quantized_image(&_color_quantized_prev, _image_prev, _connected_regions_prev);
	get_color_quantized_image(&_color_quantized_current, _image_current, _connected_regions_current);
	get_color_quantized_image(&_color_quantized_next, _image_next, _connected_regions_next);
}




/* Get connected region list
 *
 * If region_map[i] < 0 then it means the pixel is neglected.
 * So usually region_map include only the integer n > 0.
 */
template <class T>
void
BlockMatching<T>::get_connected_region_list(std::vector<std::list<VECTOR_2D<int> > >* connected_regions, const ImgVector<size_t>& region_map)
{
	const VECTOR_2D<int> adjacent[8] = {
	    VECTOR_2D<int>(-1, -1), VECTOR_2D<int>(0, -1), VECTOR_2D<int>(1, -1),
	    VECTOR_2D<int>(-1, 0), VECTOR_2D<int>(1, 0),
	    VECTOR_2D<int>(-1, 1), VECTOR_2D<int>(0, 1), VECTOR_2D<int>(1, 1)};
	std::list<std::list<VECTOR_2D<int> > > tmp_list;
	ImgVector<bool> collected(_width, _height, false);
	// Clear the vector
	connected_regions->clear();

	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			if (collected.get(x, y) == false) {
				size_t num = region_map.get(x, y);
				collected.at(x, y) = true;
				tmp_list.push_back(std::list<VECTOR_2D<int> >(0)); // Add new region pixel list
				VECTOR_2D<int> r(x, y);
				tmp_list.back().push_back(r); // Add first element
				for (std::list<VECTOR_2D<int> >::const_iterator ite = tmp_list.back().begin();
				    ite != tmp_list.back().end();
				    ++ite) {
					for (int k = 0; k < 8; k++) {
						r.x = ite->x + adjacent[k].x;
						r.y = ite->y + adjacent[k].y;
						if (0 <= r.x && r.x < _width && 0 <= r.y && r.y < _height
						    && collected.get(r.x, r.y) == false
						    && region_map.get(r.x, r.y) == num) {
							collected.at(r.x, r.y) = true;
							tmp_list.back().push_back(r);
						}
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
	double max_int = _image_prev.max();
	if (max_int > 1.0) {
		_image_prev /= max_int;
	}
	max_int = _image_current.max();
	if (max_int > 1.0) {
		_image_current /= max_int;
	}
	max_int = _image_next.max();
	if (max_int > 1.0) {
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

