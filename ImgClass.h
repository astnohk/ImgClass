#ifndef LIB_ImgClass
#define LIB_ImgClass

#include <cfloat>
#include <cstddef>
/* Macro for compatibility where the C++11 not supported
#ifndef nullptr
#define nullptr 0
#endif
*/

template <class T>
class ImgVector
{
	private:
		T *_data;
		size_t _reserved_size;
		int _width;
		int _height;

	public:
		ImgVector(void);
		ImgVector(const int Width, const int Height, const T& value = T());
		ImgVector(const int Width, const int Height, const T* array);
		explicit ImgVector(const ImgVector<T>& copy); // Copy constructor

		virtual ~ImgVector(void);

		void clear(void);
		void reserve(const int Width, const int Height);
		void reset(const int Width, const int Height, const T& value = T()); // Delete current data and resize the array
		void reset(const int Width, const int Height, const T* array); // Delete current data and resize the array
		void resize(const int Width, const int Height, const T& value = T()); // Do only resizing

		ImgVector<T>& copy(const ImgVector<T>& vector); // Assign vector to *this
		ImgVector<T>& operator=(const ImgVector<T>& vector); // Assign vector to *this
		template<class RT> ImgVector<T>& cast_copy(const ImgVector<RT>& vector);

		// Get Properties
		size_t reserved_size(void) const;
		int width(void) const;
		int height(void) const;
		size_t size(void) const;
		bool isNULL(void) const;

		// Data access
		T* data(void) const;
		// Reference to the pixel
		T& operator[](const size_t n);
		T& at(const size_t n);
		T& at(const int x, const int y);
		T& at_repeat(const int x, const int y);
		T& at_mirror(const int x, const int y);
		// Get image intencity
		const T get(const size_t n) const; // return const to avoid to mistake get() for at()
		const T get(const int x, const int y) const;
		// Get intencity with boundary treatment
		const T get_zeropad(const int x, const int y) const;
		const T get_repeat(const int x, const int y) const;
		const T get_mirror(const int x, const int y) const;
		// Get intensity interpolated by bicubic
		const T get_zeropad_cubic(const double& x, const double& y, const double& B = 0.0, const double& C = (1.0 / 2.0)) const;
		const T get_repeat_cubic(const double& x, const double& y, const double& B = 0.0, const double& C = (1.0 / 2.0)) const;
		const T get_mirror_cubic(const double& x, const double& y, const double& B = 0.0, const double& C = (1.0 / 2.0)) const;

		// Get statistical value
		const T min(void) const;
		const T max(void) const;
		const T min(const int top_left_x, const int top_left_y, const int crop_width, const int crop_height) const;
		const T max(const int top_left_x, const int top_left_y, const int crop_width, const int crop_height) const;
		const T variance(void) const;
		const T variance(const int top_left_x, const int top_left_y, const int crop_width, const int crop_height) const;

		// Cropping (Change the data of *this)
		ImgVector<T>* crop(const int top_left_x, const int top_left_y, const int crop_width, const int crop_height) const;

		// Simple image processing (Change the data of *this)
		void contrast_stretching(const T& Min, const T& Max);
		void map(T (*func)(T &value));

		// Resampling (Change the data of *this)
		void resample_zerohold(const int Width, const int Height);
		void resample_bicubic(const int Width, const int Height, const bool saturate = false, const double min = 0.0, const double max = 0.0, T (*Nearest_Integer_Method)(double &d) = nullptr, const double B = (0.0 / 3.0), const double C = (1.0 / 2.0));

		// Operators
		template<class RT> ImgVector<T>& operator+=(const RT& rvalue);
		template<class RT> ImgVector<T>& operator-=(const RT& rvalue);
		template<class RT> ImgVector<T>& operator*=(const RT& rvalue);
		template<class RT> ImgVector<T>& operator/=(const RT& rvalue);

		template<class RT> ImgVector<T>& operator+=(const ImgVector<RT>& rvector);
		template<class RT> ImgVector<T>& operator-=(const ImgVector<RT>& rvector);
		template<class RT> ImgVector<T>& operator*=(const ImgVector<RT>& rvector);
		template<class RT> ImgVector<T>& operator/=(const ImgVector<RT>& rvector);

	protected:
		double cubic(const double x, const double B, const double C) const;
};

#include "ImgClass_private.h"

#endif

