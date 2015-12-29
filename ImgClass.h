#ifndef LIB_ImgClass
#define LIB_ImgClass

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
		int _width;
		int _height;

	public:
		ImgVector(void);
		ImgVector(const int W, const int H, const T& value = T());
		ImgVector(const int W, const int H, const T* array);
		explicit ImgVector(const ImgVector<T>& copy); // Copy constructor

		virtual ~ImgVector(void);

		void clear(void);
		void reset(const int W, const int H, const T& value = T());
		void reset(const int W, const int H, const T* array);
		void resize(const int W, const int H, const T& value = T());

		ImgVector<T>& copy(const ImgVector<T>& vector);
		template<class RT> ImgVector<T>& cast_copy(const ImgVector<RT>& vector);

		ImgVector<T>& operator=(const ImgVector<T>& vector);

		// Properties
		int width(void) const;
		int height(void) const;
		int size(void) const;
		bool isNULL(void) const;

		// Data access
		T* data(void) const;

		// Reference to the pixel
		T& operator[](const int n);
		T& at(const int n);
		T& at(const int x, const int y);
		T& at_repeat(const int x, const int y);
		T& at_mirror(const int x, const int y);

		// Get image intencity
		const T get(const int n) const; // return const to avoid to mistake get() for at()
		const T get(const int x, const int y) const;

		const T get_zeropad(const int x, const int y) const;
		const T get_repeat(const int x, const int y) const;
		const T get_mirror(const int x, const int y) const;

		const T get_zeropad_cubic(const double& x, const double& y, const double& B = 0.0, const double& C = (1.0 / 2.0)) const;
		const T get_repeat_cubic(const double& x, const double& y, const double& B = 0.0, const double& C = (1.0 / 2.0)) const;
		const T get_mirror_cubic(const double& x, const double& y, const double& B = 0.0, const double& C = (1.0 / 2.0)) const;

		const T min(void) const;
		const T max(void) const;
		const T min(const int top_left_x, const int top_left_y, const int crop_width, const int crop_height) const;
		const T max(const int top_left_x, const int top_left_y, const int crop_width, const int crop_height) const;

		const T variance(void) const;
		const T variance(const int top_left_x, const int top_left_y, const int crop_width, const int crop_height) const;

		// Cropping
		ImgVector<T>* crop(const int top_left_x, const int top_left_y, const int crop_width, const int crop_height) const;

		// Simple image processing
		void contrast_stretching(const T& Min, const T& Max);
		void map(T (*func)(T &value));

		// Resampling
		void resize_zerohold(const int W, const int H);
		void resize_bicubic(const int W, const int H, const bool saturate = false, const double min = 0.0, const double max = 0.0, T (*Nearest_Integer_Method)(double &d) = nullptr, const double B = (0.0 / 3.0), const double C = (1.0 / 2.0));

		// Operators
		template<class RT> ImgVector<T>& operator+=(const RT& val);
		template<class RT> ImgVector<T>& operator-=(const RT& val);
		template<class RT> ImgVector<T>& operator*=(const RT& val);
		template<class RT> ImgVector<T>& operator/=(const RT& val);

	protected:
		double cubic(const double x, const double B, const double C) const;
};

#include "ImgClass_private.h"

#endif

