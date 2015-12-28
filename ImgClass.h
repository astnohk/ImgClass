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
		ImgVector(int W, int H, const T& value = T());
		ImgVector(int W, int H, const T* array);
		explicit ImgVector(const ImgVector<T>& copy); // Copy constructor

		virtual ~ImgVector(void);

		void clear(void);
		void reset(int W, int H, const T& value = T());
		void reset(int W, int H, const T* array);

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
		T& operator[](int n);
		T& at(int n);
		T& at(int x, int y);
		T& at_repeat(int x, int y);
		T& at_mirror(int x, int y);

		// Get image intencity
		const T get(int n) const; // return const to avoid to mistake get() for at()
		const T get(int x, int y) const;

		const T get_zeropad(int x, int y) const;
		const T get_repeat(int x, int y) const;
		const T get_mirror(int x, int y) const;

		const T get_zeropad_cubic(double x, double y, double B = 0.0, double C = (1.0 / 2.0)) const;
		const T get_repeat_cubic(double x, double y, double B = 0.0, double C = (1.0 / 2.0)) const;
		const T get_mirror_cubic(double x, double y, double B = 0.0, double C = (1.0 / 2.0)) const;

		const T min(void) const;
		const T max(void) const;
		const T min(int top_left_x, int top_left_y, int crop_width, int crop_height) const;
		const T max(int top_left_x, int top_left_y, int crop_width, int crop_height) const;

		const T variance(void) const;
		const T variance(const int top_left_x, const int top_left_y, const int crop_width, const int crop_height) const;

		// Cropping
		ImgVector<T>* crop(int top_left_x, int top_left_y, int crop_width, int crop_height) const;

		// Simple image processing
		void contrast_stretching(const T& Min, const T& Max);
		void map(T (*func)(T &value));

		// Resampling
		void resize_zerohold(int W, int H);
		void resize_bicubic(int W, int H, bool saturate = false, double min = 0.0, double max = 0.0, T (*Nearest_Integer_Method)(double &d) = nullptr, double B = (0.0 / 3.0), double C = (1.0 / 2.0));

		double cubic(double x, double B, double C) const;

		// Operators
		template<class RT> ImgVector<T>& operator+=(const RT& val);
		template<class RT> ImgVector<T>& operator-=(const RT& val);
		template<class RT> ImgVector<T>& operator*=(const RT& val);
		template<class RT> ImgVector<T>& operator/=(const RT& val);
};

#include "ImgClass_private.h"

#endif

