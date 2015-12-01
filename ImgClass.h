#ifndef nullptr
#define nullptr 0
#endif



#ifndef LIB_ImgClass
#define LIB_ImgClass

template <class T>
class ImgVector
{
	private:
		T *_data;
		int _width;
		int _height;
	public:
		ImgVector(void);
		ImgVector(const ImgVector<T> &copy);
		ImgVector(int W, int H);
		ImgVector(int W, int H, const T &value);
		ImgVector(int W, int H, const T *array);
		virtual ~ImgVector(void);
		void reset(int W, int H);
		void reset(int W, int H, const T &value);
		void reset(int W, int H, const T *array);
		ImgVector<T>& copy(const ImgVector<T> &vector);
		ImgVector<T>& copy(const ImgVector<T> *vector);
		ImgVector<T>& operator=(const ImgVector<T> &vector);
		void set(int x, int y, const T &value);

		// Properties
		int width(void) const;
		int height(void) const;
		int size(void) const;
		bool isNULL(void) const;

		// Data access
		T* data(void) const;

		// Reference to the pixel
		T& operator[](int n);
		T& at(int x, int y);
		T& at_repeat(int x, int y);
		T& at_mirror(int x, int y);

		// Get image intencity
		const T get(int n) const; // return const to avoid to mistake get() for at()
		const T get(int x, int y) const;

		const T get_zeropad(int x, int y) const;
		const T get_repeat(int x, int y) const;
		const T get_mirror(int x, int y) const;

		const T get_zeropad(double x, double y, double B = 0.0, double C = (1.0 / 2.0)) const;
		const T get_repeat(double x, double y, double B = 0.0, double C = (1.0 / 2.0)) const;
		const T get_mirror(double x, double y, double B = 0.0, double C = (1.0 / 2.0)) const;

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
		//bool resize_bicubic(int W, int H, double min = 0.0, double max = 0.0, T (*Nearest_Integer_Method)(double &d) = nullptr, double B = (1.0 / 3.0), double C = (1.0 / 3.0));
		void resize_bicubic(int W, int H, double min = 0.0, double max = 0.0, T (*Nearest_Integer_Method)(double &d) = nullptr, double B = (0.0 / 3.0), double C = (1.0 / 2.0));

		double cubic(double x, double B, double C) const;
};

#include "ImgClass_private.h"

#endif

