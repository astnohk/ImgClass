#ifndef LIB_ImgClass_CrossCorrelation
#define LIB_ImgClass_CrossCorrelation

#include "ImgStatistics.h"

class CrossCorrelation
{
	private:
		int _width;
		int _height;
		ImgStatistics _img0;
		ImgStatistics _img1;
	public:
		CrossCorrelation(void);
		CrossCorrelation(const CrossCorrelation& copy);
		CrossCorrelation(const ImgStatistics& img0, const ImgStatistics& img1);
		~CrossCorrelation(void);
		CrossCorrelation& copy(const CrossCorrelation& copy);
		CrossCorrelation& operator=(const CrossCorrelation& copy);
		int width(void);
		int height(void);
		double NCC(const int x, const int y, const int window_width, const int window_height);
		double TruncatedNCC(const int x, const int y, const int window_width, const int window_height);
};

#endif

