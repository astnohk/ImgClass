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
		explicit CrossCorrelation(const CrossCorrelation &copy);
		CrossCorrelation(ImgStatistics &img0, ImgStatistics &img1);
		~CrossCorrelation(void);
		CrossCorrelation& copy(const CrossCorrelation &copy);
		CrossCorrelation& operator=(const CrossCorrelation &copy);
		int width(void);
		int height(void);
		double NCC(int x, int y, int window_width, int window_height);
		double TruncatedNCC(int x, int y, int window_width, int window_height);
};

#endif

