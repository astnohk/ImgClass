#ifndef LIB_ImgClass_ImgStatistics
#define LIB_ImgClass_ImgStatistics

/* Macro for compatibility where the C++11 is not supported
#ifndef nullptr
#define nullptr 0
#endif
*/


class ImgStatistics
{
	private:
		int _width;
		int _height;
		double *_data;
	public:
		ImgStatistics(void);
		ImgStatistics(const ImgStatistics &copy);
		ImgStatistics(int W, int H, double *Img);

		virtual ~ImgStatistics(void);

		void set(int W, int H, double *Img);
		ImgStatistics& copy(const ImgStatistics &copy);
		ImgStatistics& operator=(const ImgStatistics &copy);

		double& image(int x, int y) const;

		int width(void) const;
		int height(void) const;

		double mean() const;
		double mean(int start_x, int start_y, int end_x, int end_y) const;
		double variance() const;
		double variance(int x, int y, int window_width, int window_height) const;
		double std_deviation(void) const;
		double std_deviation(int x, int y, int window_width, int window_height) const;
};


class Histogram
{
	private:
		int _bins;
		double *_hist;
	public:
		Histogram(void);
		Histogram(const Histogram &copy);
		explicit Histogram(int init_bins);
		Histogram& copy(const Histogram &copy);
		Histogram& reset(int init_bins);
		~Histogram(void);
		void free(void);

		// Read
		const double* data(void) const;
		int bins(void) const;
		double get(int bin) const;

		// Control
		bool add(int bin, double val);
};

#endif

