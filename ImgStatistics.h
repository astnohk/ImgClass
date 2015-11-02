#ifndef nullptr
#define nullptr 0
#endif


#ifndef LIB_ImgStatistics
#define LIB_ImgStatistics

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

