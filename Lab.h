#ifndef nullptr
#define nullptr 0
#endif



#ifndef LIB_ImgClass_Lab
#define LIB_ImgClass_Lab

/* L*a*b* color space
 * constructor will map RGB to L*a*b* color space
 */
namespace ImgClass {
	class Lab
	{
		public:
			double L;
			double a;
			double b;

			// Constructor
			Lab(void);
			Lab(const double& red, const double& green, const double& blue);
			Lab(const double& value);

			// Operators
			Lab& operator=(const Lab& value);
			Lab& operator=(const double& value);

			Lab& operator+=(const Lab& color);
			Lab& operator+=(const double& value);

			Lab& operator-=(const Lab& color);
			Lab& operator-=(const double& value);

			Lab& operator*=(const Lab& color);
			Lab& operator*=(const double& value);

			Lab& operator/=(const Lab& color);
			Lab& operator/=(const double& value);

			// friend Global Operators
			// Arithmetic
			friend Lab& operator+(Lab color);
			friend Lab& operator-(Lab color);

			friend Lab& operator+(Lab lcolor, const Lab& rcolor);
			friend Lab& operator+(Lab lcolor, const double& rvalue);
			friend Lab& operator+(const double& lvalue, Lab rcolor);

			friend Lab& operator-(Lab lcolor, const Lab& rcolor);
			friend Lab& operator-(Lab lcolor, const double& rvalue);
			friend Lab& operator-(const double& lvalue, Lab rcolor);

			friend Lab& operator*(Lab lcolor, const Lab& rcolor);
			friend Lab& operator*(Lab lcolor, const double& rvalue);
			friend Lab& operator*(const double& lvalue, Lab rcolor);

			friend Lab& operator/(Lab lcolor, const Lab& rcolor);
			friend Lab& operator/(Lab lcolor, const double& rvalue);
			friend Lab& operator/(const double& lvalue, Lab rcolor);

			// Comparator
			friend bool operator==(const Lab& lcolor, const Lab& rcolor);
			friend bool operator!=(const Lab& lcolor, const Lab& rcolor);
	};
}

// Global Operators
// Arithmetic
ImgClass::Lab& operator+(ImgClass::Lab color);
ImgClass::Lab& operator-(ImgClass::Lab color);

ImgClass::Lab& operator+(ImgClass::Lab lcolor, const ImgClass::Lab& rcolor);
ImgClass::Lab& operator+(ImgClass::Lab lcolor, const double& rvalue);
ImgClass::Lab& operator+(const double& lvalue, ImgClass::Lab rcolor);

ImgClass::Lab& operator-(ImgClass::Lab lcolor, const ImgClass::Lab& rcolor);
ImgClass::Lab& operator-(ImgClass::Lab lcolor, const double& rvalue);
ImgClass::Lab& operator-(const double& lvalue, ImgClass::Lab rcolor);

ImgClass::Lab& operator*(ImgClass::Lab lcolor, const ImgClass::Lab& rcolor);
ImgClass::Lab& operator*(ImgClass::Lab lcolor, const double& rvalue);
ImgClass::Lab& operator*(const double& lvalue, ImgClass::Lab rcolor);

ImgClass::Lab& operator/(ImgClass::Lab lcolor, const ImgClass::Lab& rcolor);
ImgClass::Lab& operator/(ImgClass::Lab lcolor, const double& rvalue);
ImgClass::Lab& operator/(const double& lvalue, ImgClass::Lab rcolor);

// Comparator
bool operator==(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor);
bool operator!=(const ImgClass::Lab& lcolor, const ImgClass::Lab& rcolor);

#endif

