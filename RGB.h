#ifndef LIB_ImgClass_RGB
#define LIB_ImgClass_RGB

namespace ImgClass {
	class RGB
	{
		public:
			double R;
			double G;
			double B;

			// Constructor
			RGB(void);
			RGB(const double& red, const double& green, const double& blue, const double& gamma = 0.0);
			RGB(const RGB& color); // Copy constructor

			RGB& set(const double& red, const double& green, const double& blue, const double& gamma = 0.0);

			// Operators
			explicit operator double() const; // return intensity

			RGB& operator=(const RGB& rvalue);
			RGB& operator=(const double& rvalue);

			RGB& operator+=(const RGB& rcolor);
			RGB& operator+=(const double& rvalue);

			RGB& operator-=(const RGB& rcolor);
			RGB& operator-=(const double& rvalue);

			RGB& operator*=(const RGB& rcolor);
			RGB& operator*=(const double& rvalue);

			RGB& operator/=(const RGB& rcolor);
			RGB& operator/=(const double& rvalue);

			// Non substituting
			const RGB operator+(const RGB& rcolor) const;
			const RGB operator+(const double& rvalue) const;

			const RGB operator-(const RGB& rcolor) const;
			const RGB operator-(const double& rvalue) const;
	};
}


const ImgClass::RGB operator+(ImgClass::RGB rcolor);
const ImgClass::RGB operator-(ImgClass::RGB rcolor);

const ImgClass::RGB operator*(const ImgClass::RGB& lcolor, const ImgClass::RGB& rcolor);
const ImgClass::RGB operator*(const ImgClass::RGB& lcolor, const double& rvalue);
const ImgClass::RGB operator*(const double& lvalue, const ImgClass::RGB& rcolor);

const ImgClass::RGB operator/(const ImgClass::RGB& lcolor, const ImgClass::RGB& rcolor);
const ImgClass::RGB operator/(const ImgClass::RGB& lcolor, const double& rvalue);
const ImgClass::RGB operator/(const double& lvalue, const ImgClass::RGB& rcolor);

// Comparator
bool operator==(const ImgClass::RGB& lcolor, const ImgClass::RGB& rcolor);
bool operator!=(const ImgClass::RGB& lcolor, const ImgClass::RGB& rcolor);

#endif

